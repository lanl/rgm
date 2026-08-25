#!/usr/bin/env python3
"""
Generate the iso_c_binding shim module (src/module_geological_model_c.f90)
by parsing the parameter components of the rgm2_curved and rgm3_curved
derived types. Rerun this script after adding/removing type parameters:

    python3 python/generate_shim.py

The Fortran types remain the single source of truth; the shim maps
parameter names (strings) to type components.
"""

import re
import os

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Output arrays exposed through getters (name -> rank); present in both
# types unless listed in only3d/only2d
outputs_common = {
    'image': None, 'rgt': None, 'facies': None, 'fault': None,
    'fault_dip': None, 'fault_disp': None, 'salt': None, 'karst': None,
    'vp': None, 'vs': None, 'rho': None, 'psf': None,
    'image_pp': None, 'image_ps': None, 'image_sp': None, 'image_ss': None,
}
outputs_only3d = {'fault_strike': None, 'fault_rake': None}

# Components that are inputs but not settable through the generic numeric
# interface (multi-dimensional custom input arrays)
skip = {'refl', 'refl_top', 'psf', 'image', 'rgt', 'facies', 'fault',
        'fault_dip', 'fault_strike', 'fault_rake', 'fault_disp', 'salt',
        'karst', 'vp', 'vs', 'rho', 'image_pp', 'image_ps', 'image_sp',
        'image_ss'}


def parse_type(src, type_name):
    """Extract (name, kind) parameter list from a derived type definition.
    kind is one of: real, int, logical, string, real1 (fixed or allocatable
    1D real array)."""
    m = re.search(r'type %s\n(.*?)\n    contains' % type_name, src, re.S)
    body = m.group(1)
    params = []
    for line in body.split('\n'):
        line = line.strip()
        if not line or line.startswith('!'):
            continue
        dm = re.match(r'(real|integer|logical|character)\s*(\(len=\d+\))?\s*(,\s*dimension\s*\(([^)]*)\))?\s*(,\s*allocatable)?(,\s*dimension\s*\(([^)]*)\))?\s*::\s*(.*)', line)
        if not dm:
            continue
        ftype = dm.group(1)
        dims = dm.group(4) or dm.group(7)
        names = dm.group(8).split('!')[0]
        # strip default initializers
        names = re.sub(r'=\s*\[[^\]]*\]', '', names)
        names = re.sub(r'=[^,]*', '', names)
        for name in [n.strip() for n in names.split(',') if n.strip()]:
            if name in skip:
                continue
            if ftype == 'character':
                params.append((name, 'string'))
            elif ftype == 'logical':
                params.append((name, 'logical'))
            elif ftype == 'integer':
                params.append((name, 'int'))
            elif ftype == 'real' and dims is None and 'allocatable' not in line:
                params.append((name, 'real'))
            elif ftype == 'real':
                params.append((name, 'real1'))
    return params


def set_cases(params, var):
    out = []
    for name, kind in params:
        if kind == 'real':
            out.append(f"            case ('{name}')\n                {var}%{name} = real(val)")
        elif kind == 'int':
            out.append(f"            case ('{name}')\n                {var}%{name} = nint(val)")
        elif kind == 'logical':
            out.append(f"            case ('{name}')\n                {var}%{name} = nint(val) /= 0")
    return '\n'.join(out)


def set_array_cases(params, var):
    out = []
    for name, kind in params:
        if kind == 'real1':
            out.append(f"            case ('{name}')\n                {var}%{name} = real(vals)")
    return '\n'.join(out)


def set_string_cases(params, var):
    out = []
    for name, kind in params:
        if kind == 'string':
            out.append(f"            case ('{name}')\n                {var}%{name} = trim(str)")
    return '\n'.join(out)


def get_shape_cases(outputs, var, rank):
    out = []
    for name in outputs:
        dims_set = '; '.join(f'dims({i + 1}) = size({var}%{name}, {i + 1})' for i in range(rank))
        out.append(f"""            case ('{name}')
                if (allocated({var}%{name})) then
                    {dims_set}
                    ok = 1
                end if""")
    return '\n'.join(out)


def get_cases(outputs, var, rank):
    out = []
    shape_expr = ', '.join(f'size(p, {i + 1})' for i in range(rank))
    for name in outputs:
        out.append(f"""            case ('{name}')
                if (allocated({var}%{name})) then
                    call copy_out_{rank}d({var}%{name}, dims, data_ptr)
                    ok = 1
                end if""")
    return '\n'.join(out)


src2 = open(os.path.join(root, 'src/module_geological_model_2d_curved.f90')).read()
src3 = open(os.path.join(root, 'src/module_geological_model_3d_curved.f90')).read()
p2 = parse_type(src2, 'rgm2_curved')
p3 = parse_type(src3, 'rgm3_curved')

out2 = dict(outputs_common)
out3 = dict(outputs_common)
out3.update(outputs_only3d)

template = open(os.path.join(root, 'python/shim_template.f90')).read()
code = template
code = code.replace('!SET_CASES_2D', set_cases(p2, 'pool2(h)%p'))
code = code.replace('!SET_CASES_3D', set_cases(p3, 'pool3(h)%p'))
code = code.replace('!SET_ARRAY_CASES_2D', set_array_cases(p2, 'pool2(h)%p'))
code = code.replace('!SET_ARRAY_CASES_3D', set_array_cases(p3, 'pool3(h)%p'))
code = code.replace('!SET_STRING_CASES_2D', set_string_cases(p2, 'pool2(h)%p'))
code = code.replace('!SET_STRING_CASES_3D', set_string_cases(p3, 'pool3(h)%p'))
code = code.replace('!GET_CASES_2D', get_cases(out2, 'pool2(h)%p', 2))
code = code.replace('!GET_CASES_3D', get_cases(out3, 'pool3(h)%p', 3))
code = code.replace('!GET_SHAPE_CASES_2D', get_shape_cases(out2, 'pool2(h)%p', 2))
code = code.replace('!GET_SHAPE_CASES_3D', get_shape_cases(out3, 'pool3(h)%p', 3))

with open(os.path.join(root, 'src/module_geological_model_c.f90'), 'w') as f:
    f.write(code)

print(f'shim generated: {len(p2)} 2D params, {len(p3)} 3D params')
