#
# © 2024-2026. Triad National Security, LLC. All rights reserved.
#
# Python interface to RGM (random geological model generation).
#
# The heavy lifting is done by the Fortran library (lib/librgm.so); this
# package is a thin ctypes wrapper and does not reimplement any algorithm.
#
# Usage:
#
#     import rgm
#
#     p = rgm.rgm3(n1=128, n2=128, n3=128, nf=4, yn_rgt=True, seed=1234)
#     p.delta_strike = [15, 25]
#     p.yn_vary_disp = True
#     p.generate()
#
#     vp = p.vp          # numpy array, shape (n1, n2, n3)
#     image = p.image
#     fault = p.fault
#
# The parameter names are identical to the Fortran derived-type components
# of rgm2_curved/rgm3_curved (see doc/README.md). Output arrays: vp, vs,
# rho, image, image_pp/ps/sp/ss, rgt, facies, fault, fault_dip, fault_disp,
# (3D only: fault_strike, fault_rake), salt, karst, psf.
#

import ctypes
import os

import numpy


__all__ = ['rgm2', 'rgm3']


def _find_library():
    candidates = []
    if 'RGM_LIB' in os.environ:
        candidates.append(os.environ['RGM_LIB'])
    here = os.path.dirname(os.path.abspath(__file__))
    candidates.append(os.path.join(here, 'librgm.so'))
    candidates.append(os.path.abspath(os.path.join(here, '..', '..', 'lib', 'librgm.so')))
    for c in candidates:
        if os.path.isfile(c):
            return c
    raise OSError(
        'librgm.so not found. Build it with "make so" in rgm/src, and/or '
        'set the environment variable RGM_LIB to its full path. Searched: '
        + '; '.join(candidates))


_lib = ctypes.CDLL(_find_library())

for _pre in ('rgm2', 'rgm3'):
    getattr(_lib, _pre + '_create').restype = ctypes.c_int
    getattr(_lib, _pre + '_create').argtypes = []
    getattr(_lib, _pre + '_free').restype = None
    getattr(_lib, _pre + '_free').argtypes = [ctypes.c_int]
    getattr(_lib, _pre + '_generate').restype = None
    getattr(_lib, _pre + '_generate').argtypes = [ctypes.c_int]
    getattr(_lib, _pre + '_set_num').restype = ctypes.c_int
    getattr(_lib, _pre + '_set_num').argtypes = [
        ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_double]
    getattr(_lib, _pre + '_set_num_array').restype = ctypes.c_int
    getattr(_lib, _pre + '_set_num_array').argtypes = [
        ctypes.c_int, ctypes.c_char_p, ctypes.c_int,
        ctypes.POINTER(ctypes.c_double), ctypes.c_int]
    getattr(_lib, _pre + '_set_string').restype = ctypes.c_int
    getattr(_lib, _pre + '_set_string').argtypes = [
        ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_int]
    getattr(_lib, _pre + '_get_shape').restype = ctypes.c_int
    getattr(_lib, _pre + '_get_shape').argtypes = [
        ctypes.c_int, ctypes.c_char_p, ctypes.c_int,
        ctypes.POINTER(ctypes.c_int)]
    getattr(_lib, _pre + '_get_array').restype = ctypes.c_int
    getattr(_lib, _pre + '_get_array').argtypes = [
        ctypes.c_int, ctypes.c_char_p, ctypes.c_int,
        ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_float)]


class _RgmBase:

    _prefix = None
    _ndim = None

    def __init__(self, **kwargs):
        h = getattr(_lib, self._prefix + '_create')()
        if h < 0:
            raise RuntimeError('too many live rgm instances; free some first')
        object.__setattr__(self, '_h', h)
        for k, v in kwargs.items():
            setattr(self, k, v)

    def __setattr__(self, name, value):
        if name.startswith('_'):
            object.__setattr__(self, name, value)
            return
        c = getattr(_lib, self._prefix + '_set_num')
        ca = getattr(_lib, self._prefix + '_set_num_array')
        cs = getattr(_lib, self._prefix + '_set_string')
        bname = name.encode()
        if isinstance(value, str):
            ok = cs(self._h, bname, len(bname), value.encode(), len(value))
        elif isinstance(value, (bool, int, float, numpy.integer, numpy.floating)):
            ok = c(self._h, bname, len(bname), float(value))
        else:
            arr = numpy.asarray(value, dtype=numpy.float64).ravel()
            ok = ca(self._h, bname, len(bname),
                    arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), arr.size)
        if ok == 0:
            raise AttributeError(
                f'{type(self).__name__} has no settable parameter "{name}" '
                f'(or the value type does not match)')

    def __getattr__(self, name):
        # only called when normal lookup fails: treat as output array request
        if name.startswith('_'):
            raise AttributeError(name)
        gs = getattr(_lib, self._prefix + '_get_shape')
        ga = getattr(_lib, self._prefix + '_get_array')
        bname = name.encode()
        dims = (ctypes.c_int*3)()
        ok = gs(self._h, bname, len(bname), dims)
        if ok == 0:
            raise AttributeError(
                f'output array "{name}" is not available (not an output '
                f'name, or generate() has not been called with the '
                f'corresponding option enabled)')
        shape = [dims[i] for i in range(self._ndim)]
        buf = numpy.empty(int(numpy.prod(shape)), dtype=numpy.float32)
        ga(self._h, bname, len(bname), dims,
           buf.ctypes.data_as(ctypes.POINTER(ctypes.c_float)))
        return buf.reshape(shape, order='F')

    def generate(self):
        getattr(_lib, self._prefix + '_generate')(self._h)

    def free(self):
        if getattr(self, '_h', -1) >= 0:
            getattr(_lib, self._prefix + '_free')(self._h)
            object.__setattr__(self, '_h', -1)

    def __del__(self):
        try:
            self.free()
        except Exception:
            pass


class rgm2(_RgmBase):
    """2D random geological model (Fortran type rgm2_curved)."""

    _prefix = 'rgm2'
    _ndim = 2


class rgm3(_RgmBase):
    """3D random geological model (Fortran type rgm3_curved)."""

    _prefix = 'rgm3'
    _ndim = 3
