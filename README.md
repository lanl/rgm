# Description
**RGM: Random Geological Model Generation package**

![Example random geological models](github.com/lanl/rgm/example/rgm.png =512x)

We develop a Fortran package for generating realistic random geological models, including seismic images, faults, relative geological time, and salt bodies. 

The work was supported by Los Alamos National Laboratory (LANL) Laboratory Directory Research and Development (LDRD) project 20240322ER. LANL is operated by Triad National Security, LLC, for the National Nuclear Security Administration (NNSA) of the U.S. Department of Energy (DOE) under Contract No. 89233218CNA000001. The research used high-performance computing resources provided by LANL's Institutional Computing program. 

The work is under LANL open source approval reference O4778.

# Reference
A manuscript associated with the code is under preparation and will be released pending separate approval.  

# Requirement
`RGM` relies on [FLIT](github.com/lanl/flit). You need to install `FLIT` first. 

# Use
The code is written in Fortran. 

To install `RGM`, 

```
cd src
make
```

The compiled library file (with a name `librgm.a`) and module files will be at [lib](lib). 

We include several simple examples to use `RGM` in `example`. To try these examples,

```
cd example
make
```

The compiled executables will be at [example/bin](example/bin). Running these executables will generate images/faults/salt bodies in the directory. All the generated files will be in little-endian single-precision raw binary format, with dimensions specificed in the respective codes. 

The [Makefile](example/Makefile) in the [example](example) directory can serve as an example on how to use `RGM` in your code, including path inclusion and link of the compiled library. 

# License
&copy; 2024. Triad National Security, LLC. All rights reserved. 

This program is Open-Source under the BSD-3 License.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
- Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
- Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# Author
Kai Gao, <kaigao@lanl.gov>

We welcome feedbacks, bug reports, and improvement ideas on `RGM`. 

If you use this package in your research and find it useful, please cite it as

* Kai Gao, Ting Chen, 2024, RGM:  Random Geological Model Generation package, url: [github.com/lanl/rgm](github.com/lanl/rgm)
