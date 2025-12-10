
# Short Tutorial for RGM

This is short tutorial on how to install and use RGM.

## Table of Contents
- [Installation](#installation)
- [Parameters](#parameters)
- [Examples](#examples)


## Installation

Installation of `RGM` is straightforward:

```bash
git clone https://github.com/lanl/rgm.git
cd rgm
cd src
make
```

You need to install [`FLIT`](https://github.com/lanl/flit) before installing `RGM`, and set the path `flitdir = ...` in the [Makefile](https://github.com/lanl/rgm/tree/main/src/Makefile). 


## Parameters

Below is an incomplete list of important parameters for using `RGM` to generate random geological models. 

Here, I use `m` to represent a `rgm2_curved` or `rgm3_curved` class, i.e., 
```fortran
    type(rgm2_curved) :: m
```
or 
```fortran
    type(rgm3_curved) :: m
```

### Randomness

The random geological model generation process contains multiple steps and components (e.g., number of faults, layer shapes, fault geometry attributes, etc.), and each of these steps and components can be randomized through some parameters. 

> **`seed` (integer)** 
- **Description**: Seed for controlling ranodmness reproducibility. By defualt, `m%seed = -1`, meaning that every time the generated random model is different. For a value that is not `-1` (e.g., `m%seed = 123`), all randomness is reproducible and the resulting random model will be same over different realizations at different times. 
- **Default**: `-1`

### Dimension
> **`n1` (integer)** 
- **Description**: Number of grid points along axis-1 of the generated random model.
- **Default**: `128`

> **`n2` (integer)** 
- **Description**: Number of grid points along axis-2 of the generated random model.
- **Default**: `128`

> **`n3` (integer)** 
- **Description**: Number of grid points along axis-3 of the generated random model.
- **Default**: `128`

For instance, the user can set
```fortran
    m%n1 = 128
    m%n2 = 256
    m%n3 = 256
```
to generate random geological models with a size of `(n1, n2, n3) = (128, 256, 256)`. 

### Layer and reflectors

> **`nl` (integer)** 
- **Description**: Number of layers. Note that due to the limitation of the current algorithm, the resulting model may not have exactly `nl` layers, but should be very close in most cases. 
- **Default**: `20`

> **`refl_height` (real, dimension(1:2))** 
- **Description**: Height range of the reflectors at the bottom of the model, e.g., `m%refl_height =  = [0.0, 10.0]` sets the maximum height of reflectors at the bottom of the model to be `10`. Higher value of the second number results in larger deformation of layers. 
- **Default**: `[0.0, 10.0]`

> **`refl_height_top` (real, dimension(1:2))** 
- **Description**: Height range of the reflectors at the top of the model, e.g., `m%refl_height_top =  = [0.0, 5.0]` sets the maximum height of reflectors at the top of the model to be `5`. Higher value of the second number results in larger deformation of layers. 
- **Default**: `[0.0, 5.0]`

> **`lwv` (real)** 
- **Description**: Degree of variation of layer thickness in the vertical direction, i.e., how layer thickness for different layers varies. Valid range is `[0, 1]`. 
- **Default**: `0.25`

> **`lwh` (real)** 
- **Description**: Degree of variation of layer thickness in the horizontal directions, i.e., how layer thickness for the same layers varies in the horizontal directions. Valid range is `[0, 1]`. 
- **Default**: `0.1`

> **`refl_shape` (character(len=24))** 
- **Description**: Shape of the bottom reflectors. Valid options are `random`, `gaussian`, `cauchy`, `perlin`, and `custom`. 
- **Default**: `random`

> **`refl_shape_top` (character(len=24))**
- **Description**: Shape of the top reflectors. Valid options are `random`, `gaussian`, `cauchy`, `perlin`, and `custom`. In between the top and the bottom reflectors, the shapes of the reflectors will be interpolated depending on depth. 
- **Default**: `random`

To use `custom` reflector shapes, the user must provide the following arrays:

> **`refl` (real, allocatable, dimension(:, :))** 
- **Description**: A 2D array (or 1D array in the 2D case) storing the shape of the bottom reflector shape. Its size must be `(n2, n3)` for 3D and `n2` in 2D. 
- **Default**: `None`

> **`refl_top` (real, allocatable, dimension(:, :))** 
- **Description**: A 2D array (or 1D array in the 2D case) storing the shape of the top reflector shape. Its size must be `(n2, n3)` for 3D and `n2` in 2D. 
- **Default**: `None`

### Fault

> **`yn_fault` (logical)** 
- **Description**: Insert fault or not. When this option is `.false.`, no faults will be inserted. 
- **Default**: `.true.`

> **`nf` (integer)** 
- **Description**: Number of faults. 
- **Default**: `4`

> **`dip` (real, allocatable, dimension(:))** 
- **Description**: Range for defining the dip of the faults. This and the following three fault geometry parameters are defined in a similar way. Basically, one can specify `dip` with `2*nsf` values, where each pair of the values define a bounding range of the faults. For instance, `m%dip = [30.0, 40.0, 60.0, 120.0]` will define two sets of faults, where the first set of faults (`nf/2`) have a dip ranging from `[30, 40]` degrees, while the second set of faults (`nf/s`) have a dip ranging from `[60, 120]` degrees. 
- **Default**: `[70.0, 110.0]` (degrees)

> **`strike` (real, allocatable, dimension(:))** 
- **Description**: Range for defining the strikes of the faults. 
- **Default**: `[0.0, 180.0]` (degrees)

> **`rake` (real, allocatable, dimension(:))** 
- **Description**: Range for defining the rake of the faults. 
- **Default**: `[0.0, 180.0]` (degrees)

> **`disp` (real, allocatable, dimension(:))** 
- **Description**: Range for defining the displacement of the faults. 
- **Default**: `[5.0, 30.0]` (grid points)

> **`yn_regular_fault` (logical)** 
- **Description**: Generate two sets of faults rather than random-strike faults, each set has a domaint strike (i.e., quasi-parallel within each group). When this option is `.true.`, specify the fault geometry parameters using only two values like `m%dip = [50.0, 110.0]`, and the resulting dip will be `[0.95, 1.05]*m%dip(1)` for the first group, and `[0.95, 1.05]*m%dip(2)` for the second group. 
- **Default**: `.false.`

> **`fwidth` (real)** 
- **Description**: Fault thickness. 
- **Default**: `2.0` (grid points)

> **`delta_dip` (real, dimension(1:2))** 
- **Description**: Dip decrease of the faults at the bottom compared with the portion at the top (near-surface region). Larger values will result in less steep dips at the bottom. Assume a fault has a dip of `60` degrees at the top, then if `delta_dip` is not `[0, 0]`, the dip at the bottom could be `60 + r` where `r` is a random value drawn from the range of `delta_dip`. When `delta_dip` is negative, the dip increases with depth. 
- **Default**: `[15.0, 30.0]`

<!-- 
logical :: yn_group_faults = .false. -->

### Unconformity

> **`unconf` (integer)** 
- **Description**: Number of unconformity surfaces. 
- **Default**: `0`

> **`unconf_z` (real, dimension(1:2))** 
- **Description**: Relative depth range of the unconformity surfaces in the vertical direction. 
- **Default**: `[0.0, 0.5]`

> **`unconf_height` (real, dimension(1:2))** 
- **Description**: Real height of the unconformity surfaces in the vertical direction. 
- **Default**: `[5.0, 15.0]`

> **`unconf_smooth` (real)** 
- **Description**: Gaussian smooth sigma for applied to the unconformity surfaces. 
- **Default**: `0.0`

<!-- !==============================================================================================
integer :: unconf_nl = 99999
real, dimension(1:2) :: unconf_refl_height = [0.0, 5.0]
real :: unconf_refl_slope = -2.5
real :: unconf_refl_smooth = 10.0
character(len=12) :: unconf_refl_shape = 'random' -->

### Salt 

> **`yn_salt` (logical)** 
- **Description**: Insert salt body or not. 
- **Default**: `.false.`

> **`nsalt` (integer)** 
- **Description**: Number of salt bodies. 
- **Default**: `1` (only useful when `yn_salt = .true.`). 

> **`salt_radius` (real, dimension(1:2))** 
- **Description**: Salt radius relative to the size of the model. Larger values results in broader salt bodies. 
- **Default**: `[0.05*0.5*(n2 + n3), 0.2*0.5*(n2 + n3)]`

<!-- !==============================================================================================
real :: salt_radius_variation = 0.7
real :: salt_path_variation = 5.0
integer :: salt_nnode = 10
real, dimension(1:2) :: salt_top_z = [0.5, 0.8]
real :: salt_vp = 5000.0
real :: salt_rho = 2150.0
real, allocatable, dimension(:, :, :) :: salt
real :: salt_top_height = 20.0 -->
<!-- !> Salt body's Vs
real :: salt_vs = 4400.0
!> Is salt before or after unconformity?
logical :: salt_before_unconf = .true. -->

### Elastic model

> **`yn_elastic` (logical)** 
- **Description**: Generate elastic or just acoustic model. When this option is `.true.`, the generated models will contain not only `vp`, but also `vs`, `rho`, and the image arrays will contain `image_pp`, `image_ps`, `image_sp`, and `image_ss`. 
- **Default**: `.false.`

<!-- !==============================================================================================
!> Vp/Vs ratios
real, dimension(1:2) :: vpvsratio = [1.5, 1.8]
!> Elastic models
real, allocatable, dimension(:, :, :) :: vp, vs, rho
real :: rho_a = 310.0, rho_b = 0.25, rho_c = 0.0
!> Elastic images
real, allocatable, dimension(:, :, :) :: image_pp, image_ps, image_sp, image_ss-->


### RGT

> **`yn_rgt` (logical)** 
- **Description**: Generate RGT volume. 
- **Default**: `.false.`

<!-- 
integer :: nl = 20
real :: refl_smooth = 20.0
real :: refl_smooth_top = 40.0
real :: dt = 1.0e-3
real :: f0 = 150.0
real, dimension(1:2) :: refl_slope = [0.0, 0.0]
real, dimension(1:2) :: refl_slope_top = [0.0, 0.0]
real, dimension(1:3) :: noise_smooth = [1.0, 1.0, 1.0]
real :: noise_level = 0.05
character(len=24) :: wave = 'ricker'

integer :: ng = 2
!> Range of Gaussian standard devision along x2 for refl_shape = gaussian
real, dimension(1:2) :: refl_sigma2 = [0.0, 0.0]
!> Range of Gaussian mean along x2 for refl_shape = gaussian
real, dimension(1:2) :: refl_mu2 = [0.0, 0.0]
!> Range of Gaussian standard devision along x3 for refl_shape = gaussian
real, dimension(1:2) :: refl_sigma3 = [0.0, 0.0]
!> Range of Gaussian mean along x3 for refl_shape = gaussian
real, dimension(1:2) :: refl_mu3 = [0.0, 0.0]

!> Secondary reflector smoothing
real :: secondary_refl_smooth = 10.0
!> For Gaussian, Cauchy surface, whether to rotate
logical :: rotate_fold = .false.

logical :: yn_facies = .false.
real, allocatable, dimension(:, :, :) :: image, rgt, facies, fault
real, allocatable, dimension(:, :, :) :: fault_dip, fault_strike, fault_rake, fault_disp
real, dimension(1:3) :: psf_sigma = [5.0, 2.5, 2.5]
real, allocatable, dimension(:, :, :) :: psf
logical :: custom_psf = .false.
real :: facies_threshold = 0.0
character(len=12) :: noise_type = 'normal'
logical :: yn_conv_noise = .false.

real, allocatable, dimension(:) :: wave_filt_freqs, wave_filt_amps
!> Min value for scaling the facies
real :: vmin = 2000.0
!> Max value for scaling the facies; after scaling, the facies will fall in [vmin, vmax]
real :: vmax = 4000.0
!> Velocity perturbation of layers
real :: delta_v = 500.0

-->

## Examples

We place several examples in the [example](https://github.com/lanl/rgm/tree/main/example) directory. Some of the generated figures display in the front page. 