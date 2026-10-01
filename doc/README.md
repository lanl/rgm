
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
- **Description**: Shape of the bottom reflectors. Valid options are `random`, `gaussian`, `cauchy`, `perlin`, `fold`, and `custom`. 
- **Default**: `random`

> **`refl_shape_top` (character(len=24))**
- **Description**: Shape of the top reflectors. Valid options are `random`, `gaussian`, `cauchy`, `perlin`, `fold`, `custom`, and `same`. With `same`, the top reflectors have the same shape as the bottom reflectors, scaled to the height range `refl_height_top`. In between the top and the bottom reflectors, the shapes of the reflectors will be interpolated depending on depth. 
- **Default**: `random`

To use `custom` reflector shapes, the user must provide the following arrays:

> **`refl` (real, allocatable, dimension(:, :))** 
- **Description**: A 2D array (or 1D array in the 2D case) storing the shape of the bottom reflector shape. Its size must be `(n2, n3)` for 3D and `n2` in 2D. 
- **Default**: `None`

> **`refl_top` (real, allocatable, dimension(:, :))** 
- **Description**: A 2D array (or 1D array in the 2D case) storing the shape of the top reflector shape. Its size must be `(n2, n3)` for 3D and `n2` in 2D. 
- **Default**: `None`

For `refl_shape = 'fold'` (or `refl_shape_top = 'fold'`), the reflector is a fold train $r(x) = A(x)\left[-\cos\psi + c\cos 2\psi\right]$, where $\psi$ is the fold phase measured from a trough, $c$ is the crest asymmetry, and $A$ is the local amplitude. The wavelength and the amplitude drift from fold to fold, so the train is not periodic. In 3D, the fold train runs across fold axes with a random strike; along strike, the anticlines rise and plunge out, and the fold axes bend.

> **`refl_fold_lambda` (real, dimension(1:2))**
- **Description**: Range of the mean fold wavelength in grid points. When `[0, 0]`, it is set to `[0.2, 0.5]*n2` in 2D and `[0.2, 0.5]*min(n2, n3)` in 3D.
- **Default**: `[0.0, 0.0]`

> **`refl_fold_crest` (real, dimension(1:2))**
- **Description**: Range of the crest asymmetry, clipped to `[-0.25, 0.25]`. Positive values give sharp anticlines and broad synclines; negative values give broad (box) anticlines and sharp synclines.
- **Default**: `[-0.25, 0.25]`

> **`refl_fold_vergence` (real, dimension(1:2))**
- **Description**: Range of the magnitude of the limb asymmetry (vergence), clipped to `[0, 0.9]`. The sense of vergence is random per model and common to the bottom and top reflectors.
- **Default**: `[0.0, 0.5]`

> **`refl_fold_lambda_drift` (real)**
- **Description**: Fold-to-fold drift of the wavelength; the local wavelength varies by a factor between `exp(-d)` and `exp(d)`.
- **Default**: `0.4`

> **`refl_fold_amp_drift` (real)**
- **Description**: Fold-to-fold drift of the amplitude; the local amplitude varies by a factor between `1 - d` and `1 + d`.
- **Default**: `0.5`

> **`refl_fold_mode` (character(len=12))**
- **Description**: Form of the vergence. With `smooth`, the fold phase is warped smoothly, and a vergence above about 0.7 develops a shoulder on the gentle limb. With `limb`, each limb is uniformly steep or gentle, and a vergence above about 0.5 makes the steep limb approach a step.
- **Default**: `smooth`

> **`refl_fold_strike` (real, dimension(1:2)); only for `rgm3_curved`**
- **Description**: Range of the strike of the fold axes in degrees, measured from x2 toward x3; common to the bottom and top reflectors.
- **Default**: `[0.0, 180.0]`

> **`refl_fold_plunge` (real); only for `rgm3_curved`**
- **Description**: Relative variation of the fold uplift along strike. Values near or above 1 make anticlines rise and plunge out along strike (periclines); `0` gives cylindrical folds.
- **Default**: `0.6`

> **`refl_fold_plunge_length` (real); only for `rgm3_curved`**
- **Description**: Along-strike correlation length of the uplift variation, in fold wavelengths.
- **Default**: `2.5`

> **`refl_fold_wobble` (real); only for `rgm3_curved`**
- **Description**: Lateral shift of the fold axes along strike, in fold wavelengths, which bends the fold axes.
- **Default**: `0.15`

For `refl_shape = 'gaussian'` or `'cauchy'`, the bumps can be asymmetric:

> **`refl_skew` (real, dimension(1:2))**
- **Description**: Range of the skew magnitude `s`, drawn per bump and clipped to `[0, 0.95]`. A bump has the scale `sigma*(1 - s)` on one side and `sigma*(1 + s)` on the other, which gives a steep limb and a gentle limb, as in forced folds and fault-propagation folds. In 3D, the skew acts along the (rotated) x2 axis of the bump. When `[0, 0]`, the bumps are symmetric.
- **Default**: `[0.0, 0.0]`

> **`refl_skew_common` (logical)**
- **Description**: Whether all bumps of a model are steep on the same side (common vergence, as in a thrust belt).
- **Default**: `.false.`

> **`refl_background` (real)**
- **Description**: Amplitude of a smooth random background added to the bumps, relative to the maximum bump height, so the reflectors between the bumps are not flat.
- **Default**: `0.0`

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

> **`delta_strike` (real, dimension(1:2)); only for `rgm3_curved`** 
- **Description**: Range of the maximum strike deviation along a fault. When not `[0, 0]`, the strike of each fault spatially varies along the fault, so in map view a fault is no longer a straight line but a smooth curve. Each fault draws a maximum deviation `r` from the range of `delta_strike`, and its local strike varies smoothly within `strike ± r` along the fault. The along-fault strike variation is deliberately low-wavenumber (smooth, gentle bends), as high-wavenumber oscillation of fault strike is geologically uncommon. In this case, `fault_strike` stores the spatially varying local strike angles, and `fault_dip` stores the local true dip angles that account for the strike deviation, i.e., `tan(dip_local) = tan(dip)/cos(strike_deviation)`. 
- **Default**: `[0.0, 0.0]` (degrees; constant-strike faults)
- **Note**: For strike-varying faults, the slip direction of the block shifting follows the local strike of the fault. 
- **See also**: [varying_strike_and_displacement.md](varying_strike_and_displacement.md) for the underlying formulation.

> **`strike_nperiod` (integer); only for `rgm3_curved`** 
- **Description**: Number of undulation periods of the along-fault strike deviation when `delta_strike` is not `[0, 0]`. A larger value results in more oscillatory faults in map view. 
- **Default**: `2`

> **`yn_vary_disp` (logical); for `rgm2_curved` and `rgm3_curved`** 
- **Description**: Whether to spatially vary the displacement on each fault surface. When `.true.`, the displacement of each fault follows an elliptical slip patch defined in the (along-strike, depth) coordinates of the fault (in 2D, the patch is an interval along the fault): the displacement is maximum at the patch center and diminishes to zero at the patch boundary (the fault tip line), so a fault can die out within the model rather than always cutting through the entire volume, and the hanging wall deforms gently near the fault tips. In this case, `fault_disp` stores the spatially varying displacement magnitudes on the fault surfaces, and the fault is only marked where the local displacement is `>= 0.5` grid points (below that, the offset is not visible). 
- **Default**: `.false.` (displacement is constant on each fault surface)
- **See also**: [varying_strike_and_displacement.md](varying_strike_and_displacement.md) for the underlying formulation.

> **`disp_radius_strike` (real, dimension(1:2)); only for `rgm3_curved`** 
- **Description**: Range of the along-strike semi-axis of the elliptical slip patch when `yn_vary_disp = .true.`, relative to the average lateral model size `0.5*(n2 + n3)`. Smaller values result in faults that die out laterally within the model. 
- **Default**: `[0.6, 1.2]`

> **`disp_radius_dip` (real, dimension(1:2)); for `rgm2_curved` and `rgm3_curved`** 
- **Description**: Range of the along-depth semi-axis of the elliptical slip patch when `yn_vary_disp = .true.`, relative to the model depth `n1`. Smaller values result in faults that die out vertically within the model. 
- **Default**: `[0.6, 1.2]`

> **`disp_center_dip` (real, dimension(1:2)); for `rgm2_curved` and `rgm3_curved`** 
- **Description**: Range of the center depth of the slip patch when `yn_vary_disp = .true.`, relative to the model depth `n1`. All the displacement tapering (and therefore the deformation of layers) occurs between the patch center and the patch boundary, so a shallower center together with a smaller `disp_radius_dip` keeps the fault tips and the associated deformation away from the deep, high-velocity section of the model. 
- **Default**: `[0.25, 0.75]`

> **`disp_center_strike` (real, dimension(1:2)); only for `rgm3_curved`** 
- **Description**: Range of the along-strike center of the slip patch when `yn_vary_disp = .true.`, relative to the average lateral model size `0.5*(n2 + n3)` and to the fault center. 
- **Default**: `[-0.2, 0.2]`

> **`yn_disp_decay` (logical); for `rgm2_curved` and `rgm3_curved`** 
- **Description**: Whether to decay the fault displacement away from the fault along the fault normal. When `.false.` (default), the displaced block shifts rigidly, which is the correct kinematics for through-going faults (a marker bed stays offset at any distance from the fault). When `.true.`, the displacement of the shifted block decays with a Gaussian profile of the distance to the fault, mimicking the deformation halo of a finite slip patch: fault drag, rollover (reverse drag), and folding above blind fault tips. The decay is applied only within the displaced block, and the displacement at the fault surface itself is unaffected, so the fault throw at the fault and the `fault_disp` label are unchanged. Can be combined freely with `yn_vary_disp` and `delta_strike`. 
- **Default**: `.false.`

> **`disp_decay_width` (real, dimension(1:2)); for `rgm2_curved` and `rgm3_curved`** 
- **Description**: Range of the displacement decay width (the standard deviation of the Gaussian decay profile) when `yn_disp_decay = .true.`, relative to the model depth `n1`. Smaller values concentrate the deformation near the fault (drag-fold style), while larger values approach the rigid-block end member. 
- **Default**: `[0.5, 1.0]`

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

> **`unconf_shape` (character(len=24)); RGM v2.0** 
- **Description**: Shape of the unconformity surfaces. Can be one of `random` (Perlin-noise topography, the classic behavior), `meander_channel` (meandering river channels), `meander_canyon` (meandering incised canyon), `drainage_channel` (dendritic drainage channel network), or `drainage_canyon` (canyon carved by a dendritic drainage network). For the channel/canyon shapes, the unconformity topography is the erosional depth map of a geomorphological simulation (map-view simulation in 3D; a random cross-flow section of the map-view simulation in 2D), and `unconf_height` still bounds the total erosional relief. 
- **Default**: `random`

> **`unconf_channel_width` (real, dimension(1:2)); RGM v2.0** 
- **Description**: Channel/canyon width as an (approximate) fraction of the lateral model extent, drawn per unconformity surface. Channels are rendered relatively thin and densely distributed; canyons are major and few. 
- **Default**: `[0.03, 0.08]`

> **`unconf_channel_sinuosity` (real); RGM v2.0** 
- **Description**: Meander maturity for `meander_*` shapes; scales the number of channel migration iterations. Larger values produce more strongly meandering channels. 
- **Default**: `1.0`

> **`unconf_channel_length` (real); RGM v2.0** 
- **Description**: Total centerline length of the meandering channel/canyon for the `meander_*` shapes, in the internal units of the migration simulation. Because the simulated map is fit to the model grid, this effectively controls how many meander bends span the model laterally: the default gives roughly 10 well-developed gooseneck bends, and doubling the length roughly doubles the bend count. The migration width is internally capped at its calibrated value, so longer channels develop more bends of the same shape rather than degenerating; values below the default give fewer, larger bends (clamped from below at `10000` to keep enough centerline nodes for realistic meandering to develop). Note that the loop amplitude of a meandering river is physically capped by neck cutoffs, so for values beyond roughly `2&ndash;3x` the default, the channel gradually reads as a wiggly band rather than distinct goosenecks &mdash; the same look real rivers have when mapped at regional scale. The practically useful range is about `10000&ndash;60000`. 
- **Default**: `25000.0`

> **`unconf_channel_density` (real, dimension(1:2)); RGM v2.0** 
- **Description**: Drainage density (fraction of cells covered by channels) for the `drainage_*` shapes, drawn per unconformity surface. The `drainage_canyon` shape internally reduces the density so that canyons remain major and few. 
- **Default**: `[0.02, 0.08]`

> **`unconf_topo` (real); RGM v2.0** 
- **Description**: Amplitude of the low-wavenumber background (interfluve) topography relative to the channel/canyon incision depth; `0` = flat background between channels. 
- **Default**: `0.25`

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

### Karst (RGM v2.0)

Karst cave systems are tube-network cave geometries inserted into the medium
parameter models like salt bodies: the caves overwrite `vp`/`vs`/`rho` with
`karst_vp`/`karst_vs`/`karst_rho`, the mask is returned in `%karst`, and the
labels (`fault*`, `rgt`, `facies`) are zeroed inside the caves. Karst can be
freely combined with faults, salt (karst carves salt where they overlap),
unconformities (`karst_before_unconf`, analogous to `salt_before_unconf`),
and elastic models.

> **`yn_karst` (logical)** 
- **Description**: Insert a karst cave system or not. 
- **Default**: `.false.`

> **`karst_z` (real, dimension(1:2))** 
- **Description**: Depth window (relative to `n1`) containing the karst passages. 
- **Default**: `[0.4, 0.9]`

> **`karst_npassage` (integer)** 
- **Description**: Number of cave passages. 
- **Default**: `10`

> **`karst_nctrl` (integer)** 
- **Description**: Number of control points per passage; more points produce longer and more complex passages. 
- **Default**: `40`

> **`karst_radius` (real, dimension(1:2))** 
- **Description**: Range of the mean passage radius in grid points; `[0, 0]` = automatically set to `[0.015, 0.03]*n1`. 
- **Default**: `[0.0, 0.0]`

> **`karst_radius_variation` (real)** 
- **Description**: Fractional standard deviation of the radius along a passage. 
- **Default**: `0.35`

> **`karst_tortuosity` (real)** 
- **Description**: Sinuosity of the passages: 0 = straight, 1 = fully random walk. 
- **Default**: `0.6`

> **`karst_connect` (real)** 
- **Description**: Probability that a new passage branches off an existing one. 
- **Default**: `0.7`

> **`karst_vp`, `karst_vs`, `karst_rho` (real)** 
- **Description**: Medium properties filling the caves. The defaults represent water/sediment-filled voids (low-velocity anomalies producing the classic "string-of-beads" seismic response). 
- **Default**: `2500.0`, `1200.0`, `2000.0`

> **`karst_before_unconf` (logical)** 
- **Description**: Whether the karst formed before the unconformities (and is therefore eroded above them). 
- **Default**: `.true.`

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

### Image

By default, the image is the vertical reflectivity convolved with a separable point spread function (PSF): a vertical wavelet multiplied by lateral Gaussians (`psf_sigma`). With this PSF, steeply dipping reflectors fade; with the default parameters, a reflector keeps about 38% of its amplitude at a dip of 30 degrees and 6% at 60 degrees, so steep fold limbs can look like faults. The following options make the image dip-independent and give it the character of migrated images. All of them are off by default, and the default images are unchanged. With or without them, the PSF is centered on the sample that the convolution aligns with the input, so the image is aligned with the model and the labels for both even and odd model sizes.

> **`yn_dip_independent` (logical)**
- **Description**: Dip-independent imaging. The reflectivity is normalized by the cosine of the layer dip, estimated from the RGT with faults, salt and karst excluded, and convolved with an isotropic PSF whose projection onto any direction is the zero-phase source wavelet. Every reflector then keeps its wavelet and amplitude whatever its dip, and each fault shows as reflector offsets plus a weak fault-plane reflection from the juxtaposed layers. The RGT and the fault labels are computed internally even when `yn_rgt` or `yn_fault` is `.false.`, and released afterwards. For the isotropic PSF, `psf_sigma(2:3)` are not used.
- **Default**: `.false.`

> **`f0_bottom` (real)**
- **Description**: Center frequency of the source wavelet at the bottom of the model. When `> 0`, RGM computes two images, one with a wavelet of center frequency `f0` and one with a wavelet of center frequency `f0_bottom`. It blends them with weights that change linearly with depth, so the top of the image uses only the `f0` wavelet and the bottom uses only the `f0_bottom` wavelet. This mimics the loss of high frequencies with depth, and the amplitude of a flat reflector stays constant.
- **Default**: `0.0`

> **`illum_level` (real)**
- **Description**: Standard deviation of the logarithm of a smooth multiplicative gain applied to the image, so the reflectors fade and strengthen along their length, as with uneven illumination. `0` means no gain.
- **Default**: `0.0`

> **`illum_smooth` (real)**
- **Description**: Correlation length of the illumination gain in grid points.
- **Default**: `30.0`

> **`noise_type` (character(len=12))**
- **Description**: Type of noise: `normal`, `uniform`, `exp`, `wavenumber`, or `migration`. The `migration` noise has the character of migrated images. It is a weighted sum of three components, each band-limited by the PSF: worm noise (short reflector-parallel segments in patches, which follow the layers through the RGT), swing noise (steeply dipping streaks that grow with depth and, by default, cross in both directions; see `noise_swing_direction`), and background noise. It is added after the PSF convolution regardless of `yn_conv_noise`, and `noise_level` is then the ratio of the noise RMS to the image RMS.
- **Default**: `normal`

> **`noise_mix` (real, dimension(1:3))**
- **Description**: Relative weights of the worm, swing and background components of the `migration` noise.
- **Default**: `[0.35, 0.30, 0.15]`

> **`noise_worm_length` (real)**
- **Description**: Along-layer correlation length of the worm noise in grid points.
- **Default**: `10.0`

> **`noise_swing_dip` (real, dimension(1:2))**
- **Description**: Range of the dips of the swing noise in degrees.
- **Default**: `[25.0, 75.0]`

> **`noise_swing_nfam` (integer)**
- **Description**: Number of dip families of the swing noise, each with a random dip from `noise_swing_dip`. With `noise_swing_direction = 'dual'`, each family contains both senses of dip, mirrored about the vertical (in 3D, a random azimuth and its opposite), with random relative strengths, and `1` gives a single crossing pattern.
- **Default**: `2`

> **`noise_swing_direction` (character(len=12))**
- **Description**: Direction of the swing noise. With `dual`, each dip family contains both senses of dip, so the streaks cross in both directions. With `single`, all dip families dip in the same, randomly chosen sense (in 3D, they share one random azimuth), so the streaks dip in one direction.
- **Default**: `dual`

> **`jitter_shift` (real, dimension(1:2))**
- **Description**: Trace-to-trace jitter: the maximum smooth static shift and the standard deviation of the random static shift of the traces, in grid points. It roughens the reflectors slightly, as residual statics do.
- **Default**: `[0.0, 0.0]`

> **`jitter_gain` (real)**
- **Description**: Standard deviation of a random gain applied to each trace, which gives faint vertical striping.
- **Default**: `0.0`

For example, the following settings give images that resemble migrated field data:
```fortran
    m%yn_dip_independent = .true.
    m%noise_type = 'migration'
    m%noise_level = 0.5
    m%f0_bottom = 100.0
    m%illum_level = 0.35
    m%jitter_shift = [0.8, 0.25]
    m%jitter_gain = 0.05
```

## Examples

We place several examples in the [example](https://github.com/lanl/rgm/tree/main/example) directory. Some of the generated figures display in the front page. 