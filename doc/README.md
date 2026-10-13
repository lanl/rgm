# Parameters for RGM

This document lists the parameters of the RGM model generators `rgm2_curved` (2D) and `rgm3_curved` (3D), with their types, default values, and meanings.

`m` is an `rgm2_curved` or `rgm3_curved` object, e.g., `type(rgm3_curved) :: m`. The Python classes `rgm.rgm2` and `rgm.rgm3` take the same parameter names, except the custom arrays `refl`, `refl_top`, and `psf`. Most two-value parameters `[a, b]` are ranges, from which each model draws a random value. Lengths are in grid points and depths are relative to `n1`, unless noted. (2D) or (3D) marks a parameter of only one type.

## Randomness

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `seed` | integer | `-1` | Random seed. `-1` gives a different model at each run; any other value makes the model reproducible. |

## Dimension

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `n1` | integer | `128` | Grid points along x1 (depth). |
| `n2` | integer | `128` | Grid points along x2. |
| `n3` (3D) | integer | `128` | Grid points along x3. |

## Layers and reflectors

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `nl` | integer | `20` | Number of layers (approximate). |
| `lwv` | real | `0.25` | Thickness variation between layers: each layer is `[1 - lwv, 1 + lwv]` times the mean thickness. |
| `lwh` | real | `0.1` | Lateral thickness variation along each layer, relative to its thickness. |
| `secondary_refl_smooth` | real | `10.0` | Smoothing length of the lateral thickness variation. |
| `vmin`, `vmax` | real | `2000.0`, `4000.0` | Range of the layer Vp (m/s). |
| `delta_v` | real | `500.0` | Random Vp perturbation of each layer (m/s), before the rescaling to `[vmin, vmax]`. |
| `rho_a`, `rho_b`, `rho_c` | real | `310.0`, `0.25`, `0.0` | Density from Gardner's relation, `rho = rho_a*vp**rho_b + rho_c`. |
| `refl_shape` | character | `'random'` | Shape of the bottom reflector: `random`, `gaussian`, `cauchy`, `perlin`, `fold`, or `custom`. |
| `refl_shape_top` | character | `'random'` | Shape of the top reflector: the same options, or `same` (the bottom shape). The layers in between blend the two. |
| `refl_height` | real(2) | 2D: `[0.0, 20.0]`<br>3D: `[0.0, 10.0]` | Relief of the bottom reflector within the model, `b - a`; `a = b` gives a flat reflector. For `gaussian` and `cauchy`, the bump heights are also drawn from `[a, b]`. |
| `refl_height_top` | real(2) | 2D: `[0.0, 10.0]`<br>3D: `[0.0, 5.0]` | Relief of the top reflector, as for `refl_height`. |
| `refl_slope` | 2D: real<br>3D: real(2) | 2D: `0.0`<br>3D: `[0.0, 0.0]` | Linear rise of the bottom reflector across the model (3D: along x2 and x3). |
| `refl_slope_top` | 2D: real<br>3D: real(2) | 2D: `0.0`<br>3D: `[0.0, 0.0]` | Linear rise of the top reflector across the model. |
| `refl_smooth` | real | `20.0` | Smoothing length of the `random` and `perlin` bottom reflectors. |
| `refl_smooth_top` | real | 2D: `20.0`<br>3D: `40.0` | Smoothing length of the `random` and `perlin` top reflectors. |
| `refl`, `refl_top` | 2D: real(:)<br>3D: real(:, :) | | Bottom and top reflector shapes for `custom`, of size `n2` (2D) or `(n2, n3)` (3D). Larger values are shallower, and `refl_height` and `refl_height_top` set the relief. |

Gaussian and Cauchy bumps (`refl_shape = 'gaussian'` or `'cauchy'`):

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `ng` | integer | `2` | Number of bumps. |
| `refl_mu2`, `refl_mu3` (3D) | real(2) | `[0.0, 0.0]` | Range of the bump centers along x2 and x3. `[0, 0]`: across the model. |
| `refl_sigma2`, `refl_sigma3` (3D) | real(2) | `[0.0, 0.0]` | Range of the bump widths along x2 and x3. `[0, 0]`: `[0.05, 0.15]` of the model size, including the margin for fault displacement. |
| `rotate_fold` (3D) | logical | `.false.` | Rotate the bumps by random angles in map view. |
| `refl_skew` | real(2) | `[0.0, 0.0]` | Range of the bump skew, clipped to `[0, 0.95]`; a skewed bump has a steep limb and a gentle limb. |
| `refl_skew_common` | logical | `.false.` | Give all bumps the same steep side; with `rotate_fold`, along each bump's own axis. |
| `refl_background` | real | `0.0` | Amplitude of a smooth random background added to the bumps, relative to the bump height. |

Fold trains (`refl_shape = 'fold'`):

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `refl_fold_lambda` | real(2) | `[0.0, 0.0]` | Range of the mean fold wavelength. `[0, 0]`: `[0.2, 0.5]*n2` (3D: `[0.2, 0.5]*min(n2, n3)`). |
| `refl_fold_crest` | real(2) | `[-0.25, 0.25]` | Range of the crest asymmetry, in `[-0.25, 0.25]`. Positive: sharp anticlines and broad synclines; negative: the reverse. |
| `refl_fold_vergence` | real(2) | `[0.0, 0.5]` | Range of the limb asymmetry (vergence), in `[0, 0.9]`. |
| `refl_fold_lambda_drift` | real | `0.4` | Fold-to-fold wavelength drift `d`: factors from `exp(-d)` to `exp(d)`. |
| `refl_fold_amp_drift` | real | `0.5` | Fold-to-fold amplitude drift `d`: factors from `1 - d` to `1 + d`. |
| `refl_fold_mode` | character | `'smooth'` | Form of the vergence: `smooth` (smooth phase warp) or `limb` (each limb uniformly steep or gentle). |
| `refl_fold_strike` (3D) | real(2) | `[0.0, 180.0]` | Range of the strike of the fold axes (degrees, from x2 toward x3). |
| `refl_fold_plunge` (3D) | real | `0.6` | Variation of the fold uplift along strike. Near or above `1`: periclines; `0`: no variation. |
| `refl_fold_plunge_length` (3D) | real | `2.5` | Along-strike correlation length of the uplift variation and of the axis wobble, in fold wavelengths. |
| `refl_fold_wobble` (3D) | real | `0.15` | Lateral shift of the fold axes along strike, in fold wavelengths. |

## Faults

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `nf` | integer | `4` | Number of faults. |
| `yn_fault` | logical | `.true.` | Fill the fault labels. `fault` holds `1` to `nf` for the faults, `nf + 1` to `nf + nthrust` for the thrust faults, and `0` elsewhere. It does not control the insertion of faults, and is set to `.false.` when `nf = 0` and `nthrust = 0`. |
| `fwidth` | real | `2.0` | Width of the fault labels. |
| `dip` | real(:) | `[70.0, 110.0]` | Ranges of the fault dip (degrees), one pair per fault set, e.g., `[30.0, 40.0, 60.0, 120.0]` for two sets. |
| `strike` (3D) | real(:) | `[0.0, 180.0]` | Ranges of the fault strike (degrees), one pair per fault set. |
| `rake` (3D) | real(:) | `[0.0, 180.0]` | Ranges of the fault rake (degrees), one pair per fault set. `90` is dip slip; `0` and `180` are strike slip, in opposite directions. |
| `disp` | real(:) | `[5.0, 30.0]` | Ranges of the fault displacement, one pair per fault set. Positive moves the hanging wall down (normal fault), and negative moves it up (reverse fault); in 3D, this holds for a rake between `0` and `180`. Unless `yn_regular_fault` is `.true.`, faults 2, 4, 6, ... take the opposite sign, so the faults alternate between normal and reverse. |
| `delta_dip` | real(2) | `[15.0, 30.0]` | Range of the dip decrease from the top to the bottom of a fault (degrees), which gives listric faults. |
| `yn_regular_fault` | logical | `.false.` | Two sets of quasi-regularly spaced faults (needs `nf >= 2`). `dip` and `disp` (and in 3D `strike` and `rake`) then give one value per set. |
| `yn_group_faults` | logical | `.false.` | With `yn_regular_fault`, keep the faults of each set next to each other. |
| `delta_strike` (3D) | real(2) | `[0.0, 0.0]` | Range of the maximum strike deviation along a fault (degrees), which curves the faults in map view; the slip follows the local strike. |
| `strike_nperiod` (3D) | integer | `2` | Number of undulation periods of the strike along a fault, including the margin outside the model. |
| `yn_vary_disp` | logical | `.false.` | Vary the displacement over an elliptical slip patch, so that faults can die out within the model. |
| `disp_radius_dip` | real(2) | `[0.6, 1.2]` | Range of the semi-axis of the slip patch along depth, relative to `n1`. |
| `disp_radius_strike` (3D) | real(2) | `[0.6, 1.2]` | Range of the semi-axis of the slip patch along strike, relative to `0.5*(n2 + n3)`. |
| `disp_center_dip` | real(2) | `[0.25, 0.75]` | Range of the depth of the slip-patch center. |
| `disp_center_strike` (3D) | real(2) | `[-0.2, 0.2]` | Range of the along-strike offset of the slip-patch center from the fault center, relative to `0.5*(n2 + n3)`. |
| `yn_disp_decay` | logical | `.false.` | Decay the displacement away from the fault (drag and rollover). |
| `disp_decay_width` | real(2) | `[0.5, 1.0]` | Range of the decay width (Gaussian standard deviation), relative to `n1`. |

The hanging wall is the block above the fault, and it is the only block that moves. In 2D, a dip below 90 degrees makes the fault dip toward -x2, and a dip above 90 degrees toward +x2; in 3D, the two cases dip toward opposite sides of the strike. The sign of the displacement sets the fault type in 2D and 3D, and `fault_disp` keeps that sign:

| Displacement | Hanging wall | Fault | `fault_disp` |
| --- | --- | --- | --- |
| `disp > 0` | Moves down | Normal | Positive |
| `disp < 0` | Moves up | Reverse | Negative |

Thrust faults (`nthrust`) are reverse faults, so their `fault_disp` is negative. In 3D, the rake splits the slip into dip slip and strike slip: a rake of 90 degrees gives pure dip slip, and 0 or 180 degrees gives pure strike slip, which moves the hanging wall only sideways. The sign rule holds for a rake between 0 and 180 degrees.

Unless `yn_regular_fault` is `.true.`, faults 2, 4, 6, ... take the opposite sign of their drawn displacement. A positive range therefore gives both fault types: `disp = [5.0, 10.0]` makes faults 1, 3, 5, ... normal and faults 2, 4, 6, ... reverse. A single fault (`nf = 1`) keeps the sign of the range. With `yn_regular_fault = .true.`, the faults do not alternate, and each fault set takes the sign of its value; for example, `disp = [5.0, -5.0]` gives one set of normal faults and one set of reverse faults.

A range that includes zero, such as `disp = [-10.0, 10.0]`, adds no fault type but can give faults with almost no offset. Without `yn_vary_disp`, the labels still mark such a fault along its whole trace, although the model shows no offset.

## Thrust faults

Thrust faults are older than the other faults. In a model with unconformities, only the unit below the deepest unconformity has thrusts. All thrusts of a model share a lower flat (decollement) and, optionally, an upper flat (roof); the flats follow horizons, except that a decollement deep below the model is horizontal. The fault labels mark a thrust where its stratigraphic offset is at least 0.5 grid points, and `fault_disp` stores the slip as a negative value, as for a reverse fault; in 3D, `fault_rake` is 90 degrees. The RGT decreases across a thrust, where older layers lie above younger layers. With `nthrust = 0`, no thrust is inserted. Typical settings:

| Structure | Parameters |
| --- | --- |
| Fault-propagation fold | `nthrust = 1`, `thrust_tip_z = [0.3, 0.5]` |
| Fault-bend fold | `nthrust = 1`, `thrust_roof = 1.0` |
| Imbricate fan | `nthrust = 3` |
| Duplex | `nthrust = 4`, `thrust_roof = 1.0`, `thrust_roof_z = [0.5, 0.55]`, `thrust_z = [0.85, 0.9]` |
| Pop-up | `nthrust = 2`, `thrust_back = 1.0` |

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `nthrust` | integer | `0` | Number of thrust faults. |
| `thrust_dip` | real(2) | `[20.0, 35.0]` | Range of the ramp dip (degrees). |
| `thrust_disp` | real(2) | `[0.0, 0.0]` | Range of the slip. `[0, 0]`: `[0.05, 0.15]` of the model size along the transport direction. |
| `thrust_z` | real(2) | `[0.85, 1.1]` | Range of the decollement depth; above `1`, the decollement is below the model. |
| `thrust_tip_z` | real(2) | `[-0.2, 0.5]` | Range of the final tip depth of a thrust without a roof; below `0`, the tip is above the model. |
| `thrust_bend` | real(2) | `[0.0, 0.0]` | Range of the length of the rounded fault bends. `[0, 0]`: `[0.15, 0.3]*n1`. |
| `thrust_widen` | real(2) | `[0.2, 0.4]` | Range of the growth of the fold-hinge width with the distance from the fault. |
| `thrust_ps` | real(2) | `[1.5, 3.0]` | Range of the ratio of tip propagation to slip (P/S). |
| `thrust_trishear` | real(2) | `[40.0, 100.0]` | Range of the apical angle of the trishear zone (degrees). |
| `thrust_roof` | real | `0.0` | Probability of a roof shared by the thrusts, which gives fault-bend folds and duplexes. |
| `thrust_roof_z` | real(2) | `[0.3, 0.6]` | Range of the roof depth. |
| `thrust_back` | real | `0.0` | Probability that a thrust after the first one is a backthrust, which gives a pop-up. |
| `thrust_sequence` | character | `'forward'` | Order of thrusting: `forward` (piggyback), `backward`, or `random`. |
| `thrust_vergence` (2D) | character | `'random'` | Transport direction: `right` (toward +x2), `left`, or `random`. |
| `thrust_azimuth` (3D) | real(2) | `[0.0, 360.0]` | Range of the transport azimuth (degrees, from x2 toward x3). |
| `thrust_delta_strike` (3D) | real(2) | `[0.0, 0.0]` | Range of the maximum strike deviation along the ramp traces (degrees). |
| `thrust_radius` (3D) | real(2) | `[0.0, 0.0]` | Range of the half length of the slip patch along strike, relative to `0.5*(n2 + n3)`. `[0, 0]`: constant slip. |

## Unconformity

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `unconf` | integer | `0` | Number of unconformities. |
| `unconf_z` | real(2) | `[0.0, 0.5]` | Range of the unconformity depth. |
| `unconf_height` | real(2) | `[5.0, 15.0]` | Range of the relief of the unconformity surfaces. |
| `unconf_smooth` | real | `0.0` | Smoothing length of the unconformity surfaces. |
| `unconf_nl` | integer | `99999` | `0`: the units above the unconformities are homogeneous (e.g., seawater); otherwise layered. |
| `unconf_refl_shape` | character | `'random'` | Reflector shape above the unconformities. |
| `unconf_refl_height` | real(2) | `[0.0, 5.0]` | Reflector height range above the unconformities. |
| `unconf_refl_slope` | real | `-2.5` | Reflector slope above the unconformities. |
| `unconf_refl_smooth` | real | `10.0` | Reflector smoothing length above the unconformities. |
| `unconf_shape` | character | `'random'` | Unconformity surface: `random` (Perlin-noise topography), `meander_channel`, `meander_canyon`, `drainage_channel`, or `drainage_canyon`. |
| `unconf_channel_width` | real(2) | `[0.03, 0.08]` | Range of the channel or canyon width, approximately relative to the lateral model size. |
| `unconf_channel_sinuosity` | real | `1.0` | Meander maturity (`meander_*`); larger values give stronger meanders. |
| `unconf_channel_length` | real | `25000.0` | Channel centerline length in simulation units (`meander_*`), which sets the number of bends. |
| `unconf_channel_density` | real(2) | `[0.02, 0.08]` | Range of the drainage density (`drainage_*`). |
| `unconf_topo` | real | `0.25` | Background topography relative to the incision depth; `0`: flat. |

## Salt

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `yn_salt` | logical | `.false.` | Insert salt bodies. |
| `nsalt` | integer | `1` | Number of salt bodies. |
| `salt_radius` | real(2) | `[0.0, 0.0]` | Range of the maximum salt radius: several bodies span this range, and a single body uses the lower bound. `[0, 0]`: `[0.05, 0.2]*n2` (3D: `[0.05, 0.2]*0.5*(n2 + n3)`). In 3D, the radius can reach twice this value. |
| `salt_radius_variation` | real | `0.7` | Relative variation of the salt radius with depth, from `0.75v` to `v`. |
| `salt_path_variation` | real | `5.0` | Scale of the lateral deviation of the salt center from a vertical line; the bodies use `0.75` to `1` times this value, growing with their radius. |
| `salt_nnode` | integer | `10` | Number of nodes of the salt boundaries along depth. |
| `salt_top_z` | real(2) | `[0.5, 0.8]` | Range of the depth of the salt top. |
| `salt_top_height` | real | `20.0` | Maximum relief of the salt top. |
| `salt_vp`, `salt_vs`, `salt_rho` | real | `5000.0`, `2850.0`, `2150.0` | Salt Vp, Vs (elastic models only), and density. |
| `salt_before_unconf` | logical | `.true.` | Salt formed before the unconformities, so that they erode it. |

## Karst

Karst caves overwrite `vp`, `vs`, and `rho` as salt does, and the labels are zero inside them.

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `yn_karst` | logical | `.false.` | Insert a karst cave system. |
| `karst_z` | real(2) | `[0.4, 0.9]` | Depth window of the passages. |
| `karst_npassage` | integer | `10` | Number of passages. |
| `karst_nctrl` | integer | `40` | Number of control points per passage; more points give longer passages. |
| `karst_radius` | real(2) | `[0.0, 0.0]` | Range of the mean passage radius. `[0, 0]`: `[0.015, 0.03]*n1`. |
| `karst_radius_variation` | real | `0.35` | Relative standard deviation of the radius along a passage, reduced to 0.4 times for branching passages. |
| `karst_tortuosity` | real | `0.6` | Sinuosity of the passages: `0` straight, `1` random walk. |
| `karst_connect` | real | `0.7` | Probability that a passage branches off an existing one. |
| `karst_vp`, `karst_vs`, `karst_rho` | real | `2500.0`, `1200.0`, `2000.0` | Cave Vp, Vs (elastic models only), and density. |
| `karst_before_unconf` | logical | `.true.` | Karst formed before the unconformities, so that they erode it. |

## Elastic model

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `yn_elastic` | logical | `.false.` | Generate `vs` and the PP, PS, SP, and SS images instead of `image`. |
| `vpvsratio` | real(2) | `[1.5, 1.8]` | Range of the Vp/Vs ratio; each layer draws its own ratio. |
| `angle_ps` | real(2) | `[10.0, 40.0]` | Range of the P-wave incidence angles (degrees) stacked in the PS and SP images. PP and SS use the ray parameters of P-wave incidence angles from 0° to 15°. |

## RGT and facies

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `yn_rgt` | logical | `.false.` | Output the relative geologic time `rgt`. |
| `yn_facies` | logical | `.false.` | Output the layer index `facies`. |

## Image

The image is the reflectivity convolved with a point spread function (PSF). By default, the PSF is a vertical wavelet times Gaussian tapers (`psf_sigma`), so steep reflectors and high-angle faults image weakly, as in field data: with the default parameters, a reflector keeps about 38% of its amplitude at a dip of 30 degrees and 6% at 60 degrees. A flat reflector is imaged with the amplitude of its reflection coefficient. A custom PSF is used as given.

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `wave` | character | `'ricker'` | Wavelet: `ricker`, `ricker_deriv`, `gaussian`, `gaussian_deriv`, `sinc`, `delta`, or `''` (no convolution). |
| `f0` | real | `150.0` | Frequency of the wavelet (Hz); the peak frequency for `ricker`. |
| `dt` | real | `1.0e-3` | Sampling interval of the wavelet along x1 (s). |
| `f0_bottom` | real | `0.0` | Center frequency at the model bottom (Hz); the image blends from `f0` at the top. `0`: off. |
| `wave_phase` | real | `0.0` | Constant phase rotation of the wavelet (degrees). |
| `wave_filt_freqs`, `wave_filt_amps` | real(:) | | Frequencies (Hz) and amplitudes that filter the `delta` wavelet, e.g., `[0, 10, 30, 40]` and `[0, 1, 1, 0]` for a band-pass. |
| `psf_sigma` | 2D: real(2)<br>3D: real(3) | 2D: `[5.0, 2.5]`<br>3D: `[5.0, 2.5, 2.5]` | Widths of the Gaussian tapers of the PSF along x1 and the lateral axes. |
| `custom_psf` | logical | `.false.` | Use the custom PSF in `psf`, of the model size. |
| `yn_dip_independent` | logical | `.false.` | Dip-independent imaging: reflectors keep their amplitude up to dips of about 60 degrees. |
| `noise_type` | character | `'normal'` | Noise: `normal`, `uniform`, or `exp` (smoothed random noise), `wavenumber` (random filtering of the image spectrum), or `migration` (noise like that of migrated images). |
| `noise_level` | real | `0.05` | Noise level: the peak ratio to the image for `normal`, `uniform`, and `exp`; the RMS ratio for `migration`; for `wavenumber`, the threshold of the random spectral mask, where larger values remove more of the spectrum. With `yn_conv_noise`, the ratio is to the reflectivity. |
| `noise_smooth` | 2D: real(2)<br>3D: real(3) | `[1.0, 1.0]`<br>3D: `[1.0, 1.0, 1.0]` | Smoothing lengths of the noise (for `wavenumber`, of the spectral mask), except `migration`. |
| `yn_conv_noise` | logical | `.false.` | Add the noise before the PSF convolution, except `migration`. |
| `noise_mix` | real(3) | `[0.35, 0.30, 0.15]` | Weights of the worm, swing, and background components of the `migration` noise. |
| `noise_worm_length` | real | `10.0` | Along-layer length of the worm noise. |
| `noise_swing_dip` | real(2) | `[25.0, 75.0]` | Range of the dips of the swing noise (degrees). |
| `noise_swing_nfam` | integer | `2` | Number of dip families of the swing noise. |
| `noise_swing_direction` | character | `'dual'` | `dual`: the swings dip both ways; `single`: one way. |
| `illum_level` | real | `0.0` | Standard deviation of the logarithm of a smooth illumination gain. `0`: off. |
| `illum_smooth` | real | `30.0` | Correlation length of the illumination gain. |
| `jitter_shift` | real(2) | `[0.0, 0.0]` | Maximum smooth static shift and standard deviation of the random static shift of the traces. |
| `jitter_gain` | real | `0.0` | Standard deviation of a random gain of each trace. |

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

## Outputs

The output arrays have the size `(n1, n2)` in 2D and `(n1, n2, n3)` in 3D.

| Array | Description |
| --- | --- |
| `vp`, `rho` | Vp (m/s) and density (kg/m³). |
| `vs` | Vs (m/s), with `yn_elastic`. |
| `image` | Seismic image of acoustic models. |
| `image_pp`, `image_ps`, `image_sp`, `image_ss` | Seismic images of elastic models, with `yn_elastic`. |
| `rgt` | Relative geologic time, larger for older layers, with `yn_rgt`. |
| `facies` | Layer index, with `yn_facies`. |
| `fault`, `fault_dip`, `fault_disp`, `fault_strike` (3D), `fault_rake` (3D) | Fault index, dip, displacement (positive for normal and negative for reverse faults), strike, and rake (angles in degrees), with `yn_fault`. |
| `salt`, `karst` | `1` inside salt or caves, with `yn_salt` or `yn_karst`. |
| `psf` | Point spread function of the image; for elastic models, of the SS image, and with `f0_bottom`, at the top. |
