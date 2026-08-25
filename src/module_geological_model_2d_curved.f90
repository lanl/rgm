!
! © 2024-2026. Triad National Security, LLC. All rights reserved.
!
! This program was produced under U.S. Government contract 89233218CNA000001
! for Los Alamos National Laboratory (LANL), which is operated by
! Triad National Security, LLC for the U.S. Department of Energy/National Nuclear
! Security Administration. All rights in the program are reserved by
! Triad National Security, LLC, and the U.S. Department of Energy/National
! Nuclear Security Administration. The Government is granted for itself and
! others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
! license in this material to reproduce, prepare derivative works,
! distribute copies to the public, perform publicly and display publicly,
! and to permit others to do so.
!
! Author:
!    Kai Gao, kaigao@lanl.gov
!

module geological_model_2d_curved

    !
    ! 2D random geological model generator (type rgm2_curved).
    !
    ! The generator follows a geologically motivated workflow:
    ! build layered velocity/density (and optionally elastic) models between two
    ! bounding curves -> deform the model with faults -> optionally insert
    ! salt bodies and unconformities -> compute reflectivity -> convolve with a
    ! wavelet to produce a synthetic migration-like image, together with
    ! pixel-wise labels (fault index/dip/displacement, relative geological
    ! time, facies, salt).
    !
    ! Faults are curved: the dip can vary with depth (listric faults via
    ! delta_dip).
    !
    ! The displacement of a fault can spatially vary along the fault when
    ! yn_vary_disp = .true. In this case, the displacement follows a slip patch
    ! that is maximum at the patch center and diminishes to zero at the patch
    ! boundary (fault tips), so a fault can die out within the model;
    ! fault_disp stores the spatially varying displacements.
    ! Additionally, the displacement can decay away from the fault when
    ! yn_disp_decay = .true., mimicking the deformation halo of a finite slip
    ! patch (fault drag, rollover, and blind-fault folding).
    !

    use libflit
    use geological_model_utility
    use geological_model_meander
    use geological_model_drainage
    use geological_model_karst

    implicit none

    ! 2D random geological model with curved faults
    type rgm2_curved

        !==============================================================================================
        !> Number of grid points along x1
        integer :: n1 = 128
        !> Number of grid points along x2
        integer :: n2 = 128
        !> Number of faults
        integer :: nf = 4
        !> Number of reflectors
        integer :: nl = 20
        !> Random seed
        integer :: seed = -1
        !> sigma (in terms of grid number) of Gaussian filter for smoothing reflectors along x2
        real :: refl_smooth = 20.0
        real :: refl_smooth_top = 20.0
        !> Linear slope (in terms of grid number) added to reflectors
        real :: refl_slope = 0.0
        real :: refl_slope_top = 0.0
        !> Center frequency of source wavelet for convolving reflectors
        real :: dt = 1.0e-3
        real :: f0 = 150.0
        !> Fracture width (in terms of grid number)
        real :: fwidth = 2.0
        !> Range of fault dip angles and displacements; for displacement/throw,
        !> + -> normal fault and - -> reverse fault
        real, allocatable, dimension(:) :: dip, disp
        !> Dip increase/descrease at the top compared with the top
        real, dimension(1:2) :: delta_dip = [15.0, 30.0]

        !> Whether to spatially vary the displacement along each fault.
        !> When = .true., the displacement follows an elliptical slip patch
        !> along the fault: maximum at the patch center and diminishing to zero
        !> at the patch boundary (fault tips), so a fault can die out within
        !> the model rather than always cutting through the entire section, and
        !> the block deforms gently near the fault tips.
        !> In this case, fault_disp stores the spatially varying displacements.
        logical :: yn_vary_disp = .false.
        !> Range of the along-depth semi-axis of the slip patch,
        !> relative to the model depth n1
        real, dimension(1:2) :: disp_radius_dip = [0.6, 1.2]
        !> Range of the center depth of the slip patch, relative to the model
        !> depth n1. Note that all the displacement tapering (and therefore the
        !> deformation of layers) occurs between the patch center and the patch
        !> boundary; setting a shallower center together with a smaller
        !> disp_radius_dip keeps the fault tips and the associated deformation
        !> away from the deep, high-velocity section
        real, dimension(1:2) :: disp_center_dip = [0.25, 0.75]

        !> Whether to decay the fault displacement away from the fault.
        !> When = .false. (default), the displaced block shifts rigidly
        !> (correct for through-going faults). When = .true., the displacement
        !> of the shifted block decays with a Gaussian profile of the distance
        !> to the fault, mimicking the deformation halo of a finite slip patch
        !> (fault drag, rollover, blind-fault folding); the decay is applied
        !> only within the displaced block, and the displacement at the fault
        !> itself is unaffected.
        logical :: yn_disp_decay = .false.
        !> Range of the displacement decay width (the Gaussian standard deviation
        !> of the decay profile), relative to the model depth n1
        real, dimension(1:2) :: disp_decay_width = [0.5, 1.0]
        !> sigmas of Gaussian filter for smoothing the generated noise
        real, dimension(1:2) :: noise_smooth = [1.0, 1.0]
        !> Level of noise in terms of max(abs(noise))/max(abs(image))
        real :: noise_level = 0.05
        !> Source wavelet for reflector convolution, can be one of
        !> ricker, ricker_deriv, gaussian, gaussian_deriv, sinc, delta
        character(len=24) :: wave = 'ricker'
        !> Shape of reflectors, can be one of
        !> random, gaussian, cauchy, perlin, custom
        !> When = custom, must specifiy the array refl
        character(len=24) :: refl_shape = 'random'
        character(len=24) :: refl_shape_top = 'random'
        real, allocatable, dimension(:) :: refl
        real, allocatable, dimension(:) :: refl_top
        !> Number of Gaussians for refl_shape = gaussian
        integer :: ng = 2
        !> Range of reflector's heights
        real, dimension(1:2) :: refl_height = 20.0
        real, dimension(1:2) :: refl_height_top = 10.0
        !> Range of Gaussian standard devision for refl_shape = gaussian
        real, dimension(1:2) :: refl_sigma2 = [0.0, 0.0]
        !> Range of Gaussian mean for refl_shape = gaussian
        real, dimension(1:2) :: refl_mu2 = [0.0, 0.0]
        !> The vertical thickness of layers varies from [1 - lwv, 1 + lwv] of average layer thickness
        real :: lwv = 0.25
        !> The horizontal thickness variation
        real :: lwh = 0.1
        !> Secondary reflector smoothing
        real :: secondary_refl_smooth = 10.0

        !> Whether or not to compute relative geological time
        !> if = .true., the array rgt will be filled with computed RGT
        logical :: yn_rgt = .false.
        !> Whether or not to compute facies (or piecewise constant random perturbations)
        !> if = .true., the array facies will be filled with computed facies
        logical :: yn_facies = .false.
        !> Whether or not to output faults
        !> if = .true., the array fault will be filled with 1, 2, 3, ..., nf, indicating
        !> the numbered faults within the model; if only fault probability is needed, then
        !> clip(fault, 0, 1) will do the work. Meanwhile, fault_dip and fault_disp will be
        !> filled with the dip angles and displacements associated with the faults.
        !> Note that this value does not affect the insertion of faults into the model set by nf, ...
        !> i.e., if = .false., then it will not fill the fault arrays, but the model will still be
        !> a faulted model.
        logical :: yn_fault = .true.
        !> Arrays for holding results
        real, allocatable, dimension(:, :) :: image, rgt, facies, fault, fault_dip, fault_disp
        !> sigmas of Gaussians for tapering the source wavelet the vertical and horizontal directions
        real, dimension(1:2) :: psf_sigma = [5.0, 2.5]
        !> Array for holding a custom point spread function
        real, allocatable, dimension(:, :) :: psf
        !> Whether or not to set a custom point spread function
        !> If yes then psf must be given with a dimension of (n1, n2)
        logical :: custom_psf = .false.
        !> Type of random noise, can be normal, uniform, or exp
        character(len=12) :: noise_type = 'normal'
        !> Convolve the psf with noise as well in addition to reflector image
        !> If = .false., then noise will be added after reflector-psf convolution
        logical :: yn_conv_noise = .false.
        !> Set faults to be with (quasi) regular spacing and dips
        logical :: yn_regular_fault = .false.
        !> For regularly spaced faults, whether to group faults with distinct dips into a group in space
        logical :: yn_group_faults = .false.
        !> Source wavelet filtering frequencies and coefficients, e.g.,
        !> amps = [0, 1, 1, 0] and freqs = [0, 10, 30, 40]
        !> is a band-pass filtering of 0, 10, 30, 40 Hz.
        real, allocatable, dimension(:) :: wave_filt_freqs, wave_filt_amps
        !> Min value for scaling the facies
        real :: vmin = 2000.0
        !> Max value for scaling the facies; after scaling, the facies will fall in [vmin, vmax]
        real :: vmax = 4000.0
        !> Velocity perturbation of layers
        real :: delta_v = 500.0

        !==============================================================================================
        !> Number of unconformity interfaces
        integer :: unconf = 0
        !> Range of depth of unconformity interfaces in terms of fraction of the entire depth,
        !> which must fall in [0, 1]; smaller values represent shallower unconformity interfaces
        real, dimension(1:2) :: unconf_z = [0.0, 0.5]
        !> Range of height of unconformity interfaces in depth, smaller values represent flatter unconformity interfaces
        real, dimension(1:2) :: unconf_height = [5.0, 15.0]
        !> Number of reflectors above the unconformity interfaces; when
        !> = 1 the region above the unconformity represents water
        integer :: unconf_nl = 99999
        !> Smoothing unconformity surfaces
        real :: unconf_smooth = 0.0
        !> Reflector height above the unconformity
        real, dimension(1:2) :: unconf_refl_height = [0.0, 5.0]
        !> Relector slope above the unconformity
        real :: unconf_refl_slope = -2.5
        !> Reflector smoothing
        real :: unconf_refl_smooth = 10.0
        !> Reflector shape
        character(len=12) :: unconf_refl_shape = 'random'
        !> Shape of the unconformity surfaces, can be one of
        !>      random - random (Perlin-noise) topography (default)
        !>      meander_channel - meandering river channels
        !>      meander_canyon - meandering incised canyon
        !>      drainage_channel - dendritic drainage channel network
        !>      drainage_canyon - canyon carved by a dendritic drainage network
        !> For the channel/canyon shapes, a map-view geomorphological simulation
        !> is run on an auxiliary grid and a random cross-flow section of its
        !> erosional depth map is used as the unconformity topography;
        !> unconf_height still bounds the total erosional relief
        character(len=24) :: unconf_shape = 'random'
        !> Channel/canyon width as an (approximate) fraction of the lateral
        !> model extent, drawn per unconformity surface
        real, dimension(1:2) :: unconf_channel_width = [0.03, 0.08]
        !> Meander maturity; scales the number of migration iterations
        real :: unconf_channel_sinuosity = 1.0
        !> Total centerline length of the meandering channel/canyon in the
        !> internal units of the migration simulation; since the simulated
        !> map is fit to the auxiliary grid, this effectively controls how
        !> many meander bends appear across the map (and thus how many
        !> channel crossings the sampled 2D section can contain). Values
        !> larger than the default add more bends of the same shape; smaller
        !> values give fewer, larger bends (clamped from below to keep enough
        !> centerline nodes for the migration to develop realistic meanders)
        real :: unconf_channel_length = 25000.0
        !> Drainage density (fraction of cells covered by channels) for the
        !> drainage shapes, drawn per unconformity surface
        real, dimension(1:2) :: unconf_channel_density = [0.02, 0.08]
        !> Amplitude of the background (interfluve) topography relative to the
        !> channel/canyon incision depth; 0 = flat background between channels
        real :: unconf_topo = 0.25

        !==============================================================================================
        !> Salt body
        logical :: yn_salt = .false.
        !> Number of salt bodies
        integer :: nsalt = 1
        !> Range of max horizontal radius of salt bodies in grid number
        real, dimension(1:2) :: salt_radius = [0.0, 0.0]
        !> Range of max depth of salt body top interfaces in fraction of n1
        real, dimension(1:2) :: salt_top_z = [0.5, 0.8]
        !> Salt body velocity
        real :: salt_vp = 5000.0
        !> Salt body density
        real :: salt_rho = 2150.0
        !> Array for holding the salt bodies
        real, allocatable, dimension(:, :) :: salt
        !> Max height of the random salt body top in grid number
        real :: salt_top_height = 20.0
        !> Degree of salt radius variation
        real :: salt_radius_variation = 0.7
        !> Maximum devaition of salt center from a vertical line
        real :: salt_path_variation = 5.0
        !> Number of nodes to form the salt vertical profile
        integer :: salt_nnode = 10

        !==============================================================================================
        !> Whether to insert a karst cave system into the model
        logical :: yn_karst = .false.
        !> Depth window (relative to n1) containing the karst passages
        real, dimension(1:2) :: karst_z = [0.4, 0.9]
        !> Number of karst passages
        integer :: karst_npassage = 10
        !> Number of control points per passage; more points -> longer passages
        integer :: karst_nctrl = 40
        !> Range of the mean passage radius (in grid points);
        !> [0, 0] (default) = automatically set to [0.015, 0.03]*n1
        real, dimension(1:2) :: karst_radius = [0.0, 0.0]
        !> Fractional standard deviation of the radius along a passage
        real :: karst_radius_variation = 0.35
        !> Sinuosity of the passages: 0 = straight, 1 = fully random walk
        real :: karst_tortuosity = 0.6
        !> Probability that a new passage branches off an existing one
        real :: karst_connect = 0.7
        !> Medium properties filling the karst caves; the defaults represent
        !> water/sediment-filled voids (low-velocity anomalies)
        real :: karst_vp = 2500.0
        real :: karst_vs = 1200.0
        real :: karst_rho = 2000.0
        !> Is karst before or after unconformity?
        logical :: karst_before_unconf = .true.
        !> Output mask: 1 = karst
        real, allocatable, dimension(:, :) :: karst

        !==============================================================================================
        !> Elastic
        logical :: yn_elastic = .false.
        !> Vp/Vs ratios
        real, dimension(1:2) :: vpvsratio = [1.5, 1.8]
        !> Elastic models
        real, allocatable, dimension(:, :) :: vp, vs, rho
        real :: rho_a = 310.0, rho_b = 0.25, rho_c = 0.0
        !> Elastic images
        real, allocatable, dimension(:, :) :: image_pp, image_ps, image_sp, image_ss
        !> Salt body's Vs
        real :: salt_vs = 4400.0
        !> Is salt before or after unconformity?
        logical :: salt_before_unconf = .true.

    contains

        procedure, private :: create_psf
        procedure, private :: generate_image
        procedure, private :: generate_image_elastic
        procedure, public :: generate => generate_2d

    end type rgm2_curved

    private
    public :: rgm2_curved

contains

    subroutine create_psf(this, n1, n2, freq)

        class(rgm2_curved), intent(inout) :: this
        integer :: n1, n2
        real, intent(in), optional :: freq

        real, allocatable, dimension(:) :: wavelet, psf1, psf2
        real, allocatable, dimension(:, :) :: psf
        real :: f0, wt
        integer :: i, j

        if (present(freq)) then
            f0 = freq
        else
            f0 = this%f0
        end if

        wavelet = zeros(n1)
        !$omp parallel do private(i, wt)
        do i = 1, n1
            wt = (i - 1.0 - (n1 - 1.0)/2.0)*this%dt
            select case (this%wave)
                case ('ricker')
                    wavelet(i) = ricker_wavelet(wt, f0)
                case ('ricker_deriv')
                    wavelet(i) = ricker_deriv_wavelet(wt, f0)
                case ('gaussian_deriv')
                    wavelet(i) = gaussian_deriv_wavelet(wt, f0)
                case ('gaussian')
                    wavelet(i) = gaussian_wavelet(wt, f0)
                case ('sinc')
                    wavelet(i) = sinc_wavelet(wt, f0)
                case default
                    wavelet(i) = ricker_wavelet(wt, f0)
            end select
        end do
        !$omp end parallel do

        if (this%wave == 'delta') then
            wavelet = 0
            wavelet((n1 + 1)/2) = 1.0
            if (allocated(this%wave_filt_freqs)) then
                call assert(size(this%wave_filt_freqs) == size(this%wave_filt_amps))
                wavelet = fourier_filt(wavelet, this%dt, this%wave_filt_freqs, this%wave_filt_amps)
            end if
        end if

        wavelet = wavelet/norm2(wavelet)

        if (.not. this%custom_psf) then
            psf = zeros(n1, n2)
            psf1 = zeros(n1)
            psf2 = zeros(n2)
            if (this%psf_sigma(1) == 0) then
                !$omp parallel do private(i)
                do i = 1, n1
                    psf1(i) = exp(-0.5*(i - 1.0 - (n1 - 1.0)/2.0)**2)
                end do
                !$omp end parallel do
                where (psf1 < maxval(psf1))
                    psf1 = 0.0
                end where
            else
                !$omp parallel do private(i)
                do i = 1, n1
                    psf1(i) = exp(-0.5*(i - 1.0 - (n1 - 1.0)/2.0)**2/this%psf_sigma(1)**2)
                end do
                !$omp end parallel do
            end if
            if (this%psf_sigma(2) == 0) then
                !$omp parallel do private(j)
                do j = 1, n2
                    psf2(j) = exp(-0.5*(j - 1.0 - (n2 - 1.0)/2.0)**2)
                end do
                !$omp end parallel do
                where (psf2 < maxval(psf2))
                    psf2 = 0.0
                end where
            else
                !$omp parallel do private(j)
                do j = 1, n2
                    psf2(j) = exp(-0.5*(j - 1.0 - (n2 - 1.0)/2.0)**2/this%psf_sigma(2)**2)
                end do
                !$omp end parallel do
            end if
            !$omp parallel do private(i, j) collapse(2)
            do j = 1, n2
                do i = 1, n1
                    psf(i, j) = wavelet(i)*psf1(i)*psf2(j)
                end do
            end do
            !$omp end parallel do
            this%psf = psf/norm2(psf)
        else
            call assert(size(this%psf, 1) == this%n1 .and. size(this%psf, 2) == this%n2, &
                ' Error: Shape of custom psf must be n1 x n2')
        end if

    end subroutine create_psf

    subroutine generate_image_elastic(this, vp, vs, rho)

        class(rgm2_curved), intent(inout) :: this
        real, dimension(:, :), intent(in) :: vp, vs, rho

        integer :: n1, n2, i, j, l
        real, allocatable, dimension(:) :: rfc
        real, allocatable, dimension(:, :) :: ww

        n1 = this%n1
        n2 = this%n2

        this%image_pp = zeros(n1, n2)
        this%image_ps = zeros(n1, n2)
        this%image_sp = zeros(n1, n2)
        this%image_ss = zeros(n1, n2)

        rfc = zeros(4)
        !$omp parallel do private(i, j, l, rfc)
        do j = 1, n2
            do i = 1, n1
                rfc = 0
                do l = 0, 5
                    rfc = rfc + elastic_reflection_coefs(real(sin((l*3.0)*const_deg2rad)/vp(i, j)), &
                        vp(i, j), vs(i, j), rho(i, j), vp(i + 1, j), vs(i + 1, j), rho(i + 1, j))
                end do
                rfc = rfc/6.0
                this%image_pp(i, j) = rfc(1)
                this%image_ps(i, j) = rfc(2)
                this%image_sp(i, j) = rfc(3)
                this%image_ss(i, j) = rfc(4)
            end do
        end do
        !$omp end parallel do

        ! Add random noise
        if (this%noise_level /= 0 .and. this%yn_conv_noise) then
            select case (this%noise_type)

                case ('normal', 'gaussian', 'uniform', 'exp')
                    ww = gauss_filt(random(n1, n2, dist=this%noise_type, seed=safe_seed(int(this%seed, 8)*17)), this%noise_smooth)
                    ww = ww - mean(ww)
                    ww = ww/maxval(abs(ww))*this%noise_level*maxval(abs(this%image_pp))
                    this%image_pp = this%image_pp + ww

                    ww = gauss_filt(random(n1, n2, dist=this%noise_type, seed=safe_seed(int(this%seed, 8)*17 + 1)), this%noise_smooth)
                    ww = ww - mean(ww)
                    ww = ww/maxval(abs(ww))*this%noise_level*maxval(abs(this%image_ps))
                    this%image_ps = this%image_ps + ww

                    ww = gauss_filt(random(n1, n2, dist=this%noise_type, seed=safe_seed(int(this%seed, 8)*17 + 2)), this%noise_smooth)
                    ww = ww - mean(ww)
                    ww = ww/maxval(abs(ww))*this%noise_level*maxval(abs(this%image_sp))
                    this%image_sp = this%image_sp + ww

                    ww = gauss_filt(random(n1, n2, dist=this%noise_type, seed=safe_seed(int(this%seed, 8)*17 + 3)), this%noise_smooth)
                    ww = ww - mean(ww)
                    ww = ww/maxval(abs(ww))*this%noise_level*maxval(abs(this%image_ss))
                    this%image_ss = this%image_ss + ww

                case ('wavenumber')
                    this%image_pp = this%image_pp + noise_wavenumber(this%image_pp, this%noise_level, this%noise_smooth, safe_seed(int(this%seed, 8)*17))
                    this%image_ps = this%image_ps + noise_wavenumber(this%image_ps, this%noise_level, this%noise_smooth, safe_seed(int(this%seed, 8)*17 + 1))
                    this%image_sp = this%image_sp + noise_wavenumber(this%image_sp, this%noise_level, this%noise_smooth, safe_seed(int(this%seed, 8)*17 + 2))
                    this%image_ss = this%image_ss + noise_wavenumber(this%image_ss, this%noise_level, this%noise_smooth, safe_seed(int(this%seed, 8)*17 + 3))

            end select
        end if

        ! Source wavelet
        if (this%wave /= '') then

            call this%create_psf(n1, n2, this%f0)
            this%image_pp = conv(this%image_pp, this%psf, 'same')
            call this%create_psf(n1, n2, this%f0*mean(0.5*(ones(n1, n2) + vp(1:this%n1, :)/vs(1:this%n1, :))))
            this%image_ps = conv(this%image_ps, this%psf, 'same')
            this%image_sp = conv(this%image_sp, this%psf, 'same')
            call this%create_psf(n1, n2, this%f0*mean(vp(1:this%n1, :)/vs(1:this%n1, :)))
            this%image_ss = conv(this%image_ss, this%psf, 'same')

        end if

        ! Add random noise
        if (this%noise_level /= 0 .and. (.not. this%yn_conv_noise)) then
            select case (this%noise_type)

                case ('normal', 'gaussian', 'uniform', 'exp')
                    ww = gauss_filt(random(n1, n2, dist=this%noise_type, seed=safe_seed(int(this%seed, 8)*17)), this%noise_smooth)
                    ww = ww - mean(ww)
                    ww = ww/maxval(abs(ww))*this%noise_level*maxval(abs(this%image_pp))
                    this%image_pp = this%image_pp + ww

                    ww = gauss_filt(random(n1, n2, dist=this%noise_type, seed=safe_seed(int(this%seed, 8)*17 + 1)), this%noise_smooth)
                    ww = ww - mean(ww)
                    ww = ww/maxval(abs(ww))*this%noise_level*maxval(abs(this%image_ps))
                    this%image_ps = this%image_ps + ww

                    ww = gauss_filt(random(n1, n2, dist=this%noise_type, seed=safe_seed(int(this%seed, 8)*17 + 2)), this%noise_smooth)
                    ww = ww - mean(ww)
                    ww = ww/maxval(abs(ww))*this%noise_level*maxval(abs(this%image_sp))
                    this%image_sp = this%image_sp + ww

                    ww = gauss_filt(random(n1, n2, dist=this%noise_type, seed=safe_seed(int(this%seed, 8)*17 + 3)), this%noise_smooth)
                    ww = ww - mean(ww)
                    ww = ww/maxval(abs(ww))*this%noise_level*maxval(abs(this%image_ss))
                    this%image_ss = this%image_ss + ww

                case ('wavenumber')
                    this%image_pp = this%image_pp + noise_wavenumber(this%image_pp, this%noise_level, this%noise_smooth, safe_seed(int(this%seed, 8)*17))
                    this%image_ps = this%image_ps + noise_wavenumber(this%image_ps, this%noise_level, this%noise_smooth, safe_seed(int(this%seed, 8)*17 - 1))
                    this%image_sp = this%image_sp + noise_wavenumber(this%image_sp, this%noise_level, this%noise_smooth, safe_seed(int(this%seed, 8)*17 - 2))
                    this%image_ss = this%image_ss + noise_wavenumber(this%image_ss, this%noise_level, this%noise_smooth, safe_seed(int(this%seed, 8)*17 - 3))

            end select
        end if

    end subroutine generate_image_elastic

    subroutine generate_image(this, vp, rho)

        class(rgm2_curved), intent(inout) :: this
        real, dimension(:, :), intent(in) :: vp, rho

        integer :: n1, n2, i, j
        real, allocatable, dimension(:, :) :: ww

        n1 = this%n1
        n2 = this%n2

        this%image = zeros(n1, n2)

        !$omp parallel do private(i, j)
        do j = 1, n2
            do i = 1, n1
                this%image(i, j) = (vp(i + 1, j)*rho(i + 1, j) - vp(i, j)*rho(i, j)) &
                    /(vp(i + 1, j)*rho(i + 1, j) + vp(i, j)*rho(i, j))
            end do
        end do
        !$omp end parallel do

        ! Add random noise
        if (this%noise_level /= 0 .and. this%yn_conv_noise) then
            select case (this%noise_type)

                case ('normal', 'gaussian', 'uniform', 'exp')
                    ww = gauss_filt(random(n1, n2, dist=this%noise_type, seed=safe_seed(int(this%seed, 8)*17)), this%noise_smooth)
                    ww = ww - mean(ww)
                    ww = ww/maxval(abs(ww))*this%noise_level*maxval(abs(this%image))
                    this%image = this%image + ww

                case ('wavenumber')
                    this%image = this%image + noise_wavenumber(this%image, this%noise_level, this%noise_smooth, safe_seed(int(this%seed, 8)*17))

            end select
        end if

        ! Source wavelet
        if (this%wave /= '') then
            call this%create_psf(n1, n2, this%f0)
            this%image = conv(this%image, this%psf, 'same')
        end if

        if (this%noise_level /= 0 .and. (.not. this%yn_conv_noise)) then
            select case (this%noise_type)

                case ('normal', 'gaussian', 'uniform', 'exp')
                    ww = gauss_filt(random(n1, n2, dist=this%noise_type, seed=safe_seed(int(this%seed, 8)*17)), this%noise_smooth)
                    ww = ww - mean(ww)
                    ww = ww/maxval(abs(ww))*this%noise_level*maxval(abs(this%image))
                    this%image = this%image + ww

                case ('wavenumber')
                    this%image = this%image + noise_wavenumber(this%image, this%noise_level, this%noise_smooth, safe_seed(int(this%seed, 8)*17))

            end select
        end if

    end subroutine generate_image

    subroutine generate_2d(this)

        class(rgm2_curved), intent(inout) :: this

        if (this%nf == 0) then
            this%yn_fault = .false.
        end if

        if (this%unconf == 0) then
            call generate_2d_geological_model(this)
        else
            call generate_2d_unconformal_geological_model(this)
        end if

    end subroutine generate_2d

    subroutine generate_2d_geological_model(this)

        type(rgm2_curved), intent(inout) :: this

        real, allocatable, dimension(:, :) :: w, f, t, m
        real, allocatable, dimension(:, :) :: ww, ff, cf, tt, mm
        integer :: nf, nl, n1, n2, i, j, fi, ne1, ne2, newi, newj
        real, allocatable, dimension(:) :: disp, dip, f2, r, rt, vv
        real :: fwidth
        real, allocatable, dimension(:) :: mu, sigma, height
        real, dimension(1:2) :: mu2, sigma2
        real, allocatable, dimension(:) :: sumdisp
        real, allocatable, dimension(:, :) :: lz
        real, allocatable, dimension(:) :: plw, delta_dip
        real, allocatable, dimension(:, :) :: fdip, fdisp, ffdip, ffdisp
        integer :: l, nsf
        real :: b, b_prev, dxys, dist, dist_prev, theta, x0
        real, allocatable, dimension(:) :: disp_az, disp_zc, decay_w, bs
        real :: alpha, dloc
        integer :: ia, ja
        real :: wa, wb
        logical :: yn_mark
        logical, allocatable, dimension(:, :) :: fblock
        real, allocatable, dimension(:, :) :: dips
        real, allocatable, dimension(:) :: x1, x2, rds, topz, rc, vds, pds, xs
        integer :: nd, isalt
        real, allocatable, dimension(:) :: gmax, salt_radius
        real :: tp
        type(fractal_noise_1d) :: pn
        real, allocatable, dimension(:, :) :: pxy
        real :: thick
        real, allocatable, dimension(:, :) :: vp, vs, rho
        real, allocatable, dimension(:) :: rfc
        real :: m1, m2, m3

        fwidth = this%fwidth

        n1 = this%n1
        n2 = this%n2
        nf = this%nf

        if (this%nf /= 0 .and. this%yn_regular_fault) then
            call assert(this%nf >= 2, ' <generate_3d_geological_model> Error: for yn_regular_fault = .true., nf must >= 2')
        end if

        if (nf >= 1) then

            if (allocated(this%dip)) then
                call assert(size(this%dip) >= 2 .and. mod(size(this%dip), 2) == 0, &
                    ' <generate_2d_geological_model> Error: fault dip is not specified properly. ')
            else
                this%dip = [70.0, 110.0]
            end if
            if (allocated(this%disp)) then
                call assert(size(this%disp) >= 2 .and. mod(size(this%disp), 2) == 0 .and. size(this%disp) == size(this%dip), &
                    ' <generate_2d_geological_model> Error: fault displacement is not specified properly. ')
            else
                this%disp = tile([5.0, 30.0], size(this%dip)/2)
            end if

            if (this%yn_regular_fault) then

                dip = [random(nint(nf/2.0), range=[0.95, 1.05]*this%dip(1), seed=this%seed)*const_deg2rad, &
                    random(nf - nint(nf/2.0), range=[0.95, 1.05]*this%dip(2), seed=this%seed)*const_deg2rad]
                disp = [random(nint(nf/2.0), range=[0.9, 1.1]*this%disp(1), seed=safe_seed(int(this%seed, 8)*2)), &
                    random(nf - nint(nf/2.0), range=[0.9, 1.1]*this%disp(2), seed=safe_seed(int(this%seed, 8)*2))]

            else

                ! Dip angles and fault displacements
                nsf = size(this%dip)/2
                dip = random(ceiling(nf*1.0/nsf), range=this%dip(1:2), seed=this%seed)
                disp = random(ceiling(nf*1.0/nsf), range=this%disp(1:2), seed=safe_seed(int(this%seed, 8)*2))
                do i = 2, nsf
                    dip = [dip, random(ceiling(nf*1.0/nsf), range=this%dip((i - 1)*2 + 1:i*2), seed=this%seed)]
                    disp = [disp, random(ceiling(nf*1.0/nsf), range=this%disp((i - 1)*2 + 1:i*2), seed=safe_seed(int(this%seed, 8)*2))]
                end do
                dip = dip(1:nf)
                dip = dip*const_deg2rad

                disp = disp(1:nf)
                disp(1:nf:2) = -disp(1:nf:2)

            end if
        else
            dip = zeros(1)
            disp = zeros(1)
        end if

        ! Compute the number of padding grid points
        sumdisp = disp*sin(dip)
        m1 = max(sum(sumdisp, mask=sumdisp > 0), -sum(sumdisp, mask=sumdisp < 0))
        m2 = max(abs(this%refl_slope), abs(this%refl_slope_top))
        m3 = max(maxval(abs(this%refl_height)), maxval(abs(this%refl_height_top)))
        ne1 = ceiling(max(m1, m2) + m3)
        n1 = n1 + 2*ne1

        sumdisp = disp*cos(dip)
        ne2 = ceiling(max(sum(sumdisp, mask=sumdisp > 0), -sum(sumdisp, mask=sumdisp < 0)))
        n2 = n2 + 2*ne2

        ! Compute the top and bottom fault dip angles
        ! Note that to ensure the faults have proper dip angle range within the final cropped image,
        ! here I add some extra degrees to the begin and end dip angles
        dips = zeros(n1, nf)
        delta_dip = random(nf, range=this%delta_dip, seed=safe_seed((int(this%seed, 8)*5 + 1)/2))*const_deg2rad
        do i = 1, nf
            if (dip(i) <= const_pi_half) then
                dips(:, i) = [ones(ne1)*dip(i), &
                    linspace(dip(i), dip(i) - delta_dip(i)*(1.0 + ne1*1.0/n1), n1 - ne1)]
                dips(:, i) = clip(dips(:, i), 0.0, real(const_pi_half))
            else
                dips(:, i) = [ones(ne1)*dip(i), &
                    linspace(dip(i), dip(i) + delta_dip(i)*(1.0 + ne1*1.0/n1), n1 - ne1)]
                dips(:, i) = clip(dips(:, i), real(const_pi_half), real(const_pi))
            end if
        end do

        ! Reflector's shape at the bottom (r) and at the top (rt)
        select case (this%refl_shape)

            case default
                r = random(n2, dist='normal', seed=safe_seed(int(this%seed, 8)*3))
                r = gauss_filt(r, this%refl_smooth)

            case ('random')
                r = random(n2, dist='normal', seed=safe_seed(int(this%seed, 8)*3))
                r = gauss_filt(r, this%refl_smooth)

            case ('gaussian', 'cauchy')
                if (maxval(this%refl_mu2) == 0) then
                    mu2 = [1.0, this%n2 - 1.0]
                else
                    mu2 = this%refl_mu2
                end if

                if (maxval(this%refl_sigma2) == 0) then
                    sigma2 = [0.05, 0.15]*n2
                else
                    sigma2 = this%refl_sigma2
                end if

                mu = random(this%ng, range=mu2, seed=safe_seed(int(this%seed, 8)*3))
                sigma = random(this%ng, range=sigma2, seed=safe_seed(int(this%seed, 8)*4))
                height = random(this%ng, range=this%refl_height, seed=safe_seed(int(this%seed, 8)*5))

                r = zeros(n2)
                do i = 1, this%ng
                    select case (this%refl_shape)
                        case ('gaussian')
                            r = r + rescale(gaussian(linspace(0.0, n2 - 1.0, n2), mu(i) + ne2, sigma(i)), [0.0, height(i)])
                        case ('cauchy')
                            r = r + rescale(cauchy(linspace(0.0, n2 - 1.0, n2), mu(i) + ne2, sigma(i)), [0.0, height(i)])
                    end select
                end do

            case ('perlin')
                pn%n1 = n2
                pn%seed = safe_seed(int(this%seed, 8)*3)
                r = gauss_filt(pn%generate(), this%refl_smooth)

            case ('custom')
                call assert(allocated(this%refl), &
                    ' <generate_2d_geological_model> Error: refl must be initialized. ')
                call assert(size(this%refl) == this%n2, &
                    '<generate_2d_geological_model> Error: size(refl) must = n2')
                r = pad(this%refl, [ne2, ne2], ['edge', 'edge'])

        end select

        ! If the shape of the top reflector is not 'same' with the bottom reflector
        if (this%refl_shape_top /= 'same') then

            select case (this%refl_shape_top)

                case default
                    rt = random(n2, dist='normal', seed=safe_seed(int(this%seed, 8)*3 - 1))
                    rt = gauss_filt(rt, this%refl_smooth_top)

                case ('random')
                    rt = random(n2, dist='normal', seed=safe_seed(int(this%seed, 8)*3 - 1))
                    rt = gauss_filt(rt, this%refl_smooth_top)

                case ('gaussian', 'cauchy')
                    if (maxval(this%refl_mu2) == 0) then
                        mu2 = [1.0, this%n2 - 1.0]
                    else
                        mu2 = this%refl_mu2
                    end if

                    if (maxval(this%refl_sigma2) == 0) then
                        sigma2 = [0.05, 0.15]*n2
                    else
                        sigma2 = this%refl_sigma2
                    end if

                    mu = random(this%ng, range=mu2, seed=safe_seed(int(this%seed, 8)*3 - 1))
                    sigma = random(this%ng, range=sigma2, seed=safe_seed(int(this%seed, 8)*4 - 1))
                    height = random(this%ng, range=this%refl_height_top, seed=safe_seed(int(this%seed, 8)*5 - 1))

                    rt = zeros(n2)
                    do i = 1, this%ng
                        select case (this%refl_shape)
                            case ('gaussian')
                                rt = rt + rescale(gaussian(linspace(0.0, n2 - 1.0, n2), mu(i) + ne2, sigma(i)), [0.0, height(i)])
                            case ('cauchy')
                                rt = rt + rescale(cauchy(linspace(0.0, n2 - 1.0, n2), mu(i) + ne2, sigma(i)), [0.0, height(i)])
                        end select
                    end do

                case ('perlin')
                    pn%n1 = n2
                    pn%seed = safe_seed(int(this%seed, 8)*3 - 1)
                    rt = gauss_filt(pn%generate(), this%refl_smooth_top)

                case ('custom')
                    call assert(allocated(this%refl_top), &
                        ' <generate_2d_geological_model> Error: refl_top must be initialized. ')
                    call assert(size(this%refl_top) == this%n2, &
                        '<generate_2d_geological_model> Error: size(refl_top) must = n2')
                    rt = pad(this%refl_top, [ne2 + 1, ne2 + 1], ['edge', 'edge'])

            end select

        else

            rt = r

        end if

        ! Rescale reflectors to their height
        r = rescale(r, this%refl_height*rov(r)/(rov(r(ne2 + 1:ne2 + this%n2)) + float_tiny))
        rt = rescale(rt, this%refl_height_top*rov(rt)/(rov(rt(ne2 + 1:ne2 + this%n2)) + float_tiny))

        ! Add slope
        !$omp parallel do private(j)
        do j = 1, n2
            r(j) = r(j) + (j - 1.0)*this%refl_slope/this%n2
            rt(j) = rt(j) + (j - 1.0)*this%refl_slope_top/this%n2
        end do
        !$omp end parallel do

        ! ... and get the final positions of top/bottom reflectors
        r = r - mean(r)
        rt = rt - mean(rt) + n1

        nl = nint(this%nl + this%nl*2.0*ne1/(n1 - 2*ne1))

        vv = random(nl - 1, seed=safe_seed(int(this%seed, 8)*6))*this%delta_v
        vv = linspace(this%vmax*(1.0 + ne1*1.0/this%n1), this%vmin*(1.0 - ne1*1.0/this%n1), nl - 1) + vv

        lz = zeros(nl, n2)
        call assert(abs(this%lwv) <= 1, ' Error: lwv must be 0 <= lwv <= 1. ')
        plw = random(nl, range=[1 - this%lwv, 1 + this%lwv], seed=safe_seed(int(this%seed, 8)*7))

        !$omp parallel do private(j)
        do j = 1, n2
            lz(:, j) = ginterp([0.0, (ne1 + 1.0)/n1, (n1 - ne1 - 1.0)/n1, 1.0], &
                [r(j), r(j) + ne1, rt(j) - ne1, rt(j)], &
                linspace(0.0, 1.0, nl), 'linear')
        end do
        !$omp end parallel do

        ! Compute layer interface positions in the vertical direction
        lz = deriv(lz, dim=1)
        !$omp parallel do private(j)
        do j = 1, n2
            lz(:, j) = r(j) + rescale(cumsum(lz(:, j)*plw), [0.0, rt(j) - r(j)])
        end do
        !$omp end parallel do

        ! Add horizontal layer thickness variation
        if (this%lwh > 0) then

            pxy = zeros(nl, n2)

            !$omp parallel do private(i, thick)
            do i = 1, nl - 1
                thick = (lz(i + 1, 1) - lz(i, 1))*this%lwh*2.0
                ! The seed here must be sign-consistent, e.g., when = -1, all must < 0; otherwise, all must > 0
                pxy(i, :) = gauss_filt(random(n2, seed=safe_seed(int(this%seed, 8)*nl - i)), this%secondary_refl_smooth)
                pxy(i, :) = rescale(pxy(i, :), [-thick, thick]*(0.5 + (nl - i + 0.0)/(nl - 1.0)*0.5))
            end do
            !$omp end parallel do

            lz = deriv(lz + pxy, dim=1)
            where (lz <= 0)
                lz = 1.0e-9
            end where
            !$omp parallel do private(j)
            do j = 1, n2
                lz(:, j) = r(j) + rescale(cumsum(lz(:, j)), [0.0, rt(j) - r(j)])
            end do
            !$omp end parallel do

        end if

        ! Stack layers to create a layer-constant velocity model
        w = zeros(n1, n2)
        if (this%yn_facies) then
            m = zeros(n1, n2)
        end if
        if (this%yn_rgt) then
            t = zeros(n1, n2)
        end if

        ! On the contrary, the following will generate a smoother velocity model by interpolation,
        ! but there will be in-between values
        nd = nint(n1*2.0/nl)
        rc = zeros((nl - 1)*nd)
        rfc = zeros((nl - 1)*nd)

        !$omp parallel do private(i, j, l, rc, rfc, tp)
        do j = 1, n2

            ! Velocity model is interpolated vertically based on reflectors
            ! Otherwise there will be staircases in the generated velocity model
            ! And correspondingly the images can has notable staircases which are unrealistic
            do l = 1, nl - 1
                tp = lz(l + 1, j) - lz(l, j)
                rc((l - 1)*nd + 1:l*nd) = linspace(lz(l, j) + tp/(nd + 1), lz(l + 1, j) - tp/(nd + 1), nd)
                rfc((l - 1)*nd + 1:l*nd) = vv(l)
            end do
            w(:, j) = ginterp(rc, rfc, linspace(n1 - 1.0, 0.0, n1), 'linear')

            ! Facies is piecewise constant and has distinct values for each layer
            if (this%yn_facies) then
                do i = 1, n1
                    loop_layer: do l = 1, nl - 1
                        if (n1 - i >= lz(l, j) .and. n1 - i < lz(l + 1, j)) then
                            m(i, j) = l
                            exit loop_layer
                        end if
                    end do loop_layer
                end do
            end if

            ! RGT is linearly interpolated based on the location of reflectors.
            if (this%yn_rgt) then
                t(:, j) = ginterp(n1 - 1.0 - lz(:, j), linspace(1.0, 0.0, nl), linspace(0.0, n1 - 1.0, n1), 'linear')
            end if

        end do
        !$omp end parallel do

        ! Add faults
        ww = zeros(n1, n2)
        f = zeros(n1, n2)
        ff = zeros(n1, n2)
        cf = zeros(n1, n2)
        fdip = zeros(n1, n2)
        fdisp = zeros(n1, n2)
        ffdip = zeros(n1, n2)
        ffdisp = zeros(n1, n2)
        if (this%yn_facies) then
            mm = zeros(n1, n2)
        end if
        if (this%yn_rgt) then
            tt = zeros(n1, n2)
        end if

        if (nf >= 1) then

            if (this%yn_regular_fault) then

                rc = random(nf - 1, range=[0.75, 1.25]*0.9*this%n2/(nf - 1.0), seed=safe_seed(int(this%seed, 8)*12))
                rc = rc*0.9*this%n2/sum(rc)
                rc = sort(rc)
                f2 = zeros(nf)
                f2(1) = ne2 + rand(range=[0.05, 0.1]*this%n2, seed=safe_seed(int(this%seed, 8)*12 - 1))
                f2(2:) = f2(1) + cumsum(rc)
                if (.not. this%yn_group_faults) then
                    f2 = random_permute(f2, seed=safe_seed(int(this%seed, 8)*12 - 2))
                end if

            else

                f2 = random(nf, range=[ne2 + 0.1*this%n2, n2 - ne2 - 0.1*this%n2], &
                    seed=safe_seed(int(this%seed, 8)*12), spacing=0.75*(n2 - 2*ne2 - 0.2*this%n2)/nf)

            end if

            ! Elliptical slip patches for spatially varying fault displacement;
            ! in 2D, the patch is an interval along the fault parameterized by depth
            disp_az = random(nf, range=this%disp_radius_dip, seed=safe_seed((int(this%seed, 8)*15 + 1)/2))*this%n1
            disp_zc = ne1 + random(nf, range=this%disp_center_dip, seed=safe_seed((int(this%seed, 8)*19 + 1)/2))*this%n1

            ! Widths of the Gaussian decay of displacement away from each fault
            decay_w = random(nf, range=this%disp_decay_width, seed=safe_seed((int(this%seed, 8)*21 + 1)/2))*this%n1

            ! Lateral positions of the fault at all depths
            bs = zeros(n1)

            do fi = 1, nf

                ww = w
                if (this%yn_facies) then
                    mm = m
                end if
                if (this%yn_rgt) then
                    tt = t
                end if
                ff = f
                ffdip = fdip
                ffdisp = fdisp
                cf = 0

                b_prev = 0
                fblock = falses(n1, n2)

                tp = (tan(dip(fi) - const_pi_half))*n1*0.5

                do i = 1, n1

                    theta = dips(i, fi)

                    if (i == 1) then
                        b = f2(fi) - tp
                        dxys = 0.0
                    else
                        dxys = -1.0/tan(dips(i - 1, fi))
                        b = b_prev + dxys
                    end if
                    bs(i) = b

                    ! Local displacement scaling from the slip patch; in 2D it
                    ! only depends on depth. Where the local displacement
                    ! diminishes below half a grid point, the fault causes no
                    ! visible offset, and the depth is beyond the fault tip
                    ! and therefore not marked as fault
                    if (this%yn_vary_disp) then
                        alpha = max(0.0, 1.0 - ((i - disp_zc(fi))/disp_az(fi))**2)
                        yn_mark = abs(disp(fi))*alpha >= 0.5
                    else
                        alpha = 1.0
                        yn_mark = .true.
                    end if

                    if (yn_mark) then

                        !$omp parallel do private(j, x0, dist)
                        do j = 1, n2

                            x0 = j - 1.0

                            dist = abs(x0 - b)
                            if (dist < 0.5*fwidth/sin(theta)) then
                                f(i, j) = fi
                                cf(i, j) = 1.0
                                fdip(i, j) = theta
                                ! Store the local displacement magnitude
                                fdisp(i, j) = abs(disp(fi))*alpha
                            end if

                        end do
                        !$omp end parallel do

                        if (abs(dxys) >= 1.0) then

                            !$omp parallel do private(j, x0, dist_prev, dist)
                            do j = 1, n2

                                x0 = j - 1.0

                                dist_prev = abs(x0 - b_prev)
                                dist = abs(x0 - b)
                                if (dist_prev <= abs(dxys) .and. dist <= abs(dxys)) then
                                    f(i, j) = fi
                                    cf(i, j) = 1.0
                                    fdip(i, j) = theta
                                    fdisp(i, j) = abs(disp(fi))*alpha
                                end if
                            end do
                            !$omp end parallel do

                        end if

                    end if

                    !$omp parallel do private(j, x0)
                    do j = 1, n2
                        x0 = j - 1.0
                        if (theta < const_pi_half) then
                            if (x0 - b <= 0) then
                                fblock(i, j) = .true.
                            end if
                        else
                            if (x0 - b >= 0) then
                                fblock(i, j) = .true.
                            end if
                        end if
                    end do
                    !$omp end parallel do

                    b_prev = b

                end do

                ! Shift the fault block; here for each target point within the fault
                ! block, the values are gathered from its source point, which is
                ! equivalent to the previous shift-source-to-target (scatter) approach
                ! for a constant fault displacement, while it naturally supports
                ! spatially varying (e.g., depth-varying) fault displacement without
                ! leaving holes in the model
                !$omp parallel do private(i, j, x0, alpha, dloc, newi, newj, ia, ja, wa, wb) collapse(2)
                do j = 1, n2
                    do i = 1, n1

                        if (fblock(i, j)) then

                            ! Local displacement from the slip patch, maximum at
                            ! the patch center and diminishing to zero at the
                            ! patch boundary (fault tips)
                            if (this%yn_vary_disp) then
                                alpha = max(0.0, 1.0 - ((i - disp_zc(fi))/disp_az(fi))**2)
                            else
                                alpha = 1.0
                            end if

                            ! Decay the displacement away from the fault with a Gaussian
                            ! profile of the distance to the fault, mimicking the deformation
                            ! halo of a finite slip patch (fault drag, rollover, and
                            ! blind-fault folding); the displacement at the fault
                            ! itself is unaffected
                            if (this%yn_disp_decay) then
                                x0 = j - 1.0
                                alpha = alpha*exp(-0.5*((x0 - bs(i))/(decay_w(fi) + float_tiny))**2)
                            end if

                            dloc = disp(fi)*alpha

                            ! The source point of the target point
                            newi = nint(i - dloc*sin(dip(fi)))
                            newj = nint(j + dloc*cos(dip(fi)))

                            if (newi >= 1 .and. newi <= n1 .and. newj >= 1 .and. newj <= n2) then

                                ! Move velocity model
                                w(i, j) = ww(newi, newj)

                                ! Move facies
                                if (this%yn_facies) then
                                    m(i, j) = mm(newi, newj)
                                end if

                                ! Move RGT
                                if (this%yn_rgt) then
                                    t(i, j) = tt(newi, newj)
                                end if

                                ! Move existing faults
                                if (cf(i, j) == 0) then
                                    f(i, j) = ff(newi, newj)
                                    fdip(i, j) = ffdip(newi, newj)
                                    fdisp(i, j) = ffdisp(newi, newj)
                                end if

                                ! For a spatially varying displacement, rounding the source
                                ! point to the nearest grid point quantizes the warp into
                                ! integer shifts, which produces terrace (staircase) artifacts
                                ! along the iso-displacement contours; therefore, the continuous
                                ! fields (velocity and RGT) are gathered with bilinear
                                ! interpolation instead, while the categorical fields (facies
                                ! and fault attributes) keep the nearest-neighbor sampling
                                if (this%yn_vary_disp .or. this%yn_disp_decay) then

                                    ia = clip(floor(i - dloc*sin(dip(fi))), 1, n1 - 1)
                                    ja = clip(floor(j + dloc*cos(dip(fi))), 1, n2 - 1)
                                    wa = clip(i - dloc*sin(dip(fi)) - ia, 0.0, 1.0)
                                    wb = clip(j + dloc*cos(dip(fi)) - ja, 0.0, 1.0)

                                    w(i, j) = &
                                        ww(ia, ja)*(1 - wa)*(1 - wb) &
                                        + ww(ia + 1, ja)*wa*(1 - wb) &
                                        + ww(ia, ja + 1)*(1 - wa)*wb &
                                        + ww(ia + 1, ja + 1)*wa*wb

                                    if (this%yn_rgt) then
                                        t(i, j) = &
                                            tt(ia, ja)*(1 - wa)*(1 - wb) &
                                            + tt(ia + 1, ja)*wa*(1 - wb) &
                                            + tt(ia, ja + 1)*(1 - wa)*wb &
                                            + tt(ia + 1, ja + 1)*wa*wb
                                    end if

                                end if

                            end if

                        end if

                    end do
                end do
                !$omp end parallel do

            end do

        end if

        ! Select model before adding salt bodies
        ! Note that n1 = this%n1 + 1 for computing reflectivity coefficients
        vp = w(ne1 + 1:ne1 + this%n1 + 1, ne2 + 1:ne2 + this%n2)
        vp = rescale(vp, [this%vmin, this%vmax])
        rho = this%rho_a*vp**this%rho_b + this%rho_c
        if (this%yn_elastic) then
            vs = vp/rescale(vp, this%vpvsratio)
        end if

        n1 = size(vp, 1)
        n2 = size(vp, 2)

        ! Add salt
        if (this%yn_salt) then

            if (maxval(this%salt_radius) == 0) then
                salt_radius = [0.05*n2, 0.2*n2]
            else
                salt_radius = this%salt_radius
            end if
            tp = mean(salt_radius)

            select case (this%refl_shape)
                case ('random', 'perlin', 'custom')
                    gmax = random(this%nsalt, range=[tp, n2 - tp], seed=safe_seed(int(this%seed, 8)*5), spacing=0.5*tp)
                case ('gaussian', 'cauchy')
                    if (this%nsalt > this%ng) then
                        gmax = [mu, random(this%nsalt - this%ng, range=[tp, n2 - tp], seed=safe_seed(int(this%seed, 8)*5), spacing=0.5*tp)]
                    else
                        gmax = mu
                    end if
            end select

            this%salt = zeros(n1, n2)
            rds = random(this%nsalt, seed=safe_seed(int(this%seed, 8)*14 - 1))
            rds = rescale(rds, salt_radius)
            vds = random(2*this%nsalt, seed=safe_seed(int(this%seed, 8)*15 - 1))
            vds = rescale(vds, [0.75*this%salt_radius_variation, this%salt_radius_variation])
            pds = rescale(rds, [0.75*this%salt_path_variation, this%salt_path_variation])

            ! Build top interface of the salt
            pn%n1 = n2
            pn%octaves = 4
            pn%seed = safe_seed(int(this%seed, 8)*14 - 2)
            topz = pn%generate()
            topz = rescale(topz, [0.0, this%salt_top_height])

            ! Insert salt bodies
            do isalt = 1, this%nsalt

                nd = nint((1.0 - rand(range=this%salt_top_z, seed=safe_seed(int(this%seed, 8)*15*isalt - 2)))*this%n1)

                ! Salt body boundaries
                pn%n1 = this%salt_nnode
                pn%octaves = 4
                pn%seed = safe_seed(int(this%seed, 8)*15*isalt - 1)
                x1 = pn%generate()

                pn%n1 = this%salt_nnode
                pn%octaves = 4
                pn%seed = safe_seed(int(this%seed, 8)*16*isalt - 1)
                x2 = pn%generate()

                x1 = interp_to(x1, n1, 'pchip')
                x2 = interp_to(x2, n1, 'pchip')

                x1 = rescale(x1, range=[1.0 - vds(isalt), 1.0]*rds(isalt))
                x2 = rescale(x2, range=[1.0 - vds(isalt), 1.0]*rds(isalt))

                ! Salt body path deviations
                pn%n1 = n1
                pn%octaves = 4
                pn%seed = safe_seed(int(this%seed, 8)*17*isalt - 1)
                xs = pn%generate()
                xs = xs - mean(xs)
                xs = median_filt(xs, 2)/maxval(xs)*pds(isalt)

                !$omp parallel do private(i, j)
                do j = 1, n2
                    do i = max(n1 - nd - ceiling(maxval(topz)), 1), n1
                        if (j >= gmax(isalt) - x1(i) + xs(i) .and. j <= gmax(isalt) + x2(i) + xs(i) &
                                .and. i >= n1 - nd - topz(j)) then
                            vp(i, j) = this%salt_vp
                            rho(i, j) = this%salt_rho
                            if (this%yn_elastic) then
                                vs(i, j) = this%salt_vs
                            end if
                            this%salt(i, j) = 1.0
                        end if
                    end do
                end do
                !$omp end parallel do

            end do

        end if

        ! Add karst
        if (this%yn_karst) then

            block

                type(karst_2d) :: kg
                real, allocatable, dimension(:, :) :: kvol
                integer :: nzk

                ! The karst passages are confined to the depth window karst_z;
                ! generate on a sub-grid spanning [0, karst_z(2)]*n1 with the
                ! caprock fraction covering [0, karst_z(1)]*n1
                nzk = clip(nint(this%karst_z(2)*n1), 2, n1)
                kg%nx = n2
                kg%nz = nzk
                kg%n_passages = this%karst_npassage
                kg%n_ctrl = this%karst_nctrl
                if (maxval(this%karst_radius) == 0) then
                    kg%r_mean = rand(range=[0.015, 0.03]*this%n1, seed=safe_seed(int(this%seed, 8)*27 - 1))
                else
                    kg%r_mean = rand(range=this%karst_radius, seed=safe_seed(int(this%seed, 8)*27 - 1))
                end if
                kg%r_std = this%karst_radius_variation
                kg%tortuosity = this%karst_tortuosity
                kg%p_connect = this%karst_connect
                kg%z_top = 1.0 - this%karst_z(1)/this%karst_z(2)
                kg%seed = safe_seed(int(this%seed, 8)*29 - 1)
                kvol = kg%generate()

                this%karst = zeros(n1, n2)
                !$omp parallel do private(i, j) collapse(2)
                do j = 1, n2
                    do i = 1, nzk
                        if (kvol(i, j) > 0.5) then
                            vp(i, j) = this%karst_vp
                            rho(i, j) = this%karst_rho
                            if (this%yn_elastic) then
                                vs(i, j) = this%karst_vs
                            end if
                            this%karst(i, j) = 1.0
                            ! Karst carves salt where they overlap
                            if (this%yn_salt) then
                                this%salt(i, j) = 0.0
                            end if
                        end if
                    end do
                end do
                !$omp end parallel do

            end block

        end if

        ! Generate images
        if (this%unconf == 0) then
            if (this%yn_elastic) then
                call this%generate_image_elastic(vp, vs, rho)
            else
                call this%generate_image(vp, rho)
            end if
        end if

        ! Output
        ! Fault and fault attributes
        if (this%yn_fault) then
            this%fault = f(ne1 + 1:ne1 + this%n1, ne2 + 1:ne2 + this%n2)
            this%fault_dip = fdip(ne1 + 1:ne1 + this%n1, ne2 + 1:ne2 + this%n2)*const_rad2deg
            this%fault_disp = fdisp(ne1 + 1:ne1 + this%n1, ne2 + 1:ne2 + this%n2)
        end if

        ! RGT
        if (this%yn_rgt) then
            this%rgt = t(ne1 + 1:ne1 + this%n1, ne2 + 1:ne2 + this%n2)
            this%rgt = rescale(this%rgt, [0.0, 1.0])
        end if

        ! Facies
        if (this%yn_facies) then
            this%facies = m(ne1 + 1:ne1 + this%n1, ne2 + 1:ne2 + this%n2)
            this%facies = this%facies - minval(this%facies) + 1
            if (this%unconf == 0) then
                this%facies = maxval(this%facies) - this%facies + 1
            end if
        end if

        ! Salt
        if (this%yn_salt) then

            this%salt = this%salt(1:this%n1, :)

            if (this%unconf == 0) then

                if (this%yn_fault) then
                    where (this%salt == 1)
                        this%fault = 0
                        this%fault_dip = 0
                        this%fault_disp = 0
                    end where
                end if

                if (this%yn_rgt) then
                    where (this%salt == 1)
                        this%rgt = 0
                    end where
                end if

                if (this%yn_facies) then
                    where (this%salt == 1)
                        this%facies = 0
                    end where
                end if

            end if

        end if

        ! Karst
        if (this%yn_karst) then

            this%karst = this%karst(1:this%n1, :)

            if (this%unconf == 0) then

                if (this%yn_fault) then
                    where (this%karst == 1)
                        this%fault = 0
                        this%fault_dip = 0
                        this%fault_disp = 0
                    end where
                end if

                if (this%yn_rgt) then
                    where (this%karst == 1)
                        this%rgt = 0
                    end where
                end if

                if (this%yn_facies) then
                    where (this%karst == 1)
                        this%facies = 0
                    end where
                end if

            end if

        end if

        ! Velocity models
        this%vp = vp(1:this%n1, :)
        this%rho = rho(1:this%n1, :)
        if (this%yn_elastic) then
            this%vs = vs(1:this%n1, :)
        end if

        ! Unconformity
        if (this%unconf > 0 .and. this%unconf_nl == 0) then
            if (this%yn_elastic) then
                this%vp = minval(this%vp)
                this%vs = minval(this%vs)
                this%rho = minval(this%rho)
            else
                this%vp = minval(this%vp)
                this%rho = minval(this%rho)
            end if
            if (this%yn_rgt) then
                this%rgt = 0
            end if
            if (this%yn_fault) then
                this%fault = 0
                this%fault_dip = 0
                this%fault_disp = 0
            end if
            if (this%yn_facies) then
                this%facies = 1
            end if
            if (this%yn_salt) then
                this%salt = 0
            end if
            if (this%yn_karst) then
                this%karst = 0
            end if
        end if

    end subroutine generate_2d_geological_model

    !
    !> Generate geological models with one or multiple unconformity surfaces
    !
    subroutine generate_2d_unconformal_geological_model(this)

        type(rgm2_curved), intent(inout) :: this

        type(rgm2_curved), allocatable, dimension(:) :: g
        integer :: iconf, i, j
        type(meta_array1_real), allocatable, dimension(:) :: uff
        real, allocatable, dimension(:) :: ufz
        real, allocatable, dimension(:, :) :: rgt_above, rgt_below
        real, allocatable, dimension(:, :) :: facies_above, facies_below
        real :: tmin, tmax
        type(fractal_noise_1d) :: q
        real, allocatable, dimension(:, :) :: vp, vs, rho

        allocate (g(1:this%unconf + 1))

        if (this%yn_elastic) then
            this%vp = zeros(this%n1 + 1, this%n2)
            this%vs = zeros(this%n1 + 1, this%n2)
            this%rho = zeros(this%n1 + 1, this%n2)
        else
            this%vp = zeros(this%n1 + 1, this%n2)
            this%rho = zeros(this%n1 + 1, this%n2)
        end if

        if (this%yn_fault) then
            this%fault = zeros(this%n1 + 1, this%n2)
            this%fault_dip = zeros(this%n1 + 1, this%n2)
            this%fault_disp = zeros(this%n1 + 1, this%n2)
        end if

        if (this%yn_rgt) then
            this%rgt = zeros(this%n1 + 1, this%n2)
        end if

        if (this%yn_facies) then
            this%facies = zeros(this%n1 + 1, this%n2)
        end if

        if (this%yn_salt) then
            this%salt = zeros(this%n1 + 1, this%n2)
        end if
        if (this%yn_karst) then
            this%karst = zeros(this%n1 + 1, this%n2)
        end if

        do i = 1, this%unconf + 1

            g(i) = this
            g(i)%n1 = this%n1 + 1

            g(i)%yn_salt = .false.
            if (this%yn_salt .and. i == this%unconf + 1) then
                g(i)%yn_salt = .true.
            end if

            g(i)%yn_karst = .false.
            if (this%yn_karst .and. i == this%unconf + 1) then
                g(i)%yn_karst = .true.
            end if

            if (this%seed /= -1) then
                g(i)%seed = safe_seed(int(g(i)%seed, 8)*i)
            end if

            if (i <= this%unconf) then
                g(i)%nf = 0
                g(i)%lwv = abs(g(i)%lwv/2.0)
                g(i)%lwh = abs(g(i)%lwh/2.0)
                g(i)%refl_height = g(i)%unconf_refl_height
                g(i)%refl_slope = g(i)%unconf_refl_slope
                g(i)%refl_height_top = g(i)%unconf_refl_height
                g(i)%refl_slope_top = g(i)%unconf_refl_slope
                g(i)%refl_shape = g(i)%unconf_refl_shape
                g(i)%refl_shape_top = g(i)%unconf_refl_shape
                g(i)%refl_smooth = g(i)%unconf_refl_smooth
                g(i)%refl_smooth_top = g(i)%unconf_refl_smooth
            end if

            if (i < this%unconf + 1 .and. this%unconf_nl == 0) then
                g(i)%unconf_nl = 0
            else
                g(i)%unconf_nl = g(i)%nl
            end if

            call generate_2d_geological_model(g(i))

        end do

        allocate (uff(1:this%unconf + 1))
        ufz = random(this%unconf, range=this%unconf_z, seed=safe_seed(int(this%seed, 8)*31))
        ufz = sort(ufz, order=1)
        do i = 1, this%unconf

            select case (this%unconf_shape)

                case ('meander_channel', 'meander_canyon', 'drainage_channel', 'drainage_canyon')
                    ! The unconformity surface is a cross-flow section of the
                    ! erosional topography of a geomorphological simulation
                    call unconformity_topography_2d(this, g(i)%seed, i, uff(i)%array)

                case default
                    ! Use Perlin noise to generate unconformity surfaces
                    q%n1 = this%n2
                    q%octaves = 5
                    q%seed = safe_seed(int(g(i)%seed, 8)*41*i)
                    uff(i)%array = q%generate()

            end select

            if (this%unconf_smooth > 0) then
                uff(i)%array = gauss_filt(uff(i)%array, this%unconf_smooth)
            end if

            uff(i)%array = rescale(uff(i)%array, &
                [0.0, rand(range=this%unconf_height, seed=safe_seed(int(g(i)%seed, 8)*51*i))]) + ufz(i)*this%n1

        end do

        ! Merge sedimentary units
        if (this%yn_elastic) then
            this%vp = g(this%unconf + 1)%vp
            this%vs = g(this%unconf + 1)%vs
            this%rho = g(this%unconf + 1)%rho
        else
            this%vp = g(this%unconf + 1)%vp
            this%rho = g(this%unconf + 1)%rho
        end if
        if (this%yn_fault) then
            this%fault = g(this%unconf + 1)%fault
            this%fault_dip = g(this%unconf + 1)%fault_dip
            this%fault_disp = g(this%unconf + 1)%fault_disp
        end if
        if (this%yn_rgt) then
            this%rgt = g(this%unconf + 1)%rgt
            rgt_above = zeros(this%n1, this%n2)
            rgt_below = zeros(this%n1, this%n2)
        end if
        if (this%yn_facies) then
            this%facies = g(this%unconf + 1)%facies
            facies_above = zeros(this%n1, this%n2)
            facies_below = zeros(this%n1, this%n2)
        end if
        if (this%yn_salt) then
            this%salt = g(this%unconf + 1)%salt
        end if
        if (this%yn_karst) then
            this%karst = g(this%unconf + 1)%karst
        end if

        do iconf = this%unconf, 1, -1

            if (this%yn_rgt) then
                rgt_above = 0
                rgt_below = 0
            end if

            if (this%yn_facies) then
                facies_above = 0
                facies_below = 0
            end if

            !$omp parallel do private(i, j)
            do j = 1, this%n2

                ! Image by soft merging
                if (this%yn_elastic) then
                    do i = 1, this%n1
                        if (i < uff(iconf)%array(j)) then
                            this%vp(i, j) = g(iconf)%vp(i, j)
                            this%vs(i, j) = g(iconf)%vs(i, j)
                            this%rho(i, j) = g(iconf)%rho(i, j)
                            if (this%yn_salt .and. this%salt_before_unconf) then
                                this%salt(i, j) = 0.0
                            end if
                            if (this%yn_karst .and. this%karst_before_unconf) then
                                this%karst(i, j) = 0.0
                            end if
                        end if
                    end do
                else
                    do i = 1, this%n1
                        if (i < uff(iconf)%array(j)) then
                            this%vp(i, j) = g(iconf)%vp(i, j)
                            this%rho(i, j) = g(iconf)%rho(i, j)
                            if (this%yn_salt .and. this%salt_before_unconf) then
                                this%salt(i, j) = 0.0
                            end if
                            if (this%yn_karst .and. this%karst_before_unconf) then
                                this%karst(i, j) = 0.0
                            end if
                        end if
                    end do
                end if

                ! Fault by hard merging
                if (this%yn_fault) then
                    do i = 1, this%n1
                        if (i < uff(iconf)%array(j)) then
                            this%fault(i, j) = g(iconf)%fault(i, j)
                            this%fault_dip(i, j) = g(iconf)%fault_dip(i, j)
                            this%fault_disp(i, j) = g(iconf)%fault_disp(i, j)
                        end if
                    end do
                end if

                ! RGT by history-consistent merging
                if (this%yn_rgt) then
                    do i = 1, this%n1
                        if (i < uff(iconf)%array(j)) then
                            rgt_above(i, j) = g(iconf)%rgt(i, j)
                        else
                            rgt_below(i, j) = this%rgt(i, j)
                        end if
                    end do
                end if

                ! Facies by hard merging
                if (this%yn_facies) then
                    do i = 1, this%n1
                        if (i < uff(iconf)%array(j)) then
                            facies_above(i, j) = g(iconf)%facies(i, j)
                        else
                            facies_below(i, j) = this%facies(i, j)
                        end if
                    end do
                end if

            end do
            !$omp end parallel do

            ! Correct RGT
            if (this%yn_rgt) then
                tmax = maxval(rgt_above)
                tmin = minval(rgt_below, mask=(rgt_below /= 0))
                !$omp parallel do private(i, j) collapse(2)
                do j = 1, this%n2
                    do i = 1, this%n1
                        if (rgt_above(i, j) /= 0) then
                            rgt_above(i, j) = rgt_above(i, j) - tmax + tmin
                        end if
                        this%rgt(i, j) = rgt_above(i, j) + rgt_below(i, j)
                    end do
                end do
                !$omp end parallel do
            end if

            ! Correct facies
            if (this%yn_facies) then
                tmin = minval(facies_above, mask=(facies_above /= 0))
                tmax = maxval(facies_below)
                !$omp parallel do private(i, j) collapse(2)
                do j = 1, this%n2
                    do i = 1, this%n1
                        if (facies_above(i, j) /= 0) then
                            facies_above(i, j) = facies_above(i, j) - tmin + tmax + 1
                        end if
                        this%facies(i, j) = facies_above(i, j) + facies_below(i, j)
                    end do
                end do
                !$omp end parallel do
            end if

        end do

        if (this%yn_salt) then
            where (this%salt == 1)
                this%vp = this%salt_vp
                this%rho = this%salt_rho
            end where
            if (this%yn_elastic) then
                where (this%salt == 1)
                    this%vs = this%salt_vs
                end where
            end if
        end if

        if (this%yn_karst) then
            where (this%karst == 1)
                this%vp = this%karst_vp
                this%rho = this%karst_rho
            end where
            if (this%yn_elastic) then
                where (this%karst == 1)
                    this%vs = this%karst_vs
                end where
            end if
        end if

        vp = this%vp
        rho = this%rho
        if (this%yn_elastic) then
            vs = this%vs
        end if

        ! Final processing
        this%vp = this%vp(1:this%n1, :)
        this%rho = this%rho(1:this%n1, :)
        if (this%yn_elastic) then
            this%vs = this%vs(1:this%n1, :)
        end if

        if (this%yn_fault) then
            this%fault = this%fault(1:this%n1, :)
            this%fault_dip = this%fault_dip(1:this%n1, :)
            this%fault_disp = this%fault_disp(1:this%n1, :)
        end if

        if (this%yn_rgt) then
            this%rgt = this%rgt(1:this%n1, :)
        end if

        if (this%yn_facies) then
            this%facies = this%facies(1:this%n1, :)
        end if

        if (this%yn_salt) then
            this%salt = this%salt(1:this%n1, :)
        end if

        if (this%yn_karst) then
            this%karst = this%karst(1:this%n1, :)
        end if

        ! Rescale RGT to [0, 1]
        if (this%yn_rgt) then

            this%rgt = rescale(this%rgt, [0.0, 1.0])

            ! When the unconformity represents seafloor, reprocess the RGT
            ! so that the RGT of the seawater is zero
            if (this%unconf_nl == 0) then
                tmin = minval(this%rgt, mask=(this%rgt > 0))
                where (this%rgt > 0)
                    this%rgt = this%rgt - tmin
                end where
                this%rgt = rescale(this%rgt, [0.0, 1.0])
            end if

        end if

        ! Process facies
        if (this%yn_facies) then
            this%facies = maxval(this%facies) - this%facies + 1
        end if

        ! Process salt
        if (this%yn_salt) then

            if (this%yn_fault) then
                where (this%salt == 1)
                    this%fault = 0
                    this%fault_dip = 0
                    this%fault_disp = 0
                end where
            end if

            if (this%yn_rgt) then
                where (this%salt == 1)
                    this%rgt = 0
                end where
            end if

            if (this%yn_facies) then
                where (this%salt == 1)
                    this%facies = 0
                end where
            end if

        end if

        ! Process karst
        if (this%yn_karst) then

            if (this%yn_fault) then
                where (this%karst == 1)
                    this%fault = 0
                    this%fault_dip = 0
                    this%fault_disp = 0
                end where
            end if

            if (this%yn_rgt) then
                where (this%karst == 1)
                    this%rgt = 0
                end where
            end if

            if (this%yn_facies) then
                where (this%karst == 1)
                    this%facies = 0
                end where
            end if

        end if

        ! Finally, generate image
        if (this%yn_elastic) then
            call this%generate_image_elastic(vp, vs, rho)
        else
            call this%generate_image(vp, rho)
        end if

        if (.not. this%custom_psf .and. this%wave /= '') then
            this%psf = g(1)%psf
        end if

    end subroutine generate_2d_unconformal_geological_model

    !
    !> Generate the topography of a channel/canyon/drainage-type unconformity
    !> curve for a 2D model. A map-view simulation is run on an auxiliary grid
    !> with the flow forced along the auxiliary x-axis, and a random cross-flow
    !> column of the erosional depth map, which typically crosses the channel
    !> belt multiple times, is used as the 1D unconformity topography.
    !
    subroutine unconformity_topography_2d(this, seed, isurf, surf)

        type(rgm2_curved), intent(in) :: this
        integer, intent(in) :: seed, isurf
        real, allocatable, dimension(:), intent(out) :: surf

        type(meandering_channel) :: mc
        type(meandering_canyon) :: my
        type(drainage_channel) :: dc
        type(drainage_canyon) :: dy
        real, allocatable, dimension(:, :) :: dm
        real :: wf
        integer :: nlev, naux, ix

        ! Reference (calibrated) configurations of the meandering channel and
        ! canyon: at length = *_length_ref with *_nbends_ref seed bends, the
        ! migration develops well-formed gooseneck meanders. The mappings
        ! below scale unconf_channel_length relative to these single
        ! calibration points
        real, parameter :: mc_length_ref = 25000.0
        real, parameter :: mc_nbends_ref = 30.0
        real, parameter :: my_length_ref = 5000.0
        real, parameter :: my_nbends_ref = 20.0
        ! Minimum accepted channel length; below this the centerline has too
        ! few nodes for the migration to develop realistic meanders
        real, parameter :: length_min = 10000.0

        wf = rand(range=this%unconf_channel_width, seed=safe_seed(int(seed, 8)*91 + isurf))

        ! Auxiliary along-flow grid size
        naux = max(this%n2, 128)

        select case (this%unconf_shape)

            case ('meander_channel')
                mc%nx = naux
                mc%ny = this%n2
                mc%nz = 64
                ! The length/W ratio sets how many meander bends span the
                ! map. The migration width W is capped at its calibrated
                ! value so that longer channels develop more bends of the
                ! same shape (the migration dynamics are width-invariant),
                ! while shorter channels shrink W proportionally to keep
                ! >= ~80 centerline nodes; the seed sine's wavelength
                ! (length/n_bends) is held at its calibrated value; the
                ! channel is rendered with W_render so that the on-grid width
                ! is approximately the requested fraction of the lateral extent
                mc%length = max(this%unconf_channel_length, length_min)
                mc%W = 0.025*min(mc%length, mc_length_ref)
                mc%W_render = wf*mc%length
                mc%W_shape = 3
                mc%n_bends = max(5, nint(mc_nbends_ref*mc%length/mc_length_ref))
                mc%n_iter = nint(1500*this%unconf_channel_sinuosity)
                mc%terrain_bg = this%unconf_topo*(mc%nz - 1.0)
                mc%seed = safe_seed(int(seed, 8)*61 + isurf)
                mc%orient = 0
                mc%yn_depth_only = .true.
                call mc%generate
                dm = mc%depth_map

            case ('meander_canyon')
                my%nx = naux
                my%ny = this%n2
                my%nz = 0
                ! As for the channel, W is capped at its calibrated value so
                ! that a longer canyon adds bends rather than changing the
                ! migration dynamics
                my%length = (my_length_ref/mc_length_ref)*max(this%unconf_channel_length, length_min)
                my%W = 0.02*min(my%length, my_length_ref)
                my%W_canyon = 0.5*wf*my%length
                my%n_bends = max(5, nint(my_nbends_ref*my%length/my_length_ref))
                my%n_iter = nint(clip(4000*this%unconf_channel_sinuosity, 2000.0, 6000.0))
                nlev = max(2, nint(my%n_iter/(my%save_every*1.0)))
                my%terrain_bg = this%unconf_topo*(nlev - 1.0)
                my%seed = safe_seed(int(seed, 8)*61 + isurf)
                my%orient = 0
                my%yn_depth_only = .true.
                call my%generate
                dm = my%depth_map

            case ('drainage_channel')
                dc%nx = naux
                dc%ny = this%n2
                dc%nz = 64
                ! Channels: relatively thin and densely distributed
                dc%W_max = 0.5*wf*this%n2
                dc%D_max = 48.0
                dc%channel_frac = rand(range=this%unconf_channel_density, seed=safe_seed(int(seed, 8)*81 + isurf))
                dc%terrain_bg = this%unconf_topo*(dc%nz - 1.0)
                dc%seed = safe_seed(int(seed, 8)*61 + isurf)
                dc%orient = 2
                dc%yn_depth_only = .true.
                call dc%generate
                dm = dc%depth_map

            case ('drainage_canyon')
                dy%nx = naux
                dy%ny = this%n2
                dy%nz = 64
                ! Canyons: major and few - full width but sparse network
                dy%W_max = wf*this%n2
                dy%channel_frac = 0.25*rand(range=this%unconf_channel_density, seed=safe_seed(int(seed, 8)*81 + isurf))
                dy%terrain_bg = this%unconf_topo*(dy%nz - 1.0)
                dy%seed = safe_seed(int(seed, 8)*61 + isurf)
                dy%orient = 2
                dy%yn_depth_only = .true.
                call dy%generate
                dm = dy%depth_map

        end select

        ! Take a random cross-flow column as the unconformity topography
        ix = clip(nint(rand(range=[0.1, 0.9], seed=safe_seed(int(seed, 8)*101 + isurf))*naux), 1, naux)
        surf = dm(:, ix)

    end subroutine unconformity_topography_2d

end module geological_model_2d_curved

