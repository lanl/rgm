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

module geological_model_3d_curved

    use libflit
    use geological_model_utility
    use geological_model_meander
    use geological_model_drainage
    use geological_model_karst

    implicit none

    !
    ! 3D random geological model generator (type rgm3_curved).
    !
    ! The generator follows a geologically motivated workflow:
    ! build layered velocity/density (and optionally elastic) models between two
    ! bounding surfaces -> deform the model with faults -> optionally insert
    ! salt bodies and unconformities -> compute reflectivity -> convolve with a
    ! wavelet to produce a synthetic migration-like image, together with
    ! voxel-wise labels (fault index/dip/strike/rake/displacement, relative
    ! geological time, facies, salt).
    !
    ! Faults are curved: the dip can vary with depth (listric faults via
    ! delta_dip), and the strike can vary along the fault (delta_strike).
    !
    ! The strike of a fault can spatially vary along the fault when
    ! delta_strike > 0. In this case, the fault surface is no longer planar
    ! in the horizontal directions; the local strike angle stored in fault_strike
    ! varies along the fault, and the local dip angle stored in fault_dip
    ! is corrected for the strike deviation.
    !
    ! The displacement of a fault can spatially vary on the fault surface when
    ! yn_vary_disp = .true. In this case, the displacement follows an elliptical
    ! slip patch that is maximum at the patch center and diminishes to zero at
    ! the patch boundary (fault tip line), so a fault can die out within the
    ! model; fault_disp stores the spatially varying displacements.
    !
    ! Additionally, the displacement can decay away from the fault along the
    ! fault normal when yn_disp_decay = .true., mimicking the deformation halo
    ! of a finite slip patch (fault drag, rollover, and blind-fault folding);
    ! and for strike-varying faults, the slip direction follows the local
    ! strike of the fault.
    !

    ! 3D random geological model with curved faults
    ! Meanings of the parameters are similar with those in the 2D case
    type rgm3_curved

        !==============================================================================================
        integer :: n1 = 128
        integer :: n2 = 128
        integer :: n3 = 128
        integer :: nf = 4
        integer :: nl = 20
        integer :: seed = -1
        real :: refl_smooth = 20.0
        real :: refl_smooth_top = 40.0
        real :: dt = 1.0e-3
        real :: f0 = 150.0
        real :: fwidth = 2.0
        real, allocatable, dimension(:) :: dip, strike, rake, disp
        real, dimension(1:2) :: refl_slope = [0.0, 0.0]
        real, dimension(1:2) :: refl_slope_top = [0.0, 0.0]
        real, dimension(1:3) :: noise_smooth = [1.0, 1.0, 1.0]
        real :: noise_level = 0.05
        character(len=24) :: wave = 'ricker'
        character(len=24) :: refl_shape = 'random'
        character(len=24) :: refl_shape_top = 'random'
        real, allocatable, dimension(:, :) :: refl
        real, allocatable, dimension(:, :) :: refl_top
        integer :: ng = 2
        real, dimension(1:2) :: refl_height = [0.0, 10.0]
        real, dimension(1:2) :: refl_height_top = [0.0, 5.0]
        !> Range of Gaussian standard devision along x2 for refl_shape = gaussian
        real, dimension(1:2) :: refl_sigma2 = [0.0, 0.0]
        !> Range of Gaussian mean along x2 for refl_shape = gaussian
        real, dimension(1:2) :: refl_mu2 = [0.0, 0.0]
        !> Range of Gaussian standard devision along x3 for refl_shape = gaussian
        real, dimension(1:2) :: refl_sigma3 = [0.0, 0.0]
        !> Range of Gaussian mean along x3 for refl_shape = gaussian
        real, dimension(1:2) :: refl_mu3 = [0.0, 0.0]
        real :: lwv = 0.25
        !> Horizontal layer thickness variation
        real :: lwh = 0.1
        !> Secondary reflector smoothing
        real :: secondary_refl_smooth = 10.0
        !> For Gaussian, Cauchy surface, whether to rotate
        logical :: rotate_fold = .false.

        logical :: yn_rgt = .false.
        logical :: yn_facies = .false.
        logical :: yn_fault = .true.
        real, allocatable, dimension(:, :, :) :: image, rgt, facies, fault
        real, allocatable, dimension(:, :, :) :: fault_dip, fault_strike, fault_rake, fault_disp
        real, dimension(1:3) :: psf_sigma = [5.0, 2.5, 2.5]
        real, allocatable, dimension(:, :, :) :: psf
        logical :: custom_psf = .false.
        real :: facies_threshold = 0.0
        character(len=12) :: noise_type = 'normal'
        logical :: yn_conv_noise = .false.
        logical :: yn_regular_fault = .false.
        logical :: yn_group_faults = .false.
        real, allocatable, dimension(:) :: wave_filt_freqs, wave_filt_amps
        !> Min value for scaling the facies
        real :: vmin = 2000.0
        !> Max value for scaling the facies; after scaling, the facies will fall in [vmin, vmax]
        real :: vmax = 4000.0
        !> Velocity perturbation of layers
        real :: delta_v = 500.0

        !> Dip increase/descrease at the top compared with the top
        real, dimension(1:2) :: delta_dip = [15.0, 30.0]

        !> Range of the maximum strike deviation (in degrees) along a fault.
        !> When > 0, the strike of each fault spatially varies along the fault,
        !> so a fault is no longer a straight line in map view but a smooth curve;
        !> the local strike angles are stored in fault_strike, and fault_dip
        !> stores the local true dip angles that account for the strike deviation.
        !> When = [0, 0] (default), the faults have constant strikes.
        real, dimension(1:2) :: delta_strike = [0.0, 0.0]
        !> Number of undulation periods of the strike deviation along a fault;
        !> a larger value results in more oscillatory faults in map view.
        !> The strike variation is deliberately low-wavenumber (smooth, gentle bends)
        !> as high-wavenumber oscillation of fault strike is geologically uncommon.
        integer :: strike_nperiod = 2

        !> Whether to spatially vary the displacement on each fault surface.
        !> When = .true., the displacement follows an elliptical slip patch:
        !> maximum at the patch center and diminishing to zero at the patch
        !> boundary (fault tip line), so a fault can die out within the model
        !> rather than always cutting through the entire volume, and the
        !> hanging wall deforms gently near the fault tips.
        !> In this case, fault_disp stores the spatially varying displacements.
        logical :: yn_vary_disp = .false.
        !> Range of the along-strike semi-axis of the elliptical slip patch,
        !> relative to the average lateral model size 0.5*(n2 + n3)
        real, dimension(1:2) :: disp_radius_strike = [0.6, 1.2]
        !> Range of the along-depth semi-axis of the elliptical slip patch,
        !> relative to the model depth n1
        real, dimension(1:2) :: disp_radius_dip = [0.6, 1.2]
        !> Range of the center depth of the slip patch, relative to the model
        !> depth n1. Note that all the displacement tapering (and therefore the
        !> deformation of layers) occurs between the patch center and the patch
        !> boundary; setting a shallower center together with a smaller
        !> disp_radius_dip keeps the fault tips and the associated deformation
        !> away from the deep, high-velocity section
        real, dimension(1:2) :: disp_center_dip = [0.25, 0.75]
        !> Range of the along-strike center of the slip patch, relative to the
        !> average lateral model size 0.5*(n2 + n3) and to the fault center
        real, dimension(1:2) :: disp_center_strike = [-0.2, 0.2]

        !> Whether to decay the fault displacement away from the fault
        !> along the fault normal. When = .false. (default), the displaced block
        !> shifts rigidly (correct for through-going faults). When = .true.,
        !> the displacement of the shifted block decays with a Gaussian profile
        !> of the distance to the fault, mimicking the deformation halo of a
        !> finite slip patch (fault drag, rollover, blind-fault folding);
        !> the decay is applied only within the displaced block, and the
        !> displacement at the fault surface itself is unaffected.
        logical :: yn_disp_decay = .false.
        !> Range of the displacement decay width (the Gaussian standard deviation
        !> of the decay profile), relative to the model depth n1; smaller values
        !> concentrate the deformation near the fault (drag-fold style), while
        !> larger values approach the rigid-block end member.
        real, dimension(1:2) :: disp_decay_width = [0.5, 1.0]

        !==============================================================================================
        integer :: unconf = 0
        real, dimension(1:2) :: unconf_z = [0.0, 0.5]
        real, dimension(1:2) :: unconf_height = [5.0, 15.0]
        integer :: unconf_nl = 99999
        real :: unconf_smooth = 0.0
        real, dimension(1:2) :: unconf_refl_height = [0.0, 5.0]
        real :: unconf_refl_slope = -2.5
        real :: unconf_refl_smooth = 10.0
        character(len=12) :: unconf_refl_shape = 'random'
        !> Shape of the unconformity surfaces, can be one of
        !>      random - random (Perlin-noise) topography (default)
        !>      meander_channel - meandering river channels
        !>      meander_canyon - meandering incised canyon
        !>      drainage_channel - dendritic drainage channel network
        !>      drainage_canyon - canyon carved by a dendritic drainage network
        !> For the channel/canyon shapes, the topography is the erosional
        !> depth map of a geomorphological simulation, and unconf_height
        !> still bounds the total erosional relief
        character(len=24) :: unconf_shape = 'random'
        !> Channel/canyon width as an (approximate) fraction of the lateral
        !> model extent, drawn per unconformity surface
        real, dimension(1:2) :: unconf_channel_width = [0.03, 0.08]
        !> Meander maturity; scales the number of migration iterations
        real :: unconf_channel_sinuosity = 1.0
        !> Total centerline length of the meandering channel/canyon in the
        !> internal units of the migration simulation; since the simulated
        !> map is fit to the model grid, this effectively controls how many
        !> meander bends span the model laterally. Values larger than the
        !> default add more bends of the same shape; smaller values give
        !> fewer, larger bends (clamped from below to keep enough centerline
        !> nodes for the migration to develop realistic meanders)
        real :: unconf_channel_length = 25000.0
        !> Drainage density (fraction of cells covered by channels) for the
        !> drainage shapes, drawn per unconformity surface
        real, dimension(1:2) :: unconf_channel_density = [0.02, 0.08]
        !> Amplitude of the background (interfluve) topography relative to the
        !> channel/canyon incision depth; 0 = flat background between channels
        real :: unconf_topo = 0.25

        !==============================================================================================
        logical :: yn_salt = .false.
        integer :: nsalt = 1
        real, dimension(1:2) :: salt_radius = [0.0, 0.0]
        real :: salt_radius_variation = 0.7
        real :: salt_path_variation = 5.0
        integer :: salt_nnode = 10
        real, dimension(1:2) :: salt_top_z = [0.5, 0.8]
        real :: salt_vp = 5000.0
        real :: salt_rho = 2150.0
        real, allocatable, dimension(:, :, :) :: salt
        real :: salt_top_height = 20.0

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
        real, allocatable, dimension(:, :, :) :: karst

        !==============================================================================================
        !> Elastic
        logical :: yn_elastic = .false.
        !> Vp/Vs ratios
        real, dimension(1:2) :: vpvsratio = [1.5, 1.8]
        !> Elastic models
        real, allocatable, dimension(:, :, :) :: vp, vs, rho
        real :: rho_a = 310.0, rho_b = 0.25, rho_c = 0.0
        !> Elastic images
        real, allocatable, dimension(:, :, :) :: image_pp, image_ps, image_sp, image_ss
        !> Salt body's Vs
        real :: salt_vs = 4400.0
        !> Is salt before or after unconformity?
        logical :: salt_before_unconf = .true.

    contains

        procedure, private :: create_psf
        procedure, private :: generate_image
        procedure, private :: generate_image_elastic
        procedure, public :: generate => generate_3d

    end type rgm3_curved

    private
    public :: rgm3_curved

contains

    subroutine create_psf(this, n1, n2, n3, freq)

        class(rgm3_curved), intent(inout) :: this
        integer :: n1, n2, n3
        real, intent(in), optional :: freq

        real, allocatable, dimension(:) :: wavelet, psf1, psf2, psf3
        real, allocatable, dimension(:, :, :) :: psf
        real :: f0, wt
        integer :: i, j, k

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
                wavelet = fourier_filt(wavelet, this%dt, this%wave_filt_freqs, this%wave_filt_amps)
            end if
        end if

        wavelet = wavelet/norm2(wavelet)

        if (.not. this%custom_psf) then
            psf = zeros(n1, n2, n3)
            psf1 = zeros(n1)
            psf2 = zeros(n2)
            psf3 = zeros(n3)
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
            if (this%psf_sigma(3) == 0) then
                !$omp parallel do private(k)
                do k = 1, n3
                    psf3(k) = exp(-0.5*(k - 1.0 - (n3 - 1.0)/2.0)**2)
                end do
                !$omp end parallel do
                where (psf3 < maxval(psf3))
                    psf3 = 0.0
                end where
            else
                !$omp parallel do private(k)
                do k = 1, n3
                    psf3(k) = exp(-0.5*(k - 1.0 - (n3 - 1.0)/2.0)**2/this%psf_sigma(3)**2)
                end do
                !$omp end parallel do
            end if
            !$omp parallel do private(i, j, k) collapse(3)
            do k = 1, n3
                do j = 1, n2
                    do i = 1, n1
                        psf(i, j, k) = wavelet(i)*psf1(i)*psf2(j)*psf3(k)
                    end do
                end do
            end do
            !$omp end parallel do
            this%psf = psf/norm2(psf)
        else
            call assert(size(this%psf, 1) == this%n1 .and. size(this%psf, 2) == this%n2 &
                .and. size(this%psf, 3) == this%n3, ' Error: shape of psf must be n1 x n2 x n3')
        end if

    end subroutine create_psf

    subroutine generate_image_elastic(this, vp, vs, rho)

        class(rgm3_curved), intent(inout) :: this
        real, dimension(:, :, :), intent(in) :: vp, vs, rho

        integer :: n1, n2, n3, i, j, k, l
        real, allocatable, dimension(:) :: rfc
        real, allocatable, dimension(:, :, :) :: ww

        n1 = this%n1
        n2 = this%n2
        n3 = this%n3

        this%image_pp = zeros(this%n1, this%n2, this%n3)
        this%image_ps = zeros(this%n1, this%n2, this%n3)
        this%image_sp = zeros(this%n1, this%n2, this%n3)
        this%image_ss = zeros(this%n1, this%n2, this%n3)

        rfc = zeros(4)
        !$omp parallel do private(i, j, k, l, rfc)
        do k = 1, n3
            do j = 1, n2
                do i = 1, n1
                    rfc = 0
                    do l = 0, 5
                        rfc = rfc + elastic_reflection_coefs(real(sin((l*3.0)*const_deg2rad)/vp(i, j, k)), &
                            vp(i, j, k), vs(i, j, k), rho(i, j, k), vp(i + 1, j, k), vs(i + 1, j, k), rho(i + 1, j, k))
                    end do
                    rfc = rfc/6.0
                    this%image_pp(i, j, k) = rfc(1)
                    this%image_ps(i, j, k) = rfc(2)
                    this%image_sp(i, j, k) = rfc(3)
                    this%image_ss(i, j, k) = rfc(4)
                end do
            end do
        end do
        !$omp end parallel do

        ! Add random noise
        if (this%noise_level /= 0 .and. this%yn_conv_noise) then
            select case (this%noise_type)

                case ('normal', 'gaussian', 'uniform', 'exp')
                    ww = gauss_filt(random(n1, n2, n3, dist=this%noise_type, seed=safe_seed(int(this%seed, 8)*23)), this%noise_smooth)
                    ww = ww - mean(ww)
                    ww = ww/maxval(abs(ww))*this%noise_level*maxval(abs(this%image_pp))
                    this%image_pp = this%image_pp + ww

                    ww = gauss_filt(random(n1, n2, n3, dist=this%noise_type, seed=safe_seed(int(this%seed, 8)*23 + 1)), this%noise_smooth)
                    ww = ww - mean(ww)
                    ww = ww/maxval(abs(ww))*this%noise_level*maxval(abs(this%image_ps))
                    this%image_ps = this%image_ps + ww

                    ww = gauss_filt(random(n1, n2, n3, dist=this%noise_type, seed=safe_seed(int(this%seed, 8)*23 + 2)), this%noise_smooth)
                    ww = ww - mean(ww)
                    ww = ww/maxval(abs(ww))*this%noise_level*maxval(abs(this%image_sp))
                    this%image_sp = this%image_sp + ww

                    ww = gauss_filt(random(n1, n2, n3, dist=this%noise_type, seed=safe_seed(int(this%seed, 8)*23 + 3)), this%noise_smooth)
                    ww = ww - mean(ww)
                    ww = ww/maxval(abs(ww))*this%noise_level*maxval(abs(this%image_ss))
                    this%image_ss = this%image_ss + ww

                case ('wavenumber')
                    this%image_pp = this%image_pp + noise_wavenumber(this%image_pp, this%noise_level, this%noise_smooth, safe_seed(int(this%seed, 8)*23))
                    this%image_ps = this%image_ps + noise_wavenumber(this%image_ps, this%noise_level, this%noise_smooth, safe_seed(int(this%seed, 8)*23 + 1))
                    this%image_sp = this%image_sp + noise_wavenumber(this%image_sp, this%noise_level, this%noise_smooth, safe_seed(int(this%seed, 8)*23 + 2))
                    this%image_ss = this%image_ss + noise_wavenumber(this%image_ss, this%noise_level, this%noise_smooth, safe_seed(int(this%seed, 8)*23 + 3))

            end select
        end if

        if (this%wave /= '') then

            call this%create_psf(n1, n2, n3, this%f0)
            this%image_pp = conv(this%image_pp, this%psf, 'same')
            call this%create_psf(n1, n2, n3, this%f0*mean(0.5*(ones(n1, n2, n3) + vp(1:this%n1, :, :)/vs(1:this%n1, :, :))))
            this%image_ps = conv(this%image_ps, this%psf, 'same')
            this%image_sp = conv(this%image_sp, this%psf, 'same')
            call this%create_psf(n1, n2, n3, this%f0*mean(vp(1:this%n1, :, :)/vs(1:this%n1, :, :)))
            this%image_ss = conv(this%image_ss, this%psf, 'same')

        end if

        ! Add random noise
        if (this%noise_level /= 0 .and. (.not. this%yn_conv_noise)) then
            select case (this%noise_type)

                case ('normal', 'gaussian', 'uniform', 'exp')
                    ww = gauss_filt(random(n1, n2, n3, dist=this%noise_type, seed=safe_seed(int(this%seed, 8)*23)), this%noise_smooth)
                    ww = ww - mean(ww)
                    ww = ww/maxval(abs(ww))*this%noise_level*maxval(abs(this%image_pp))
                    this%image_pp = this%image_pp + ww

                    ww = gauss_filt(random(n1, n2, n3, dist=this%noise_type, seed=safe_seed(int(this%seed, 8)*23 + 1)), this%noise_smooth)
                    ww = ww - mean(ww)
                    ww = ww/maxval(abs(ww))*this%noise_level*maxval(abs(this%image_ps))
                    this%image_ps = this%image_ps + ww

                    ww = gauss_filt(random(n1, n2, n3, dist=this%noise_type, seed=safe_seed(int(this%seed, 8)*23 + 2)), this%noise_smooth)
                    ww = ww - mean(ww)
                    ww = ww/maxval(abs(ww))*this%noise_level*maxval(abs(this%image_sp))
                    this%image_sp = this%image_sp + ww

                    ww = gauss_filt(random(n1, n2, n3, dist=this%noise_type, seed=safe_seed(int(this%seed, 8)*23 + 3)), this%noise_smooth)
                    ww = ww - mean(ww)
                    ww = ww/maxval(abs(ww))*this%noise_level*maxval(abs(this%image_ss))
                    this%image_ss = this%image_ss + ww

                case ('wavenumber')
                    this%image_pp = this%image_pp + noise_wavenumber(this%image_pp, this%noise_level, this%noise_smooth, safe_seed(int(this%seed, 8)*23))
                    this%image_ps = this%image_ps + noise_wavenumber(this%image_ps, this%noise_level, this%noise_smooth, safe_seed(int(this%seed, 8)*23 + 1))
                    this%image_sp = this%image_sp + noise_wavenumber(this%image_sp, this%noise_level, this%noise_smooth, safe_seed(int(this%seed, 8)*23 + 2))
                    this%image_ss = this%image_ss + noise_wavenumber(this%image_ss, this%noise_level, this%noise_smooth, safe_seed(int(this%seed, 8)*23 + 3))

            end select

        end if

    end subroutine generate_image_elastic

    subroutine generate_image(this, vp, rho)

        class(rgm3_curved), intent(inout) :: this
        real, dimension(:, :, :), intent(in) :: vp, rho

        integer :: n1, n2, n3, i, j, k
        real, allocatable, dimension(:, :, :) :: ww

        n1 = this%n1
        n2 = this%n2
        n3 = this%n3

        this%image = zeros(n1, n2, n3)

        !$omp parallel do private(i, j, k)
        do k = 1, this%n3
            do j = 1, this%n2
                do i = 1, this%n1
                    this%image(i, j, k) = (vp(i + 1, j, k)*rho(i + 1, j, k) - vp(i, j, k)*rho(i, j, k)) &
                        /(vp(i + 1, j, k)*rho(i + 1, j, k) + vp(i, j, k)*rho(i, j, k))
                end do
            end do
        end do
        !$omp end parallel do

        ! Add random noise
        if (this%noise_level /= 0 .and. this%yn_conv_noise) then
            select case (this%noise_type)

                case ('normal', 'gaussian', 'uniform', 'exp')
                    ww = gauss_filt(random(n1, n2, n3, dist=this%noise_type, seed=safe_seed(int(this%seed, 8)*23)), this%noise_smooth)
                    ww = ww - mean(ww)
                    ww = ww/maxval(abs(ww))*this%noise_level*maxval(abs(this%image))
                    this%image = this%image + ww

                case ('wavenumber')
                    this%image = this%image + noise_wavenumber(this%image, this%noise_level, this%noise_smooth, safe_seed(int(this%seed, 8)*23))

            end select
        end if

        if (this%wave /= '') then

            call this%create_psf(n1, n2, n3, this%f0)
            this%image = conv(this%image, this%psf, 'same')

        end if

        ! Add random noise
        if (this%noise_level /= 0 .and. (.not. this%yn_conv_noise)) then
            select case (this%noise_type)

                case ('normal', 'gaussian', 'uniform', 'exp')
                    ww = gauss_filt(random(n1, n2, n3, dist=this%noise_type, seed=safe_seed(int(this%seed, 8)*23)), this%noise_smooth)
                    ww = ww - mean(ww)
                    ww = ww/maxval(abs(ww))*this%noise_level*maxval(abs(this%image))
                    this%image = this%image + ww

                case ('wavenumber')
                    this%image = this%image + noise_wavenumber(this%image, this%noise_level, this%noise_smooth, safe_seed(int(this%seed, 8)*23))

            end select
        end if

    end subroutine generate_image

    subroutine generate_3d(this)

        class(rgm3_curved), intent(inout) :: this

        if (this%nf == 0) then
            this%yn_fault = .false.
        end if

        if (this%unconf == 0) then
            call generate_3d_geological_model(this)
        else
            call generate_3d_unconformal_geological_model(this)
        end if

    end subroutine generate_3d

    subroutine generate_3d_geological_model(this)

        type(rgm3_curved), intent(inout) :: this

        real, allocatable, dimension(:, :, :) :: w, f, t, m
        real, allocatable, dimension(:, :, :) :: ww, ff, cf, tt, lz, mm
        integer :: nf, nl, n1, n2, n3
        real, allocatable, dimension(:) :: disp, dip, strike, rake, f2, f3, vv
        integer :: newi, newj, newk
        real, allocatable, dimension(:, :) :: r, rt
        integer :: i, j, k, fi, ne1, ne2, ne3, l
        real, allocatable, dimension(:) :: plw
        real, allocatable, dimension(:, :, :) :: fdip, fstrike, frake, fdisp
        real, allocatable, dimension(:, :, :) :: ffdip, ffstrike, ffrake, ffdisp
        real :: x0, y0, fwidth, dt
        real, allocatable, dimension(:) :: mu2, sigma2, mu3, sigma3, height, gtheta
        real, dimension(1:2) :: gmu2, gmu3, gsigma2, gsigma3
        real, allocatable, dimension(:) :: sumdisp, rc
        real :: theta
        real, allocatable, dimension(:) :: delta_dip, delta_strike
        real, dimension(1:2) :: pt
        real, allocatable, dimension(:, :) :: dips
        integer :: nsf
        real :: xc, yc, cphi, sphi, u0, v0, umin, umax
        real :: dev, dev_prev, ddev, theta_loc, dp, gc, wu, vdist, vdist_prev, sl
        integer :: iu, nu
        logical :: yn_mark
        real, allocatable, dimension(:) :: dphi, gcurve, devs
        real, allocatable, dimension(:) :: disp_au, disp_az, disp_uc, disp_zc, decay_w
        real :: alpha, dloc, s1, s2, s3, cr, sr, cd, s2c, s3c
        logical :: yn_pointwise
        integer :: ia, ja, ka
        real :: wa, wb, wc
        type(fractal_noise_1d) :: sn
        logical, allocatable, dimension(:, :, :) :: fblock
        real, allocatable, dimension(:) :: x1, x2, rds, pds, vds, salt_radius
        integer :: nd, isalt
        real, allocatable, dimension(:, :) :: slice, topz
        real :: dist, tp
        integer :: ag
        integer, allocatable, dimension(:) :: gmax, hmax
        type(fractal_noise_2d) :: pn
        type(fractal_noise_1d) :: qn
        real, allocatable, dimension(:, :, :) :: pxy
        real :: thick
        real, allocatable, dimension(:, :, :) :: vp, vs, rho
        real, allocatable, dimension(:) :: rfc
        real :: xcenter, ycenter
        real :: m1, m2, m3

        fwidth = this%fwidth
        dt = 0.001

        n1 = this%n1
        n2 = this%n2
        n3 = this%n3
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
            if (allocated(this%strike)) then
                call assert(size(this%strike) >= 2 .and. mod(size(this%strike), 2) == 0 .and. size(this%strike) == size(this%dip), &
                    ' <generate_2d_geological_model> Error: fault strike is not specified properly. ')
            else
                this%strike = tile([0.0, 180.0], size(this%dip)/2)
            end if
            if (allocated(this%rake)) then
                call assert(size(this%rake) >= 2 .and. mod(size(this%rake), 2) == 0 .and. size(this%rake) == size(this%dip), &
                    ' <generate_2d_geological_model> Error: fault rake is not specified properly. ')
            else
                this%rake = tile([0.0, 180.0], size(this%dip)/2)
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
                strike = [random(nint(nf/2.0), range=[0.975, 1.025]*this%strike(1), seed=safe_seed(int(this%seed, 8)*2))*const_deg2rad, &
                    random(nf - nint(nf/2.0), range=[0.975, 1.025]*this%strike(2), seed=safe_seed(int(this%seed, 8)*2))*const_deg2rad]
                rake = [random(nint(nf/2.0), range=[0.95, 1.05]*this%rake(1), seed=safe_seed(int(this%seed, 8)*3))*const_deg2rad, &
                    random(nf - nint(nf/2.0), range=[0.95, 1.05]*this%rake(2), seed=safe_seed(int(this%seed, 8)*3))*const_deg2rad]
                disp = [random(nint(nf/2.0), range=[0.9, 1.1]*this%disp(1), seed=safe_seed(int(this%seed, 8)*4)), &
                    random(nf - nint(nf/2.0), range=[0.9, 1.1]*this%disp(2), seed=safe_seed(int(this%seed, 8)*4))]

                strike = clip(strike, 0.0, real(const_pi))
                rake = clip(rake, 0.0, real(const_pi))

            else

                ! Dip, strike, rake angles, and fault displacements
                nsf = size(this%dip)/2
                dip = random(ceiling(nf*1.0/nsf), range=this%dip(1:2), seed=this%seed)
                strike = random(ceiling(nf*1.0/nsf), range=this%strike(1:2), seed=safe_seed(int(this%seed, 8)*2))
                rake = random(ceiling(nf*1.0/nsf), range=this%rake(1:2), seed=safe_seed(int(this%seed, 8)*3))
                disp = random(ceiling(nf*1.0/nsf), range=this%disp(1:2), seed=safe_seed(int(this%seed, 8)*4))
                do i = 2, nsf
                    dip = [dip, random(ceiling(nf*1.0/nsf), range=this%dip((i - 1)*2 + 1:i*2), seed=this%seed)]
                    strike = [strike, random(ceiling(nf*1.0/nsf), range=this%strike((i - 1)*2 + 1:i*2), seed=safe_seed(int(this%seed, 8)*2))]
                    rake = [rake, random(ceiling(nf*1.0/nsf), range=this%rake((i - 1)*2 + 1:i*2), seed=safe_seed(int(this%seed, 8)*2))]
                    disp = [disp, random(ceiling(nf*1.0/nsf), range=this%disp((i - 1)*2 + 1:i*2), seed=safe_seed(int(this%seed, 8)*4))]
                end do
                dip = dip(1:nf)
                dip = dip*const_deg2rad

                strike = strike(1:nf)
                strike = strike*const_deg2rad

                rake = rake(1:nf)
                rake = rake*const_deg2rad

                disp = disp(1:nf)
                disp(1:nf:2) = -disp(1:nf:2)

            end if
        else
            dip = zeros(1)
            strike = zeros(1)
            rake = zeros(1)
            disp = zeros(1)
        end if

        ! The dimensions of the padded model
        sumdisp = disp*(-sin(rake)*sin(dip))
        m1 = max(sum(sumdisp, mask=sumdisp > 0), -sum(sumdisp, mask=sumdisp < 0))
        m2 = max(maxval(abs(this%refl_slope)), maxval(abs(this%refl_slope_top)))
        m3 = max(maxval(abs(this%refl_height)), maxval(abs(this%refl_height_top)))
        ne1 = ceiling(max(m1, m2) + m3)
        n1 = n1 + 2*ne1

        if (maxval(abs(this%delta_strike)) > 0) then
            ! For strike-varying faults, the lateral slip direction rotates with
            ! the local strike, so pad with the strike-independent bound of the
            ! lateral slip components
            sumdisp = disp*sqrt(cos(rake)**2 + (sin(rake)*cos(dip))**2)
            ne2 = ceiling(max(sum(sumdisp, mask=sumdisp > 0), -sum(sumdisp, mask=sumdisp < 0)))
            ne3 = ne2
        else
            sumdisp = disp*(cos(rake)*sin(strike) - sin(rake)*cos(dip)*cos(strike))
            ne2 = ceiling(max(sum(sumdisp, mask=sumdisp > 0), -sum(sumdisp, mask=sumdisp < 0)))

            sumdisp = disp*(cos(rake)*cos(strike) + sin(rake)*cos(dip)*sin(strike))
            ne3 = ceiling(max(sum(sumdisp, mask=sumdisp > 0), -sum(sumdisp, mask=sumdisp < 0)))
        end if
        n2 = n2 + 2*ne2
        n3 = n3 + 2*ne3

        ! Compute the top and bottom fault dip angles
        ! Note that to ensure the faults have proper dip angle range within the final cropped image,
        ! here I add some extra degrees to the begin and end dip angles
        dips = zeros(n1, nf)
        delta_dip = random(nf, range=this%delta_dip, seed=safe_seed((int(this%seed, 8)*9 + 1)/2))*const_deg2rad
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

        ! Reflector's shape
        select case (this%refl_shape)

            case ('random')
                r = random(n2, n3, dist='normal', seed=safe_seed(int(this%seed, 8)*5))
                r = gauss_filt(r, [this%refl_smooth, this%refl_smooth])

            case ('gaussian', 'cauchy')
                if (maxval(this%refl_mu2) == 0) then
                    gmu2 = [1.0, this%n2 - 1.0]
                else
                    gmu2 = this%refl_mu2
                end if

                if (maxval(this%refl_sigma2) == 0) then
                    gsigma2 = [0.05, 0.15]*n2
                else
                    gsigma2 = this%refl_sigma2
                end if

                if (maxval(this%refl_mu3) == 0) then
                    gmu3 = [1.0, this%n3 - 1.0]
                else
                    gmu3 = this%refl_mu3
                end if

                if (maxval(this%refl_sigma3) == 0) then
                    gsigma3 = [0.05, 0.15]*n3
                else
                    gsigma3 = this%refl_sigma3
                end if

                mu2 = random(this%ng, range=gmu2, seed=safe_seed(int(this%seed, 8)*5))
                sigma2 = random(this%ng, range=gsigma2, seed=safe_seed(int(this%seed, 8)*6))
                mu3 = random(this%ng, range=gmu3, seed=safe_seed(int(this%seed, 8)*7))
                sigma3 = random(this%ng, range=gsigma3, seed=safe_seed(int(this%seed, 8)*8))
                height = random(this%ng, range=this%refl_height, seed=safe_seed(int(this%seed, 8)*9))
                gtheta = random(this%ng, range=[0.0, 180.0], seed=safe_seed(int(this%seed, 8)*9 - 1))*ifelse(this%rotate_fold, real(const_deg2rad), 0.0)

                r = zeros(n2, n3)
                do i = 1, this%ng
                    select case (this%refl_shape)
                        case ('gaussian')
                            r = r + rescale(gaussian(linspace(0.0, n2 - 1.0, n2), linspace(0.0, n3 - 1.0, n3), &
                                [mu2(i) + ne2, mu3(i) + ne3], [sigma2(i), sigma3(i)], gtheta(i)), [0.0, height(i)])
                        case ('cauchy')
                            r = r + rescale(cauchy(linspace(0.0, n2 - 1.0, n2), linspace(0.0, n3 - 1.0, n3), &
                                [mu2(i) + ne2, mu3(i) + ne3], [sigma2(i), sigma3(i)], gtheta(i)), [0.0, height(i)])
                    end select
                end do

            case ('perlin')
                pn%n1 = n2
                pn%n2 = n3
                pn%seed = safe_seed(int(this%seed, 8)*5)
                r = gauss_filt(pn%generate(), [this%refl_smooth, this%refl_smooth])

            case ('custom')
                call assert(allocated(this%refl), &
                    ' <generate_2d_geological_model> Error: refl must be initialized. ')
                call assert(size(this%refl, 1) == this%n2 .and. size(this%refl, 2) == this%n3, &
                    '<generate_2d_geological_model> Error: size(refl) must = (n2, n3)')
                r = pad(this%refl, [ne2 + 1, ne2 + 1, n3 + 1, n3 + 1], ['edge', 'edge', 'edge', 'edge'])

        end select

        if (this%refl_shape_top /= 'same') then

            select case (this%refl_shape_top)

                case default
                    rt = random(n2, n3, dist='normal', seed=safe_seed(int(this%seed, 8)*5 - 1))
                    rt = gauss_filt(rt, [this%refl_smooth_top, this%refl_smooth_top])

                case ('random')
                    rt = random(n2, n3, dist='normal', seed=safe_seed(int(this%seed, 8)*5 - 1))
                    rt = gauss_filt(rt, [this%refl_smooth_top, this%refl_smooth_top])

                case ('gaussian', 'cauchy')
                    if (maxval(this%refl_mu2) == 0) then
                        gmu2 = [1.0, this%n2 - 1.0]
                    else
                        gmu2 = this%refl_mu2
                    end if

                    if (maxval(this%refl_sigma2) == 0) then
                        gsigma2 = [0.05, 0.15]*n2
                    else
                        gsigma2 = this%refl_sigma2
                    end if

                    if (maxval(this%refl_mu3) == 0) then
                        gmu3 = [1.0, this%n3 - 1.0]
                    else
                        gmu3 = this%refl_mu3
                    end if

                    if (maxval(this%refl_sigma3) == 0) then
                        gsigma3 = [0.05, 0.15]*n3
                    else
                        gsigma3 = this%refl_sigma3
                    end if

                    mu2 = random(this%ng, range=gmu2, seed=safe_seed(int(this%seed, 8)*5 - 1))
                    sigma2 = random(this%ng, range=gsigma2, seed=safe_seed(int(this%seed, 8)*6 - 1))
                    mu3 = random(this%ng, range=gmu3, seed=safe_seed(int(this%seed, 8)*7 - 1))
                    sigma3 = random(this%ng, range=gsigma3, seed=safe_seed(int(this%seed, 8)*8 - 1))
                    height = random(this%ng, range=this%refl_height, seed=safe_seed(int(this%seed, 8)*9 - 1))
                    gtheta = random(this%ng, range=[0.0, 180.0], seed=safe_seed(int(this%seed, 8)*9 - 2))*ifelse(this%rotate_fold, real(const_deg2rad), 0.0)

                    rt = zeros(n2, n3)
                    do i = 1, this%ng
                        select case (this%refl_shape)
                            case ('gaussian')
                                rt = rt + rescale(gaussian(linspace(0.0, n2 - 1.0, n2), linspace(0.0, n3 - 1.0, n3), &
                                    [mu2(i) + ne2, mu3(i) + ne3], [sigma2(i), sigma3(i)], gtheta(i)), [0.0, height(i)])
                            case ('cauchy')
                                rt = rt + rescale(cauchy(linspace(0.0, n2 - 1.0, n2), linspace(0.0, n3 - 1.0, n3), &
                                    [mu2(i) + ne2, mu3(i) + ne3], [sigma2(i), sigma3(i)], gtheta(i)), [0.0, height(i)])
                        end select
                    end do

                case ('perlin')
                    pn%n1 = n2
                    pn%n2 = n3
                    pn%seed = safe_seed(int(this%seed, 8)*5 - 1)
                    rt = gauss_filt(pn%generate(), [this%refl_smooth_top, this%refl_smooth_top])

                case ('custom')
                    call assert(allocated(this%refl_top), &
                        ' <generate_2d_geological_model> Error: refl_top must be initialized. ')
                    call assert(size(this%refl_top, 1) == this%n2 .and. size(this%refl_top, 2) == this%n3, &
                        '<generate_2d_geological_model> Error: size(refl_top) must = (n2, n3)')
                    rt = pad(this%refl_top, [ne2 + 1, ne2 + 1, n3 + 1, n3 + 1], ['edge', 'edge', 'edge', 'edge'])

            end select

        else

            rt = r

        end if

        ! Rescale reflectors to their height
        r = rescale(r, this%refl_height*rov(r)/(rov(r(ne2 + 1:ne2 + this%n2, ne3 + 1:ne3 + this%n3)) + float_tiny))
        rt = rescale(rt, this%refl_height_top*rov(rt)/(rov(rt(ne2 + 1:ne2 + this%n2, ne3 + 1:ne3 + this%n3)) + float_tiny))

        ! Add slopes
        !$omp parallel do private(j, k) collapse(2)
        do k = 1, n3
            do j = 1, n2
                r(j, k) = r(j, k) + (j - 1.0)*this%refl_slope(1)/this%n2 &
                    + (k - 1.0)*this%refl_slope(2)/this%n3
                rt(j, k) = rt(j, k) + (j - 1.0)*this%refl_slope_top(1)/this%n2 &
                    + (k - 1.0)*this%refl_slope_top(2)/this%n3
            end do
        end do
        !$omp end parallel do

        ! ... and get the final positions of top/bottom reflectors
        r = r - mean(r)
        rt = rt - mean(rt) + n1

        nl = nint(this%nl + this%nl*2.0*ne1/(n1 - 2*ne1))

        vv = random(nl - 1, seed=safe_seed(int(this%seed, 8)*6))*this%delta_v
        vv = linspace(this%vmax*(1.0 + ne1*1.0/this%n1), this%vmin*(1.0 - ne1*1.0/this%n1), nl - 1) + vv

        lz = zeros(nl, n2, n3)
        plw = random(nl, range=[1 - this%lwv, 1 + this%lwv], seed=safe_seed(int(this%seed, 8)*7))

        !$omp parallel do private(j, k)
        do k = 1, n3
            do j = 1, n2
                lz(:, j, k) = ginterp([0.0, (ne1 + 1.0)/n1, (n1 - ne1 - 1.0)/n1, 1.0], &
                    [r(j, k), r(j, k) + ne1, rt(j, k) - ne1, rt(j, k)], &
                    linspace(0.0, 1.0, nl), 'linear')
            end do
        end do
        !$omp end parallel do

        lz = deriv(lz, dim=1)
        !$omp parallel do private(j, k, l, tp)
        do k = 1, n3
            do j = 1, n2
                lz(:, j, k) = r(j, k) + rescale(cumsum(lz(:, j, k)*plw), [0.0, rt(j, k) - r(j, k)])
            end do
        end do
        !$omp end parallel do

        ! Add horizontal layer thickness variation
        if (this%lwh > 0) then

            pxy = zeros(nl, n2, n3)

            !$omp parallel do private(i, thick)
            do i = 1, nl - 1
                thick = (lz(i + 1, 1, 1) - lz(i, 1, 1))*this%lwh*2.0
                ! The seed here must be sign-consistent, e.g., when = -1, all must < 0; otherwise, all must > 0
                pxy(i, :, :) = gauss_filt(random(n2, n3, seed=safe_seed(int(this%seed, 8)*nl - i)), [1, 1]*this%secondary_refl_smooth)
                pxy(i, :, :) = rescale(pxy(i, :, :), [-thick, thick]*(0.5 + (nl - i + 0.0)/(nl - 1.0)*0.5))
            end do
            !$omp end parallel do

            lz = deriv(lz + pxy, dim=1)
            where (lz <= 0)
                lz = 1.0e-9
            end where
            !$omp parallel do private(j, k)
            do k = 1, n3
                do j = 1, n2
                    lz(:, j, k) = r(j, k) + rescale(cumsum(lz(:, j, k)), [0.0, rt(j, k) - r(j, k)])
                end do
            end do
            !$omp end parallel do

        end if

        ! Stack layers to create a layer-constant velocity model
        w = zeros(n1, n2, n3)
        if (this%yn_facies) then
            m = zeros(n1, n2, n3)
        end if
        if (this%yn_rgt) then
            t = zeros(n1, n2, n3)
        end if

        nd = nint(n1*2.0/nl)
        rc = zeros((nl - 1)*nd)
        rfc = zeros((nl - 1)*nd)

        !$omp parallel do private(i, j, k, l, rc, rfc, tp)
        do k = 1, n3
            do j = 1, n2

                ! Velocity model is interpolated vertically based on reflectors
                ! Otherwise there will be staircases in the generated velocity model
                ! And correspondingly the images can has notable staircases which are unrealistic
                do l = 1, nl - 1
                    tp = lz(l + 1, j, k) - lz(l, j, k)
                    rc((l - 1)*nd + 1:l*nd) = linspace(lz(l, j, k) + tp/(nd + 1), lz(l + 1, j, k) - tp/(nd + 1), nd)
                    rfc((l - 1)*nd + 1:l*nd) = vv(l)
                end do
                w(:, j, k) = ginterp(rc, rfc, linspace(n1 - 1.0, 0.0, n1), 'linear')

                ! Facies is piecewise constant and has distinct values for each layer
                if (this%yn_facies) then
                    do i = 1, n1
                        loop_layer: do l = 1, nl - 1
                            if (n1 - i >= lz(l, j, k) .and. n1 - i < lz(l + 1, j, k)) then
                                m(i, j, k) = l
                                exit loop_layer
                            end if
                        end do loop_layer
                    end do
                end if

                ! RGT is linearly interpolated based on the location of reflectors
                if (this%yn_rgt) then
                    t(:, j, k) = ginterp(n1 - 1.0 - lz(:, j, k), linspace(1.0, 0.0, nl), linspace(0.0, n1 - 1.0, n1), 'linear')
                end if

            end do
        end do
        !$omp end parallel do

        ! Add faults
        f = zeros(n1, n2, n3)
        ff = zeros(n1, n2, n3)
        cf = zeros(n1, n2, n3)
        ww = zeros(n1, n2, n3)

        fdip = zeros(n1, n2, n3)
        fstrike = zeros(n1, n2, n3)
        frake = zeros(n1, n2, n3)
        fdisp = zeros(n1, n2, n3)

        ffdip = zeros(n1, n2, n3)
        ffstrike = zeros(n1, n2, n3)
        ffrake = zeros(n1, n2, n3)
        ffdisp = zeros(n1, n2, n3)

        if (this%yn_facies) then
            mm = zeros(n1, n2, n3)
        end if
        if (this%yn_rgt) then
            tt = zeros(n1, n2, n3)
        end if

        if (nf >= 1) then

            if (this%yn_regular_fault) then

                rc = random(nf - 1, range=[0.75, 1.25]*0.9*this%n2/(nf - 1.0), seed=safe_seed(int(this%seed, 8)*12))
                rc = rc*0.9*this%n2/sum(rc)
                rc = sort(rc)
                f2 = zeros(nf)
                f2(1) = ne2 + rand(range=[0.2, 0.8]*(this%n2 - sum(rc)), seed=safe_seed(int(this%seed, 8)*12 - 1))
                f2(2:) = f2(1) + cumsum(rc)
                if (.not. this%yn_group_faults) then
                    f2 = random_permute(f2, seed=safe_seed(int(this%seed, 8)*12 - 2))
                end if

                rc = random(nf - 1, range=[0.75, 1.25]*0.9*this%n3/(nf - 1.0), seed=safe_seed(int(this%seed, 8)*13))
                rc = rc*0.9*this%n3/sum(rc)
                rc = sort(rc)
                f3 = zeros(nf)
                f3(1) = ne3 + rand(range=[0.2, 0.8]*(this%n3 - sum(rc)), seed=safe_seed(int(this%seed, 8)*13 - 1))
                f3(2:) = f3(1) + cumsum(rc)
                if (.not. this%yn_group_faults) then
                    f3 = random_permute(f3, seed=safe_seed(int(this%seed, 8)*13 - 2))
                end if

                ! Rotate the points of fault centers
                do i = 1, nf
                    ! Correct f2, f3, so that the quasi-parallel fault pattern occurs in the depth center rather than at the surface
                    f2(i) = f2(i) - tan(dip(i) - const_pi_half)*n1*0.5*cos(strike(i))
                    f3(i) = f3(i) - tan(dip(i) - const_pi_half)*n1*0.5*sin(strike(i))
                    pt = rotate_point([f2(i), f3(i)], real(const_pi_half) - strike(i), [0.5*(n2 - 1.0), 0.5*(n3 - 1.0)])
                    f2(i) = pt(1)
                    f3(i) = pt(2)
                end do

            else

                f2 = random(nf, range=[ne2 + 0.1*this%n2, n2 - ne2 - 0.1*this%n2], &
                    seed=safe_seed(int(this%seed, 8)*16), spacing=0.75*(n2 - 2*ne2 - 0.2*this%n2)/nf)
                f3 = random(nf, range=[ne3 + 0.1*this%n3, n3 - ne3 - 0.1*this%n3], &
                    seed=safe_seed(int(this%seed, 8)*17), spacing=0.75*(n3 - 2*ne3 - 0.2*this%n3)/nf)

            end if

            ! The maximum strike deviation along each fault
            delta_strike = random(nf, range=this%delta_strike, seed=safe_seed((int(this%seed, 8)*11 + 1)/2))*const_deg2rad

            ! Elliptical slip patches for spatially varying fault displacement;
            ! the semi-axes and centers are defined in the (along-strike, depth)
            ! coordinates of each fault
            disp_au = random(nf, range=this%disp_radius_strike, seed=safe_seed((int(this%seed, 8)*13 + 1)/2))*0.5*(this%n2 + this%n3)
            disp_az = random(nf, range=this%disp_radius_dip, seed=safe_seed((int(this%seed, 8)*15 + 1)/2))*this%n1
            disp_uc = random(nf, range=this%disp_center_strike, seed=safe_seed((int(this%seed, 8)*17 + 1)/2))*0.5*(this%n2 + this%n3)
            disp_zc = ne1 + random(nf, range=this%disp_center_dip, seed=safe_seed((int(this%seed, 8)*19 + 1)/2))*this%n1

            ! Widths of the Gaussian decay of displacement away from each fault
            decay_w = random(nf, range=this%disp_decay_width, seed=safe_seed((int(this%seed, 8)*21 + 1)/2))*this%n1

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
                ffstrike = fstrike
                ffrake = frake
                ffdisp = fdisp
                cf = 0

                ! Local coordinates: u is along the base strike direction and
                ! v is normal to it, with the origin at the fault center
                xc = f3(fi)
                yc = f2(fi)
                cphi = cos(strike(fi))
                sphi = sin(strike(fi))

                ! The u-axis range that covers the entire (padded) model
                umin = float_huge
                umax = -float_huge
                do k = 1, 2
                    do j = 1, 2
                        tp = ((k - 1)*(n3 - 1.0) - xc)*cphi + ((j - 1)*(n2 - 1.0) - yc)*sphi
                        umin = min(umin, tp)
                        umax = max(umax, tp)
                    end do
                end do
                nu = ceiling(umax - umin) + 2

                ! The local strike of the fault is strike + dphi(u), where dphi
                ! is a smooth random deviation angle along the fault, and the
                ! resulting lateral deviation of the fault trace in map view
                ! w.r.t. the straight trace is the along-u integration of tan(dphi(u))
                if (delta_strike(fi) /= 0) then
                    ! Low-wavenumber (smooth) variation of the strike along the fault;
                    ! a small number of octaves avoids geologically uncommon
                    ! high-wavenumber oscillations of the fault surface
                    sn%n1 = nu
                    sn%octaves = 2
                    sn%periods1 = max(1, this%strike_nperiod)
                    sn%seed = safe_seed(int(this%seed, 8)*27*fi - 1)
                    dphi = sn%generate()
                    dphi = dphi/maxval(abs(dphi))*delta_strike(fi)
                else
                    dphi = zeros(nu)
                end if
                gcurve = cumsum(tan(dphi))
                ! Let the curved trace pass through the fault center
                gcurve = gcurve - gcurve(clip(nint(-umin) + 1, 1, nu))

                fblock = falses(n1, n2, n3)
                dev_prev = 0.0

                do i = 1, n1

                    theta = dips(i, fi)

                    ! The lateral deviation of the curved (listric) fault
                    ! w.r.t. the topmost fault trace, along the base strike normal
                    dev = 1.0/tan(theta)*(i - 1.0)
                    if (i == 1) then
                        dev_prev = dev
                    end if
                    ddev = dev - dev_prev

                    !$omp parallel do private(j, k, x0, y0, u0, v0, iu, wu, dp, gc, theta_loc, sl, vdist, vdist_prev, yn_mark, alpha) collapse(2)
                    do k = 1, n3
                        do j = 1, n2

                            x0 = k - 1.0
                            y0 = j - 1.0

                            u0 = (x0 - xc)*cphi + (y0 - yc)*sphi
                            v0 = -(x0 - xc)*sphi + (y0 - yc)*cphi

                            if (delta_strike(fi) /= 0) then

                                ! Local strike deviation angle and trace deviation by linear interpolation
                                iu = clip(floor(u0 - umin) + 1, 1, nu - 1)
                                wu = clip(u0 - umin - (iu - 1.0), 0.0, 1.0)
                                dp = dphi(iu)*(1.0 - wu) + dphi(iu + 1)*wu
                                gc = gcurve(iu)*(1.0 - wu) + gcurve(iu + 1)*wu

                                ! Local true dip angle satisfying tan(theta_loc) = tan(theta)/cos(dp),
                                ! where the atan2 form is robust for theta = pi/2 and
                                ! preserves theta_loc in (0, pi)
                                theta_loc = atan2(sin(theta), cos(theta)*cos(dp))

                            else

                                ! Constant strike: the local geometry is the base geometry
                                dp = 0.0
                                gc = 0.0
                                theta_loc = theta

                            end if

                            ! Signed horizontal distance to the fault, measured
                            ! along the local strike normal
                            vdist = (v0 - gc - dev)*cos(dp)

                            ! Local displacement scaling from the elliptical slip patch
                            if (this%yn_vary_disp) then
                                alpha = max(0.0, 1.0 - ((u0 - disp_uc(fi))/disp_au(fi))**2 - ((i - disp_zc(fi))/disp_az(fi))**2)
                            else
                                alpha = 1.0
                            end if

                            yn_mark = abs(vdist) < 0.5*fwidth/sin(theta_loc)

                            ! Fill possible gaps caused by fast lateral deviation of a low-dip fault
                            if (.not. yn_mark .and. abs(ddev) >= sqrt(2.0)) then
                                vdist_prev = (v0 - gc - dev_prev)*cos(dp)
                                yn_mark = abs(vdist_prev) <= sqrt(2.0)*abs(ddev) .and. abs(vdist) <= sqrt(2.0)*abs(ddev)
                            end if

                            ! Where the local displacement diminishes below half a grid point,
                            ! the fault causes no visible offset, and the point is beyond
                            ! the fault tip line and therefore not marked as fault
                            if (this%yn_vary_disp) then
                                yn_mark = yn_mark .and. abs(disp(fi))*alpha >= 0.5
                            end if

                            if (yn_mark) then
                                ! Wrap the local strike angle into [0, pi)
                                sl = strike(fi) + dp
                                if (sl < 0) then
                                    sl = sl + const_pi
                                else if (sl >= const_pi) then
                                    sl = sl - const_pi
                                end if
                                f(i, j, k) = fi
                                cf(i, j, k) = 1.0
                                fdip(i, j, k) = theta_loc
                                fstrike(i, j, k) = sl
                                frake(i, j, k) = rake(fi)
                                ! Store the local displacement magnitude
                                fdisp(i, j, k) = abs(disp(fi))*alpha
                            end if

                            ! Mark the fault block (hanging wall side) to be shifted
                            if (theta < const_pi_half) then
                                fblock(i, j, k) = vdist > 0
                            else
                                fblock(i, j, k) = vdist < 0
                            end if

                        end do
                    end do
                    !$omp end parallel do

                    ! next iteration
                    dev_prev = dev

                end do

                ! Shift the fault block; here for each target point within the fault
                ! block, the values are gathered from its source point, which is
                ! equivalent to the previous shift-source-to-target (scatter) approach
                ! for a constant fault displacement, while it naturally supports
                ! spatially varying fault displacement without leaving holes in the model
                s1 = -sin(rake(fi))*sin(dip(fi))
                cr = cos(rake(fi))
                sr = sin(rake(fi))
                cd = cos(dip(fi))

                ! Fast path: when the strike does not vary along this fault and the
                ! displacement is constant, the slip vector and displacement are
                ! per-fault constants, and the per-point local-strike evaluation
                ! below is unnecessary
                yn_pointwise = delta_strike(fi) /= 0 .or. this%yn_vary_disp .or. this%yn_disp_decay
                s2c = cr*sin(strike(fi)) - sr*cd*cos(strike(fi))
                s3c = cr*cos(strike(fi)) + sr*cd*sin(strike(fi))

                ! Lateral deviations of the fault surface at all depths, for
                ! computing the distance to the fault when yn_disp_decay = .true.
                devs = zeros(n1)
                do i = 1, n1
                    devs(i) = 1.0/tan(dips(i, fi))*(i - 1.0)
                end do

                !$omp parallel do private(i, j, k, x0, y0, u0, v0, iu, wu, dp, gc, vdist, alpha, dloc, s2, s3, newi, newj, newk, ia, ja, ka, wa, wb, wc) collapse(3)
                do k = 1, n3
                    do j = 1, n2
                        do i = 1, n1

                            if (fblock(i, j, k)) then

                                if (yn_pointwise) then

                                    x0 = k - 1.0
                                    y0 = j - 1.0
                                    u0 = (x0 - xc)*cphi + (y0 - yc)*sphi

                                    ! Local strike deviation angle at the target point
                                    iu = clip(floor(u0 - umin) + 1, 1, nu - 1)
                                    wu = clip(u0 - umin - (iu - 1.0), 0.0, 1.0)
                                    dp = dphi(iu)*(1.0 - wu) + dphi(iu + 1)*wu

                                else

                                    dp = 0.0

                                end if

                                ! Local displacement from the elliptical slip patch,
                                ! maximum at the patch center and diminishing to zero
                                ! at the patch boundary (fault tip line)
                                if (this%yn_vary_disp) then
                                    alpha = max(0.0, 1.0 - ((u0 - disp_uc(fi))/disp_au(fi))**2 - ((i - disp_zc(fi))/disp_az(fi))**2)
                                else
                                    alpha = 1.0
                                end if

                                ! Decay the displacement away from the fault with a Gaussian
                                ! profile of the distance to the fault, mimicking the deformation
                                ! halo of a finite slip patch (fault drag, rollover, and
                                ! blind-fault folding); the displacement at the fault surface
                                ! itself is unaffected
                                if (this%yn_disp_decay) then
                                    v0 = -(x0 - xc)*sphi + (y0 - yc)*cphi
                                    gc = gcurve(iu)*(1.0 - wu) + gcurve(iu + 1)*wu
                                    vdist = (v0 - gc - devs(i))*cos(dp)
                                    alpha = alpha*exp(-0.5*(vdist/(decay_w(fi) + float_tiny))**2)
                                end if

                                dloc = disp(fi)*alpha

                                ! The slip direction follows the local strike of the fault
                                if (yn_pointwise) then
                                    s2 = cr*sin(strike(fi) + dp) - sr*cd*cos(strike(fi) + dp)
                                    s3 = cr*cos(strike(fi) + dp) + sr*cd*sin(strike(fi) + dp)
                                else
                                    s2 = s2c
                                    s3 = s3c
                                end if

                                ! The source point of the target point
                                newi = nint(i - s1*dloc)
                                newj = nint(j - s2*dloc)
                                newk = nint(k - s3*dloc)

                                if (newi >= 1 .and. newi <= n1 .and. newj >= 1 .and. newj <= n2 .and. newk >= 1 .and. newk <= n3) then

                                    ! Move velocity
                                    w(i, j, k) = ww(newi, newj, newk)

                                    ! Move facies
                                    if (this%yn_facies) then
                                        m(i, j, k) = mm(newi, newj, newk)
                                    end if

                                    ! Move RGT
                                    if (this%yn_rgt) then
                                        t(i, j, k) = tt(newi, newj, newk)
                                    end if

                                    ! Move existing faults
                                    if (cf(i, j, k) == 0) then
                                        f(i, j, k) = ff(newi, newj, newk)
                                        fdip(i, j, k) = ffdip(newi, newj, newk)
                                        fstrike(i, j, k) = ffstrike(newi, newj, newk)
                                        frake(i, j, k) = ffrake(newi, newj, newk)
                                        fdisp(i, j, k) = ffdisp(newi, newj, newk)
                                    end if

                                    ! For a spatially varying displacement, rounding the source
                                    ! point to the nearest grid point quantizes the warp into
                                    ! integer shifts, which produces terrace (ring-like) artifacts
                                    ! along the iso-displacement contours, especially visible on
                                    ! horizontal slices; therefore, the continuous fields
                                    ! (velocity and RGT) are gathered with trilinear interpolation
                                    ! instead, while the categorical fields (facies and fault
                                    ! attributes) keep the nearest-neighbor sampling
                                    if (this%yn_vary_disp .or. this%yn_disp_decay) then

                                        ia = clip(floor(i - s1*dloc), 1, n1 - 1)
                                        ja = clip(floor(j - s2*dloc), 1, n2 - 1)
                                        ka = clip(floor(k - s3*dloc), 1, n3 - 1)
                                        wa = clip(i - s1*dloc - ia, 0.0, 1.0)
                                        wb = clip(j - s2*dloc - ja, 0.0, 1.0)
                                        wc = clip(k - s3*dloc - ka, 0.0, 1.0)

                                        w(i, j, k) = &
                                            ww(ia, ja, ka)*(1 - wa)*(1 - wb)*(1 - wc) &
                                            + ww(ia + 1, ja, ka)*wa*(1 - wb)*(1 - wc) &
                                            + ww(ia, ja + 1, ka)*(1 - wa)*wb*(1 - wc) &
                                            + ww(ia, ja, ka + 1)*(1 - wa)*(1 - wb)*wc &
                                            + ww(ia + 1, ja + 1, ka)*wa*wb*(1 - wc) &
                                            + ww(ia + 1, ja, ka + 1)*wa*(1 - wb)*wc &
                                            + ww(ia, ja + 1, ka + 1)*(1 - wa)*wb*wc &
                                            + ww(ia + 1, ja + 1, ka + 1)*wa*wb*wc

                                        if (this%yn_rgt) then
                                            t(i, j, k) = &
                                                tt(ia, ja, ka)*(1 - wa)*(1 - wb)*(1 - wc) &
                                                + tt(ia + 1, ja, ka)*wa*(1 - wb)*(1 - wc) &
                                                + tt(ia, ja + 1, ka)*(1 - wa)*wb*(1 - wc) &
                                                + tt(ia, ja, ka + 1)*(1 - wa)*(1 - wb)*wc &
                                                + tt(ia + 1, ja + 1, ka)*wa*wb*(1 - wc) &
                                                + tt(ia + 1, ja, ka + 1)*wa*(1 - wb)*wc &
                                                + tt(ia, ja + 1, ka + 1)*(1 - wa)*wb*wc &
                                                + tt(ia + 1, ja + 1, ka + 1)*wa*wb*wc
                                        end if

                                    end if

                                end if

                            end if

                        end do
                    end do
                end do
                !$omp end parallel do

            end do

        end if

        ! Select model before adding salt bodies
        ! Note that n1 = this%n1 + 1 for computing reflectivity coefficients
        vp = w(ne1 + 1:ne1 + this%n1 + 1, ne2 + 1:ne2 + this%n2, ne3 + 1:ne3 + this%n3)
        vp = rescale(vp, [this%vmin, this%vmax])
        rho = this%rho_a*vp**this%rho_b + this%rho_c
        if (this%yn_elastic) then
            vs = vp/rescale(vp, this%vpvsratio)
        end if

        n1 = size(vp, 1)
        n2 = size(vp, 2)
        n3 = size(vp, 3)

        ! Add salt
        if (this%yn_salt) then

            if (maxval(this%salt_radius) == 0) then
                salt_radius = [0.05*0.5*(this%n2 + this%n3), 0.2*0.5*(this%n2 + this%n3)]
            else
                salt_radius = this%salt_radius
            end if
            tp = mean(salt_radius)

            select case (this%refl_shape)
                case ('random', 'perlin', 'custom')
                    gmax = random(this%nsalt, range=[tp, n2 - tp], seed=safe_seed(int(this%seed, 8)*5), spacing=0.5*tp)
                    hmax = random(this%nsalt, range=[tp, n3 - tp], seed=safe_seed(int(this%seed, 8)*6), spacing=0.5*tp)
                case ('gaussian', 'cauchy')
                    if (this%nsalt > this%ng) then
                        gmax = [mu2, random(this%nsalt - this%ng, range=[tp, n2 - tp], seed=safe_seed(int(this%seed, 8)*5), spacing=0.5*tp)]
                        hmax = [mu3, random(this%nsalt - this%ng, range=[tp, n3 - tp], seed=safe_seed(int(this%seed, 8)*6), spacing=0.5*tp)]
                    else
                        gmax = mu2
                        hmax = mu3
                    end if
            end select

            this%salt = zeros(n1, n2, n3)
            rds = random(this%nsalt, seed=safe_seed(int(this%seed, 8)*15 - 1))
            rds = rescale(rds, salt_radius)
            vds = random(this%nsalt, seed=safe_seed(int(this%seed, 8)*15 - 2))
            vds = rescale(vds, [0.75*this%salt_radius_variation, this%salt_radius_variation])
            pds = rescale(rds, [0.75*this%salt_path_variation, this%salt_path_variation])

            ! Define the top surface of the salt bodies
            pn%n1 = n2
            pn%n2 = n3
            pn%octaves = 4
            pn%seed = safe_seed(int(this%seed, 8)*15 - 3)
            topz = pn%generate()
            topz = rescale(topz, [0.0, this%salt_top_height])

            ! Iterate over all salt bodies
            do isalt = 1, this%nsalt

                nd = nint((1.0 - rand(range=this%salt_top_z, seed=safe_seed(int(this%seed, 8)*14*isalt - 1)))*this%n1)

                ! Salt body boundaries --> salt radius at control surface depths
                qn%n1 = this%salt_nnode
                qn%octaves = 4
                qn%seed = safe_seed(int(this%seed, 8)*15*isalt - 1)
                x1 = qn%generate()

                qn%n1 = this%salt_nnode
                qn%octaves = 4
                qn%seed = safe_seed(int(this%seed, 8)*16*isalt - 1)
                x2 = qn%generate()

                x1 = rescale(x1, range=[1.0 - vds(isalt), 1.0]*rds(isalt))
                x2 = rescale(x2, range=[1.0 - vds(isalt), 1.0]*rds(isalt))

                rc = x2 + x1

                ! Define path deviation curves
                qn%n1 = n1
                qn%octaves = 4
                qn%seed = safe_seed(int(this%seed, 8)*19*isalt - 1)
                x1 = qn%generate()
                x1 = x1 - mean(x1)
                x1 = median_filt(x1, 2)/maxval(x1)*pds(isalt)

                qn%n1 = n1
                qn%octaves = 4
                qn%seed = safe_seed(int(this%seed, 8)*20*isalt - 1)
                x2 = qn%generate()
                x2 = x2 - mean(x2)
                x2 = median_filt(x2, 2)/maxval(x2)*pds(isalt)

                ! Define controlling closed curves
                slice = zeros(360, this%salt_nnode)
                !$omp parallel do private(i)
                do i = 1, this%salt_nnode
                    slice(:, i) = random_circular(360, 0.75*rc(i), 0.25*rc(i), 0.3, 10.0, seed=safe_seed(int(this%seed, 8)*22*isalt - i))
                end do
                !$omp end parallel do
                slice = interp_to(slice, [360, n1], ['', 'pchip'])

                ! Fill enclosed region to get salt body
                !$omp parallel do private(i, j, k, dist, ag, xcenter, ycenter)
                do k = 1, n3
                    do j = 1, n2
                        do i = max(n1 - nd - ceiling(maxval(topz)), 1), n1

                            xcenter = hmax(isalt) + x1(i)
                            ycenter = gmax(isalt) + x2(i)

                            dist = sqrt((j - ycenter)**2 + (k - xcenter)**2)
                            ag = clip(nint(atan2(k - xcenter, j - ycenter + float_tiny)*const_rad2deg) + 181, 1, 360)

                            if (dist <= slice(ag, i) .and. i >= n1 - nd - topz(j, k)) then
                                vp(i, j, k) = this%salt_vp
                                rho(i, j, k) = this%salt_rho
                                if (this%yn_elastic) then
                                    vs(i, j, k) = this%salt_vs
                                end if
                                this%salt(i, j, k) = 1.0
                            end if

                        end do
                    end do

                end do
                !$omp end parallel do

            end do

        end if

        ! Add karst
        if (this%yn_karst) then

            block

                type(karst_3d) :: kg
                real, allocatable, dimension(:, :, :) :: kvol
                integer :: nzk

                ! The karst passages are confined to the depth window karst_z;
                ! generate on a sub-grid spanning [0, karst_z(2)]*n1 with the
                ! caprock fraction covering [0, karst_z(1)]*n1
                nzk = clip(nint(this%karst_z(2)*n1), 2, n1)
                kg%nx = n3
                kg%ny = n2
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

                this%karst = zeros(n1, n2, n3)
                !$omp parallel do private(i, j, k) collapse(3)
                do k = 1, n3
                    do j = 1, n2
                        do i = 1, nzk
                            if (kvol(i, j, k) > 0.5) then
                                vp(i, j, k) = this%karst_vp
                                rho(i, j, k) = this%karst_rho
                                if (this%yn_elastic) then
                                    vs(i, j, k) = this%karst_vs
                                end if
                                this%karst(i, j, k) = 1.0
                                ! Karst carves salt where they overlap
                                if (this%yn_salt) then
                                    this%salt(i, j, k) = 0.0
                                end if
                            end if
                        end do
                    end do
                end do
                !$omp end parallel do

            end block

        end if

        ! Final processing

        if (this%yn_elastic) then
            call this%generate_image_elastic(vp, vs, rho)
        else
            call this%generate_image(vp, rho)
        end if

        ! Output
        ! Fault and fault attributes
        if (this%yn_fault) then
            this%fault = f(ne1 + 1:ne1 + this%n1, ne2 + 1:ne2 + this%n2, ne3 + 1:ne3 + this%n3)
            this%fault_dip = fdip(ne1 + 1:ne1 + this%n1, ne2 + 1:ne2 + this%n2, ne3 + 1:ne3 + this%n3)*const_rad2deg
            this%fault_strike = fstrike(ne1 + 1:ne1 + this%n1, ne2 + 1:ne2 + this%n2, ne3 + 1:ne3 + this%n3)*const_rad2deg
            this%fault_rake = frake(ne1 + 1:ne1 + this%n1, ne2 + 1:ne2 + this%n2, ne3 + 1:ne3 + this%n3)*const_rad2deg
            this%fault_disp = fdisp(ne1 + 1:ne1 + this%n1, ne2 + 1:ne2 + this%n2, ne3 + 1:ne3 + this%n3)
        end if

        ! RGT
        if (this%yn_rgt) then
            this%rgt = t(ne1 + 1:ne1 + this%n1, ne2 + 1:ne2 + this%n2, ne3 + 1:ne3 + this%n3)
            this%rgt = rescale(this%rgt, [0.0, 1.0])
        end if

        ! Facies
        if (this%yn_facies) then
            this%facies = m(ne1 + 1:ne1 + this%n1, ne2 + 1:ne2 + this%n2, ne3 + 1:ne3 + this%n3)
            this%facies = this%facies - minval(this%facies) + 1
            if (this%unconf == 0) then
                this%facies = maxval(this%facies) - this%facies + 1
            end if
        end if

        ! Salt
        if (this%yn_salt) then

            this%salt = this%salt(1:this%n1, :, :)

            if (this%unconf == 0) then

                if (this%yn_fault) then
                    where (this%salt == 1)
                        this%fault = 0
                        this%fault_dip = 0
                        this%fault_strike = 0
                        this%fault_rake = 0
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

            this%karst = this%karst(1:this%n1, :, :)

            if (this%unconf == 0) then

                if (this%yn_fault) then
                    where (this%karst == 1)
                        this%fault = 0
                        this%fault_dip = 0
                        this%fault_strike = 0
                        this%fault_rake = 0
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
        this%vp = vp(1:this%n1, :, :)
        this%rho = rho(1:this%n1, :, :)
        if (this%yn_elastic) then
            this%vs = vs(1:this%n1, :, :)
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
                this%fault_strike = 0
                this%fault_rake = 0
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

    end subroutine generate_3d_geological_model

    !
    !> Generate geological models with one or multiple unconformity surfaces
    !
    subroutine generate_3d_unconformal_geological_model(this)

        type(rgm3_curved), intent(inout) :: this

        type(rgm3_curved), allocatable, dimension(:) :: g
        integer :: iconf, i, j, k
        type(meta_array2_real), allocatable, dimension(:) :: uff
        real, allocatable, dimension(:) :: ufz
        real, allocatable, dimension(:, :, :) :: rgt_above, rgt_below
        real, allocatable, dimension(:, :, :) :: facies_above, facies_below
        real :: tmin, tmax
        type(fractal_noise_2d) :: q
        real, allocatable, dimension(:, :, :) :: vp, vs, rho

        allocate (g(1:this%unconf + 1))

        if (this%yn_elastic) then
            this%vp = zeros(this%n1 + 1, this%n2, this%n3)
            this%vs = zeros(this%n1 + 1, this%n2, this%n3)
            this%rho = zeros(this%n1 + 1, this%n2, this%n3)
        else
            this%vp = zeros(this%n1 + 1, this%n2, this%n3)
            this%rho = zeros(this%n1 + 1, this%n2, this%n3)
        end if

        if (this%yn_fault) then
            this%fault = zeros(this%n1 + 1, this%n2, this%n3)
            this%fault_dip = zeros(this%n1 + 1, this%n2, this%n3)
            this%fault_strike = zeros(this%n1 + 1, this%n2, this%n3)
            this%fault_rake = zeros(this%n1 + 1, this%n2, this%n3)
            this%fault_disp = zeros(this%n1 + 1, this%n2, this%n3)
        end if
        if (this%yn_rgt) then
            this%rgt = zeros(this%n1 + 1, this%n2, this%n3)
        end if
        if (this%yn_facies) then
            this%facies = zeros(this%n1 + 1, this%n2, this%n3)
        end if
        if (this%yn_salt) then
            this%salt = zeros(this%n1 + 1, this%n2, this%n3)
        end if
        if (this%yn_karst) then
            this%karst = zeros(this%n1 + 1, this%n2, this%n3)
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

            call generate_3d_geological_model(g(i))

        end do

        allocate (uff(1:this%unconf + 1))
        ufz = random(this%unconf, range=this%unconf_z, seed=safe_seed(int(this%seed, 8)*31))
        ufz = sort(ufz, order=1)
        do i = 1, this%unconf

            select case (this%unconf_shape)

                case ('meander_channel', 'meander_canyon', 'drainage_channel', 'drainage_canyon')
                    ! The unconformity surface is the erosional topography of a
                    ! geomorphological (channel/canyon/drainage) simulation
                    call unconformity_topography_3d(this, g(i)%seed, i, uff(i)%array)

                case default
                    ! Use Perlin noise to generate unconformity surfaces
                    q%n1 = this%n2
                    q%n2 = this%n3
                    q%octaves = 5
                    q%seed = safe_seed(int(g(i)%seed, 8)*41*i)
                    uff(i)%array = q%generate()

            end select

            if (this%unconf_smooth > 0) then
                uff(i)%array = gauss_filt(uff(i)%array, [this%unconf_smooth, this%unconf_smooth])
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
            this%fault_strike = g(this%unconf + 1)%fault_strike
            this%fault_rake = g(this%unconf + 1)%fault_rake
            this%fault_disp = g(this%unconf + 1)%fault_disp
        end if
        if (this%yn_rgt) then
            this%rgt = g(this%unconf + 1)%rgt
            rgt_above = zeros(this%n1, this%n2, this%n3)
            rgt_below = zeros(this%n1, this%n2, this%n3)
        end if
        if (this%yn_facies) then
            this%facies = g(this%unconf + 1)%facies
            facies_above = zeros(this%n1, this%n2, this%n3)
            facies_below = zeros(this%n1, this%n2, this%n3)
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

            !$omp parallel do private(i, j, k)
            do k = 1, this%n3
                do j = 1, this%n2

                    ! Image by soft merging
                    if (this%yn_elastic) then
                        do i = 1, this%n1
                            if (i < uff(iconf)%array(j, k)) then
                                this%vp(i, j, k) = g(iconf)%vp(i, j, k)
                                this%vs(i, j, k) = g(iconf)%vs(i, j, k)
                                this%rho(i, j, k) = g(iconf)%rho(i, j, k)
                                if (this%yn_salt .and. this%salt_before_unconf) then
                                    this%salt(i, j, k) = 0.0
                                end if
                                if (this%yn_karst .and. this%karst_before_unconf) then
                                    this%karst(i, j, k) = 0.0
                                end if
                            end if
                        end do
                    else
                        do i = 1, this%n1
                            if (i < uff(iconf)%array(j, k)) then
                                this%vp(i, j, k) = g(iconf)%vp(i, j, k)
                                this%rho(i, j, k) = g(iconf)%rho(i, j, k)
                                if (this%yn_salt .and. this%salt_before_unconf) then
                                    this%salt(i, j, k) = 0.0
                                end if
                                if (this%yn_karst .and. this%karst_before_unconf) then
                                    this%karst(i, j, k) = 0.0
                                end if
                            end if
                        end do
                    end if

                    ! Fault by hard merging
                    if (this%yn_fault) then
                        do i = 1, this%n1
                            if (i < uff(iconf)%array(j, k)) then
                                this%fault(i, j, k) = g(iconf)%fault(i, j, k)
                                this%fault_dip(i, j, k) = g(iconf)%fault_dip(i, j, k)
                                this%fault_strike(i, j, k) = g(iconf)%fault_strike(i, j, k)
                                this%fault_rake(i, j, k) = g(iconf)%fault_rake(i, j, k)
                                this%fault_disp(i, j, k) = g(iconf)%fault_disp(i, j, k)
                            end if
                        end do
                    end if

                    ! RGT by history-consistent merging
                    if (this%yn_rgt) then
                        do i = 1, this%n1
                            if (i < uff(iconf)%array(j, k)) then
                                rgt_above(i, j, k) = g(iconf)%rgt(i, j, k)
                            else
                                rgt_below(i, j, k) = this%rgt(i, j, k)
                            end if
                        end do
                    end if

                    ! Facies by hard merging
                    if (this%yn_facies) then
                        do i = 1, this%n1
                            if (i < uff(iconf)%array(j, k)) then
                                facies_above(i, j, k) = g(iconf)%facies(i, j, k)
                            else
                                facies_below(i, j, k) = this%facies(i, j, k)
                            end if
                        end do
                    end if

                end do
            end do
            !$omp end parallel do

            ! Correct RGT
            if (this%yn_rgt) then
                tmax = maxval(rgt_above)
                tmin = minval(rgt_below, mask=(rgt_below /= 0))
                !$omp parallel do private(i, j, k) collapse(3)
                do k = 1, this%n3
                    do j = 1, this%n2
                        do i = 1, this%n1
                            if (rgt_above(i, j, k) /= 0) then
                                rgt_above(i, j, k) = rgt_above(i, j, k) - tmax + tmin
                            end if
                            this%rgt(i, j, k) = rgt_above(i, j, k) + rgt_below(i, j, k)
                        end do
                    end do
                end do
                !$omp end parallel do
            end if

            ! Correct facies
            if (this%yn_facies) then
                tmin = minval(facies_above, mask=(facies_above /= 0))
                tmax = maxval(facies_below)
                !$omp parallel do private(i, j, k) collapse(3)
                do k = 1, this%n3
                    do j = 1, this%n2
                        do i = 1, this%n1
                            if (facies_above(i, j, k) /= 0) then
                                facies_above(i, j, k) = facies_above(i, j, k) - tmin + tmax + 1
                            end if
                            this%facies(i, j, k) = facies_above(i, j, k) + facies_below(i, j, k)
                        end do
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
        this%vp = this%vp(1:this%n1, :, :)
        this%rho = this%rho(1:this%n1, :, :)
        if (this%yn_elastic) then
            this%vs = this%vs(1:this%n1, :, :)
        end if

        if (this%yn_fault) then
            this%fault = this%fault(1:this%n1, :, :)
            this%fault_dip = this%fault_dip(1:this%n1, :, :)
            this%fault_strike = this%fault_strike(1:this%n1, :, :)
            this%fault_rake = this%fault_rake(1:this%n1, :, :)
            this%fault_disp = this%fault_disp(1:this%n1, :, :)
        end if

        if (this%yn_rgt) then
            this%rgt = this%rgt(1:this%n1, :, :)
        end if

        if (this%yn_facies) then
            this%facies = this%facies(1:this%n1, :, :)
        end if

        if (this%yn_salt) then
            this%salt = this%salt(1:this%n1, :, :)
        end if

        if (this%yn_karst) then
            this%karst = this%karst(1:this%n1, :, :)
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
                    this%fault_strike = 0
                    this%fault_rake = 0
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
                    this%fault_strike = 0
                    this%fault_rake = 0
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

    end subroutine generate_3d_unconformal_geological_model

    !
    !> Generate the topography of a channel/canyon/drainage-type unconformity
    !> surface as the erosional depth map of a geomorphological simulation.
    !> The returned surface has shape (n2, n3), is downward-positive, and is
    !> subsequently rescaled by unconf_height like the default random surface.
    !
    subroutine unconformity_topography_3d(this, seed, isurf, surf)

        type(rgm3_curved), intent(in) :: this
        integer, intent(in) :: seed, isurf
        real, allocatable, dimension(:, :), intent(out) :: surf

        type(meandering_channel) :: mc
        type(meandering_canyon) :: my
        type(drainage_channel) :: dc
        type(drainage_canyon) :: dy
        real :: wf
        integer :: nlev

        ! Reference (calibrated) configurations of the meandering channel and
        ! canyon: at length = *_length_ref with *_nbends_ref seed bends, the
        ! migration develops well-formed gooseneck meanders. The mappings
        ! below scale unconf_channel_length relative to these single
        ! calibration points
        real, parameter :: mc_length_ref = 25000.0
        real, parameter :: mc_nbends_ref = 30.0
        real, parameter :: my_length_ref = 15000.0
        real, parameter :: my_nbends_ref = 20.0
        ! Minimum accepted channel length; below this the centerline has too
        ! few nodes for the migration to develop realistic meanders
        real, parameter :: length_min = 10000.0

        ! Width fraction of this surface's channel/canyon
        wf = rand(range=this%unconf_channel_width, seed=safe_seed(int(seed, 8)*91 + isurf))

        select case (this%unconf_shape)

            case ('meander_channel')
                mc%nx = this%n3
                mc%ny = this%n2
                mc%nz = 64
                ! The length/W ratio sets how many meander bends span the
                ! domain. The migration width W is capped at its calibrated
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
                mc%orient = irand(range=[0, 3], seed=safe_seed(int(seed, 8)*71 + isurf))
                mc%yn_depth_only = .true.
                call mc%generate
                surf = mc%depth_map

            case ('meander_canyon')
                my%nx = this%n3
                my%ny = this%n2
                my%nz = 0
                ! Keep length/W large so that the centerline has enough nodes to
                ! survive the migration/cutoff process across all snapshots;
                ! the canyon aperture is controlled by W_canyon instead. As for
                ! the channel, W is capped at its calibrated value so that a
                ! longer canyon adds bends rather than changing the dynamics
                my%length = (my_length_ref/mc_length_ref)*max(this%unconf_channel_length, length_min)
                my%W = 0.03*min(my%length, my_length_ref)
                my%W_canyon = 0.5*wf*my%length
                my%n_bends = max(5, nint(my_nbends_ref*my%length/my_length_ref))
                my%n_iter = nint(clip(4000*this%unconf_channel_sinuosity, 2000.0, 6000.0))
                nlev = max(2, nint(my%n_iter/(my%save_every*1.0)))
                my%terrain_bg = this%unconf_topo*(nlev - 1.0)
                my%seed = safe_seed(int(seed, 8)*61 + isurf)
                my%orient = irand(range=[0, 3], seed=safe_seed(int(seed, 8)*71 + isurf))
                my%yn_depth_only = .true.
                call my%generate
                surf = my%depth_map

            case ('drainage_channel')
                dc%nx = this%n3
                dc%ny = this%n2
                dc%nz = 64
                ! Channels: relatively thin and densely distributed
                dc%W_max = 0.5*wf*0.5*(this%n2 + this%n3)
                dc%D_max = 48.0
                dc%channel_frac = rand(range=this%unconf_channel_density, seed=safe_seed(int(seed, 8)*81 + isurf))
                dc%terrain_bg = this%unconf_topo*(dc%nz - 1.0)
                dc%seed = safe_seed(int(seed, 8)*61 + isurf)
                dc%orient = irand(range=[0, 3], seed=safe_seed(int(seed, 8)*71 + isurf))
                dc%yn_depth_only = .true.
                call dc%generate
                surf = dc%depth_map

            case ('drainage_canyon')
                dy%nx = this%n3
                dy%ny = this%n2
                dy%nz = 64
                ! Canyons: major and few - full width but sparse network
                dy%W_max = wf*0.5*(this%n2 + this%n3)
                dy%channel_frac = 0.25*rand(range=this%unconf_channel_density, seed=safe_seed(int(seed, 8)*81 + isurf))
                dy%terrain_bg = this%unconf_topo*(dy%nz - 1.0)
                dy%seed = safe_seed(int(seed, 8)*61 + isurf)
                dy%orient = irand(range=[0, 3], seed=safe_seed(int(seed, 8)*71 + isurf))
                dy%yn_depth_only = .true.
                call dy%generate
                surf = dy%depth_map

        end select

    end subroutine unconformity_topography_3d

end module geological_model_3d_curved

