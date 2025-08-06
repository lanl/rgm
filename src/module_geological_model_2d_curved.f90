!
! © 2024-2025. Triad National Security, LLC. All rights reserved.
!
! This program was produced under U.S. Government contract 89233218CNA000001
! for Los Alamos National Laboratory (LANL), which is operated by
! Triad National Security, LLC for the U.S. Department of Energy/National Nuclear
! Security Administration. All rights in the program are reserved by
! Triad National Security, LLC, and the U.S. Department of Energy/National
! Nuclear Security Administration. The Government is granted for itself and
! others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
! license in this material to reproduce, prepare. derivative works,
! distribute copies to the public, perform publicly and display publicly,
! and to permit others to do so.
!
! Author:
!    Kai Gao, kaigao@lanl.gov
!

module geological_model_2d_curved

    !
    ! The module is mostly a simplified version of rgm2 and rgm2_elastic,
    ! with a few changes:
    !   - The way that generates layers and reflectors changes from
    !       generate reflectivity series -> add faults -> convolve with wavelet
    !       generate velocity model -> add faults -> compute reflectivity -> convolve with wavelet
    !       The process is closer to realistic geological sedimentation and deformation.
    !       Note that in this case, the image of faults may not be zeros, especially for faults with
    !       moderate angles.
    !   - Support generating curved faults to further enhance reality where
    !       the fault dip at the top can be different from that of the bottom
    !   - Simplifies salt body insertion and reflectivity computation,
    !       so the reflectivity image of salt may be more accurate now.
    !   - Support generating both acoustic and elastic velocity model and migration images.
    ! As such, there are some changes on the parameters of rgm2_curved
    ! compared with rgm2/rgm2_elastic, but the changes are minimized to keep consistency
    !

    use libflit
    use geological_model_utility

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
                    ww = gauss_filt(random(n1, n2, dist=this%noise_type, seed=this%seed*17), this%noise_smooth)
                    ww = ww - mean(ww)
                    ww = ww/maxval(abs(ww))*this%noise_level*maxval(abs(this%image_pp))
                    this%image_pp = this%image_pp + ww

                    ww = gauss_filt(random(n1, n2, dist=this%noise_type, seed=this%seed*17 + 1), this%noise_smooth)
                    ww = ww - mean(ww)
                    ww = ww/maxval(abs(ww))*this%noise_level*maxval(abs(this%image_ps))
                    this%image_ps = this%image_ps + ww

                    ww = gauss_filt(random(n1, n2, dist=this%noise_type, seed=this%seed*17 + 2), this%noise_smooth)
                    ww = ww - mean(ww)
                    ww = ww/maxval(abs(ww))*this%noise_level*maxval(abs(this%image_sp))
                    this%image_sp = this%image_sp + ww

                    ww = gauss_filt(random(n1, n2, dist=this%noise_type, seed=this%seed*17 + 3), this%noise_smooth)
                    ww = ww - mean(ww)
                    ww = ww/maxval(abs(ww))*this%noise_level*maxval(abs(this%image_ss))
                    this%image_ss = this%image_ss + ww

                case ('wavenumber')
                    this%image_pp = this%image_pp + noise_wavenumber(this%image_pp, this%noise_level, this%noise_smooth, this%seed*17)
                    this%image_ps = this%image_ps + noise_wavenumber(this%image_ps, this%noise_level, this%noise_smooth, this%seed*17 + 1)
                    this%image_sp = this%image_sp + noise_wavenumber(this%image_sp, this%noise_level, this%noise_smooth, this%seed*17 + 2)
                    this%image_ss = this%image_ss + noise_wavenumber(this%image_ss, this%noise_level, this%noise_smooth, this%seed*17 + 3)

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
                    ww = gauss_filt(random(n1, n2, dist=this%noise_type, seed=this%seed*17), this%noise_smooth)
                    ww = ww - mean(ww)
                    ww = ww/maxval(abs(ww))*this%noise_level*maxval(abs(this%image_pp))
                    this%image_pp = this%image_pp + ww

                    ww = gauss_filt(random(n1, n2, dist=this%noise_type, seed=this%seed*17 + 1), this%noise_smooth)
                    ww = ww - mean(ww)
                    ww = ww/maxval(abs(ww))*this%noise_level*maxval(abs(this%image_ps))
                    this%image_ps = this%image_ps + ww

                    ww = gauss_filt(random(n1, n2, dist=this%noise_type, seed=this%seed*17 + 2), this%noise_smooth)
                    ww = ww - mean(ww)
                    ww = ww/maxval(abs(ww))*this%noise_level*maxval(abs(this%image_sp))
                    this%image_sp = this%image_sp + ww

                    ww = gauss_filt(random(n1, n2, dist=this%noise_type, seed=this%seed*17 + 3), this%noise_smooth)
                    ww = ww - mean(ww)
                    ww = ww/maxval(abs(ww))*this%noise_level*maxval(abs(this%image_ss))
                    this%image_ss = this%image_ss + ww

                case ('wavenumber')
                    this%image_pp = this%image_pp + noise_wavenumber(this%image_pp, this%noise_level, this%noise_smooth, this%seed*17)
                    this%image_ps = this%image_ps + noise_wavenumber(this%image_ps, this%noise_level, this%noise_smooth, this%seed*17 - 1)
                    this%image_sp = this%image_sp + noise_wavenumber(this%image_sp, this%noise_level, this%noise_smooth, this%seed*17 - 2)
                    this%image_ss = this%image_ss + noise_wavenumber(this%image_ss, this%noise_level, this%noise_smooth, this%seed*17 - 3)

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
                    ww = gauss_filt(random(n1, n2, dist=this%noise_type, seed=this%seed*17), this%noise_smooth)
                    ww = ww - mean(ww)
                    ww = ww/maxval(abs(ww))*this%noise_level*maxval(abs(this%image))
                    this%image = this%image + ww

                case ('wavenumber')
                    this%image = this%image + noise_wavenumber(this%image, this%noise_level, this%noise_smooth, this%seed*17)

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
                    ww = gauss_filt(random(n1, n2, dist=this%noise_type, seed=this%seed*17), this%noise_smooth)
                    ww = ww - mean(ww)
                    ww = ww/maxval(abs(ww))*this%noise_level*maxval(abs(this%image))
                    this%image = this%image + ww

                case ('wavenumber')
                    this%image = this%image + noise_wavenumber(this%image, this%noise_level, this%noise_smooth, this%seed*17)

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
                disp = [random(nint(nf/2.0), range=[0.9, 1.1]*this%disp(1), seed=this%seed*2), &
                    random(nf - nint(nf/2.0), range=[0.9, 1.1]*this%disp(2), seed=this%seed*2)]

            else

                ! Dip angles and fault displacements
                nsf = size(this%dip)/2
                dip = random(ceiling(nf*1.0/nsf), range=this%dip(1:2), seed=this%seed)
                disp = random(ceiling(nf*1.0/nsf), range=this%disp(1:2), seed=this%seed*2)
                do i = 2, nsf
                    dip = [dip, random(ceiling(nf*1.0/nsf), range=this%dip((i - 1)*2 + 1:i*2), seed=this%seed)]
                    disp = [disp, random(ceiling(nf*1.0/nsf), range=this%disp((i - 1)*2 + 1:i*2), seed=this%seed*2)]
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
        delta_dip = random(nf, range=this%delta_dip, seed=nint(this%seed*2.5))*const_deg2rad
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
                r = random(n2, dist='normal', seed=this%seed*3)
                r = gauss_filt(r, this%refl_smooth)

            case ('random')
                r = random(n2, dist='normal', seed=this%seed*3)
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

                mu = random(this%ng, range=mu2, seed=this%seed*3)
                sigma = random(this%ng, range=sigma2, seed=this%seed*4)
                height = random(this%ng, range=this%refl_height, seed=this%seed*5)

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
                pn%seed = this%seed*3
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
                    rt = random(n2, dist='normal', seed=this%seed*3 - 1)
                    rt = gauss_filt(rt, this%refl_smooth_top)

                case ('random')
                    rt = random(n2, dist='normal', seed=this%seed*3 - 1)
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

                    mu = random(this%ng, range=mu2, seed=this%seed*3 - 1)
                    sigma = random(this%ng, range=sigma2, seed=this%seed*4 - 1)
                    height = random(this%ng, range=this%refl_height_top, seed=this%seed*5 - 1)

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
                    pn%seed = this%seed*3 - 1
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

        vv = random(nl - 1, seed=this%seed*6)*this%delta_v
        vv = linspace(this%vmax*(1.0 + ne1*1.0/this%n1), this%vmin*(1.0 - ne1*1.0/this%n1), nl - 1) + vv

        lz = zeros(nl, n2)
        call assert(abs(this%lwv) <= 1, ' Error: lwv must be 0 <= lwv <= 1. ')
        plw = random(nl, range=[1 - this%lwv, 1 + this%lwv], seed=this%seed*7)

        !$omp parallel do private(j)
        do j = 1, n2
            lz(:, j) = ginterp([0.0, (ne1 + 1.0)/n1, (n1 - ne1 - 1.0)/n1, 1.0], &
                [r(j), r(j) + ne1, rt(j) - ne1, rt(j)], &
                linspace(0.0, 1.0, nl))
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
                pxy(i, :) = gauss_filt(random(n2, seed=this%seed*7 - i), this%secondary_refl_smooth)
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
                t(:, j) = ginterp(n1 - 1.0 - lz(:, j), linspace(1.0, 0.0, nl), linspace(0.0, n1 - 1.0, n1))
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

                rc = random(nf - 1, range=[0.75, 1.25]*0.9*this%n2/(nf - 1.0), seed=this%seed*12)
                rc = rc*0.9*this%n2/sum(rc)
                rc = sort(rc)
                f2 = zeros(nf)
                f2(1) = ne2 + rand(range=[0.05, 0.1]*this%n2, seed=this%seed*12 - 1)
                f2(2:) = f2(1) + cumsum(rc)
                if (.not. this%yn_group_faults) then
                    f2 = random_permute(f2, seed=this%seed*12 - 2)
                end if

            else

                f2 = random(nf, range=[ne2 + 0.1*this%n2, n2 - ne2 - 0.1*this%n2], &
                    seed=this%seed*12, spacing=0.75*(n2 - 2*ne2 - 0.2*this%n2)/nf)

            end if

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

                    !$omp parallel do private(j, x0, dist)
                    do j = 1, n2

                        x0 = j - 1.0

                        dist = abs(x0 - b)
                        if (dist < 0.5*fwidth/sin(theta)) then
                            f(i, j) = fi
                            cf(i, j) = 1.0
                            fdip(i, j) = theta
                            fdisp(i, j) = 1.0
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
                                fdisp(i, j) = 1.0
                            end if
                        end do
                        !$omp end parallel do

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

                ! Shift block
                ! Here I use constant throw in both x and z, rather than depth-varying throw
                ! Modeling depth-varying throw seems very difficult in a regular-grid setting
                ! Another reason is that if using depth-varying throw, particularly in x,
                ! then the dip of previous faults will have to change accordingly since
                ! the points of faults have shifted nonuniformly
                !$omp parallel do private(i, j, newi, newj) collapse(2)
                do j = 1, n2
                    do i = 1, n1

                        newi = nint(i + disp(fi)*sin(dip(fi)))
                        newj = nint(j - disp(fi)*cos(dip(fi)))

                        if (newi >= 1 .and. newi <= n1 .and. newj >= 1 .and. newj <= n2) then
                            if (fblock(newi, newj)) then

                                ! Move velocity model
                                w(newi, newj) = ww(i, j)

                                ! Move facies
                                if (this%yn_facies) then
                                    m(newi, newj) = mm(i, j)
                                end if

                                ! Move RGT
                                if (this%yn_rgt) then
                                    t(newi, newj) = tt(i, j)
                                end if

                                ! Move existing faults
                                if (cf(newi, newj) == 0) then
                                    f(newi, newj) = ff(i, j)
                                    fdip(newi, newj) = ffdip(i, j)
                                    fdisp(newi, newj) = ffdisp(i, j)
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
                    gmax = random(this%nsalt, range=[tp, n2 - tp], seed=this%seed*5, spacing=0.5*tp)
                case ('gaussian', 'cauchy')
                    if (this%nsalt > this%ng) then
                        gmax = [mu, random(this%nsalt - this%ng, range=[tp, n2 - tp], seed=this%seed*5, spacing=0.5*tp)]
                    else
                        gmax = mu
                    end if
            end select

            this%salt = zeros(n1, n2)
            rds = random(this%nsalt, seed=this%seed*14 - 1)
            rds = rescale(rds, salt_radius)
            vds = random(2*this%nsalt, seed=this%seed*15 - 1)
            vds = rescale(vds, [0.75*this%salt_radius_variation, this%salt_radius_variation])
            pds = rescale(rds, [0.75*this%salt_path_variation, this%salt_path_variation])

            ! Build top interface of the salt
            pn%n1 = n2
            pn%octaves = 4
            pn%seed = this%seed*14 - 2
            topz = pn%generate()
            topz = rescale(topz, [0.0, this%salt_top_height])

            ! Insert salt bodies
            do isalt = 1, this%nsalt

                nd = nint((1.0 - rand(range=this%salt_top_z, seed=this%seed*15*isalt - 2))*this%n1)

                ! Salt body boundaries
                pn%n1 = this%salt_nnode
                pn%octaves = 4
                pn%seed = this%seed*15*isalt - 1
                x1 = pn%generate()

                pn%n1 = this%salt_nnode
                pn%octaves = 4
                pn%seed = this%seed*16*isalt - 1
                x2 = pn%generate()

                x1 = interp_to(x1, n1, 'pchip')
                x2 = interp_to(x2, n1, 'pchip')

                x1 = rescale(x1, range=[1.0 - vds(isalt), 1.0]*rds(isalt))
                x2 = rescale(x2, range=[1.0 - vds(isalt), 1.0]*rds(isalt))

                ! Salt body path deviations
                pn%n1 = n1
                pn%octaves = 4
                pn%seed = this%seed*17*isalt - 1
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

        do i = 1, this%unconf + 1

            g(i) = this
            g(i)%n1 = this%n1 + 1

            g(i)%yn_salt = .false.
            if (this%yn_salt .and. i == this%unconf + 1) then
                g(i)%yn_salt = .true.
            end if

            if (this%seed /= -1) then
                g(i)%seed = g(i)%seed*i
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
        ufz = random(this%unconf, range=this%unconf_z, seed=this%seed*31)
        ufz = sort(ufz, order=1)
        do i = 1, this%unconf

            ! Use Perlin noise to generate unconformity surfaces
            q%n1 = this%n2
            q%octaves = 5
            q%seed = g(i)%seed*41*i
            uff(i)%array = q%generate()
            if (this%unconf_smooth > 0) then
                uff(i)%array = gauss_filt(uff(i)%array, this%unconf_smooth)
            end if

            uff(i)%array = rescale(uff(i)%array, &
                [0.0, rand(range=this%unconf_height, seed=g(i)%seed*51*i)]) + ufz(i)*this%n1

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

        where (this%salt == 1)
            this%vp = this%salt_vp
            this%rho = this%salt_rho
        end where
        if (this%yn_elastic) then
            where (this%salt == 1)
                this%vs = this%salt_vs
            end where
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

end module geological_model_2d_curved

