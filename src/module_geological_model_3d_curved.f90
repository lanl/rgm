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

module geological_model_3d_curved

    use libflit
    use geological_model_utility

    implicit none

    !
    ! The module is mostly a simplified version of rgm3 and rgm3_elastic,
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
    ! compared with rgm3/rgm3_elastic, but the changes are minimized to keep consistency
    !
    ! Currently, the strike of a fault is still constant; future
    ! version may introduce faults with spatially varying strikes.
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
                    ww = gauss_filt(random(n1, n2, n3, dist=this%noise_type, seed=this%seed*23), this%noise_smooth)
                    ww = ww - mean(ww)
                    ww = ww/maxval(abs(ww))*this%noise_level*maxval(abs(this%image_pp))
                    this%image_pp = this%image_pp + ww

                    ww = gauss_filt(random(n1, n2, n3, dist=this%noise_type, seed=this%seed*23 + 1), this%noise_smooth)
                    ww = ww - mean(ww)
                    ww = ww/maxval(abs(ww))*this%noise_level*maxval(abs(this%image_ps))
                    this%image_ps = this%image_ps + ww

                    ww = gauss_filt(random(n1, n2, n3, dist=this%noise_type, seed=this%seed*23 + 2), this%noise_smooth)
                    ww = ww - mean(ww)
                    ww = ww/maxval(abs(ww))*this%noise_level*maxval(abs(this%image_sp))
                    this%image_sp = this%image_sp + ww

                    ww = gauss_filt(random(n1, n2, n3, dist=this%noise_type, seed=this%seed*23 + 3), this%noise_smooth)
                    ww = ww - mean(ww)
                    ww = ww/maxval(abs(ww))*this%noise_level*maxval(abs(this%image_ss))
                    this%image_ss = this%image_ss + ww

                case ('wavenumber')
                    this%image_pp = this%image_pp + noise_wavenumber(this%image_pp, this%noise_level, this%noise_smooth, this%seed*23)
                    this%image_ps = this%image_ps + noise_wavenumber(this%image_ps, this%noise_level, this%noise_smooth, this%seed*23 + 1)
                    this%image_sp = this%image_sp + noise_wavenumber(this%image_sp, this%noise_level, this%noise_smooth, this%seed*23 + 2)
                    this%image_ss = this%image_ss + noise_wavenumber(this%image_ss, this%noise_level, this%noise_smooth, this%seed*23 + 3)

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
                    ww = gauss_filt(random(n1, n2, n3, dist=this%noise_type, seed=this%seed*23), this%noise_smooth)
                    ww = ww - mean(ww)
                    ww = ww/maxval(abs(ww))*this%noise_level*maxval(abs(this%image_pp))
                    this%image_pp = this%image_pp + ww

                    ww = gauss_filt(random(n1, n2, n3, dist=this%noise_type, seed=this%seed*23 + 1), this%noise_smooth)
                    ww = ww - mean(ww)
                    ww = ww/maxval(abs(ww))*this%noise_level*maxval(abs(this%image_ps))
                    this%image_ps = this%image_ps + ww

                    ww = gauss_filt(random(n1, n2, n3, dist=this%noise_type, seed=this%seed*23 + 2), this%noise_smooth)
                    ww = ww - mean(ww)
                    ww = ww/maxval(abs(ww))*this%noise_level*maxval(abs(this%image_sp))
                    this%image_sp = this%image_sp + ww

                    ww = gauss_filt(random(n1, n2, n3, dist=this%noise_type, seed=this%seed*23 + 3), this%noise_smooth)
                    ww = ww - mean(ww)
                    ww = ww/maxval(abs(ww))*this%noise_level*maxval(abs(this%image_ss))
                    this%image_ss = this%image_ss + ww

                case ('wavenumber')
                    this%image_pp = this%image_pp + noise_wavenumber(this%image_pp, this%noise_level, this%noise_smooth, this%seed*23)
                    this%image_ps = this%image_ps + noise_wavenumber(this%image_ps, this%noise_level, this%noise_smooth, this%seed*23 + 1)
                    this%image_sp = this%image_sp + noise_wavenumber(this%image_sp, this%noise_level, this%noise_smooth, this%seed*23 + 2)
                    this%image_ss = this%image_ss + noise_wavenumber(this%image_ss, this%noise_level, this%noise_smooth, this%seed*23 + 3)

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
                    ww = gauss_filt(random(n1, n2, n3, dist=this%noise_type, seed=this%seed*23), this%noise_smooth)
                    ww = ww - mean(ww)
                    ww = ww/maxval(abs(ww))*this%noise_level*maxval(abs(this%image))
                    this%image = this%image + ww

                case ('wavenumber')
                    this%image = this%image + noise_wavenumber(this%image, this%noise_level, this%noise_smooth, this%seed*23)

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
                    ww = gauss_filt(random(n1, n2, n3, dist=this%noise_type, seed=this%seed*23), this%noise_smooth)
                    ww = ww - mean(ww)
                    ww = ww/maxval(abs(ww))*this%noise_level*maxval(abs(this%image))
                    this%image = this%image + ww

                case ('wavenumber')
                    this%image = this%image + noise_wavenumber(this%image, this%noise_level, this%noise_smooth, this%seed*23)

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
        real, allocatable, dimension(:) :: delta_dip
        real, dimension(1:2) :: pt
        real, allocatable, dimension(:, :) :: dips
        integer :: nsf
        real :: a, b, xs, ys, xys, xys_prev, dxys, b_prev, dist_prev
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
                strike = [random(nint(nf/2.0), range=[0.975, 1.025]*this%strike(1), seed=this%seed*2)*const_deg2rad, &
                    random(nf - nint(nf/2.0), range=[0.975, 1.025]*this%strike(2), seed=this%seed*2)*const_deg2rad]
                rake = [random(nint(nf/2.0), range=[0.95, 1.05]*this%rake(1), seed=this%seed*3)*const_deg2rad, &
                    random(nf - nint(nf/2.0), range=[0.95, 1.05]*this%rake(2), seed=this%seed*3)*const_deg2rad]
                disp = [random(nint(nf/2.0), range=[0.9, 1.1]*this%disp(1), seed=this%seed*4), &
                    random(nf - nint(nf/2.0), range=[0.9, 1.1]*this%disp(2), seed=this%seed*4)]

                strike = clip(strike, 0.0, real(const_pi))
                rake = clip(rake, 0.0, real(const_pi))

            else

                ! Dip, strike, rake angles, and fault displacements
                nsf = size(this%dip)/2
                dip = random(ceiling(nf*1.0/nsf), range=this%dip(1:2), seed=this%seed)
                strike = random(ceiling(nf*1.0/nsf), range=this%strike(1:2), seed=this%seed*2)
                rake = random(ceiling(nf*1.0/nsf), range=this%rake(1:2), seed=this%seed*3)
                disp = random(ceiling(nf*1.0/nsf), range=this%disp(1:2), seed=this%seed*4)
                do i = 2, nsf
                    dip = [dip, random(ceiling(nf*1.0/nsf), range=this%dip((i - 1)*2 + 1:i*2), seed=this%seed)]
                    strike = [strike, random(ceiling(nf*1.0/nsf), range=this%strike((i - 1)*2 + 1:i*2), seed=this%seed*2)]
                    rake = [rake, random(ceiling(nf*1.0/nsf), range=this%rake((i - 1)*2 + 1:i*2), seed=this%seed*2)]
                    disp = [disp, random(ceiling(nf*1.0/nsf), range=this%disp((i - 1)*2 + 1:i*2), seed=this%seed*4)]
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

        sumdisp = disp*(cos(rake)*sin(strike) - sin(rake)*cos(dip)*cos(strike))
        ne2 = ceiling(max(sum(sumdisp, mask=sumdisp > 0), -sum(sumdisp, mask=sumdisp < 0)))
        n2 = n2 + 2*ne2

        sumdisp = disp*(cos(rake)*cos(strike) + sin(rake)*cos(dip)*sin(strike))
        ne3 = ceiling(max(sum(sumdisp, mask=sumdisp > 0), -sum(sumdisp, mask=sumdisp < 0)))
        n3 = n3 + 2*ne3

        ! Compute the top and bottom fault dip angles
        ! Note that to ensure the faults have proper dip angle range within the final cropped image,
        ! here I add some extra degrees to the begin and end dip angles
        dips = zeros(n1, nf)
        delta_dip = random(nf, range=this%delta_dip, seed=nint(this%seed*4.5))*const_deg2rad
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
                r = random(n2, n3, dist='normal', seed=this%seed*5)
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

                mu2 = random(this%ng, range=gmu2, seed=this%seed*5)
                sigma2 = random(this%ng, range=gsigma2, seed=this%seed*6)
                mu3 = random(this%ng, range=gmu3, seed=this%seed*7)
                sigma3 = random(this%ng, range=gsigma3, seed=this%seed*8)
                height = random(this%ng, range=this%refl_height, seed=this%seed*9)
                gtheta = random(this%ng, range=[0.0, 180.0], seed=this%seed*9 - 1)*ifelse(this%rotate_fold, real(const_deg2rad), 0.0)

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
                pn%seed = this%seed*5
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
                    rt = random(n2, n3, dist='normal', seed=this%seed*5 - 1)
                    rt = gauss_filt(rt, [this%refl_smooth_top, this%refl_smooth_top])

                case ('random')
                    rt = random(n2, n3, dist='normal', seed=this%seed*5 - 1)
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

                    mu2 = random(this%ng, range=gmu2, seed=this%seed*5 - 1)
                    sigma2 = random(this%ng, range=gsigma2, seed=this%seed*6 - 1)
                    mu3 = random(this%ng, range=gmu3, seed=this%seed*7 - 1)
                    sigma3 = random(this%ng, range=gsigma3, seed=this%seed*8 - 1)
                    height = random(this%ng, range=this%refl_height, seed=this%seed*9 - 1)
                    gtheta = random(this%ng, range=[0.0, 180.0], seed=this%seed*9 - 2)*ifelse(this%rotate_fold, real(const_deg2rad), 0.0)

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
                    pn%seed = this%seed*5 - 1
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

        vv = random(nl - 1, seed=this%seed*6)*this%delta_v
        vv = linspace(this%vmax*(1.0 + ne1*1.0/this%n1), this%vmin*(1.0 - ne1*1.0/this%n1), nl - 1) + vv

        lz = zeros(nl, n2, n3)
        plw = random(nl, range=[1 - this%lwv, 1 + this%lwv], seed=this%seed*7)

        !$omp parallel do private(j, k)
        do k = 1, n3
            do j = 1, n2
                lz(:, j, k) = ginterp([0.0, (ne1 + 1.0)/n1, (n1 - ne1 - 1.0)/n1, 1.0], &
                    [r(j, k), r(j, k) + ne1, rt(j, k) - ne1, rt(j, k)], &
                    linspace(0.0, 1.0, nl))
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
                pxy(i, :, :) = gauss_filt(random(n2, n3, seed=this%seed*7 - i), [1, 1]*this%secondary_refl_smooth)
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
                ! As such, it may be slightly different from what obtained with rgm2
                if (this%yn_rgt) then
                    t(:, j, k) = ginterp(n1 - 1.0 - lz(:, j, k), linspace(1.0, 0.0, nl), linspace(0.0, n1 - 1.0, n1))
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

                rc = random(nf - 1, range=[0.75, 1.25]*0.9*this%n2/(nf - 1.0), seed=this%seed*12)
                rc = rc*0.9*this%n2/sum(rc)
                rc = sort(rc)
                f2 = zeros(nf)
                f2(1) = ne2 + rand(range=[0.2, 0.8]*(this%n2 - sum(rc)), seed=this%seed*12 - 1)
                f2(2:) = f2(1) + cumsum(rc)
                if (.not. this%yn_group_faults) then
                    f2 = random_permute(f2, seed=this%seed*12 - 2)
                end if

                rc = random(nf - 1, range=[0.75, 1.25]*0.9*this%n3/(nf - 1.0), seed=this%seed*13)
                rc = rc*0.9*this%n3/sum(rc)
                rc = sort(rc)
                f3 = zeros(nf)
                f3(1) = ne3 + rand(range=[0.2, 0.8]*(this%n3 - sum(rc)), seed=this%seed*13 - 1)
                f3(2:) = f3(1) + cumsum(rc)
                if (.not. this%yn_group_faults) then
                    f3 = random_permute(f3, seed=this%seed*13 - 2)
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
                    seed=this%seed*16, spacing=0.75*(n2 - 2*ne2 - 0.2*this%n2)/nf)
                f3 = random(nf, range=[ne3 + 0.1*this%n3, n3 - ne3 - 0.1*this%n3], &
                    seed=this%seed*17, spacing=0.75*(n3 - 2*ne3 - 0.2*this%n3)/nf)

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
                ffstrike = fstrike
                ffrake = frake
                ffdisp = fdisp
                cf = 0

                xys_prev = 0
                xys = 0
                b_prev = 0
                fblock = falses(n1, n2, n3)

                do i = 1, n1

                    theta = dips(i, fi)

                    !! The deviation of the curved fault w.r.t. the straight fault
                    ! xs = -(1.0/tan(theta) - 1.0/tan(dip(fi)))*sin(strike(fi))*(i - 1.0)
                    ! ys = +(1.0/tan(theta) - 1.0/tan(dip(fi)))*cos(strike(fi))*(i - 1.0)

                    ! The location of the curved fault
                    xs = -1.0/tan(theta)*sin(strike(fi))*(i - 1.0)
                    ys = +1.0/tan(theta)*cos(strike(fi))*(i - 1.0)

                    xys = sqrt(xs**2 + ys**2)
                    if (i == 1) then
                        xys_prev = xys
                    end if
                    dxys = xys - xys_prev

                    if (abs(strike(fi) - const_pi_half) >= const_pi/4.0) then

                        a = tan(strike(fi))
                        b = (f2(fi) + ys) - a*(f3(fi) + xs)
                        if (i == 1) then
                            b_prev = b
                        end if

                        !$omp parallel do private(j, k, x0, y0, dist) collapse(2)
                        do k = 1, n3
                            do j = 1, n2

                                x0 = k - 1.0
                                y0 = j - 1.0

                                dist = abs(y0 - (a*x0 + b))/sqrt(1.0 + a**2)
                                if (dist < 0.5*fwidth/sin(theta)) then
                                    f(i, j, k) = fi
                                    cf(i, j, k) = 1.0
                                    fdip(i, j, k) = theta
                                    fstrike(i, j, k) = strike(fi)
                                    frake(i, j, k) = rake(fi)
                                end if

                            end do
                        end do
                        !$omp end parallel do

                        if (abs(dxys) >= sqrt(2.0)) then

                            !$omp parallel do private(j, k, x0, y0, dist_prev, dist) collapse(2)
                            do k = 1, n3
                                do j = 1, n2

                                    x0 = k - 1.0
                                    y0 = j - 1.0

                                    dist_prev = abs(y0 - (a*x0 + b_prev))/sqrt(1.0 + a**2)
                                    dist = abs(y0 - (a*x0 + b))/sqrt(1.0 + a**2)
                                    if (dist_prev <= sqrt(2.0)*abs(dxys) .and. dist <= sqrt(2.0)*abs(dxys)) then
                                        f(i, j, k) = fi
                                        cf(i, j, k) = 1.0
                                        fdip(i, j, k) = theta
                                        fstrike(i, j, k) = strike(fi)
                                        frake(i, j, k) = rake(fi)
                                    end if
                                end do
                            end do
                            !$omp end parallel do

                        end if

                        !$omp parallel do private(j, k, x0, y0) collapse(2)
                        do k = 1, n3
                            do j = 1, n2
                                x0 = k - 1.0
                                y0 = j - 1.0
                                if (theta < const_pi_half) then
                                    if (y0 - (a*x0 + b) < 0) then
                                        fblock(i, j, k) = .true.
                                    end if
                                else
                                    if (y0 - (a*x0 + b) > 0) then
                                        fblock(i, j, k) = .true.
                                    end if
                                end if
                            end do
                        end do
                        !$omp end parallel do

                    else

                        a = tan(const_pi_half - strike(fi))
                        b = (f3(fi) + xs) - a*(f2(fi) + ys)
                        if (i == 1) then
                            b_prev = b
                        end if

                        !$omp parallel do private(j, k, x0, y0, dist) collapse(2)
                        do k = 1, n3
                            do j = 1, n2

                                x0 = k - 1.0
                                y0 = j - 1.0

                                dist = abs(x0 - (a*y0 + b))/sqrt(1.0 + a**2)
                                if (dist < 0.5*fwidth/sin(theta)) then
                                    f(i, j, k) = fi
                                    cf(i, j, k) = 1.0
                                    fdip(i, j, k) = theta
                                    fstrike(i, j, k) = strike(fi)
                                    frake(i, j, k) = rake(fi)
                                end if

                            end do
                        end do
                        !$omp end parallel do

                        if (abs(dxys) >= sqrt(2.0)) then

                            !$omp parallel do private(j, k, x0, y0, dist_prev, dist) collapse(2)
                            do k = 1, n3
                                do j = 1, n2

                                    x0 = k - 1.0
                                    y0 = j - 1.0

                                    dist_prev = abs(x0 - (a*y0 + b_prev))/sqrt(1.0 + a**2)
                                    dist = abs(x0 - (a*y0 + b))/sqrt(1.0 + a**2)
                                    if (dist_prev <= sqrt(2.0)*abs(dxys) .and. dist <= sqrt(2.0)*abs(dxys)) then
                                        f(i, j, k) = fi
                                        cf(i, j, k) = 1.0
                                        fdip(i, j, k) = theta
                                        fstrike(i, j, k) = strike(fi)
                                        frake(i, j, k) = rake(fi)
                                    end if
                                end do
                            end do
                            !$omp end parallel do

                        end if

                        !$omp parallel do private(j, k, x0, y0) collapse(2)
                        do k = 1, n3
                            do j = 1, n2
                                x0 = k - 1.0
                                y0 = j - 1.0
                                if (theta < const_pi_half) then
                                    if (x0 - (a*y0 + b_prev) <= 0) then
                                        fblock(i, j, k) = .true.
                                    end if
                                else
                                    if (x0 - (a*y0 + b_prev) >= 0) then
                                        fblock(i, j, k) = .true.
                                    end if
                                end if
                            end do
                        end do
                        !$omp end parallel do

                    end if

                    ! next iteration
                    xys_prev = xys
                    b_prev = b

                end do

                ! Shift block
                !$omp parallel do private(i, j, k, newi, newj, newk) collapse(3)
                do k = 1, n3
                    do j = 1, n2
                        do i = 1, n1

                            newi = nint(i + (-sin(rake(fi))*sin(dip(fi)))*disp(fi))
                            newj = nint(j + (cos(rake(fi))*sin(strike(fi)) - sin(rake(fi))*cos(dip(fi))*cos(strike(fi)))*disp(fi))
                            newk = nint(k + (cos(rake(fi))*cos(strike(fi)) + sin(rake(fi))*cos(dip(fi))*sin(strike(fi)))*disp(fi))

                            if (newi >= 1 .and. newi <= n1 .and. newj >= 1 .and. newj <= n2 .and. newk >= 1 .and. newk <= n3) then
                                if (fblock(newi, newj, newk)) then

                                    ! Move velocity
                                    w(newi, newj, newk) = ww(i, j, k)

                                    ! Move facies
                                    if (this%yn_facies) then
                                        m(newi, newj, newk) = mm(i, j, k)
                                    end if

                                    ! Move RGT
                                    if (this%yn_rgt) then
                                        t(newi, newj, newk) = tt(i, j, k)
                                    end if

                                    ! Move existing faults
                                    if (cf(newi, newj, newk) == 0) then
                                        f(newi, newj, newk) = ff(i, j, k)
                                        fdip(newi, newj, newk) = ffdip(i, j, k)
                                        fstrike(newi, newj, newk) = ffstrike(i, j, k)
                                        frake(newi, newj, newk) = ffrake(i, j, k)
                                        fdisp(newi, newj, newk) = ffdisp(i, j, k)
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
                    gmax = random(this%nsalt, range=[tp, n2 - tp], seed=this%seed*5, spacing=0.5*tp)
                    hmax = random(this%nsalt, range=[tp, n3 - tp], seed=this%seed*6, spacing=0.5*tp)
                case ('gaussian', 'cauchy')
                    if (this%nsalt > this%ng) then
                        gmax = [mu2, random(this%nsalt - this%ng, range=[tp, n2 - tp], seed=this%seed*5, spacing=0.5*tp)]
                        hmax = [mu3, random(this%nsalt - this%ng, range=[tp, n3 - tp], seed=this%seed*6, spacing=0.5*tp)]
                    else
                        gmax = mu2
                        hmax = mu3
                    end if
            end select

            this%salt = zeros(n1, n2, n3)
            rds = random(this%nsalt, seed=this%seed*15 - 1)
            rds = rescale(rds, salt_radius)
            vds = random(this%nsalt, seed=this%seed*15 - 2)
            vds = rescale(vds, [0.75*this%salt_radius_variation, this%salt_radius_variation])
            pds = rescale(rds, [0.75*this%salt_path_variation, this%salt_path_variation])

            ! Define the top surface of the salt bodies
            pn%n1 = n2
            pn%n2 = n3
            pn%octaves = 4
            pn%seed = this%seed*15 - 3
            topz = pn%generate()
            topz = rescale(topz, [0.0, this%salt_top_height])

            ! Iterate over all salt bodies
            do isalt = 1, this%nsalt

                nd = nint((1.0 - rand(range=this%salt_top_z, seed=this%seed*14*isalt - 1))*this%n1)

                ! Salt body boundaries --> salt radius at control surface depths
                qn%n1 = this%salt_nnode
                qn%octaves = 4
                qn%seed = this%seed*15*isalt - 1
                x1 = qn%generate()

                qn%n1 = this%salt_nnode
                qn%octaves = 4
                qn%seed = this%seed*16*isalt - 1
                x2 = qn%generate()

                x1 = rescale(x1, range=[1.0 - vds(isalt), 1.0]*rds(isalt))
                x2 = rescale(x2, range=[1.0 - vds(isalt), 1.0]*rds(isalt))

                rc = x2 + x1

                ! Define path deviation curves
                qn%n1 = n1
                qn%octaves = 4
                qn%seed = this%seed*19*isalt - 1
                x1 = qn%generate()
                x1 = x1 - mean(x1)
                x1 = median_filt(x1, 2)/maxval(x1)*pds(isalt)

                qn%n1 = n1
                qn%octaves = 4
                qn%seed = this%seed*20*isalt - 1
                x2 = qn%generate()
                x2 = x2 - mean(x2)
                x2 = median_filt(x2, 2)/maxval(x2)*pds(isalt)

                ! Define controlling closed curves
                slice = zeros(360, this%salt_nnode)
                !$omp parallel do private(i)
                do i = 1, this%salt_nnode
                    slice(:, i) = random_circular(360, 0.75*rc(i), 0.25*rc(i), 0.3, 10.0, seed=this%seed*22*isalt - i)
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

            call generate_3d_geological_model(g(i))

        end do

        allocate (uff(1:this%unconf + 1))
        ufz = random(this%unconf, range=this%unconf_z, seed=this%seed*31)
        ufz = sort(ufz, order=1)
        do i = 1, this%unconf

            ! Use Perlin noise to generate unconformity surfaces
            q%n1 = this%n2
            q%n2 = this%n3
            q%octaves = 5
            q%seed = g(i)%seed*41*i
            uff(i)%array = q%generate()
            if (this%unconf_smooth > 0) then
                uff(i)%array = gauss_filt(uff(i)%array, [this%unconf_smooth, this%unconf_smooth])
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

end module geological_model_3d_curved

