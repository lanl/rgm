!
! © 2024. Triad National Security, LLC. All rights reserved.
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
    ! The module is mostly a simplified version of geological_model_3d
    ! with a few changes:
    !   - The way that generates layers and reflectors changes from
    !       generate reflectivity series -> add faults -> convolve with wavelet
    !       generate velocity model -> add faults -> compute reflectivity -> convolve with wavelet
    !       The process is closer to realistic geological sedimentary and deformation process.
    !       Note that in this case, the image of faults may not be zeros, especially for faults with
    !       moderate angles.
    !   - Allows for curved faults to further enhance reality where
    !       the fault dip at the top can be different from that of the bottom
    !   - Simplifies salt body insertion and reflectivity computation,
    !       so the reflectivity image of salt may be more accurate now.
    ! As such, there are some changes on the parameters of rgm3_curved
    ! compared with rgm3, but the changes are minimized to keep consistency
    !
    ! Currently, the strike of a fault is still constant; future
    ! version may introduce varying-strike faults.
    !

    ! 3D random geological model -- acoustic
    type rgm3_curved
        !> Meanings of the parameters are similar with those in the 2D case

        !==============================================================================================
        integer :: n1 = 128
        integer :: n2 = 128
        integer :: n3 = 128
        integer :: nf = 4
        integer :: nl = 20
        integer :: seed = -1
        real :: refl_smooth = 20.0
        real :: refl_smooth_top = 40.0
        real :: f0 = 150.0
        real :: fwidth = 2.0
        real, allocatable, dimension(:) :: dip, strike, rake, disp
        real, dimension(1:2) :: refl_slope = [0.0, 0.0]
        real, dimension(1:2) :: refl_slope_top = [0.0, 0.0]
        real, dimension(1:3) :: noise_smooth = [1.0, 1.0, 1.0]
        real :: noise_level = 0.25
        character(len=24) :: wave = 'ricker'
        character(len=24) :: refl_shape = 'random'
        real, allocatable, dimension(:, :) :: refl
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
        real :: lwv = 0.0

        logical :: yn_rgt = .false.
        logical :: yn_facies = .false.
        logical :: yn_fault = .true.
        real, allocatable, dimension(:, :, :) :: image, rgt, facies, fault
        real, allocatable, dimension(:, :, :) :: fault_dip, fault_strike, fault_rake, fault_disp
        real, dimension(1:3) :: psf_sigma = [5.0, 2.5, 2.5]
        real, allocatable, dimension(:, :, :) :: psf
        logical :: custom_psf = .false.
        real :: facies_threshold = 0.0
        !        character(len=12) :: refl_amp_dist = 'normal'
        character(len=12) :: noise_type = 'normal'
        logical :: yn_conv_noise = .false.
        real :: secondary_refl_height_ratio = 0.0
        logical :: yn_regular_fault = .false.
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
        real, dimension(1:2) :: unconf_amp = [0.05, 0.15]
        integer :: unconf_nl = 99999

        !==============================================================================================
        logical :: yn_salt = .false.
        integer :: nsalt = 1
        integer :: nstem = 5
        real, dimension(1:2) :: salt_max_radius = [0.05, 0.15]
        real, dimension(1:2) :: salt_top_max = [0.2, 0.4]
        real :: salt_top_smooth = 10.0
        character(len=24) :: salt_type = 'dome'
        real :: salt_vel = 5000.0
        real, allocatable, dimension(:, :, :) :: salt
        real :: salt_top_amp = 20.0
        real :: salt_noise = 0.0
        real :: salt_amp = 2.0

    contains
        procedure :: generate => generate_3d
    end type rgm3_curved

    private
    public :: rgm3_curved

contains

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
        real, allocatable, dimension(:, :, :) :: ww, ff, cf, tt, mm, lz
        integer :: nf, nl, n1, n2, n3
        real, allocatable, dimension(:) :: disp, dip, strike, rake, f2, f3, vv
        integer :: newi, newj, newk
        real, allocatable, dimension(:, :) :: r, rt, sr
        integer :: i, j, k, fi, ne1, ne2, ne3, l
        real, allocatable, dimension(:) :: plw
        real, allocatable, dimension(:, :, :) :: fdip, fstrike, frake, fdisp
        real, allocatable, dimension(:, :, :) :: ffdip, ffstrike, ffrake, ffdisp
        real, allocatable, dimension(:) :: wavelet
        real :: f0, x0, y0, fwidth, dt, wt
        real, allocatable, dimension(:) :: mu2, sigma2, mu3, sigma3, height
        real, dimension(1:2) :: gmu2, gmu3, gsigma2, gsigma3
        real, allocatable, dimension(:, :, :) :: psf
        real, allocatable, dimension(:) :: psf1, psf2, psf3
        real, allocatable, dimension(:) :: sumdisp, rc
        real :: theta
        real, allocatable, dimension(:) :: delta_dip
        real, dimension(1:2) :: pt
        real, allocatable, dimension(:, :) :: dips
        integer :: nsf
        real :: a, b, xs, ys, xys, xys_prev, dxys, b_prev, dist_prev
        logical, allocatable, dimension(:, :, :) :: fblock
        real, allocatable, dimension(:) :: x1, x2, x3, x4, rds
        integer :: nd, isalt
        real, allocatable, dimension(:, :) :: slice, topz
        real :: dist, tp
        integer :: ag, nex
        integer, allocatable, dimension(:) :: gmax, hmax
        type(fractal_noise_2d) :: pn

        fwidth = this%fwidth
        dt = 0.001

        n1 = this%n1
        n2 = this%n2
        n3 = this%n3
        nf = this%nf

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

        ! The following value may be too big (therefore too slow) for 3D
        ! ne = nint(sum(abs(disp))) + ceiling(sum(abs(this%refl_slope))) + ceiling(maxval(abs(this%refl_height)))
        ! The following one seems to be sufficient
        sumdisp = disp*(-sin(rake)*sin(dip))
        ne1 = max(nint(sum(sumdisp, mask=sumdisp > 0)), -nint(sum(sumdisp, mask=sumdisp < 0)))
        ne1 = max(ne1, ceiling(maxval(abs(this%refl_slope))))
        ne1 = ne1 + maxval(this%refl_height*(1.0 + this%secondary_refl_height_ratio))
        n1 = n1 + 2*ne1 + 2
        if (mod(n1, 2) == 1 .and. this%wave == 'delta') then
            n1 = n1 + 1
        end if

        ! Compute the top and bottom fault dip angles
        ! Note that to ensure the faults have proper dip angle range within the final cropped image,
        ! here I add some extra degrees to the begin and end dip angles
        dips = zeros(n1, nf)
        delta_dip = random(nf, range=this%delta_dip, seed=nint(this%seed*4.5))
        do i = 1, nf
            if (dip(i) <= const_pi_half) then
                dips(:, i) = linspace(dip(i) + delta_dip(i)*const_deg2rad*ne1*1.0/this%n1, &
                    dip(i) - delta_dip(i)*const_deg2rad*(1.0 + ne1*1.0/n1), n1)
                where (dips(:, i) >= const_pi_half)
                    dips(:, i) = const_pi_half*0.99
                end where
                where (dips(:, i) < 0)
                    dips(:, i) = 0
                end where
            else
                dips(:, i) = linspace(dip(i) - delta_dip(i)*const_deg2rad*ne1*1.0/this%n1, &
                    dip(i) + delta_dip(i)*const_deg2rad*(1.0 + ne1*1.0/n1), n1)
                where (dips(:, i) <= const_pi_half)
                    dips(:, i) = 1.01*const_pi_half
                end where
                where (dips(:, i) > const_pi)
                    dips(:, i) = const_pi
                end where
            end if
        end do

        ! Compute the number of padding grid points along x2
        sumdisp = disp*(cos(rake)*sin(strike) - sin(rake)*cos(dips(n1, :))*cos(strike))
        ne2 = max(nint(sum(sumdisp, mask=sumdisp > 0)), -nint(sum(sumdisp, mask=sumdisp < 0)))
        n2 = n2 + 2*ne2 + 2

        sumdisp = disp*(cos(rake)*cos(strike) + sin(rake)*cos(dips(n1, :))*sin(strike))
        ne3 = max(nint(sum(sumdisp, mask=sumdisp > 0)), -nint(sum(sumdisp, mask=sumdisp < 0)))
        n3 = n3 + 2*ne3 + 2

        ! Reflector's shape
        select case (this%refl_shape)

            case ('random')
                r = random(n2, n3, dist='normal', seed=this%seed*5)
                r = gauss_filt(r, [this%refl_smooth, this%refl_smooth])
                r = rescale(r, this%refl_height)

                rt = random(n2, n3, dist='normal', seed=this%seed*5 - 1)
                rt = gauss_filt(rt, [this%refl_smooth_top, this%refl_smooth_top])
                rt = rescale(rt, this%refl_height_top)

            case ('gaussian')

                if (sum(abs(this%refl_mu2)) == 0) then
                    gmu2 = [1.0, this%n2 - 1.0]
                    gsigma2 = [0.05, 0.15]*n2
                else
                    gmu2 = this%refl_mu2
                    gsigma2 = this%refl_sigma2
                end if

                if (sum(abs(this%refl_mu3)) == 0) then
                    gmu3 = [1.0, this%n3 - 1.0]
                    gsigma3 = [0.05, 0.15]*n3
                else
                    gmu3 = this%refl_mu3
                    gsigma3 = this%refl_sigma3
                end if

                mu2 = random(this%ng, range=gmu2, seed=this%seed*5)
                sigma2 = random(this%ng, range=gsigma2, seed=this%seed*6)
                mu3 = random(this%ng, range=gmu3, seed=this%seed*7)
                sigma3 = random(this%ng, range=gsigma3, seed=this%seed*8)
                height = random(this%ng, range=this%refl_height, seed=this%seed*9)

                r = zeros(n2, n3)
                do i = 1, this%ng
                    r = r + rescale(gaussian([mu2(i) + ne2, mu3(i) + ne3], &
                        [sigma2(i), sigma3(i)], &
                        linspace(0.0, n2 - 1.0, n2), linspace(0.0, n3 - 1.0, n3)), [0.0, height(i)])
                end do
                r = rescale(r, this%refl_height)
                rt = rescale(r, this%refl_height_top)

            case ('perlin')
                pn%n1 = n2
                pn%n2 = n3
                pn%seed = this%seed*5
                r = gauss_filt(pn%generate(), [this%refl_smooth, this%refl_smooth])
                r = rescale(r, this%refl_height)

                pn%n1 = n2
                pn%n2 = n3
                pn%seed = this%seed*5 - 1
                rt = gauss_filt(pn%generate(), [this%refl_smooth, this%refl_smooth])
                rt = rescale(rt, this%refl_height)

            case ('custom')
                call assert(size(this%refl, 1) == this%n2 - 2*ne2 .and. &
                    size(this%refl, 2) == this%n3 - 2*ne3, &
                    '<generate_3d_geological_model> Error: size(refl) must = (n2, n3)')
                r = this%refl
                rt = rescale(r, this%refl_height_top)

        end select

        ! Add secondary height fluctuation to the reflectors
        if (this%secondary_refl_height_ratio > 0) then

            sr = random(n2, n3, dist='normal', seed=this%seed*5 - 1)
            sr = gauss_filt(sr, [this%refl_smooth, this%refl_smooth]/5.0)
            sr = rescale(sr, [minval(r), maxval(r)]*this%secondary_refl_height_ratio)
            r = r + sr

            sr = random(n2, n3, dist='normal', seed=this%seed*5 - 2)
            sr = gauss_filt(sr, [this%refl_smooth_top, this%refl_smooth_top]/5.0)
            sr = rescale(sr, [minval(rt), maxval(rt)]*this%secondary_refl_height_ratio)
            rt = rt + sr

        end if

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
        r = r - maxval(r)
        rt = rt - minval(rt) + n1

        nl = nint(this%nl + this%nl*2.0*ne1/(n1 - 2*ne1))

        vv = random(nl - 1, seed=this%seed*6)*this%delta_v
        vv = linspace(this%vmax*(1.0 + ne1*1.0/this%n1), this%vmin*(1.0 - ne1*1.0/this%n1), nl) + vv

        lz = zeros(nl, n2, n3)
        plw = random(nl, range=[1 - this%lwv, 1 + this%lwv], seed=this%seed*7)

        !$omp parallel do private(j, k)
        do k = 1, n3
            do j = 1, n2
                lz(:, j, k) = linspace(r(j, k), rt(j, k), nl)
            end do
        end do
        !$omp end parallel do
        lz = deriv(lz, dim=1)

        !$omp parallel do private(j, k)
        do k = 1, n3
            do j = 1, n2
                lz(:, j, k) = r(j, k) + rescale(integ(lz(:, j, k)*plw), [0.0, rt(j, k) - r(j, k)])
            end do
        end do
        !$omp end parallel do

        ! Stack layers to create a layer-constant velocity model
        w = zeros(n1, n2, n3)
        if (this%yn_rgt) then
            t = zeros(n1, n2, n3)
        end if
        !$omp parallel do private(i, j, k, l)
        do k = 1, n3
            do j = 1, n2
                do i = 1, n1

                    loop_layer: do l = 1, nl - 1
                        if (n1 - i >= lz(l, j, k) .and. n1 - i < lz(l + 1, j, k)) then
                            w(i, j, k) = vv(l)
                            exit loop_layer
                        end if
                    end do loop_layer

                end do

                ! The RGT is linearly interpolated based on the location of reflectors
                ! As such, it may be slightly different from what obtained with rgm2
                if (this%yn_rgt) then
                    t(:, j, k) = ginterp(n1 - 1.0 - lz(:, j, k), linspace(1.0, 0.0, nl), regspace(0.0, 1.0, n1 - 1.0))
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

        if (this%yn_rgt) then
            tt = zeros(n1, n2, n3)
        end if
        if (this%yn_facies) then
            mm = zeros(n1, n2, n3)
        end if

        if (nf >= 1) then

            if (this%yn_regular_fault) then

                rc = random(nf - 1, range=[0.75, 1.25]*(n2 - 2*ne2)/(nf - 1.0), seed=this%seed*12)
                rc = rc*(n2 - 2*ne2)/sum(rc)
                rc = sort(rc)
                f2 = zeros(nf)
                f2(1) = ne2
                f2(2:) = ne2 + integ(rc)

                rc = random(nf - 1, range=[0.75, 1.25]*(n3 - 2*ne3)/(nf - 1.0), seed=this%seed*13)
                rc = rc*(n3 - 2*ne3)/sum(rc)
                rc = sort(rc)
                f3 = zeros(nf)
                f3(1) = ne3
                f3(2:) = ne3 + integ(rc)

                ! Rotate the points of fault centers
                do i = 1, nf
                    pt = rotate_point([f2(i), f3(i)], real(const_pi_half) - strike(i), [0.5*(n2 - 1.0), 0.5*(n3 - 1.0)])
                    f2(i) = pt(1)
                    f3(i) = pt(2)
                end do

            else

                f2 = random(nf, range=[ne2 + 0.1*n2, n2 - ne2 - 0.1*n2], seed=this%seed*16, spacing=0.75*(n2 - 2*ne2 - 0.2*n2)/nf)
                f3 = random(nf, range=[ne3 + 0.1*n3, n3 - ne3 - 0.1*n3], seed=this%seed*17, spacing=0.75*(n3 - 2*ne3 - 0.2*n3)/nf)

            end if

            do fi = 1, nf

                ww = w
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

                                dist = abs(y0 - (a*x0 + b))
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

                                    dist_prev = abs(y0 - (a*x0 + b_prev))
                                    dist = abs(y0 - (a*x0 + b))
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

                                dist = abs(x0 - (a*y0 + b))
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

                                    dist_prev = abs(x0 - (a*y0 + b_prev))
                                    dist = abs(x0 - (a*y0 + b))
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

                ! shift block
                !$omp parallel do private(i, j, k, newi, newj, newk) collapse(3)
                do k = 1, n3
                    do j = 1, n2
                        do i = 1, n1

                            newi = nint(i + (-sin(rake(fi))*sin(dip(fi)))*disp(fi))
                            newj = nint(j + (cos(rake(fi))*sin(strike(fi)) - sin(rake(fi))*cos(dip(fi))*cos(strike(fi)))*disp(fi))
                            newk = nint(k + (cos(rake(fi))*cos(strike(fi)) + sin(rake(fi))*cos(dip(fi))*sin(strike(fi)))*disp(fi))

                            if (newi >= 1 .and. newi <= n1 .and. newj >= 1 .and. newj <= n2 .and. newk >= 1 .and. newk <= n3) then
                                if (fblock(newi, newj, newk)) then
                                    w(newi, newj, newk) = ww(i, j, k)
                                    if (this%yn_rgt) then
                                        t(newi, newj, newk) = tt(i, j, k)
                                    end if
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

        ! Add salt
        if (this%yn_salt) then

            tp = mean(this%salt_max_radius)

            select case (this%refl_shape)
                case ('random')
                    gmax = random(this%nsalt, range=[ne2 + tp*n2/2.0, n2 - ne2 - tp*n2/2.0], seed=this%seed*5)
                    hmax = random(this%nsalt, range=[ne3 + tp*n3/2.0, n3 - ne3 - tp*n3/2.0], seed=this%seed*6)
                case ('gaussian')
                    if (this%nsalt > this%ng) then
                        gmax = [mu2, random(this%nsalt - this%ng, range=[ne2 + tp*n2/2.0, n2 - ne2 - tp*n2/2.0], seed=this%seed*5)]
                        hmax = [mu3, random(this%nsalt - this%ng, range=[ne3 + tp*n3/2.0, n3 - ne3 - tp*n3/2.0], seed=this%seed*6)]
                    else
                        gmax = mu2
                        hmax = mu3
                    end if
            end select

            nex = 11
            this%salt = zeros(n1 + nex, n2, n3)
            rds = random(this%nsalt, range=this%salt_max_radius, seed=this%seed*15 - 1)

            topz = random(n2, n3, seed=this%seed*15 - 2, dist='normal')
            topz = gauss_filt(topz, [this%salt_top_smooth, this%salt_top_smooth])
            topz = rescale(topz, [0.0, this%salt_top_amp])

            do isalt = 1, this%nsalt

                x1 = random(this%nstem + 1, seed=this%seed*15*isalt)
                x2 = random(this%nstem + 1, seed=this%seed*16*isalt)
                x3 = random(this%nstem + 1, seed=this%seed*17*isalt)
                x4 = random(this%nstem + 1, seed=this%seed*18*isalt)

                nd = nint((1.0 - rand(range=this%salt_top_max, seed=this%seed*14*isalt - 1))*this%n1 + ne1)

                x1 = ginterp([0.0, rand(range=[0.05*nd, 0.15*nd], seed=this%seed*19*isalt - 1), &
                    sort(random(this%nstem - 2, range=[0.2*nd, nd - 0.1*nd], seed=this%seed*19*isalt)), nd - 1.0], &
                    x1, linspace(0.0, nd + 10.0, nd + 11), method='pchip')
                x2 = ginterp([0.0, rand(range=[0.05*nd, 0.15*nd], seed=this%seed*20*isalt - 1), &
                    sort(random(this%nstem - 2, range=[0.2*nd, nd - 0.1*nd], seed=this%seed*20*isalt)), nd - 1.0], &
                    x2, linspace(0.0, nd + 10.0, nd + 11), method='pchip')
                x3 = ginterp([0.0, rand(range=[0.05*nd, 0.15*nd], seed=this%seed*21*isalt - 1), &
                    sort(random(this%nstem - 2, range=[0.2*nd, nd - 0.1*nd], seed=this%seed*21*isalt)), nd - 1.0], &
                    x3, linspace(0.0, nd + 10.0, nd + 11), method='pchip')
                x4 = ginterp([0.0, rand(range=[0.05*nd, 0.15*nd], seed=this%seed*22*isalt - 1), &
                    sort(random(this%nstem - 2, range=[0.2*nd, nd - 0.1*nd], seed=this%seed*22*isalt)), nd - 1.0], &
                    x4, linspace(0.0, nd + 10.0, nd + 11), method='pchip')

                x1 = rescale(x1, range=[0.2, 1.0]*rds(isalt)*max(this%n2, this%n3))
                x2 = rescale(x2, range=[0.2, 1.0]*rds(isalt)*max(this%n2, this%n3))
                x3 = rescale(x3, range=[0.2, 1.0]*rds(isalt)*max(this%n2, this%n3))
                x4 = rescale(x4, range=[0.2, 1.0]*rds(isalt)*max(this%n2, this%n3))
                x1 = pad(x1, [n1 - nd, 0], method=['edge', 'edge'])
                x2 = pad(x2, [n1 - nd, 0], method=['edge', 'edge'])
                x3 = pad(x3, [n1 - nd, 0], method=['edge', 'edge'])
                x4 = pad(x4, [n1 - nd, 0], method=['edge', 'edge'])

                slice = zeros(360, n1 + nex)
                !$omp parallel do private(i, pt)
                do i = 1, n1 + nex
                    pt = maxval(abs([x1(i), x2(i), x3(i), x4(i)]))*0.2
                    slice(:, i) = interp([x1(i), x2(i), x3(i), x4(i), x1(i)], 5, 1.0/4, 0.0, &
                        360, 1.0/359.0, 0.0, method='pchip')
                    slice(:, i) = slice(:, i) + rescale(gauss_filt(random(360, seed=this%seed*22*isalt + 1), 6.0), range=[-pt, pt])
                end do
                !$omp end parallel do

                !$omp parallel do private(i, j, k, dist, ag)
                do k = 1, n3
                    do j = 1, n2
                        dist = sqrt((j - gmax(isalt))**2 + (k - hmax(isalt))**2 + 0.0)
                        ag = clip(nint(atan2(k - hmax(isalt) + 0.0, j - gmax(isalt) + float_tiny)*const_rad2deg) + 181, 1, 360)
                        do i = n1 - nd - ceiling(maxval(topz)), n1 + nex
                            if (dist <= slice(ag, i) .and. i >= n1 - nd - topz(j, k)) then
                                this%salt(i, j, k) = 1.0
                                if (i >= 1 .and. i <= n1) then
                                    w(i, j, k) = this%salt_vel
                                    f(i, j, k) = 0.0
                                    fdip(i, j, k) = 0.0
                                    fstrike(i, j, k) = 0.0
                                    frake(i, j, k) = 0.0
                                    fdisp(i, j, k) = 0.0
                                    if (this%yn_rgt) then
                                        t(i, j, k) = 1.0
                                    end if
                                end if
                            end if
                        end do
                    end do
                end do
                !$omp end parallel do
            end do

            this%salt = this%salt(1:n1, 1:n2, 1:n3)

        end if

        ! Final processing

        ! After fault block shifting, the model is w, while RGT is t
        if (this%yn_facies) then
            m = w
        end if

        ww = w
        !$omp parallel do private(i, j, k)
        do k = 1, n3
            do j = 1, n2
                do i = 1, n1 - 1
                    w(i, j, k) = (ww(i + 1, j, k) - ww(i, j, k))/(ww(i + 1, j, k) + ww(i, j, k))
                end do
            end do
        end do
        !$omp end parallel do
        where (f /= 0)
            w = w*0.5
        end where

        this%image = w(ne1 + 1:ne1 + this%n1, ne2 + 1:ne2 + this%n2, ne3 + 1:ne3 + this%n3)
        n1 = this%n1
        n2 = this%n2
        n3 = this%n3

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

            wavelet = zeros(n1)
            f0 = this%f0
            !$omp parallel do private(i, wt)
            do i = 1, n1
                wt = (i - 1.0 - (n1 - 1.0)/2.0)*dt
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
                    wavelet = fourier_filt(wavelet, dt, this%wave_filt_freqs, this%wave_filt_amps)
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

        ! Output
        this%image = this%image(1:this%n1, 1:this%n2, 1:this%n3)

        if (this%yn_fault) then
            this%fault = f(ne1 + 1:ne1 + this%n1, ne2 + 1:ne2 + this%n2, ne3 + 1:ne3 + this%n3)
            this%fault_dip = fdip(ne1 + 1:ne1 + this%n1, ne2 + 1:ne2 + this%n2, ne3 + 1:ne3 + this%n3)*const_rad2deg
            this%fault_strike = fstrike(ne1 + 1:ne1 + this%n1, ne2 + 1:ne2 + this%n2, ne3 + 1:ne3 + this%n3)*const_rad2deg
            this%fault_rake = frake(ne1 + 1:ne1 + this%n1, ne2 + 1:ne2 + this%n2, ne3 + 1:ne3 + this%n3)*const_rad2deg
            this%fault_disp = fdisp(ne1 + 1:ne1 + this%n1, ne2 + 1:ne2 + this%n2, ne3 + 1:ne3 + this%n3)
        end if

        if (this%yn_salt) then
            this%salt = this%salt(ne1 + 1:ne1 + this%n1, ne2 + 1:ne2 + this%n2, ne3 + 1:ne3 + this%n3)
        end if

        if (this%yn_rgt) then
            this%rgt = t(ne1 + 1:ne1 + this%n1, ne2 + 1:ne2 + this%n2, ne3 + 1:ne3 + this%n3)
            this%rgt = rescale(this%rgt, [0.0, 1.0])
        end if

        if (this%yn_facies) then
            this%facies = m(ne1 + 1:ne1 + this%n1, ne2 + 1:ne2 + this%n2, ne3 + 1:ne3 + this%n3)
            if (this%yn_salt) then
                this%facies = rescale(this%facies, [this%vmin, max(this%vmax, this%salt_vel)])
            else
                this%facies = rescale(this%facies, [this%vmin, this%vmax])
            end if
        end if

        if (this%unconf > 0 .and. this%unconf_nl == 0) then
            this%image = 0
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
                this%facies = 0
            end if
            if (this%yn_salt) then
                this%salt = 0
            end if
        end if

    end subroutine generate_3d_geological_model

    subroutine generate_3d_unconformal_geological_model(this)

        type(rgm3_curved), intent(inout) :: this

        type(rgm3_curved), allocatable, dimension(:) :: g
        integer :: iconf, i, j, k
        type(meta_array2_real), allocatable, dimension(:) :: uff
        real, allocatable, dimension(:) :: ufz, aufz
        real, allocatable, dimension(:, :, :) :: rgt_above, rgt_below
        real :: tmin, tmax
        integer, allocatable, dimension(:, :) :: st

        allocate (g(1:this%unconf + 1))

        this%image = zeros(this%n1, this%n2, this%n3)

        if (this%yn_fault) then
            this%fault = zeros(this%n1, this%n2, this%n3)
            this%fault_dip = zeros(this%n1, this%n2, this%n3)
            this%fault_strike = zeros(this%n1, this%n2, this%n3)
            this%fault_rake = zeros(this%n1, this%n2, this%n3)
            this%fault_disp = zeros(this%n1, this%n2, this%n3)
        end if
        if (this%yn_rgt) then
            this%rgt = zeros(this%n1, this%n2, this%n3)
        end if
        if (this%yn_facies) then
            this%facies = zeros(this%n1, this%n2, this%n3)
        end if
        if (this%yn_salt) then
            this%salt = zeros(this%n1, this%n2, this%n3)
        end if

        do i = 1, this%unconf + 1

            if (this%wave == 'delta' .and. allocated(this%wave_filt_freqs)) then
                g(i)%wave_filt_freqs = zeros(size(this%wave_filt_freqs))
                g(i)%wave_filt_amps = zeros(size(this%wave_filt_freqs))
            end if

            g(i) = this

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
                g(i)%refl_height = g(i)%refl_height/2.0
                g(i)%refl_slope = -g(i)%refl_slope/4.0
                g(i)%refl_height_top = g(i)%refl_height_top/2.0
                g(i)%refl_slope_top = -g(i)%refl_slope_top/4.0
            end if

            if (i < this%unconf + 1 .and. this%unconf_nl == 0) then
                g(i)%unconf_nl = 0
            else
                g(i)%unconf_nl = g(i)%nl
            end if

            call generate_3d_geological_model(g(i))

            g(i)%image = g(i)%image/norm2(g(i)%image)

        end do

        if (g(this%unconf + 1)%yn_salt) then
            st = maxloc(g(this%unconf + 1)%salt, 1)
            aufz = clip(this%unconf_z, minval(this%unconf_z), &
                minval(st, mask=(st /= 1))*1.0/this%n1 - maxval(this%unconf_amp))
        else
            aufz = this%unconf_z
        end if

        allocate (uff(1:this%unconf + 1))
        ufz = random(this%unconf, range=aufz, seed=this%seed*31)
        ufz = sort(ufz, order=1)
        do i = 1, this%unconf
            uff(i)%array = random(this%n2, this%n3, seed=g(i)%seed*41*i)
            uff(i)%array = gauss_filt(uff(i)%array, [this%n2*0.2, this%n3*0.2])
            uff(i)%array = rescale(uff(i)%array, &
                [0.0, rand(range=this%unconf_amp, seed=g(i)%seed*51*i)*this%n1]) + ufz(i)*this%n1
        end do

        ! Merge sedimentary units
        this%image = g(this%unconf + 1)%image
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
        end if
        if (this%yn_salt) then
            this%salt = g(this%unconf + 1)%salt
        end if

        do iconf = this%unconf, 1, -1

            rgt_above = 0
            rgt_below = 0

            !$omp parallel do private(i, j, k)
            do k = 1, this%n3
                do j = 1, this%n2

                    ! Image by soft merging
                    this%image(:, j, k) = taper(this%image(:, j, k), protect=[nint(uff(iconf)%array(j, k)) + 5, this%n1], len=[10, 0]) &
                        + taper(g(iconf)%image(:, j, k), protect=[1, nint(uff(iconf)%array(j, k)) - 5], len=[0, 10])

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

                    ! Model perturbation by hard merging
                    if (this%yn_facies) then
                        do i = 1, this%n1
                            if (i < uff(iconf)%array(j, k)) then
                                this%facies(i, j, k) = g(iconf)%facies(i, j, k)
                            end if
                        end do
                    end if

                end do
            end do
            !$omp end parallel do

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

        end do

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

        if (.not. this%custom_psf .and. this%wave /= '') then
            this%psf = g(1)%psf
        end if

    end subroutine generate_3d_unconformal_geological_model

end module geological_model_3d_curved

