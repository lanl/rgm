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


module geological_model_2d_elastic

    use libflit
    use geological_model_utility

    implicit none

    ! 2D random geological model -- elastic
    type rgm2_elastic
        integer :: n1 = 128
        integer :: n2 = 128
        integer :: nf = 4
        integer :: nl = 20
        integer :: seed = -1
        real :: refl_smooth = 20.0
        real :: refl_slope = 0.0
        real :: f0 = 150.0
        real :: fwidth = 2.0
        real, allocatable, dimension(:) :: dip, disp
        real, dimension(1:2) :: noise_smooth = [1.0, 1.0]
        real :: noise_level = 0.25
        real, dimension(1:2) :: image_smooth = [0.0, 0.0]
        character(len=24) :: wave = 'ricker'
        character(len=24) :: refl_shape = 'random'
        integer :: ng = 2
        real, dimension(1:2) :: refl_height = [0.0, 10.0]
        real, dimension(1:2) :: refl_amp = [0.5, 1.0]
        real, dimension(1:2) :: refl_sigma2 = [0.0, 0.0]
        real, dimension(1:2) :: refl_mu2 = [0.0, 0.0]
        real, allocatable, dimension(:) :: refl
        real :: lwv = 0.0
        integer :: unconf = 0
        real, dimension(1:2) :: unconf_z = [0.0, 0.5]
        real, dimension(1:2) :: unconf_amp = [0.05, 0.15]
        integer :: unconf_nl = 99999
        real :: lwmin = 0.5
        logical :: yn_rgt = .false.
        logical :: yn_facies = .false.
        logical :: yn_fault = .true.
        real, allocatable, dimension(:, :) :: image_pp, image_ps, image_sp, image_ss
        real, allocatable, dimension(:, :) :: vp, vs
        real, allocatable, dimension(:, :) :: rgt, facies, fault, fault_dip, fault_disp
        real, dimension(1:2) :: psf_sigma = [5.0, 2.5]
        real, allocatable, dimension(:, :) :: psf
        logical :: custom_psf = .false.
        real :: facies_threshold = 0.0
        character(len=12) :: refl_amp_dist = 'normal'
        character(len=12) :: noise_type = 'normal'
        logical :: yn_conv_noise = .true.
        real :: secondary_refl_amp = 0.0
        logical :: yn_regular_fault = .false.
        real, allocatable, dimension(:) :: wave_filt_freqs, wave_filt_amps

        ! Salt generation is not implemented in the elastic version yet
        logical :: yn_salt = .false.
        integer :: nsalt = 1
        integer :: nstem = 5
        real, dimension(1:2) :: salt_max_radius = [0.05, 0.15]
        real, dimension(1:2) :: salt_top_max = [0.2, 0.4]
        real :: salt_top_smooth = 10.0
        character(len=24) :: salt_type = 'dome'
        real :: vpmin = 2000.0
        real :: vpmax = 4000.0
        real :: vpvsratiomin = const_sqrt2
        real :: vpvsratiomax = const_sqrt3
        real :: perturb_max = 0.2
        real :: salt_vel = 5000.0
        real, allocatable, dimension(:, :) :: salt
        real :: salt_top_amp = 20.0
        real :: salt_noise = 0.0
        real :: salt_amp = 2.0

    contains
        procedure, private :: create_psf => create_psf_2d
        procedure, public :: generate => generate_2d_elastic
    end type rgm2_elastic

    private
    public :: rgm2_elastic

contains

    subroutine generate_2d_elastic(this)

        class(rgm2_elastic), intent(inout) :: this

        if (this%nf == 0) then
            this%yn_fault = .false.
        end if

        if (this%unconf == 0) then
            call generate_2d_geological_model_elastic(this)
        else
            call generate_2d_unconformal_geological_model_elastic(this)
        end if

    end subroutine generate_2d_elastic

    subroutine generate_2d_geological_model_elastic(this)

        type(rgm2_elastic), intent(inout) :: this

        real, allocatable, dimension(:, :) :: w_pp, w_ps, w_sp, w_ss, f, t, mp, ms
        real, allocatable, dimension(:, :) :: ww, ww_pp, ww_ps, ww_sp, ww_ss, ff, cf, tt, mmp, mms
        integer :: nf, nl, n1, n2, i, j, fi, ne1, ne2, depth, newi, newj
        real, allocatable, dimension(:) :: disp, dip, refl, reflf, phi, f1, f2, r, sr
        real, allocatable, dimension(:) :: refl_pp, refl_ps, refl_sp, refl_ss, rr_pp, rr_ps, rr_sp, rr_ss, pp
        real, allocatable, dimension(:) :: qqp, qqs, imp, impp, imps
        integer, allocatable, dimension(:) :: ri
        real :: z0, x0, z1, x1, z2, x2, dt, fwidth
        real, allocatable, dimension(:) :: mu, sigma, height, layer_thick
        real :: lw
        integer :: nlt
        real, dimension(1:2) :: mu2, sigma2
        real, allocatable, dimension(:) :: sumdisp, rc1, rc2, rc
        real, allocatable, dimension(:) :: rt
        real, allocatable, dimension(:) :: wp, uv, ratio, lr
        integer :: nuv, nsf

        fwidth = this%fwidth
        dt = 0.001

        n1 = this%n1
        n2 = this%n2
        nf = this%nf

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

                !                ! Dip angles
                !                dip = random(nf, range=this%dip, seed=this%seed)*const_deg2rad
                !
                !                ! Fault displacements
                !                disp = random(nf, range=this%disp, seed=this%seed*2)
                !                disp(1:nf:2) = -disp(1:nf:2)

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

        ! The following seems to be sufficient
        sumdisp = disp*sin(dip)
        ne1 = max(nint(sum(sumdisp, mask=sumdisp > 0)), -nint(sum(sumdisp, mask=sumdisp < 0)))
        ne1 = max(ne1, ceiling(abs(this%refl_slope)))
        ne1 = ne1 + maxval(this%refl_height*(1.0 + this%secondary_refl_amp))
        n1 = n1 + 2*ne1 + 2
        if (mod(n1, 2) == 1 .and. this%wave == 'delta') then
            n1 = n1 + 1
        end if

        sumdisp = disp*cos(dip)
        ne2 = max(nint(sum(sumdisp, mask=sumdisp > 0)), -nint(sum(sumdisp, mask=sumdisp < 0)))
        n2 = n2 + 2*ne2 + 2

        ! Reflector's shape
        select case (this%refl_shape)

            case ('random')
                r = random(n2, dist='normal', seed=this%seed*3)
                r = gauss_filt(r, this%refl_smooth)
                r = rescale(r, this%refl_height)

            case ('gaussian')

                if (sum(abs(this%refl_mu2)) == 0) then
                    mu2 = [1.0, this%n2 - 1.0]
                    sigma2 = [0.05, 0.15]*n2
                else
                    mu2 = this%refl_mu2
                    sigma2 = this%refl_sigma2
                end if

                mu = random(this%ng, range=mu2, seed=this%seed*3)
                sigma = random(this%ng, range=sigma2, seed=this%seed*4)
                height = random(this%ng, range=this%refl_height, seed=this%seed*5)

                r = zeros(n2)
                do i = 1, this%ng
                    r = r + rescale(gaussian(mu(i) + ne2, sigma(i), linspace(0.0, n2 - 1.0, n2)), [0.0, height(i)])
                end do
                if (this%lwv >= 0) then
                    r = -r
                end if

            case ('custom')
                call assert(size(this%refl) == this%n2 - 2*ne2, &
                    '<generate_2d_geological_model> Error: size(refl) must = n2')
                r = this%refl
                r = -r

        end select

        if (this%secondary_refl_amp > 0) then
            sr = random(n2, dist='normal', seed=this%seed*3 - 1)
            sr = gauss_filt(sr, this%refl_smooth/5.0)
            sr = rescale(sr, [minval(r), maxval(r)]*this%secondary_refl_amp)
            r = r + sr
        end if

        !$omp parallel do private(j)
        do j = 1, n2
            r(j) = r(j) + (j - 1.0)*this%refl_slope/this%n2
        end do
        !$omp end parallel do
        r = r - mean(r)

        ! Reflectivity
        nl = nint(this%nl + this%nl*2.0*ne1/(n1 - 2*ne1))

        rc1 = random(nint(nl/2.0), seed=this%seed*6, range=this%refl_amp, dist=this%refl_amp_dist)
        if (this%refl_amp_dist == 'normal' .or. this%refl_amp_dist == 'gaussian') then
            rc1 = rescale(rc1, this%refl_amp)
        end if
        rc2 = random(nl - nint(nl/2.0), seed=this%seed*7, range=this%refl_amp, dist=this%refl_amp_dist)
        if (this%refl_amp_dist == 'normal' .or. this%refl_amp_dist == 'gaussian') then
            rc2 = rescale(rc2, this%refl_amp)
        end if
        rc2 = -rc2

        lw = (n1 - 1.0)/nl
        ri = nint(linspace(1.0, n1 - lw*2.0/3.0, nl) + random(nl, range=[0.0, lw*2.0/3.0], seed=this%seed*8))
        ri = random_permute(ri, seed=this%seed*9)

        rc = zeros(nl)
        rc = [rc1, rc2]
        rc = rc(random_order(nl, seed=this%seed*10))

        refl = zeros(n1)
        j = 1
        do i = 1, n1
            if (any(i == ri)) then
                refl(i) = rc(j)
                j = j + 1
            end if
        end do

        ! Model and reflectivity
        reflf = refl
        where (abs(reflf) < this%facies_threshold*maxval(abs(reflf)))
            reflf = 0
        end where
        imp = integ(reflf)
        imp = imp - mean(imp)

        uv = unique(imp)
        nuv = size(uv)
        wp = linspace(0.0, 1.0, nuv)
        ! Vp/Vs ratio
        ratio = random(nuv, range=[this%vpvsratiomin, this%vpvsratiomax])

        reflf = 0
        lr = zeros_like(reflf)
        do i = 1, nuv
            where (imp == uv(i))
                reflf = wp(i)
                lr = ratio(i)
            end where
        end do
        reflf = rescale(reflf, [0.0, 1.0])
        imp = imp/maxval(imp)
        imp = reflf + this%perturb_max*imp

        ! Create Vp, Vs models using imp
        impp = rescale(imp, [this%vpmin, this%vpmax])
        imps = impp/lr

        ! Compute PP, PS, SP, SS reflectivities
        refl_pp = zeros(n1)
        refl_ps = zeros(n1)
        refl_sp = zeros(n1)
        refl_ss = zeros(n1)
        rt = zeros(4)
        !$omp parallel do private(i, j, rt)
        do i = 1, n1 - 1
            rt = 0
            do j = 0, 5
                rt = rt + elastic_reflection_coefs(real(sin((j*3.0)*const_deg2rad)/impp(i)), &
                    impp(i), imps(i), 1.0, impp(i + 1), imps(i + 1), 1.0)

            end do
            rt = rt/6.0

            refl_pp(i) = rt(1)
            refl_ps(i) = rt(2)
            refl_sp(i) = rt(3)
            refl_ss(i) = rt(4)
        end do
        !$omp end parallel do

        if (this%yn_rgt) then
            phi = linspace(0.0, 1.0, n1)
        end if

        ! Create non-fault image and RGT
        w_pp = zeros(n1, n2)
        w_ps = zeros(n1, n2)
        w_sp = zeros(n1, n2)
        w_ss = zeros(n1, n2)

        if (this%yn_rgt) then
            t = zeros(n1, n2)
        end if

        mp = zeros(n1, n2)
        ms = zeros(n1, n2)

        ! Spatially variant layer thickness
        layer_thick = rescale(r, [1.0, 1.0 + abs(this%lwv)])

        nlt = nint(maxval(layer_thick*n1))
        rr_pp = zeros(nlt)
        rr_ps = zeros(nlt)
        rr_sp = zeros(nlt)
        rr_ss = zeros(nlt)
        pp = zeros(nlt)
        qqp = zeros(nlt)
        qqs = zeros(nlt)

        !$omp parallel do private(i, j, nlt, depth, rr_pp, rr_ps, rr_sp, rr_ss, pp, qqp, qqs)
        do j = 1, n2

            depth = nint(r(j))
            nlt = nint(n1*layer_thick(j))

            ! Reflectivity
            rr_pp(1:nlt) = interp(refl_pp, n1, 1.0, r(j), nlt, 1.0/layer_thick(j), depth*1.0, 'linear')
            rr_ps(1:nlt) = interp(refl_ps, n1, 1.0, r(j), nlt, 1.0/layer_thick(j), depth*1.0, 'linear')
            rr_sp(1:nlt) = interp(refl_sp, n1, 1.0, r(j), nlt, 1.0/layer_thick(j), depth*1.0, 'linear')
            rr_ss(1:nlt) = interp(refl_ss, n1, 1.0, r(j), nlt, 1.0/layer_thick(j), depth*1.0, 'linear')
            ! RGT
            if (this%yn_rgt) then
                pp(1:nlt) = interp(phi, n1, 1.0, r(j), nlt, 1.0/layer_thick(j), depth*1.0, 'linear')
            end if
            ! Model perturbation
            if (this%yn_facies) then
                qqp(1:nlt) = interp(impp, n1, 1.0, r(j), nlt, 1.0/layer_thick(j), depth*1.0, 'nearest')
                qqs(1:nlt) = interp(imps, n1, 1.0, r(j), nlt, 1.0/layer_thick(j), depth*1.0, 'nearest')
            end if

            do i = 1, n1
                if (depth + i >= 1 .and. depth + i <= n1) then
                    ! Reflectivity
                    w_pp(depth + i, j) = rr_pp(i)
                    w_ps(depth + i, j) = rr_ps(i)
                    w_sp(depth + i, j) = rr_sp(i)
                    w_ss(depth + i, j) = rr_ss(i)
                    ! RGT
                    if (this%yn_rgt) then
                        t(depth + i, j) = pp(i)
                    end if
                    ! Model perturbation
                    if (this%yn_facies) then
                        mp(depth + i, j) = qqp(i)
                        ms(depth + i, j) = qqs(i)
                    end if
                end if
            end do

        end do
        !$omp end parallel do

        ! Add faults
        ww_pp = zeros(n1, n2)
        ww_ps = zeros(n1, n2)
        ww_sp = zeros(n1, n2)
        ww_ss = zeros(n1, n2)
        f = zeros(n1, n2)
        ff = zeros(n1, n2)
        cf = zeros(n1, n2)

        if (this%yn_rgt) then
            tt = zeros(n1, n2)
        end if

        if (this%yn_facies) then
            mmp = zeros(n1, n2)
            mms = zeros(n1, n2)
        end if

        if (nf >= 1) then
            if (this%yn_regular_fault) then

                f1 = 0.5*n1*ones(nf)

                rc = random(nf - 1, range=[0.75, 1.25]*(n2 - 2*ne2)/(nf - 1.0), seed=this%seed*12)
                rc = rc*(n2 - 2*ne2)/sum(rc)
                rc = sort(rc)
                f2 = zeros(nf)
                f2(1) = ne2
                f2(2:) = ne2 + integ(rc)

            else

                f1 = random(nf, range=[ne1 + 5.0, n1 - ne1 - 5.0], seed=this%seed*11)
                f2 = random(nf, range=[ne2 + 0.05*this%n2, n2 - ne2 - 0.05*this%n2], seed=this%seed*12)

            end if

            do fi = 1, nf

                x1 = f2(fi)
                z1 = f1(fi)
                x2 = x1 - fwidth*sin(dip(fi))
                z2 = z1 - fwidth*cos(dip(fi))

                ww_pp = w_pp
                ww_ps = w_ps
                ww_sp = w_sp
                ww_ss = w_ss
                ff = f
                cf = 0

                if (this%yn_rgt) then
                    tt = t
                end if

                if (this%yn_facies) then
                    mmp = mp
                    mms = ms
                end if

                ! Make randomly oriented faults with random displacements
                !$omp parallel do private(i, j, x0, z0) collapse(2)
                do j = 1, n2
                    do i = 1, n1

                        ! Coordinates of a grid point
                        x0 = j - 1.0
                        z0 = i - 1.0
                        if ((z0 - z1)*(z2 - z1) + (x0 - x1)*(x2 - x1) >= 0 .and. &
                                (z0 - z1)*(z2 - z1) + (x0 - x1)*(x2 - x1) &
                                <= (z2 - z1)**2 + (x2 - x1)**2) then
                            f(i, j) = fi*1.0
                            cf(i, j) = 1.0
                        end if

                    end do
                end do
                !$omp end parallel do

                ! Shift blocks
                !$omp parallel do private(i, j, x0, z0, newi, newj) collapse(2)
                do j = 1, n2
                    do i = 1, n1

                        ! Shift center to ensure alignment
                        x0 = j - 1.0 + 0.5*sign(1.0, x2 - x1)
                        z0 = i - 1.0 + 0.5*sign(1.0, z2 - z1)

                        if ((z0 - z2)*(z2 - z1) + (x0 - x2)*(x2 - x1) >= 0) then

                            newi = clip(nint(i - disp(fi)*sin(dip(fi))), 1, n1)
                            newj = clip(nint(j + disp(fi)*cos(dip(fi))), 1, n2)

                            if (newi >= 1 .and. newi <= n1 .and. newj >= 1 .and. newj <= n2) then
                                ! Reflectivity
                                w_pp(newi, newj) = ww_pp(i, j)
                                w_ps(newi, newj) = ww_ps(i, j)
                                w_sp(newi, newj) = ww_sp(i, j)
                                w_ss(newi, newj) = ww_ss(i, j)
                                ! RGT
                                if (this%yn_rgt) then
                                    t(newi, newj) = tt(i, j)
                                end if
                                ! Model perturbation
                                if (this%yn_facies) then
                                    mp(newi, newj) = mmp(i, j)
                                    ms(newi, newj) = mms(i, j)
                                end if
                                ! Fault
                                if (cf(newi, newj) == 0) then
                                    f(newi, newj) = ff(i, j)
                                end if
                            end if

                        end if
                    end do
                end do
                !$omp end parallel do

            end do
        end if

        this%image_pp = w_pp(ne1 + 1:ne1 + this%n1 + 2, ne2 + 1:ne2 + this%n2 + 2)
        this%image_ps = w_ps(ne1 + 1:ne1 + this%n1 + 2, ne2 + 1:ne2 + this%n2 + 2)
        this%image_sp = w_sp(ne1 + 1:ne1 + this%n1 + 2, ne2 + 1:ne2 + this%n2 + 2)
        this%image_ss = w_ss(ne1 + 1:ne1 + this%n1 + 2, ne2 + 1:ne2 + this%n2 + 2)
        n1 = size(this%image_pp, 1)
        n2 = size(this%image_pp, 2)

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

            call this%create_psf(n1, n2, dt, this%f0)
            this%image_pp = conv(this%image_pp, this%psf, 'same')
            call this%create_psf(n1, n2, dt, this%f0*mean(0.5*(lr + ones_like(lr))))
            this%image_ps = conv(this%image_ps, this%psf, 'same')
            this%image_sp = conv(this%image_sp, this%psf, 'same')
            call this%create_psf(n1, n2, dt, this%f0*mean(lr))
            this%image_ss = conv(this%image_ss, this%psf, 'same')

        end if

        ! Image smooth
        call assert(all(this%image_smooth >= 0), ' <generate_2d_geological_model> Error: image_smooth must >= 0')
        if (any(this%image_smooth > 0)) then
            this%image_pp = gauss_filt(this%image_pp, this%image_smooth)
            this%image_ps = gauss_filt(this%image_ps, this%image_smooth)
            this%image_sp = gauss_filt(this%image_sp, this%image_smooth)
            this%image_ss = gauss_filt(this%image_ss, this%image_smooth)
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

        ! Output
        this%image_pp = this%image_pp(1:this%n1, 1:this%n2)
        this%image_ps = this%image_ps(1:this%n1, 1:this%n2)
        this%image_sp = this%image_sp(1:this%n1, 1:this%n2)
        this%image_ss = this%image_ss(1:this%n1, 1:this%n2)

        if (this%yn_fault) then
            this%fault = f(ne1 + 1:ne1 + this%n1, ne2 + 1:ne2 + this%n2)
            ! Fault attributes
            this%fault_dip = this%fault
            this%fault_disp = this%fault
            do fi = 1, this%nf
                !$omp parallel do private(i, j) collapse(2)
                do j = 1, this%n2
                    do i = 1, this%n1
                        if (nint(this%fault(i, j)) == fi) then
                            this%fault_dip(i, j) = dip(fi)*const_rad2deg
                            this%fault_disp(i, j) = disp(fi)
                        end if
                    end do
                end do
                !$omp end parallel do
            end do
        end if

        if (this%yn_rgt) then
            this%rgt = t(ne1 + 1:ne1 + this%n1, ne2 + 1:ne2 + this%n2)
            this%rgt = rescale(this%rgt, [0.0, 1.0])
        end if

        if (this%yn_facies) then
            this%vp = mp(ne1 + 1:ne1 + this%n1, ne2 + 1:ne2 + this%n2)
            this%vs = ms(ne1 + 1:ne1 + this%n1, ne2 + 1:ne2 + this%n2)
        end if

        if (this%unconf_nl == 0) then
            this%image_pp = 0
            this%image_ps = 0
            this%image_sp = 0
            this%image_ss = 0
            if (this%yn_rgt) then
                this%rgt = 0
            end if
            if (this%yn_fault) then
                this%fault = 0
                this%fault_dip = 0
                this%fault_disp = 0
            end if
            if (this%yn_facies) then
                this%vp = 0
                this%vs = 0
            end if
            if (this%yn_salt) then
                this%salt = 0
            end if
        end if

    end subroutine generate_2d_geological_model_elastic

    subroutine create_psf_2d(this, n1, n2, dt, freq)

        class(rgm2_elastic), intent(inout) :: this
        integer :: n1, n2
        real :: dt
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
                call assert(size(this%wave_filt_freqs) == size(this%wave_filt_amps))
                wavelet = fourier_filt(wavelet, dt, this%wave_filt_freqs, this%wave_filt_amps)
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

    end subroutine create_psf_2d

    subroutine generate_2d_unconformal_geological_model_elastic(this)

        type(rgm2_elastic), intent(inout) :: this

        type(rgm2_elastic), allocatable, dimension(:) :: g
        integer :: iconf, i, j
        type(meta_array1_real), allocatable, dimension(:) :: uff
        real, allocatable, dimension(:) :: ufz, aufz
        real, allocatable, dimension(:, :) :: rgt_above, rgt_below
        real :: tmin, tmax
        integer, allocatable, dimension(:) :: st

        allocate (g(1:this%unconf + 1))

        this%image_pp = zeros(this%n1, this%n2)
        this%image_ps = zeros(this%n1, this%n2)
        this%image_sp = zeros(this%n1, this%n2)
        this%image_ss = zeros(this%n1, this%n2)

        if (this%yn_fault) then
            this%fault = zeros(this%n1, this%n2)
            this%fault_dip = zeros(this%n1, this%n2)
            this%fault_disp = zeros(this%n1, this%n2)
        end if
        if (this%yn_rgt) then
            this%rgt = zeros(this%n1, this%n2)
        end if
        if (this%yn_facies) then
            this%vp = zeros(this%n1, this%n2)
            this%vs = zeros(this%n1, this%n2)
        end if
        if (this%yn_salt) then
            this%salt = zeros(this%n1, this%n2)
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
            end if

            if (i < this%unconf + 1 .and. this%unconf_nl == 0) then
                g(i)%unconf_nl = 0
            else
                g(i)%unconf_nl = g(i)%nl
            end if

            call generate_2d_geological_model_elastic(g(i))

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

            uff(i)%array = random(this%n2, seed=g(i)%seed*41*i)
            uff(i)%array = gauss_filt(uff(i)%array, this%n2*0.2)
            uff(i)%array = rescale(uff(i)%array, &
                [0.0, rand(range=this%unconf_amp, seed=g(i)%seed*51*i)*this%n1]) + ufz(i)*this%n1

        end do

        ! Merge sedimentary units
        this%image_pp = g(this%unconf + 1)%image_pp
        this%image_ps = g(this%unconf + 1)%image_ps
        this%image_sp = g(this%unconf + 1)%image_sp
        this%image_ss = g(this%unconf + 1)%image_ss

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
            this%vp = g(this%unconf + 1)%vp
            this%vs = g(this%unconf + 1)%vs
        end if
        if (this%yn_salt) then
            this%salt = g(this%unconf + 1)%salt
        end if

        do iconf = this%unconf, 1, -1

            rgt_above = 0
            rgt_below = 0

            !$omp parallel do private(i, j)
            do j = 1, this%n2

                ! Image by soft merging
                this%image_pp(:, j) = &
                    taper(this%image_pp(:, j), protect=[nint(uff(iconf)%array(j)) + 5, this%n1], len=[10, 0]) &
                    + taper(g(iconf)%image_pp(:, j), protect=[1, nint(uff(iconf)%array(j)) - 5], len=[0, 10])
                this%image_ps(:, j) = &
                    taper(this%image_ps(:, j), protect=[nint(uff(iconf)%array(j)) + 5, this%n1], len=[10, 0]) &
                    + taper(g(iconf)%image_ps(:, j), protect=[1, nint(uff(iconf)%array(j)) - 5], len=[0, 10])
                this%image_sp(:, j) = &
                    taper(this%image_sp(:, j), protect=[nint(uff(iconf)%array(j)) + 5, this%n1], len=[10, 0]) &
                    + taper(g(iconf)%image_sp(:, j), protect=[1, nint(uff(iconf)%array(j)) - 5], len=[0, 10])
                this%image_ss(:, j) = &
                    taper(this%image_ss(:, j), protect=[nint(uff(iconf)%array(j)) + 5, this%n1], len=[10, 0]) &
                    + taper(g(iconf)%image_ss(:, j), protect=[1, nint(uff(iconf)%array(j)) - 5], len=[0, 10])

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

                ! Model perturbation by hard merging
                if (this%yn_facies) then
                    do i = 1, this%n1
                        if (i < uff(iconf)%array(j)) then
                            this%vp(i, j) = g(iconf)%vp(i, j)
                            this%vs(i, j) = g(iconf)%vs(i, j)
                        end if
                    end do
                end if

            end do
            !$omp end parallel do

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

    end subroutine generate_2d_unconformal_geological_model_elastic

end module geological_model_2d_elastic
