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
    ! The module is mostly a simplified version of geological_model_2d
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
    ! As such, there are some changes on the parameters of rgm2_curved
    ! compared with rgm2, but the changes are minimized to keep consistency
    !

    use libflit
    use geological_model_utility

    implicit none

    ! 2D random geological model
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
        real :: noise_level = 0.25
        !> Source wavelet for reflector convolution, can be one of
        !> ricker, ricker_deriv, gaussian, gaussian_deriv, sinc, delta
        character(len=24) :: wave = 'ricker'
        !> Shape of reflectors, can be one of
        !> random, gaussian, custom
        !> When = custom, must specifiy the array refl
        character(len=24) :: refl_shape = 'random'
        real, allocatable, dimension(:) :: refl
        !> Number of Gaussians for refl_shape = gaussian
        integer :: ng = 2
        !> Range of reflector's heights
        real, dimension(1:2) :: refl_height = [0.0, 10.0]
        real, dimension(1:2) :: refl_height_top = [0.0, 10.0]
        !> Range of Gaussian standard devision for refl_shape = gaussian
        real, dimension(1:2) :: refl_sigma2 = [0.0, 0.0]
        !> Range of Gaussian mean for refl_shape = gaussian
        real, dimension(1:2) :: refl_mu2 = [0.0, 0.0]
        !> The  thickness of layers varies from [1 - lwv, 1 + lwv] of average layer thickness
        real :: lwv = 0.1

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
        !> Max height of secondary random fluctuations added to reflectors in terms of max height of base reflectors
        real :: secondary_refl_height_ratio = 0.0
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
        !> Range of height of unconformity interfaces in depth, measured in terms of fraction of the entire depth,
        !> which must fall in [0, 1]; smaller values represent flatter unconformity interfaces
        real, dimension(1:2) :: unconf_amp = [0.05, 0.15]
        !> Number of reflectors above the unconformity interfaces; when
        !> = 1 the region above the unconformity represents water
        integer :: unconf_nl = 99999

        !==============================================================================================
        !> Whether or not to insert salt bodies into the model
        logical :: yn_salt = .false.
        !> Number of salt bodies
        integer :: nsalt = 1
        !> Number of random control points for creating the vertical shape of salt bodies
        integer :: nstem = 5
        !> Range of max horizontal radius of salt bodies in terms of fraction of n2
        real, dimension(1:2) :: salt_max_radius = [0.05, 0.15]
        !> Range of max depth of salt body top interfaces in terms of fraction of n1
        real, dimension(1:2) :: salt_top_max = [0.2, 0.4]
        !> sigma of Gaussian filter for smoothing the top of salt body
        real :: salt_top_smooth = 10.0
        !> Reserved parameter not in use at this momoment
        character(len=24) :: salt_type = 'dome'
        !> Salt body velocity
        real :: salt_vel = 5000.0
        !> Array for holding the salt bodies
        real, allocatable, dimension(:, :) :: salt
        !> Max height of the random salt body top
        real :: salt_top_amp = 20.0

    contains
        procedure :: generate => generate_2d
    end type rgm2_curved

    private
    public :: rgm2_curved

contains

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
        real, allocatable, dimension(:) :: disp, dip, f2, r, rt, sr, wavelet, vv
        real :: f0, dt, fwidth, wt
        real, allocatable, dimension(:) :: mu, sigma, height
        real, dimension(1:2) :: mu2, sigma2
        real, allocatable, dimension(:, :) :: psf
        real, allocatable, dimension(:) :: psf1, psf2
        real, allocatable, dimension(:) :: sumdisp
        real, allocatable, dimension(:, :) :: lz
        real, allocatable, dimension(:) :: plw, delta_dip
        real, allocatable, dimension(:, :) :: fdip, fdisp, ffdip, ffdisp
        integer :: l, nsf
        real :: b, b_prev, dxys, dist, dist_prev, theta, x0
        logical, allocatable, dimension(:, :) :: fblock
        real, allocatable, dimension(:, :) :: dips
        real, allocatable, dimension(:) :: x1, x2, rds, topz, rc
        integer :: nd, isalt
        real, allocatable, dimension(:) :: gmax
        real :: tp
        integer :: nex
        type(fractal_noise_1d) :: pn

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

        ! Compute the number of padding grid points along x1
        sumdisp = disp*sin(dip)
        ne1 = max(nint(sum(sumdisp, mask=sumdisp > 0)), -nint(sum(sumdisp, mask=sumdisp < 0)))
        ne1 = max(ne1, ceiling(abs(this%refl_slope)), ceiling(abs(this%refl_slope_top)))
        ne1 = ne1 + max(maxval(this%refl_height*(1.0 + this%secondary_refl_height_ratio)), &
            maxval(this%refl_height_top*(1.0 + this%secondary_refl_height_ratio)))
        n1 = n1 + 2*ne1 + 2
        if (mod(n1, 2) == 1 .and. this%wave == 'delta') then
            n1 = n1 + 1
        end if

        ! Compute the top and bottom fault dip angles
        ! Note that to ensure the faults have proper dip angle range within the final cropped image,
        ! here I add some extra degrees to the begin and end dip angles
        dips = zeros(n1, nf)
        delta_dip = random(nf, range=this%delta_dip, seed=nint(this%seed*2.5))
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
        sumdisp = disp*cos(dips(n1, :))
        ne2 = max(nint(sum(sumdisp, mask=sumdisp > 0)), -nint(sum(sumdisp, mask=sumdisp < 0)))
        n2 = n2 + 2*ne2 + 2

        ! Reflector's shape at the bottom (r) and at the top (rt)
        select case (this%refl_shape)

            case ('random')
                r = random(n2, dist='normal', seed=this%seed*3)
                r = gauss_filt(r, this%refl_smooth)
                r = rescale(r, this%refl_height)

                rt = random(n2, dist='normal', seed=this%seed*3 - 1)
                rt = gauss_filt(rt, this%refl_smooth_top)
                rt = rescale(rt, this%refl_height_top)

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
                r = rescale(r, this%refl_height)
                rt = rescale(r, this%refl_height_top)

            case ('perlin')
                pn%n1 = n2
                pn%seed = this%seed*3
                r = gauss_filt(pn%generate(), this%refl_smooth)
                r = rescale(r, this%refl_height)

                pn%n1 = n2
                pn%seed = this%seed*3 - 1
                rt = gauss_filt(pn%generate(), this%refl_smooth)
                rt = rescale(rt, this%refl_height)

            case ('custom')
                call assert(size(this%refl) == this%n2 - 2*ne2, &
                    '<generate_2d_geological_model> Error: size(refl) must = n2')
                r = this%refl
                rt = rescale(r, this%refl_height_top)

        end select

        ! Add secondary height fluctuation to the reflectors
        if (this%secondary_refl_height_ratio > 0) then

            sr = random(n2, dist='normal', seed=this%seed*3 - 1)
            sr = gauss_filt(sr, this%refl_smooth/5.0)
            sr = rescale(sr, [minval(r), maxval(r)]*this%secondary_refl_height_ratio)
            r = r + sr

            sr = random(n2, dist='normal', seed=this%seed*3 - 2)
            sr = gauss_filt(sr, this%refl_smooth_top/5.0)
            sr = rescale(sr, [minval(r), maxval(r)]*this%secondary_refl_height_ratio)
            rt = rt + sr

        end if

        !$omp parallel do private(j)
        do j = 1, n2
            r(j) = r(j) + (j - 1.0)*this%refl_slope/this%n2
            rt(j) = rt(j) + (j - 1.0)*this%refl_slope_top/this%n2
        end do
        !$omp end parallel do
        r = r - maxval(r)
        rt = rt - minval(rt) + n1

        nl = nint(this%nl + this%nl*2.0*ne1/(n1 - 2*ne1))

        vv = random(nl - 1, seed=this%seed*6)*this%delta_v
        vv = linspace(this%vmax*(1.0 + ne1*1.0/this%n1), this%vmin*(1.0 - ne1*1.0/this%n1), nl) + vv

        lz = zeros(nl, n2)
        plw = random(nl, range=[1 - this%lwv, 1 + this%lwv], seed=this%seed*7)

        !$omp parallel do private(j)
        do j = 1, n2
            lz(:, j) = linspace(r(j), rt(j), nl)
        end do
        !$omp end parallel do
        lz = deriv(lz, dim=1)
        !$omp parallel do private(j)
        do j = 1, n2
            lz(:, j) = r(j) + rescale(cumsum(lz(:, j)*plw), [0.0, rt(j) - r(j)])
        end do
        !$omp end parallel do

        ! Stack layers to create a layer-constant velocity model
        w = zeros(n1, n2)
        if (this%yn_rgt) then
            t = zeros(n1, n2)
        end if

        !$omp parallel do private(i, j, l)
        do j = 1, n2
            do i = 1, n1

                loop_layer: do l = 1, nl - 1
                    if (n1 - i >= lz(l, j) .and. n1 - i < lz(l + 1, j)) then
                        w(i, j) = vv(l)
                        exit loop_layer
                    end if
                end do loop_layer

            end do

            ! The RGT is linearly interpolated based on the location of reflectors
            ! As such, it may be slightly different from what obtained with rgm2
            if (this%yn_rgt) then
                t(:, j) = ginterp(n1 - 1.0 - lz(:, j), linspace(1.0, 0.0, nl), regspace(0.0, 1.0, n1 - 1.0))
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
        if (this%yn_rgt) then
            tt = zeros(n1, n2)
        end if
        if (this%yn_facies) then
            mm = zeros(n1, n2)
        end if

        if (nf >= 1) then

            if (this%yn_regular_fault) then

                rc = random(nf - 1, range=[0.75, 1.25]*(n2 - 2*ne2)/(nf - 1.0), seed=this%seed*12)
                rc = rc*(n2 - 2*ne2)/sum(rc)
                rc = sort(rc)
                f2 = zeros(nf)
                f2(1) = ne2
                f2(2:) = ne2 + cumsum(rc)
                if (.not. this%yn_group_faults) then
                    f2 = random_permute(f2, seed=this%seed*12 + 1)
                end if

            else

                f2 = random(nf, range=[ne2 + 0.1*n2, n2 - ne2 - 0.1*n2], seed=this%seed*12, spacing=0.75*(n2 - 2*ne2 - 0.2*n2)/nf)

            end if

            do fi = 1, nf

                ww = w
                if (this%yn_rgt) then
                    tt = t
                end if
                ff = f
                ffdip = fdip
                ffdisp = fdisp
                cf = 0

                b_prev = 0
                fblock = falses(n1, n2)

                do i = 1, n1

                    theta = dips(i, fi)

                    if (i == 1) then
                        b = f2(fi)
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
                                w(newi, newj) = ww(i, j)
                                if (this%yn_rgt) then
                                    t(newi, newj) = tt(i, j)
                                end if
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

        ! Add salt
        if (this%yn_salt) then

            tp = mean(this%salt_max_radius)

            select case (this%refl_shape)
                case ('random', 'perlin', 'custom')
                    gmax = random(this%nsalt, range=[ne2 + tp*this%n2/2.0, n2 - ne2 - tp*this%n2/2.0], seed=this%seed*5)
                case ('gaussian')
                    if (this%nsalt > this%ng) then
                        gmax = [mu, random(this%nsalt - this%ng, range=[ne2 + tp*n2/2.0, n2 - ne2 - tp*n2/2.0], seed=this%seed*5)]
                    else
                        gmax = mu
                    end if
            end select

            nex = 11
            this%salt = zeros(n1 + nex, n2)
            rds = random(this%nsalt, range=this%salt_max_radius, seed=this%seed*14 - 1)

            ! Build top interface of the salt
            topz = random(n2, seed=this%seed*14 - 2, dist='normal')
            topz = gauss_filt(topz, this%salt_top_smooth)
            topz = rescale(topz, [0.0, this%salt_top_amp])

            do isalt = 1, this%nsalt

                x1 = random(this%nstem + 1, seed=this%seed*13*isalt)
                x2 = random(this%nstem + 1, seed=this%seed*14*isalt)

                nd = nint((1.0 - rand(range=this%salt_top_max, seed=this%seed*15*isalt - 2))*this%n1 + ne1)

                x1 = ginterp([0.0, rand(range=[0.05*nd, 0.15*nd], seed=this%seed*15*isalt - 1), &
                    sort(random(this%nstem - 2, range=[0.2*nd, nd - 0.1*nd], seed=this%seed*15*isalt)), nd - 1.0], &
                    x1, linspace(0.0, nd + nex - 1.0, nd + nex), method='pchip')
                x2 = ginterp([0.0, rand(range=[0.05*nd, 0.15*nd], seed=this%seed*16*isalt - 1), &
                    sort(random(this%nstem - 2, range=[0.2*nd, nd - 0.1*nd], seed=this%seed*16*isalt)), nd - 1.0], &
                    x2, linspace(0.0, nd + nex - 1.0, nd + nex), method='pchip')

                x1 = rescale(x1, range=[0.1, 1.0]*rds(isalt)*this%n2)
                x2 = rescale(x2, range=[0.1, 1.0]*rds(isalt)*this%n2)
                x1 = pad(x1, [n1 - nd, 0], method=['edge', 'edge'])
                x2 = pad(x2, [n1 - nd, 0], method=['edge', 'edge'])

                !$omp parallel do private(i, j)
                do j = 1, n2
                    do i = n1 - nd - ceiling(maxval(topz)), n1 + nex
                        if (j >= gmax(isalt) - x1(i) .and. j <= gmax(isalt) + x2(i) .and. i >= n1 - nd - topz(j)) then
                            this%salt(i, j) = 1.0
                            if (i >= 1 .and. i <= n1) then
                                w(i, j) = this%salt_vel
                                f(i, j) = 0.0
                                fdip(i, j) = 0.0
                                fdisp(i, j) = 0.0
                                if (this%yn_rgt) then
                                    t(i, j) = 1.0
                                end if
                            end if
                        end if
                    end do
                end do
                !$omp end parallel do

            end do

            this%salt = this%salt(1:n1, 1:n2)

        end if

        ! Final processing

        ! After fault block shifting, the model is w, while RGT is t
        if (this%yn_facies) then
            m = w
        end if

        ww = w
        !$omp parallel do private(i, j)
        do j = 1, n2
            do i = 1, n1 - 1
                w(i, j) = (ww(i + 1, j) - ww(i, j))/(ww(i + 1, j) + ww(i, j))
            end do
        end do
        !$omp end parallel do
        where (f /= 0)
            w = w*0.5
        end where

        this%image = w(ne1 + 1:ne1 + this%n1, ne2 + 1:ne2 + this%n2)
        n1 = this%n1
        n2 = this%n2

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

        ! Output
        if (this%yn_fault) then
            this%fault = f(ne1 + 1:ne1 + this%n1, ne2 + 1:ne2 + this%n2)
            this%fault_dip = fdip(ne1 + 1:ne1 + this%n1, ne2 + 1:ne2 + this%n2)*const_rad2deg
            this%fault_disp = fdisp(ne1 + 1:ne1 + this%n1, ne2 + 1:ne2 + this%n2)
        end if

        if (this%yn_salt) then
            this%salt = this%salt(ne1 + 1:ne1 + this%n1, ne2 + 1:ne2 + this%n2)
        end if

        if (this%yn_rgt) then
            this%rgt = t(ne1 + 1:ne1 + this%n1, ne2 + 1:ne2 + this%n2)
            this%rgt = rescale(this%rgt, [0.0, 1.0])
        end if

        if (this%yn_facies) then
            this%facies = m(ne1 + 1:ne1 + this%n1, ne2 + 1:ne2 + this%n2)
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
                this%fault_disp = 0
            end if
            if (this%yn_facies) then
                this%facies = 0
            end if
            if (this%yn_salt) then
                this%salt = 0
            end if
        end if

    end subroutine generate_2d_geological_model

    subroutine generate_2d_unconformal_geological_model(this)

        type(rgm2_curved), intent(inout) :: this

        type(rgm2_curved), allocatable, dimension(:) :: g
        integer :: iconf, i, j
        type(meta_array1_real), allocatable, dimension(:) :: uff
        real, allocatable, dimension(:) :: ufz, aufz
        real, allocatable, dimension(:, :) :: rgt_above, rgt_below
        real :: tmin, tmax
        integer, allocatable, dimension(:) :: st

        allocate (g(1:this%unconf + 1))

        this%image = zeros(this%n1, this%n2)

        if (this%yn_fault) then
            this%fault = zeros(this%n1, this%n2)
            this%fault_dip = zeros(this%n1, this%n2)
            this%fault_disp = zeros(this%n1, this%n2)
        end if
        if (this%yn_rgt) then
            this%rgt = zeros(this%n1, this%n2)
        end if
        if (this%yn_facies) then
            this%facies = zeros(this%n1, this%n2)
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
                g(i)%refl_slope = -g(i)%refl_slope/4.0
                g(i)%refl_height_top = g(i)%refl_height_top/2.0
                g(i)%refl_slope_top = -g(i)%refl_slope_top/4.0
            end if

            if (i < this%unconf + 1 .and. this%unconf_nl == 0) then
                g(i)%unconf_nl = 0
            else
                g(i)%unconf_nl = g(i)%nl
            end if

            call generate_2d_geological_model(g(i))

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

            uff(i)%array = random(this%n2, seed=g(i)%seed*41*i)
            uff(i)%array = gauss_filt(uff(i)%array, this%n2*0.2)
            uff(i)%array = rescale(uff(i)%array, &
                [0.0, rand(range=this%unconf_amp, seed=g(i)%seed*51*i)*this%n1]) + ufz(i)*this%n1

        end do

        ! Merge sedimentary units
        this%image = g(this%unconf + 1)%image
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
                this%image(:, j) = &
                    taper(this%image(:, j), protect=[nint(uff(iconf)%array(j)) + 5, this%n1], len=[10, 0]) &
                    + taper(g(iconf)%image(:, j), protect=[1, nint(uff(iconf)%array(j)) - 5], len=[0, 10])

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
                            this%facies(i, j) = g(iconf)%facies(i, j)
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

    end subroutine generate_2d_unconformal_geological_model

end module geological_model_2d_curved

