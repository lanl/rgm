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


!
! Thrust faults of rgm2_curved and rgm3_curved: fault-bend folds,
! fault-propagation folds, imbricate fans, duplexes and pop-ups.
!
! A thrust is a fault path in vertical sections along its transport direction
! (xi; z positive downward). The lower flat (decollement) and the optional
! upper flat (roof) follow horizons of the pre-thrust model, and the ramp cuts
! up-section at the ramp dip; the bends are rounded. Behind the fault tip, the
! hanging wall moves along flow lines: fault-parallel flow by inclined shear
! along the bisector of the flat-ramp bend, where the flow line at distance d
! from the fault is the fault smoothed by a Gaussian of width widen*d, so the
! fold hinges widen away from the fault. Ahead of the tip, a trishear zone
! (Erslev, 1991; Zehnder and Allmendinger, 2000) distributes the shear while
! the tip propagates along the path at P/S times the slip rate. Each grid
! point is traced backward in time to its source point (gather), so repeated
! sections and overturned limbs need no special treatment. In 3D, each
! section is deformed in plane strain, while the ramp trace, the slip and the
! horizons vary from section to section.
!
module geological_model_thrust

    use libflit
    use geological_model_utility
    use geological_model_realism

    implicit none

    !
    !> A thrust fault and, once built for a model, its path tables in nsec
    !> sections along strike (nsec = 1 in 2D)
    !
    type thrust_fault

        !> Sense relative to the transport direction of the thrust system
        !> (1 = forethrust, -1 = backthrust)
        integer :: sense = 1
        !> Ramp dip (radians), maximum slip, propagation-to-slip ratio, tangent
        !> of the half apical angle of the trishear zone, arc length of the
        !> rounded bends, and growth of the hinge width with distance from the fault
        real :: dip = 0.5, slip = 0.0, ps = 2.0, m = 0.5, bend = 45.0, widen = 0.3
        !> Depths of the decollement and of the ramp top (tip or roof) below the
        !> model top, and position of the ramp top along the transport direction
        !> relative to the model center
        real :: zdec = 0.0, ztop = 0.0, xtop = 0.0
        !> Whether the ramp flattens into a roof (upper flat)
        logical :: roof = .false.
        !> 3D: maximum strike deviation of the ramp trace (radians), and half
        !> length and center of the slip patch along strike (radius = 0: through-going)
        real :: dstrike = 0.0, radius = 0.0, center = 0.0
        integer :: seed = -1

        !> Frame: transport azimuth (radians, from x2 toward x3) and origin
        real :: azimuth = 0.0, c2 = 0.0, c3 = 0.0
        !> Inclined shear: angle from the vertical, its tangent and cosine
        real :: alpha = 0.0, ta = 0.0, ca = 1.0
        !> Tables: g = xi - ta*z on [g0, g0 + (ng - 1)*dg]; arc length on uniform
        !> grids of spacing du (fault) and da (flow lines); smoothing levels sig;
        !> sections along strike at eta0 + (is - 1)*deta; visible separation
        !> along the fault at spacing dsep; time step ds
        integer :: ng = 0, na = 0, nu = 0, nlev = 1, nsec = 1, nsep = 0
        real :: g0 = 0.0, dg = 0.5, da = 0.5, du = 0.5, eta0 = 0.0, deta = 1.0, dsep = 1.0, ds = 0.5
        real, allocatable, dimension(:) :: sig
        real, allocatable, dimension(:, :) :: z0, th, ut, gu, sep
        real, allocatable, dimension(:, :, :) :: zt, arc, ga
        real, allocatable, dimension(:) :: sl, lf, l0, pse, gl0, dev, ddev
        logical, allocatable, dimension(:) :: straight

    end type thrust_fault

    private
    public :: thrust_fault
    public :: thrust_system
    public :: thrust_padding
    public :: thrust_apply_2d
    public :: thrust_apply_3d

contains

    !
    !> Linear interpolation in a table sampled at x0 + (i - 1)*dx, clamped at the ends
    !
    pure function lin(a, x0, dx, x) result(v)

        real, dimension(:), intent(in) :: a
        real, intent(in) :: x0, dx, x
        real :: v

        real :: f, w
        integer :: i, n

        n = size(a)
        f = (x - x0)/dx
        if (f <= 0.0) then
            v = a(1)
        else if (f >= n - 1.0) then
            v = a(n)
        else
            i = int(f)
            w = f - i
            v = a(i + 1)*(1.0 - w) + a(i + 2)*w
        end if

    end function lin

    !
    !> Inverse of an increasing table y sampled at x0 + (i - 1)*dx, on the grid
    !> (i - 1)*dy, i = 1, ..., ny (clamped outside the range of y)
    !
    pure function inverse_table(y, x0, dx, dy, ny) result(x)

        real, dimension(:), intent(in) :: y
        real, intent(in) :: x0, dx, dy
        integer, intent(in) :: ny
        real, dimension(1:ny) :: x

        integer :: i, j, n
        real :: yv

        n = size(y)
        j = 1
        do i = 1, ny
            yv = (i - 1)*dy
            if (yv <= y(1)) then
                x(i) = x0
            else if (yv >= y(n)) then
                x(i) = x0 + (n - 1)*dx
            else
                do while (j < n - 1 .and. y(j + 1) < yv)
                    j = j + 1
                end do
                x(i) = x0 + (j - 1 + (yv - y(j))/max(y(j + 1) - y(j), float_tiny))*dx
            end if
        end do

    end function inverse_table

    !
    !> Polynomial smooth minimum and maximum with blending width k
    !
    elemental function smin(a, b, k) result(c)

        real, intent(in) :: a, b, k
        real :: c

        real :: h

        if (k <= 0) then
            c = min(a, b)
        else
            h = max(k - abs(a - b), 0.0)/k
            c = min(a, b) - h*h*k*0.25
        end if

    end function smin

    elemental function smax(a, b, k) result(c)

        real, intent(in) :: a, b, k
        real :: c

        c = -smin(-a, -b, k)

    end function smax

    !
    !> Gaussian smoothing with standard deviation sigma (in samples), with the
    !> ends extended by their edge values; a direct convolution for small
    !> sigma and three box filters for large sigma
    !
    function smooth_edge(a, sigma) result(b)

        real, dimension(:), intent(in) :: a
        real, intent(in) :: sigma
        real, allocatable, dimension(:) :: b

        integer :: n, r, p, i, j, pass
        real, allocatable, dimension(:) :: e, kern
        double precision, allocatable, dimension(:) :: c

        n = size(a)
        if (sigma <= 0 .or. n < 2) then
            b = a
            return
        end if

        if (sigma <= 6.0) then
            r = ceiling(3.0*sigma)
            kern = [(exp(-0.5*(i/sigma)**2), i = -r, r)]
            kern = kern/sum(kern)
            allocate (b(1:n))
            do i = 1, n
                b(i) = 0.0
                do j = -r, r
                    b(i) = b(i) + kern(j + r + 1)*a(min(max(i + j, 1), n))
                end do
            end do
        else
            ! Three boxes of width 2r + 1 have the variance ((2r + 1)^2 - 1)/4
            r = max(nint(0.5*(sqrt(4.0*sigma**2 + 1.0) - 1.0)), 1)
            p = 3*r + 1
            e = [(a(1), i = 1, p), a, (a(n), i = 1, p)]
            allocate (c(1:size(e) + 1))
            do pass = 1, 3
                c(1) = 0.0d0
                do i = 1, size(e)
                    c(i + 1) = c(i) + e(i)
                end do
                do i = 1 + r, size(e) - r
                    e(i) = real((c(i + r + 1) - c(i - r))/(2*r + 1.0d0))
                end do
            end do
            b = e(p + 1:p + n)
        end if

    end function smooth_edge

    !
    !> Depth on the uniform g grid of a curve sampled at (xh(i), zh(i)), xh
    !> increasing. A horizon cut by an earlier thrust can jump; g is kept
    !> increasing, so a jump collapses into a near-vertical segment.
    !
    subroutine curve_in_g(xh, zh, ta, g0, dg, ng, zg)

        real, dimension(:), intent(in) :: xh, zh
        real, intent(in) :: ta, g0, dg
        integer, intent(in) :: ng
        real, allocatable, dimension(:), intent(out) :: zg

        real, allocatable, dimension(:) :: gh
        integer :: n, i, j
        real :: g

        n = size(xh)
        gh = xh - ta*zh
        do i = 2, n
            gh(i) = max(gh(i), gh(i - 1) + 1.0e-3)
        end do

        allocate (zg(1:ng))
        i = 1
        do j = 1, ng
            g = g0 + (j - 1)*dg
            if (g <= gh(1)) then
                zg(j) = zh(1)
            else if (g >= gh(n)) then
                zg(j) = zh(n)
            else
                do while (i < n - 1 .and. gh(i + 1) < g)
                    i = i + 1
                end do
                zg(j) = zh(i) + (zh(i + 1) - zh(i))*(g - gh(i))/(gh(i + 1) - gh(i))
            end if
        end do

    end subroutine curve_in_g

    !
    !> Deepest depth (0-based, fractional) where a column c, which increases with
    !> depth, crosses tau; clamped to the column ends
    !
    pure function column_crossing(c, tau) result(h)

        real, dimension(:), intent(in) :: c
        real, intent(in) :: tau
        real :: h

        integer :: i, n

        n = size(c)
        do i = n - 1, 1, -1
            if (c(i) <= tau .and. c(i + 1) > tau) then
                h = (i - 1.0) + (tau - c(i))/(c(i + 1) - c(i))
                return
            end if
        end do
        if (tau <= minval(c)) then
            h = 0.0
        else
            h = n - 1.0
        end if

    end function column_crossing

    !
    !> Value of a column at the fractional 0-based depth z
    !
    pure function column_value(c, z) result(v)

        real, dimension(:), intent(in) :: c
        real, intent(in) :: z
        real :: v

        v = lin(c, 0.0, 1.0, z)

    end function column_value

    !
    !> Draw the thrusts of a thrust system. The thrusts share a decollement,
    !> and with probability roof also a roof (upper flat), which makes a
    !> fault-bend fold for one thrust and a duplex for several. A thrust after
    !> the first is a backthrust with probability back, placed in the hanging
    !> wall of the previous forethrust (a pop-up). For the forward (piggyback)
    !> sequence, the forethrusts form from the hinterland toward the foreland,
    !> and each is placed back by the slip of the younger thrusts that carry it.
    !> Depths and positions are relative to the model top and center; nlat is
    !> the model size along the transport direction.
    !
    function thrust_system(nt, n1, nlat, seed, dip, disp, zdec, ztip, bend, widen, ps, apex, &
            roof, zroof, back, sequence) result(th)

        integer, intent(in) :: nt, n1, seed
        real, intent(in) :: nlat
        real, dimension(1:2), intent(in) :: dip, disp, zdec, ztip, bend, widen, ps, apex, zroof
        real, intent(in) :: roof, back
        character(len=*), intent(in) :: sequence
        type(thrust_fault), allocatable, dimension(:) :: th

        real, allocatable, dimension(:) :: vdip, vdisp, vtip, vbend, vwiden, vps, vapex, rb, pos, rs, lr
        real, dimension(1:2) :: rdisp, rbend
        real :: zd, zr, carry, za, zc, xc
        logical :: yn_centered
        logical :: yn_roof
        integer :: i, j, nfore, ib
        integer, allocatable, dimension(:) :: order
        type(thrust_fault), allocatable, dimension(:) :: tmp

        allocate (th(1:nt))
        if (nt == 0) then
            return
        end if

        ! Slip and bend: [0, 0] = automatic, relative to the model size
        if (maxval(disp) == 0) then
            rdisp = [0.05, 0.15]*nlat
        else
            rdisp = disp
        end if
        if (maxval(bend) == 0) then
            rbend = [0.15, 0.3]*n1
        else
            rbend = bend
        end if

        vdip = random(nt, range=dip, seed=safe_seed(int(seed, 8)*139))*const_deg2rad
        vdisp = random(nt, range=rdisp, seed=safe_seed(int(seed, 8)*149))
        vtip = random(nt, range=ztip, seed=safe_seed(int(seed, 8)*151))*n1
        vbend = random(nt, range=rbend, seed=safe_seed(int(seed, 8)*157))
        vwiden = random(nt, range=widen, seed=safe_seed(int(seed, 8)*163))
        vps = random(nt, range=ps, seed=safe_seed(int(seed, 8)*167))
        vapex = random(nt, range=apex, seed=safe_seed(int(seed, 8)*173))*0.5*const_deg2rad
        rb = random(nt, range=[0.0, 1.0], seed=safe_seed(int(seed, 8)*179))
        rs = random(nt, range=[0.0, 1.0], seed=safe_seed(int(seed, 8)*197))

        ! Shared decollement and roof
        zd = rand(range=zdec, seed=safe_seed(int(seed, 8)*181))*n1
        yn_roof = rand(range=[0.0, 1.0], seed=safe_seed(int(seed, 8)*191)) < roof
        zr = min(rand(range=zroof, seed=safe_seed(int(seed, 8)*193))*n1, zd - 0.1*n1)

        do i = 1, nt
            th(i)%dip = vdip(i)
            th(i)%slip = vdisp(i)
            th(i)%ps = vps(i)
            th(i)%m = tan(vapex(i))
            th(i)%bend = vbend(i)
            th(i)%widen = vwiden(i)
            th(i)%zdec = zd
            th(i)%roof = yn_roof
            if (yn_roof) then
                th(i)%ztop = zr
            else
                ! The tip stays above the decollement
                th(i)%ztop = min(vtip(i), zd - 0.1*n1)
            end if
            th(i)%seed = safe_seed(int(seed, 8)*199 + i)
        end do

        ! Backthrusts
        do i = 2, nt
            if (rb(i) < back) then
                th(i)%sense = -1
            end if
        end do
        nfore = count(th(:)%sense == 1)

        ! Horizontal ramp lengths
        allocate (lr(1:nt))
        do i = 1, nt
            lr(i) = (th(i)%zdec - th(i)%ztop)/tan(th(i)%dip)
        end do

        ! Forethrusts from the hinterland toward the foreland. Each ramp is
        ! anchored a third of the way down its depth range within the model, at
        ! a position spread over the model; the ramps of a duplex are spaced by
        ! a fraction of their horizontal lengths, so that the horses stack
        allocate (pos(1:nfore))
        if (yn_roof .and. nfore > 1) then
            j = 0
            do i = 1, nt
                if (th(i)%sense == 1) then
                    j = j + 1
                    if (j == 1) then
                        pos(j) = 0.0
                    else
                        pos(j) = pos(j - 1) + (0.3 + 0.3*rs(i))*lr(i)
                    end if
                end if
            end do
            pos = pos - 0.5*(pos(1) + pos(nfore))
        else if (nfore == 1) then
            pos = [rand(range=[-0.15, 0.15]*nlat, seed=safe_seed(int(seed, 8)*211))]
        else
            pos = sort(random(nfore, range=[-0.3, 0.3]*nlat, seed=safe_seed(int(seed, 8)*211), &
                spacing=0.6*0.6*nlat/nfore))
        end if

        j = 0
        do i = 1, nt
            if (th(i)%sense == 1) then
                j = j + 1
                za = max(th(i)%ztop, 0.0) + (min(th(i)%zdec, real(n1)) - max(th(i)%ztop, 0.0))/3.0
                th(i)%xtop = pos(j) + (za - th(i)%ztop)/tan(th(i)%dip)
            end if
        end do

        ! A backthrust rises toward the hinterland from a point at depth on the
        ! ramp of the previous forethrust, so the block between them pops up;
        ! the first such pair is centered at the anchor of the forethrust
        j = 0
        ib = 0
        do i = 1, nt
            if (th(i)%sense == 1) then
                j = j + 1
                ib = i
                yn_centered = .false.
                cycle
            end if
            zc = th(ib)%zdec - (0.2 + 0.3*rs(i))*(th(ib)%zdec - max(th(ib)%ztop, 0.0))
            xc = th(ib)%xtop - (zc - th(ib)%ztop)/tan(th(ib)%dip)
            th(i)%ztop = min(th(i)%ztop, zc - 0.1*n1)
            th(i)%xtop = xc - (zc - th(i)%ztop)/tan(th(i)%dip)
            if (.not. yn_centered) then
                xc = pos(j) - 0.5*(th(i)%xtop + th(ib)%xtop)
                th(ib)%xtop = th(ib)%xtop + xc
                th(i)%xtop = th(i)%xtop + xc
                yn_centered = .true.
            end if
        end do

        ! Piggyback: older thrusts are carried toward the foreland by the
        ! younger forethrusts, so they are placed back by their slips
        if (sequence /= 'backward' .and. sequence /= 'random') then
            do i = 1, nt
                carry = 0.0
                do j = i + 1, nt
                    if (th(j)%sense == 1) then
                        carry = carry + th(j)%slip*cos(th(j)%dip)
                    end if
                end do
                th(i)%xtop = th(i)%xtop - carry
            end do
        end if

        ! A backthrust is positioned along its own (opposite) transport direction
        do i = 1, nt
            if (th(i)%sense == -1) then
                th(i)%xtop = -th(i)%xtop
            end if
        end do

        ! Sequence of activity
        select case (sequence)
            case ('backward')
                ! Break-back: the forethrusts form from the foreland toward the hinterland
                allocate (order(1:nt))
                j = 0
                do i = nt, 1, -1
                    if (th(i)%sense == 1) then
                        j = j + 1
                        order(j) = i
                    end if
                end do
                do i = 1, nt
                    if (th(i)%sense == -1) then
                        j = j + 1
                        order(j) = i
                    end if
                end do
                tmp = th(order)
                th = tmp
            case ('random')
                order = nint(random_permute(real([(i, i = 1, nt)]), seed=safe_seed(int(seed, 8)*223)))
                tmp = th(order)
                th = tmp
        end select

    end function thrust_system

    !
    !> Padding needed by the thrusts: the hanging wall is lifted from below and
    !> moved from the hinterland
    !
    pure subroutine thrust_padding(th, pe1, pe2)

        type(thrust_fault), dimension(:), intent(in) :: th
        integer, intent(out) :: pe1, pe2

        integer :: i
        real :: up, lat

        up = 0.0
        lat = 0.0
        do i = 1, size(th)
            up = up + th(i)%slip*sin(th(i)%dip)
            lat = lat + th(i)%slip
        end do
        pe1 = ceiling(up) + 2
        pe2 = ceiling(lat) + 2

    end subroutine thrust_padding

    !
    !> Margins of the thrust padding that thrust k does not need to compute: the
    !> younger thrusts read only the part of the model that their slips reach
    !
    pure subroutine unused_padding(th, k, u1, u2)

        type(thrust_fault), dimension(:), intent(in) :: th
        integer, intent(in) :: k
        integer, intent(out) :: u1, u2

        integer :: pe1, pe2, j
        real :: up, lat

        call thrust_padding(th, pe1, pe2)
        up = 0.0
        lat = 0.0
        do j = k + 1, size(th)
            up = up + th(j)%slip*sin(th(j)%dip)
            lat = lat + th(j)%slip
        end do
        u1 = max(pe1 - ceiling(up) - 2, 0)
        u2 = max(pe2 - ceiling(lat) - 2, 0)

    end subroutine unused_padding

    !
    !> Smoothing levels, table sizes and allocation of the tables for a model
    !> of depth n1 and a section frame covering xi in [xmin, xmax]
    !
    subroutine thrust_init(t, n1, xmin, xmax, nsec)

        type(thrust_fault), intent(inout) :: t
        integer, intent(in) :: n1, nsec
        real, intent(in) :: xmin, xmax

        real, dimension(1:13), parameter :: levels = [0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 16.0, 24.0, 32.0, 48.0, 64.0, 96.0, 128.0]
        real :: smax_, margin, gmin, gmax
        integer :: i

        t%alpha = 0.5*t%dip
        t%ta = tan(t%alpha)
        t%ca = cos(t%alpha)

        ! Smoothing levels up to the width at the distance n1 from the fault
        smax_ = min(t%widen*n1, levels(13))
        if (smax_ <= 0.5) then
            t%sig = [0.0]
        else
            t%sig = [0.0]
            do i = 2, 13
                if (levels(i) < smax_) then
                    t%sig = [t%sig, levels(i)]
                end if
            end do
            t%sig = [t%sig, smax_]
        end if
        t%nlev = size(t%sig)

        ! g range: the frame, the shear offset, the slip, and the smoothing support
        margin = 4.0*t%sig(t%nlev) + t%bend + 50.0
        gmin = xmin - t%ta*(n1 + 0.0) - t%slip - margin
        gmax = xmax + t%ta*(n1 + 0.0) + margin
        t%dg = 0.5
        t%g0 = gmin
        t%ng = ceiling((gmax - gmin)/t%dg) + 1
        t%du = 0.5
        t%da = 0.5
        t%nu = ceiling(1.5*(gmax - gmin)/t%du) + 2
        t%na = t%nu
        t%dsep = 1.0
        t%nsep = ceiling(1.5*(gmax - gmin)/t%dsep) + 2
        t%nsec = nsec

        allocate (t%z0(1:t%ng, 1:nsec), t%th(1:t%ng, 1:nsec), t%ut(1:t%ng, 1:nsec))
        allocate (t%gu(1:t%nu, 1:nsec), t%sep(1:t%nsep, 1:nsec))
        allocate (t%zt(1:t%ng, 1:t%nlev, 1:nsec), t%arc(1:t%ng, 1:t%nlev, 1:nsec), t%ga(1:t%na, 1:t%nlev, 1:nsec))
        allocate (t%sl(1:nsec), t%lf(1:nsec), t%l0(1:nsec), t%pse(1:nsec), t%gl0(1:nsec), t%straight(1:nsec))
        if (.not. allocated(t%dev)) then
            allocate (t%dev(1:nsec), t%ddev(1:nsec))
            t%dev = 0.0
            t%ddev = 0.0
        end if
        t%sep = 0.0

    end subroutine thrust_init

    !
    !> Path tables of section is: the decollement and roof horizons hdec and hroof
    !> sampled at the uniform positions xh, the ramp through the anchor (xa, za)
    !> at the ramp dip, and the slip of the section
    !
    subroutine thrust_section(t, is, xh, hdec, hroof, xa, za, slip)

        type(thrust_fault), intent(inout) :: t
        integer, intent(in) :: is
        real, dimension(:), intent(in) :: xh, hdec, hroof
        real, intent(in) :: xa, za, slip

        real, allocatable, dimension(:) :: zdec, zroof, zramp, z, dz, dxx, dudg, zs
        real :: tt, kb, lfoot, lfc, l0, dmin, dmax, a
        integer :: i, lev, n

        n = t%ng
        tt = tan(t%dip)
        kb = 0.5*t%bend*tt

        call curve_in_g(xh, hdec, t%ta, t%g0, t%dg, n, zdec)
        allocate (zramp(1:n))
        do i = 1, n
            zramp(i) = (za - (t%g0 + (i - 1)*t%dg - xa)*tt)/(1.0 + t%ta*tt)
        end do
        z = smin(zramp, zdec, kb)
        if (t%roof) then
            call curve_in_g(xh, hroof, t%ta, t%g0, t%dg, n, zroof)
            z = smax(z, zroof, kb)
        end if
        t%z0(:, is) = z

        ! Dip and arc length along the fault
        allocate (dz(1:n), dxx(1:n), dudg(1:n))
        dz(1) = (z(2) - z(1))/t%dg
        dz(n) = (z(n) - z(n - 1))/t%dg
        do i = 2, n - 1
            dz(i) = (z(i + 1) - z(i - 1))/(2.0*t%dg)
        end do
        dxx = 1.0 + t%ta*dz
        t%th(:, is) = atan2(-dz, dxx)
        dudg = sqrt(dxx**2 + dz**2)
        t%ut(1, is) = 0.0
        do i = 2, n
            t%ut(i, is) = t%ut(i - 1, is) + 0.5*(dudg(i) + dudg(i - 1))*t%dg
        end do
        t%gu(:, is) = inverse_table(t%ut(:, is), t%g0, t%dg, t%du, t%nu)

        ! Flow lines: the fault smoothed at each level, their arc lengths and inverses
        do lev = 1, t%nlev
            zs = smooth_edge(z, t%sig(lev)/t%dg)
            t%zt(:, lev, is) = zs
            dz(1) = (zs(2) - zs(1))/t%dg
            dz(n) = (zs(n) - zs(n - 1))/t%dg
            do i = 2, n - 1
                dz(i) = (zs(i + 1) - zs(i - 1))/(2.0*t%dg)
            end do
            dudg = sqrt((1.0 + t%ta*dz)**2 + dz**2)
            t%arc(1, lev, is) = 0.0
            do i = 2, n
                t%arc(i, lev, is) = t%arc(i - 1, lev, is) + 0.5*(dudg(i) + dudg(i - 1))*t%dg
            end do
            t%ga(:, lev, is) = inverse_table(t%arc(:, lev, is), t%g0, t%dg, t%da, t%na)
        end do

        ! Ramp foot: where the ramp leaves the decollement
        lfoot = t%ut(1, is)
        do i = 1, n
            if (zramp(i) <= zdec(i)) then
                lfoot = t%ut(i, is)
                exit
            end if
        end do

        ! Final tip at the maximum slip (the anchor; a roof thrust has no tip
        ! in the model), initial tip and the effective P/S
        if (t%roof) then
            lfc = t%ut(n, is)
        else
            lfc = lin(t%ut(:, is), t%g0, t%dg, xa - t%ta*za)
        end if
        l0 = max(lfc - t%ps*t%slip, lfoot)
        t%l0(is) = l0
        t%pse(is) = (lfc - l0)/max(t%slip, float_tiny)
        t%sl(is) = slip
        t%lf(is) = l0 + t%pse(is)*slip
        t%gl0(is) = lin(t%gu(:, is), 0.0, t%du, l0)

        ! Whether the tip propagates along a straight part of the fault
        dmin = float_huge
        dmax = -float_huge
        do i = 0, 20
            a = lin(t%th(:, is), t%g0, t%dg, lin(t%gu(:, is), 0.0, t%du, l0 + (t%lf(is) - l0)*i/20.0))
            dmin = min(dmin, a)
            dmax = max(dmax, a)
        end do
        t%straight(is) = dmax - dmin < 0.01

    end subroutine thrust_section

    !
    !> Smoothing level and weight for the distance d from the fault
    !
    pure subroutine flow_level(t, d, k, w)

        type(thrust_fault), intent(in) :: t
        real, intent(in) :: d
        integer, intent(out) :: k
        real, intent(out) :: w

        real :: s

        if (t%nlev == 1) then
            k = 1
            w = 0.0
            return
        end if
        s = min(max(t%widen*d, 0.0), t%sig(t%nlev))
        k = 1
        do while (k < t%nlev - 1)
            if (s < t%sig(k + 1)) then
                exit
            end if
            k = k + 1
        end do
        w = min(max((s - t%sig(k))/(t%sig(k + 1) - t%sig(k)), 0.0), 1.0)

    end subroutine flow_level

    !
    !> Depth of the flow line of level (k, w) at g
    !
    pure function flow_depth(t, is, g, k, w) result(z)

        type(thrust_fault), intent(in) :: t
        integer, intent(in) :: is, k
        real, intent(in) :: g, w
        real :: z

        if (t%nlev == 1) then
            z = lin(t%zt(:, 1, is), t%g0, t%dg, g)
        else
            z = (1.0 - w)*lin(t%zt(:, k, is), t%g0, t%dg, g) + w*lin(t%zt(:, k + 1, is), t%g0, t%dg, g)
        end if

    end function flow_depth

    !
    !> Flow-line coordinates of the point (x, z): g, the distance d along the
    !> shear direction (> 0 in the hanging wall) and the smoothing level
    !
    pure subroutine flow_coord(t, is, x, z, g, d, k, w)

        type(thrust_fault), intent(in) :: t
        integer, intent(in) :: is
        real, intent(in) :: x, z
        real, intent(out) :: g, d, w
        integer, intent(out) :: k

        integer :: it

        g = x - t%ta*z
        d = (lin(t%z0(:, is), t%g0, t%dg, g) - z)/t%ca
        call flow_level(t, d, k, w)
        if (t%nlev > 1) then
            do it = 1, 4
                d = (flow_depth(t, is, g, k, w) - z)/t%ca
                call flow_level(t, d, k, w)
            end do
        end if

    end subroutine flow_coord

    !
    !> Move the point (g, d) back along its flow line by the arc length ds
    !
    pure subroutine flow_move(t, is, g, d, k, w, ds, x, z)

        type(thrust_fault), intent(in) :: t
        integer, intent(in) :: is, k
        real, intent(in) :: g, d, w, ds
        real, intent(out) :: x, z

        real :: a0, a1, gn

        if (t%nlev == 1) then
            a0 = lin(t%arc(:, 1, is), t%g0, t%dg, g) - ds
            gn = lin(t%ga(:, 1, is), 0.0, t%da, a0)
        else
            a0 = lin(t%arc(:, k, is), t%g0, t%dg, g) - ds
            a1 = lin(t%arc(:, k + 1, is), t%g0, t%dg, g) - ds
            gn = (1.0 - w)*lin(t%ga(:, k, is), 0.0, t%da, a0) + w*lin(t%ga(:, k + 1, is), 0.0, t%da, a1)
        end if
        z = flow_depth(t, is, gn, k, w) - d*t%ca
        x = gn + t%ta*z

    end subroutine flow_move

    !
    !> Velocity ahead of the tip at the arc length ltip: rigid translation along
    !> the fault above the trishear zone, a smoothstep profile across the zone,
    !> and no motion below it
    !
    pure subroutine tip_velocity(t, is, x, z, ltip, vx, vz)

        type(thrust_fault), intent(in) :: t
        integer, intent(in) :: is
        real, intent(in) :: x, z, ltip
        real, intent(out) :: vx, vz

        real :: g, xt, zt, tt, tx, tz, nx, nz, xi, eta, vxi, veta, q

        g = lin(t%gu(:, is), 0.0, t%du, ltip)
        zt = lin(t%z0(:, is), t%g0, t%dg, g)
        xt = g + t%ta*zt
        tt = lin(t%th(:, is), t%g0, t%dg, g)
        tx = cos(tt)
        tz = -sin(tt)
        nx = -sin(tt)
        nz = -cos(tt)
        xi = (x - xt)*tx + (z - zt)*tz
        eta = (x - xt)*nx + (z - zt)*nz
        vxi = 0.0
        veta = 0.0
        if (xi > 0) then
            if (eta >= t%m*xi) then
                vxi = 1.0
            else if (eta > -t%m*xi) then
                q = 0.5*(eta/(t%m*xi) + 1.0)
                vxi = q*q*(3.0 - 2.0*q)
                veta = -3.0*t%m*q*q*(1.0 - q)**2
            end if
        else if (eta >= 0) then
            vxi = 1.0
        end if
        vx = vxi*tx + veta*nx
        vz = vxi*tz + veta*nz

    end subroutine tip_velocity

    !
    !> Source (pre-thrust) position (x, z) of the point (x0, z0) of section is,
    !> traced backward in time. A point behind the initial tip moves along its
    !> flow line by the whole slip at once (hanging wall) or not at all
    !> (footwall). When the tip propagates along a straight ramp, a point ahead
    !> of the tip that stays below the trishear zone does not move, and one that
    !> stays above it translates along the ramp. Other points take time steps,
    !> along their flow lines behind the moving tip and with the tip velocity
    !> ahead of it.
    !
    pure subroutine thrust_backward(t, is, x0, z0, x, z)

        type(thrust_fault), intent(in) :: t
        integer, intent(in) :: is
        real, intent(in) :: x0, z0
        real, intent(out) :: x, z

        real :: s, ds, g, d, w, ltip, gtip, vx, vz, xm, zm
        real :: gf, xt, zt, tt, tx, tz, xi, eta
        integer :: k

        x = x0
        z = z0
        s = t%sl(is)
        if (s <= 0) then
            return
        end if

        g = x - t%ta*z
        if (g < t%gl0(is)) then
            call flow_coord(t, is, x, z, g, d, k, w)
            if (d >= 0) then
                call flow_move(t, is, g, d, k, w, s, x, z)
            end if
            return
        end if

        if (t%straight(is)) then
            gf = lin(t%gu(:, is), 0.0, t%du, t%lf(is))
            if (g >= gf) then
                zt = lin(t%z0(:, is), t%g0, t%dg, gf)
                xt = gf + t%ta*zt
                tt = lin(t%th(:, is), t%g0, t%dg, gf)
                tx = cos(tt)
                tz = -sin(tt)
                xi = (x - xt)*tx + (z - zt)*tz
                eta = -(x - xt)*sin(tt) - (z - zt)*cos(tt)
                ! Below the zone at all times: the tip retreats by lf - l0
                if (eta < 0 .and. eta <= -t%m*max(xi + t%lf(is) - t%l0(is), 0.0)) then
                    return
                end if
                ! Above the zone at all times: the point moves along the ramp at
                ! the slip rate while the tip moves at P/S times it
                if (eta >= 0 .and. eta >= t%m*max(xi + max(t%pse(is) - 1.0, 0.0)*s, 0.0)) then
                    x = x - s*tx
                    z = z - s*tz
                    return
                end if
            end if
        end if

        ! Time steps. Behind the tip, a point on a flow line moves exactly by
        ! any step, so it advances until the retreating tip could reach it, and
        ! a footwall point stays until the tip retreats past it; ahead of the
        ! tip, the step grows with the distance from the tip, where the
        ! velocity varies slowly
        do while (s > 0)
            ltip = t%l0(is) + t%pse(is)*s
            gtip = lin(t%gu(:, is), 0.0, t%du, ltip)
            g = x - t%ta*z
            if (g < gtip) then
                call flow_coord(t, is, x, z, g, d, k, w)
                if (g < t%gl0(is)) then
                    ! Behind the initial tip: the rest at once
                    if (d >= 0) then
                        call flow_move(t, is, g, d, k, w, s, x, z)
                    end if
                    return
                end if
                if (d >= 0) then
                    if (t%pse(is) <= 1.0) then
                        call flow_move(t, is, g, d, k, w, s, x, z)
                        return
                    end if
                    ds = min(s, max(t%ds, 0.5*(gtip - g)/(t%pse(is) - 1.0)))
                    call flow_move(t, is, g, d, k, w, ds, x, z)
                else
                    ! The tip passes the point when l0 + pse*s = u(g)
                    ds = s - (lin(t%ut(:, is), t%g0, t%dg, g) - t%l0(is))/max(t%pse(is), float_tiny)
                    ds = min(s, max(ds, t%ds))
                end if
            else
                zt = lin(t%z0(:, is), t%g0, t%dg, gtip)
                xt = gtip + t%ta*zt
                ds = min(s, max(t%ds, min(0.15*sqrt((x - xt)**2 + (z - zt)**2), 8.0*t%ds)))
                call tip_velocity(t, is, x, z, ltip, vx, vz)
                xm = x - 0.5*ds*vx
                zm = z - 0.5*ds*vz
                call tip_velocity(t, is, xm, zm, t%l0(is) + t%pse(is)*(s - 0.5*ds), vx, vz)
                x = x - ds*vx
                z = z - ds*vz
            end if
            s = s - ds
        end do

    end subroutine thrust_backward

    !
    !> Arc length of the foot, perpendicular distance to the fault and local
    !> dip at the point (x, z) of section is, in the post-thrust frame
    !
    pure subroutine thrust_point(t, is, x, z, u, dperp, dip)

        type(thrust_fault), intent(in) :: t
        integer, intent(in) :: is
        real, intent(in) :: x, z
        real, intent(out) :: u, dperp, dip

        real :: g

        g = x - t%ta*z
        u = lin(t%ut(:, is), t%g0, t%dg, g)
        dip = lin(t%th(:, is), t%g0, t%dg, g)
        dperp = (lin(t%z0(:, is), t%g0, t%dg, g) - z)/t%ca*cos(dip - t%alpha)

    end subroutine thrust_point

    !
    !> Visibility of a fault along its length, from the separation sep (in grid
    !> points) of the samples spaced by dsep: 1 where the separation is at least
    !> 0.5 grid points, and 0 elsewhere. On flats, the separation is near zero
    !> but noisy, for instance where an older fault runs along the decollement,
    !> so gaps shorter than gmin samples are filled, and visible runs shorter
    !> than lmin samples are removed.
    !
    pure subroutine clean_visibility(sep, gmin, lmin)

        real, dimension(:), intent(inout) :: sep
        integer, intent(in) :: gmin, lmin

        logical, dimension(size(sep)) :: v
        integer :: n, i, j

        n = size(sep)
        v = sep >= 0.5

        ! Fill short gaps between visible runs
        i = 1
        do while (i <= n)
            if (v(i)) then
                i = i + 1
                cycle
            end if
            j = i
            do while (j < n)
                if (v(j + 1)) then
                    exit
                end if
                j = j + 1
            end do
            if (i > 1 .and. j < n .and. j - i + 1 < gmin) then
                v(i:j) = .true.
            end if
            i = j + 1
        end do

        ! Remove short visible runs
        i = 1
        do while (i <= n)
            if (.not. v(i)) then
                i = i + 1
                cycle
            end if
            j = i
            do while (j < n)
                if (.not. v(j + 1)) then
                    exit
                end if
                j = j + 1
            end do
            if (j - i + 1 < lmin) then
                v(i:j) = .false.
            end if
            i = j + 1
        end do

        sep = merge(1.0, 0.0, v)

    end subroutine clean_visibility

    !
    !> Points just above (a) and below (b) the fault of section is, at the arc
    !> lengths (i - 1)*dsep up to the final tip
    !
    subroutine side_points(t, is, xa, za, xb, zb, nn)

        type(thrust_fault), intent(in) :: t
        integer, intent(in) :: is
        real, allocatable, dimension(:), intent(out) :: xa, za, xb, zb
        integer, intent(out) :: nn

        integer :: i
        real :: g, x, z, a

        nn = min(t%nsep, floor(t%lf(is)/t%dsep) + 1)
        allocate (xa(1:nn), za(1:nn), xb(1:nn), zb(1:nn))
        do i = 1, nn
            g = lin(t%gu(:, is), 0.0, t%du, (i - 1)*t%dsep)
            z = lin(t%z0(:, is), t%g0, t%dg, g)
            x = g + t%ta*z
            a = lin(t%th(:, is), t%g0, t%dg, g)
            xa(i) = x - sin(a)
            za(i) = z - cos(a)
            xb(i) = x + sin(a)
            zb(i) = z + cos(a)
        end do

    end subroutine side_points

    !
    !> Bilinear and nearest sampling of a 2D array at 0-based fractional positions
    !
    pure function sample_2d(a, z, x) result(v)

        real, dimension(:, :), intent(in) :: a
        real, intent(in) :: z, x
        real :: v

        integer :: i, j, n1, n2
        real :: wz, wx, zc, xc

        n1 = size(a, 1)
        n2 = size(a, 2)
        zc = min(max(z, 0.0), n1 - 1.001)
        xc = min(max(x, 0.0), n2 - 1.001)
        i = int(zc) + 1
        j = int(xc) + 1
        wz = zc - (i - 1)
        wx = xc - (j - 1)
        v = a(i, j)*(1 - wz)*(1 - wx) + a(i + 1, j)*wz*(1 - wx) &
            + a(i, j + 1)*(1 - wz)*wx + a(i + 1, j + 1)*wz*wx

    end function sample_2d

    pure function nearest_2d(a, z, x) result(v)

        real, dimension(:, :), intent(in) :: a
        real, intent(in) :: z, x
        real :: v

        v = a(min(max(nint(z), 0), size(a, 1) - 1) + 1, min(max(nint(x), 0), size(a, 2) - 1) + 1)

    end function nearest_2d

    pure function sample_3d(a, z, x, y) result(v)

        real, dimension(:, :, :), intent(in) :: a
        real, intent(in) :: z, x, y
        real :: v

        integer :: i, j, k, n1, n2, n3
        real :: wz, wx, wy, zc, xc, yc

        n1 = size(a, 1)
        n2 = size(a, 2)
        n3 = size(a, 3)
        zc = min(max(z, 0.0), n1 - 1.001)
        xc = min(max(x, 0.0), n2 - 1.001)
        yc = min(max(y, 0.0), n3 - 1.001)
        i = int(zc) + 1
        j = int(xc) + 1
        k = int(yc) + 1
        wz = zc - (i - 1)
        wx = xc - (j - 1)
        wy = yc - (k - 1)
        v = (a(i, j, k)*(1 - wz)*(1 - wx) + a(i + 1, j, k)*wz*(1 - wx) &
            + a(i, j + 1, k)*(1 - wz)*wx + a(i + 1, j + 1, k)*wz*wx)*(1 - wy) &
            + (a(i, j, k + 1)*(1 - wz)*(1 - wx) + a(i + 1, j, k + 1)*wz*(1 - wx) &
            + a(i, j + 1, k + 1)*(1 - wz)*wx + a(i + 1, j + 1, k + 1)*wz*wx)*wy

    end function sample_3d

    pure function nearest_3d(a, z, x, y) result(v)

        real, dimension(:, :, :), intent(in) :: a
        real, intent(in) :: z, x, y
        real :: v

        v = a(min(max(nint(z), 0), size(a, 1) - 1) + 1, min(max(nint(x), 0), size(a, 2) - 1) + 1, &
            min(max(nint(y), 0), size(a, 3) - 1) + 1)

    end function nearest_3d

    !
    !> Bilinear sampling of p at a 0-based fractional position from the neighbors
    !> that belong to fault idv (id == idv); far where idv = 0 or no neighbor does
    !
    pure function sample_id_2d(p, id, idv, far, z, x) result(v)

        real, dimension(:, :), intent(in) :: p, id
        real, intent(in) :: idv, far, z, x
        real :: v

        integer :: i, j, a, b
        real :: wz, wx, zc, xc, w, sw

        v = far
        if (idv == 0) then
            return
        end if
        zc = min(max(z, 0.0), size(p, 1) - 1.001)
        xc = min(max(x, 0.0), size(p, 2) - 1.001)
        i = int(zc) + 1
        j = int(xc) + 1
        wz = zc - (i - 1)
        wx = xc - (j - 1)
        sw = 0.0
        v = 0.0
        do b = 0, 1
            do a = 0, 1
                if (id(i + a, j + b) == idv) then
                    w = merge(wz, 1.0 - wz, a == 1)*merge(wx, 1.0 - wx, b == 1)
                    v = v + w*p(i + a, j + b)
                    sw = sw + w
                end if
            end do
        end do
        if (sw > 0) then
            v = v/sw
        else
            v = far
        end if

    end function sample_id_2d

    pure function sample_id_3d(p, id, idv, far, z, x, y) result(v)

        real, dimension(:, :, :), intent(in) :: p, id
        real, intent(in) :: idv, far, z, x, y
        real :: v

        integer :: i, j, k, a, b, c
        real :: wz, wx, wy, zc, xc, yc, w, sw

        v = far
        if (idv == 0) then
            return
        end if
        zc = min(max(z, 0.0), size(p, 1) - 1.001)
        xc = min(max(x, 0.0), size(p, 2) - 1.001)
        yc = min(max(y, 0.0), size(p, 3) - 1.001)
        i = int(zc) + 1
        j = int(xc) + 1
        k = int(yc) + 1
        wz = zc - (i - 1)
        wx = xc - (j - 1)
        wy = yc - (k - 1)
        sw = 0.0
        v = 0.0
        do c = 0, 1
            do b = 0, 1
                do a = 0, 1
                    if (id(i + a, j + b, k + c) == idv) then
                        w = merge(wz, 1.0 - wz, a == 1)*merge(wx, 1.0 - wx, b == 1)*merge(wy, 1.0 - wy, c == 1)
                        v = v + w*p(i + a, j + b, k + c)
                        sw = sw + w
                    end if
                end do
            end do
        end do
        if (sw > 0) then
            v = v/sw
        else
            v = far
        end if

    end function sample_id_3d

    !
    !> Derivative of p along one axis at a grid point, from the neighbors that
    !> belong to the same fault; pm and pp are the neighbors, lm and lp whether
    !> they belong to it
    !
    pure function deriv_id(p0, pm, pp, lm, lp) result(g)

        real, intent(in) :: p0, pm, pp
        logical, intent(in) :: lm, lp
        real :: g

        if (lm .and. lp) then
            g = 0.5*(pp - pm)
        else if (lp) then
            g = pp - p0
        else if (lm) then
            g = p0 - pm
        else
            g = 0.0
        end if

    end function deriv_id

    !
    !> Keep a thrust label (index above nf) where the distance to its fault is
    !> below half the fault width. The distance is the signed distance phi
    !> divided by the norm of its gradient, because a younger thrust stretches
    !> or shortens phi when it folds the fault; the labels then have the same
    !> width everywhere.
    !
    subroutine thrust_rasterize_2d(phi, f, fdip, fdisp, nf, fwidth)

        real, dimension(:, :), intent(in) :: phi
        real, dimension(:, :), intent(inout) :: f, fdip, fdisp
        integer, intent(in) :: nf
        real, intent(in) :: fwidth

        logical, allocatable, dimension(:, :) :: keep
        integer :: n1, n2, i, j
        real :: idv, g1, g2, pm, pp
        logical :: lm, lp

        n1 = size(f, 1)
        n2 = size(f, 2)
        allocate (keep(1:n1, 1:n2))
        keep = .true.

        !$omp parallel do private(i, j, idv, g1, g2, pm, pp, lm, lp) collapse(2)
        do j = 1, n2
            do i = 1, n1
                idv = f(i, j)
                if (idv > nf) then
                    lm = .false.
                    lp = .false.
                    pm = 0.0
                    pp = 0.0
                    if (i > 1) then
                        lm = f(i - 1, j) == idv
                        pm = phi(i - 1, j)
                    end if
                    if (i < n1) then
                        lp = f(i + 1, j) == idv
                        pp = phi(i + 1, j)
                    end if
                    g1 = deriv_id(phi(i, j), pm, pp, lm, lp)
                    lm = .false.
                    lp = .false.
                    if (j > 1) then
                        lm = f(i, j - 1) == idv
                        pm = phi(i, j - 1)
                    end if
                    if (j < n2) then
                        lp = f(i, j + 1) == idv
                        pp = phi(i, j + 1)
                    end if
                    g2 = deriv_id(phi(i, j), pm, pp, lm, lp)
                    keep(i, j) = abs(phi(i, j)) < 0.5*fwidth*max(sqrt(g1**2 + g2**2), 0.25)
                end if
            end do
        end do
        !$omp end parallel do

        where (.not. keep)
            f = 0.0
            fdip = 0.0
            fdisp = 0.0
        end where

    end subroutine thrust_rasterize_2d

    subroutine thrust_rasterize_3d(phi, f, fdip, fstrike, frake, fdisp, nf, fwidth)

        real, dimension(:, :, :), intent(in) :: phi
        real, dimension(:, :, :), intent(inout) :: f, fdip, fstrike, frake, fdisp
        integer, intent(in) :: nf
        real, intent(in) :: fwidth

        logical, allocatable, dimension(:, :, :) :: keep
        integer :: n1, n2, n3, i, j, l
        real :: idv, g1, g2, g3, pm, pp
        logical :: lm, lp

        n1 = size(f, 1)
        n2 = size(f, 2)
        n3 = size(f, 3)
        allocate (keep(1:n1, 1:n2, 1:n3))
        keep = .true.

        !$omp parallel do private(i, j, l, idv, g1, g2, g3, pm, pp, lm, lp) collapse(3)
        do l = 1, n3
            do j = 1, n2
                do i = 1, n1
                    idv = f(i, j, l)
                    if (idv > nf) then
                        lm = .false.
                        lp = .false.
                        pm = 0.0
                        pp = 0.0
                        if (i > 1) then
                            lm = f(i - 1, j, l) == idv
                            pm = phi(i - 1, j, l)
                        end if
                        if (i < n1) then
                            lp = f(i + 1, j, l) == idv
                            pp = phi(i + 1, j, l)
                        end if
                        g1 = deriv_id(phi(i, j, l), pm, pp, lm, lp)
                        lm = .false.
                        lp = .false.
                        if (j > 1) then
                            lm = f(i, j - 1, l) == idv
                            pm = phi(i, j - 1, l)
                        end if
                        if (j < n2) then
                            lp = f(i, j + 1, l) == idv
                            pp = phi(i, j + 1, l)
                        end if
                        g2 = deriv_id(phi(i, j, l), pm, pp, lm, lp)
                        lm = .false.
                        lp = .false.
                        if (l > 1) then
                            lm = f(i, j, l - 1) == idv
                            pm = phi(i, j, l - 1)
                        end if
                        if (l < n3) then
                            lp = f(i, j, l + 1) == idv
                            pp = phi(i, j, l + 1)
                        end if
                        g3 = deriv_id(phi(i, j, l), pm, pp, lm, lp)
                        keep(i, j, l) = abs(phi(i, j, l)) < 0.5*fwidth*max(sqrt(g1**2 + g2**2 + g3**2), 0.25)
                    end if
                end do
            end do
        end do
        !$omp end parallel do

        where (.not. keep)
            f = 0.0
            fdip = 0.0
            fstrike = 0.0
            frake = 0.0
            fdisp = 0.0
        end where

    end subroutine thrust_rasterize_3d

    !
    !> Apply the thrusts to the padded 2D model arrays (velocity w, RGT t,
    !> facies m) and mark them in the fault arrays as faults nf + 1, nf + 2, ...
    !> The thrusts are applied in the order of the array (older first); ne1 is
    !> the padding above the model, to which the depths of the thrusts refer;
    !> th(:)%azimuth is 0 for transport toward +x2 and pi toward -x2.
    !
    subroutine thrust_apply_2d(th, ne1, fwidth, nf, w, t, m, f, fdip, fdisp, yn_facies)

        type(thrust_fault), dimension(:), intent(inout) :: th
        integer, intent(in) :: ne1, nf
        real, intent(in) :: fwidth
        real, dimension(:, :), intent(inout) :: w, t, f, fdip, fdisp
        real, allocatable, dimension(:, :), intent(inout) :: m
        logical, intent(in) :: yn_facies

        real, allocatable, dimension(:, :) :: ww, tt, mm, ff, ffdip, ffdisp, phi, pp
        real, allocatable, dimension(:) :: xh, hdec, hroof, cdec, croof, xa, za, xb, zb, sep
        real :: c2, dir, tau_dec, tau_roof, xi, z, xs, zs, u, dperp, dip, dtdz, x
        real :: xsa, zsa, xsb, zsb, xq, zq, idv, dband
        integer :: n1, n2, k, i, j, nh, nn, u1, u2
        logical :: yn_dec_flat

        n1 = size(w, 1)
        n2 = size(w, 2)
        c2 = 0.5*(n2 - 1.0)

        ! Half width of the band over which a fault travels as fields
        dband = 0.5*fwidth + 2.0
        phi = zeros(n1, n2) + dband

        ! Horizons followed by the flats, identified by the RGT of the pre-thrust
        ! model at the center column; a decollement below the padded model is
        ! horizontal
        yn_dec_flat = th(1)%zdec + ne1 >= n1 - 2.0
        tau_dec = column_value(t(:, nint(c2) + 1), th(1)%zdec + ne1)
        tau_roof = column_value(t(:, nint(c2) + 1), th(1)%ztop + ne1)

        do k = 1, size(th)

            dir = sign(1.0, cos(th(k)%azimuth))
            th(k)%c2 = c2
            call thrust_init(th(k), n1, -c2, c2, 1)

            ! Horizons in the current model, sampled along the section
            allocate (cdec(1:n2), croof(1:n2))
            do j = 1, n2
                if (yn_dec_flat) then
                    cdec(j) = th(k)%zdec + ne1
                else
                    cdec(j) = column_crossing(t(:, j), tau_dec)
                end if
                croof(j) = column_crossing(t(:, j), tau_roof)
            end do
            nh = ceiling((th(k)%ng*th(k)%dg + 2*th(k)%ta*n1 + 200.0)/0.5)
            allocate (xh(1:nh), hdec(1:nh), hroof(1:nh))
            do i = 1, nh
                xh(i) = th(k)%g0 - 100.0 + (i - 1)*0.5
                x = c2 + dir*xh(i)
                hdec(i) = lin(cdec, 0.0, 1.0, x)
                hroof(i) = lin(croof, 0.0, 1.0, x)
            end do
            call thrust_section(th(k), 1, xh, hdec, hroof, th(k)%xtop, th(k)%ztop + ne1, th(k)%slip)
            deallocate (cdec, croof, xh, hdec, hroof)

            ! Visible separation along the fault: the stratigraphic shift of the
            ! hanging wall relative to the footwall, from the RGT
            dtdz = mean(abs(t(2:n1, :) - t(1:n1 - 1, :))) + float_tiny
            call side_points(th(k), 1, xa, za, xb, zb, nn)
            allocate (sep(1:nn))
            !$omp parallel do private(i, xsa, zsa, xsb, zsb)
            do i = 1, nn
                call thrust_backward(th(k), 1, xa(i), za(i), xsa, zsa)
                call thrust_backward(th(k), 1, xb(i), zb(i), xsb, zsb)
                sep(i) = abs((sample_2d(t, zsa, c2 + dir*xsa) - sample_2d(t, za(i), c2 + dir*xa(i))) &
                    - (sample_2d(t, zsb, c2 + dir*xsb) - sample_2d(t, zb(i), c2 + dir*xb(i))))/dtdz
            end do
            !$omp end parallel do
            th(k)%sep(1:nn, 1) = sep
            call clean_visibility(th(k)%sep(1:nn, 1), 6, 10)
            deallocate (sep)

            ! Gather from the source points. The faults travel as fields over a band
            ! around each fault: f is the index of the nearest fault (fdip and fdisp
            ! its attributes), and phi is the signed distance to it, positive in its
            ! hanging wall; phi is interpolated, so that the labels keep their width
            ! when a younger thrust folds the fault. In the footwall next to this
            ! fault, the fields continue across this fault from its hanging wall, so
            ! that an older fault that ran along this fault (a shared flat) moves as
            ! a whole instead of leaving its lower half behind
            ww = w
            tt = t
            ff = f
            ffdip = fdip
            ffdisp = fdisp
            pp = phi
            if (yn_facies) then
                mm = m
            end if
            call unused_padding(th, k, u1, u2)
            !$omp parallel do private(i, j, xi, z, xs, zs, u, dperp, dip, xq, zq, idv) collapse(2)
            do j = 1 + u2, n2 - u2
                do i = 1 + u1, n1 - u1
                    xi = dir*(j - 1.0 - c2)
                    z = i - 1.0
                    call thrust_backward(th(k), 1, xi, z, xs, zs)
                    xs = c2 + dir*xs
                    w(i, j) = sample_2d(ww, zs, xs)
                    t(i, j) = sample_2d(tt, zs, xs)
                    if (yn_facies) then
                        m(i, j) = nearest_2d(mm, zs, xs)
                    end if
                    call thrust_point(th(k), 1, xi, z, u, dperp, dip)
                    if (dperp <= 0 .and. dperp > -dband .and. u <= th(k)%lf(1) .and. th(k)%sl(1) >= 0.5) then
                        call thrust_backward(th(k), 1, xi - (0.5 - dperp)*sin(dip), z - (0.5 - dperp)*cos(dip), xq, zq)
                        xq = c2 + dir*xq
                        idv = nearest_2d(ff, zq, xq)
                        phi(i, j) = max(sample_id_2d(pp, ff, idv, dband, zq, xq) - (0.5 - dperp), -dband)
                    else
                        xq = xs
                        zq = zs
                        idv = nearest_2d(ff, zq, xq)
                        phi(i, j) = sample_id_2d(pp, ff, idv, dband, zq, xq)
                    end if
                    f(i, j) = idv
                    fdip(i, j) = nearest_2d(ffdip, zq, xq)
                    fdisp(i, j) = nearest_2d(ffdisp, zq, xq)
                    ! This fault where its offset is visible: within half its width, and
                    ! within the band where it is nearer than the older faults
                    if (abs(dperp) < dband .and. u <= th(k)%lf(1) .and. th(k)%sl(1) >= 0.5) then
                        if ((abs(dperp) < 0.5*fwidth .or. abs(dperp) < abs(phi(i, j))) &
                                .and. lin(th(k)%sep(:, 1), 0.0, th(k)%dsep, u) >= 0.5) then
                            f(i, j) = nf + k
                            ! Dip in (0, pi): below pi/2 the fault dips toward -x2
                            fdip(i, j) = modulo(dir*dip, real(const_pi))
                            fdisp(i, j) = th(k)%sl(1)
                            phi(i, j) = dperp
                        end if
                    end if
                end do
            end do
            !$omp end parallel do

        end do

        ! The labels: within half the fault width of the faults
        call thrust_rasterize_2d(phi, f, fdip, fdisp, nf, fwidth)

    end subroutine thrust_apply_2d

    !
    !> Apply the thrusts to the padded 3D model arrays, as thrust_apply_2d. Each
    !> thrust is computed in vertical sections along its transport azimuth,
    !> spaced by 2 grid points along strike, and a point takes the linear blend
    !> of the two nearest sections. th(:)%azimuth must be set.
    !
    subroutine thrust_apply_3d(th, ne1, fwidth, nf, w, t, m, f, fdip, fstrike, frake, fdisp, yn_facies)

        type(thrust_fault), dimension(:), intent(inout) :: th
        integer, intent(in) :: ne1, nf
        real, intent(in) :: fwidth
        real, dimension(:, :, :), intent(inout) :: w, t, f, fdip, fstrike, frake, fdisp
        real, allocatable, dimension(:, :, :), intent(inout) :: m
        logical, intent(in) :: yn_facies

        real, allocatable, dimension(:, :, :) :: ww, tt, mm, ff, ffdip, ffstrike, ffrake, ffdisp, phi, pp
        real, allocatable, dimension(:, :) :: cdec, croof
        real, allocatable, dimension(:) :: xh, hdec, hroof, xa, za, xb, zb, dphi
        real :: c2, c3, cp, sp, tau_dec, tau_roof, xi, eta, z, xs, zs, xs2, zs2, e, ws
        real :: u, dperp, dip, u2, dperp2, dip2, sepv, dtdz, emin, emax, xmax_, x2, x3, slip, taper
        real :: xsa, zsa, xsb, zsb, sl, dl, beta, tdip, idv, dband
        integer :: n1, n2, n3, k, i, j, l, is, nh, nn, isec, i0, m1, m2
        logical :: yn_dec_flat

        n1 = size(w, 1)
        n2 = size(w, 2)
        n3 = size(w, 3)
        c2 = 0.5*(n2 - 1.0)
        c3 = 0.5*(n3 - 1.0)

        ! Half width of the band over which a fault travels as fields
        dband = 0.5*fwidth + 2.0
        phi = zeros(n1, n2, n3) + dband

        yn_dec_flat = th(1)%zdec + ne1 >= n1 - 2.0
        tau_dec = column_value(t(:, nint(c2) + 1, nint(c3) + 1), th(1)%zdec + ne1)
        tau_roof = column_value(t(:, nint(c2) + 1, nint(c3) + 1), th(1)%ztop + ne1)

        allocate (cdec(1:n2, 1:n3), croof(1:n2, 1:n3))

        do k = 1, size(th)

            th(k)%c2 = c2
            th(k)%c3 = c3
            cp = cos(th(k)%azimuth)
            sp = sin(th(k)%azimuth)

            ! Sections along strike covering the padded model
            xmax_ = c2*abs(cp) + c3*abs(sp)
            emax = c2*abs(sp) + c3*abs(cp)
            emin = -emax
            th(k)%deta = 2.0
            th(k)%eta0 = emin - th(k)%deta
            isec = ceiling((emax - emin)/th(k)%deta) + 3
            if (allocated(th(k)%dev)) then
                deallocate (th(k)%dev, th(k)%ddev)
            end if
            allocate (th(k)%dev(1:isec), th(k)%ddev(1:isec))

            ! Ramp trace: a smooth strike deviation, integrated along strike and
            ! zero at the model center
            if (th(k)%dstrike > 0) then
                dphi = smooth_noise_1d(isec, 0.25*(emax - emin)/th(k)%deta, th(k)%seed)
                dphi = dphi/(maxval(abs(dphi)) + float_tiny)*th(k)%dstrike
            else
                dphi = zeros(isec)
            end if
            th(k)%ddev = tan(dphi)
            th(k)%dev(1) = 0.0
            do is = 2, isec
                th(k)%dev(is) = th(k)%dev(is - 1) + 0.5*(th(k)%ddev(is) + th(k)%ddev(is - 1))*th(k)%deta
            end do
            i0 = nint(-th(k)%eta0/th(k)%deta) + 1
            th(k)%dev = th(k)%dev - th(k)%dev(min(max(i0, 1), isec))

            call thrust_init(th(k), n1, -xmax_, xmax_, isec)

            ! Horizons in the current model
            !$omp parallel do private(i, j) collapse(2)
            do j = 1, n3
                do i = 1, n2
                    if (yn_dec_flat) then
                        cdec(i, j) = th(k)%zdec + ne1
                    else
                        cdec(i, j) = column_crossing(t(:, i, j), tau_dec)
                    end if
                    croof(i, j) = column_crossing(t(:, i, j), tau_roof)
                end do
            end do
            !$omp end parallel do

            ! Section tables
            nh = ceiling((th(k)%ng*th(k)%dg + 2*th(k)%ta*n1 + 200.0)/0.5)
            !$omp parallel do private(is, eta, xh, hdec, hroof, i, x2, x3, taper, slip)
            do is = 1, isec
                eta = th(k)%eta0 + (is - 1)*th(k)%deta
                allocate (xh(1:nh), hdec(1:nh), hroof(1:nh))
                do i = 1, nh
                    xh(i) = th(k)%g0 - 100.0 + (i - 1)*0.5
                    x2 = c2 + xh(i)*cp - eta*sp
                    x3 = c3 + xh(i)*sp + eta*cp
                    hdec(i) = sample_2d(cdec, x2, x3)
                    hroof(i) = sample_2d(croof, x2, x3)
                end do
                ! Slip along strike: a bell-shaped profile within the slip patch,
                ! which vanishes with zero slope at the lateral tips, so the
                ! sections do not tear apart there
                if (th(k)%radius > 0) then
                    taper = max(0.0, 1.0 - ((eta - th(k)%center)/th(k)%radius)**2)**2
                else
                    taper = 1.0
                end if
                slip = th(k)%slip*taper
                call thrust_section(th(k), is, xh, hdec, hroof, th(k)%xtop + th(k)%dev(is), th(k)%ztop + ne1, slip)
                deallocate (xh, hdec, hroof)
            end do
            !$omp end parallel do

            ! Visible separation along the fault in each section
            dtdz = mean(abs(t(2:n1, :, :) - t(1:n1 - 1, :, :))) + float_tiny
            !$omp parallel do private(is, eta, xa, za, xb, zb, nn, i, xsa, zsa, xsb, zsb)
            do is = 1, isec
                eta = th(k)%eta0 + (is - 1)*th(k)%deta
                if (th(k)%sl(is) <= 0) then
                    cycle
                end if
                call side_points(th(k), is, xa, za, xb, zb, nn)
                do i = 1, nn
                    call thrust_backward(th(k), is, xa(i), za(i), xsa, zsa)
                    call thrust_backward(th(k), is, xb(i), zb(i), xsb, zsb)
                    th(k)%sep(i, is) = abs( &
                        (sample_3d(t, zsa, c2 + xsa*cp - eta*sp, c3 + xsa*sp + eta*cp) &
                        - sample_3d(t, za(i), c2 + xa(i)*cp - eta*sp, c3 + xa(i)*sp + eta*cp)) &
                        - (sample_3d(t, zsb, c2 + xsb*cp - eta*sp, c3 + xsb*sp + eta*cp) &
                        - sample_3d(t, zb(i), c2 + xb(i)*cp - eta*sp, c3 + xb(i)*sp + eta*cp)))/dtdz
                end do
                call clean_visibility(th(k)%sep(1:nn, is), 6, 10)
            end do
            !$omp end parallel do

            ! Gather, with the faults traveling as fields as in 2D
            ww = w
            tt = t
            ff = f
            ffdip = fdip
            ffstrike = fstrike
            ffrake = frake
            ffdisp = fdisp
            pp = phi
            if (yn_facies) then
                mm = m
            end if
            call unused_padding(th, k, m1, m2)
            !$omp parallel do private(i, j, l, xi, eta, z, e, is, ws, xs, zs, xs2, zs2, x2, x3, &
            !$omp u, dperp, dip, u2, dperp2, dip2, sepv, sl, dl, beta, tdip, idv) collapse(3)
            do l = 1 + m2, n3 - m2
                do j = 1 + m2, n2 - m2
                    do i = 1 + m1, n1 - m1
                        xi = (j - 1.0 - c2)*cp + (l - 1.0 - c3)*sp
                        eta = -(j - 1.0 - c2)*sp + (l - 1.0 - c3)*cp
                        z = i - 1.0
                        e = (eta - th(k)%eta0)/th(k)%deta
                        is = min(max(int(e) + 1, 1), th(k)%nsec - 1)
                        ws = min(max(e - (is - 1), 0.0), 1.0)
                        call thrust_backward(th(k), is, xi, z, xs, zs)
                        call thrust_backward(th(k), is + 1, xi, z, xs2, zs2)
                        xs = (1.0 - ws)*xs + ws*xs2
                        zs = (1.0 - ws)*zs + ws*zs2
                        x2 = c2 + xs*cp - eta*sp
                        x3 = c3 + xs*sp + eta*cp
                        w(i, j, l) = sample_3d(ww, zs, x2, x3)
                        t(i, j, l) = sample_3d(tt, zs, x2, x3)
                        if (yn_facies) then
                            m(i, j, l) = nearest_3d(mm, zs, x2, x3)
                        end if
                        call thrust_point(th(k), is, xi, z, u, dperp, dip)
                        call thrust_point(th(k), is + 1, xi, z, u2, dperp2, dip2)
                        u = (1.0 - ws)*u + ws*u2
                        dperp = (1.0 - ws)*dperp + ws*dperp2
                        dip = (1.0 - ws)*dip + ws*dip2
                        sepv = (1.0 - ws)*lin(th(k)%sep(:, is), 0.0, th(k)%dsep, u) &
                            + ws*lin(th(k)%sep(:, is + 1), 0.0, th(k)%dsep, u)
                        sl = (1.0 - ws)*th(k)%sl(is) + ws*th(k)%sl(is + 1)
                        dl = (1.0 - ws)*th(k)%lf(is) + ws*th(k)%lf(is + 1)
                        if (dperp <= 0 .and. dperp > -dband .and. u <= dl .and. sl >= 0.5) then
                            call thrust_backward(th(k), is, xi - (0.5 - dperp)*sin(dip), z - (0.5 - dperp)*cos(dip), xs, zs)
                            call thrust_backward(th(k), is + 1, xi - (0.5 - dperp)*sin(dip), z - (0.5 - dperp)*cos(dip), xs2, zs2)
                            xs = (1.0 - ws)*xs + ws*xs2
                            zs = (1.0 - ws)*zs + ws*zs2
                            x2 = c2 + xs*cp - eta*sp
                            x3 = c3 + xs*sp + eta*cp
                            idv = nearest_3d(ff, zs, x2, x3)
                            phi(i, j, l) = max(sample_id_3d(pp, ff, idv, dband, zs, x2, x3) - (0.5 - dperp), -dband)
                        else
                            idv = nearest_3d(ff, zs, x2, x3)
                            phi(i, j, l) = sample_id_3d(pp, ff, idv, dband, zs, x2, x3)
                        end if
                        f(i, j, l) = idv
                        fdip(i, j, l) = nearest_3d(ffdip, zs, x2, x3)
                        fstrike(i, j, l) = nearest_3d(ffstrike, zs, x2, x3)
                        frake(i, j, l) = nearest_3d(ffrake, zs, x2, x3)
                        fdisp(i, j, l) = nearest_3d(ffdisp, zs, x2, x3)
                        ! This fault where its offset is visible
                        if (abs(dperp) < dband .and. u <= dl .and. sl >= 0.5) then
                            if ((abs(dperp) < 0.5*fwidth .or. abs(dperp) < abs(phi(i, j, l))) .and. sepv >= 0.5) then
                                f(i, j, l) = nf + k
                                ! The fault dips toward the hinterland; where the ramp
                                ! trace turns by beta, the dip direction turns by -beta
                                ! and the true dip is steeper than the dip in the section
                                beta = atan((1.0 - ws)*th(k)%ddev(is) + ws*th(k)%ddev(is + 1))
                                tdip = atan(tan(dip)/cos(beta))
                                call strike_dip(th(k)%azimuth - beta, tdip, fstrike(i, j, l), fdip(i, j, l))
                                frake(i, j, l) = 0.5*const_pi
                                fdisp(i, j, l) = sl
                                phi(i, j, l) = dperp
                            end if
                        end if
                    end do
                end do
            end do
            !$omp end parallel do

        end do

        ! The labels: within half the fault width of the faults
        call thrust_rasterize_3d(phi, f, fdip, fstrike, frake, fdisp, nf, fwidth)

    end subroutine thrust_apply_3d

    !
    !> Strike in [0, pi) and dip in (0, pi) in the convention of rgm3_curved
    !> (strike vector (sin(strike), cos(strike)) in (x2, x3); a dip below pi/2
    !> dips toward (cos(strike), -sin(strike))) for a fault that dips by dip
    !> toward the hinterland of a thrust with transport azimuth psi; a negative
    !> dip dips toward the foreland
    !
    pure subroutine strike_dip(psi, dip, strike, dipl)

        real, intent(in) :: psi, dip
        real, intent(out) :: strike, dipl

        if (dip >= 0) then
            strike = modulo(real(const_pi) - psi, 2.0*real(const_pi))
            dipl = dip
        else
            strike = modulo(-psi, 2.0*real(const_pi))
            dipl = -dip
        end if
        if (strike >= const_pi) then
            strike = strike - real(const_pi)
            dipl = real(const_pi) - dipl
        end if

    end subroutine strike_dip

end module geological_model_thrust
