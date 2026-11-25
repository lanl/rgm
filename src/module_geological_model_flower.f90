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
! Flower-like branched fault systems of rgm2_curved and rgm3_curved, called
! flowers here for short: closely aligned listric faults with their bottoms
! chopped out. A flower has a steep master fault and two sets of branches, one
! leaning toward each side. The branches of a set share one listric shape
! (steep at the top, flatter with depth, as the faults of rgm2_curved), so
! they are copies of one curve shifted sideways at a quasi regular spacing.
! Followed downward, each branch ends where it meets the master fault or the
! first branch of the other set (its bottom is chopped out there). The master
! fault bends into the first branch of one set at the split point. Each fault
! shifts its hanging wall along its dip, as the faults of rgm2_curved do; the
! displacement decays away from the fault, so the blocks tilt toward the
! faults, and fades to zero toward the start of a branch. Normal faults make
! a negative flower (a sag), reverse faults a positive flower (a pop-up). In a
! 3D model, a flower is this section across the strike, carried along the
! trace of the master fault in map view. The faults are carried as chains of
! points (meshes in 3D) that move with the rocks, and the labels are marked
! from them at the end, so that they keep their width.
!
module geological_model_flower

    use libflit
    use geological_model_utility
    use geological_model_label

    implicit none

    ! Degrees to radians and pi in single precision
    real, parameter :: d2r = real(const_deg2rad)
    real, parameter :: rpi = real(const_pi)

    !
    !> A fault of a flower: its trace x(i) (0-based column at the 0-based row
    !> i - 1 of the padded model), the row where it starts (its bottom, on its
    !> parent), its displacement (positive for a normal fault), the width of
    !> the decay of the displacement away from the fault (0: no decay) and the
    !> side of its hanging wall (+1: toward +x2); the master fault starts
    !> below the model
    !
    type flower_fault
        real :: zroot = 0.0, disp = 0.0, decay = 0.0, hside = 0.0
        logical :: master = .false.
        real, allocatable, dimension(:) :: x
    end type flower_fault

    !
    !> A flower: its faults and the order in which they act
    !
    type flower
        type(flower_fault), allocatable, dimension(:) :: fault
        integer, allocatable, dimension(:) :: order
        ! Column and row of the split point
        real :: xs = 0.0, zs = 0.0
    end type flower

    !
    !> A flower of a 3D model: the flower sec of the section across the strike
    !> (the traces of its faults count from the column sec%xs of the split
    !> point), carried along the trace of the master fault in map view. The
    !> trace passes through the point (c2 along x2, c3 along x3) and strikes at
    !> the angle strike (from x3 toward x2); along it (from u0 on, every grid
    !> point), g is the deviation of the trace across the strike, s a scale of
    !> the section, and a(:, k) the weight of the displacement of the fault k.
    !
    type flower_3d
        type(flower) :: sec
        real :: c2 = 0.0, c3 = 0.0, strike = 0.0, u0 = 0.0
        real, allocatable, dimension(:) :: g, s
        real, allocatable, dimension(:, :) :: a
    end type flower_3d

    !
    !> A fault of a flower of a 3D model, in the map coordinates (u, v) around
    !> the point (c2 along x2, c3 along x3): u runs along the strike (an angle
    !> from x3 toward x2, with its cosine cphi and sine sphi) and v across it.
    !> The fault surface is v = g(u) + s(u)*x(row): x is the trace of the fault
    !> in the section across the strike, g the trace of the master fault in map
    !> view, and s the scale of the section; a(u) is the weight of the
    !> displacement along u. g, s and a are sampled every grid point from u0
    !> on. The hanging wall is on the side hside of the surface in v, where the
    !> weights alpha(row) and a(u) are positive. Its rocks shift along the
    !> local dip by disp*alpha*a (positive: down the dip), times a Gaussian
    !> decay of the width decay away from the surface (0: none).
    !
    type flower_surface
        real :: c2 = 0.0, c3 = 0.0, cphi = 1.0, sphi = 0.0, u0 = 0.0
        real, allocatable, dimension(:) :: g, s, a, x, alpha
        real :: disp = 0.0, decay = 0.0, hside = 1.0
    end type flower_surface

    private
    public :: flower, flower_fault, flower_3d
    public :: flower_draw_2d
    public :: flower_apply_2d
    public :: flower_draw_3d
    public :: flower_apply_3d

contains

    !
    !> Lean of a listric branch from the vertical at the row z: phi1 at the
    !> model top (row ztm) and above it, growing by dphi over the model
    !> depth m1
    !
    pure function branch_lean(z, phi1, dphi, ztm, m1) result(phi)

        real, intent(in) :: z, phi1, dphi, ztm, m1
        real :: phi

        phi = min(phi1 + dphi*min(max((z - ztm)/m1, 0.0), 1.25), 65.0*d2r)

    end function branch_lean

    !
    !> Trace of a branch from the point (xr, zr) upward, leaning toward the
    !> side (+1: toward +x2) by branch_lean; below zr, the trace keeps its
    !> values
    !
    subroutine branch_trace(x, xr, zr, side, phi1, dphi, ztm, m1)

        real, dimension(:), intent(inout) :: x
        real, intent(in) :: xr, zr, side, phi1, dphi, ztm, m1

        integer :: i, i1, n

        n = size(x)
        i1 = min(max(floor(zr), 0), n - 1) + 1
        x(i1) = xr + side*(zr - (i1 - 1.0))*tan(branch_lean(zr, phi1, dphi, ztm, m1))
        do i = i1 - 1, 1, -1
            x(i) = x(i + 1) + side*tan(branch_lean(i - 0.5, phi1, dphi, ztm, m1))
        end do

    end subroutine branch_trace

    !
    !> Value of a trace x (at the rows 0, 1, ...) at the fractional row z
    !
    pure function row_value(x, z) result(v)

        real, dimension(:), intent(in) :: x
        real, intent(in) :: z
        real :: v

        integer :: i
        real :: w

        i = min(max(int(z), 0), size(x) - 2)
        w = min(max(z - i, 0.0), 1.0)
        v = (1.0 - w)*x(i + 1) + w*x(i + 2)

    end function row_value

    !
    !> Draw nfl flowers in a padded 2D model of n1 rows, with the model top at
    !> the row ne1 and the model columns from ne2 on; m1 and m2 are the model
    !> sizes. The master fault (dip from fdip) rises to the split point (depth
    !> from fz, relative to m1) and bends into the first branch of one set;
    !> the first branch of the other set starts at the split point. A branch
    !> dips by an angle from bdip at the model top, and its dip decreases by
    !> an angle from ddip over the model depth (a listric fault). Each set has
    !> a number of branches from nbr (drawn for each flower), at the horizontal
    !> spacing from spacing (relative to m1, with a jitter of 25%): the
    !> branches outside the first one start on the master fault below the
    !> split point, and those inside it cross above the split point and start
    !> on the first branch of the other set. disp: range of the displacement
    !> ([0, 0] = [0.02, 0.05]*m1); mdisp: range of the displacement of the
    !> master fault relative to it; decay: range of the width of the decay of
    !> the displacement away from each fault (relative to m1). ftype:
    !> 'negative' (normal faults), 'positive' (reverse faults), 'hybrid'
    !> (reverse master fault and first branches, the other branches normal
    !> and smaller) or 'random' (negative or positive).
    !
    function flower_draw_2d(nfl, n1, ne1, ne2, m1, m2, seed, fdip, fz, nbr, bdip, ddip, disp, mdisp, spacing, decay, ftype) result(fl)

        integer, intent(in) :: nfl, n1, ne1, ne2, m1, m2, seed
        integer, dimension(1:2), intent(in) :: nbr
        real, dimension(1:2), intent(in) :: fdip, fz, bdip, ddip, disp, mdisp, spacing, decay
        character(len=*), intent(in) :: ftype
        type(flower), allocatable, dimension(:) :: fl

        type(flower_fault), allocatable, dimension(:) :: ff
        real, allocatable, dimension(:) :: pos
        real, dimension(1:2) :: rdisp, phi1, dphi
        real :: r(1:16), rr(1:6), xs, zs, ztm, sgn, side, sd, mlean, kr, z, q, lam, dec, dx, zn, zi, zo, p1, dp, rm1
        integer :: f, i, j, s, so, nf, ia, nin, is, nins(1:2), nouts(1:2), nb
        integer, allocatable, dimension(:) :: kind
        logical :: yn_done

        allocate (fl(1:nfl))
        if (nfl == 0) then
            return
        end if
        if (maxval(abs(disp)) == 0) then
            rdisp = [0.02, 0.05]*m1
        else
            rdisp = disp
        end if
        ztm = real(ne1)
        rm1 = real(m1)
        ! One flower lies near the center of the model; several flowers spread
        ! evenly across it, with a jitter
        pos = random(nfl, range=[-0.15, 0.15]*m2/nfl, seed=safe_seed(int(seed, 8)*331))
        if (nfl == 1) then
            pos = pos*2.0/3.0
        else
            do f = 1, nfl
                pos(f) = pos(f) + ((f - 0.5)/nfl - 0.5)*m2
            end do
        end if

        do f = 1, nfl

            r = random(16, range=[0.0, 1.0], seed=safe_seed(int(seed, 8)*337 + f))
            nb = min(nbr(1) + int(r(15)*(nbr(2) - nbr(1) + 1)), nbr(2))
            select case (ftype)
                case ('negative')
                    sgn = 1.0
                case ('positive', 'hybrid')
                    sgn = -1.0
                case default
                    sgn = merge(1.0, -1.0, r(1) < 0.5)
            end select
            zs = ne1 + (fz(1) + (fz(2) - fz(1))*r(2))*m1
            xs = ne2 + 0.5*(m2 - 1.0) + pos(f)
            ! The set whose first branch continues the master fault (1: the
            ! branches lean toward -x2, 2: toward +x2)
            is = merge(2, 1, r(3) < 0.5)
            side = merge(1.0, -1.0, is == 2)
            mlean = (90.0 - (fdip(1) + (fdip(2) - fdip(1))*r(4)))*d2r
            ! The listric shape of each set
            do s = 1, 2
                phi1(s) = (90.0 - (bdip(1) + (bdip(2) - bdip(1))*r(4 + s)))*d2r
                dphi(s) = (ddip(1) + (ddip(2) - ddip(1))*r(6 + s))*d2r
            end do
            dec = (decay(1) + (decay(2) - decay(1))*r(9))*m1

            allocate (ff(1:2*max(nb, 1)), kind(1:2*max(nb, 1)))

            ! The master fault: steep below the split point, where it leans
            ! slightly toward the side of its branch (so that its hanging wall
            ! stays on one side), and bending into that branch over 2*kr rows
            kr = 0.03*m1
            ff(1)%master = .true.
            ff(1)%zroot = n1 + 1.0
            ff(1)%disp = sgn*(rdisp(1) + (rdisp(2) - rdisp(1))*r(10))*(mdisp(1) + (mdisp(2) - mdisp(1))*r(16))
            ff(1)%decay = max(dec, 0.4*m1)
            ff(1)%hside = -side
            kind(1) = 0
            allocate (ff(1)%x(1:n1))
            do i = n1, 1, -1
                z = i - 1.0
                if (z >= zs + kr .or. i == n1) then
                    ff(1)%x(i) = xs - side*(z - zs)*tan(mlean)
                else
                    q = min(max((zs + kr - (z + 0.5))/(2.0*kr), 0.0), 1.0)
                    q = q*q*(3.0 - 2.0*q)
                    lam = mlean + (branch_lean(z + 0.5, phi1(is), dphi(is), ztm, rm1) - mlean)*q
                    ff(1)%x(i) = ff(1)%x(i + 1) + side*tan(lam)
                end if
            end do

            ! The first branch of the other set starts at the split point
            so = 3 - is
            ff(2)%x = ff(1)%x
            call branch_trace(ff(2)%x, row_value(ff(1)%x, zs), zs, -side, phi1(so), dphi(so), ztm, rm1)
            ff(2)%zroot = zs
            ff(2)%disp = sgn*(rdisp(1) + (rdisp(2) - rdisp(1))*r(11))
            ff(2)%decay = dec
            ff(2)%hside = side
            kind(2) = 1
            nf = 2

            ! The other branches of each set: copies of the first branch
            ! shifted sideways. An outer branch starts on the master fault below
            ! the split point; an inner branch crosses above the split point
            ! and starts on the first branch of the other set (fault ia)
            ! Numbers of outer branches (up to two in each set) and of inner
            ! branches (the rest in one set; in the other set up to two, and
            ! none where the first set has more than two, so that the inner
            ! branches of the two sets do not cross in a dense lattice)
            s = merge(1, 2, r(12) < 0.5)
            nouts(s) = min(int(r(13)*3.0), 2, nb - 1)
            nins(s) = nb - 1 - nouts(s)
            nouts(3 - s) = min(int(r(14)*3.0), 2, nb - 1)
            nins(3 - s) = merge(0, min(nb - 1 - nouts(3 - s), 2), nins(s) > 2)
            do s = 1, 2
                sd = merge(1.0, -1.0, s == 2)
                so = 3 - s
                ia = merge(2, 1, s == is)
                nin = nins(s)
                zi = zs
                zo = zs
                do j = 1, nin + nouts(s)
                    rr = random(6, range=[0.0, 1.0], seed=safe_seed(int(seed, 8)*349 + 1000*f + 100*s + j))
                    dx = (spacing(1) + (spacing(2) - spacing(1))*rr(1))*m1*(0.75 + 0.5*rr(2))
                    ! The shape of this branch, close to that of its set
                    p1 = max(phi1(s) + (2.0*rr(3) - 1.0)*d2r, 2.0*d2r)
                    dp = dphi(s) + 2.0*(2.0*rr(4) - 1.0)*d2r
                    yn_done = .false.
                    if (j <= nin) then
                        zn = zi - dx/(tan(branch_lean(zi, phi1(s), dphi(s), ztm, rm1)) &
                            + tan(branch_lean(zi, phi1(so), dphi(so), ztm, rm1)))
                        if (zn >= ztm + 0.3*m1) then
                            zi = zn
                            nf = nf + 1
                            ff(nf)%x = ff(ia)%x
                            call branch_trace(ff(nf)%x, row_value(ff(ia)%x, zn), zn, sd, p1, dp, ztm, rm1)
                            yn_done = .true.
                        end if
                    else
                        zn = zo + dx/tan(branch_lean(zo, phi1(s), dphi(s), ztm, rm1))
                        if (zn <= ztm + 0.95*m1) then
                            zo = zn
                            nf = nf + 1
                            ff(nf)%x = ff(1)%x
                            call branch_trace(ff(nf)%x, row_value(ff(1)%x, zn), zn, sd, p1, dp, ztm, rm1)
                            yn_done = .true.
                        end if
                    end if
                    if (yn_done) then
                        ff(nf)%zroot = zn
                        ff(nf)%disp = sgn*(rdisp(1) + (rdisp(2) - rdisp(1))*rr(5))
                        if (ftype == 'hybrid') then
                            ff(nf)%disp = 0.7*abs(ff(nf)%disp)
                        end if
                        ff(nf)%decay = dec*(0.8 + 0.4*rr(6))
                        ff(nf)%hside = -sd
                        kind(nf) = merge(merge(3, 4, s == is), 2, j <= nin)
                    end if
                end do
            end do

            ! The order of action. A fault must not move the start of an older
            ! branch off its parent, and a branch must start on its parent where
            ! the parent is: the outer branches act first, from the deepest start
            ! up (a fault does not move the rocks below its start); then the
            ! inner branches that start on the first branch of the other set,
            ! from the lowest start up (each lies in the footwall of the next),
            ! and that first branch; then the inner branches that start on the
            ! master fault, and the master fault, whose hanging wall carries
            ! the starts of its branches along it
            allocate (fl(f)%order(1:nf))
            j = 0
            do
                i = 0
                do s = 3, nf
                    if (kind(s) == 2) then
                        if (i == 0) then
                            i = s
                        else if (ff(s)%zroot > ff(i)%zroot) then
                            i = s
                        end if
                    end if
                end do
                if (i == 0) then
                    exit
                end if
                j = j + 1
                fl(f)%order(j) = i
                kind(i) = -1
            end do
            do s = 3, nf
                if (kind(s) == 3) then
                    j = j + 1
                    fl(f)%order(j) = s
                end if
            end do
            j = j + 1
            fl(f)%order(j) = 2
            do s = 3, nf
                if (kind(s) == 4) then
                    j = j + 1
                    fl(f)%order(j) = s
                end if
            end do
            fl(f)%order(nf) = 1
            fl(f)%fault = ff(1:nf)
            fl(f)%xs = xs
            fl(f)%zs = zs
            deallocate (ff, kind)

        end do

    end function flower_draw_2d

    !
    !> Apply the flowers to the padded 2D model arrays (velocity w, RGT t,
    !> facies m), each fault as a fault of rgm2_curved: the hanging wall shifts
    !> along the local dip by the displacement. The displacement decays away
    !> from the fault with a Gaussian profile (as with yn_disp_decay of
    !> rgm2_curved), and fades to zero toward the start of a branch over the
    !> length lfade, or over three times the displacement if that is longer.
    !>
    !> The faults are carried as chains of points in pc (labels id0 + 1,
    !> id0 + 2, ...): each fault moves the points of the older faults (those
    !> of the thrusts included) as it moves the rocks, and then adds its own
    !> points, from its start upward. The labels are marked from the chains
    !> afterwards (chain_label_2d), so that they keep their width however the
    !> faults are bent or cut.
    !
    subroutine flower_apply_2d(fl, id0, lfade, w, t, m, pc, yn_facies)

        type(flower), dimension(:), intent(in) :: fl
        integer, intent(in) :: id0
        real, intent(in) :: lfade
        real, dimension(:, :), intent(inout) :: w, t
        real, allocatable, dimension(:, :), intent(inout) :: m
        type(fault_chain), intent(inout) :: pc
        logical, intent(in) :: yn_facies

        real, allocatable, dimension(:, :) :: ww, tt, mm
        real, allocatable, dimension(:) :: b, th, alpha, qx, qz, qd
        integer, allocatable, dimension(:) :: qid, qch
        real :: x0, dloc, zi, xj, q, dsp, dec, lf, hs, zs, xs, zp, xp, zn, xn, z1
        integer :: n1, n2, fi, k, io, i, j, l, id, ni, nj, it, nq
        logical :: yn_master
        ! Spacing of the points along a fault (rows)
        real, parameter :: dp = 0.5

        n1 = size(w, 1)
        n2 = size(w, 2)
        id = id0
        allocate (b(1:n1), th(1:n1), alpha(1:n1))
        call chain_append(pc, [real ::], [real ::], [real ::], [integer ::], [integer ::])

        do fi = 1, size(fl)
            do io = 1, size(fl(fi)%fault)

                k = fl(fi)%order(io)
                b = fl(fi)%fault(k)%x
                dsp = fl(fi)%fault(k)%disp
                dec = fl(fi)%fault(k)%decay
                hs = fl(fi)%fault(k)%hside
                yn_master = fl(fi)%fault(k)%master
                lf = max(lfade, 3.0*abs(dsp))
                if (dsp < 0) then
                    lf = max(2.5*lfade, 6.0*abs(dsp))
                end if

                ! Local dip (below pi/2 the fault dips toward -x2, as in
                ! rgm2_curved) and the weight of the displacement: zero below the
                ! start of a branch, rising smoothly over lf
                do i = 1, n1
                    th(i) = 0.5*const_pi + atan(0.5*(b(min(i + 1, n1)) - b(max(i - 1, 1))))
                    if (yn_master) then
                        alpha(i) = 1.0
                    else
                        q = min(max((fl(fi)%fault(k)%zroot - (i - 1.0))/lf, 0.0), 1.0)
                        alpha(i) = q*q*(3.0 - 2.0*q)
                    end if
                end do

                ! Shift the hanging wall: each point takes the values of its source
                ! point, up the dip from it for a normal fault
                ww = w
                tt = t
                if (yn_facies) then
                    mm = m
                end if
                !$omp parallel do private(i, j, x0, dloc, zi, xj, ni, nj) collapse(2)
                do j = 1, n2
                    do i = 1, n1
                        x0 = j - 1.0
                        if (alpha(i) <= 0 .or. (x0 - b(i))*hs < 0) then
                            cycle
                        end if
                        dloc = dsp*alpha(i)
                        if (dec > 0) then
                            dloc = dloc*exp(-0.5*((x0 - b(i))/dec)**2)
                        end if
                        zi = i - dloc*sin(th(i))
                        xj = j + dloc*cos(th(i))
                        w(i, j) = bilinear(ww, zi - 1.0, xj - 1.0)
                        t(i, j) = bilinear(tt, zi - 1.0, xj - 1.0)
                        if (yn_facies) then
                            ni = min(max(nint(zi), 1), n1)
                            nj = min(max(nint(xj), 1), n2)
                            m(i, j) = mm(ni, nj)
                        end if
                    end do
                end do
                !$omp end parallel do

                ! Move the points of the older faults with the rocks: the point
                ! that the shift took from the source s is now at p = s + u(p),
                ! where u is the displacement at p (a fixed-point iteration)
                !$omp parallel do private(l, zs, xs, zp, xp, zn, xn, it, dloc)
                do l = 1, pc%n
                    zs = pc%z(l)
                    xs = pc%x(l)
                    if ((xs - row_value(b, zs))*hs < 0) then
                        cycle
                    end if
                    zp = zs
                    xp = xs
                    do it = 1, 40
                        dloc = dsp*row_value(alpha, zp)
                        if (dec > 0) then
                            dloc = dloc*exp(-0.5*((xp - row_value(b, zp))/dec)**2)
                        end if
                        zn = zs + dloc*sin(row_value(th, zp))
                        xn = xs - dloc*cos(row_value(th, zp))
                        if (abs(zn - zp) + abs(xn - xp) < 1.0e-3) then
                            zp = zn
                            xp = xn
                            exit
                        end if
                        zp = zn
                        xp = xn
                    end do
                    pc%z(l) = zp
                    pc%x(l) = xp
                end do
                !$omp end parallel do

                ! The points of this fault, from just above its start upward
                z1 = merge(n1 - 1.0, fl(fi)%fault(k)%zroot - 0.5*dp, yn_master)
                nq = max(floor(z1/dp), 0) + 1
                allocate (qx(1:nq), qz(1:nq), qd(1:nq), qid(1:nq), qch(1:nq))
                do i = 1, nq
                    qz(i) = z1 - (i - 1)*dp
                    qx(i) = row_value(b, qz(i))
                    qd(i) = dsp*row_value(alpha, qz(i))
                end do
                qid = id + k
                qch = 1
                call chain_append(pc, qx, qz, qd, qid, qch)
                deallocate (qx, qz, qd, qid, qch)

            end do
            id = id + size(fl(fi)%fault)
        end do

    contains

        pure function bilinear(a, z, x) result(val)

            real, dimension(:, :), intent(in) :: a
            real, intent(in) :: z, x
            real :: val

            integer :: i1, j1, na, nb2
            real :: wz, wx, zc, xc

            na = size(a, 1)
            nb2 = size(a, 2)
            zc = min(max(z, 0.0), na - 1.001)
            xc = min(max(x, 0.0), nb2 - 1.001)
            i1 = int(zc) + 1
            j1 = int(xc) + 1
            wz = zc - (i1 - 1)
            wx = xc - (j1 - 1)
            val = a(i1, j1)*(1 - wz)*(1 - wx) + a(i1 + 1, j1)*wz*(1 - wx) &
                + a(i1, j1 + 1)*(1 - wz)*wx + a(i1 + 1, j1 + 1)*wz*wx

        end function bilinear

    end subroutine flower_apply_2d

    !
    !> Draw nfl flowers in a padded 3D model (n2 by n3 columns; the model of
    !> m1 by m2 by m3 starts at the row ne1 and the columns ne2 and ne3). Each
    !> flower is a flower of flower_draw_2d in the section across its strike,
    !> carried along the trace of its master fault. fstrike: range of the
    !> strike of the master fault (degrees, from x3 toward x2); dstrike: range
    !> of the largest deviation of the local strike along the trace (degrees);
    !> radius: range of the half length along strike over which the
    !> displacement of a branch falls to zero (relative to (m2 + m3)/2), so a
    !> branch can die out along strike; the master fault and the two first
    !> branches run through. The section also widens and narrows along strike
    !> by up to 25%. The flowers spread evenly across the mean strike.
    !
    function flower_draw_3d(nfl, n1, n2, n3, ne1, ne2, ne3, m1, m2, m3, seed, fdip, fz, nbr, bdip, ddip, disp, mdisp, &
            spacing, decay, ftype, fstrike, dstrike, radius) result(fl3)

        integer, intent(in) :: nfl, n1, n2, n3, ne1, ne2, ne3, m1, m2, m3, seed
        integer, dimension(1:2), intent(in) :: nbr
        real, dimension(1:2), intent(in) :: fdip, fz, bdip, ddip, disp, mdisp, spacing, decay, fstrike, dstrike, radius
        character(len=*), intent(in) :: ftype
        type(flower_3d), allocatable, dimension(:) :: fl3

        type(flower), allocatable, dimension(:) :: fl
        type(fractal_noise_1d) :: sn
        real, allocatable, dimension(:) :: dphi
        real :: r(1:4), rr(1:2), mw, off, cphi, sphi, umin, umax, tp, uc, ul, u
        integer :: f, j, k, nu, nk, iu

        allocate (fl3(1:nfl))
        if (nfl == 0) then
            return
        end if
        mw = 0.5*(m2 + m3)
        fl = flower_draw_2d(nfl, n1, ne1, 0, m1, nint(mw), seed, fdip, fz, nbr, bdip, ddip, disp, mdisp, spacing, decay, ftype)

        do f = 1, nfl

            r = random(4, range=[0.0, 1.0], seed=safe_seed(int(seed, 8)*353 + f))
            fl3(f)%sec = fl(f)
            fl3(f)%strike = (fstrike(1) + (fstrike(2) - fstrike(1))*r(1))*d2r
            cphi = cos(fl3(f)%strike)
            sphi = sin(fl3(f)%strike)
            ! The split point lies off the model center across the strike
            off = fl(f)%xs - 0.5*(nint(mw) - 1.0)
            fl3(f)%c3 = ne3 + 0.5*(m3 - 1.0) - off*sphi
            fl3(f)%c2 = ne2 + 0.5*(m2 - 1.0) + off*cphi

            ! The along-strike axis covers the padded model
            umin = float_huge
            umax = -float_huge
            do k = 1, 2
                do j = 1, 2
                    tp = ((k - 1)*(n3 - 1.0) - fl3(f)%c3)*cphi + ((j - 1)*(n2 - 1.0) - fl3(f)%c2)*sphi
                    umin = min(umin, tp)
                    umax = max(umax, tp)
                end do
            end do
            nu = ceiling(umax - umin) + 2
            fl3(f)%u0 = umin

            ! The trace of the master fault: a smooth deviation of the strike
            sn%n1 = nu
            sn%octaves = 2
            sn%periods1 = 1
            sn%seed = safe_seed(int(seed, 8)*359 + f)
            dphi = sn%generate()
            dphi = dphi/maxval(abs(dphi))*(dstrike(1) + (dstrike(2) - dstrike(1))*r(2))*d2r
            fl3(f)%g = cumsum(tan(dphi))
            fl3(f)%g = fl3(f)%g - fl3(f)%g(min(max(nint(-umin) + 1, 1), nu))

            ! The section widens and narrows along strike
            sn%seed = safe_seed(int(seed, 8)*367 + f)
            dphi = sn%generate()
            fl3(f)%s = 1.0 + 0.25*dphi/maxval(abs(dphi))

            ! The displacement of each branch along strike: an elliptical patch
            nk = size(fl(f)%fault)
            allocate (fl3(f)%a(1:nu, 1:nk))
            fl3(f)%a = 1.0
            do k = 3, nk
                rr = random(2, range=[0.0, 1.0], seed=safe_seed(int(seed, 8)*373 + 100*f + k))
                uc = (0.6*rr(1) - 0.3)*mw
                ul = (radius(1) + (radius(2) - radius(1))*rr(2))*mw
                do iu = 1, nu
                    u = umin + (iu - 1.0)
                    fl3(f)%a(iu, k) = max(0.0, 1.0 - ((u - uc)/ul)**2)
                end do
            end do

        end do

    end function flower_draw_3d

    !
    !> Apply the flowers to the padded 3D model arrays (velocity w, RGT t,
    !> facies m), as flower_apply_2d: each fault shifts its hanging wall along
    !> the local dip, with the displacement decaying away from the fault and
    !> fading toward the start of a branch. The faults are carried as meshes
    !> of points in pm (labels id0 + 1, id0 + 2, ...): each fault moves the
    !> meshes of the older faults (those of the thrusts included) as it moves
    !> the rocks, and then adds its own mesh. The labels are marked from the
    !> meshes afterwards (mesh_label_3d).
    !
    subroutine flower_apply_3d(fl3, id0, lfade, w, t, m, pm, yn_facies)

        type(flower_3d), dimension(:), intent(in) :: fl3
        integer, intent(in) :: id0
        real, intent(in) :: lfade
        real, dimension(:, :, :), intent(inout) :: w, t
        real, allocatable, dimension(:, :, :), intent(inout) :: m
        type(fault_mesh), intent(inout) :: pm
        logical, intent(in) :: yn_facies

        type(flower_surface) :: sf
        integer :: n1, fi, io, k, i, id
        real :: lf, q

        n1 = size(w, 1)
        id = id0
        allocate (sf%alpha(1:n1))

        do fi = 1, size(fl3)
            sf%c2 = fl3(fi)%c2
            sf%c3 = fl3(fi)%c3
            sf%cphi = cos(fl3(fi)%strike)
            sf%sphi = sin(fl3(fi)%strike)
            sf%u0 = fl3(fi)%u0
            sf%g = fl3(fi)%g
            sf%s = fl3(fi)%s
            do io = 1, size(fl3(fi)%sec%fault)

                k = fl3(fi)%sec%order(io)
                sf%a = fl3(fi)%a(:, k)
                sf%x = fl3(fi)%sec%fault(k)%x - fl3(fi)%sec%xs
                sf%disp = fl3(fi)%sec%fault(k)%disp
                sf%decay = fl3(fi)%sec%fault(k)%decay
                sf%hside = fl3(fi)%sec%fault(k)%hside

                ! The weight of the displacement with depth: zero below the start
                ! of a branch, rising smoothly (over a longer distance for a
                ! reverse fault, which stretches the layers there)
                lf = max(lfade, 3.0*abs(sf%disp))
                if (sf%disp < 0) then
                    lf = max(2.5*lfade, 6.0*abs(sf%disp))
                end if
                do i = 1, n1
                    if (fl3(fi)%sec%fault(k)%master) then
                        sf%alpha(i) = 1.0
                    else
                        q = min(max((fl3(fi)%sec%fault(k)%zroot - (i - 1.0))/lf, 0.0), 1.0)
                        sf%alpha(i) = q*q*(3.0 - 2.0*q)
                    end if
                end do

                call surface_shift_3d(sf, w, t, m, pm, yn_facies)
                call surface_mesh_3d(sf, pm, id + k, size(w, 2), size(w, 3))

            end do
            id = id + size(fl3(fi)%sec%fault)
        end do

    end subroutine flower_apply_3d

    !
    !> Whether the point (z, y, x) (0-based depth and columns along x2 and x3)
    !> lies in the hanging wall of the fault sf, and the displacement
    !> (dz, dy, dx) of the rocks at it: along the local dip, that is, down the
    !> surface and normal to its trace in map view
    !
    pure subroutine surface_kinematics(sf, z, y, x, yn, dz, dy, dx)

        type(flower_surface), intent(in) :: sf
        real, intent(in) :: z, y, x
        logical, intent(out) :: yn
        real, intent(out) :: dz, dy, dx

        integer :: i, iu, n1, nu
        real :: w, wu, xr, az, slope, d3, d2, u, v, uq, g, s, a, gp, sp, b, q, nn, th, d, dh

        n1 = size(sf%x)
        nu = size(sf%g)
        d3 = x - sf%c3
        d2 = y - sf%c2
        u = d3*sf%cphi + d2*sf%sphi
        v = -d3*sf%sphi + d2*sf%cphi
        uq = u - sf%u0
        iu = min(max(floor(uq) + 1, 1), nu - 1)
        wu = min(max(uq - (iu - 1.0), 0.0), 1.0)
        g = (1.0 - wu)*sf%g(iu) + wu*sf%g(iu + 1)
        s = (1.0 - wu)*sf%s(iu) + wu*sf%s(iu + 1)
        a = (1.0 - wu)*sf%a(iu) + wu*sf%a(iu + 1)
        gp = sf%g(iu + 1) - sf%g(iu)
        sp = sf%s(iu + 1) - sf%s(iu)
        i = min(max(floor(z) + 1, 1), n1 - 1)
        w = min(max(z + 1.0 - i, 0.0), 1.0)
        xr = (1.0 - w)*sf%x(i) + w*sf%x(i + 1)
        az = (1.0 - w)*sf%alpha(i) + w*sf%alpha(i + 1)
        slope = (1.0 - w)*0.5*(sf%x(i + 1) - sf%x(max(i - 1, 1))) + w*0.5*(sf%x(min(i + 2, n1)) - sf%x(i))
        b = g + s*xr

        ! The normal of the trace in map view (toward +v) and the dip
        q = gp + sp*xr
        nn = sqrt(1.0 + q**2)
        th = 0.5*rpi + atan(s*slope/nn)

        yn = az*a > 0 .and. (v - b)*sf%hside >= 0
        d = sf%disp*az*a
        if (sf%decay > 0) then
            d = d*exp(-0.5*((v - b)/(nn*sf%decay))**2)
        end if
        dz = d*sin(th)
        dh = -d*cos(th)
        dx = dh*(-sf%sphi - q*sf%cphi)/nn
        dy = dh*(sf%cphi - q*sf%sphi)/nn

    end subroutine surface_kinematics

    !
    !> Shift the hanging wall of the fault sf in the padded 3D model arrays
    !> (velocity w, RGT t, facies m): each grid point of it takes the values
    !> of its source point. Then move the meshes pm with the rocks: a point s
    !> stays where it is if it lies outside the hanging wall, and it appears
    !> at the point p of the hanging wall with p = s + u(p), where u is the
    !> displacement at p (a fixed-point iteration), if that p lies in the
    !> hanging wall.
    !
    subroutine surface_shift_3d(sf, w, t, m, pm, yn_facies)

        type(flower_surface), intent(in) :: sf
        real, dimension(:, :, :), intent(inout) :: w, t
        real, allocatable, dimension(:, :, :), intent(inout) :: m
        type(fault_mesh), intent(inout) :: pm
        logical, intent(in) :: yn_facies

        real, allocatable, dimension(:, :, :) :: ww, tt, mm
        real, allocatable, dimension(:) :: t2, t3, tz
        logical, allocatable, dimension(:) :: stay, move
        integer :: n1, n2, n3, i, j, k, l, it, np
        real :: dz, dy, dx, z0, y0, x0, z, y, x, zn, yn_, xn
        logical :: yn

        n1 = size(w, 1)
        n2 = size(w, 2)
        n3 = size(w, 3)
        ww = w
        tt = t
        if (yn_facies) then
            mm = m
        end if
        !$omp parallel do private(i, j, k, yn, dz, dy, dx) collapse(3)
        do k = 1, n3
            do j = 1, n2
                do i = 1, n1
                    call surface_kinematics(sf, i - 1.0, j - 1.0, k - 1.0, yn, dz, dy, dx)
                    if (.not. yn) then
                        cycle
                    end if
                    if (abs(dz) + abs(dy) + abs(dx) < 1.0e-4) then
                        cycle
                    end if
                    w(i, j, k) = trilinear(ww, i - 1.0 - dz, j - 1.0 - dy, k - 1.0 - dx)
                    t(i, j, k) = trilinear(tt, i - 1.0 - dz, j - 1.0 - dy, k - 1.0 - dx)
                    if (yn_facies) then
                        m(i, j, k) = mm(min(max(nint(i - dz), 1), n1), min(max(nint(j - dy), 1), n2), &
                            min(max(nint(k - dx), 1), n3))
                    end if
                end do
            end do
        end do
        !$omp end parallel do

        np = pm%n
        if (np == 0) then
            return
        end if
        allocate (t2(1:np), t3(1:np), tz(1:np), stay(1:np), move(1:np))
        !$omp parallel do private(l, it, z0, y0, x0, z, y, x, zn, yn_, xn, yn, dz, dy, dx)
        do l = 1, np
            t2(l) = pm%x2(l)
            t3(l) = pm%x3(l)
            tz(l) = pm%z(l)
            stay(l) = .false.
            move(l) = .false.
            if (.not. pm%v(l)) then
                cycle
            end if
            z0 = pm%z(l)
            y0 = pm%x2(l)
            x0 = pm%x3(l)
            call surface_kinematics(sf, z0, y0, x0, yn, dz, dy, dx)
            stay(l) = .not. yn
            z = z0
            y = y0
            x = x0
            do it = 1, 50
                call surface_kinematics(sf, z, y, x, yn, dz, dy, dx)
                zn = z0 + dz
                yn_ = y0 + dy
                xn = x0 + dx
                if (abs(zn - z) + abs(yn_ - y) + abs(xn - x) < 1.0e-3) then
                    z = zn
                    y = yn_
                    x = xn
                    exit
                end if
                z = zn
                y = yn_
                x = xn
            end do
            call surface_kinematics(sf, z, y, x, yn, dz, dy, dx)
            move(l) = yn
            tz(l) = z
            t2(l) = y
            t3(l) = x
        end do
        !$omp end parallel do
        call mesh_split_3d(pm, t2, t3, tz, stay, move)

    contains

        pure function trilinear(a, z, y, x) result(val)

            real, dimension(:, :, :), intent(in) :: a
            real, intent(in) :: z, y, x
            real :: val

            integer :: i1, j1, k1
            real :: wz, wy, wx, zc, yc, xc

            zc = min(max(z, 0.0), size(a, 1) - 1.001)
            yc = min(max(y, 0.0), size(a, 2) - 1.001)
            xc = min(max(x, 0.0), size(a, 3) - 1.001)
            i1 = int(zc) + 1
            j1 = int(yc) + 1
            k1 = int(xc) + 1
            wz = zc - (i1 - 1)
            wy = yc - (j1 - 1)
            wx = xc - (k1 - 1)
            val = a(i1, j1, k1)*(1 - wz)*(1 - wy)*(1 - wx) + a(i1 + 1, j1, k1)*wz*(1 - wy)*(1 - wx) &
                + a(i1, j1 + 1, k1)*(1 - wz)*wy*(1 - wx) + a(i1, j1, k1 + 1)*(1 - wz)*(1 - wy)*wx &
                + a(i1 + 1, j1 + 1, k1)*wz*wy*(1 - wx) + a(i1 + 1, j1, k1 + 1)*wz*(1 - wy)*wx &
                + a(i1, j1 + 1, k1 + 1)*(1 - wz)*wy*wx + a(i1 + 1, j1 + 1, k1 + 1)*wz*wy*wx

        end function trilinear

    end subroutine surface_shift_3d

    !
    !> Add to the meshes pm the mesh of the fault sf with the label id: a
    !> section at each sample of u, with a point at each row, within the padded
    !> model of n2 by n3 columns. A point is valid where the weight alpha is
    !> positive and the displacement along u is at least half a grid point (the
    !> fault causes no visible offset beyond).
    !
    subroutine surface_mesh_3d(sf, pm, id, n2, n3)

        type(flower_surface), intent(in) :: sf
        type(fault_mesh), intent(inout) :: pm
        integer, intent(in) :: id, n2, n3

        real, allocatable, dimension(:) :: q2, q3, qz, qd
        logical, allocatable, dimension(:) :: qv
        integer :: n1, nu, iu, i, j
        real :: u, v, du

        n1 = size(sf%x)
        nu = size(sf%g)
        allocate (q2(1:nu*n1), q3(1:nu*n1), qz(1:nu*n1), qd(1:nu*n1), qv(1:nu*n1))

        !$omp parallel do private(iu, i, j, u, v)
        do iu = 1, nu
            u = sf%u0 + (iu - 1.0)
            do i = 1, n1
                j = (iu - 1)*n1 + i
                v = sf%g(iu) + sf%s(iu)*sf%x(i)
                q3(j) = sf%c3 + u*sf%cphi - v*sf%sphi
                q2(j) = sf%c2 + u*sf%sphi + v*sf%cphi
                qz(j) = i - 1.0
                qv(j) = sf%alpha(i) > 0 .and. abs(sf%disp)*sf%a(iu) >= 0.5 .and. q2(j) >= -2.0 .and. q2(j) <= n2 + 1.0 &
                    .and. q3(j) >= -2.0 .and. q3(j) <= n3 + 1.0
                ! The displacement of the hanging wall, positive for a normal fault
                qd(j) = sf%disp*sf%alpha(i)*sf%a(iu)
            end do
        end do
        !$omp end parallel do

        ! The points of a section are one row apart: farther apart along a
        ! gently dipping part of the surface
        du = 1.0
        do i = 1, n1 - 1
            du = max(du, sqrt(1.0 + (maxval(abs(sf%s))*(sf%x(i + 1) - sf%x(i)))**2))
        end do
        call mesh_append(pm, q2, q3, qz, qd, qv, nu, n1, id, 0.5*rpi, 0.0, du, 1.0)

    end subroutine surface_mesh_3d

end module geological_model_flower
