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
! Fault labels of rgm2_curved and rgm3_curved carried as geometry.
!
! A fault whose displacement varies in space (a slip patch, a decay away from
! the fault) moves the rocks non-rigidly. Labels of the older faults that are
! moved from grid point to grid point then become irregular in width and
! break into scattered points. Here a fault is instead a chain of points in a
! 2D model and a mesh of points in a 3D model. Each younger fault moves the
! points with the rocks, and the labels are marked once at the end: the grid
! points within half the fault width of a chain or a mesh. The labels then
! keep their width however the faults are bent or cut.
!
module geological_model_label

    use libflit

    implicit none

    !
    !> Fault labels of a 2D model carried as chains of points: column x and
    !> depth z of each point (0-based, in the padded model), its displacement,
    !> its fault label and the index of its chain; n points in nchain chains
    !
    type fault_chain
        real, allocatable, dimension(:) :: x, z, disp
        integer, allocatable, dimension(:) :: id, ch
        integer :: n = 0, nchain = 0
    end type fault_chain

    !
    !> Fault labels of a 3D model carried as meshes of points. Mesh m has
    !> ns(m) sections of nu(m) points each, stored section by section from the
    !> index off(m) + 1 on: the columns x2 and x3 and the depth z of each point
    !> (0-based, in the padded model), its displacement, and whether it is
    !> valid (v: the fault is there and visible). id is the fault label, and du
    !> and ds are the spacings of the points along a section and between the
    !> sections when the mesh was made. rake is the rake angle for the strike
    !> direction on which the fault has a dip below (side = 1) or above
    !> (side = -1) pi/2; where a triangle of the mesh gets its strike label in
    !> the opposite direction (the strike label lies in [0, pi)), its dip label
    !> is on the other side of pi/2, and its rake label is pi - rake. n points
    !> in nm meshes.
    !
    type fault_mesh
        real, allocatable, dimension(:) :: x2, x3, z, disp
        logical, allocatable, dimension(:) :: v
        integer, allocatable, dimension(:) :: off, ns, nu, id
        real, allocatable, dimension(:) :: rake, side, du, ds
        integer :: n = 0, nm = 0
    end type fault_mesh

    !
    !> A fault of rgm3_curved, for moving fault meshes through its block shift
    !> and for making its own mesh. In map view, u runs along the base strike
    !> and v across it, from the fault center (xc along x3, yc along x2). The
    !> fault surface is v = gcurve(u) + devs(row): dphi and gcurve are the
    !> deviations of the strike angle and of the trace along u, sampled every
    !> grid point from umin on, and devs is the lateral deviation of the
    !> listric surface at each row. The block that shifts is on the positive
    !> side of the surface at the rows where pos is true. The displacement disp
    !> follows the slip patch (yn_patch: center uc, zc and half lengths au, az)
    !> and decays away from the surface (yn_decay: width wd).
    !
    type regular_fault_3d
        real :: xc = 0.0, yc = 0.0, strike = 0.0, dip = 0.0, rake = 0.0, disp = 0.0, umin = 0.0
        real, allocatable, dimension(:) :: dphi, gcurve, devs
        logical, allocatable, dimension(:) :: pos
        logical :: yn_patch = .false., yn_decay = .false.
        real :: uc = 0.0, au = 1.0, zc = 0.0, az = 1.0, wd = 1.0
    end type regular_fault_3d

    private
    public :: fault_chain
    public :: fault_mesh
    public :: regular_fault_3d
    public :: chain_append
    public :: chain_shift_2d
    public :: chain_fault_2d
    public :: chain_label_2d
    public :: mesh_append
    public :: mesh_split_3d
    public :: mesh_shift_3d
    public :: mesh_fault_3d
    public :: mesh_label_3d

contains

    !
    !> Add the points (x, z) with their displacements, labels and chain indices
    !> (1, 2, ..., counted within these points) to the chains pc
    !
    subroutine chain_append(pc, x, z, disp, id, ch)

        type(fault_chain), intent(inout) :: pc
        real, dimension(:), intent(in) :: x, z, disp
        integer, dimension(:), intent(in) :: id, ch

        if (.not. allocated(pc%x)) then
            allocate (pc%x(1:0), pc%z(1:0), pc%disp(1:0), pc%id(1:0), pc%ch(1:0))
            pc%n = 0
            pc%nchain = 0
        end if
        if (size(x) == 0) then
            return
        end if
        pc%x = [pc%x, x]
        pc%z = [pc%z, z]
        pc%disp = [pc%disp, disp]
        pc%id = [pc%id, id]
        pc%ch = [pc%ch, ch + pc%nchain]
        pc%n = pc%n + size(x)
        pc%nchain = pc%nchain + maxval(ch)

    end subroutine chain_append

    !
    !> Move the chains pc through the shift of a fault block of rgm2_curved.
    !> The block lies on one side of the fault trace b (0-based column at each
    !> row; left(i) = .true.: the block is toward -x2 at the row i). A grid
    !> point p in the block takes its rocks from the source s = p - u(p), where
    !> u is the displacement disp*alpha(p) along the dip angle dipn; alpha is
    !> the slip patch along depth (yn_patch: center row zc and half length az)
    !> times the Gaussian decay away from the trace (yn_decay: width wd).
    !>
    !> A point s of a chain therefore appears at the point p of the block with
    !> p = s + u(p), found by a fixed-point iteration, if that p lies in the
    !> block; and it stays where it is if it lies outside the block. Along a
    !> curved trace, the shift is not parallel to the trace everywhere, so a
    !> point near the trace can do both (the rocks there are duplicated) or
    !> neither (the rocks there are overridden). The points that appear in the
    !> block form new chains.
    !
    subroutine chain_shift_2d(pc, b, left, disp, dipn, yn_patch, zc, az, yn_decay, wd)

        type(fault_chain), intent(inout) :: pc
        real, dimension(:), intent(in) :: b
        logical, dimension(:), intent(in) :: left
        real, intent(in) :: disp, dipn, zc, az, wd
        logical, intent(in) :: yn_patch, yn_decay

        real, allocatable, dimension(:) :: tx, tz
        logical, allocatable, dimension(:) :: stay, move
        integer :: n1, np, l, it
        real :: z0, x0, z, x, zn, xn, a

        if (pc%n == 0) then
            return
        end if
        n1 = size(b)
        np = pc%n
        allocate (tx(1:np), tz(1:np), stay(1:np), move(1:np))

        !$omp parallel do private(l, z0, x0, z, x, zn, xn, it, a)
        do l = 1, np
            z0 = pc%z(l)
            x0 = pc%x(l)
            z = z0
            x = x0
            do it = 1, 50
                a = 1.0
                if (yn_patch) then
                    a = max(0.0, 1.0 - ((z + 1.0 - zc)/az)**2)
                end if
                if (yn_decay) then
                    a = a*exp(-0.5*((x - trace(z))/(wd + float_tiny))**2)
                end if
                zn = z0 + disp*a*sin(dipn)
                xn = x0 - disp*a*cos(dipn)
                if (abs(zn - z) + abs(xn - x) < 1.0e-3) then
                    z = zn
                    x = xn
                    exit
                end if
                z = zn
                x = xn
            end do
            tz(l) = z
            tx(l) = x
            move(l) = in_block(z, x)
            stay(l) = .not. in_block(z0, x0)
        end do
        !$omp end parallel do

        call chain_split_2d(pc, tx, tz, stay, move)

    contains

        ! Column of the trace at the fractional row z (0-based)
        pure function trace(z) result(v)

            real, intent(in) :: z
            real :: v

            integer :: i
            real :: w

            i = min(max(floor(z) + 1, 1), n1 - 1)
            w = min(max(z + 1.0 - i, 0.0), 1.0)
            v = (1.0 - w)*b(i) + w*b(i + 1)

        end function trace

        ! Whether the point (z, x) lies in the block
        pure function in_block(z, x) result(yn)

            real, intent(in) :: z, x
            logical :: yn

            if (left(min(max(nint(z) + 1, 1), n1))) then
                yn = x - trace(z) <= 0
            else
                yn = x - trace(z) >= 0
            end if

        end function in_block

    end subroutine chain_shift_2d

    !
    !> Move the chains pc through the shift of a fault block. The caller gives,
    !> for each point of pc, whether it stays where it is (it lies outside the
    !> block, whose rocks are replaced) and whether it appears in the block
    !> (move), at the point (tx, tz). A point near a curved fault can do both
    !> (the rocks there are duplicated) or neither (the rocks there are
    !> overridden). The points that stay keep their chains; the points that
    !> appear in the block form a new chain for each chain that they come from.
    !
    subroutine chain_split_2d(pc, tx, tz, stay, move)

        type(fault_chain), intent(inout) :: pc
        real, dimension(:), intent(in) :: tx, tz
        logical, dimension(:), intent(in) :: stay, move

        real, allocatable, dimension(:) :: nx, nz, nd
        integer, allocatable, dimension(:) :: nid, nch, newch
        integer :: np, l, nk, nm, j, nchain

        np = pc%n
        nk = count(stay(1:np))
        nm = count(move(1:np))
        allocate (nx(1:nk + nm), nz(1:nk + nm), nd(1:nk + nm), nid(1:nk + nm), nch(1:nk + nm))
        allocate (newch(1:max(pc%nchain, 1)))
        newch = 0
        nchain = pc%nchain
        j = 0
        do l = 1, np
            if (stay(l)) then
                j = j + 1
                nx(j) = pc%x(l)
                nz(j) = pc%z(l)
                nd(j) = pc%disp(l)
                nid(j) = pc%id(l)
                nch(j) = pc%ch(l)
            end if
        end do
        do l = 1, np
            if (move(l)) then
                if (newch(pc%ch(l)) == 0) then
                    nchain = nchain + 1
                    newch(pc%ch(l)) = nchain
                end if
                j = j + 1
                nx(j) = tx(l)
                nz(j) = tz(l)
                nd(j) = pc%disp(l)
                nid(j) = pc%id(l)
                nch(j) = newch(pc%ch(l))
            end if
        end do
        pc%x = nx
        pc%z = nz
        pc%disp = nd
        pc%id = nid
        pc%ch = nch
        pc%n = nk + nm
        pc%nchain = nchain

    end subroutine chain_split_2d

    !
    !> Add to the chains pc the points of a fault of rgm2_curved with the label
    !> id along its trace b (0-based column at each row), every half row. With
    !> yn_patch, the displacement disp follows the slip patch along depth
    !> (center row zc, half length az), and only the points where it is at
    !> least half a grid point are added: the fault causes no visible offset
    !> beyond. Each run of such points is a chain.
    !
    subroutine chain_fault_2d(pc, b, id, disp, yn_patch, zc, az)

        type(fault_chain), intent(inout) :: pc
        real, dimension(:), intent(in) :: b
        integer, intent(in) :: id
        real, intent(in) :: disp, zc, az
        logical, intent(in) :: yn_patch

        real, allocatable, dimension(:) :: qx, qz, qd
        integer, allocatable, dimension(:) :: qid, qch
        integer :: n1, nq, i, j, k, irun
        real :: z, a, w
        logical :: prev

        n1 = size(b)
        nq = 2*(n1 - 1) + 1
        allocate (qx(1:nq), qz(1:nq), qd(1:nq), qid(1:nq), qch(1:nq))
        j = 0
        irun = 0
        prev = .false.
        do i = 1, nq
            z = (i - 1)*0.5
            a = 1.0
            if (yn_patch) then
                a = max(0.0, 1.0 - ((z + 1.0 - zc)/az)**2)
            end if
            if (abs(disp)*a >= 0.5 .or. .not. yn_patch) then
                if (.not. prev) then
                    irun = irun + 1
                end if
                j = j + 1
                k = min(max(floor(z) + 1, 1), n1 - 1)
                w = min(max(z + 1.0 - k, 0.0), 1.0)
                qx(j) = (1.0 - w)*b(k) + w*b(k + 1)
                qz(j) = z
                qd(j) = disp*a
                qid(j) = id
                qch(j) = irun
                prev = .true.
            else
                prev = .false.
            end if
        end do
        call chain_append(pc, qx(1:j), qz(1:j), qd(1:j), qid(1:j), qch(1:j))

    end subroutine chain_fault_2d

    !
    !> Mark the labels of the chains pc in the fault arrays: the grid points
    !> within half the fault width fwidth of a chain take its label, dip and
    !> displacement
    !
    subroutine chain_label_2d(pc, fwidth, f, fdip, fdisp)

        type(fault_chain), intent(in) :: pc
        real, intent(in) :: fwidth
        real, dimension(:, :), intent(inout) :: f, fdip, fdisp

        if (pc%n > 0) then
            call chain_rasterize_2d(pc%x, pc%z, pc%disp, pc%id, pc%ch, fwidth, f, fdip, fdisp)
        end if

    end subroutine chain_label_2d

    !
    !> Mark the grid points within half the fault width fwidth of the chains of
    !> points (x, z) (0-based column and depth) with the label id of the nearest
    !> chain, its dip in (0, pi) (below pi/2 the fault dips toward -x2), and
    !> its displacement disp. Consecutive points of a chain (same ch) farther
    !> apart than 2.5 grid points are not connected: a younger fault cut the
    !> chain there.
    !
    subroutine chain_rasterize_2d(x, z, disp, id, ch, fwidth, f, fdip, fdisp)

        real, dimension(:), intent(in) :: x, z, disp
        integer, dimension(:), intent(in) :: id, ch
        real, intent(in) :: fwidth
        real, dimension(:, :), intent(inout) :: f, fdip, fdisp

        real, allocatable, dimension(:, :) :: best
        integer :: n1, n2, i, j, l, i0, i1, j0, j1
        real :: r, dx, dz, len2, a, d, gx, gz

        n1 = size(f, 1)
        n2 = size(f, 2)
        r = 0.5*fwidth
        best = zeros(n1, n2) + r

        do l = 1, size(x) - 1
            if (ch(l + 1) /= ch(l)) then
                cycle
            end if
            dx = x(l + 1) - x(l)
            dz = z(l + 1) - z(l)
            len2 = dx**2 + dz**2
            if (len2 > 2.5**2 .or. len2 <= 0) then
                cycle
            end if
            j0 = max(floor(min(x(l), x(l + 1)) - r), 0) + 1
            j1 = min(ceiling(max(x(l), x(l + 1)) + r), n2 - 1) + 1
            i0 = max(floor(min(z(l), z(l + 1)) - r), 0) + 1
            i1 = min(ceiling(max(z(l), z(l + 1)) + r), n1 - 1) + 1
            do j = j0, j1
                do i = i0, i1
                    gx = j - 1.0 - x(l)
                    gz = i - 1.0 - z(l)
                    a = min(max((gx*dx + gz*dz)/len2, 0.0), 1.0)
                    d = sqrt((gx - a*dx)**2 + (gz - a*dz)**2)
                    if (d < best(i, j)) then
                        best(i, j) = d
                        f(i, j) = id(l)
                        fdip(i, j) = modulo(atan2(-dz, dx), real(const_pi))
                        fdisp(i, j) = (1.0 - a)*disp(l) + a*disp(l + 1)
                    end if
                end do
            end do
        end do

    end subroutine chain_rasterize_2d

    !
    !> Add a mesh of ns sections of nu points (x2, x3, z, displacement, valid
    !> or not) to the meshes pm, with its label, rake (and the dip side that the
    !> rake refers to; 0 for a rake of pi/2), and point spacings
    !
    subroutine mesh_append(pm, x2, x3, z, disp, v, ns, nu, id, rake, side, du, ds)

        type(fault_mesh), intent(inout) :: pm
        real, dimension(:), intent(in) :: x2, x3, z, disp
        logical, dimension(:), intent(in) :: v
        integer, intent(in) :: ns, nu, id
        real, intent(in) :: rake, side, du, ds

        if (.not. allocated(pm%x2)) then
            allocate (pm%x2(1:0), pm%x3(1:0), pm%z(1:0), pm%disp(1:0), pm%v(1:0))
            allocate (pm%off(1:0), pm%ns(1:0), pm%nu(1:0), pm%id(1:0), pm%rake(1:0), pm%side(1:0), pm%du(1:0), pm%ds(1:0))
            pm%n = 0
            pm%nm = 0
        end if
        pm%off = [pm%off, pm%n]
        pm%ns = [pm%ns, ns]
        pm%nu = [pm%nu, nu]
        pm%id = [pm%id, id]
        pm%rake = [pm%rake, rake]
        pm%side = [pm%side, side]
        pm%du = [pm%du, du]
        pm%ds = [pm%ds, ds]
        pm%x2 = [pm%x2, x2]
        pm%x3 = [pm%x3, x3]
        pm%z = [pm%z, z]
        pm%disp = [pm%disp, disp]
        pm%v = [pm%v, v]
        pm%n = pm%n + ns*nu
        pm%nm = pm%nm + 1

    end subroutine mesh_append

    !
    !> Move the meshes pm through the shift of a fault block. The caller gives,
    !> for each point of pm, whether it stays where it is (it lies outside the
    !> block, whose rocks are replaced) and whether it appears in the block
    !> (move), at the point (t2, t3, tz). Along a curved fault a point near the
    !> fault can do both (the rocks there are duplicated) or neither (the rocks
    !> there are overridden). A mesh with points of both kinds is split: the
    !> points that stay keep the mesh, and the points that appear in the block
    !> form a new mesh with the same label.
    !
    subroutine mesh_split_3d(pm, t2, t3, tz, stay, move)

        type(fault_mesh), intent(inout) :: pm
        real, dimension(:), intent(in) :: t2, t3, tz
        logical, dimension(:), intent(in) :: stay, move

        integer :: m, nm0, o, n
        logical :: yn_stay, yn_move

        nm0 = pm%nm
        do m = 1, nm0
            o = pm%off(m)
            n = pm%ns(m)*pm%nu(m)
            yn_stay = any(pm%v(o + 1:o + n) .and. stay(o + 1:o + n))
            yn_move = any(pm%v(o + 1:o + n) .and. move(o + 1:o + n))
            if (yn_move .and. yn_stay) then
                call mesh_append(pm, t2(o + 1:o + n), t3(o + 1:o + n), tz(o + 1:o + n), pm%disp(o + 1:o + n), &
                    pm%v(o + 1:o + n) .and. move(o + 1:o + n), pm%ns(m), pm%nu(m), pm%id(m), pm%rake(m), pm%side(m), &
                    pm%du(m), pm%ds(m))
                pm%v(o + 1:o + n) = pm%v(o + 1:o + n) .and. stay(o + 1:o + n)
            else if (yn_move) then
                pm%x2(o + 1:o + n) = t2(o + 1:o + n)
                pm%x3(o + 1:o + n) = t3(o + 1:o + n)
                pm%z(o + 1:o + n) = tz(o + 1:o + n)
                pm%v(o + 1:o + n) = pm%v(o + 1:o + n) .and. move(o + 1:o + n)
            else
                pm%v(o + 1:o + n) = pm%v(o + 1:o + n) .and. stay(o + 1:o + n)
            end if
        end do

    end subroutine mesh_split_3d

    !
    !> Move the meshes pm through the block shift of the fault rf of
    !> rgm3_curved. A grid point p in the block takes its rocks from the source
    !> s = p - u(p), where u is the local displacement along the slip vector
    !> (which follows the local strike). A point s of a mesh therefore appears
    !> at the point p of the block with p = s + u(p), found by a fixed-point
    !> iteration, if that p lies in the block; and it stays where it is if it
    !> lies outside the block (see mesh_split_3d).
    !
    subroutine mesh_shift_3d(pm, rf)

        type(fault_mesh), intent(inout) :: pm
        type(regular_fault_3d), intent(in) :: rf

        real, allocatable, dimension(:) :: t2, t3, tz
        logical, allocatable, dimension(:) :: stay, move
        integer :: np, l, it, n1, nu
        real :: cphi, sphi, cr, sr, cd, s1, z0, y0, x0, z, y, x, zn, yn, xn, dp, vdist, a

        if (pm%n == 0) then
            return
        end if
        np = pm%n
        n1 = size(rf%devs)
        nu = size(rf%dphi)
        cphi = cos(rf%strike)
        sphi = sin(rf%strike)
        cr = cos(rf%rake)
        sr = sin(rf%rake)
        cd = cos(rf%dip)
        s1 = sr*sin(rf%dip)
        allocate (t2(1:np), t3(1:np), tz(1:np), stay(1:np), move(1:np))

        !$omp parallel do private(l, it, z0, y0, x0, z, y, x, zn, yn, xn, dp, vdist, a)
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
            call local(z0, y0, x0, dp, vdist, a)
            stay(l) = .not. in_block(z0, vdist)
            z = z0
            y = y0
            x = x0
            do it = 1, 50
                call local(z, y, x, dp, vdist, a)
                zn = z0 + s1*rf%disp*a
                yn = y0 - (cr*sin(rf%strike + dp) - sr*cd*cos(rf%strike + dp))*rf%disp*a
                xn = x0 - (cr*cos(rf%strike + dp) + sr*cd*sin(rf%strike + dp))*rf%disp*a
                if (abs(zn - z) + abs(yn - y) + abs(xn - x) < 1.0e-3) then
                    z = zn
                    y = yn
                    x = xn
                    exit
                end if
                z = zn
                y = yn
                x = xn
            end do
            call local(z, y, x, dp, vdist, a)
            move(l) = in_block(z, vdist)
            tz(l) = z
            t2(l) = y
            t3(l) = x
        end do
        !$omp end parallel do

        call mesh_split_3d(pm, t2, t3, tz, stay, move)

    contains

        ! Strike deviation, signed distance to the fault surface and
        ! displacement weight at the point (z, y along x2, x along x3)
        pure subroutine local(z, y, x, dp, vdist, a)

            real, intent(in) :: z, y, x
            real, intent(out) :: dp, vdist, a

            real :: u0, v0, wu, gc, w, dev
            integer :: iu, i

            u0 = (x - rf%xc)*cphi + (y - rf%yc)*sphi
            v0 = -(x - rf%xc)*sphi + (y - rf%yc)*cphi
            iu = min(max(floor(u0 - rf%umin) + 1, 1), nu - 1)
            wu = min(max(u0 - rf%umin - (iu - 1.0), 0.0), 1.0)
            dp = rf%dphi(iu)*(1.0 - wu) + rf%dphi(iu + 1)*wu
            gc = rf%gcurve(iu)*(1.0 - wu) + rf%gcurve(iu + 1)*wu
            i = min(max(floor(z) + 1, 1), n1 - 1)
            w = min(max(z + 1.0 - i, 0.0), 1.0)
            dev = rf%devs(i)*(1.0 - w) + rf%devs(i + 1)*w
            vdist = (v0 - gc - dev)*cos(dp)
            a = 1.0
            if (rf%yn_patch) then
                a = max(0.0, 1.0 - ((u0 - rf%uc)/rf%au)**2 - ((z + 1.0 - rf%zc)/rf%az)**2)
            end if
            if (rf%yn_decay) then
                a = a*exp(-0.5*(vdist/(rf%wd + float_tiny))**2)
            end if

        end subroutine local

        ! Whether a point at the depth z and the signed distance vdist lies in the block
        pure function in_block(z, vdist) result(yn)

            real, intent(in) :: z, vdist
            logical :: yn

            if (rf%pos(min(max(nint(z) + 1, 1), n1))) then
                yn = vdist > 0
            else
                yn = vdist < 0
            end if

        end function in_block

    end subroutine mesh_shift_3d

    !
    !> Add to the meshes pm the mesh of the fault rf of rgm3_curved with the
    !> label id: a point on the fault surface every grid point along the base
    !> strike (a section) and every row, within the padded model of n2 by n3
    !> columns. With a slip patch, a point is valid where the local
    !> displacement is at least half a grid point: the fault causes no visible
    !> offset beyond.
    !
    subroutine mesh_fault_3d(pm, rf, id, n2, n3)

        type(fault_mesh), intent(inout) :: pm
        type(regular_fault_3d), intent(in) :: rf
        integer, intent(in) :: id, n2, n3

        real, allocatable, dimension(:) :: q2, q3, qz, qd
        logical, allocatable, dimension(:) :: qv
        integer :: n1, nu, iu, i, j
        real :: cphi, sphi, u, v, a, du

        n1 = size(rf%devs)
        nu = size(rf%dphi)
        cphi = cos(rf%strike)
        sphi = sin(rf%strike)
        allocate (q2(1:nu*n1), q3(1:nu*n1), qz(1:nu*n1), qd(1:nu*n1), qv(1:nu*n1))

        !$omp parallel do private(iu, i, j, u, v, a)
        do iu = 1, nu
            u = rf%umin + (iu - 1.0)
            do i = 1, n1
                j = (iu - 1)*n1 + i
                v = rf%gcurve(iu) + rf%devs(i)
                q3(j) = rf%xc + u*cphi - v*sphi
                q2(j) = rf%yc + u*sphi + v*cphi
                qz(j) = i - 1.0
                a = 1.0
                if (rf%yn_patch) then
                    a = max(0.0, 1.0 - ((u - rf%uc)/rf%au)**2 - ((i - rf%zc)/rf%az)**2)
                end if
                qd(j) = rf%disp*a
                qv(j) = q2(j) >= -2.0 .and. q2(j) <= n2 + 1.0 .and. q3(j) >= -2.0 .and. q3(j) <= n3 + 1.0
                if (rf%yn_patch) then
                    qv(j) = qv(j) .and. abs(rf%disp)*a >= 0.5
                end if
            end do
        end do
        !$omp end parallel do

        ! The points of a section are one row apart: farther apart along a
        ! gently dipping part of the surface
        du = 1.0
        do i = 1, n1 - 1
            du = max(du, sqrt(1.0 + (rf%devs(i + 1) - rf%devs(i))**2))
        end do
        call mesh_append(pm, q2, q3, qz, qd, qv, nu, n1, id, rf%rake, sign(1.0, real(0.5*const_pi) - rf%dip), du, 1.0)

    end subroutine mesh_fault_3d

    !
    !> Mark the labels of the meshes pm in the fault arrays: the grid points
    !> within half the fault width fwidth of a mesh take its label and rake,
    !> and the strike, dip and displacement of the nearest triangle
    !
    subroutine mesh_label_3d(pm, fwidth, f, fdip, fstrike, frake, fdisp)

        type(fault_mesh), intent(in) :: pm
        real, intent(in) :: fwidth
        real, dimension(:, :, :), intent(inout) :: f, fdip, fstrike, frake, fdisp

        if (pm%nm > 0) then
            call mesh_rasterize_3d(pm%x2, pm%x3, pm%z, pm%disp, pm%v, pm%off, pm%ns, pm%nu, pm%id, pm%du, pm%ds, pm%rake, &
                pm%side, fwidth, f, fdip, fstrike, frake, fdisp)
        end if

    end subroutine mesh_label_3d

    !
    !> Mark the grid points within half the fault width fwidth of the fault
    !> meshes with the label mid of the nearest mesh, its rake (mrake and mside,
    !> see fault_mesh), and the strike, dip and displacement of the nearest
    !> triangle. Mesh m has mns(m) sections of mnu(m) points (x2, x3, z:
    !> 0-based columns and depth) from the index moff(m) + 1 on, spaced by
    !> mdu(m) along a section and by mds(m) between the sections; a point is
    !> valid (v) where the fault is visible. A quad of four valid points is
    !> split into two triangles, unless a younger fault cut it: an edge along a
    !> section longer than 4*mdu, or between the sections longer than 3*mds.
    !
    subroutine mesh_rasterize_3d(x2, x3, z, disp, v, moff, mns, mnu, mid, mdu, mds, mrake, mside, fwidth, &
            f, fdip, fstrike, frake, fdisp)

        real, dimension(:), intent(in) :: x2, x3, z, disp
        logical, dimension(:), intent(in) :: v
        integer, dimension(:), intent(in) :: moff, mns, mnu, mid
        real, dimension(:), intent(in) :: mdu, mds, mrake, mside
        real, intent(in) :: fwidth
        real, dimension(:, :, :), intent(inout) :: f, fdip, fstrike, frake, fdisp

        real, allocatable, dimension(:, :, :) :: best
        integer :: n1, n2, n3, mm, is, iu, a, b, c, d
        real :: r

        n1 = size(f, 1)
        n2 = size(f, 2)
        n3 = size(f, 3)
        r = 0.5*fwidth
        best = zeros(n1, n2, n3) + r

        do mm = 1, size(moff)
            do is = 1, mns(mm) - 1
                do iu = 1, mnu(mm) - 1
                    a = moff(mm) + (is - 1)*mnu(mm) + iu
                    b = a + 1
                    c = a + mnu(mm)
                    d = c + 1
                    if (.not. (v(a) .and. v(b) .and. v(c) .and. v(d))) then
                        cycle
                    end if
                    if (edge(a, b) > 4.0*mdu(mm) .or. edge(c, d) > 4.0*mdu(mm) .or. edge(a, c) > 3.0*mds(mm) &
                            .or. edge(b, d) > 3.0*mds(mm)) then
                        cycle
                    end if
                    call triangle(a, b, d, mid(mm), mrake(mm), mside(mm))
                    call triangle(a, d, c, mid(mm), mrake(mm), mside(mm))
                end do
            end do
        end do

    contains

        pure function edge(i, j) result(e)

            integer, intent(in) :: i, j
            real :: e

            e = sqrt((x2(i) - x2(j))**2 + (x3(i) - x3(j))**2 + (z(i) - z(j))**2)

        end function edge

        !
        !> Rasterize the triangle of the points ia, ib, ic: for each grid point
        !> near it, the closest point on the triangle (by the regions of its
        !> barycentric coordinates)
        !
        subroutine triangle(ia, ib, ic, id, rake, side)

            integer, intent(in) :: ia, ib, ic, id
            real, intent(in) :: rake, side

            real, dimension(1:3) :: pa, pb, pc, ab, ac, ap, bp, cp, q, nrm
            real :: d1, d2, d3, d4, d5, d6, va, vb, vc, wa, wb, wc, den, dd, strike, dipl
            integer :: i, j, l, i0, i1, j0, j1, l0, l1

            pa = [z(ia), x2(ia), x3(ia)]
            pb = [z(ib), x2(ib), x3(ib)]
            pc = [z(ic), x2(ic), x3(ic)]
            ab = pb - pa
            ac = pc - pa
            nrm = [ab(2)*ac(3) - ab(3)*ac(2), ab(3)*ac(1) - ab(1)*ac(3), ab(1)*ac(2) - ab(2)*ac(1)]
            if (sum(nrm**2) <= 1.0e-12) then
                return
            end if
            call normal_strike_dip(nrm(1), nrm(2), nrm(3), strike, dipl)

            i0 = max(floor(min(pa(1), pb(1), pc(1)) - r), 0) + 1
            i1 = min(ceiling(max(pa(1), pb(1), pc(1)) + r), n1 - 1) + 1
            j0 = max(floor(min(pa(2), pb(2), pc(2)) - r), 0) + 1
            j1 = min(ceiling(max(pa(2), pb(2), pc(2)) + r), n2 - 1) + 1
            l0 = max(floor(min(pa(3), pb(3), pc(3)) - r), 0) + 1
            l1 = min(ceiling(max(pa(3), pb(3), pc(3)) + r), n3 - 1) + 1

            do l = l0, l1
                do j = j0, j1
                    do i = i0, i1
                        q = [i - 1.0, j - 1.0, l - 1.0]
                        ap = q - pa
                        d1 = dot_product(ab, ap)
                        d2 = dot_product(ac, ap)
                        bp = q - pb
                        d3 = dot_product(ab, bp)
                        d4 = dot_product(ac, bp)
                        cp = q - pc
                        d5 = dot_product(ab, cp)
                        d6 = dot_product(ac, cp)
                        vc = d1*d4 - d3*d2
                        vb = d5*d2 - d1*d6
                        va = d3*d6 - d5*d4
                        if (d1 <= 0 .and. d2 <= 0) then
                            wa = 1.0
                            wb = 0.0
                            wc = 0.0
                        else if (d3 >= 0 .and. d4 <= d3) then
                            wa = 0.0
                            wb = 1.0
                            wc = 0.0
                        else if (vc <= 0 .and. d1 >= 0 .and. d3 <= 0) then
                            wb = d1/(d1 - d3)
                            wa = 1.0 - wb
                            wc = 0.0
                        else if (d6 >= 0 .and. d5 <= d6) then
                            wa = 0.0
                            wb = 0.0
                            wc = 1.0
                        else if (vb <= 0 .and. d2 >= 0 .and. d6 <= 0) then
                            wc = d2/(d2 - d6)
                            wa = 1.0 - wc
                            wb = 0.0
                        else if (va <= 0 .and. d4 - d3 >= 0 .and. d5 - d6 >= 0) then
                            wc = (d4 - d3)/((d4 - d3) + (d5 - d6))
                            wb = 1.0 - wc
                            wa = 0.0
                        else
                            den = 1.0/(va + vb + vc)
                            wb = vb*den
                            wc = vc*den
                            wa = 1.0 - wb - wc
                        end if
                        dd = sqrt(sum((q - (wa*pa + wb*pb + wc*pc))**2))
                        if (dd < best(i, j, l)) then
                            best(i, j, l) = dd
                            f(i, j, l) = id
                            fstrike(i, j, l) = strike
                            fdip(i, j, l) = dipl
                            if ((0.5*const_pi - dipl)*side < 0) then
                                frake(i, j, l) = const_pi - rake
                            else
                                frake(i, j, l) = rake
                            end if
                            fdisp(i, j, l) = wa*disp(ia) + wb*disp(ib) + wc*disp(ic)
                        end if
                    end do
                end do
            end do

        end subroutine triangle

    end subroutine mesh_rasterize_3d

    !
    !> Strike in [0, pi) and dip in [0, pi) in the convention of rgm3_curved
    !> (strike vector (sin(strike), cos(strike)) in (x2, x3); a dip below pi/2
    !> dips toward (cos(strike), -sin(strike))) of a plane with the normal
    !> (a1, a2, a3) along (x1, x2, x3), x1 positive downward
    !
    pure subroutine normal_strike_dip(a1, a2, a3, strike, dipl)

        real, intent(in) :: a1, a2, a3
        real, intent(out) :: strike, dipl

        real :: b1, b2, b3, h, del

        ! The upward normal leans toward the dip direction
        b1 = a1
        b2 = a2
        b3 = a3
        if (b1 > 0) then
            b1 = -b1
            b2 = -b2
            b3 = -b3
        end if
        h = sqrt(b2**2 + b3**2)
        del = atan2(h, -b1)
        if (h < 1.0e-6) then
            strike = 0.0
            dipl = 0.0
            return
        end if
        ! The strike vector (sin(strike), cos(strike)) is normal to the dip
        ! direction; a dip below pi/2 dips toward (cos(strike), -sin(strike))
        strike = modulo(atan2(-b3, b2), real(const_pi))
        if (cos(strike)*b2 - sin(strike)*b3 > 0) then
            dipl = del
        else
            dipl = real(const_pi) - del
        end if

    end subroutine normal_strike_dip

end module geological_model_label
