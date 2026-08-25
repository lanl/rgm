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

!==============================================================================
! Flow-accumulation drainage-network generator
!
! Algorithm
! ---------
! 1. Terrain: fractal Brownian motion (fBm) surface built by summing octaves
!    of bilinear-upsampled random noise.  Octave k has grid size min(nx,4·2^k)
!    and amplitude 2^{−k·hurst}.  Higher hurst → smoother, fewer tributaries;
!    lower hurst → rougher, denser branching.  A linear regional tilt is added
!    so water drains consistently toward one side (iy = 1).
!
! 2. Pit filling: Planchon-Darboux iterative forward + backward sweeps raise
!    interior depressions by eps = 1e-6 per pass until every cell can drain to
!    a boundary cell.  This eliminates flat-floored pits that would stall D8
!    routing.
!
! 3. D8 flow routing: each cell is assigned the direction (out of 8 compass
!    directions) of its steepest downslope neighbour on the filled DEM.
!    Diagonal neighbours use distance √2 in the slope calculation.
!
! 4. Topological accumulation: cells are sorted high→low by elevation; each
!    cell passes its accumulated flow count to its D8 downstream neighbour.
!    Result: acc(ix,iy) = number of cells draining through (ix,iy), so the
!    main trunk has the highest acc value.
!
! 5. Threshold: a geometric binary search finds the smallest accumulation
!    threshold such that exactly channel_frac · nx · ny cells qualify as
!    channel cells.
!
! 6. Reach decomposition: the network is split into reaches at headwaters
!    (no upstream channel neighbours) and confluences (≥ 2 upstream channel
!    neighbours).  Each reach is smoothed with a 3-point weighted average
!    (endpoints fixed at their exact grid positions so junctions are gapless),
!    then rendered as a sequence of capsule segments.
!    Width (and depth in 3-D) scale as W_max · √(acc / acc_max), so the trunk
!    is widest and tributaries taper toward their headwaters.
!
! 7. Canyon mode: instead of rendering capsule tubes, the distance from each
!    XY cell to the nearest channel axis is computed.  A depth fraction
!      df = (1 − dist / r_local)^canyon_exp
!    is derived for each cell (r_local ∝ √acc).  The cell is then filled
!    erosionally from iz_min = nz − round((nz−1) · df^wall_exp) up to iz = nz,
!    producing a body that is deepest at the thalweg and tapers to zero at
!    the canyon rim.
!
! Public OOP API
! --------------
!   type(drainage_channel) :: dc
!   call dc%generate     ! fills dc%depth_map(nx,ny) and dc%vol(nx,ny,nz)
!
!   type(drainage_canyon) :: dc
!   call dc%generate     ! fills dc%depth_map(nx,ny) and dc%vol(nx,ny,nz)
!
!   depth_map: 0 = background; value = depth of channel floor in iz-units (0..nz-1)
!   vol:       1.0 = channel / canyon body, 0.0 = host rock
!==============================================================================

module geological_model_drainage

    use libflit
    use geological_model_utility
    use, intrinsic :: iso_fortran_env, only: dp => real64, sp => real32
    implicit none
    private

    public :: drainage_channel, drainage_canyon

    ! D8 neighbour offsets (E NE N NW W SW S SE)
    integer, parameter :: DX8(8) = [1, 1, 0, -1, -1, -1, 0, 1]
    integer, parameter :: DY8(8) = [0, 1, 1, 1, 0, -1, -1, -1]
    real(dp), parameter :: DIST8(8) = [1.0_dp, 1.41421356_dp, 1.0_dp, 1.41421356_dp, &
        1.0_dp, 1.41421356_dp, 1.0_dp, 1.41421356_dp]

    !> Dendritic drainage channel network: plan-view depth map + 3-D incised body.
    !> generate() fills depth_map(ny,nx) (downward-positive) and vol/obj(nnz,ny,nx).
    type drainage_channel
        integer  :: nx = 256, ny = 256  ! output grid dimensions (map view)
        integer  :: nz = 64             ! vertical layers (incision depth)
        integer  :: nnz = 0
        real(dp) :: W_max = 10.0_dp     ! trunk channel width in grid cells (tributaries scale as √acc/acc_max)
        real(dp) :: D_max = 32.0_dp     ! trunk channel depth in layers (tributaries scale as √acc/acc_max)
        real(dp) :: channel_frac = 0.05_dp ! fraction of grid cells rendered as channel (0.02=sparse trunk, 0.1=dense network)
        real(dp) :: hurst = 0.75_dp     ! fBm Hurst exponent: 0.5=rough/branchy terrain, 0.85=smooth/few tributaries
        real(dp) :: tilt = 3.0_dp       ! regional slope magnitude: high (3–5) → single trunk, low (1–2) → competing trunks
        real(dp) :: terrain_bg = 0.0_dp ! if > 0, fBm terrain (normalised to [0, terrain_bg]) fills non-channel cells in depth_map
        integer  :: seed = -1           ! random seed; also controls drain direction (one of 4 cardinal edges) unless orient is set
        integer  :: orient = -1         ! drain direction: -1 = derived from seed; 0-3 = fixed (2 = toward +x, i.e., flow along x)
        logical  :: yn_depth_only = .false. ! if .true., skip rebuilding 3D vol/obj; only depth_map is produced
        real(sp), allocatable :: depth_map(:, :)  ! 0 = background (or terrain if terrain_bg > 0); value = channel-floor depth (0..nz-1)
        real(sp), allocatable :: vol(:, :, :)     ! 1.0 = channel sand body, 0.0 = host rock
        real(sp), allocatable :: obj(:, :, :)
    contains
        procedure :: generate => generate_drainage_channel
    end type drainage_channel

    !> Submarine / incised canyon carved by a dendritic drainage network.
    !> generate() fills depth_map(ny,nx) (downward-positive) and vol/obj(nnz,ny,nx).
    type drainage_canyon
        integer  :: nx = 256, ny = 256  ! output grid dimensions (map view)
        integer  :: nz = 64             ! vertical layers (canyon depth)
        real(dp) :: W_max = 20.0_dp     ! canyon full-width at surface for the main trunk (grid cells)
        real(dp) :: channel_frac = 0.015_dp ! channel fraction: 0.01–0.02 → trunk + 1–2 tributaries; 0.05 → more branching
        real(dp) :: hurst = 0.75_dp     ! fBm Hurst exponent controlling terrain roughness
        real(dp) :: tilt = 4.0_dp       ! regional slope (higher → single dominant trunk)
        real(dp) :: canyon_exp = 1.0_dp ! cross-section shape exponent: 1=linear V, 0.5=U-shape, 2=narrow-V
        real(dp) :: wall_exp = 1.0_dp   ! vertical wall curvature: 1=linear, <1=concave/flared, >1=convex/steep-stem
        real(dp) :: terrain_bg = 0.0_dp ! if > 0, fBm terrain fills non-canyon cells; canyon cells show terrain − canyon_depth
        integer  :: seed = -1           ! random seed; also controls drain direction (one of 4 cardinal edges) unless orient is set
        integer  :: orient = -1         ! drain direction: -1 = derived from seed; 0-3 = fixed (2 = toward +x, i.e., flow along x)
        logical  :: yn_depth_only = .false. ! if .true., skip rebuilding 3D vol/obj; only depth_map is produced
        integer  :: nnz = 0
        real(sp), allocatable :: depth_map(:, :)  ! downward-positive depth map
        real(sp), allocatable :: vol(:, :, :)     ! 1.0 = solid below surface, 0.0 = void above
        real(sp), allocatable :: obj(:, :, :)
    contains
        procedure :: generate => generate_drainage_canyon
    end type drainage_canyon

contains

    ! ── bilinear upsample: coarse (nx_s × ny_s) → fine (nx_t × ny_t) ──────────

    subroutine bilinear_up(src, nx_s, ny_s, nx_t, ny_t, dst)
        integer, intent(in)  :: nx_s, ny_s, nx_t, ny_t
        real(dp), intent(in)  :: src(nx_s, ny_s)
        real(dp), intent(out) :: dst(nx_t, ny_t)

        real(dp) :: fx, fy, wx, wy
        integer  :: ix, iy, ix0, iy0

        do iy = 1, ny_t
            fy = 1.0_dp + real(iy - 1, dp)/real(ny_t - 1, dp)*real(ny_s - 1, dp)
            iy0 = max(1, min(ny_s - 1, int(fy)))
            wy = fy - real(iy0, dp)
            do ix = 1, nx_t
                fx = 1.0_dp + real(ix - 1, dp)/real(nx_t - 1, dp)*real(nx_s - 1, dp)
                ix0 = max(1, min(nx_s - 1, int(fx)))
                wx = fx - real(ix0, dp)
                dst(ix, iy) = (1 - wx)*(1 - wy)*src(ix0, iy0) &
                    + wx*(1 - wy)*src(ix0 + 1, iy0) &
                    + (1 - wx)*wy*src(ix0, iy0 + 1) &
                    + wx*wy*src(ix0 + 1, iy0 + 1)
            end do
        end do
    end subroutine bilinear_up

    ! ── fBm terrain: sum of random octaves, coarse-to-fine ────────────────────
    ! Octave k uses a nc = min(nx, 4·2^k) random grid, amplitude = 2^{-k·hurst}.
    ! A regional tilt drives drainage toward one of 4 cardinal edges.
    ! drain_dir: 1=toward iy=ny (+y), 2=toward iy=1 (-y),
    !            3=toward ix=nx (+x), 4=toward ix=1 (-x).

    subroutine gen_terrain(nx, ny, hurst, n_oct, tilt, seed, dem, drain_dir)
        integer, intent(in)           :: nx, ny, n_oct, seed
        integer, intent(in), optional :: drain_dir
        real(dp), intent(in)          :: hurst, tilt
        real(dp), intent(out)         :: dem(nx, ny)

        real(dp), allocatable :: rn(:), coarse(:, :), fine(:, :)
        real(dp) :: amp, drange
        integer  :: oct, nc, s, ix, iy, dir_val

        dem = 0.0_dp
        amp = 1.0_dp
        s = mod(abs(seed), 2147483647)

        do oct = 0, n_oct - 1
            nc = min(nx, max(4, 4*(2**oct)))    ! 4, 8, 16, 32, … up to nx
            allocate (coarse(nc, nc), fine(nx, ny))
            rn = random(nc*nc, dist='uniform', &
                seed=safe_seed(abs(int(s, 8)*(oct + 1)*1009 + 7)))
            coarse = reshape(rn, [nc, nc])
            call bilinear_up(coarse, nc, nc, nx, ny, fine)
            dem = dem + amp*fine
            amp = amp*(0.5_dp**hurst)
            deallocate (coarse, fine, rn)
        end do

        ! Normalise to [0,1]
        drange = maxval(dem) - minval(dem)
        if (drange > 1.0e-12_dp) then
            dem = (dem - minval(dem))/drange
        end if

        dir_val = 1
        if (present(drain_dir)) dir_val = drain_dir
        select case (dir_val)
            case (1)  ! drain toward iy=ny; source at iy=1
                do iy = 1, ny
                    dem(:, iy) = dem(:, iy) + tilt*real(ny - iy, dp)/real(max(1, ny - 1), dp)
                end do
            case (2)  ! drain toward iy=1; source at iy=ny
                do iy = 1, ny
                    dem(:, iy) = dem(:, iy) + tilt*real(iy - 1, dp)/real(max(1, ny - 1), dp)
                end do
            case (3)  ! drain toward ix=nx; source at ix=1
                do ix = 1, nx
                    dem(ix, :) = dem(ix, :) + tilt*real(nx - ix, dp)/real(max(1, nx - 1), dp)
                end do
            case default  ! drain toward ix=1; source at ix=nx
                do ix = 1, nx
                    dem(ix, :) = dem(ix, :) + tilt*real(ix - 1, dp)/real(max(1, nx - 1), dp)
                end do
        end select
    end subroutine gen_terrain

    ! ── quicksort on index array, descending by key ────────────────────────────

    recursive subroutine qsort(idx, key, lo, hi)
        integer, intent(inout) :: idx(:)
        real(dp), intent(in)    :: key(:)
        integer, intent(in)    :: lo, hi

        integer  :: i, j, tmp
        real(dp) :: pval

        if (lo >= hi) then
            return
        end if
        pval = key(idx((lo + hi)/2))
        i = lo
        j = hi
        do while (i <= j)
            do while (key(idx(i)) > pval)
                i = i + 1
            end do
            do while (key(idx(j)) < pval)
                j = j - 1
            end do
            if (i <= j) then
                tmp = idx(i)
                idx(i) = idx(j)
                idx(j) = tmp
                i = i + 1
                j = j - 1
            end if
        end do
        if (lo < j) then
            call qsort(idx, key, lo, j)
        end if
        if (i < hi) then
            call qsort(idx, key, i, hi)
        end if
    end subroutine qsort

    ! ── Planchon-Darboux pit-filling ──────────────────────────────────────────
    ! Raises each pit just enough (by eps per step) so every cell drains to a
    ! boundary cell.  Forward + backward sweeps converge in O(ny) iterations.

    subroutine pitfill(nx, ny, dem, filled)
        integer, intent(in)  :: nx, ny
        real(dp), intent(in)  :: dem(nx, ny)
        real(dp), intent(out) :: filled(nx, ny)

        real(dp), parameter :: eps = 1.0e-6_dp
        logical  :: changed
        integer  :: ix, iy, jx, jy, d
        real(dp) :: new_val

        ! Borders keep their true elevation; interior starts at huge
        filled = huge(0.0_dp)
        do iy = 1, ny
            do ix = 1, nx
                if (ix == 1 .or. ix == nx .or. iy == 1 .or. iy == ny) then
                    filled(ix, iy) = dem(ix, iy)
                end if
            end do
        end do

        do
            changed = .false.
            ! Forward sweep
            do iy = 1, ny
                do ix = 1, nx
                    if (filled(ix, iy) <= dem(ix, iy) + eps) then
                        cycle
                    end if
                    do d = 1, 8
                        jx = ix + DX8(d)
                        jy = iy + DY8(d)
                        if (jx < 1 .or. jx > nx .or. jy < 1 .or. jy > ny) then
                            cycle
                        end if
                        new_val = max(dem(ix, iy), filled(jx, jy) + eps)
                        if (new_val < filled(ix, iy)) then
                            filled(ix, iy) = new_val
                            changed = .true.
                        end if
                    end do
                end do
            end do
            ! Backward sweep
            do iy = ny, 1, -1
                do ix = nx, 1, -1
                    if (filled(ix, iy) <= dem(ix, iy) + eps) then
                        cycle
                    end if
                    do d = 1, 8
                        jx = ix + DX8(d)
                        jy = iy + DY8(d)
                        if (jx < 1 .or. jx > nx .or. jy < 1 .or. jy > ny) then
                            cycle
                        end if
                        new_val = max(dem(ix, iy), filled(jx, jy) + eps)
                        if (new_val < filled(ix, iy)) then
                            filled(ix, iy) = new_val
                            changed = .true.
                        end if
                    end do
                end do
            end do
            if (.not. changed) then
                exit
            end if
        end do
    end subroutine pitfill

    ! ── D8 flow direction + topological accumulation ──────────────────────────

    subroutine flow_accum_2d(nx, ny, dem, dir8, acc)
        integer, intent(in)  :: nx, ny
        real(dp), intent(in)  :: dem(nx, ny)
        integer, intent(out) :: dir8(nx, ny)   ! D8 direction (1–8) or 0 at outlets
        real(dp), intent(out) :: acc(nx, ny)

        integer, allocatable :: idx(:)
        real(dp), allocatable :: flat_dem(:)
        real(dp) :: slp, max_slp
        integer  :: ix, iy, jx, jy, d, n, i

        ! Steepest-descent D8 direction
        dir8 = 0
        do iy = 1, ny
            do ix = 1, nx
                max_slp = 0.0_dp
                do d = 1, 8
                    jx = ix + DX8(d)
                    jy = iy + DY8(d)
                    if (jx < 1 .or. jx > nx .or. jy < 1 .or. jy > ny) then
                        cycle
                    end if
                    slp = (dem(ix, iy) - dem(jx, jy))/DIST8(d)
                    if (slp > max_slp) then
                        max_slp = slp
                        dir8(ix, iy) = d
                    end if
                end do
            end do
        end do

        ! Sort cells high→low by elevation
        n = nx*ny
        allocate (idx(n), flat_dem(n))
        flat_dem = reshape(dem, [n])
        do i = 1, n
            idx(i) = i
        end do
        call qsort(idx, flat_dem, 1, n)
        deallocate (flat_dem)

        ! Accumulate: each cell passes its flow to its D8 downstream neighbour
        acc = 1.0_dp
        do i = 1, n
            ix = mod(idx(i) - 1, nx) + 1
            iy = (idx(i) - 1)/nx + 1
            d = dir8(ix, iy)
            if (d == 0) then
                cycle
            end if
            jx = ix + DX8(d)
            jy = iy + DY8(d)
            if (jx < 1 .or. jx > nx .or. jy < 1 .or. jy > ny) then
                cycle
            end if
            acc(jx, jy) = acc(jx, jy) + acc(ix, iy)
        end do
        deallocate (idx)
    end subroutine flow_accum_2d

    ! ── 2-D capsule-segment renderer ──────────────────────────────────────────

    subroutine render_seg2(ax, ay, bx, by, W, nx, ny, grid)
        real(dp), intent(in)    :: ax, ay, bx, by, W
        integer, intent(in)    :: nx, ny
        real(sp), intent(inout) :: grid(nx, ny)

        real(dp) :: ex, ey, len2, px, py, r, t
        integer  :: jx, jy, jx0, jx1, jy0, jy1

        r = max(0.5_dp, W/2.0_dp)
        ex = bx - ax
        ey = by - ay
        len2 = ex*ex + ey*ey
        jx0 = max(1, int(min(ax, bx) - r))
        jx1 = min(nx, int(max(ax, bx) + r) + 1)
        jy0 = max(1, int(min(ay, by) - r))
        jy1 = min(ny, int(max(ay, by) + r) + 1)
        do jy = jy0, jy1
            py = real(jy, dp)
            do jx = jx0, jx1
                px = real(jx, dp)
                t = 0.0_dp
                if (len2 > 1.0e-20_dp) then
                    t = max(0.0_dp, min(1.0_dp, ((px - ax)*ex + (py - ay)*ey)/len2))
                end if
                if (hypot(px - ax - t*ex, py - ay - t*ey) <= r) then
                    grid(jx, jy) = 1.0_sp
                end if
            end do
        end do
    end subroutine render_seg2

    ! ── 3-D incised capsule-segment renderer (U-shaped cross-section) ─────────

    subroutine render_seg3(ax, ay, bx, by, W, D, nx, ny, nz, vol)
        real(dp), intent(in)    :: ax, ay, bx, by, W, D
        integer, intent(in)    :: nx, ny, nz
        real(sp), intent(inout) :: vol(nx, ny, nz)

        real(dp) :: ex, ey, len2, px, py, r, t, dist_h, dz_max
        integer  :: jx, jy, jz, jx0, jx1, jy0, jy1, iz_bot

        r = max(0.5_dp, W/2.0_dp)
        ex = bx - ax
        ey = by - ay
        len2 = ex*ex + ey*ey
        jx0 = max(1, int(min(ax, bx) - r))
        jx1 = min(nx, int(max(ax, bx) + r) + 1)
        jy0 = max(1, int(min(ay, by) - r))
        jy1 = min(ny, int(max(ay, by) + r) + 1)
        do jy = jy0, jy1
            py = real(jy, dp)
            do jx = jx0, jx1
                px = real(jx, dp)
                t = 0.0_dp
                if (len2 > 1.0e-20_dp) then
                    t = max(0.0_dp, min(1.0_dp, ((px - ax)*ex + (py - ay)*ey)/len2))
                end if
                dist_h = hypot(px - ax - t*ex, py - ay - t*ey)
                if (dist_h > r) then
                    cycle
                end if
                dz_max = D*sqrt(max(0.0_dp, 1.0_dp - dist_h/r))
                iz_bot = max(1, nz - max(1, nint(dz_max)) + 1)
                do jz = iz_bot, nz
                    vol(jx, jy, jz) = 1.0_sp
                end do
            end do
        end do
    end subroutine render_seg3

    ! ── smooth a path in-place (weighted 3-point moving average) ─────────────

    subroutine smooth_path(px, py, n, n_pass)
        integer, intent(in)    :: n, n_pass
        real(dp), intent(inout) :: px(n), py(n)
        real(dp), allocatable   :: tx(:), ty(:)
        integer :: i, p
        allocate (tx(n), ty(n))
        do p = 1, n_pass
            tx = px
            ty = py
            do i = 2, n - 1
                px(i) = 0.25_dp*tx(i - 1) + 0.50_dp*tx(i) + 0.25_dp*tx(i + 1)
                py(i) = 0.25_dp*ty(i - 1) + 0.50_dp*ty(i) + 0.25_dp*ty(i + 1)
            end do
        end do
        deallocate (tx, ty)
    end subroutine smooth_path

    ! ── trace reaches (break-point to break-point), smooth, render (2-D) ────────
    !
    ! Decomposes the channel tree into "reaches": segments between headwaters,
    ! confluences (n_in >= 2), and the outlet.  Each reach is smoothed with its
    ! two endpoints FIXED at exact grid positions so that every junction is a
    ! pixel-accurate shared point — no gaps at any scale.

    subroutine render_network_2d(nx, ny, dir8, acc, threshold, acc_max, W_max, grid)
        integer, intent(in)    :: nx, ny, dir8(nx, ny)
        real(dp), intent(in)    :: acc(nx, ny), threshold, acc_max, W_max
        real(sp), intent(inout) :: grid(nx, ny)

        logical, allocatable :: is_chan(:, :)
        integer, allocatable :: n_in(:, :)
        real(dp), allocatable :: px(:), py(:), pw(:)
        integer :: ix, iy, jx, jy, kx, ky, d, n_path, i

        allocate (is_chan(nx, ny), n_in(nx, ny))
        is_chan = (acc >= threshold)

        ! Count upstream channel cells that route into each channel cell
        n_in = 0
        do iy = 1, ny
            do ix = 1, nx
                if (.not. is_chan(ix, iy)) then
                    cycle
                end if
                d = dir8(ix, iy)
                if (d == 0) then
                    cycle
                end if
                jx = ix + DX8(d)
                jy = iy + DY8(d)
                if (jx < 1 .or. jx > nx .or. jy < 1 .or. jy > ny) then
                    cycle
                end if
                if (is_chan(jx, jy)) then
                    n_in(jx, jy) = n_in(jx, jy) + 1
                end if
            end do
        end do

        allocate (px(nx*ny), py(nx*ny), pw(nx*ny))

        do iy = 1, ny
            do ix = 1, nx
                if (.not. is_chan(ix, iy)) then
                    cycle
                end if
                if (n_in(ix, iy) == 1) then
                    cycle   ! interior — start of reach handled upstream
                end if

                ! Break point: headwater (n_in==0) or confluence (n_in>=2)
                n_path = 0
                jx = ix
                jy = iy
                do
                    n_path = n_path + 1
                    px(n_path) = real(jx, dp)
                    py(n_path) = real(jy, dp)
                    pw(n_path) = max(1.0_dp, W_max*sqrt(acc(jx, jy)/acc_max))
                    d = dir8(jx, jy)
                    if (d == 0) then
                        exit
                    end if
                    kx = jx + DX8(d)
                    ky = jy + DY8(d)
                    if (kx < 1 .or. kx > nx .or. ky < 1 .or. ky > ny) then
                        exit
                    end if
                    if (.not. is_chan(kx, ky)) then
                        exit
                    end if
                    jx = kx
                    jy = ky
                    ! Stop at next confluence (add it as fixed endpoint, then exit)
                    if (n_in(jx, jy) >= 2) then
                        n_path = n_path + 1
                        px(n_path) = real(jx, dp)
                        py(n_path) = real(jy, dp)
                        pw(n_path) = max(1.0_dp, W_max*sqrt(acc(jx, jy)/acc_max))
                        exit
                    end if
                end do

                if (n_path >= 2) then
                    call smooth_path(px, py, n_path, 3)
                    do i = 1, n_path - 1
                        call render_seg2(px(i), py(i), px(i + 1), py(i + 1), pw(i), nx, ny, grid)
                    end do
                end if
            end do
        end do
        deallocate (is_chan, n_in, px, py, pw)
    end subroutine render_network_2d

    ! ── trace reaches, smooth, render (3-D) ───────────────────────────────────

    subroutine render_network_3d(nx, ny, nz, dir8, acc, threshold, acc_max, W_max, D_max, vol)
        integer, intent(in)    :: nx, ny, nz, dir8(nx, ny)
        real(dp), intent(in)    :: acc(nx, ny), threshold, acc_max, W_max, D_max
        real(sp), intent(inout) :: vol(nx, ny, nz)

        logical, allocatable :: is_chan(:, :)
        integer, allocatable :: n_in(:, :)
        real(dp), allocatable :: px(:), py(:), pw(:), pd(:)
        integer :: ix, iy, jx, jy, kx, ky, d, n_path, i
        real(dp) :: ratio

        allocate (is_chan(nx, ny), n_in(nx, ny))
        is_chan = (acc >= threshold)

        n_in = 0
        do iy = 1, ny
            do ix = 1, nx
                if (.not. is_chan(ix, iy)) then
                    cycle
                end if
                d = dir8(ix, iy)
                if (d == 0) then
                    cycle
                end if
                jx = ix + DX8(d)
                jy = iy + DY8(d)
                if (jx < 1 .or. jx > nx .or. jy < 1 .or. jy > ny) then
                    cycle
                end if
                if (is_chan(jx, jy)) then
                    n_in(jx, jy) = n_in(jx, jy) + 1
                end if
            end do
        end do

        allocate (px(nx*ny), py(nx*ny), pw(nx*ny), pd(nx*ny))

        do iy = 1, ny
            do ix = 1, nx
                if (.not. is_chan(ix, iy)) then
                    cycle
                end if
                if (n_in(ix, iy) == 1) then
                    cycle
                end if

                n_path = 0
                jx = ix
                jy = iy
                do
                    n_path = n_path + 1
                    px(n_path) = real(jx, dp)
                    py(n_path) = real(jy, dp)
                    ratio = sqrt(acc(jx, jy)/acc_max)
                    pw(n_path) = max(1.0_dp, W_max*ratio)
                    pd(n_path) = max(1.0_dp, D_max*ratio)
                    d = dir8(jx, jy)
                    if (d == 0) then
                        exit
                    end if
                    kx = jx + DX8(d)
                    ky = jy + DY8(d)
                    if (kx < 1 .or. kx > nx .or. ky < 1 .or. ky > ny) then
                        exit
                    end if
                    if (.not. is_chan(kx, ky)) then
                        exit
                    end if
                    jx = kx
                    jy = ky
                    if (n_in(jx, jy) >= 2) then
                        n_path = n_path + 1
                        px(n_path) = real(jx, dp)
                        py(n_path) = real(jy, dp)
                        ratio = sqrt(acc(jx, jy)/acc_max)
                        pw(n_path) = max(1.0_dp, W_max*ratio)
                        pd(n_path) = max(1.0_dp, D_max*ratio)
                        exit
                    end if
                end do

                if (n_path >= 2) then
                    call smooth_path(px, py, n_path, 3)
                    do i = 1, n_path - 1
                        call render_seg3(px(i), py(i), px(i + 1), py(i + 1), pw(i), pd(i), nx, ny, nz, vol)
                    end do
                end if
            end do
        end do
        deallocate (is_chan, n_in, px, py, pw, pd)
    end subroutine render_network_3d

    ! ── accumulation threshold via geometric binary search ────────────────────
    ! Finds the smallest threshold such that at least n_target cells have acc >= threshold.
    ! Geometric bisection is appropriate for the power-law accumulation distribution.

    subroutine acc_threshold(acc, nx, ny, channel_frac, threshold, acc_max_out)
        real(dp), intent(in)  :: acc(nx, ny)
        integer, intent(in)  :: nx, ny
        real(dp), intent(in)  :: channel_frac
        real(dp), intent(out) :: threshold, acc_max_out

        real(dp) :: lo_t, hi_t, mid_t
        integer  :: n_target, i

        n_target = max(1, nint(channel_frac*real(nx*ny, dp)))
        acc_max_out = maxval(acc)
        lo_t = 1.0_dp
        hi_t = acc_max_out*2.0_dp + 1.0_dp

        do i = 1, 128
            mid_t = (lo_t + hi_t)*0.5_dp
            if (count(acc >= mid_t) >= n_target) then
                lo_t = mid_t    ! threshold still low enough; try higher
            else
                hi_t = mid_t    ! threshold too high; try lower
            end if
            if ((hi_t - lo_t) < 0.5_dp) then
                exit
            end if
        end do
        threshold = lo_t
    end subroutine acc_threshold

    ! ── public: 2-D drainage network ─────────────────────────────────────────

    !>   W_max        : maximum channel width in grid cells (trunk)
    !>   channel_frac : fraction of cells to render as channels (e.g. 0.05)
    !>   hurst        : fBm Hurst exponent — controls sinuosity
    !>                  0.5=rough/sinuous, 0.85=smooth/straight
    !>   seed         : random seed (different values → different patterns)
    !>   tilt         : regional slope strength (optional, default 3.0)
    !>                  high (3–5) → 1 dominant trunk
    !>                  moderate (1–2) → 2–3 competing trunks
    subroutine drainage_slice(nx, ny, grid, W_max, channel_frac, hurst, seed, tilt)
        integer, intent(in)           :: nx, ny, seed
        real(sp), allocatable, intent(out) :: grid(:, :)
        real(dp), intent(in)           :: W_max, channel_frac, hurst
        real(dp), intent(in), optional :: tilt

        real(dp), allocatable :: dem(:, :), filled(:, :), acc(:, :)
        integer, allocatable :: dir8(:, :)
        real(dp) :: threshold, acc_max, tilt_val
        integer  :: s

        tilt_val = 3.0_dp
        if (present(tilt)) then
            tilt_val = tilt
        end if
        s = mod(abs(seed), 2147483647)

        allocate (dem(nx, ny), filled(nx, ny), acc(nx, ny), dir8(nx, ny))
        call gen_terrain(nx, ny, hurst, 7, tilt_val, s, dem)
        call pitfill(nx, ny, dem, filled)
        call flow_accum_2d(nx, ny, filled, dir8, acc)
        call acc_threshold(acc, nx, ny, channel_frac, threshold, acc_max)
        deallocate (dem, filled)

        allocate (grid(nx, ny))
        grid = 0.0_sp
        call render_network_2d(nx, ny, dir8, acc, threshold, acc_max, W_max, grid)
        deallocate (acc, dir8)
    end subroutine drainage_slice

    ! ── public: 3-D drainage volume ───────────────────────────────────────────

    !>   D_max : maximum channel depth in layers (incised from iz=nz downward)
    !>   tilt  : optional regional slope (see drainage_slice for guidance)
    subroutine drainage_volume(nx, ny, nz, vol, W_max, D_max, channel_frac, hurst, seed, tilt, drain_dir_in)
        integer, intent(in)           :: nx, ny, nz, seed
        real(sp), allocatable, intent(out) :: vol(:, :, :)
        real(dp), intent(in)           :: W_max, D_max, channel_frac, hurst
        real(dp), intent(in), optional :: tilt
        integer, intent(in), optional  :: drain_dir_in

        real(dp), allocatable :: dem(:, :), filled(:, :), acc(:, :)
        integer, allocatable :: dir8(:, :)
        real(dp) :: threshold, acc_max, tilt_val
        integer  :: s, drain_dir

        tilt_val = 3.0_dp
        if (present(tilt)) then
            tilt_val = tilt
        end if
        s = mod(abs(seed), 2147483647)
        if (present(drain_dir_in)) then
            drain_dir = drain_dir_in
        else
            drain_dir = mod(s, 4) + 1
        end if

        allocate (dem(nx, ny), filled(nx, ny), acc(nx, ny), dir8(nx, ny))
        call gen_terrain(nx, ny, hurst, 7, tilt_val, s, dem, drain_dir=drain_dir)
        call pitfill(nx, ny, dem, filled)
        call flow_accum_2d(nx, ny, filled, dir8, acc)
        call acc_threshold(acc, nx, ny, channel_frac, threshold, acc_max)
        deallocate (dem, filled)

        allocate (vol(nx, ny, nz))
        vol = 0.0_sp
        call render_network_3d(nx, ny, nz, dir8, acc, threshold, acc_max, W_max, D_max, vol)
        deallocate (acc, dir8)
    end subroutine drainage_volume

    ! ── canyon rendering: distance-to-channel → erosional depth fill ─────────────
    !
    ! For every XY cell, finds the nearest channel segment and computes a
    ! normalised depth fraction:
    !
    !   d_norm = dist / r_seg          (r_seg = half-width ∝ √acc)
    !   df     = (1 − d_norm)^p        (p = canyon_exp, controls cross-section shape)
    !            p = 1 → linear V,  p < 1 → U-shape,  p > 1 → narrow V
    !
    ! Then fills vol(iix,iiy, iz_min:nz) = 1, where
    !   iz_min = nz − round((nz−1) × df)
    ! giving a body that is DEEPEST at the thalweg and tapers to zero at the rim.

    subroutine render_canyon_3d(nx, ny, nz, dir8, acc, threshold, acc_max, &
            W_max, canyon_exp, vol, wall_exp)
        integer, intent(in)    :: nx, ny, nz, dir8(nx, ny)
        real(dp), intent(in)    :: acc(nx, ny), threshold, acc_max, W_max, canyon_exp
        real(dp), intent(in), optional :: wall_exp
        real(sp), intent(inout) :: vol(nx, ny, nz)

        real(dp), allocatable :: depth_norm(:, :)
        logical, allocatable :: is_chan(:, :)
        integer, allocatable :: n_in(:, :)
        real(dp), allocatable :: px(:), py(:), pw(:)

        real(dp) :: ex, ey, len2, ptx, pty, t_proj, dist_h, r_seg, df
        integer  :: ix, iy, jx, jy, kx, ky, d, n_path, i
        integer  :: jx0, jx1, jy0, jy1, iix, iiy, iz_min

        allocate (depth_norm(nx, ny))
        depth_norm = 0.0_dp

        ! ── reach decomposition (same logic as render_network_2d) ─────────────────
        allocate (is_chan(nx, ny), n_in(nx, ny))
        is_chan = (acc >= threshold)
        n_in = 0
        do iy = 1, ny
            do ix = 1, nx
                if (.not. is_chan(ix, iy)) then
                    cycle
                end if
                d = dir8(ix, iy)
                if (d == 0) then
                    cycle
                end if
                jx = ix + DX8(d)
                jy = iy + DY8(d)
                if (jx < 1 .or. jx > nx .or. jy < 1 .or. jy > ny) then
                    cycle
                end if
                if (is_chan(jx, jy)) then
                    n_in(jx, jy) = n_in(jx, jy) + 1
                end if
            end do
        end do

        allocate (px(nx*ny), py(nx*ny), pw(nx*ny))

        do iy = 1, ny
            do ix = 1, nx
                if (.not. is_chan(ix, iy)) then
                    cycle
                end if
                if (n_in(ix, iy) == 1) then
                    cycle
                end if

                n_path = 0
                jx = ix
                jy = iy
                do
                    n_path = n_path + 1
                    px(n_path) = real(jx, dp)
                    py(n_path) = real(jy, dp)
                    pw(n_path) = max(2.0_dp, W_max*sqrt(acc(jx, jy)/acc_max))
                    d = dir8(jx, jy)
                    if (d == 0) then
                        exit
                    end if
                    kx = jx + DX8(d)
                    ky = jy + DY8(d)
                    if (kx < 1 .or. kx > nx .or. ky < 1 .or. ky > ny) then
                        exit
                    end if
                    if (.not. is_chan(kx, ky)) then
                        exit
                    end if
                    jx = kx
                    jy = ky
                    if (n_in(jx, jy) >= 2) then
                        n_path = n_path + 1
                        px(n_path) = real(jx, dp)
                        py(n_path) = real(jy, dp)
                        pw(n_path) = max(2.0_dp, W_max*sqrt(acc(jx, jy)/acc_max))
                        exit
                    end if
                end do

                if (n_path < 2) then
                    cycle
                end if
                call smooth_path(px, py, n_path, 3)

                do i = 1, n_path - 1
                    r_seg = 0.5_dp*max(pw(i), pw(i + 1))
                    ex = px(i + 1) - px(i)
                    ey = py(i + 1) - py(i)
                    len2 = ex*ex + ey*ey
                    jx0 = max(1, int(min(px(i), px(i + 1)) - r_seg))
                    jx1 = min(nx, int(max(px(i), px(i + 1)) + r_seg) + 1)
                    jy0 = max(1, int(min(py(i), py(i + 1)) - r_seg))
                    jy1 = min(ny, int(max(py(i), py(i + 1)) + r_seg) + 1)
                    do iiy = jy0, jy1
                        pty = real(iiy, dp)
                        do iix = jx0, jx1
                            ptx = real(iix, dp)
                            t_proj = 0.0_dp
                            if (len2 > 1.0e-20_dp) then
                                t_proj = max(0.0_dp, min(1.0_dp, &
                                    ((ptx - px(i))*ex + (pty - py(i))*ey)/len2))
                            end if
                            dist_h = hypot(ptx - px(i) - t_proj*ex, pty - py(i) - t_proj*ey)
                            if (dist_h >= r_seg) then
                                cycle
                            end if
                            df = (1.0_dp - dist_h/r_seg)**canyon_exp
                            if (df > depth_norm(iix, iiy)) then
                                depth_norm(iix, iiy) = df
                            end if
                        end do
                    end do
                end do
            end do
        end do
        deallocate (is_chan, n_in, px, py, pw)

        ! ── erosional fill: each cell open from iz_min up to the seabed (iz=nz) ───
        ! wall_exp remaps depth_norm before the iz_min conversion:
        !   wall_exp = 1.0 : linear walls (default)
        !   wall_exp < 1.0 : concave – canyon flares open quickly from the thalweg
        !   wall_exp > 1.0 : convex – narrow stem, rapid widening near surface
        block
            real(dp) :: wexp, dn
            wexp = 1.0_dp
            if (present(wall_exp)) then
                wexp = max(0.1_dp, wall_exp)
            end if
            do iiy = 1, ny
                do iix = 1, nx
                    dn = depth_norm(iix, iiy)
                    if (dn <= 0.0_dp) then
                        cycle
                    end if
                    dn = dn**wexp
                    iz_min = max(1, nz - nint(real(nz - 1, dp)*dn))
                    vol(iix, iiy, iz_min:nz) = 1.0_sp
                end do
            end do
        end block
        deallocate (depth_norm)

    end subroutine render_canyon_3d

    ! ── public: 3-D submarine canyon volume ──────────────────────────────────────
    !
    !> Generate a 3-D submarine canyon body from a dendritic drainage network.
    !>
    !> The canyon depth at each XY cell is determined by how close that cell is to
    !> the nearest channel segment, scaled by local channel width (∝ √accumulation):
    !>   depth_frac = (1 − dist/r_local)^canyon_exp,  r_local ∝ √acc / acc_max
    !>
    !> The volume is filled erosionally: vol(ix,iy, iz_min:nz) = 1, so the body
    !> is deepest at the thalweg and tapers to zero at the canyon rim.
    !>
    !>   W_max        : canyon full-width at surface for the main trunk (grid cells)
    !>   channel_frac : fraction of cells that are channels
    !>                  0.01–0.02 → main trunk only + 1–2 tributaries (best for
    !>                  canyon simulation);  0.05 → more branching
    !>   canyon_exp   : cross-section exponent  (1=V-shape, 0.5=U-shape, 2=narrow-V)
    !>   wall_exp     : vertical wall curvature (1=linear, <1=concave/flared, >1=convex/steep-stem)
    !>   tilt         : regional slope (3–5 → single dominant trunk)
    subroutine build_drainage_canyon(nx, ny, nz, vol, W_max, channel_frac, hurst, &
            seed, tilt, canyon_exp, wall_exp, drain_dir_in)
        integer, intent(in)           :: nx, ny, nz, seed
        real(sp), allocatable, intent(out) :: vol(:, :, :)
        real(dp), intent(in)           :: W_max, channel_frac, hurst
        real(dp), intent(in), optional :: tilt, canyon_exp, wall_exp
        integer, intent(in), optional  :: drain_dir_in

        real(dp), allocatable :: dem(:, :), filled(:, :), acc(:, :)
        integer, allocatable :: dir8(:, :)
        real(dp) :: threshold, acc_max, tilt_val, cexp, wexp
        integer  :: s, drain_dir

        tilt_val = 3.0_dp
        if (present(tilt)) then
            tilt_val = tilt
        end if
        cexp = 1.0_dp
        if (present(canyon_exp)) then
            cexp = canyon_exp
        end if
        wexp = 1.0_dp
        if (present(wall_exp)) then
            wexp = wall_exp
        end if
        s = mod(abs(seed), 2147483647)
        if (present(drain_dir_in)) then
            drain_dir = drain_dir_in
        else
            drain_dir = mod(s, 4) + 1
        end if

        allocate (dem(nx, ny), filled(nx, ny), acc(nx, ny), dir8(nx, ny))
        call gen_terrain(nx, ny, hurst, 7, tilt_val, s, dem, drain_dir=drain_dir)
        call pitfill(nx, ny, dem, filled)
        call flow_accum_2d(nx, ny, filled, dir8, acc)
        call acc_threshold(acc, nx, ny, channel_frac, threshold, acc_max)
        deallocate (dem, filled)

        allocate (vol(nx, ny, nz))
        vol = 0.0_sp
        call render_canyon_3d(nx, ny, nz, dir8, acc, threshold, acc_max, W_max, cexp, vol, wexp)
        deallocate (acc, dir8)

    end subroutine build_drainage_canyon

    ! ── OOP generate subroutines ─────────────────────────────────────────────────

    subroutine generate_drainage_channel(self)
        class(drainage_channel), intent(inout) :: self
        integer :: iix, iiy, iiz, nz_g, eff_seed, s_val, drain_dir_val
        real(dp), allocatable :: dem_bg(:, :)

        eff_seed = safe_seed(int(self%seed, 8))

        ! Drain direction: fixed by orient if given, otherwise derived from the seed
        if (self%orient >= 0) then
            drain_dir_val = mod(self%orient, 4) + 1
        else
            drain_dir_val = mod(mod(abs(eff_seed), 2147483647), 4) + 1
        end if

        call drainage_volume(self%nx, self%ny, self%nz, self%vol, &
            self%W_max, self%D_max, self%channel_frac, self%hurst, eff_seed, &
            tilt=self%tilt, drain_dir_in=drain_dir_val)

        nz_g = size(self%vol, 3)
        if (allocated(self%depth_map)) then
            deallocate (self%depth_map)
        end if
        allocate (self%depth_map(self%nx, self%ny))
        self%depth_map = 0.0_sp
        do iiy = 1, self%ny
            do iix = 1, self%nx
                do iiz = 1, nz_g
                    if (self%vol(iix, iiy, iiz) > 0.5_sp) then
                        self%depth_map(iix, iiy) = real(nz_g - iiz + 1, sp)/real(nz_g, sp)
                        exit
                    end if
                end do
            end do
        end do
        self%depth_map = self%depth_map*(nz_g - 1.0_sp)
        self%depth_map = transpose(self%depth_map)
        if (.not. self%yn_depth_only) then
            self%vol = flip(permute(self%vol, 321), [1])
        end if

        ! DEM
        allocate (dem_bg(self%nx, self%ny))
        dem_bg = 0.0
        if (self%terrain_bg > 0.0_dp) then
            s_val = mod(abs(eff_seed), 2147483647)
            call gen_terrain(self%nx, self%ny, self%hurst, 7, self%tilt, s_val, dem_bg, &
                drain_dir=drain_dir_val)
            dem_bg = real(transpose(dem_bg)*self%terrain_bg, sp)
            ! Background cells get terrain value [0, terrain_bg].
            ! Channel cells get terrain − channel_depth (negative = incised below surface).
            where (self%depth_map == 0.0_sp)
                self%depth_map = dem_bg
            elsewhere
                self%depth_map = dem_bg - self%depth_map
            end where
        else
            self%depth_map = -self%depth_map
            dem_bg = real(transpose(dem_bg), sp)
        end if
        dem_bg = nint(dem_bg)
        dem_bg = maxval(dem_bg) - dem_bg

        ! Convert to downward-positive coordinates
        self%depth_map = nint(self%depth_map)
        self%depth_map = maxval(self%depth_map) - self%depth_map

        self%nnz = int(abs(maxval(self%depth_map) - minval(self%depth_map))) + 1
        deallocate (self%vol)
        if (.not. self%yn_depth_only) then
            allocate (self%vol(self%nnz, self%ny, self%nx))
            self%vol = 1.0
            if (allocated(self%obj)) then
                deallocate (self%obj)
            end if
            allocate (self%obj(self%nnz, self%ny, self%nx))
            self%obj = 0.0
            !$omp parallel do private(iix, iiy)
            do iiy = 1, self%ny
                do iix = 1, self%nx
                    self%vol(1:self%depth_map(iiy, iix) + 1, iiy, iix) = 0.0
                    self%obj(nint(dem_bg(iiy, iix)) + 2:self%depth_map(iiy, iix) + 1, iiy, iix) = 1.0
                end do
            end do
            !$omp end parallel do
        end if

    end subroutine generate_drainage_channel

    subroutine generate_drainage_canyon(self)
        class(drainage_canyon), intent(inout) :: self
        integer :: iix, iiy, iiz, nz_g, eff_seed, s_val, drain_dir_val
        real(dp), allocatable :: dem_bg(:, :)

        eff_seed = safe_seed(int(self%seed, 8))

        ! Drain direction: fixed by orient if given, otherwise derived from the seed
        if (self%orient >= 0) then
            drain_dir_val = mod(self%orient, 4) + 1
        else
            drain_dir_val = mod(mod(abs(eff_seed), 2147483647), 4) + 1
        end if

        call build_drainage_canyon(self%nx, self%ny, self%nz, self%vol, &
            self%W_max, self%channel_frac, self%hurst, &
            eff_seed, self%tilt, self%canyon_exp, self%wall_exp, drain_dir_in=drain_dir_val)

        nz_g = size(self%vol, 3)
        if (allocated(self%depth_map)) then
            deallocate (self%depth_map)
        end if
        allocate (self%depth_map(self%nx, self%ny))
        self%depth_map = 0.0_sp
        do iiy = 1, self%ny
            do iix = 1, self%nx
                do iiz = 1, nz_g
                    if (self%vol(iix, iiy, iiz) > 0.5_sp) then
                        self%depth_map(iix, iiy) = real(nz_g - iiz + 1, sp)/real(nz_g, sp)
                        exit
                    end if
                end do
            end do
        end do
        self%depth_map = self%depth_map*(nz_g - 1.0_sp)
        self%depth_map = transpose(self%depth_map)
        if (.not. self%yn_depth_only) then
            self%vol = flip(permute(self%vol, 321), [1])
        end if

        ! DEM
        allocate (dem_bg(self%nx, self%ny))
        dem_bg = 0.0
        if (self%terrain_bg > 0.0_dp) then
            s_val = mod(abs(eff_seed), 2147483647)
            call gen_terrain(self%nx, self%ny, self%hurst, 7, self%tilt, s_val, dem_bg, &
                drain_dir=drain_dir_val)
            dem_bg = real(transpose(dem_bg)*self%terrain_bg, sp)
            where (self%depth_map == 0.0_sp)
                self%depth_map = dem_bg
            elsewhere
                self%depth_map = dem_bg - self%depth_map
            end where
        else
            self%depth_map = -self%depth_map
            dem_bg = real(transpose(dem_bg), sp)
        end if
        dem_bg = nint(dem_bg)
        dem_bg = maxval(dem_bg) - dem_bg

        ! Convert to downward-positive coordinates
        self%depth_map = nint(self%depth_map)
        self%depth_map = maxval(self%depth_map) - self%depth_map

        self%nnz = int(abs(maxval(self%depth_map) - minval(self%depth_map))) + 1
        deallocate (self%vol)
        if (.not. self%yn_depth_only) then
            allocate (self%vol(self%nnz, self%ny, self%nx))
            self%vol = 1.0
            if (allocated(self%obj)) then
                deallocate (self%obj)
            end if
            allocate (self%obj(self%nnz, self%ny, self%nx))
            self%obj = 0.0
            !$omp parallel do private(iix, iiy)
            do iiy = 1, self%ny
                do iix = 1, self%nx
                    self%vol(1:self%depth_map(iiy, iix) + 1, iiy, iix) = 0.0
                    self%obj(nint(dem_bg(iiy, iix)) + 2:self%depth_map(iiy, iix) + 1, iiy, iix) = 1.0
                end do
            end do
            !$omp end parallel do
        end if

    end subroutine generate_drainage_canyon

end module
