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
! Meandering river / submarine canyon generator
!
! Algorithm — Howard-Knutson (1983) planform migration
! -----------------------------------------------------
! 1. Initialise: a sine-generated centerline of ~n_bends half-wavelengths is
!    laid out; nodes are spaced uniformly at ds = W/2.
!
! 2. Curvature: the signed curvature C(s) at each node is estimated from the
!    cross product of consecutive unit-tangent vectors, divided by ds.
!
! 3. Transverse velocity u*(s): computed by upstream integration of C(s),
!    following Ikeda, Parker & Sawai (1981).  The integrand decays with a
!    spatial lag scale k, and is modulated by the secondary-flow coefficient
!    omega and the ratio gamma = width / depth.
!
! 4. Migration: each node is displaced normal to the centreline by
!      Δn = kl · u* · dt
!    (Howard & Knutson 1983).  Nodes near the ends of the domain are tapered
!    to zero displacement to suppress boundary reflections.
!
! 5. Cutoffs: after every migration step, pairs of non-adjacent nodes closer
!    than crdist = W are candidates for oxbow cutoffs.  When a cutoff is
!    detected the intervening loop is removed and the two endpoints are joined.
!
! 6. Resampling: after migration (and after each cutoff) the centerline is
!    resampled to restore uniform node spacing at ds.
!
! 7. Steps 2–6 repeat for n_iter time steps.
!
! 3-D channel volume (simulate_meander_3d)
! ----------------------------------------
! Boundary padding: the first/last npad = max(10, round(0.05*N0)) nodes are
! held fixed (zero displacement) to suppress edge artefacts.  Before each
! snapshot is projected, those pinned nodes are trimmed so only freely-migrating
! nodes reach the grid edges.
!
! OBB projection: the trimmed centreline is mapped to the nx×ny output grid
! via an Oriented Bounding Box (endpoint unit vector ê).  The u-axis (along ê)
! maps to grid columns 1..nx; the v-axis is padded by W/2 on each side and maps
! to rows 1..ny, keeping the channel body away from the non-start/end edges.
!
! Cross-section shape (W_shape): at depth layer iz (t = (iz-1)/(nz-1)),
! half-width W_iz = W_bot + (W - W_bot)*f(t) with four profile options:
!   0  power law  f = t^W_exp  (W_exp<1 → U, =1 → V, >1 → flared)
!   1  smoothstep f = 3t²-2t³  (flat floor, steep walls, gentle shoulders)
!   2  Gaussian CDF  f = normalised erf, steepness sigma = W_exp
!   3  Gaussian bell (upside-down bell), sigma = W_exp
! Every snapshot accumulates into a single volume; deepest iz per (x,y) column
! is tracked and the body fills upward to iz = nz (composite channel sand body).
!
! 3-D canyon volume (simulate_meander_canyon_3d)
! -----------------------------------------------
! Snapshots are trimmed of npad boundary nodes before saving, and the output
! grid is internally padded by pad = round(0.1*nx) columns on each side so
! that the canyon exits the grid edges naturally (padded columns are cropped
! from the output after rendering).
!
! z-assignment follows a V-trajectory: t=0 and t=1 → iz ≈ nz (surface);
! t=0.5 → iz = 1 (maximum incision).  n_stages independently random longitudinal
! offset profiles add stratigraphic roughness (amplitude dz_sigma*nz).
! W_eff widens linearly from W (thalweg) to W_canyon (surface), giving the
! characteristic erosional flare of a submarine canyon.
!
! Public OOP API
! --------------
!   type(meandering_channel) :: chan
!   call chan%generate          ! fills chan%depth_map(nx,ny) and chan%vol(nx,ny,nz)
!
!   type(meandering_canyon) :: canyon
!   call canyon%generate        ! fills canyon%depth_map(nx,ny) and canyon%vol(nx,ny,nz)
!
!   depth_map: 0 = background; value = depth of channel floor in iz-units (0..nz-1)
!   vol:       1.0 = channel sand / erosional body, 0.0 = host rock
!==============================================================================

module geological_model_meander

    use libflit
    use geological_model_utility

    use, intrinsic :: iso_fortran_env, only: dp => real64, sp => real32
    implicit none
    private

    !   real(dp), parameter :: PI = acos(-1.0_dp)

    !    public :: simulate
    !    public :: simulate_meander_2d
    !    public :: simulate_meander_3d
    !    public :: simulate_meander_canyon_3d
    !    public :: project_centerline_to_grid
    !    public :: project_channel_to_grid
    !> Meandering river channel sand body.
    !> call generate() fills depth_map(nx,ny) and vol(nx,ny,nz).
    type meandering_channel
        integer  :: nx = 256, ny = 256  ! output grid dimensions (map view)
        integer  :: nz = 32             ! vertical layers (stratigraphic thickness)
        real(dp) :: length = 5000.0_dp  ! initial channel length (same units as W)
        real(dp) :: W = 200.0_dp        ! bankfull channel width
        real(dp) :: D = 10.0_dp         ! bankfull channel depth
        real(dp) :: Cf = 0.011_dp       ! friction coefficient (controls u* magnitude)
        real(dp) :: k = 1.0_dp          ! upstream-weighting decay constant for u*
        real(dp) :: kl = 60.0_dp        ! lateral migration rate coefficient
        real(dp) :: omega = -1.0_dp     ! secondary-flow weight (negative → outward bank migration)
        real(dp) :: gamma = 2.5_dp      ! width-to-depth ratio W/D used in u* integral
        integer  :: n_bends = 20        ! number of initial half-wavelength bends
        real(dp) :: dt = 0.05_dp        ! dimensionless time step (fraction of W/kl)
        integer  :: n_iter = 3000       ! total migration steps
        integer  :: seed = -1           ! random seed; -1 = clock-based
        real(dp) :: W_bot = 1.0_dp      ! channel-floor width; ≤ 0 → use W/10
        real(dp) :: W_exp = 0.5_dp      ! width-vs-depth exponent (power law); sigma for Gaussian CDF
        integer  :: W_shape = 0         ! cross-section profile: 0=power law, 1=smoothstep, 2=Gaussian CDF, 3=Gaussian bell
        integer  :: nnz = 0
        real(dp) :: terrain_bg = 0.0_dp ! if > 0, terrain (normalised to [0, terrain_bg]) fills non-channel cells in depth_map
        integer  :: orient = -1         ! output orientation: -1 = derived from seed; 0-3 = fixed (0 = flow along x)
        logical  :: yn_depth_only = .false. ! if .true., skip rebuilding 3D vol/obj; only depth_map is produced
        real(dp) :: W_render = 0.0_dp   ! rendered channel width; <= 0 -> use W; allows a wide channel to be
        ! rendered while the migration dynamics keep a healthy node count
        real(sp), allocatable, dimension(:, :) :: depth_map
        real(sp), allocatable, dimension(:, :, :) :: vol
        real(sp), allocatable, dimension(:, :, :) :: obj
    contains
        procedure :: generate => generate_channel
    end type meandering_channel

    !> Meandering submarine canyon (stacked incision history).
    !> call generate() fills depth_map(nx,ny) and vol(nx,ny,nz).
    type meandering_canyon
        integer  :: nx = 256, ny = 256  ! output grid dimensions (map view)
        integer  :: nz = 0              ! z-levels; 0 = auto-set to n_snapshots
        real(dp) :: length = 5000.0_dp  ! initial channel length
        real(dp) :: W = 200.0_dp        ! channel width at maximum incision (thalweg)
        real(dp) :: D = 10.0_dp         ! channel depth
        real(dp) :: Cf = 0.011_dp       ! friction coefficient
        real(dp) :: k = 1.0_dp          ! upstream-weighting decay for u*
        real(dp) :: kl = 60.0_dp        ! lateral migration rate coefficient
        real(dp) :: omega = -1.0_dp     ! secondary-flow weight
        real(dp) :: gamma = 2.5_dp      ! width-to-depth ratio used in u* integral
        integer  :: n_bends = 20        ! initial number of half-wavelength bends
        real(dp) :: dt = 0.05_dp        ! dimensionless time step
        integer  :: n_iter = 3000       ! total migration steps
        integer  :: seed = -1           ! random seed; -1 = clock-based
        integer  :: save_every = 50     ! migration steps between depth-plane snapshots
        real(dp) :: W_canyon = 1.0_dp   ! canyon aperture at surface; ≤ W → same as W
        integer  :: n_stages = 10       ! number of z-trajectory stages (controls sinuosity of incision path)
        real(dp) :: dz_sigma = 0.05_dp  ! Gaussian jitter on z-position as fraction of nz (stratigraphic roughness)
        integer  :: nnz = 0
        real(dp) :: terrain_bg = 0.0_dp ! if > 0, terrain (normalised to [0, terrain_bg]) fills non-channel cells in depth_map
        integer  :: orient = -1         ! output orientation: -1 = derived from seed; 0-3 = fixed (0 = flow along x)
        logical  :: yn_depth_only = .false. ! if .true., skip rebuilding 3D vol/obj; only depth_map is produced
        real(sp), allocatable, dimension(:, :) :: depth_map
        real(sp), allocatable, dimension(:, :, :) :: vol
        real(sp), allocatable, dimension(:, :, :) :: obj
    contains
        procedure :: generate => generate_canyon
    end type meandering_canyon

    public :: meandering_channel
    public :: meandering_canyon

contains

    ! ── numerical utilities ─────────────────────────────────────────────────────

    !> Central-difference gradient with unit spacing (same as np.gradient for 1D).
    !> da(1) = a(2)-a(1),  da(i) = 0.5*(a(i+1)-a(i-1)),  da(n) = a(n)-a(n-1)
    subroutine grad1d(a, da)
        real(dp), intent(in)  :: a(:)
        real(dp), intent(out) :: da(:)
        integer :: i, n
        n = size(a)
        if (n < 2) then
            da(1) = 0.0_dp
            return
        end if
        da(1) = a(2) - a(1)
        do i = 2, n - 1
            da(i) = 0.5_dp*(a(i + 1) - a(i - 1))
        end do
        da(n) = a(n) - a(n - 1)
    end subroutine grad1d

    !> Causal single-pole IIR filter matching scipy.signal.lfilter([1-d],[1,-d],x).
    !> y(1) = (1-d)*x(1),   y(i) = (1-d)*x(i) + d*y(i-1)
    subroutine iir_causal(x, d, y)
        real(dp), intent(in)  :: x(:), d
        real(dp), intent(out) :: y(:)
        integer :: i
        y(1) = (1.0_dp - d)*x(1)
        do i = 2, size(x)
            y(i) = (1.0_dp - d)*x(i) + d*y(i - 1)
        end do
    end subroutine iir_causal

    ! ── arc-length resampling ────────────────────────────────────────────────────

    !> Resample (x,y) to uniform arc-length spacing deltas using cubic splines.
    !> x and y are reallocated in place.
    subroutine resample_chan(x, y, deltas)
        real(dp), allocatable, intent(inout) :: x(:), y(:)
        real(dp), intent(in) :: deltas

        real(dp), allocatable :: s(:), sd(:), xd(:), yd(:), s_new(:)
        integer :: i, n, nd, n_new

        n = size(x)
        allocate (s(n))

        ! cumulative arc length
        s(1) = 0.0_dp
        do i = 2, n
            s(i) = s(i - 1) + hypot(x(i) - x(i - 1), y(i) - y(i - 1))
        end do

        if (s(n) < 3.0_dp*deltas .or. n < 4) then
            return
        end if

        ! remove consecutive duplicate arc-length values (zero-length steps)
        allocate (sd(n), xd(n), yd(n))
        nd = 1
        sd(1) = s(1)
        xd(1) = x(1)
        yd(1) = y(1)
        do i = 2, n
            if (s(i) > sd(nd) + 1.0e-10_dp) then
                nd = nd + 1
                sd(nd) = s(i)
                xd(nd) = x(i)
                yd(nd) = y(i)
            end if
        end do

        if (nd < 4) then
            return
        end if

        ! uniform query locations
        n_new = max(4, nint(s(n)/deltas) + 1)
        s_new = linspace(0.0_dp, s(n), n_new)

        ! cubic spline interpolation via flit
        deallocate (x, y)
        x = ginterp(sd(1:nd), xd(1:nd), s_new, method='cubic')
        y = ginterp(sd(1:nd), yd(1:nd), s_new, method='cubic')

    end subroutine resample_chan

    ! ── cutoff detection ─────────────────────────────────────────────────────────

    !> Detect and remove neck cutoffs.  Scan upstream-to-downstream for the first
    !> pair of non-adjacent nodes closer than crdist; excise the loop; repeat.
    subroutine cut_off_cutoffs(x, y, crdist, deltas, n_ox)
        real(dp), allocatable, intent(inout) :: x(:), y(:)
        real(dp), intent(in)  :: crdist, deltas
        integer, intent(out) :: n_ox

        real(dp), allocatable :: xt(:), yt(:)
        integer :: n, diag_blank, i, j, i1, j1, n_new
        logical :: found

        n_ox = 0
        diag_blank = max(5, nint(2.0_dp*crdist/deltas))

        do
            n = size(x)
            if (n < 2*diag_blank + 4) exit

            ! scan for upstream-most cutoff pair (Python: argmin(minimum(i,j)) = min i)
            found = .false.
            outer: do i = 1, n - diag_blank
                do j = i + diag_blank, n
                    if (hypot(x(i) - x(j), y(i) - y(j)) < crdist) then
                        i1 = i
                        j1 = j
                        found = .true.
                        exit outer
                    end if
                end do
            end do outer

            if (.not. found) then
                exit
            end if

            n_ox = n_ox + 1
            ! excise loop [i1+1 .. j1-1]: keep x(1:i1) and x(j1:n)
            n_new = i1 + (n - j1 + 1)
            allocate (xt(n_new), yt(n_new))
            xt(1:i1) = x(1:i1)
            xt(i1 + 1:n_new) = x(j1:n)
            yt(1:i1) = y(1:i1)
            yt(i1 + 1:n_new) = y(j1:n)
            call move_alloc(xt, x)
            call move_alloc(yt, y)
            call resample_chan(x, y, deltas)
        end do

    end subroutine cut_off_cutoffs

    ! ── public simulate ──────────────────────────────────────────────────────────

    !> Simulate meandering channel migration (Howard-Knutson 1983 / meanderpy).
    !!
    !! Parameters match meander2d.py::simulate() exactly.
    !!   kl   : lateral migration coefficient (primary knob for curve size)
    !!   dt   : time step — use 1.0–2.0 for well-developed meanders
    !!
    !! On return x_out/y_out hold the final centreline node coordinates.
    subroutine simulate_meander_2d(length, W, D, Cf, k, kl, omega, gamma, &
            n_bends, dt, n_iter, crdist, seed, &
            x_out, y_out, n_oxbows_out, pad_out)

        real(dp), intent(in)  :: length, W, D, Cf, k, kl, omega, gamma, dt
        integer, intent(in)  :: n_bends, n_iter, seed
        real(dp), intent(in), optional :: crdist
        real(dp), allocatable, intent(out) :: x_out(:), y_out(:)
        integer, intent(out), optional   :: n_oxbows_out, pad_out

        ! working arrays
        real(dp), allocatable :: x(:), y(:), raw(:), rd(:)
        real(dp), allocatable :: dx(:), dy(:), ddx(:), ddy(:)
        real(dp), allocatable :: den(:), curv(:), R0(:), Omg(:), R1(:)
        real(dp), allocatable :: ds(:), tx(:), ty(:)
        real(dp), allocatable :: dpx(:), dpy(:), mag(:), sc(:)

        real(dp) :: crdist_val, deltas, alpha, decay, phi_k
        real(dp) :: kap, th, pos, slope_val, rms_val, arc, chord
        integer  :: n0, n, it, j, npad, n_ox, n_ox_total

        ! ── setup ─────────────────────────────────────────────────────────────────
        if (present(crdist)) then
            crdist_val = crdist
        else
            crdist_val = 1.5_dp*W
        end if

        deltas = W/2.0_dp
        alpha = k*2.0_dp*Cf/D
        decay = exp(-alpha*deltas)

        ! call init_rng(seed)

        ! ── OU random initial centreline ──────────────────────────────────────────
        n0 = max(12, nint(length/deltas) + 1)
        npad = max(10, nint(n0*0.05_dp))   ! ~5 % of nodes, min 10
        phi_k = exp(-deltas/(length/max(1, n_bends)*0.4_dp))

        allocate (raw(n0))
        kap = 0.0_dp
        th = 0.0_dp
        pos = 0.0_dp
        rd = random(n0, dist='normal', seed=seed)
        do j = 1, n0
            kap = phi_k*kap + sqrt(1.0_dp - phi_k**2)*rd(j)
            th = th + kap
            th = max(-1.2_dp, min(1.2_dp, th))
            pos = pos + sin(th)
            raw(j) = pos
        end do

        ! detrend and normalise RMS to 2W
        slope_val = (raw(n0) - raw(1))/real(n0 - 1, dp)
        do j = 1, n0
            raw(j) = raw(j) - slope_val*real(j - 1, dp)
        end do
        raw = raw - sum(raw)/real(n0, dp)
        rms_val = sqrt(sum(raw**2)/real(n0, dp))
        if (rms_val > 1.0e-10_dp) then
            raw = raw*(2.0_dp*W/rms_val)
        end if

        allocate (x(n0), y(n0))
        do j = 1, n0
            x(j) = real(j - 1, dp)*length/real(n0 - 1, dp)
            y(j) = raw(j)
        end do
        deallocate (raw)
        call resample_chan(x, y, deltas)

        n_ox_total = 0

        ! ── main migration loop ────────────────────────────────────────────────────
        do it = 1, n_iter
            n = size(x)
            if (n < 2*npad + 4) then
                exit
            end if

            ! curvature  κ = (x'y'' - y'x'') / (x'^2 + y'^2)^(3/2)
            allocate (dx(n), dy(n), ddx(n), ddy(n), den(n), curv(n))
            call grad1d(x, dx)
            call grad1d(y, dy)
            call grad1d(dx, ddx)
            call grad1d(dy, ddy)
            den = (dx**2 + dy**2)**1.5_dp
            where (den > 1.0e-12_dp)
                curv = (dx*ddy - dy*ddx)/den
            elsewhere
                curv = 0.0_dp
            end where

            ! migration rate  R1 = omega*R0 + gamma*Omega
            ! Omega = lfilter([1-decay],[1,-decay], R0)  (causal IIR, O(n))
            allocate (R0(n), Omg(n), R1(n))
            R0 = kl*W*curv
            call iir_causal(R0, decay, Omg)
            R1 = omega*R0 + gamma*Omg

            ! sinuosity correction: slow migration as channel becomes more sinuous
            allocate (ds(n))
            ds = hypot(dx, dy)
            arc = sum(ds)
            chord = max(hypot(x(n) - x(1), y(n) - y(1)), 1.0e-10_dp)
            R1 = R1*(chord/arc)**(2.0_dp/3.0_dp)

            ! displacement along right normal (ty, -tx)
            ds = max(ds, 1.0e-12_dp)
            allocate (tx(n), ty(n))
            tx = dx/ds
            ty = dy/ds

            allocate (dpx(n - 2*npad), dpy(n - 2*npad))
            dpx = R1(npad + 1:n - npad)*ty(npad + 1:n - npad)*dt
            dpy = -R1(npad + 1:n - npad)*tx(npad + 1:n - npad)*dt

            ! CFL clamp: max displacement = 0.5*deltas per step
            allocate (mag(n - 2*npad), sc(n - 2*npad))
            mag = hypot(dpx, dpy)
            where (mag > 0.5_dp*deltas)
                sc = 0.5_dp*deltas/(mag + 1.0e-10_dp)
            elsewhere
                sc = 1.0_dp
            end where

            x(npad + 1:n - npad) = x(npad + 1:n - npad) + dpx*sc
            y(npad + 1:n - npad) = y(npad + 1:n - npad) + dpy*sc

            deallocate (dx, dy, ddx, ddy, den, curv, R0, Omg, R1, ds, tx, ty, dpx, dpy, mag, sc)

            ! resample + cutoff detection
            call resample_chan(x, y, deltas)
            call cut_off_cutoffs(x, y, crdist_val, deltas, n_ox)
            n_ox_total = n_ox_total + n_ox

        end do

        x_out = x
        y_out = y
        if (present(n_oxbows_out)) then
            n_oxbows_out = n_ox_total
        end if
        if (present(pad_out)) then
            pad_out = npad
        end if

    end subroutine

    ! ── grid projection ──────────────────────────────────────────────────────────

    !> Rasterize centreline segments as a thin (1-cell-wide) line.
    !> Matches Python project_centerline_to_grid.
    !> grid(nx,ny): 1.0=centreline, 0.0=background.
    subroutine project_centerline_to_grid(x, y, nx, ny, grid, &
            xmin_opt, xmax_opt, ymin_opt, ymax_opt)
        real(dp), intent(in)  :: x(:), y(:)
        integer, intent(in)  :: nx, ny
        real(sp), allocatable, intent(out) :: grid(:, :)
        real(dp), intent(in), optional :: xmin_opt, xmax_opt, ymin_opt, ymax_opt

        real(dp) :: xmin, xmax, ymin, ymax, ddx, ddy, gpad
        real(dp) :: gx0, gy0, gx1, gy1
        integer  :: i, n

        n = size(x)
        if (present(xmin_opt)) then
            xmin = xmin_opt
            xmax = xmax_opt
            ymin = ymin_opt
            ymax = ymax_opt
        else
            xmin = minval(x)
            xmax = maxval(x)
            ymin = minval(y)
            ymax = maxval(y)
            gpad = max(xmax - xmin, ymax - ymin)*0.01_dp
            xmin = xmin - gpad
            xmax = xmax + gpad
            ymin = ymin - gpad
            ymax = ymax + gpad
        end if

        ddx = (xmax - xmin)/nx
        ddy = (ymax - ymin)/ny

        allocate (grid(nx, ny))
        grid = 0.0_sp

        do i = 1, n - 1
            gx0 = (x(i) - xmin)/ddx
            gy0 = (y(i) - ymin)/ddy
            gx1 = (x(i + 1) - xmin)/ddx
            gy1 = (y(i + 1) - ymin)/ddy
            call bresenham_seg(gx0, gy0, gx1, gy1, nx, ny, grid)
        end do
    end subroutine project_centerline_to_grid

    !> Walk one segment with Bresenham and mark intersected grid cells.
    subroutine bresenham_seg(gx0, gy0, gx1, gy1, nx, ny, grid)
        real(dp), intent(in)    :: gx0, gy0, gx1, gy1
        integer, intent(in)    :: nx, ny
        real(sp), intent(inout) :: grid(nx, ny)

        integer :: ix0, iy0, ix1, iy1, iix, iiy, absdx, absdy, sx, sy, berr, e2

        ix0 = nint(gx0)
        iy0 = nint(gy0)
        ix1 = nint(gx1)
        iy1 = nint(gy1)
        absdx = abs(ix1 - ix0)
        absdy = abs(iy1 - iy0)
        if (ix1 >= ix0) then
            sx = 1
        else
            sx = -1
        end if
        if (iy1 >= iy0) then
            sy = 1
        else
            sy = -1
        end if

        iix = ix0
        iiy = iy0
        berr = absdx - absdy
        do
            if (iix >= 1 .and. iix <= nx .and. iiy >= 1 .and. iiy <= ny) &
                grid(iix, iiy) = 1.0_sp
            if (iix == ix1 .and. iiy == iy1) exit
            e2 = 2*berr
            if (e2 > -absdy) then
                berr = berr - absdy
                iix = iix + sx
            end if
            if (e2 < absdx) then
                berr = berr + absdx
                iiy = iiy + sy
            end if
        end do
    end subroutine bresenham_seg

    !> Project channel (with width W) to grid using distance-to-segment.
    !> Matches Python project_channel_to_grid.
    !> grid(nx,ny): 1.0=channel, 0.0=background.
    subroutine project_channel_to_grid(x, y, W, nx, ny, grid, &
            xmin_opt, xmax_opt, ymin_opt, ymax_opt)
        real(dp), intent(in)  :: x(:), y(:), W
        integer, intent(in)  :: nx, ny
        real(sp), allocatable, intent(out) :: grid(:, :)
        real(dp), intent(in), optional :: xmin_opt, xmax_opt, ymin_opt, ymax_opt

        real(dp) :: xmin, xmax, ymin, ymax, ddx, ddy, gpad
        real(dp) :: xc, yc, dmin, d
        integer  :: i, n, iix, iiy

        n = size(x)
        if (present(xmin_opt)) then
            xmin = xmin_opt
            xmax = xmax_opt
            ymin = ymin_opt
            ymax = ymax_opt
        else
            xmin = minval(x)
            xmax = maxval(x)
            ymin = minval(y)
            ymax = maxval(y)
            gpad = max(xmax - xmin, ymax - ymin)*0.01_dp
            xmin = xmin - gpad
            xmax = xmax + gpad
            ymin = ymin - gpad
            ymax = ymax + gpad
        end if

        ddx = (xmax - xmin)/nx
        ddy = (ymax - ymin)/ny

        allocate (grid(nx, ny))
        grid = 0.0_sp

        do iiy = 1, ny
            yc = ymin + (iiy - 0.5_dp)*ddy
            do iix = 1, nx
                xc = xmin + (iix - 0.5_dp)*ddx
                dmin = huge(dmin)
                do i = 1, n - 1
                    d = seg_dist(xc, yc, x(i), y(i), x(i + 1), y(i + 1))
                    if (d < dmin) dmin = d
                end do
                if (dmin < 0.5_dp*W) grid(iix, iiy) = 1.0_sp
            end do
        end do
    end subroutine project_channel_to_grid

    !> Minimum distance from point (px,py) to segment (ax,ay)-(bx,by).
    pure function seg_dist(px, py, ax, ay, bx, by) result(d)
        real(dp), intent(in) :: px, py, ax, ay, bx, by
        real(dp) :: d, t, ex, ey, len2
        ex = bx - ax
        ey = by - ay
        len2 = ex*ex + ey*ey
        if (len2 < 1.0e-20_dp) then
            d = hypot(px - ax, py - ay)
        else
            t = max(0.0_dp, min(1.0_dp, ((px - ax)*ex + (py - ay)*ey)/len2))
            d = hypot(px - ax - t*ex, py - ay - t*ey)
        end if
    end function seg_dist

    ! ── 3-D channel volume ───────────────────────────────────────────────────────

    !> Simulate one centreline and project it to a 3-D grid where channel width
    !> increases from W_bot (bottom, iz=1) to W (top, iz=nz) via a power law:
    !>   W_iz = W_bot + (W - W_bot) * t**W_exp,   t = (iz-1)/(nz-1)
    !> W_exp < 1 (default 0.5 = sqrt) gives a rounded U-shape;
    !> W_exp = 1 gives a linear V-shape; W_exp > 1 gives a flared bell.
    !>
    !> vol(nx_g, ny_g, nz_g): 1.0 = channel, 0.0 = background.
    !> Fortran storage: x fastest, z slowest (column-major).
    !> When read in Python as order='F', reshape to (nz_g, ny_g, nx_g).
    subroutine simulate_meander_3d(length, W, D, Cf, k, kl, omega, gamma, &
            n_bends, dt, n_iter, seed, &
            nx_g, ny_g, nz_g, vol, &
            W_bot, W_exp, W_shape, W_render)

        real(dp), intent(in) :: length, W, D, Cf, k, kl, omega, gamma, dt
        integer, intent(in) :: n_bends, n_iter, seed, nx_g, ny_g, nz_g
        real(sp), allocatable, intent(out) :: vol(:, :, :)
        real(dp), intent(in), optional :: W_bot, W_exp, W_render
        integer, intent(in), optional :: W_shape

        real(dp), allocatable :: x(:), y(:)
        real(dp) :: umin, umax, vmin, vmax, gpad_v
        real(dp) :: u_c, v_c, dxe, dye, dlene
        real(dp) :: xc, yc, dmin, dseg, W_iz, W_bottom, W_exponent, t, W_r
        real(dp) :: gauss_lo, gauss_hi, gauss_denom, fval
        integer  :: i, nn, iix, iiy, iiz, npad, W_shape_val

        call simulate_meander_2d(length=length, W=W, D=D, Cf=Cf, k=k, kl=kl, &
            omega=omega, gamma=gamma, n_bends=n_bends, &
            dt=dt, n_iter=n_iter, seed=seed, &
            x_out=x, y_out=y, pad_out=npad)

        ! Trim fixed boundary nodes; OBB maps the free endpoints to the grid edges.
        nn = size(x)
        x = x(npad + 1:nn - npad)
        y = y(npad + 1:nn - npad)
        nn = size(x)

        dxe = x(nn) - x(1)
        dye = y(nn) - y(1)
        dlene = sqrt(dxe**2 + dye**2)
        if (dlene < 1.0e-10_dp) then
            dxe = 1.0_dp
            dye = 0.0_dp
        else
            dxe = dxe/dlene
            dye = dye/dlene
        end if
        umin = x(1)*dxe + y(1)*dye
        umax = x(nn)*dxe + y(nn)*dye
        vmin = minval(-x*dye + y*dxe)
        vmax = maxval(-x*dye + y*dxe)
        ! The rendered channel width may differ from the migration width so that
        ! wide channels can be rendered without degrading the centerline node count
        W_r = W
        if (present(W_render)) then
            if (W_render > 0.0_dp) then
                W_r = W_render
            end if
        end if

        gpad_v = 0.5_dp*W_r
        vmin = vmin - gpad_v
        vmax = vmax + gpad_v

        W_bottom = merge(W_bot, W_r/10.0_dp, present(W_bot))
        W_exponent = merge(W_exp, 0.5_dp, present(W_exp))
        W_shape_val = 0
        if (present(W_shape)) W_shape_val = W_shape

        ! Precompute shape-dependent constants (W_exponent = sigma for shapes 2 & 3)
        gauss_lo = 0.0_dp
        gauss_hi = 1.0_dp
        gauss_denom = 1.0_dp
        if (W_shape_val == 2) then
            gauss_lo = erf(-0.5_dp/(max(1.0e-6_dp, W_exponent)*sqrt(2.0_dp)))
            gauss_hi = erf(0.5_dp/(max(1.0e-6_dp, W_exponent)*sqrt(2.0_dp)))
            gauss_denom = max(1.0e-10_dp, gauss_hi - gauss_lo)
        else if (W_shape_val == 3) then
            ! gauss_lo repurposed as exp(-1/sigma²) for the bell cross-section
            gauss_lo = exp(-1.0_dp/max(1.0e-6_dp, W_exponent)**2)
        end if

        allocate (vol(nx_g, ny_g, nz_g))
        vol = 0.0_sp

        do iiy = 1, ny_g
            do iix = 1, nx_g
                u_c = umin + (real(iix, dp) - 0.5_dp)/real(nx_g, dp)*(umax - umin)
                v_c = vmin + (real(iiy, dp) - 0.5_dp)/real(ny_g, dp)*(vmax - vmin)
                xc = u_c*dxe - v_c*dye
                yc = u_c*dye + v_c*dxe
                dmin = huge(dmin)
                do i = 1, nn - 1
                    dseg = seg_dist(xc, yc, x(i), y(i), x(i + 1), y(i + 1))
                    if (dseg < dmin) dmin = dseg
                end do
                do iiz = 1, nz_g
                    if (nz_g > 1) then
                        t = (iiz - 1.0_dp)/(nz_g - 1.0_dp)
                        select case (W_shape_val)
                            case (1)  ! smoothstep: flat bottom, steep walls, gentle shoulders
                                W_iz = W_bottom + (W_r - W_bottom)*(t*t*(3.0_dp - 2.0_dp*t))
                            case (2)  ! Gaussian CDF: steepness set by W_exp as sigma
                                fval = (erf((t - 0.5_dp)/(max(1.0e-6_dp, W_exponent)*sqrt(2.0_dp))) &
                                    - gauss_lo)/gauss_denom
                                W_iz = W_bottom + (W_r - W_bottom)*fval
                            case (3)  ! Gaussian bell: z(x) profile is a Gaussian bell; W_exp = sigma (0.3-0.7)
                                fval = W_exponent*sqrt(max(0.0_dp, &
                                    -log((1.0_dp - t)*(1.0_dp - gauss_lo) + gauss_lo)))
                                W_iz = W_bottom + (W_r - W_bottom)*min(1.0_dp, fval)
                            case default  ! power law (original)
                                W_iz = W_bottom + (W_r - W_bottom)*t**W_exponent
                        end select
                    else
                        W_iz = W_r
                    end if
                    if (dmin < 0.5_dp*W_iz) vol(iix, iiy, iiz) = 1.0_sp
                end do
            end do
        end do

    end subroutine

    ! ── submarine canyon volume ──────────────────────────────────────────────────

    !> Simulate meander migration while saving channel snapshots at every
    !> save_every iterations.  Each snapshot is assigned a depth level that
    !> follows a V-shaped incision-then-aggradation trajectory:
    !>
    !>   iz = 1 + round((nz_g-1) * (1 - 2*|t - 0.5|)),   t = snap_index / n_snap
    !>
    !> so early and late snapshots land at iz = nz_g (seabed surface) and the
    !> mid-simulation snapshot (maximum meander development) lands at iz = 1
    !> (deepest incision).  The union of all painted capsules gives a canyon body
    !> that is narrow at depth and widens upward.
    subroutine simulate_meander_canyon_3d(length, W, D, Cf, k, kl, omega, gamma, &
            n_bends, dt, n_iter, seed, &
            nx_g, ny_g, nz_g, vol, &
            save_every, W_canyon, &
            n_stages, dz_sigma)

        real(dp), intent(in) :: length, W, D, Cf, k, kl, omega, gamma, dt
        integer, intent(in) :: n_bends, n_iter, seed, nx_g, ny_g
        integer, intent(in) :: nz_g   ! z-levels; pass 0 to auto-set = n_snap
        real(sp), allocatable, intent(out) :: vol(:, :, :)
        integer, intent(in), optional :: save_every
        real(dp), intent(in), optional :: W_canyon  ! canyon width at surface; default = W
        integer, intent(in), optional :: n_stages  ! number of longitudinal profile stages
        real(dp), intent(in), optional :: dz_sigma  ! z-profile amplitude as fraction of nz (e.g. 0.15)

        type :: snap_t
            real(dp), allocatable :: x(:), y(:)
        end type snap_t
        type(snap_t), allocatable :: snaps(:)

        real(dp), allocatable :: x(:), y(:), raw(:), rd(:)
        real(dp), allocatable :: dx(:), dy(:), ddx(:), ddy(:)
        real(dp), allocatable :: den(:), curv(:), R0(:), Omg(:), R1(:)
        real(dp), allocatable :: ds(:), tx(:), ty(:)
        real(dp), allocatable :: dpx(:), dpy(:), mag(:), sc(:)

        real(dp) :: crdist_val, deltas, alpha, decay, phi_k
        real(dp) :: kap, th, pos, slope_val, rms_val, arc, chord
        integer  :: n0, n, it, j, npad, n_ox, n_ox_total

        integer  :: save_every_val, n_snap_max, i_snap, i, iix, iiy, nn
        integer  :: r_ix, r_iy, ix_lo, ix_hi, iy_lo, iy_hi, nz_actual
        integer  :: n_stages_val, s1_stage, s2_stage, ic, iz_j
        integer, allocatable :: min_iz_map(:, :)
        real(dp) :: t, z_frac, xmin, xmax, ymin, ymax, ddx_g, ddy_g
        real(dp) :: xc, yc, ax, ay, bx_seg, by_seg
        real(dp) :: W_canyon_val, W_eff
        real(dp) :: dz_sigma_val, frac_j, dz_j, dz_j1, dz_j2, alpha_j
        real(dp) :: iz_mean_dp, t_stage, blend_stage
        real(dp), allocatable :: stage_profiles(:, :), rd_tmp(:)

        save_every_val = 50
        if (present(save_every)) save_every_val = max(1, save_every)
        n_snap_max = n_iter/save_every_val + 4

        ! ── identical setup as simulate ───────────────────────────────────────────
        crdist_val = 1.5_dp*W
        deltas = W/2.0_dp
        alpha = k*2.0_dp*Cf/D
        decay = exp(-alpha*deltas)

        n0 = max(12, nint(length/deltas) + 1)
        npad = max(10, nint(n0*0.05_dp))
        phi_k = exp(-deltas/(length/max(1, n_bends)*0.4_dp))

        allocate (raw(n0))
        kap = 0.0_dp
        th = 0.0_dp
        pos = 0.0_dp
        rd = random(n0, dist='normal', seed=seed)
        do j = 1, n0
            kap = phi_k*kap + sqrt(1.0_dp - phi_k**2)*rd(j)
            th = th + kap
            th = max(-1.2_dp, min(1.2_dp, th))
            pos = pos + sin(th)
            raw(j) = pos
        end do
        slope_val = (raw(n0) - raw(1))/real(n0 - 1, dp)
        do j = 1, n0
            raw(j) = raw(j) - slope_val*real(j - 1, dp)
        end do
        raw = raw - sum(raw)/real(n0, dp)
        rms_val = sqrt(sum(raw**2)/real(n0, dp))
        if (rms_val > 1.0e-10_dp) raw = raw*(2.0_dp*W/rms_val)

        allocate (x(n0), y(n0))
        do j = 1, n0
            x(j) = real(j - 1, dp)*length/real(n0 - 1, dp)
            y(j) = raw(j)
        end do
        deallocate (raw)
        call resample_chan(x, y, deltas)
        n_ox_total = 0

        allocate (snaps(n_snap_max))
        i_snap = 0

        ! save initial snapshot (t=0, will be placed at seabed surface)
        nn = size(x)
        if (nn > 2*npad + 1) then
            i_snap = i_snap + 1
            snaps(i_snap)%x = x(npad + 1:nn - npad)
            snaps(i_snap)%y = y(npad + 1:nn - npad)
        end if

        ! ── main migration loop with snapshot saving ───────────────────────────────
        do it = 1, n_iter
            n = size(x)
            if (n < 2*npad + 4) exit

            allocate (dx(n), dy(n), ddx(n), ddy(n), den(n), curv(n))
            call grad1d(x, dx)
            call grad1d(y, dy)
            call grad1d(dx, ddx)
            call grad1d(dy, ddy)
            den = (dx**2 + dy**2)**1.5_dp
            where (den > 1.0e-12_dp)
                curv = (dx*ddy - dy*ddx)/den
            elsewhere
                curv = 0.0_dp
            end where

            allocate (R0(n), Omg(n), R1(n))
            R0 = kl*W*curv
            call iir_causal(R0, decay, Omg)
            R1 = omega*R0 + gamma*Omg

            allocate (ds(n))
            ds = hypot(dx, dy)
            arc = sum(ds)
            chord = max(hypot(x(n) - x(1), y(n) - y(1)), 1.0e-10_dp)
            R1 = R1*(chord/arc)**(2.0_dp/3.0_dp)
            ds = max(ds, 1.0e-12_dp)

            allocate (tx(n), ty(n))
            tx = dx/ds
            ty = dy/ds

            allocate (dpx(n - 2*npad), dpy(n - 2*npad))
            dpx = R1(npad + 1:n - npad)*ty(npad + 1:n - npad)*dt
            dpy = -R1(npad + 1:n - npad)*tx(npad + 1:n - npad)*dt

            allocate (mag(n - 2*npad), sc(n - 2*npad))
            mag = hypot(dpx, dpy)
            where (mag > 0.5_dp*deltas)
                sc = 0.5_dp*deltas/(mag + 1.0e-10_dp)
            elsewhere
                sc = 1.0_dp
            end where

            x(npad + 1:n - npad) = x(npad + 1:n - npad) + dpx*sc
            y(npad + 1:n - npad) = y(npad + 1:n - npad) + dpy*sc

            deallocate (dx, dy, ddx, ddy, den, curv, R0, Omg, R1, ds, tx, ty, dpx, dpy, mag, sc)

            call resample_chan(x, y, deltas)
            call cut_off_cutoffs(x, y, crdist_val, deltas, n_ox)
            n_ox_total = n_ox_total + n_ox

            if (mod(it, save_every_val) == 0 .and. i_snap < n_snap_max) then
                nn = size(x)
                if (nn > 2*npad + 1) then
                    i_snap = i_snap + 1
                    snaps(i_snap)%x = x(npad + 1:nn - npad)
                    snaps(i_snap)%y = y(npad + 1:nn - npad)
                end if
            end if
        end do

        ! also save the final state
        nn = size(x)
        if (nn > 2*npad + 1 .and. i_snap < n_snap_max) then
            i_snap = i_snap + 1
            snaps(i_snap)%x = x(npad + 1:nn - npad)
            snaps(i_snap)%y = y(npad + 1:nn - npad)
        end if

        if (i_snap == 0) then
            allocate (vol(nx_g, ny_g, nz_g))
            vol = 0.0_sp
            return
        end if

        ! ── grid bounds ───────────────────────────────────────────────────────────
        xmin = huge(xmin)
        xmax = -huge(xmax)
        ymin = huge(ymin)
        ymax = -huge(ymax)
        do i = 1, i_snap
            xmin = min(xmin, minval(snaps(i)%x))
            xmax = max(xmax, maxval(snaps(i)%x))
            ymin = min(ymin, minval(snaps(i)%y))
            ymax = max(ymax, maxval(snaps(i)%y))
        end do
        ymin = ymin - 0.5_dp*W
        ymax = ymax + 0.5_dp*W

        ddx_g = (xmax - xmin)/nx_g
        ddy_g = (ymax - ymin)/ny_g

        ! canyon surface width; when omitted the width stays W at all depths and
        ! natural meander migration widens the canyon cross-section via cumulative fill
        if (present(W_canyon)) then
            W_canyon_val = max(W, W_canyon)
        else
            W_canyon_val = W
        end if
        r_ix = ceiling(0.5_dp*W_canyon_val/ddx_g) + 1
        r_iy = ceiling(0.5_dp*W_canyon_val/ddy_g) + 1

        ! When nz_g = 0 (auto), set nz_actual = i_snap
        if (nz_g <= 0) then
            nz_actual = i_snap
        else
            nz_actual = nz_g
        end if

        ! ── staged longitudinal z-profiles ───────────────────────────────────────
        ! Each stage covers a contiguous block of snapshots.  Within a stage every
        ! segment of the channel has the same mean iz (from the V-trajectory) PLUS
        ! a smooth random offset that varies along the channel arc.  The offset is
        ! interpolated from n_ctrl_pts control values drawn once per stage so that
        ! the profile is consistent across all snapshots in that stage.
        n_stages_val = 10
        if (present(n_stages)) n_stages_val = max(1, n_stages)
        dz_sigma_val = 0.0_dp
        if (present(dz_sigma)) dz_sigma_val = max(0.0_dp, min(0.5_dp, dz_sigma))

        ! n_ctrl_pts: number of control knots for the longitudinal profile (fixed)
        ! stage_profiles(ctrl, stage) = dz offset in z-cells at that knot/stage
        allocate (stage_profiles(10, n_stages_val))
        do s1_stage = 1, n_stages_val
            rd_tmp = random(10, dist='normal', seed=safe_seed(int(seed, 8) + s1_stage*7919))
            stage_profiles(:, s1_stage) = rd_tmp*(dz_sigma_val*real(nz_actual, dp))
        end do

        allocate (vol(nx_g, ny_g, nz_actual))
        vol = 0.0_sp

        ! ── erosional cumulative fill with per-segment iz ─────────────────────────
        ! V-trajectory: t=0 → iz_mean=nz_actual (surface), t=0.5 → iz_mean=1
        ! (maximum incision / most sinuous), t=1 → iz_mean=nz_actual (surface).
        ! Stage profiles add a smooth random dz along the channel arc, cross-faded
        ! between adjacent stages.  min_iz_map records the deepest iz any segment
        ! carved each XY cell; vol fills upward from there to nz_actual.
        allocate (min_iz_map(nx_g, ny_g))
        min_iz_map = nz_actual + 1  ! sentinel: not carved

        do i = 1, i_snap
            t = real(i - 1, dp)/real(max(1, i_snap - 1), dp)
            z_frac = 1.0_dp - 2.0_dp*abs(t - 0.5_dp)
            iz_mean_dp = real(nz_actual, dp) - real(nz_actual - 1, dp)*z_frac

            t_stage = t*real(n_stages_val, dp)
            s1_stage = min(n_stages_val, 1 + int(t_stage))
            s2_stage = min(n_stages_val, s1_stage + 1)
            blend_stage = t_stage - real(s1_stage - 1, dp)
            blend_stage = blend_stage*blend_stage*(3.0_dp - 2.0_dp*blend_stage)

            nn = size(snaps(i)%x)
            do j = 1, nn - 1
                frac_j = real(j, dp)/real(max(1, nn), dp)
                alpha_j = frac_j*9.0_dp
                ic = min(9, int(alpha_j)) + 1
                alpha_j = alpha_j - real(ic - 1, dp)
                dz_j1 = (1.0_dp - alpha_j)*stage_profiles(ic, s1_stage) &
                    + alpha_j*stage_profiles(ic + 1, s1_stage)
                dz_j2 = (1.0_dp - alpha_j)*stage_profiles(ic, s2_stage) &
                    + alpha_j*stage_profiles(ic + 1, s2_stage)
                dz_j = ((1.0_dp - blend_stage)*dz_j1 + blend_stage*dz_j2)*z_frac
                iz_j = max(1, min(nz_actual, nint(iz_mean_dp + dz_j)))

                W_eff = W + (W_canyon_val - W)*real(iz_j - 1, dp) &
                    /real(max(1, nz_actual - 1), dp)

                ax = (snaps(i)%x(j) - xmin)/ddx_g
                ay = (snaps(i)%y(j) - ymin)/ddy_g
                bx_seg = (snaps(i)%x(j + 1) - xmin)/ddx_g
                by_seg = (snaps(i)%y(j + 1) - ymin)/ddy_g
                ix_lo = max(1, nint(min(ax, bx_seg)) - r_ix)
                ix_hi = min(nx_g, nint(max(ax, bx_seg)) + r_ix)
                iy_lo = max(1, nint(min(ay, by_seg)) - r_iy)
                iy_hi = min(ny_g, nint(max(ay, by_seg)) + r_iy)
                do iiy = iy_lo, iy_hi
                    yc = ymin + (iiy - 0.5_dp)*ddy_g
                    do iix = ix_lo, ix_hi
                        xc = xmin + (iix - 0.5_dp)*ddx_g
                        if (seg_dist(xc, yc, snaps(i)%x(j), snaps(i)%y(j), &
                                snaps(i)%x(j + 1), snaps(i)%y(j + 1)) < 0.5_dp*W_eff) then
                            if (iz_j < min_iz_map(iix, iiy)) min_iz_map(iix, iiy) = iz_j
                        end if
                    end do
                end do
            end do
        end do
        deallocate (stage_profiles)

        do iiy = 1, ny_g
            do iix = 1, nx_g
                if (min_iz_map(iix, iiy) <= nz_actual) then
                    vol(iix, iiy, min_iz_map(iix, iiy):nz_actual) = 1.0_sp
                end if
            end do
        end do
        deallocate (min_iz_map)

        do i = 1, i_snap
            if (allocated(snaps(i)%x)) deallocate (snaps(i)%x)
            if (allocated(snaps(i)%y)) deallocate (snaps(i)%y)
        end do
        deallocate (snaps)

    end subroutine

    ! ── OOP generate functions ───────────────────────────────────────────────────

    !> Generate a 2-D depth map from a 3-D meandering channel volume.
    !> Produces depth_map(ny,nx), downward-positive: 0 = highest topography.
    !> Depth is measured from the channel bottom: centre cells score ~1.0,
    !> bank-edge cells score a small positive fraction.
    subroutine generate_channel(self)

        class(meandering_channel), intent(inout) :: self
        integer  :: iix, iiy, iiz, nz_g, eff_seed, tmp_n, ori
        real(dp) :: W_bot_val
        real(dp), allocatable :: dem_bg(:, :)

        eff_seed = safe_seed(int(self%seed, 8))

        ! Output orientation: fixed by orient if given, otherwise derived from the seed
        if (self%orient >= 0) then
            ori = mod(self%orient, 4)
        else
            ori = mod(abs(eff_seed), 4)
        end if

        ! For 90° rotations: pre-swap nx↔ny so the whole pipeline runs in the
        ! transposed frame; swapped back after the vol rebuild so output stays (ny,nx).
        select case (ori)
            case (2, 3)
                tmp_n = self%nx
                self%nx = self%ny
                self%ny = tmp_n
        end select

        W_bot_val = self%W_bot
        if (W_bot_val <= 0.0_dp) then
            if (self%W_render > 0.0_dp) then
                W_bot_val = self%W_render/2.0_dp
            else
                W_bot_val = self%w/2.0_dp
            end if
        end if

        call simulate_meander_3d(self%length, self%W, self%D, self%Cf, &
            self%k, self%kl, self%omega, self%gamma, &
            self%n_bends, self%dt, self%n_iter, eff_seed, &
            self%nx, self%ny, self%nz, self%vol, &
            W_bot=w_bot_val, W_exp=self%W_exp, W_shape=self%W_shape, &
            W_render=self%W_render)

        nz_g = size(self%vol, 3)
        if (allocated(self%depth_map)) then
            deallocate (self%depth_map)
        end if
        allocate (self%depth_map(self%nx, self%ny))
        self%depth_map = 0.0_sp

        !$omp parallel do private(iix, iiy, iiz)
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
        !$omp end parallel do

        self%depth_map = self%depth_map*(nz_g - 1.0)
        self%depth_map = transpose(self%depth_map)
        if (.not. self%yn_depth_only) then
            self%vol = flip(permute(self%vol, 321), [1])
        end if

        ! DEM
        allocate (dem_bg(self%nx, self%ny))
        dem_bg = 0.0
        if (self%terrain_bg > 0.0_dp) then
            dem_bg = rescale(gauss_filt(random(self%nx, self%ny, dist='normal', seed=safe_seed(int(eff_seed, 8)*7 + 11)), &
                [0.1, 0.1]*0.5*(self%nx + self%ny)), [0.0, 1.0])
            dem_bg = real(transpose(dem_bg)*self%terrain_bg, sp)
            ! Background cells get terrain value [0, terrain_bg].
            ! Channel cells get terrain − channel_depth (negative = incised below surface).
            where (self%depth_map == 0.0_sp)
                self%depth_map = dem_bg
            elsewhere
                ! Here, different from drainage, use the original depth map for channel pixels
                ! This will result in mesa-like topography
                ! self%depth_map = dem_bg - self%depth_map
                self%depth_map = -self%depth_map
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

        ! Rotate output to final orientation; for cases 2/3 also swap nx↔ny back.
        select case (ori)
            case (1)
                self%depth_map = flip(self%depth_map, [2])
                if (.not. self%yn_depth_only) then
                    self%vol = flip(self%vol, [3])
                    self%obj = flip(self%obj, [3])
                end if
            case (2)
                self%depth_map = transpose(self%depth_map)
                if (.not. self%yn_depth_only) then
                    self%vol = permute(self%vol, 132)
                    self%obj = permute(self%obj, 132)
                end if
                tmp_n = self%nx
                self%nx = self%ny
                self%ny = tmp_n
            case (3)
                self%depth_map = flip(transpose(self%depth_map), [2])
                if (.not. self%yn_depth_only) then
                    self%vol = permute(self%vol, 132)
                    self%vol = flip(self%vol, [3])
                    self%obj = permute(self%obj, 132)
                    self%obj = flip(self%obj, [3])
                end if
                tmp_n = self%nx
                self%nx = self%ny
                self%ny = tmp_n
        end select

    end subroutine

    !> Generate a 2-D depth map from a 3-D meander canyon volume.
    !> Produces depth_map(ny,nx), downward-positive: 0 = highest topography.
    !> Depth encodes how deep the canyon floor was carved at each (x,y) cell.
    subroutine generate_canyon(self)

        class(meandering_canyon), intent(inout) :: self
        integer  :: iix, iiy, iiz, nz_g, eff_seed, tmp_n, nx_original, npad, ori
        real(dp) :: W_canyon_val
        real(dp), allocatable :: dem_bg(:, :)

        eff_seed = safe_seed(int(self%seed, 8))

        ! Output orientation: fixed by orient if given, otherwise derived from the seed
        if (self%orient >= 0) then
            ori = mod(self%orient, 4)
        else
            ori = mod(abs(eff_seed), 4)
        end if

        ! For 90° rotations: pre-swap nx↔ny so the whole pipeline runs in the
        ! transposed frame; swapped back after the vol rebuild so output stays (ny,nx).
        select case (ori)
            case (2, 3)
                tmp_n = self%nx
                self%nx = self%ny
                self%ny = tmp_n
        end select

        ! Extend the grid in the start/end (x) direction so the canyon flows
        ! naturally off the edges; cropped back after rendering.
        nx_original = self%nx
        npad = nint(0.1_dp*real(self%nx, dp))
        self%nx = self%nx + 2*npad

        W_canyon_val = max(self%W_canyon, self%W)

        call simulate_meander_canyon_3d(self%length, self%W, self%D, self%Cf, &
            self%k, self%kl, self%omega, self%gamma, &
            self%n_bends, self%dt, self%n_iter, eff_seed, &
            self%nx, self%ny, self%nz, self%vol, &
            save_every=self%save_every, &
            W_canyon=W_canyon_val, &
            n_stages=self%n_stages, dz_sigma=self%dz_sigma)

        nz_g = size(self%vol, 3)
        if (allocated(self%depth_map)) then
            deallocate (self%depth_map)
        end if
        allocate (self%depth_map(self%nx, self%ny))
        self%depth_map = 0.0_sp

        !$omp parallel do private(iix, iiy, iiz)
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
        !$omp end parallel do

        self%depth_map = self%depth_map*(nz_g - 1.0)

        ! Crop start/end padding and restore original nx.
        self%depth_map = self%depth_map(npad + 1:npad + nx_original, :)
        self%nx = nx_original

        self%depth_map = transpose(self%depth_map)
        if (.not. self%yn_depth_only) then
            self%vol = flip(permute(self%vol, 321), [1])
            ! Crop the x-dimension (3rd dim after permute) to match.
            self%vol = self%vol(:, :, npad + 1:npad + nx_original)
        end if

        ! DEM
        allocate (dem_bg(self%nx, self%ny))
        dem_bg = 0.0
        if (self%terrain_bg > 0.0_dp) then
            dem_bg = rescale(gauss_filt(random(self%nx, self%ny, dist='normal', seed=safe_seed(int(eff_seed, 8)*7 + 11)), &
                [0.1, 0.1]*0.5*(self%nx + self%ny)), [0.0, 1.0])
            dem_bg = real(transpose(dem_bg)*self%terrain_bg, sp)
            ! Background cells get terrain value [0, terrain_bg].
            ! Channel cells get terrain − channel_depth (negative = incised below surface).
            where (self%depth_map == 0.0_sp)
                self%depth_map = dem_bg
            elsewhere
                ! self%depth_map = dem_bg - self%depth_map
                self%depth_map = -self%depth_map
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

        ! Rotate output to final orientation; for cases 2/3 also swap nx↔ny back.
        select case (ori)
            case (1)
                self%depth_map = flip(self%depth_map, [2])
                if (.not. self%yn_depth_only) then
                    self%vol = flip(self%vol, [3])
                    self%obj = flip(self%obj, [3])
                end if
            case (2)
                self%depth_map = transpose(self%depth_map)
                if (.not. self%yn_depth_only) then
                    self%vol = permute(self%vol, 132)
                    self%obj = permute(self%obj, 132)
                end if
                tmp_n = self%nx
                self%nx = self%ny
                self%ny = tmp_n
            case (3)
                self%depth_map = flip(transpose(self%depth_map), [2])
                if (.not. self%yn_depth_only) then
                    self%vol = permute(self%vol, 132)
                    self%vol = flip(self%vol, [3])
                    self%obj = permute(self%obj, 132)
                    self%obj = flip(self%obj, [3])
                end if
                tmp_n = self%nx
                self%nx = self%ny
                self%ny = tmp_n
        end select

    end subroutine

end module
