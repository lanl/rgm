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
! Supporting routines of the fold-train and skewed reflector shapes, the
! dip-independent imaging, and the migration-like image noise and effects of
! rgm2_curved and rgm3_curved
!
module geological_model_realism

    use libflit
    use geological_model_utility

    implicit none

    interface release_unless
        module procedure :: release_unless_2d
        module procedure :: release_unless_3d
    end interface release_unless

    private
    public :: derive_seed
    public :: smooth_noise_1d
    public :: smooth_noise_2d
    public :: smooth_noise_3d
    public :: fold_train
    public :: fold_surface
    public :: skewed_bump_1d
    public :: skewed_bump_2d
    public :: dip_cosine_2d
    public :: dip_cosine_3d
    public :: isotropic_psf_2d
    public :: isotropic_psf_3d
    public :: noise_band_2d
    public :: noise_band_3d
    public :: noise_worm_2d
    public :: noise_worm_3d
    public :: noise_swing_2d
    public :: noise_swing_3d
    public :: illumination_2d
    public :: illumination_3d
    public :: trace_jitter_2d
    public :: trace_jitter_3d
    public :: fold_reflector_1d
    public :: fold_reflector_2d
    public :: bump_reflector_1d
    public :: bump_reflector_2d
    public :: layer_dip_cosine_2d
    public :: layer_dip_cosine_3d
    public :: migration_noise_2d
    public :: migration_noise_3d
    public :: depth_blend_2d
    public :: depth_blend_3d
    public :: release_unless

contains

    !
    !> Derive a sub-seed from a seed and an index; a negative seed (pure
    !> randomness) stays negative
    !
    pure function derive_seed(seed, k) result(s)

        integer, intent(in) :: seed, k
        integer :: s

        if (seed < 0) then
            s = -1
        else
            s = safe_seed(int(seed, 8)*131 + k)
        end if

    end function derive_seed

    !
    !> Radius of libflit's recursive Gaussian filter (gauss_filt, method 'rsf')
    !> for sigma, plus a safety margin
    !
    pure function gauss_radius(sigma) result(r)

        real, intent(in) :: sigma
        integer :: r

        r = nint(3.0*sigma + 0.5) + 2

    end function gauss_radius

    !
    !> 1D Gaussian smoothing that is valid for any sigma. libflit's gauss_filt
    !> reflects indices only once at each end, so it reads beyond the array when
    !> its radius (about 3*sigma) exceeds the array length; the array is then
    !> extended with its edge values before filtering. An array with sigma <= 0
    !> is returned unchanged.
    !
    function gauss_smooth_1d(w, sigma) result(ws)

        real, dimension(:), intent(in) :: w
        real, intent(in) :: sigma
        real, allocatable, dimension(:) :: ws

        integer :: n, p

        n = size(w)
        if (sigma <= 0 .or. n < 2) then
            ws = w
        else if (gauss_radius(sigma) <= n - 2) then
            ws = gauss_filt(w, sigma)
        else
            p = (gauss_radius(sigma) - n)/2 + 2
            ws = gauss_filt(pad(w, [p, p], ['edge', 'edge']), sigma)
            ws = ws(p + 1:p + n)
        end if

    end function gauss_smooth_1d

    !
    !> Gaussian smoothing applied as successive 1D filters along the axes
    !> (gauss_smooth_1d), which stays valid when sigma is large relative to the
    !> array size; an axis with sigma <= 0 is not smoothed
    !
    function gauss_smooth_2d(w, sigma) result(ws)

        real, dimension(:, :), intent(in) :: w
        real, dimension(1:2), intent(in) :: sigma
        real, allocatable, dimension(:, :) :: ws

        integer :: i, j

        ws = w
        if (sigma(1) > 0) then
            !$omp parallel do private(j)
            do j = 1, size(ws, 2)
                ws(:, j) = gauss_smooth_1d(ws(:, j), sigma(1))
            end do
            !$omp end parallel do
        end if
        if (sigma(2) > 0) then
            !$omp parallel do private(i)
            do i = 1, size(ws, 1)
                ws(i, :) = gauss_smooth_1d(ws(i, :), sigma(2))
            end do
            !$omp end parallel do
        end if

    end function gauss_smooth_2d

    function gauss_smooth_3d(w, sigma) result(ws)

        real, dimension(:, :, :), intent(in) :: w
        real, dimension(1:3), intent(in) :: sigma
        real, allocatable, dimension(:, :, :) :: ws

        integer :: i, j, k

        ws = w
        if (sigma(1) > 0) then
            !$omp parallel do private(j, k) collapse(2)
            do k = 1, size(ws, 3)
                do j = 1, size(ws, 2)
                    ws(:, j, k) = gauss_smooth_1d(ws(:, j, k), sigma(1))
                end do
            end do
            !$omp end parallel do
        end if
        if (sigma(2) > 0) then
            !$omp parallel do private(i, k) collapse(2)
            do k = 1, size(ws, 3)
                do i = 1, size(ws, 1)
                    ws(i, :, k) = gauss_smooth_1d(ws(i, :, k), sigma(2))
                end do
            end do
            !$omp end parallel do
        end if
        if (sigma(3) > 0) then
            !$omp parallel do private(i, j) collapse(2)
            do j = 1, size(ws, 2)
                do i = 1, size(ws, 1)
                    ws(i, j, :) = gauss_smooth_1d(ws(i, j, :), sigma(3))
                end do
            end do
            !$omp end parallel do
        end if

    end function gauss_smooth_3d

    !
    !> Length of the axis on which white noise is generated, so that smoothing
    !> it with sigma is valid without padding; the center part of length n is used
    !
    pure function noise_length(n, sigma) result(m)

        integer, intent(in) :: n
        real, intent(in) :: sigma
        integer :: m

        if (sigma > 0) then
            m = max(n, gauss_radius(sigma) + 2)
        else
            m = n
        end if

    end function noise_length

    !
    !> Zero-mean, unit-variance smooth random function; sigma is the standard
    !> deviation (in grid points) of the Gaussian smoothing
    !
    function smooth_noise_1d(n, sigma, seed) result(w)

        integer, intent(in) :: n, seed
        real, intent(in) :: sigma
        real, allocatable, dimension(:) :: w

        integer :: m, o

        m = noise_length(n, sigma)
        o = (m - n)/2
        w = gauss_smooth_1d(random(m, dist='normal', seed=seed), sigma)
        w = w(o + 1:o + n)
        w = w - mean(w)
        w = w/(std(w) + float_tiny)

    end function smooth_noise_1d

    function smooth_noise_2d(n1, n2, sigma, seed) result(w)

        integer, intent(in) :: n1, n2, seed
        real, dimension(1:2), intent(in) :: sigma
        real, allocatable, dimension(:, :) :: w

        integer :: m1, m2, o1, o2

        m1 = noise_length(n1, sigma(1))
        m2 = noise_length(n2, sigma(2))
        o1 = (m1 - n1)/2
        o2 = (m2 - n2)/2
        w = gauss_smooth_2d(random(m1, m2, dist='normal', seed=seed), sigma)
        w = w(o1 + 1:o1 + n1, o2 + 1:o2 + n2)
        w = w - mean(w)
        w = w/(std(w) + float_tiny)

    end function smooth_noise_2d

    function smooth_noise_3d(n1, n2, n3, sigma, seed) result(w)

        integer, intent(in) :: n1, n2, n3, seed
        real, dimension(1:3), intent(in) :: sigma
        real, allocatable, dimension(:, :, :) :: w

        integer :: m1, m2, m3, o1, o2, o3

        m1 = noise_length(n1, sigma(1))
        m2 = noise_length(n2, sigma(2))
        m3 = noise_length(n3, sigma(3))
        o1 = (m1 - n1)/2
        o2 = (m2 - n2)/2
        o3 = (m3 - n3)/2
        w = gauss_smooth_3d(random(m1, m2, m3, dist='normal', seed=seed), sigma)
        w = w(o1 + 1:o1 + n1, o2 + 1:o2 + n2, o3 + 1:o3 + n3)
        w = w - mean(w)
        w = w/(std(w) + float_tiny)

    end function smooth_noise_3d

    !
    !> Fold train with a crest asymmetry and a limb asymmetry (vergence). The
    !> reflector is amp(x)*prof(x), with
    !>     prof = -cos(psi) + crest*cos(2*psi),
    !> where psi is the fold phase measured from a trough (psi = 0 at troughs and
    !> pi at crests) and d(theta)/dx = 2*pi/lambda(x). For mode = 'smooth',
    !> psi = theta + verg*sin(theta) + pi/2; for mode = 'limb', psi is a
    !> piecewise-linear warp of theta in which the rising limb takes (1 - verg)/2
    !> of each period, so that each limb is uniformly steep or gentle. crest > 0
    !> gives sharp anticlines and broad synclines, crest < 0 broad (box)
    !> anticlines and sharp synclines; |crest| <= 0.25 keeps the profile
    !> monotonic between hinges. The wavelength and the amplitude drift from fold
    !> to fold: lambda(x) = lam*exp(dlam*tanh(n1(x))) and
    !> amp(x) = 1 + damp*tanh(n2(x)), where n1 and n2 are smooth random functions
    !> with correlation lengths lam and 0.7*lam.
    !
    subroutine fold_train(n, lam, verg, crest, dlam, damp, mode, seed, prof, amp)

        integer, intent(in) :: n, seed
        real, intent(in) :: lam, verg, crest, dlam, damp
        character(len=*), intent(in) :: mode
        real, allocatable, dimension(:), intent(out) :: prof, amp

        real, allocatable, dimension(:) :: theta, psi
        real :: a, pi
        integer :: i

        pi = real(const_pi)

        theta = 2*pi*cumsum(1.0/(lam*exp(dlam*tanh(smooth_noise_1d(n, lam, derive_seed(seed, 1)))))) &
            + rand(range=[0.0, 2*pi], seed=derive_seed(seed, 3))
        amp = 1.0 + damp*tanh(smooth_noise_1d(n, 0.7*lam, derive_seed(seed, 2)))

        psi = zeros(n)
        select case (mode)
            case ('limb')
                a = pi*(1.0 - verg)
                do i = 1, n
                    psi(i) = modulo(theta(i) + 0.5*pi, 2*pi)
                    if (psi(i) < a) then
                        psi(i) = psi(i)/(1.0 - verg)
                    else
                        psi(i) = pi + (psi(i) - a)/(1.0 + verg)
                    end if
                end do
            case default
                psi = theta + verg*sin(theta) + 0.5*pi
        end select

        prof = -cos(psi) + crest*cos(2.0*psi)

    end subroutine fold_train

    !
    !> Fold-train surface for a 3D model. The fold train runs across the fold
    !> axes, whose strike is measured from x2 toward x3. Along strike, the fold
    !> uplift above the trough level is scaled by a smooth positive random field
    !> of relative variation plunge, correlated over plunge_length*lam along
    !> strike and 0.5*lam across strike, so each anticline rises and plunges out
    !> at its own place; the fold axes are shifted laterally by a smooth random
    !> function of the along-strike position with amplitude wobble*lam.
    !
    function fold_surface(n2, n3, lam, strike, verg, crest, dlam, damp, mode, &
            plunge, plunge_length, wobble, seed) result(r)

        integer, intent(in) :: n2, n3, seed
        real, intent(in) :: lam, strike, verg, crest, dlam, damp, plunge, plunge_length, wobble
        character(len=*), intent(in) :: mode
        real, allocatable, dimension(:, :) :: r

        real, allocatable, dimension(:, :) :: u, v, g
        real, allocatable, dimension(:) :: prof, amp, wob
        real :: x, y, a, b, pa, pr, am, env, trough
        integer :: j, k, nu, nv, m, l

        ! Positions across (u) and along (v) the fold axes
        u = zeros(n2, n3)
        v = zeros(n2, n3)
        do k = 1, n3
            do j = 1, n2
                u(j, k) = (j - 1.0)*sin(strike) - (k - 1.0)*cos(strike)
                v(j, k) = (j - 1.0)*cos(strike) + (k - 1.0)*sin(strike)
            end do
        end do
        v = v - minval(v)
        nv = ceiling(maxval(v)) + 2

        ! Lateral shift of the fold axes along strike
        wob = wobble*lam*smooth_noise_1d(nv, plunge_length*lam, derive_seed(seed, 11))
        do k = 1, n3
            do j = 1, n2
                y = v(j, k) + 1.0
                l = min(floor(y), nv - 1)
                b = y - l
                u(j, k) = u(j, k) + wob(l)*(1.0 - b) + wob(l + 1)*b
            end do
        end do
        u = u - minval(u)
        nu = ceiling(maxval(u)) + 2

        call fold_train(nu, lam, verg, crest, dlam, damp, mode, seed, prof, amp)
        g = smooth_noise_2d(nu, nv, [0.5*lam, plunge_length*lam], derive_seed(seed, 12))
        trough = -1.0 + crest

        r = zeros(n2, n3)
        !$omp parallel do private(j, k, x, y, m, l, a, b, pa, pr, am, env)
        do k = 1, n3
            do j = 1, n2
                x = u(j, k) + 1.0
                m = min(floor(x), nu - 1)
                a = x - m
                y = v(j, k) + 1.0
                l = min(floor(y), nv - 1)
                b = y - l
                pr = prof(m)*(1.0 - a) + prof(m + 1)*a
                am = amp(m)*(1.0 - a) + amp(m + 1)*a
                pa = g(m, l)*(1.0 - a)*(1.0 - b) + g(m + 1, l)*a*(1.0 - b) &
                    + g(m, l + 1)*(1.0 - a)*b + g(m + 1, l + 1)*a*b
                ! Smooth positive envelope (softplus) of 1 + plunge*pa
                env = 4.0*(1.0 + plunge*pa)
                env = (log(1.0 + exp(-abs(env))) + max(env, 0.0))/4.0
                r(j, k) = am*(trough + env*(pr - trough))
            end do
        end do
        !$omp end parallel do

    end function fold_surface

    !
    !> Two-piece (skewed) bump: the scale is sigma*(1 - skew) for x < mu and
    !> sigma*(1 + skew) for x >= mu (|skew| < 1); kind = 'cauchy' gives
    !> 1/(1 + q^2) and otherwise exp(-q^2/2), as libflit's gaussian and cauchy
    !
    function skewed_bump_1d(x, mu, sigma, skew, kind) result(f)

        real, dimension(:), intent(in) :: x
        real, intent(in) :: mu, sigma, skew
        character(len=*), intent(in) :: kind
        real, allocatable, dimension(:) :: f

        real, allocatable, dimension(:) :: q

        q = (x - mu)/merge(sigma*(1.0 - skew), sigma*(1.0 + skew), x < mu)
        if (kind == 'cauchy') then
            f = 1.0/(1.0 + q**2)
        else
            f = exp(-0.5*q**2)
        end if

    end function skewed_bump_1d

    !
    !> Rotated two-piece bump on the grid (x, y), with the rotation of libflit's
    !> 2D gaussian and cauchy: g1 = dx*cos(theta) - dy*sin(theta) and
    !> g2 = dx*sin(theta) + dy*cos(theta); the scale along g1 is
    !> sigma(1)*(1 - skew) for g1 < 0 and sigma(1)*(1 + skew) for g1 >= 0
    !
    function skewed_bump_2d(x, y, mu, sigma, theta, skew, kind) result(f)

        real, dimension(:), intent(in) :: x, y
        real, dimension(1:2), intent(in) :: mu, sigma
        real, intent(in) :: theta, skew
        character(len=*), intent(in) :: kind
        real, allocatable, dimension(:, :) :: f

        real :: g1, g2, s1
        integer :: i, j
        logical :: yn_cauchy

        yn_cauchy = kind == 'cauchy'
        f = zeros(size(x), size(y))
        !$omp parallel do private(i, j, g1, g2, s1)
        do j = 1, size(y)
            do i = 1, size(x)
                g1 = (x(i) - mu(1))*cos(theta) - (y(j) - mu(2))*sin(theta)
                g2 = (x(i) - mu(1))*sin(theta) + (y(j) - mu(2))*cos(theta)
                if (g1 < 0) then
                    s1 = sigma(1)*(1.0 - skew)
                else
                    s1 = sigma(1)*(1.0 + skew)
                end if
                if (yn_cauchy) then
                    f(i, j) = (1.0 + (g1/s1)**2 + (g2/sigma(2))**2)**(-1.5)
                else
                    f(i, j) = exp(-0.5*((g1/s1)**2 + (g2/sigma(2))**2))
                end if
            end do
        end do
        !$omp end parallel do

    end function skewed_bump_2d

    !
    !> Cosine of the layer dip from the gradient of the relative geological time
    !> (RGT). Voxels flagged in mask (e.g., faults, salt, karst) and their
    !> immediate neighbors do not contribute, because the RGT jumps there and its
    !> gradient does not measure the layer dip; the remaining gradients are
    !> averaged with a normalized Gaussian of standard deviation sigma.
    !
    function dip_cosine_2d(rgt, mask, sigma) result(c)

        real, dimension(:, :), intent(in) :: rgt
        logical, dimension(:, :), intent(in) :: mask
        real, intent(in) :: sigma
        real, allocatable, dimension(:, :) :: c

        real, allocatable, dimension(:, :) :: g1, g2, w, s
        integer :: n1, n2, i, j

        n1 = size(rgt, 1)
        n2 = size(rgt, 2)

        g1 = zeros(n1, n2)
        g2 = zeros(n1, n2)
        !$omp parallel do private(i, j)
        do j = 1, n2
            do i = 1, n1
                g1(i, j) = (rgt(min(i + 1, n1), j) - rgt(max(i - 1, 1), j))/(min(i + 1, n1) - max(i - 1, 1))
                g2(i, j) = (rgt(i, min(j + 1, n2)) - rgt(i, max(j - 1, 1)))/(min(j + 1, n2) - max(j - 1, 1))
            end do
        end do
        !$omp end parallel do

        w = gauss_smooth_2d(merge(1.0, 0.0, mask), [1.0, 1.0])
        w = merge(0.0, 1.0, w > 0.05)
        s = gauss_smooth_2d(w, [sigma, sigma]) + float_tiny
        g1 = gauss_smooth_2d(g1*w, [sigma, sigma])/s
        g2 = gauss_smooth_2d(g2*w, [sigma, sigma])/s
        c = abs(g1)/(sqrt(g1**2 + g2**2) + float_tiny)

    end function dip_cosine_2d

    function dip_cosine_3d(rgt, mask, sigma) result(c)

        real, dimension(:, :, :), intent(in) :: rgt
        logical, dimension(:, :, :), intent(in) :: mask
        real, intent(in) :: sigma
        real, allocatable, dimension(:, :, :) :: c

        real, allocatable, dimension(:, :, :) :: g1, g2, g3, w, s
        integer :: n1, n2, n3, i, j, k

        n1 = size(rgt, 1)
        n2 = size(rgt, 2)
        n3 = size(rgt, 3)

        g1 = zeros(n1, n2, n3)
        g2 = zeros(n1, n2, n3)
        g3 = zeros(n1, n2, n3)
        !$omp parallel do private(i, j, k) collapse(2)
        do k = 1, n3
            do j = 1, n2
                do i = 1, n1
                    g1(i, j, k) = (rgt(min(i + 1, n1), j, k) - rgt(max(i - 1, 1), j, k))/(min(i + 1, n1) - max(i - 1, 1))
                    g2(i, j, k) = (rgt(i, min(j + 1, n2), k) - rgt(i, max(j - 1, 1), k))/(min(j + 1, n2) - max(j - 1, 1))
                    g3(i, j, k) = (rgt(i, j, min(k + 1, n3)) - rgt(i, j, max(k - 1, 1)))/(min(k + 1, n3) - max(k - 1, 1))
                end do
            end do
        end do
        !$omp end parallel do

        w = gauss_smooth_3d(merge(1.0, 0.0, mask), [1.0, 1.0, 1.0])
        w = merge(0.0, 1.0, w > 0.05)
        s = gauss_smooth_3d(w, [sigma, sigma, sigma]) + float_tiny
        g1 = gauss_smooth_3d(g1*w, [sigma, sigma, sigma])/s
        g2 = gauss_smooth_3d(g2*w, [sigma, sigma, sigma])/s
        g3 = gauss_smooth_3d(g3*w, [sigma, sigma, sigma])/s
        c = abs(g1)/(sqrt(g1**2 + g2**2 + g3**2) + float_tiny)

    end function dip_cosine_3d

    !
    !> Frequency (cycles per sample) of the i-th FFT coefficient of n samples
    !
    pure function fft_freq(i, n) result(f)

        integer, intent(in) :: i, n
        real :: f

        if (i - 1 <= (n - 1)/2) then
            f = (i - 1.0)/n
        else
            f = (i - 1.0 - n)/n
        end if

    end function fft_freq

    !
    !> Isotropic point spread function whose projection onto any direction is the
    !> zero-phase version of the 1D wavelet: by the projection-slice theorem, its
    !> Fourier transform is the 1D amplitude spectrum evaluated at |k|. Convolving
    !> a reflectivity with it gives every reflector the same wavelet whatever its
    !> dip. With n1 = size(wavelet), the PSF is centered at the 0-based sample
    !> (n1/2, n2/2) (integer division), which conv(..., 'same') aligns with the
    !> input, and normalized to unit L2 norm.
    !
    function isotropic_psf_2d(wavelet, n2) result(psf)

        real, dimension(:), intent(in) :: wavelet
        integer, intent(in) :: n2
        real, allocatable, dimension(:, :) :: psf

        complex, allocatable, dimension(:, :) :: pf
        real, allocatable, dimension(:) :: amp
        real :: f1, f2, x, w, c1, c2, pi
        integer :: n1, nh, i, j, m

        pi = real(const_pi)
        n1 = size(wavelet)
        nh = n1/2
        amp = abs(fft(wavelet))
        amp = amp(1:nh + 1)
        c1 = real(n1/2)
        c2 = real(n2/2)

        allocate (pf(1:n1, 1:n2))
        !$omp parallel do private(i, j, f1, f2, x, m, w)
        do j = 1, n2
            f2 = fft_freq(j, n2)
            do i = 1, n1
                f1 = fft_freq(i, n1)
                x = sqrt(f1**2 + f2**2)*n1
                m = floor(x)
                if (m < nh) then
                    w = x - m
                    pf(i, j) = (amp(m + 1)*(1.0 - w) + amp(m + 2)*w) &
                        *exp(cmplx(0.0, -2.0*pi*(f1*c1 + f2*c2)))
                else
                    pf(i, j) = cmplx(0.0, 0.0)
                end if
            end do
        end do
        !$omp end parallel do

        psf = ifft(pf, real=.true.)
        psf = psf/norm2(psf)

    end function isotropic_psf_2d

    function isotropic_psf_3d(wavelet, n2, n3) result(psf)

        real, dimension(:), intent(in) :: wavelet
        integer, intent(in) :: n2, n3
        real, allocatable, dimension(:, :, :) :: psf

        complex, allocatable, dimension(:, :, :) :: pf
        real, allocatable, dimension(:) :: amp
        real :: f1, f2, f3, x, w, c1, c2, c3, pi
        integer :: n1, nh, i, j, k, m

        pi = real(const_pi)
        n1 = size(wavelet)
        nh = n1/2
        amp = abs(fft(wavelet))
        amp = amp(1:nh + 1)
        c1 = real(n1/2)
        c2 = real(n2/2)
        c3 = real(n3/2)

        allocate (pf(1:n1, 1:n2, 1:n3))
        !$omp parallel do private(i, j, k, f1, f2, f3, x, m, w) collapse(2)
        do k = 1, n3
            do j = 1, n2
                f3 = fft_freq(k, n3)
                f2 = fft_freq(j, n2)
                do i = 1, n1
                    f1 = fft_freq(i, n1)
                    x = sqrt(f1**2 + f2**2 + f3**2)*n1
                    m = floor(x)
                    if (m < nh) then
                        w = x - m
                        pf(i, j, k) = (amp(m + 1)*(1.0 - w) + amp(m + 2)*w) &
                            *exp(cmplx(0.0, -2.0*pi*(f1*c1 + f2*c2 + f3*c3)))
                    else
                        pf(i, j, k) = cmplx(0.0, 0.0)
                    end if
                end do
            end do
        end do
        !$omp end parallel do

        psf = ifft(pf, real=.true.)
        psf = psf/norm2(psf)

    end function isotropic_psf_3d

    !
    !> Normalize to zero mean and unit standard deviation
    !
    subroutine unit_std_2d(w)
        real, dimension(:, :), intent(inout) :: w
        w = w - mean(w)
        w = w/(std(w) + float_tiny)
    end subroutine unit_std_2d

    subroutine unit_std_3d(w)
        real, dimension(:, :, :), intent(inout) :: w
        w = w - mean(w)
        w = w/(std(w) + float_tiny)
    end subroutine unit_std_3d

    !
    !> Band-limited background noise: white noise convolved with the PSF
    !> (returned with zero mean and unit standard deviation)
    !
    function noise_band_2d(psf, seed) result(w)

        real, dimension(:, :), intent(in) :: psf
        integer, intent(in) :: seed
        real, allocatable, dimension(:, :) :: w

        w = conv(random(size(psf, 1), size(psf, 2), dist='normal', seed=seed), psf, 'same')
        call unit_std_2d(w)

    end function noise_band_2d

    function noise_band_3d(psf, seed) result(w)

        real, dimension(:, :, :), intent(in) :: psf
        integer, intent(in) :: seed
        real, allocatable, dimension(:, :, :) :: w

        w = conv(random(size(psf, 1), size(psf, 2), size(psf, 3), dist='normal', seed=seed), psf, 'same')
        call unit_std_3d(w)

    end function noise_band_3d

    !
    !> Noise conformable to the layers ("worm" noise): white noise on a grid of
    !> (RGT, lateral position) is smoothed mostly along the lateral axes
    !> (correlation length along, in grid points), mapped back to depth through
    !> the RGT, band-limited by the PSF, and switched on in patches. It appears as
    !> short reflector-parallel segments. Voxels where valid is false (e.g.,
    !> salt, karst, water) get no worm noise.
    !
    function noise_worm_2d(rgt, valid, psf, along, seed) result(w)

        real, dimension(:, :), intent(in) :: rgt, psf
        logical, dimension(:, :), intent(in) :: valid
        real, intent(in) :: along
        integer, intent(in) :: seed
        real, allocatable, dimension(:, :) :: w

        real, allocatable, dimension(:, :) :: base
        real :: x, a, tmin, tmax
        integer :: n1, n2, nt, i, j, m

        n1 = size(rgt, 1)
        n2 = size(rgt, 2)
        nt = 2*n1
        base = gauss_smooth_2d(random(nt, n2, dist='normal', seed=derive_seed(seed, 1)), [0.7, along])
        tmin = minval(rgt)
        tmax = maxval(rgt)

        w = zeros(n1, n2)
        !$omp parallel do private(i, j, x, m, a)
        do j = 1, n2
            do i = 1, n1
                if (valid(i, j)) then
                    x = (rgt(i, j) - tmin)/(tmax - tmin + float_tiny)*(nt - 1.0) + 1.0
                    m = min(max(floor(x), 1), nt - 1)
                    a = x - m
                    w(i, j) = base(m, j)*(1.0 - a) + base(m + 1, j)*a
                end if
            end do
        end do
        !$omp end parallel do

        w = conv(w, psf, 'same')
        call unit_std_2d(w)
        w = w/(1.0 + exp(-2.0*smooth_noise_2d(n1, n2, [25.0, 25.0], derive_seed(seed, 2))))
        call unit_std_2d(w)

    end function noise_worm_2d

    function noise_worm_3d(rgt, valid, psf, along, seed) result(w)

        real, dimension(:, :, :), intent(in) :: rgt, psf
        logical, dimension(:, :, :), intent(in) :: valid
        real, intent(in) :: along
        integer, intent(in) :: seed
        real, allocatable, dimension(:, :, :) :: w

        real, allocatable, dimension(:, :, :) :: base
        real :: x, a, tmin, tmax
        integer :: n1, n2, n3, nt, i, j, k, m

        n1 = size(rgt, 1)
        n2 = size(rgt, 2)
        n3 = size(rgt, 3)
        nt = 2*n1
        base = gauss_smooth_3d(random(nt, n2, n3, dist='normal', seed=derive_seed(seed, 1)), [0.7, along, along])
        tmin = minval(rgt)
        tmax = maxval(rgt)

        w = zeros(n1, n2, n3)
        !$omp parallel do private(i, j, k, x, m, a) collapse(2)
        do k = 1, n3
            do j = 1, n2
                do i = 1, n1
                    if (valid(i, j, k)) then
                        x = (rgt(i, j, k) - tmin)/(tmax - tmin + float_tiny)*(nt - 1.0) + 1.0
                        m = min(max(floor(x), 1), nt - 1)
                        a = x - m
                        w(i, j, k) = base(m, j, k)*(1.0 - a) + base(m + 1, j, k)*a
                    end if
                end do
            end do
        end do
        !$omp end parallel do

        w = conv(w, psf, 'same')
        call unit_std_3d(w)
        w = w/(1.0 + exp(-2.0*smooth_noise_3d(n1, n2, n3, [25.0, 25.0, 25.0], derive_seed(seed, 2))))
        call unit_std_3d(w)

    end function noise_worm_3d

    !
    !> Migration swing noise: white noise restricted to narrow fans (5 degree
    !> half-width) of steep dips, band-limited by the PSF, growing with depth, and
    !> switched on in patches. Each of the nfam dip families, with a dip drawn
    !> from dip_range (degrees), contains both senses of dip, mirrored about the
    !> vertical (in 3D, a random azimuth and its opposite), with random relative
    !> strengths between 0.6 and 1, so the streaks cross the reflectors and each
    !> other in both directions, like the two arms of migration smiles. When dual
    !> is present and .false., all families dip in the same sense instead (in 3D,
    !> they share one random azimuth), so the streaks dip in a single direction.
    !
    function noise_swing_2d(psf, dip_range, nfam, seed, dual) result(w)

        real, dimension(:, :), intent(in) :: psf
        real, dimension(1:2), intent(in) :: dip_range
        integer, intent(in) :: nfam, seed
        logical, intent(in), optional :: dual
        real, allocatable, dimension(:, :) :: w

        complex, allocatable, dimension(:, :) :: wf
        real, allocatable, dimension(:) :: dips, cd, sd, amps
        real :: k1, k2, kn, fan, width
        integer :: n1, n2, i, j, l

        n1 = size(psf, 1)
        n2 = size(psf, 2)
        width = 5.0*const_deg2rad
        dips = random(nfam, range=dip_range, seed=derive_seed(seed, 1))*const_deg2rad
        cd = cos(dips)
        sd = sin(dips)
        amps = random(2*nfam, range=[0.6, 1.0], seed=derive_seed(seed, 2))
        if (present(dual)) then
            if (.not. dual) then
                ! Keep one sense of dip, common to all families
                if (rand(range=[-1.0, 1.0], seed=derive_seed(seed, 6)) >= 0) then
                    amps(2::2) = 0.0
                else
                    amps(1::2) = 0.0
                end if
            end if
        end if

        wf = fft(random(n1, n2, dist='normal', seed=derive_seed(seed, 3)))
        !$omp parallel do private(i, j, l, k1, k2, kn, fan)
        do j = 1, n2
            k2 = fft_freq(j, n2)
            do i = 1, n1
                k1 = fft_freq(i, n1)
                kn = sqrt(k1**2 + k2**2)
                fan = 0.0
                if (kn > 0) then
                    do l = 1, nfam
                        fan = fan + amps(2*l - 1)*exp(-0.5*(acos(min(1.0, abs(k1*cd(l) - k2*sd(l))/kn))/width)**2) &
                            + amps(2*l)*exp(-0.5*(acos(min(1.0, abs(k1*cd(l) + k2*sd(l))/kn))/width)**2)
                    end do
                end if
                wf(i, j) = wf(i, j)*fan
            end do
        end do
        !$omp end parallel do

        w = conv(ifft(wf, real=.true.), psf, 'same')
        do i = 1, n1
            w(i, :) = w(i, :)*(0.25 + 0.75*((i - 1.0)/max(n1 - 1, 1))**1.5)
        end do
        call unit_std_2d(w)
        w = w/(1.0 + exp(-2.0*smooth_noise_2d(n1, n2, [25.0, 25.0], derive_seed(seed, 4))))
        call unit_std_2d(w)

    end function noise_swing_2d

    function noise_swing_3d(psf, dip_range, nfam, seed, dual) result(w)

        real, dimension(:, :, :), intent(in) :: psf
        real, dimension(1:2), intent(in) :: dip_range
        integer, intent(in) :: nfam, seed
        logical, intent(in), optional :: dual
        real, allocatable, dimension(:, :, :) :: w

        complex, allocatable, dimension(:, :, :) :: wf
        real, allocatable, dimension(:) :: dips, azim, n1v, n2v, n3v, amps
        real :: k1, k2, k3, kn, fan, width
        integer :: n1, n2, n3, i, j, k, l

        n1 = size(psf, 1)
        n2 = size(psf, 2)
        n3 = size(psf, 3)
        width = 5.0*const_deg2rad
        dips = random(nfam, range=dip_range, seed=derive_seed(seed, 1))*const_deg2rad
        azim = random(nfam, range=[0.0, 360.0], seed=derive_seed(seed, 2))*const_deg2rad
        if (present(dual)) then
            if (.not. dual) then
                ! One azimuth, common to all families
                azim = azim(1)
            end if
        end if
        ! Unit normals of the dip families; the mirrored sense of each family has
        ! the opposite azimuth, i.e., the opposite sign of the lateral components
        n1v = cos(dips)
        n2v = sin(dips)*cos(azim)
        n3v = sin(dips)*sin(azim)
        amps = random(2*nfam, range=[0.6, 1.0], seed=derive_seed(seed, 5))
        if (present(dual)) then
            if (.not. dual) then
                ! Keep one sense of dip
                amps(2::2) = 0.0
            end if
        end if

        wf = fft(random(n1, n2, n3, dist='normal', seed=derive_seed(seed, 3)))
        !$omp parallel do private(i, j, k, l, k1, k2, k3, kn, fan) collapse(2)
        do k = 1, n3
            do j = 1, n2
                k3 = fft_freq(k, n3)
                k2 = fft_freq(j, n2)
                do i = 1, n1
                    k1 = fft_freq(i, n1)
                    kn = sqrt(k1**2 + k2**2 + k3**2)
                    fan = 0.0
                    if (kn > 0) then
                        do l = 1, nfam
                            fan = fan &
                                + amps(2*l - 1)*exp(-0.5*(acos(min(1.0, abs(k1*n1v(l) - k2*n2v(l) - k3*n3v(l))/kn))/width)**2) &
                                + amps(2*l)*exp(-0.5*(acos(min(1.0, abs(k1*n1v(l) + k2*n2v(l) + k3*n3v(l))/kn))/width)**2)
                        end do
                    end if
                    wf(i, j, k) = wf(i, j, k)*fan
                end do
            end do
        end do
        !$omp end parallel do

        w = conv(ifft(wf, real=.true.), psf, 'same')
        do i = 1, n1
            w(i, :, :) = w(i, :, :)*(0.25 + 0.75*((i - 1.0)/max(n1 - 1, 1))**1.5)
        end do
        call unit_std_3d(w)
        w = w/(1.0 + exp(-2.0*smooth_noise_3d(n1, n2, n3, [25.0, 25.0, 25.0], derive_seed(seed, 4))))
        call unit_std_3d(w)

    end function noise_swing_3d

    !
    !> Smooth multiplicative illumination gain with unit mean; level is the
    !> standard deviation of its logarithm
    !
    function illumination_2d(n1, n2, level, sigma, seed) result(g)

        integer, intent(in) :: n1, n2, seed
        real, intent(in) :: level, sigma
        real, allocatable, dimension(:, :) :: g

        g = exp(level*smooth_noise_2d(n1, n2, [sigma, sigma], seed))
        g = g/mean(g)

    end function illumination_2d

    function illumination_3d(n1, n2, n3, level, sigma, seed) result(g)

        integer, intent(in) :: n1, n2, n3, seed
        real, intent(in) :: level, sigma
        real, allocatable, dimension(:, :, :) :: g

        g = exp(level*smooth_noise_3d(n1, n2, n3, [sigma, sigma, sigma], seed))
        g = g/mean(g)

    end function illumination_3d

    !
    !> Trace-to-trace jitter: every trace is shifted vertically by a smooth static
    !> (maximum shift(1) grid points) plus a random static (standard deviation
    !> shift(2) grid points), and scaled by a random gain 1 + gain*N(0, 1)
    !
    subroutine trace_jitter_2d(w, shift, gain, seed)

        real, dimension(:, :), intent(inout) :: w
        real, dimension(1:2), intent(in) :: shift
        real, intent(in) :: gain
        integer, intent(in) :: seed

        real, allocatable, dimension(:) :: s, g, x
        integer :: n1, n2, j

        n1 = size(w, 1)
        n2 = size(w, 2)
        s = smooth_noise_1d(n2, 3.0, derive_seed(seed, 1))
        s = shift(1)*s/(maxval(abs(s)) + float_tiny) + shift(2)*random(n2, dist='normal', seed=derive_seed(seed, 2))
        g = 1.0 + gain*random(n2, dist='normal', seed=derive_seed(seed, 3))
        x = linspace(0.0, n1 - 1.0, n1)

        !$omp parallel do private(j)
        do j = 1, n2
            w(:, j) = ginterp(x, w(:, j), clip(x + s(j), 0.0, n1 - 1.0), 'cubic')*g(j)
        end do
        !$omp end parallel do

    end subroutine trace_jitter_2d

    subroutine trace_jitter_3d(w, shift, gain, seed)

        real, dimension(:, :, :), intent(inout) :: w
        real, dimension(1:2), intent(in) :: shift
        real, intent(in) :: gain
        integer, intent(in) :: seed

        real, allocatable, dimension(:, :) :: s, g
        real, allocatable, dimension(:) :: x
        integer :: n1, n2, n3, j, k

        n1 = size(w, 1)
        n2 = size(w, 2)
        n3 = size(w, 3)
        s = smooth_noise_2d(n2, n3, [3.0, 3.0], derive_seed(seed, 1))
        s = shift(1)*s/(maxval(abs(s)) + float_tiny) + shift(2)*random(n2, n3, dist='normal', seed=derive_seed(seed, 2))
        g = 1.0 + gain*random(n2, n3, dist='normal', seed=derive_seed(seed, 3))
        x = linspace(0.0, n1 - 1.0, n1)

        !$omp parallel do private(j, k) collapse(2)
        do k = 1, n3
            do j = 1, n2
                w(:, j, k) = ginterp(x, w(:, j, k), clip(x + s(j, k), 0.0, n1 - 1.0), 'cubic')*g(j, k)
            end do
        end do
        !$omp end parallel do

    end subroutine trace_jitter_3d

    !==============================================================================================
    ! Building blocks used directly by rgm2_curved and rgm3_curved
    !==============================================================================================

    !
    !> Fold-train reflector of length n for rgm2_curved (refl_shape = 'fold'). The
    !> mean wavelength, the crest asymmetry and the vergence magnitude are drawn
    !> from their ranges (lambda = [0, 0] means [0.2, 0.5]*nref), and vsign is the
    !> sense of vergence.
    !
    function fold_reflector_1d(n, nref, lambda, crest, vergence, vsign, lambda_drift, amp_drift, mode, seed) result(r)

        integer, intent(in) :: n, nref, seed
        real, dimension(1:2), intent(in) :: lambda, crest, vergence
        real, intent(in) :: vsign, lambda_drift, amp_drift
        character(len=*), intent(in) :: mode
        real, allocatable, dimension(:) :: r

        real, allocatable, dimension(:) :: amp
        real, dimension(1:2) :: lam
        real :: verg, crst

        if (maxval(lambda) == 0) then
            lam = [0.2, 0.5]*nref
        else
            lam = lambda
        end if
        verg = min(max(rand(range=vergence, seed=derive_seed(seed, 1)), 0.0), 0.9)*vsign
        crst = min(max(rand(range=crest, seed=derive_seed(seed, 2)), -0.25), 0.25)
        call fold_train(n, rand(range=lam, seed=derive_seed(seed, 3)), verg, crst, &
            lambda_drift, amp_drift, mode, derive_seed(seed, 4), r, amp)
        r = r*amp

    end function fold_reflector_1d

    !
    !> Fold-train surface of size (n2, n3) for rgm3_curved (refl_shape = 'fold'),
    !> with fold axes striking at strike (radians); see fold_reflector_1d and
    !> fold_surface
    !
    function fold_reflector_2d(n2, n3, nref, lambda, crest, vergence, vsign, strike, &
            lambda_drift, amp_drift, mode, plunge, plunge_length, wobble, seed) result(r)

        integer, intent(in) :: n2, n3, nref, seed
        real, dimension(1:2), intent(in) :: lambda, crest, vergence
        real, intent(in) :: vsign, strike, lambda_drift, amp_drift, plunge, plunge_length, wobble
        character(len=*), intent(in) :: mode
        real, allocatable, dimension(:, :) :: r

        real, dimension(1:2) :: lam
        real :: verg, crst

        if (maxval(lambda) == 0) then
            lam = [0.2, 0.5]*nref
        else
            lam = lambda
        end if
        verg = min(max(rand(range=vergence, seed=derive_seed(seed, 1)), 0.0), 0.9)*vsign
        crst = min(max(rand(range=crest, seed=derive_seed(seed, 2)), -0.25), 0.25)
        r = fold_surface(n2, n3, rand(range=lam, seed=derive_seed(seed, 3)), strike, verg, crst, &
            lambda_drift, amp_drift, mode, plunge, plunge_length, wobble, derive_seed(seed, 4))

    end function fold_reflector_2d

    !
    !> Signed skews of ng bumps: magnitudes drawn from skew_range (clipped to
    !> [0, 0.95]) and random signs, common to all bumps if yn_common
    !
    function bump_skews(ng, skew_range, yn_common, seed) result(skew)

        integer, intent(in) :: ng, seed
        real, dimension(1:2), intent(in) :: skew_range
        logical, intent(in) :: yn_common
        real, allocatable, dimension(:) :: skew

        real, allocatable, dimension(:) :: sgn

        skew = min(max(random(ng, range=skew_range, seed=derive_seed(seed, 1)), 0.0), 0.95)
        sgn = sign(ones(ng), random(ng, range=[-1.0, 1.0], seed=derive_seed(seed, 2)))
        if (yn_common) then
            sgn = sgn(1)
        end if
        skew = skew*sgn

    end function bump_skews

    !
    !> Reflector of length n made of Gaussian or Cauchy bumps (kind) centered at
    !> mu with scales sigma and heights height, skewed when skew_range > 0 (see
    !> bump_skews and skewed_bump_1d), plus a smooth random background of
    !> relative amplitude background (correlation length 0.15*nref)
    !
    function bump_reflector_1d(n, mu, sigma, height, kind, skew_range, yn_common, background, nref, &
            skew_seed, background_seed) result(r)

        integer, intent(in) :: n, nref, skew_seed, background_seed
        real, dimension(:), intent(in) :: mu, sigma, height
        character(len=*), intent(in) :: kind
        real, dimension(1:2), intent(in) :: skew_range
        logical, intent(in) :: yn_common
        real, intent(in) :: background
        real, allocatable, dimension(:) :: r

        real, allocatable, dimension(:) :: x, skew, sb
        integer :: i

        x = linspace(0.0, n - 1.0, n)
        r = zeros(n)
        if (maxval(skew_range) > 0) then
            skew = bump_skews(size(mu), skew_range, yn_common, skew_seed)
            do i = 1, size(mu)
                r = r + rescale(skewed_bump_1d(x, mu(i), sigma(i), skew(i), kind), [0.0, height(i)])
            end do
        else
            do i = 1, size(mu)
                if (kind == 'cauchy') then
                    r = r + rescale(cauchy(x, mu(i), sigma(i)), [0.0, height(i)])
                else
                    r = r + rescale(gaussian(x, mu(i), sigma(i)), [0.0, height(i)])
                end if
            end do
        end if
        if (background > 0) then
            sb = smooth_noise_1d(n, 0.15*nref, background_seed)
            r = r + background*maxval(abs(r))*sb/maxval(abs(sb))
        end if

    end function bump_reflector_1d

    !
    !> Surface of size (n2, n3) made of rotated Gaussian or Cauchy bumps (kind);
    !> the skew acts along the rotated x2 axis of each bump (see bump_reflector_1d
    !> and skewed_bump_2d)
    !
    function bump_reflector_2d(n2, n3, mu2, mu3, sigma2, sigma3, theta, height, kind, skew_range, yn_common, &
            background, nref2, nref3, skew_seed, background_seed) result(r)

        integer, intent(in) :: n2, n3, nref2, nref3, skew_seed, background_seed
        real, dimension(:), intent(in) :: mu2, mu3, sigma2, sigma3, theta, height
        character(len=*), intent(in) :: kind
        real, dimension(1:2), intent(in) :: skew_range
        logical, intent(in) :: yn_common
        real, intent(in) :: background
        real, allocatable, dimension(:, :) :: r

        real, allocatable, dimension(:) :: x, y, skew
        real, allocatable, dimension(:, :) :: sb
        integer :: i

        x = linspace(0.0, n2 - 1.0, n2)
        y = linspace(0.0, n3 - 1.0, n3)
        r = zeros(n2, n3)
        if (maxval(skew_range) > 0) then
            skew = bump_skews(size(mu2), skew_range, yn_common, skew_seed)
            do i = 1, size(mu2)
                r = r + rescale(skewed_bump_2d(x, y, [mu2(i), mu3(i)], [sigma2(i), sigma3(i)], theta(i), skew(i), kind), &
                    [0.0, height(i)])
            end do
        else
            do i = 1, size(mu2)
                if (kind == 'cauchy') then
                    r = r + rescale(cauchy(x, y, [mu2(i), mu3(i)], [sigma2(i), sigma3(i)], theta(i)), [0.0, height(i)])
                else
                    r = r + rescale(gaussian(x, y, [mu2(i), mu3(i)], [sigma2(i), sigma3(i)], theta(i)), [0.0, height(i)])
                end if
            end do
        end if
        if (background > 0) then
            sb = smooth_noise_2d(n2, n3, [0.15*nref2, 0.15*nref3], background_seed)
            r = r + background*maxval(abs(r))*sb/maxval(abs(sb))
        end if

    end function bump_reflector_2d

    !
    !> Voxels that do not belong to continuous layers: salt, karst and, if
    !> yn_with_fault, faults; an array is used only if its flag is set, it is
    !> allocated, and its size is (n1, n2)
    !
    function nonlayer_mask_2d(n1, n2, yn_with_fault, yn_fault, fault, yn_salt, salt, yn_karst, karst) result(m)

        integer, intent(in) :: n1, n2
        logical, intent(in) :: yn_with_fault, yn_fault, yn_salt, yn_karst
        real, allocatable, dimension(:, :), intent(in) :: fault, salt, karst
        logical, allocatable, dimension(:, :) :: m

        m = falses(n1, n2)
        if (yn_with_fault .and. yn_fault .and. allocated(fault)) then
            if (size(fault, 1) == n1 .and. size(fault, 2) == n2) then
                m = m .or. fault > 0
            end if
        end if
        if (yn_salt .and. allocated(salt)) then
            if (size(salt, 1) == n1 .and. size(salt, 2) == n2) then
                m = m .or. salt == 1
            end if
        end if
        if (yn_karst .and. allocated(karst)) then
            if (size(karst, 1) == n1 .and. size(karst, 2) == n2) then
                m = m .or. karst == 1
            end if
        end if

    end function nonlayer_mask_2d

    function nonlayer_mask_3d(n1, n2, n3, yn_with_fault, yn_fault, fault, yn_salt, salt, yn_karst, karst) result(m)

        integer, intent(in) :: n1, n2, n3
        logical, intent(in) :: yn_with_fault, yn_fault, yn_salt, yn_karst
        real, allocatable, dimension(:, :, :), intent(in) :: fault, salt, karst
        logical, allocatable, dimension(:, :, :) :: m

        m = falses(n1, n2, n3)
        if (yn_with_fault .and. yn_fault .and. allocated(fault)) then
            if (size(fault, 1) == n1 .and. size(fault, 2) == n2 .and. size(fault, 3) == n3) then
                m = m .or. fault > 0
            end if
        end if
        if (yn_salt .and. allocated(salt)) then
            if (size(salt, 1) == n1 .and. size(salt, 2) == n2 .and. size(salt, 3) == n3) then
                m = m .or. salt == 1
            end if
        end if
        if (yn_karst .and. allocated(karst)) then
            if (size(karst, 1) == n1 .and. size(karst, 2) == n2 .and. size(karst, 3) == n3) then
                m = m .or. karst == 1
            end if
        end if

    end function nonlayer_mask_3d

    !
    !> Cosine of the layer dip from the RGT, with faults, salt and karst excluded
    !> (see dip_cosine_2d); ones if the RGT is not allocated
    !
    function layer_dip_cosine_2d(n1, n2, rgt, yn_fault, fault, yn_salt, salt, yn_karst, karst) result(c)

        integer, intent(in) :: n1, n2
        real, allocatable, dimension(:, :), intent(in) :: rgt, fault, salt, karst
        logical, intent(in) :: yn_fault, yn_salt, yn_karst
        real, allocatable, dimension(:, :) :: c

        if (allocated(rgt)) then
            c = dip_cosine_2d(rgt, nonlayer_mask_2d(n1, n2, .true., yn_fault, fault, yn_salt, salt, yn_karst, karst), 2.0)
        else
            c = ones(n1, n2)
        end if

    end function layer_dip_cosine_2d

    function layer_dip_cosine_3d(n1, n2, n3, rgt, yn_fault, fault, yn_salt, salt, yn_karst, karst) result(c)

        integer, intent(in) :: n1, n2, n3
        real, allocatable, dimension(:, :, :), intent(in) :: rgt, fault, salt, karst
        logical, intent(in) :: yn_fault, yn_salt, yn_karst
        real, allocatable, dimension(:, :, :) :: c

        if (allocated(rgt)) then
            c = dip_cosine_3d(rgt, nonlayer_mask_3d(n1, n2, n3, .true., yn_fault, fault, yn_salt, salt, yn_karst, karst), 2.0)
        else
            c = ones(n1, n2, n3)
        end if

    end function layer_dip_cosine_3d

    !
    !> Migration-image noise with zero mean and unit standard deviation: worm,
    !> swing and band-limited background components weighted by mix and
    !> band-limited by psf; the worm noise follows the layers through the RGT and
    !> is absent from salt, karst and voxels with zero RGT (e.g., water)
    !
    function migration_noise_2d(psf, rgt, yn_salt, salt, yn_karst, karst, mix, worm_length, &
            swing_dip, swing_nfam, yn_swing_dual, seed) result(w)

        real, dimension(:, :), intent(in) :: psf
        real, allocatable, dimension(:, :), intent(in) :: rgt, salt, karst
        logical, intent(in) :: yn_salt, yn_karst, yn_swing_dual
        real, dimension(1:3), intent(in) :: mix
        real, intent(in) :: worm_length
        real, dimension(1:2), intent(in) :: swing_dip
        integer, intent(in) :: swing_nfam, seed
        real, allocatable, dimension(:, :) :: w

        real, allocatable, dimension(:, :) :: dummy
        real, dimension(1:3) :: wt
        integer :: n1, n2

        n1 = size(psf, 1)
        n2 = size(psf, 2)
        wt = mix/(norm2(mix) + float_tiny)
        w = zeros(n1, n2)
        if (wt(1) > 0 .and. allocated(rgt)) then
            w = w + wt(1)*noise_worm_2d(rgt, (rgt > 0) .and. &
                (.not. nonlayer_mask_2d(n1, n2, .false., .false., dummy, yn_salt, salt, yn_karst, karst)), &
                psf, worm_length, derive_seed(seed, 1))
        end if
        if (wt(2) > 0) then
            w = w + wt(2)*noise_swing_2d(psf, swing_dip, swing_nfam, derive_seed(seed, 2), yn_swing_dual)
        end if
        if (wt(3) > 0) then
            w = w + wt(3)*noise_band_2d(psf, derive_seed(seed, 3))
        end if
        w = (w - mean(w))/(std(w) + float_tiny)

    end function migration_noise_2d

    function migration_noise_3d(psf, rgt, yn_salt, salt, yn_karst, karst, mix, worm_length, &
            swing_dip, swing_nfam, yn_swing_dual, seed) result(w)

        real, dimension(:, :, :), intent(in) :: psf
        real, allocatable, dimension(:, :, :), intent(in) :: rgt, salt, karst
        logical, intent(in) :: yn_salt, yn_karst, yn_swing_dual
        real, dimension(1:3), intent(in) :: mix
        real, intent(in) :: worm_length
        real, dimension(1:2), intent(in) :: swing_dip
        integer, intent(in) :: swing_nfam, seed
        real, allocatable, dimension(:, :, :) :: w

        real, allocatable, dimension(:, :, :) :: dummy
        real, dimension(1:3) :: wt
        integer :: n1, n2, n3

        n1 = size(psf, 1)
        n2 = size(psf, 2)
        n3 = size(psf, 3)
        wt = mix/(norm2(mix) + float_tiny)
        w = zeros(n1, n2, n3)
        if (wt(1) > 0 .and. allocated(rgt)) then
            w = w + wt(1)*noise_worm_3d(rgt, (rgt > 0) .and. &
                (.not. nonlayer_mask_3d(n1, n2, n3, .false., .false., dummy, yn_salt, salt, yn_karst, karst)), &
                psf, worm_length, derive_seed(seed, 1))
        end if
        if (wt(2) > 0) then
            w = w + wt(2)*noise_swing_3d(psf, swing_dip, swing_nfam, derive_seed(seed, 2), yn_swing_dual)
        end if
        if (wt(3) > 0) then
            w = w + wt(3)*noise_band_3d(psf, derive_seed(seed, 3))
        end if
        w = (w - mean(w))/(std(w) + float_tiny)

    end function migration_noise_3d

    !
    !> Blend images made with the PSFs psf_top and psf_bot linearly with depth,
    !> from wtop at the top to wbot at the bottom, with the amplitude of a flat
    !> reflector kept constant
    !
    function depth_blend_2d(wtop, wbot, psf_top, psf_bot) result(w)

        real, dimension(:, :), intent(in) :: wtop, wbot, psf_top, psf_bot
        real, allocatable, dimension(:, :) :: w

        real, allocatable, dimension(:, :) :: wb
        real :: ab, at
        integer :: i, n1

        n1 = size(wtop, 1)
        ab = maxval(abs(sum(psf_bot, dim=2)))
        at = maxval(abs(sum(psf_top, dim=2)))
        wb = wbot*at/(ab + float_tiny)
        w = wtop
        do i = 1, n1
            w(i, :) = w(i, :) + (wb(i, :) - w(i, :))*(i - 1.0)/max(n1 - 1, 1)
        end do

    end function depth_blend_2d

    function depth_blend_3d(wtop, wbot, psf_top, psf_bot) result(w)

        real, dimension(:, :, :), intent(in) :: wtop, wbot, psf_top, psf_bot
        real, allocatable, dimension(:, :, :) :: w

        real, allocatable, dimension(:, :, :) :: wb
        real :: ab, at
        integer :: i, n1

        n1 = size(wtop, 1)
        ab = maxval(abs(sum(sum(psf_bot, dim=3), dim=2)))
        at = maxval(abs(sum(sum(psf_top, dim=3), dim=2)))
        wb = wbot*at/(ab + float_tiny)
        w = wtop
        do i = 1, n1
            w(i, :, :) = w(i, :, :) + (wb(i, :, :) - w(i, :, :))*(i - 1.0)/max(n1 - 1, 1)
        end do

    end function depth_blend_3d

    !
    !> Deallocate an output array unless it was requested (keep)
    !
    subroutine release_unless_2d(keep, a)

        logical, intent(in) :: keep
        real, allocatable, dimension(:, :), intent(inout) :: a

        if (.not. keep .and. allocated(a)) then
            deallocate (a)
        end if

    end subroutine release_unless_2d

    subroutine release_unless_3d(keep, a)

        logical, intent(in) :: keep
        real, allocatable, dimension(:, :, :), intent(inout) :: a

        if (.not. keep .and. allocated(a)) then
            deallocate (a)
        end if

    end subroutine release_unless_3d

end module geological_model_realism
