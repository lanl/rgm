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


module geological_model_utility

    use libflit

    implicit none

    interface noise_wavenumber
        module procedure :: noise_wavenumber_2d
        module procedure :: noise_wavenumber_3d
    end interface

    interface noise_random_mix
        module procedure :: noise_random_mix_2d
        module procedure :: noise_random_mix_3d
    end interface

    interface random_mask_smooth
        module procedure :: random_mask_smooth_2d
        module procedure :: random_mask_smooth_3d
    end interface

    ! Fractal noise (multi-octave Perlin noise)
    ! See https://we.copernicus.org/articles/22/1/2022/
    type fractal_noise_1d
        ! Number of grid points
        integer :: n1
        ! Number of periods
        integer :: periods1 = 5
        ! Randomness seed
        integer :: seed = -1
        ! Number of octaves
        integer :: octaves = 5
        ! Weight for decaying as persistence^i P_i(x), in the range (0, 1]
        real :: persistence = 0.5
        ! In a fractal, fills space, higher lacunarity indicating more or larger gaps.
        ! In the context of Perlin noise the lacunarity parameter is a multiplier that
        ! determines the rate of change in the number of periods between successive waves.
        ! For example, when lacunarity = 2, then the number of periods doubles for
        ! each successive wave (Fig. 1a to c), or when lacunarity = 4,
        ! then the number of periods quadruples for each successive wave.
        real :: lacunarity = 2.0
    contains
        procedure, public :: generate => generate_fractal_noise_1d
    end type

    type fractal_noise_2d
        integer :: n1, n2
        integer :: periods1 = 5
        integer :: periods2 = 5
        integer :: seed = -1
        integer :: octaves = 5
        real :: persistence = 0.5
        real :: lacunarity = 2.0
    contains
        procedure, public :: generate => generate_fractal_noise_2d
    end type

    type fractal_noise_3d
        integer :: n1, n2, n3
        integer :: periods1 = 5
        integer :: periods2 = 5
        integer :: periods3 = 5
        integer :: seed = -1
        integer :: octaves = 5
        real :: persistence = 0.5
        real :: lacunarity = 2.0
    contains
        procedure, public :: generate => generate_fractal_noise_3d
    end type

    private
    public :: sinc_wavelet
    public :: gaussian_wavelet
    public :: gaussian_deriv_wavelet
    public :: gaussian_wavelet_hdur
    public :: ricker_wavelet
    public :: ricker_deriv_wavelet
    public :: ormsby_wavelet
    public :: noise_random_mix
    public :: noise_wavenumber
    public :: fractal_noise_1d
    public :: fractal_noise_2d
    public :: fractal_noise_3d
    public :: elastic_reflection_coefs
    public :: random_mask_smooth
    public :: random_circular

contains

    function random_circular(n, r, dr, alpha, smooth, seed) result(m)

        integer :: n
        real :: r, dr, alpha, smooth
        integer, optional :: seed
        real, allocatable, dimension(:) :: m

        real, allocatable, dimension(:) :: wk

        if (present(seed)) then
            m = random(n, range=[0.0, 2*real(const_pi)], seed=seed)
        else
            m = random(n, range=[0.0, 2*real(const_pi)])
        end if

        wk = fft_omega(n, 1.0)
        m = ifft(exp(const_i*m)/sqrt(1.0/(alpha*n)**2 + wk**2)*exp(-wk**2*smooth), real=.true.)

        m = m - mean(m)
        m = m/maxval(abs(m))*dr + r

    end function random_circular


    function sinc_wavelet(t, f0) result(wavelet)

        real, intent(in) :: t, f0
        real :: wavelet

        wavelet = sinc(const_pi*f0*t)

    end function sinc_wavelet

    !
    !> @brief Gaussian wavelet
    !
    function gaussian_wavelet(t, f0) result(wavelet)

        real, intent(in) :: t, f0
        real :: wavelet

        wavelet = exp(-(const_pi*f0*t)**2)

    end function gaussian_wavelet


    function gaussian_wavelet_hdur(hdur, t) result(wavelet)

        real :: t, hdur
        real :: wavelet

        real :: a

        a = 1.0/hdur**2
        wavelet = exp(-a*t**2)/(sqrt(const_pi)*hdur)

    end function gaussian_wavelet_hdur

    !
    !> @brief Gaussian first derivative wavelet
    !
    function gaussian_deriv_wavelet(t, f0) result(wavelet)

        real, intent(in) :: t, f0
        real :: wavelet

        wavelet = -t*exp(-(const_pi*f0*t)**2)

    end function gaussian_deriv_wavelet

    !
    !> @brief Gaussian second derivative wavelet, a.k.a. Ricker wavelet
    !
    function ricker_wavelet(t, f0) result(wavelet)

        real, intent(in) :: t, f0
        real :: wavelet

        wavelet = (const_pi*f0*t)**2
        wavelet = (1 - 2.0*wavelet)*exp(-wavelet)

    end function ricker_wavelet

    !
    !> @brief First-order derivative of Ricker wavelet
    !
    function ricker_deriv_wavelet(t, f0) result(wavelet)

        real, intent(in) :: t, f0
        real :: wavelet

        wavelet = (const_pi*f0*t)**2
        wavelet = (2.0*wavelet - 3.0)*exp(-wavelet)*t

    end function ricker_deriv_wavelet

    !
    !> @brief Omsby wavelet
    !
    function ormsby_wavelet(f, t) result(wavelet)

        real, dimension(1:4), intent(in) :: f
        real, intent(in) :: t
        real :: wavelet

        real :: f1, f2, f3, f4

        f1 = f(1)
        f2 = f(2)
        f3 = f(3)
        f4 = f(4)

        wavelet = const_pi &
            *(f4**2/(f4 - f3)*sinc(const_pi*f4*t)**2 &
            - f3**2/(f4 - f3)*sinc(const_pi*f3*t)**2 &
            - f2**2/(f2 - f1)*sinc(const_pi*f2*t)**2 &
            + f1**2/(f2 - f1)*sinc(const_pi*f1*t)**2)

    end function ormsby_wavelet

    function noise_wavenumber_2d(w, level, smooth, seed) result(wt)

        real, dimension(:, :) :: w
        real, optional :: level
        real, dimension(1:2), optional :: smooth
        integer, optional :: seed
        real, allocatable, dimension(:, :) :: wt

        real, allocatable, dimension(:, :) :: r
        integer :: n1, n2, sd
        real :: amp, gs(1:2)

        n1 = size(w, 1)
        n2 = size(w, 2)

        if (present(seed)) then
            sd = seed
        else
            sd = -1
        end if

        if (present(level)) then
            amp = level
        else
            amp = 0.5
        end if

        if (present(smooth)) then
            gs = smooth
        else
            gs = [2.0, 2.0]
        end if

        r = random(n1, n2, seed=sd)
        r = gauss_filt(r, gs)
        r = rescale(r, [0.0, 1.0])
        where (r < amp)
            r = 0.0
        end where

        wt = ifft(fft(w)*r, real=.true.) - w

    end function noise_wavenumber_2d

    function noise_wavenumber_3d(w, level, smooth, seed) result(wt)

        real, dimension(:, :, :) :: w
        real, optional :: level
        real, dimension(1:3), optional :: smooth
        integer, optional :: seed
        real, allocatable, dimension(:, :, :) :: wt

        real, allocatable, dimension(:, :, :) :: r
        integer :: n1, n2, n3, sd
        real :: amp, gs(1:3)

        n1 = size(w, 1)
        n2 = size(w, 2)
        n3 = size(w, 3)

        if (present(seed)) then
            sd = seed
        else
            sd = -1
        end if

        if (present(level)) then
            amp = level
        else
            amp = 0.5
        end if

        if (present(smooth)) then
            gs = smooth
        else
            gs = [2.0, 2.0, 2.0]
        end if

        r = random(n1, n2, n3, seed=sd)
        r = gauss_filt(r, gs)
        r = rescale(r, [0.0, 1.0])
        where (r < amp)
            r = 0.0
        end where

        wt = ifft(fft(w)*r, real=.true.) - w

    end function noise_wavenumber_3d

    function noise_random_mix_2d(n1, n2, level, smooth, seed) result(w)

        integer :: n1, n2
        real, dimension(:), optional :: level, smooth
        integer, optional :: seed
        real, allocatable, dimension(:, :) :: w

        integer :: i, sd
        real, allocatable, dimension(:) :: amp, gs
        real, allocatable, dimension(:, :) :: r

        if (present(level)) then
            amp = level
        else
            amp = [9.0, 3.0, 1.0]
        end if

        if (present(smooth)) then
            gs = smooth
        else
            gs = [9.0, 3.0, 0.0]
        end if

        if (present(seed)) then
            sd = seed
        else
            sd = -1
        end if

        w = zeros(n1, n2)
        do i = 1, size(amp)

            r = gauss_filt(random(n1, n2, seed=sd*i, dist='normal'), [gs(i), gs(i)])
            r = r - mean(r)
            w = w + r/maxval(abs(r))*amp(i)

        end do

        w = w/maxval(abs(w))

    end function noise_random_mix_2d

    function noise_random_mix_3d(n1, n2, n3, level, smooth, seed) result(w)

        integer :: n1, n2, n3
        real, dimension(:), optional :: level, smooth
        integer, optional :: seed
        real, allocatable, dimension(:, :, :) :: w

        integer :: i, sd
        real, allocatable, dimension(:) :: amp, gs
        real, allocatable, dimension(:, :, :) :: r

        if (present(level)) then
            amp = level
        else
            amp = [9.0, 3.0, 1.0]
        end if

        if (present(smooth)) then
            gs = smooth
        else
            gs = [9.0, 3.0, 0.0]
        end if

        if (present(seed)) then
            sd = seed
        else
            sd = -1
        end if

        w = zeros(n1, n2, n3)
        do i = 1, size(amp)

            r = gauss_filt(random(n1, n2, n3, seed=sd*i, dist='normal'), [gs(i), gs(i), gs(i)])
            r = r - mean(r)
            w = w + r/maxval(abs(r))*amp(i)

        end do

        w = w/maxval(abs(w))

    end function noise_random_mix_3d

    !
    !! Compute elastic reflection coefficients
    !
    function elastic_reflection_coefs(p, a1, b1, rho1, a2, b2, rho2) result(r)

        real :: p, a1, b1, rho1, a2, b2, rho2
        real, allocatable, dimension(:) :: r

        real :: a, b, c, d, e, f, g, h, dd
        real :: i1, i2, j1, j2
        real :: rpp, rps, rsp, rss

        a = rho2*(1 - 2*b2**2*p**2) - rho1*(1 - 2*b1**2*p**2)
        b = rho2*(1 - 2*b2**2*p**2) + 2*rho1*b1**2*p**2
        c = rho1*(1 - 2*b1**2*p**2) + 2*rho2*b2**2*p**2
        d = 2*(rho2*b2**2 - rho1*b1**2)

        i1 = asin(min(p*a1, 1.0))
        i2 = asin(min(p*a2, 1.0))
        j1 = asin(min(p*b1, 1.0))
        j2 = asin(min(p*b2, 1.0))

        e = b*cos(i1)/a1 + c*cos(i2)/a2
        f = b*cos(j1)/b1 + c*cos(j2)/b2
        g = a - d*cos(i1)/a1*cos(j2)/b2
        h = a - d*cos(i2)/a2*cos(j1)/b1
        dd = e*f + g*h*p**2

        ! To avoid division by zero
        if (dd == 0) then
            dd = 1.0e-10
        end if

        rpp = (((b*cos(i1)/a1) - c*cos(i2)/a2)*f - (a + d*cos(i1)/a1*cos(j2)/b2)*h*p**2)/dd
        rss = -(((b*cos(j1)/b1) - c*cos(j2)/b2)*e - (a + d*cos(i2)/a2*cos(j1)/b1)*g*p**2)/dd
        rps = -2*cos(i1)/b1*(a*b + c*d*cos(i2)/a2*cos(j2)/b2)*p/dd
        rsp = -2*cos(j1)/a1*(a*b + c*d*cos(i2)/a2*cos(j2)/b2)*p/dd

        r = [rpp, rps, rsp, rss]

    end function elastic_reflection_coefs

    function random_mask_smooth_2d(n1, n2, gs, mask_out, seed) result(r)

        integer :: n1, n2
        real, dimension(1:2) :: gs
        real :: mask_out
        integer, optional :: seed
        real, allocatable, dimension(:, :) :: r

        real :: ratio, bv
        integer :: rs

        if (present(seed)) then
            rs = seed
        else
            rs = -1
        end if

        r = random(n1, n2, dist='normal', seed=rs)
        r = gauss_filt(r, gs)
        r = rescale(r, [0.0, 1.0])
        ratio = 1.0
        bv = 1.0
        do while(ratio > mask_out)
            ratio = count(r <= bv)*1.0/(n1*n2)
            bv = bv - 0.001
        end do
        r = binarize(r, bv, [0.0, 1.0])

    end function random_mask_smooth_2d

    function random_mask_smooth_3d(n1, n2, n3, gs, mask_out, seed) result(r)

        integer :: n1, n2, n3
        real, dimension(1:3) :: gs
        real :: mask_out
        integer, optional :: seed
        real, allocatable, dimension(:, :, :) :: r

        real :: ratio, bv
        integer :: rs

        if (present(seed)) then
            rs = seed
        else
            rs = -1
        end if

        r = random(n1, n2, n3, dist='normal', seed=rs)
        r = gauss_filt(r, gs)
        r = rescale(r, [0.0, 1.0])
        ratio = 1.0
        bv = 1.0
        do while(ratio > mask_out)
            ratio = count(r <= bv)*1.0/(n1*n2*n3)
            bv = bv - 0.001
        end do
        r = binarize(r, bv, [0.0, 1.0])

    end function random_mask_smooth_3d

    !
    ! Fade function (cubic polynomial for smooth interpolation)
    !
    pure function fade(t) result(f)

        real, intent(in) :: t
        real :: f

        f = 6.0*t**5 - 15.0*t**4 + 10.0*t**3

    end function fade

    !
    ! Linear interpolation
    !
    pure function lerp(a, b, t) result(s)

        real, intent(in) :: a, b, t
        real :: s

        s = a + t*(b - a)

    end function lerp

    !
    ! Generate 1D Perlin noise
    !
    function perlin_1d(periodx, nx, seed) result(noise)

        real, intent(in) :: periodx
        integer, intent(in) :: nx, seed
        real, allocatable, dimension(:) :: noise
        real, allocatable, dimension(:) :: x
        integer :: x0, x1, i
        real, allocatable, dimension(:) :: gd
        real :: sx, g0, g1, dx

        noise = zeros(nx)
        x = regspace(0.0, 1.0, nx - 1.0)/nx*periodx

        gd = random_permute(binarize(random(ceiling(periodx) + 1, seed=seed), 0.5, [-1.0, 1.0]), seed=seed + 1)

        !$omp parallel do private(i, x0, x1, dx, g0, g1, sx)
        do i = 1, nx

            x0 = floor(x(i)) + 1
            x1 = x0 + 1

            dx = x(i) - (x0 - 1.0)

            g0 = gd(x0)*dx
            g1 = gd(x1)*(dx - 1.0)

            sx = fade(dx)

            noise(i) = lerp(g0, g1, sx)

        end do
        !$omp end parallel do

    end function perlin_1d

    !
    ! Generate 1D fractal noise based on Perlin noise
    !
    function generate_fractal_noise_1d(this) result(fractal_noise)

        class(fractal_noise_1d), intent(in) :: this

        integer :: n1, period1, seed, octaves
        real :: persistence, lacunarity
        real, allocatable, dimension(:) :: fractal_noise
        real, allocatable, dimension(:) :: noise
        real :: amplitude, frequency
        integer :: octave, nx, grid_check

        n1 = this%n1
        period1 = this%periods1
        seed = this%seed*20
        octaves = this%octaves
        persistence = this%persistence
        lacunarity = this%lacunarity

        grid_check = floor(lacunarity**(octaves - 1)*period1)
        nx = ceiling(n1*1.0/grid_check)*grid_check

        ! Allocate fractal noise array
        fractal_noise = zeros(nx)

        ! Start with base frequency and amplitude
        frequency = 1.0
        amplitude = 1.0

        ! Loop over octaves
        do octave = 1, octaves

            ! Generate Perlin noise for the current octave
            noise = perlin_1d(period1*frequency, nx, seed + octave)

            ! Add the current octave to the fractal noise
            fractal_noise = fractal_noise + amplitude*noise

            ! Update amplitude and frequency for the next octave
            amplitude = amplitude*persistence
            frequency = frequency*lacunarity

        end do

        ! Normalize fractal noise
        fractal_noise = fractal_noise(1:n1)
        fractal_noise = (fractal_noise - mean(fractal_noise))/std(fractal_noise)

    end function generate_fractal_noise_1d

    !
    ! Generate 2D Perlin noise
    !
    function perlin_2d(periodx, periody, nx, ny, seed) result(noise)

        real, intent(in) :: periodx, periody
        integer, intent(in) :: nx, ny, seed
        real, allocatable, dimension(:, :) :: noise
        integer :: x0, x1, y0, y1, i, j
        real :: sx, sy, g00, g01, g10, g11, dx, dy
        integer :: mx, my
        real, allocatable, dimension(:) :: x, y
        real, allocatable, dimension(:, :, :) :: gd
        real, allocatable, dimension(:, :) :: r

        noise = zeros(nx, ny)
        x = regspace(0.0, 1.0, nx - 1.0)/nx*periodx
        y = regspace(0.0, 1.0, ny - 1.0)/ny*periody

        mx = ceiling(periodx) + 1
        my = ceiling(periody) + 1
        r = random(mx, my, seed=seed)*2*const_pi
        gd = zeros(mx, my, 2)

        !$omp parallel do private(i, j)
        do j = 1, my
            do i = 1, mx
                gd(i, j, 1) = cos(r(i, j))
                gd(i, j, 2) = sin(r(i, j))
            end do
        end do
        !$omp end parallel do

        !$omp parallel do private(i, j, x0, x1, y0, y1, dx, dy, g00, g01, g10, g11, sx, sy)
        do j = 1, ny
            do i = 1, nx

                x0 = floor(x(i)) + 1
                x1 = x0 + 1
                y0 = floor(y(j)) + 1
                y1 = y0 + 1

                dx = x(i) - (x0 - 1.0)
                dy = y(j) - (y0 - 1.0)

                g00 = gd(x0, y0, 1)*dx         + gd(x0, y0, 2)*dy
                g01 = gd(x0, y1, 1)*dx         + gd(x0, y1, 2)*(dy - 1.0)
                g10 = gd(x1, y0, 1)*(dx - 1.0) + gd(x1, y0, 2)*dy
                g11 = gd(x1, y1, 1)*(dx - 1.0) + gd(x1, y1, 2)*(dy - 1.0)

                sx = fade(dx)
                sy = fade(dy)

                noise(i, j) = lerp(lerp(g00, g10, sx), lerp(g01, g11, sx), sy)

            end do
        end do
        !$omp end parallel do

    end function perlin_2d

    !
    ! Generate 2D fractal noise based on Perlin noise
    !
    function generate_fractal_noise_2d(this) result(fractal_noise)

        class(fractal_noise_2d), intent(in) :: this

        integer :: n1, n2, periods1, periods2, seed, octaves
        real :: persistence, lacunarity
        real, allocatable, dimension(:, :) :: fractal_noise
        real, allocatable, dimension(:, :) :: noise
        real :: amplitude, frequency
        integer :: octave
        integer :: grid_check, nx, ny

        n1 = this%n1
        n2 = this%n2
        periods1 = this%periods1
        periods2 = this%periods2
        seed = this%seed*20
        octaves = this%octaves
        persistence = this%persistence
        lacunarity = this%lacunarity

        grid_check = floor(lacunarity**(octaves - 1)*periods1)
        nx = ceiling(n1*1.0/grid_check)*grid_check

        grid_check = floor(lacunarity**(octaves - 1)*periods2)
        ny = ceiling(n2*1.0/grid_check)*grid_check

        ! Initialize the fractal noise array
        fractal_noise = zeros(nx, ny)

        ! Start with base frequency and amplitude
        frequency = 1.0
        amplitude = 1.0

        ! Loop over octaves
        do octave = 1, octaves

            ! Generate Perlin noise for the current octave
            noise = perlin_2d(periods1*frequency, periods2*frequency, nx, ny, seed + octave)

            ! Add the current octave to the fractal noise
            fractal_noise = fractal_noise + amplitude*noise

            ! Update amplitude and frequency for the next octave
            amplitude = amplitude*persistence
            frequency = frequency*lacunarity

        end do

        ! Normalize fractal noise to [0, 1]
        fractal_noise = fractal_noise(1:n1, 1:n2)
        fractal_noise = (fractal_noise - mean(fractal_noise))/std(fractal_noise)

    end function generate_fractal_noise_2d

    !
    ! Generate 3D Perlin noise
    !
    function perlin_3d(periodx, periody, periodz, nx, ny, nz, seed) result(noise)

        real, intent(in) :: periodx, periody, periodz
        integer, intent(in) :: nx, ny, nz, seed
        real, allocatable, dimension(:, :, :) :: noise

        integer :: x0, x1, y0, y1, z0, z1, i, j, k
        real :: sx, sy, sz, g000, g001, g010, g011, g100, g101, g110, g111, dx, dy, dz
        integer :: mx, my, mz
        real, allocatable, dimension(:) :: x, y, z
        real, allocatable, dimension(:, :, :, :) :: gd
        real, allocatable, dimension(:, :, :) :: r, t

        noise = zeros(nx, ny, nz)
        x = regspace(0.0, 1.0, nx - 1.0)/nx*periodx
        y = regspace(0.0, 1.0, ny - 1.0)/ny*periody
        z = regspace(0.0, 1.0, nz - 1.0)/nz*periodz

        mx = ceiling(periodx) + 1
        my = ceiling(periody) + 1
        mz = ceiling(periodz) + 1
        r = random(mx, my, mz, seed=seed)*const_pi
        t = random(mx, my, mz, seed=seed*2)*2*const_pi
        gd = zeros(mx, my, mz, 3)

        !$omp parallel do private(i, j, k)
        do k = 1, mz
            do j = 1, my
                do i = 1, mx
                    gd(i, j, k, 1) = sin(r(i, j, k))*cos(t(i, j, k))
                    gd(i, j, k, 2) = sin(r(i, j, k))*sin(t(i, j, k))
                    gd(i, j, k, 3) = cos(r(i, j, k))
                end do
            end do
        end do
        !$omp end parallel do

        !$omp parallel do private(i, j, k, x0, x1, y0, y1, z0, z1, dx, dy, dz, &
            !$omp   g000, g001, g010, g011, g100, g101, g110, g111, sx, sy, sz)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx

                    x0 = floor(x(i)) + 1
                    x1 = x0 + 1
                    y0 = floor(y(j)) + 1
                    y1 = y0 + 1
                    z0 = floor(z(k)) + 1
                    z1 = z0 + 1

                    dx = x(i) - (x0 - 1.0)
                    dy = y(j) - (y0 - 1.0)
                    dz = z(k) - (z0 - 1.0)

                    g000 = gd(x0, y0, z0, 1)*dx         + gd(x0, y0, z0, 2)*dy         + gd(x0, y0, z0, 3)*dz
                    g001 = gd(x0, y0, z1, 1)*dx         + gd(x0, y0, z1, 2)*dy         + gd(x0, y0, z1, 3)*(dz - 1.0)
                    g010 = gd(x0, y1, z0, 1)*dx         + gd(x0, y1, z0, 2)*(dy - 1.0) + gd(x0, y1, z0, 3)*dz
                    g011 = gd(x0, y1, z1, 1)*dx         + gd(x0, y1, z1, 2)*(dy - 1.0) + gd(x0, y1, z1, 3)*(dz - 1.0)
                    g100 = gd(x1, y0, z0, 1)*(dx - 1.0) + gd(x1, y0, z0, 2)*dy         + gd(x1, y0, z0, 3)*dz
                    g101 = gd(x1, y0, z1, 1)*(dx - 1.0) + gd(x1, y0, z1, 2)*dy         + gd(x1, y0, z1, 3)*(dz - 1.0)
                    g110 = gd(x1, y1, z0, 1)*(dx - 1.0) + gd(x1, y1, z0, 2)*(dy - 1.0) + gd(x1, y1, z0, 3)*dz
                    g111 = gd(x1, y1, z1, 1)*(dx - 1.0) + gd(x1, y1, z1, 2)*(dy - 1.0) + gd(x1, y1, z1, 3)*(dz - 1.0)

                    sx = fade(dx)
                    sy = fade(dy)
                    sz = fade(dz)

                    noise(i, j, k) = lerp( &
                        lerp(lerp(g000, g100, sx), lerp(g010, g110, sx), sy), &
                        lerp(lerp(g001, g101, sx), lerp(g011, g111, sx), sy), sz)

                end do
            end do
        end do
        !$omp end parallel do

    end function perlin_3d

    !
    ! Generate 3D fractal noise based on Perlin noise
    !
    function generate_fractal_noise_3d(this) result(fractal_noise)

        class(fractal_noise_3d), intent(in) :: this

        integer :: n1, n2, n3, periods1, periods2, periods3, seed, octaves
        real :: persistence, lacunarity
        real, allocatable, dimension(:, :, :) :: fractal_noise
        real, allocatable, dimension(:, :, :) :: noise
        real :: amplitude, frequency
        integer :: octave
        integer :: grid_check, nx, ny, nz

        n1 = this%n1
        n2 = this%n2
        n3 = this%n3
        periods1 = this%periods1
        periods2 = this%periods2
        periods3 = this%periods3
        seed = this%seed*20
        octaves = this%octaves
        persistence = this%persistence
        lacunarity = this%lacunarity

        grid_check = floor(lacunarity**(octaves - 1)*periods1)
        nx = ceiling(n1*1.0/grid_check)*grid_check

        grid_check = floor(lacunarity**(octaves - 1)*periods2)
        ny = ceiling(n2*1.0/grid_check)*grid_check

        grid_check = floor(lacunarity**(octaves - 1)*periods3)
        nz = ceiling(n3*1.0/grid_check)*grid_check

        ! Initialize the fractal noise array
        fractal_noise = zeros(nx, ny, nz)

        ! Start with base frequency and amplitude
        frequency = 1.0
        amplitude = 1.0

        ! Loop over octaves
        do octave = 1, octaves

            ! Generate Perlin noise for the current octave
            noise = perlin_3d(periods1*frequency, periods2*frequency, periods3*frequency, nx, ny, nz, seed + octave)

            ! Add the current octave to the fractal noise
            fractal_noise = fractal_noise + amplitude*noise

            ! Update amplitude and frequency for the next octave
            amplitude = amplitude*persistence
            frequency = frequency*lacunarity

        end do

        ! Normalize fractal noise to [0, 1]
        fractal_noise = fractal_noise(1:n1, 1:n2, 1:n3)
        fractal_noise = (fractal_noise - mean(fractal_noise))/std(fractal_noise)

    end function generate_fractal_noise_3d

end module geological_model_utility
