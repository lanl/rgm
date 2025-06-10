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

    type fractal_noise_1d
        integer :: n1
        integer :: res1 = 5
        integer :: seed = -1
        integer :: octaves = 5
        real :: persistence = 0.5
        real :: lacunarity = 2.0
    contains
        procedure, public :: generate => generate_fractal_noise_1d
    end type

    type fractal_noise_2d
        integer :: n1, n2
        integer :: res1 = 5
        integer :: res2 = 5
        integer :: seed = -1
        integer :: octaves = 5
        real :: persistence = 0.5
        real :: lacunarity = 2.0
    contains
        procedure, public :: generate => generate_fractal_noise_2d
    end type

    type fractal_noise_3d
        integer :: n1, n2, n3
        integer :: res1 = 5
        integer :: res2 = 5
        integer :: res3 = 5
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

contains

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
    ! Compute 1D gradient dot product
    !
    function gradient_1d(perm, x0, x) result(dot)

        integer, intent(in), dimension(:) :: perm
        integer, intent(in) :: x0
        real, intent(in) :: x
        real :: dot
        integer :: h
        real :: gradient

        ! Compute hash value from permutation table
        h = perm(mod(x0, 256))

        ! Select gradient based on hash (either +1 or -1)
        if (mod(h, 2) == 0) then
            gradient = 1.0
        else
            gradient = -1.0
        end if

        ! Compute dot product: gradient * relative position
        dot = gradient*(x - x0)

    end function gradient_1d

    !
    ! Generate 1D Perlin noise
    !
    function perlin_1d(res, nx, seed) result(noise)

        real, intent(in) :: res
        integer, intent(in) :: nx, seed
        real, allocatable, dimension(:) :: noise
        real, allocatable, dimension(:) :: x
        integer :: perm(0:255), x0, x1, i
        real :: sx, n0, n1

        ! Allocate noise array and coordinate array
        allocate (noise(nx))
        allocate (x(nx))

        ! Generate coordinates scaled by res
        x = [(i/real(nx - 1), i=0, nx - 1)]*res

        ! Initialize permutation table
        perm = random_permute(regspace(0, 1, 255), seed=seed)

        ! Compute Perlin noise
        !$omp parallel do private(i, x0, x1, sx, n0, n1)
        do i = 1, nx

            x0 = floor(x(i))
            x1 = x0 + 1

            sx = fade(x(i) - x0)

            ! Compute gradients and dot products at the two corners
            n0 = gradient_1d(perm, x0, x(i))
            n1 = gradient_1d(perm, x1, x(i))

            ! Interpolate between the two values
            noise(i) = (1.0 - sx)*n0 + sx*n1

        end do
        !$omp end parallel do

    end function perlin_1d

    !
    ! Generate 1D fractal noise based on Perlin noise
    !
    function generate_fractal_noise_1d(this) result(fractal_noise)

        class(fractal_noise_1d), intent(in) :: this

        integer :: n1, res, seed, octaves
        real :: persistence, lacunarity
        real, allocatable, dimension(:) :: fractal_noise
        real, allocatable, dimension(:) :: noise
        real :: amplitude, frequency
        integer :: octave, nx, grid_check

        n1 = this%n1
        res = this%res1
        seed = this%seed*20
        octaves = this%octaves
        persistence = this%persistence
        lacunarity = this%lacunarity

        grid_check = floor(lacunarity**(octaves - 1)*res)
        nx = ceiling(n1*1.0/grid_check)*grid_check

        ! Allocate fractal noise array
        fractal_noise = zeros(nx)

        ! Start with base frequency and amplitude
        frequency = 1.0
        amplitude = 1.0

        ! Loop over octaves
        do octave = 1, octaves

            ! Generate Perlin noise for the current octave
            noise = perlin_1d(res*frequency, nx, seed + octave)

            ! Add the current octave to the fractal noise
            fractal_noise = fractal_noise + amplitude*noise

            ! Update amplitude and frequency for the next octave
            amplitude = amplitude*persistence
            frequency = frequency*lacunarity

        end do

        ! Normalize fractal noise to [0, 1]
        fractal_noise = fractal_noise(1:n1)
        fractal_noise = (fractal_noise - mean(fractal_noise))/std(fractal_noise)

    end function generate_fractal_noise_1d

    !
    ! Compute 2D gradient dot product
    !
    pure function gradient_2d(hash, x, y) result(dot)

        integer, intent(in) :: hash
        real, intent(in) :: x, y
        real :: dot
        real, dimension(4, 2) :: gradients
        integer :: h

        ! Gradients for Perlin noise (limited set of vectors)
        gradients = reshape([ &
            0.0, 1.0, &
            0.0, -1.0, &
            1.0, 0.0, &
            -1.0, 0.0], [4, 2])

        ! Ensure index is between 1 and 4
        h = mod(hash, 4) + 1

        ! Dot product of gradient vector with relative position
        dot = gradients(h, 1)*x + gradients(h, 2)*y

    end function gradient_2d

    !
    ! Generate 2D Perlin noise
    !
    function perlin_2d(res_x, res_y, nx, ny, seed) result(noise)

        real, intent(in) :: res_x, res_y
        integer, intent(in) :: nx, ny, seed
        real, allocatable, dimension(:, :) :: noise
        integer :: x0, x1, y0, y1, i, j
        real :: sx, sy, u, v, n0, n1
        real, allocatable, dimension(:) :: x, y
        integer, allocatable, dimension(:, :) :: permutation_table
        integer, allocatable, dimension(:) :: perm

        ! Allocate noise array and coordinate arrays
        allocate (noise(nx, ny))
        allocate (x(nx))
        allocate (y(ny))

        ! Generate coordinates scaled by res
        x = [(i/real(nx - 1), i=0, nx - 1)]*res_x
        y = [(j/real(ny - 1), j=0, ny - 1)]*res_y

        ! Initialize permutation table
        perm = random_permute(regspace(0, 1, 255), seed=seed)

        allocate (permutation_table(0:255, 0:255))
        do i = 0, 255
            do j = 0, 255
                permutation_table(i, j) = mod(perm(mod(i + perm(j), 256)), 4)
            end do
        end do

        ! Compute Perlin noise
        !$omp parallel do private(i, j, x0, x1, y0, y1, sx, sy, n0, n1, u, v)
        do j = 1, ny
            do i = 1, nx

                x0 = floor(x(i))
                x1 = x0 + 1
                y0 = floor(y(j))
                y1 = y0 + 1

                sx = fade(x(i) - x0)
                sy = fade(y(j) - y0)

                n0 = gradient_2d(permutation_table(mod(x0, 256), mod(y0, 256)), x(i) - x0, y(j) - y0)
                n1 = gradient_2d(permutation_table(mod(x1, 256), mod(y0, 256)), x(i) - x1, y(j) - y0)
                u = (1.0 - sx)*n0 + sx*n1

                n0 = gradient_2d(permutation_table(mod(x0, 256), mod(y1, 256)), x(i) - x0, y(j) - y1)
                n1 = gradient_2d(permutation_table(mod(x1, 256), mod(y1, 256)), x(i) - x1, y(j) - y1)
                v = (1.0 - sx)*n0 + sx*n1

                noise(i, j) = (1.0 - sy)*u + sy*v

            end do
        end do
        !$omp end parallel do

    end function perlin_2d

    !
    ! Generate 2D fractal noise based on Perlin noise
    !
    function generate_fractal_noise_2d(this) result(fractal_noise)

        class(fractal_noise_2d), intent(in) :: this

        integer :: n1, n2, res_x, res_y, seed, octaves
        real :: persistence, lacunarity
        real, allocatable, dimension(:, :) :: fractal_noise
        real, allocatable, dimension(:, :) :: noise
        real :: amplitude, frequency
        integer :: octave
        integer :: grid_check, nx, ny

        n1 = this%n1
        n2 = this%n2
        res_x = this%res1
        res_y = this%res2
        seed = this%seed*20
        octaves = this%octaves
        persistence = this%persistence
        lacunarity = this%lacunarity

        grid_check = floor(lacunarity**(octaves - 1)*res_x)
        nx = ceiling(n1*1.0/grid_check)*grid_check

        grid_check = floor(lacunarity**(octaves - 1)*res_y)
        ny = ceiling(n2*1.0/grid_check)*grid_check

        ! Initialize the fractal noise array
        fractal_noise = zeros(nx, ny)

        ! Start with base frequency and amplitude
        frequency = 1.0
        amplitude = 1.0

        ! Loop over octaves
        do octave = 1, octaves

            ! Generate Perlin noise for the current octave
            noise = perlin_2d(res_x*frequency, res_y*frequency, nx, ny, seed + octave)

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
    ! Compute 3D gradient dot product
    !
    pure function gradient_3d(perm, gradients, x0, y0, z0, x, y, z) result(dot)

        integer, intent(in), dimension(:) :: perm
        real, intent(in), dimension(:, :) :: gradients
        integer, intent(in) :: x0, y0, z0
        real, intent(in) :: x, y, z
        real :: dot
        integer :: h

        ! Compute hash value from permutation table
        h = perm(mod(x0 + perm(mod(y0 + perm(mod(z0, 256)), 256)), 256))
        h = mod(h, 12) + 1

        ! Dot product of gradient vector with relative position
        dot = gradients(h, 1)*(x - x0) + gradients(h, 2)*(y - y0) + gradients(h, 3)*(z - z0)

    end function gradient_3d

    !
    ! Generate 3D Perlin noise
    !
    function perlin_3d(res_x, res_y, res_z, nx, ny, nz, seed) result(noise)

        real, intent(in) :: res_x, res_y, res_z
        integer, intent(in) :: nx, ny, nz, seed
        real, allocatable, dimension(:, :, :) :: noise
        integer :: x0, x1, y0, y1, z0, z1, i, j, k
        real :: sx, sy, sz, n000, n100, n010, n110, n001, n101, n011, n111
        real :: nx00, nx10, nx01, nx11, nxy0, nxy1
        real, allocatable, dimension(:) :: x, y, z
        integer, allocatable, dimension(:) :: perm
        real, dimension(12, 3) :: gradients

        ! Gradients for Perlin noise (12 predefined unit vectors in 3D)
        gradients = reshape([ &
            1.0, 1.0, 0.0, &
            -1.0, 1.0, 0.0, &
            1.0, -1.0, 0.0, &
            -1.0, -1.0, 0.0, &
            1.0, 0.0, 1.0, &
            -1.0, 0.0, 1.0, &
            1.0, 0.0, -1.0, &
            -1.0, 0.0, -1.0, &
            0.0, 1.0, 1.0, &
            0.0, -1.0, 1.0, &
            0.0, 1.0, -1.0, &
            0.0, -1.0, -1.0], [12, 3])

        ! Allocate noise array and coordinate arrays
        allocate (noise(nx, ny, nz))
        allocate (x(nx))
        allocate (y(ny))
        allocate (z(nz))
        allocate (perm(0:255))

        ! Initialize coordinates scaled by res
        x = [(i/real(nx - 1), i=0, nx - 1)]*res_x
        y = [(j/real(ny - 1), j=0, ny - 1)]*res_y
        z = [(k/real(nz - 1), k=0, nz - 1)]*res_z

        ! Initialize permutation table
        perm = random_permute(regspace(0, 1, 255), seed=seed)

        ! Compute Perlin noise
        !$omp parallel do private(i, j, k, x0, x1, y0, y1, z0, z1, sx, sy, sz, &
            !$omp & n000, n100, n010, n110, n001, n101, n011, n111, nx00, nx10, nx01, nx11, nxy0, nxy1)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx

                    x0 = floor(x(i))
                    x1 = x0 + 1
                    y0 = floor(y(j))
                    y1 = y0 + 1
                    z0 = floor(z(k))
                    z1 = z0 + 1

                    sx = fade(x(i) - x0)
                    sy = fade(y(j) - y0)
                    sz = fade(z(k) - z0)

                    ! Compute gradients and dot products at the 8 cube corners
                    n000 = gradient_3d(perm, gradients, x0, y0, z0, x(i), y(j), z(k))
                    n100 = gradient_3d(perm, gradients, x1, y0, z0, x(i), y(j), z(k))
                    n010 = gradient_3d(perm, gradients, x0, y1, z0, x(i), y(j), z(k))
                    n110 = gradient_3d(perm, gradients, x1, y1, z0, x(i), y(j), z(k))
                    n001 = gradient_3d(perm, gradients, x0, y0, z1, x(i), y(j), z(k))
                    n101 = gradient_3d(perm, gradients, x1, y0, z1, x(i), y(j), z(k))
                    n011 = gradient_3d(perm, gradients, x0, y1, z1, x(i), y(j), z(k))
                    n111 = gradient_3d(perm, gradients, x1, y1, z1, x(i), y(j), z(k))

                    ! Interpolate along x-axis
                    nx00 = (1.0 - sx)*n000 + sx*n100
                    nx10 = (1.0 - sx)*n010 + sx*n110
                    nx01 = (1.0 - sx)*n001 + sx*n101
                    nx11 = (1.0 - sx)*n011 + sx*n111

                    ! Interpolate along y-axis
                    nxy0 = (1.0 - sy)*nx00 + sy*nx10
                    nxy1 = (1.0 - sy)*nx01 + sy*nx11

                    ! Interpolate along z-axis
                    noise(i, j, k) = (1.0 - sz)*nxy0 + sz*nxy1

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

        integer :: n1, n2, n3, res_x, res_y, res_z, seed, octaves
        real :: persistence, lacunarity
        real, allocatable, dimension(:, :, :) :: fractal_noise
        real, allocatable, dimension(:, :, :) :: noise
        real :: amplitude, frequency
        integer :: octave
        integer :: grid_check, nx, ny, nz

        n1 = this%n1
        n2 = this%n2
        n3 = this%n3
        res_x = this%res1
        res_y = this%res2
        res_z = this%res3
        seed = this%seed*20
        octaves = this%octaves
        persistence = this%persistence
        lacunarity = this%lacunarity

        grid_check = floor(lacunarity**(octaves - 1)*res_x)
        nx = ceiling(n1*1.0/grid_check)*grid_check

        grid_check = floor(lacunarity**(octaves - 1)*res_y)
        ny = ceiling(n2*1.0/grid_check)*grid_check

        grid_check = floor(lacunarity**(octaves - 1)*res_z)
        nz = ceiling(n3*1.0/grid_check)*grid_check

        ! Initialize the fractal noise array
        fractal_noise = zeros(nx, ny, nz)

        ! Start with base frequency and amplitude
        frequency = 1.0
        amplitude = 1.0

        ! Loop over octaves
        do octave = 1, octaves

            ! Generate Perlin noise for the current octave
            noise = perlin_3d(res_x*frequency, res_y*frequency, res_z*frequency, nx, ny, nz, seed + octave)

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
