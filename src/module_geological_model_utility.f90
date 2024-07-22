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

end module geological_model_utility
