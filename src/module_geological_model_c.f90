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
! C-interoperable interface to rgm2_curved/rgm3_curved for the Python wrapper.
!
! *** GENERATED FILE — do not edit by hand. ***
! Regenerate with: python3 python/generate_shim.py
! (parses the parameter components of the two derived types).
!

module geological_model_c

    use, intrinsic :: iso_c_binding
    use geological_model_2d_curved
    use geological_model_3d_curved

    implicit none

    integer, parameter :: max_pool = 64

    type :: rgm2_slot
        type(rgm2_curved), allocatable :: p
    end type rgm2_slot

    type :: rgm3_slot
        type(rgm3_curved), allocatable :: p
    end type rgm3_slot

    type(rgm2_slot), dimension(1:max_pool), save :: pool2
    type(rgm3_slot), dimension(1:max_pool), save :: pool3

contains

    subroutine cstr_to_f(cstr, n, fstr)

        character(kind=c_char), dimension(*), intent(in) :: cstr
        integer(c_int), intent(in), value :: n
        character(len=:), allocatable, intent(out) :: fstr

        integer :: i

        allocate (character(len=n) :: fstr)
        do i = 1, n
            fstr(i:i) = cstr(i)
        end do

    end subroutine cstr_to_f

    subroutine copy_out_2d(p, dims, data_ptr)

        real, dimension(:, :), intent(in) :: p
        integer(c_int), dimension(*), intent(out) :: dims
        real(c_float), dimension(*), intent(out) :: data_ptr

        integer :: i, j, l

        dims(1) = size(p, 1)
        dims(2) = size(p, 2)
        dims(3) = 1
        l = 0
        do j = 1, size(p, 2)
            do i = 1, size(p, 1)
                l = l + 1
                data_ptr(l) = p(i, j)
            end do
        end do

    end subroutine copy_out_2d

    subroutine copy_out_3d(p, dims, data_ptr)

        real, dimension(:, :, :), intent(in) :: p
        integer(c_int), dimension(*), intent(out) :: dims
        real(c_float), dimension(*), intent(out) :: data_ptr

        integer :: i, j, k, l

        dims(1) = size(p, 1)
        dims(2) = size(p, 2)
        dims(3) = size(p, 3)
        l = 0
        do k = 1, size(p, 3)
            do j = 1, size(p, 2)
                do i = 1, size(p, 1)
                    l = l + 1
                    data_ptr(l) = p(i, j, k)
                end do
            end do
        end do

    end subroutine copy_out_3d

    !===============================================================================
    ! 2D interface

    function rgm2_create() result(h) bind(c, name='rgm2_create')

        integer(c_int) :: h

        integer :: i

        h = -1
        do i = 1, max_pool
            if (.not. allocated(pool2(i)%p)) then
                allocate (pool2(i)%p)
                h = i
                return
            end if
        end do

    end function rgm2_create

    subroutine rgm2_free(h) bind(c, name='rgm2_free')

        integer(c_int), intent(in), value :: h

        if (h >= 1 .and. h <= max_pool) then
            if (allocated(pool2(h)%p)) then
                deallocate (pool2(h)%p)
            end if
        end if

    end subroutine rgm2_free

    function rgm2_set_num(h, cname, n, val) result(ok) bind(c, name='rgm2_set_num')

        integer(c_int), intent(in), value :: h, n
        character(kind=c_char), dimension(*), intent(in) :: cname
        real(c_double), intent(in), value :: val
        integer(c_int) :: ok

        character(len=:), allocatable :: name

        call cstr_to_f(cname, n, name)
        ok = 1

        select case (name)
            case ('n1')
                pool2(h)%p%n1 = nint(val)
            case ('n2')
                pool2(h)%p%n2 = nint(val)
            case ('nf')
                pool2(h)%p%nf = nint(val)
            case ('nl')
                pool2(h)%p%nl = nint(val)
            case ('seed')
                pool2(h)%p%seed = nint(val)
            case ('refl_smooth')
                pool2(h)%p%refl_smooth = real(val)
            case ('refl_smooth_top')
                pool2(h)%p%refl_smooth_top = real(val)
            case ('refl_slope')
                pool2(h)%p%refl_slope = real(val)
            case ('refl_slope_top')
                pool2(h)%p%refl_slope_top = real(val)
            case ('dt')
                pool2(h)%p%dt = real(val)
            case ('f0')
                pool2(h)%p%f0 = real(val)
            case ('fwidth')
                pool2(h)%p%fwidth = real(val)
            case ('yn_vary_disp')
                pool2(h)%p%yn_vary_disp = nint(val) /= 0
            case ('yn_disp_decay')
                pool2(h)%p%yn_disp_decay = nint(val) /= 0
            case ('noise_level')
                pool2(h)%p%noise_level = real(val)
            case ('ng')
                pool2(h)%p%ng = nint(val)
            case ('lwv')
                pool2(h)%p%lwv = real(val)
            case ('lwh')
                pool2(h)%p%lwh = real(val)
            case ('secondary_refl_smooth')
                pool2(h)%p%secondary_refl_smooth = real(val)
            case ('yn_rgt')
                pool2(h)%p%yn_rgt = nint(val) /= 0
            case ('yn_facies')
                pool2(h)%p%yn_facies = nint(val) /= 0
            case ('yn_fault')
                pool2(h)%p%yn_fault = nint(val) /= 0
            case ('custom_psf')
                pool2(h)%p%custom_psf = nint(val) /= 0
            case ('yn_conv_noise')
                pool2(h)%p%yn_conv_noise = nint(val) /= 0
            case ('yn_regular_fault')
                pool2(h)%p%yn_regular_fault = nint(val) /= 0
            case ('yn_group_faults')
                pool2(h)%p%yn_group_faults = nint(val) /= 0
            case ('vmin')
                pool2(h)%p%vmin = real(val)
            case ('vmax')
                pool2(h)%p%vmax = real(val)
            case ('delta_v')
                pool2(h)%p%delta_v = real(val)
            case ('unconf')
                pool2(h)%p%unconf = nint(val)
            case ('unconf_nl')
                pool2(h)%p%unconf_nl = nint(val)
            case ('unconf_smooth')
                pool2(h)%p%unconf_smooth = real(val)
            case ('unconf_refl_slope')
                pool2(h)%p%unconf_refl_slope = real(val)
            case ('unconf_refl_smooth')
                pool2(h)%p%unconf_refl_smooth = real(val)
            case ('unconf_channel_sinuosity')
                pool2(h)%p%unconf_channel_sinuosity = real(val)
            case ('unconf_channel_length')
                pool2(h)%p%unconf_channel_length = real(val)
            case ('unconf_topo')
                pool2(h)%p%unconf_topo = real(val)
            case ('yn_salt')
                pool2(h)%p%yn_salt = nint(val) /= 0
            case ('nsalt')
                pool2(h)%p%nsalt = nint(val)
            case ('salt_vp')
                pool2(h)%p%salt_vp = real(val)
            case ('salt_rho')
                pool2(h)%p%salt_rho = real(val)
            case ('salt_top_height')
                pool2(h)%p%salt_top_height = real(val)
            case ('salt_radius_variation')
                pool2(h)%p%salt_radius_variation = real(val)
            case ('salt_path_variation')
                pool2(h)%p%salt_path_variation = real(val)
            case ('salt_nnode')
                pool2(h)%p%salt_nnode = nint(val)
            case ('yn_karst')
                pool2(h)%p%yn_karst = nint(val) /= 0
            case ('karst_npassage')
                pool2(h)%p%karst_npassage = nint(val)
            case ('karst_nctrl')
                pool2(h)%p%karst_nctrl = nint(val)
            case ('karst_radius_variation')
                pool2(h)%p%karst_radius_variation = real(val)
            case ('karst_tortuosity')
                pool2(h)%p%karst_tortuosity = real(val)
            case ('karst_connect')
                pool2(h)%p%karst_connect = real(val)
            case ('karst_vp')
                pool2(h)%p%karst_vp = real(val)
            case ('karst_vs')
                pool2(h)%p%karst_vs = real(val)
            case ('karst_rho')
                pool2(h)%p%karst_rho = real(val)
            case ('karst_before_unconf')
                pool2(h)%p%karst_before_unconf = nint(val) /= 0
            case ('yn_elastic')
                pool2(h)%p%yn_elastic = nint(val) /= 0
            case ('rho_a')
                pool2(h)%p%rho_a = real(val)
            case ('rho_b')
                pool2(h)%p%rho_b = real(val)
            case ('rho_c')
                pool2(h)%p%rho_c = real(val)
            case ('salt_vs')
                pool2(h)%p%salt_vs = real(val)
            case ('salt_before_unconf')
                pool2(h)%p%salt_before_unconf = nint(val) /= 0
            case default
                ok = 0
        end select

    end function rgm2_set_num

    function rgm2_set_num_array(h, cname, n, vals, nval) result(ok) bind(c, name='rgm2_set_num_array')

        integer(c_int), intent(in), value :: h, n, nval
        character(kind=c_char), dimension(*), intent(in) :: cname
        real(c_double), dimension(nval), intent(in) :: vals
        integer(c_int) :: ok

        character(len=:), allocatable :: name

        call cstr_to_f(cname, n, name)
        ok = 1

        select case (name)
            case ('dip')
                pool2(h)%p%dip = real(vals)
            case ('disp')
                pool2(h)%p%disp = real(vals)
            case ('delta_dip')
                pool2(h)%p%delta_dip = real(vals)
            case ('disp_radius_dip')
                pool2(h)%p%disp_radius_dip = real(vals)
            case ('disp_center_dip')
                pool2(h)%p%disp_center_dip = real(vals)
            case ('disp_decay_width')
                pool2(h)%p%disp_decay_width = real(vals)
            case ('noise_smooth')
                pool2(h)%p%noise_smooth = real(vals)
            case ('refl_height')
                pool2(h)%p%refl_height = real(vals)
            case ('refl_height_top')
                pool2(h)%p%refl_height_top = real(vals)
            case ('refl_sigma2')
                pool2(h)%p%refl_sigma2 = real(vals)
            case ('refl_mu2')
                pool2(h)%p%refl_mu2 = real(vals)
            case ('psf_sigma')
                pool2(h)%p%psf_sigma = real(vals)
            case ('wave_filt_freqs')
                pool2(h)%p%wave_filt_freqs = real(vals)
            case ('wave_filt_amps')
                pool2(h)%p%wave_filt_amps = real(vals)
            case ('unconf_z')
                pool2(h)%p%unconf_z = real(vals)
            case ('unconf_height')
                pool2(h)%p%unconf_height = real(vals)
            case ('unconf_refl_height')
                pool2(h)%p%unconf_refl_height = real(vals)
            case ('unconf_channel_width')
                pool2(h)%p%unconf_channel_width = real(vals)
            case ('unconf_channel_density')
                pool2(h)%p%unconf_channel_density = real(vals)
            case ('salt_radius')
                pool2(h)%p%salt_radius = real(vals)
            case ('salt_top_z')
                pool2(h)%p%salt_top_z = real(vals)
            case ('karst_z')
                pool2(h)%p%karst_z = real(vals)
            case ('karst_radius')
                pool2(h)%p%karst_radius = real(vals)
            case ('vpvsratio')
                pool2(h)%p%vpvsratio = real(vals)
            case default
                ok = 0
        end select

    end function rgm2_set_num_array

    function rgm2_set_string(h, cname, n, cval, nval) result(ok) bind(c, name='rgm2_set_string')

        integer(c_int), intent(in), value :: h, n, nval
        character(kind=c_char), dimension(*), intent(in) :: cname, cval
        integer(c_int) :: ok

        character(len=:), allocatable :: name, str

        call cstr_to_f(cname, n, name)
        call cstr_to_f(cval, nval, str)
        ok = 1

        select case (name)
            case ('wave')
                pool2(h)%p%wave = trim(str)
            case ('refl_shape')
                pool2(h)%p%refl_shape = trim(str)
            case ('refl_shape_top')
                pool2(h)%p%refl_shape_top = trim(str)
            case ('noise_type')
                pool2(h)%p%noise_type = trim(str)
            case ('unconf_refl_shape')
                pool2(h)%p%unconf_refl_shape = trim(str)
            case ('unconf_shape')
                pool2(h)%p%unconf_shape = trim(str)
            case default
                ok = 0
        end select

    end function rgm2_set_string

    subroutine rgm2_generate(h) bind(c, name='rgm2_generate')

        integer(c_int), intent(in), value :: h

        call pool2(h)%p%generate

    end subroutine rgm2_generate

    function rgm2_get_shape(h, cname, n, dims) result(ok) bind(c, name='rgm2_get_shape')

        integer(c_int), intent(in), value :: h, n
        character(kind=c_char), dimension(*), intent(in) :: cname
        integer(c_int), dimension(*), intent(out) :: dims
        integer(c_int) :: ok

        character(len=:), allocatable :: name

        call cstr_to_f(cname, n, name)
        ok = 0
        dims(1:3) = 0

        select case (name)
            case ('image')
                if (allocated(pool2(h)%p%image)) then
                    dims(1) = size(pool2(h)%p%image, 1); dims(2) = size(pool2(h)%p%image, 2)
                    ok = 1
                end if
            case ('rgt')
                if (allocated(pool2(h)%p%rgt)) then
                    dims(1) = size(pool2(h)%p%rgt, 1); dims(2) = size(pool2(h)%p%rgt, 2)
                    ok = 1
                end if
            case ('facies')
                if (allocated(pool2(h)%p%facies)) then
                    dims(1) = size(pool2(h)%p%facies, 1); dims(2) = size(pool2(h)%p%facies, 2)
                    ok = 1
                end if
            case ('fault')
                if (allocated(pool2(h)%p%fault)) then
                    dims(1) = size(pool2(h)%p%fault, 1); dims(2) = size(pool2(h)%p%fault, 2)
                    ok = 1
                end if
            case ('fault_dip')
                if (allocated(pool2(h)%p%fault_dip)) then
                    dims(1) = size(pool2(h)%p%fault_dip, 1); dims(2) = size(pool2(h)%p%fault_dip, 2)
                    ok = 1
                end if
            case ('fault_disp')
                if (allocated(pool2(h)%p%fault_disp)) then
                    dims(1) = size(pool2(h)%p%fault_disp, 1); dims(2) = size(pool2(h)%p%fault_disp, 2)
                    ok = 1
                end if
            case ('salt')
                if (allocated(pool2(h)%p%salt)) then
                    dims(1) = size(pool2(h)%p%salt, 1); dims(2) = size(pool2(h)%p%salt, 2)
                    ok = 1
                end if
            case ('karst')
                if (allocated(pool2(h)%p%karst)) then
                    dims(1) = size(pool2(h)%p%karst, 1); dims(2) = size(pool2(h)%p%karst, 2)
                    ok = 1
                end if
            case ('vp')
                if (allocated(pool2(h)%p%vp)) then
                    dims(1) = size(pool2(h)%p%vp, 1); dims(2) = size(pool2(h)%p%vp, 2)
                    ok = 1
                end if
            case ('vs')
                if (allocated(pool2(h)%p%vs)) then
                    dims(1) = size(pool2(h)%p%vs, 1); dims(2) = size(pool2(h)%p%vs, 2)
                    ok = 1
                end if
            case ('rho')
                if (allocated(pool2(h)%p%rho)) then
                    dims(1) = size(pool2(h)%p%rho, 1); dims(2) = size(pool2(h)%p%rho, 2)
                    ok = 1
                end if
            case ('psf')
                if (allocated(pool2(h)%p%psf)) then
                    dims(1) = size(pool2(h)%p%psf, 1); dims(2) = size(pool2(h)%p%psf, 2)
                    ok = 1
                end if
            case ('image_pp')
                if (allocated(pool2(h)%p%image_pp)) then
                    dims(1) = size(pool2(h)%p%image_pp, 1); dims(2) = size(pool2(h)%p%image_pp, 2)
                    ok = 1
                end if
            case ('image_ps')
                if (allocated(pool2(h)%p%image_ps)) then
                    dims(1) = size(pool2(h)%p%image_ps, 1); dims(2) = size(pool2(h)%p%image_ps, 2)
                    ok = 1
                end if
            case ('image_sp')
                if (allocated(pool2(h)%p%image_sp)) then
                    dims(1) = size(pool2(h)%p%image_sp, 1); dims(2) = size(pool2(h)%p%image_sp, 2)
                    ok = 1
                end if
            case ('image_ss')
                if (allocated(pool2(h)%p%image_ss)) then
                    dims(1) = size(pool2(h)%p%image_ss, 1); dims(2) = size(pool2(h)%p%image_ss, 2)
                    ok = 1
                end if
        end select

    end function rgm2_get_shape

    function rgm2_get_array(h, cname, n, dims, data_ptr) result(ok) bind(c, name='rgm2_get_array')

        integer(c_int), intent(in), value :: h, n
        character(kind=c_char), dimension(*), intent(in) :: cname
        integer(c_int), dimension(*), intent(out) :: dims
        real(c_float), dimension(*), intent(out) :: data_ptr
        integer(c_int) :: ok

        character(len=:), allocatable :: name

        call cstr_to_f(cname, n, name)
        ok = 0

        select case (name)
            case ('image')
                if (allocated(pool2(h)%p%image)) then
                    call copy_out_2d(pool2(h)%p%image, dims, data_ptr)
                    ok = 1
                end if
            case ('rgt')
                if (allocated(pool2(h)%p%rgt)) then
                    call copy_out_2d(pool2(h)%p%rgt, dims, data_ptr)
                    ok = 1
                end if
            case ('facies')
                if (allocated(pool2(h)%p%facies)) then
                    call copy_out_2d(pool2(h)%p%facies, dims, data_ptr)
                    ok = 1
                end if
            case ('fault')
                if (allocated(pool2(h)%p%fault)) then
                    call copy_out_2d(pool2(h)%p%fault, dims, data_ptr)
                    ok = 1
                end if
            case ('fault_dip')
                if (allocated(pool2(h)%p%fault_dip)) then
                    call copy_out_2d(pool2(h)%p%fault_dip, dims, data_ptr)
                    ok = 1
                end if
            case ('fault_disp')
                if (allocated(pool2(h)%p%fault_disp)) then
                    call copy_out_2d(pool2(h)%p%fault_disp, dims, data_ptr)
                    ok = 1
                end if
            case ('salt')
                if (allocated(pool2(h)%p%salt)) then
                    call copy_out_2d(pool2(h)%p%salt, dims, data_ptr)
                    ok = 1
                end if
            case ('karst')
                if (allocated(pool2(h)%p%karst)) then
                    call copy_out_2d(pool2(h)%p%karst, dims, data_ptr)
                    ok = 1
                end if
            case ('vp')
                if (allocated(pool2(h)%p%vp)) then
                    call copy_out_2d(pool2(h)%p%vp, dims, data_ptr)
                    ok = 1
                end if
            case ('vs')
                if (allocated(pool2(h)%p%vs)) then
                    call copy_out_2d(pool2(h)%p%vs, dims, data_ptr)
                    ok = 1
                end if
            case ('rho')
                if (allocated(pool2(h)%p%rho)) then
                    call copy_out_2d(pool2(h)%p%rho, dims, data_ptr)
                    ok = 1
                end if
            case ('psf')
                if (allocated(pool2(h)%p%psf)) then
                    call copy_out_2d(pool2(h)%p%psf, dims, data_ptr)
                    ok = 1
                end if
            case ('image_pp')
                if (allocated(pool2(h)%p%image_pp)) then
                    call copy_out_2d(pool2(h)%p%image_pp, dims, data_ptr)
                    ok = 1
                end if
            case ('image_ps')
                if (allocated(pool2(h)%p%image_ps)) then
                    call copy_out_2d(pool2(h)%p%image_ps, dims, data_ptr)
                    ok = 1
                end if
            case ('image_sp')
                if (allocated(pool2(h)%p%image_sp)) then
                    call copy_out_2d(pool2(h)%p%image_sp, dims, data_ptr)
                    ok = 1
                end if
            case ('image_ss')
                if (allocated(pool2(h)%p%image_ss)) then
                    call copy_out_2d(pool2(h)%p%image_ss, dims, data_ptr)
                    ok = 1
                end if
        end select

    end function rgm2_get_array

    !===============================================================================
    ! 3D interface

    function rgm3_create() result(h) bind(c, name='rgm3_create')

        integer(c_int) :: h

        integer :: i

        h = -1
        do i = 1, max_pool
            if (.not. allocated(pool3(i)%p)) then
                allocate (pool3(i)%p)
                h = i
                return
            end if
        end do

    end function rgm3_create

    subroutine rgm3_free(h) bind(c, name='rgm3_free')

        integer(c_int), intent(in), value :: h

        if (h >= 1 .and. h <= max_pool) then
            if (allocated(pool3(h)%p)) then
                deallocate (pool3(h)%p)
            end if
        end if

    end subroutine rgm3_free

    function rgm3_set_num(h, cname, n, val) result(ok) bind(c, name='rgm3_set_num')

        integer(c_int), intent(in), value :: h, n
        character(kind=c_char), dimension(*), intent(in) :: cname
        real(c_double), intent(in), value :: val
        integer(c_int) :: ok

        character(len=:), allocatable :: name

        call cstr_to_f(cname, n, name)
        ok = 1

        select case (name)
            case ('n1')
                pool3(h)%p%n1 = nint(val)
            case ('n2')
                pool3(h)%p%n2 = nint(val)
            case ('n3')
                pool3(h)%p%n3 = nint(val)
            case ('nf')
                pool3(h)%p%nf = nint(val)
            case ('nl')
                pool3(h)%p%nl = nint(val)
            case ('seed')
                pool3(h)%p%seed = nint(val)
            case ('refl_smooth')
                pool3(h)%p%refl_smooth = real(val)
            case ('refl_smooth_top')
                pool3(h)%p%refl_smooth_top = real(val)
            case ('dt')
                pool3(h)%p%dt = real(val)
            case ('f0')
                pool3(h)%p%f0 = real(val)
            case ('fwidth')
                pool3(h)%p%fwidth = real(val)
            case ('noise_level')
                pool3(h)%p%noise_level = real(val)
            case ('ng')
                pool3(h)%p%ng = nint(val)
            case ('lwv')
                pool3(h)%p%lwv = real(val)
            case ('lwh')
                pool3(h)%p%lwh = real(val)
            case ('secondary_refl_smooth')
                pool3(h)%p%secondary_refl_smooth = real(val)
            case ('rotate_fold')
                pool3(h)%p%rotate_fold = nint(val) /= 0
            case ('yn_rgt')
                pool3(h)%p%yn_rgt = nint(val) /= 0
            case ('yn_facies')
                pool3(h)%p%yn_facies = nint(val) /= 0
            case ('yn_fault')
                pool3(h)%p%yn_fault = nint(val) /= 0
            case ('custom_psf')
                pool3(h)%p%custom_psf = nint(val) /= 0
            case ('facies_threshold')
                pool3(h)%p%facies_threshold = real(val)
            case ('yn_conv_noise')
                pool3(h)%p%yn_conv_noise = nint(val) /= 0
            case ('yn_regular_fault')
                pool3(h)%p%yn_regular_fault = nint(val) /= 0
            case ('yn_group_faults')
                pool3(h)%p%yn_group_faults = nint(val) /= 0
            case ('vmin')
                pool3(h)%p%vmin = real(val)
            case ('vmax')
                pool3(h)%p%vmax = real(val)
            case ('delta_v')
                pool3(h)%p%delta_v = real(val)
            case ('strike_nperiod')
                pool3(h)%p%strike_nperiod = nint(val)
            case ('yn_vary_disp')
                pool3(h)%p%yn_vary_disp = nint(val) /= 0
            case ('yn_disp_decay')
                pool3(h)%p%yn_disp_decay = nint(val) /= 0
            case ('unconf')
                pool3(h)%p%unconf = nint(val)
            case ('unconf_nl')
                pool3(h)%p%unconf_nl = nint(val)
            case ('unconf_smooth')
                pool3(h)%p%unconf_smooth = real(val)
            case ('unconf_refl_slope')
                pool3(h)%p%unconf_refl_slope = real(val)
            case ('unconf_refl_smooth')
                pool3(h)%p%unconf_refl_smooth = real(val)
            case ('unconf_channel_sinuosity')
                pool3(h)%p%unconf_channel_sinuosity = real(val)
            case ('unconf_channel_length')
                pool3(h)%p%unconf_channel_length = real(val)
            case ('unconf_topo')
                pool3(h)%p%unconf_topo = real(val)
            case ('yn_salt')
                pool3(h)%p%yn_salt = nint(val) /= 0
            case ('nsalt')
                pool3(h)%p%nsalt = nint(val)
            case ('salt_radius_variation')
                pool3(h)%p%salt_radius_variation = real(val)
            case ('salt_path_variation')
                pool3(h)%p%salt_path_variation = real(val)
            case ('salt_nnode')
                pool3(h)%p%salt_nnode = nint(val)
            case ('salt_vp')
                pool3(h)%p%salt_vp = real(val)
            case ('salt_rho')
                pool3(h)%p%salt_rho = real(val)
            case ('salt_top_height')
                pool3(h)%p%salt_top_height = real(val)
            case ('yn_karst')
                pool3(h)%p%yn_karst = nint(val) /= 0
            case ('karst_npassage')
                pool3(h)%p%karst_npassage = nint(val)
            case ('karst_nctrl')
                pool3(h)%p%karst_nctrl = nint(val)
            case ('karst_radius_variation')
                pool3(h)%p%karst_radius_variation = real(val)
            case ('karst_tortuosity')
                pool3(h)%p%karst_tortuosity = real(val)
            case ('karst_connect')
                pool3(h)%p%karst_connect = real(val)
            case ('karst_vp')
                pool3(h)%p%karst_vp = real(val)
            case ('karst_vs')
                pool3(h)%p%karst_vs = real(val)
            case ('karst_rho')
                pool3(h)%p%karst_rho = real(val)
            case ('karst_before_unconf')
                pool3(h)%p%karst_before_unconf = nint(val) /= 0
            case ('yn_elastic')
                pool3(h)%p%yn_elastic = nint(val) /= 0
            case ('rho_a')
                pool3(h)%p%rho_a = real(val)
            case ('rho_b')
                pool3(h)%p%rho_b = real(val)
            case ('rho_c')
                pool3(h)%p%rho_c = real(val)
            case ('salt_vs')
                pool3(h)%p%salt_vs = real(val)
            case ('salt_before_unconf')
                pool3(h)%p%salt_before_unconf = nint(val) /= 0
            case default
                ok = 0
        end select

    end function rgm3_set_num

    function rgm3_set_num_array(h, cname, n, vals, nval) result(ok) bind(c, name='rgm3_set_num_array')

        integer(c_int), intent(in), value :: h, n, nval
        character(kind=c_char), dimension(*), intent(in) :: cname
        real(c_double), dimension(nval), intent(in) :: vals
        integer(c_int) :: ok

        character(len=:), allocatable :: name

        call cstr_to_f(cname, n, name)
        ok = 1

        select case (name)
            case ('dip')
                pool3(h)%p%dip = real(vals)
            case ('strike')
                pool3(h)%p%strike = real(vals)
            case ('rake')
                pool3(h)%p%rake = real(vals)
            case ('disp')
                pool3(h)%p%disp = real(vals)
            case ('refl_slope')
                pool3(h)%p%refl_slope = real(vals)
            case ('refl_slope_top')
                pool3(h)%p%refl_slope_top = real(vals)
            case ('noise_smooth')
                pool3(h)%p%noise_smooth = real(vals)
            case ('refl_height')
                pool3(h)%p%refl_height = real(vals)
            case ('refl_height_top')
                pool3(h)%p%refl_height_top = real(vals)
            case ('refl_sigma2')
                pool3(h)%p%refl_sigma2 = real(vals)
            case ('refl_mu2')
                pool3(h)%p%refl_mu2 = real(vals)
            case ('refl_sigma3')
                pool3(h)%p%refl_sigma3 = real(vals)
            case ('refl_mu3')
                pool3(h)%p%refl_mu3 = real(vals)
            case ('psf_sigma')
                pool3(h)%p%psf_sigma = real(vals)
            case ('wave_filt_freqs')
                pool3(h)%p%wave_filt_freqs = real(vals)
            case ('wave_filt_amps')
                pool3(h)%p%wave_filt_amps = real(vals)
            case ('delta_dip')
                pool3(h)%p%delta_dip = real(vals)
            case ('delta_strike')
                pool3(h)%p%delta_strike = real(vals)
            case ('disp_radius_strike')
                pool3(h)%p%disp_radius_strike = real(vals)
            case ('disp_radius_dip')
                pool3(h)%p%disp_radius_dip = real(vals)
            case ('disp_center_dip')
                pool3(h)%p%disp_center_dip = real(vals)
            case ('disp_center_strike')
                pool3(h)%p%disp_center_strike = real(vals)
            case ('disp_decay_width')
                pool3(h)%p%disp_decay_width = real(vals)
            case ('unconf_z')
                pool3(h)%p%unconf_z = real(vals)
            case ('unconf_height')
                pool3(h)%p%unconf_height = real(vals)
            case ('unconf_refl_height')
                pool3(h)%p%unconf_refl_height = real(vals)
            case ('unconf_channel_width')
                pool3(h)%p%unconf_channel_width = real(vals)
            case ('unconf_channel_density')
                pool3(h)%p%unconf_channel_density = real(vals)
            case ('salt_radius')
                pool3(h)%p%salt_radius = real(vals)
            case ('salt_top_z')
                pool3(h)%p%salt_top_z = real(vals)
            case ('karst_z')
                pool3(h)%p%karst_z = real(vals)
            case ('karst_radius')
                pool3(h)%p%karst_radius = real(vals)
            case ('vpvsratio')
                pool3(h)%p%vpvsratio = real(vals)
            case default
                ok = 0
        end select

    end function rgm3_set_num_array

    function rgm3_set_string(h, cname, n, cval, nval) result(ok) bind(c, name='rgm3_set_string')

        integer(c_int), intent(in), value :: h, n, nval
        character(kind=c_char), dimension(*), intent(in) :: cname, cval
        integer(c_int) :: ok

        character(len=:), allocatable :: name, str

        call cstr_to_f(cname, n, name)
        call cstr_to_f(cval, nval, str)
        ok = 1

        select case (name)
            case ('wave')
                pool3(h)%p%wave = trim(str)
            case ('refl_shape')
                pool3(h)%p%refl_shape = trim(str)
            case ('refl_shape_top')
                pool3(h)%p%refl_shape_top = trim(str)
            case ('noise_type')
                pool3(h)%p%noise_type = trim(str)
            case ('unconf_refl_shape')
                pool3(h)%p%unconf_refl_shape = trim(str)
            case ('unconf_shape')
                pool3(h)%p%unconf_shape = trim(str)
            case default
                ok = 0
        end select

    end function rgm3_set_string

    subroutine rgm3_generate(h) bind(c, name='rgm3_generate')

        integer(c_int), intent(in), value :: h

        call pool3(h)%p%generate

    end subroutine rgm3_generate

    function rgm3_get_shape(h, cname, n, dims) result(ok) bind(c, name='rgm3_get_shape')

        integer(c_int), intent(in), value :: h, n
        character(kind=c_char), dimension(*), intent(in) :: cname
        integer(c_int), dimension(*), intent(out) :: dims
        integer(c_int) :: ok

        character(len=:), allocatable :: name

        call cstr_to_f(cname, n, name)
        ok = 0
        dims(1:3) = 0

        select case (name)
            case ('image')
                if (allocated(pool3(h)%p%image)) then
                    dims(1) = size(pool3(h)%p%image, 1); dims(2) = size(pool3(h)%p%image, 2); dims(3) = size(pool3(h)%p%image, 3)
                    ok = 1
                end if
            case ('rgt')
                if (allocated(pool3(h)%p%rgt)) then
                    dims(1) = size(pool3(h)%p%rgt, 1); dims(2) = size(pool3(h)%p%rgt, 2); dims(3) = size(pool3(h)%p%rgt, 3)
                    ok = 1
                end if
            case ('facies')
                if (allocated(pool3(h)%p%facies)) then
                    dims(1) = size(pool3(h)%p%facies, 1); dims(2) = size(pool3(h)%p%facies, 2); dims(3) = size(pool3(h)%p%facies, 3)
                    ok = 1
                end if
            case ('fault')
                if (allocated(pool3(h)%p%fault)) then
                    dims(1) = size(pool3(h)%p%fault, 1); dims(2) = size(pool3(h)%p%fault, 2); dims(3) = size(pool3(h)%p%fault, 3)
                    ok = 1
                end if
            case ('fault_dip')
                if (allocated(pool3(h)%p%fault_dip)) then
                    dims(1) = size(pool3(h)%p%fault_dip, 1); dims(2) = size(pool3(h)%p%fault_dip, 2); dims(3) = size(pool3(h)%p%fault_dip, 3)
                    ok = 1
                end if
            case ('fault_disp')
                if (allocated(pool3(h)%p%fault_disp)) then
                    dims(1) = size(pool3(h)%p%fault_disp, 1); dims(2) = size(pool3(h)%p%fault_disp, 2); dims(3) = size(pool3(h)%p%fault_disp, 3)
                    ok = 1
                end if
            case ('salt')
                if (allocated(pool3(h)%p%salt)) then
                    dims(1) = size(pool3(h)%p%salt, 1); dims(2) = size(pool3(h)%p%salt, 2); dims(3) = size(pool3(h)%p%salt, 3)
                    ok = 1
                end if
            case ('karst')
                if (allocated(pool3(h)%p%karst)) then
                    dims(1) = size(pool3(h)%p%karst, 1); dims(2) = size(pool3(h)%p%karst, 2); dims(3) = size(pool3(h)%p%karst, 3)
                    ok = 1
                end if
            case ('vp')
                if (allocated(pool3(h)%p%vp)) then
                    dims(1) = size(pool3(h)%p%vp, 1); dims(2) = size(pool3(h)%p%vp, 2); dims(3) = size(pool3(h)%p%vp, 3)
                    ok = 1
                end if
            case ('vs')
                if (allocated(pool3(h)%p%vs)) then
                    dims(1) = size(pool3(h)%p%vs, 1); dims(2) = size(pool3(h)%p%vs, 2); dims(3) = size(pool3(h)%p%vs, 3)
                    ok = 1
                end if
            case ('rho')
                if (allocated(pool3(h)%p%rho)) then
                    dims(1) = size(pool3(h)%p%rho, 1); dims(2) = size(pool3(h)%p%rho, 2); dims(3) = size(pool3(h)%p%rho, 3)
                    ok = 1
                end if
            case ('psf')
                if (allocated(pool3(h)%p%psf)) then
                    dims(1) = size(pool3(h)%p%psf, 1); dims(2) = size(pool3(h)%p%psf, 2); dims(3) = size(pool3(h)%p%psf, 3)
                    ok = 1
                end if
            case ('image_pp')
                if (allocated(pool3(h)%p%image_pp)) then
                    dims(1) = size(pool3(h)%p%image_pp, 1); dims(2) = size(pool3(h)%p%image_pp, 2); dims(3) = size(pool3(h)%p%image_pp, 3)
                    ok = 1
                end if
            case ('image_ps')
                if (allocated(pool3(h)%p%image_ps)) then
                    dims(1) = size(pool3(h)%p%image_ps, 1); dims(2) = size(pool3(h)%p%image_ps, 2); dims(3) = size(pool3(h)%p%image_ps, 3)
                    ok = 1
                end if
            case ('image_sp')
                if (allocated(pool3(h)%p%image_sp)) then
                    dims(1) = size(pool3(h)%p%image_sp, 1); dims(2) = size(pool3(h)%p%image_sp, 2); dims(3) = size(pool3(h)%p%image_sp, 3)
                    ok = 1
                end if
            case ('image_ss')
                if (allocated(pool3(h)%p%image_ss)) then
                    dims(1) = size(pool3(h)%p%image_ss, 1); dims(2) = size(pool3(h)%p%image_ss, 2); dims(3) = size(pool3(h)%p%image_ss, 3)
                    ok = 1
                end if
            case ('fault_strike')
                if (allocated(pool3(h)%p%fault_strike)) then
                    dims(1) = size(pool3(h)%p%fault_strike, 1); dims(2) = size(pool3(h)%p%fault_strike, 2); dims(3) = size(pool3(h)%p%fault_strike, 3)
                    ok = 1
                end if
            case ('fault_rake')
                if (allocated(pool3(h)%p%fault_rake)) then
                    dims(1) = size(pool3(h)%p%fault_rake, 1); dims(2) = size(pool3(h)%p%fault_rake, 2); dims(3) = size(pool3(h)%p%fault_rake, 3)
                    ok = 1
                end if
        end select

    end function rgm3_get_shape

    function rgm3_get_array(h, cname, n, dims, data_ptr) result(ok) bind(c, name='rgm3_get_array')

        integer(c_int), intent(in), value :: h, n
        character(kind=c_char), dimension(*), intent(in) :: cname
        integer(c_int), dimension(*), intent(out) :: dims
        real(c_float), dimension(*), intent(out) :: data_ptr
        integer(c_int) :: ok

        character(len=:), allocatable :: name

        call cstr_to_f(cname, n, name)
        ok = 0

        select case (name)
            case ('image')
                if (allocated(pool3(h)%p%image)) then
                    call copy_out_3d(pool3(h)%p%image, dims, data_ptr)
                    ok = 1
                end if
            case ('rgt')
                if (allocated(pool3(h)%p%rgt)) then
                    call copy_out_3d(pool3(h)%p%rgt, dims, data_ptr)
                    ok = 1
                end if
            case ('facies')
                if (allocated(pool3(h)%p%facies)) then
                    call copy_out_3d(pool3(h)%p%facies, dims, data_ptr)
                    ok = 1
                end if
            case ('fault')
                if (allocated(pool3(h)%p%fault)) then
                    call copy_out_3d(pool3(h)%p%fault, dims, data_ptr)
                    ok = 1
                end if
            case ('fault_dip')
                if (allocated(pool3(h)%p%fault_dip)) then
                    call copy_out_3d(pool3(h)%p%fault_dip, dims, data_ptr)
                    ok = 1
                end if
            case ('fault_disp')
                if (allocated(pool3(h)%p%fault_disp)) then
                    call copy_out_3d(pool3(h)%p%fault_disp, dims, data_ptr)
                    ok = 1
                end if
            case ('salt')
                if (allocated(pool3(h)%p%salt)) then
                    call copy_out_3d(pool3(h)%p%salt, dims, data_ptr)
                    ok = 1
                end if
            case ('karst')
                if (allocated(pool3(h)%p%karst)) then
                    call copy_out_3d(pool3(h)%p%karst, dims, data_ptr)
                    ok = 1
                end if
            case ('vp')
                if (allocated(pool3(h)%p%vp)) then
                    call copy_out_3d(pool3(h)%p%vp, dims, data_ptr)
                    ok = 1
                end if
            case ('vs')
                if (allocated(pool3(h)%p%vs)) then
                    call copy_out_3d(pool3(h)%p%vs, dims, data_ptr)
                    ok = 1
                end if
            case ('rho')
                if (allocated(pool3(h)%p%rho)) then
                    call copy_out_3d(pool3(h)%p%rho, dims, data_ptr)
                    ok = 1
                end if
            case ('psf')
                if (allocated(pool3(h)%p%psf)) then
                    call copy_out_3d(pool3(h)%p%psf, dims, data_ptr)
                    ok = 1
                end if
            case ('image_pp')
                if (allocated(pool3(h)%p%image_pp)) then
                    call copy_out_3d(pool3(h)%p%image_pp, dims, data_ptr)
                    ok = 1
                end if
            case ('image_ps')
                if (allocated(pool3(h)%p%image_ps)) then
                    call copy_out_3d(pool3(h)%p%image_ps, dims, data_ptr)
                    ok = 1
                end if
            case ('image_sp')
                if (allocated(pool3(h)%p%image_sp)) then
                    call copy_out_3d(pool3(h)%p%image_sp, dims, data_ptr)
                    ok = 1
                end if
            case ('image_ss')
                if (allocated(pool3(h)%p%image_ss)) then
                    call copy_out_3d(pool3(h)%p%image_ss, dims, data_ptr)
                    ok = 1
                end if
            case ('fault_strike')
                if (allocated(pool3(h)%p%fault_strike)) then
                    call copy_out_3d(pool3(h)%p%fault_strike, dims, data_ptr)
                    ok = 1
                end if
            case ('fault_rake')
                if (allocated(pool3(h)%p%fault_rake)) then
                    call copy_out_3d(pool3(h)%p%fault_rake, dims, data_ptr)
                    ok = 1
                end if
        end select

    end function rgm3_get_array

end module geological_model_c
