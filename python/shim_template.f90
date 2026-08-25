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
!SET_CASES_2D
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
!SET_ARRAY_CASES_2D
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
!SET_STRING_CASES_2D
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
!GET_SHAPE_CASES_2D
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
!GET_CASES_2D
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
!SET_CASES_3D
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
!SET_ARRAY_CASES_3D
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
!SET_STRING_CASES_3D
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
!GET_SHAPE_CASES_3D
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
!GET_CASES_3D
        end select

    end function rgm3_get_array

end module geological_model_c
