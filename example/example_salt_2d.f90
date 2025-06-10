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


program main

    use libflit
    use librgm

    implicit none

    type(rgm2) :: p
    integer :: i, nm
    real, allocatable, dimension(:) :: psfm, fs, noise, height, slope, lwv, height2, lwv2
    integer, allocatable, dimension(:) :: nf, nl
    character(len=1024) :: dir_output

    dir_output = './dataset2_salt'
    nm = 100

    call make_directory(tidy(dir_output))

    fs = random(nm, range=[150.0, 210.0])
    nf = irandom(nm,  range=[4, 16])
    nl = nint(rescale(fs, range=[25.0, 50.0]))
    noise = random(nm,  range=[0.2, 0.4])
    height = random(nm,  range=[2.0, 12.0])
    height2 = random(nm,  range=[10.0, 20.0])
    lwv = random(nm,  range=[0.1, 0.4])
    lwv2 = random(nm,  range=[0.3, 0.8])
    slope = random(nm,  range=[-30.0, 30.0])
    p%dip = [50.0, 130.0]
    psfm = random(nm, range=[1.0, 3.0])

    do i = 1, 10

        p%n1 = 256
        p%n2 = 256
        p%nf = nf(i)
        p%refl_slope = slope(i)
        p%nl = nl(i)
        p%refl_amp = [0.1, 1.0]
        p%disp = [5.0, 30.0]
        p%fwidth = 3.0
        if (mod(irand(range=[1, nm]), 2) == 0) then
            p%refl_shape = 'gaussian'
            p%refl_mu2 = [0.0, p%n2 - 1.0]
            p%refl_sigma2 = [20.0, 50.0]
            p%ng = irand(range=[2, 6])
            p%refl_height = [0.0, height2(i)]
            p%lwv = lwv2(i)
            p%secondary_refl_height_ratio = rand(range=[0.0, 0.2])

        else
            p%refl_shape = 'random'
            p%refl_smooth = 25
            p%refl_height = [0.0, height(i)]
            p%lwv = lwv(i)
            p%secondary_refl_height_ratio = rand(range=[0.0, 0.2])
        end if
        select case (mod(irand(range=[1, nm]), 3))
            case (0)
                p%wave = 'ricker'
            case (1)
                p%wave = 'ricker_deriv'
            case (2)
                p%wave = 'gaussian_deriv'
        end select

        if (mod(irand(range=[1, nm]), 3) /= 0) then
            p%refl_amp_dist = 'normal'
        else
            p%refl_amp_dist = 'uniform'
        end if

        if (mod(irand(range=[1, nm]), 3) == 0) then
            p%unconf = 2
            p%unconf_amp = [0.05, 0.1]
            p%unconf_z = [0.1, 0.7]
            p%salt_top_max = [0.75, 0.9]
        else
            p%unconf = 0
            p%salt_top_max = [0.3, 0.9]
        end if

        if (mod(irand(range=[1, nm]), 3) == 0) then
            p%yn_regular_fault = .true.
            p%nf = irand(range=[16, 30])
            if (mod(irand(range=[1, nm]), 2) == 0) then
                p%dip = [rand(range=[100.0, 120.0]), rand(range=[60.0, 80.0])]
                p%disp = [3.0, -3.0]
            else
                p%dip = [rand(range=[60.0, 80.0]), rand(range=[100.0, 120.0])]
                p%disp = [-3.0, 3.0]
            end if
        else
            p%yn_regular_fault = .false.
            p%nf = nf(i)
            p%disp = [2.0, 15.0]
            p%dip = [50.0, 130.0]
        end if

        ! Clean image
        p%f0 = fs(i)
        p%noise_level = 0
        p%psf_sigma = [12.0, 0.0]
        p%nstem = 6
        p%nsalt = irand(range=[1, 4])
        p%salt_max_radius = [0.05, 0.4]
        p%salt_top_smooth = random(range=[2.0, 15.0])
        p%yn_fault = .true.
        p%yn_salt = .true.
        p%salt_amp = rand(range=[1.0, 3.0])
        call p%generate
        where (p%fault /= 0)
            p%fault = 1.0
        end where

        call output_array(p%image, tidy(dir_output)//'/'//num2str(i - 1)//'_img.bin')
        call output_array(p%salt, tidy(dir_output)//'/'//num2str(i - 1)//'_salt.bin')
        call output_array(p%fault, tidy(dir_output)//'/'//num2str(i - 1)//'_fsem.bin')
        call output_array(p%fault_dip/180.0, tidy(dir_output)//'/'//num2str(i - 1)//'_fdip.bin')

        print *, date_time_compact(), i

    end do

end program main
