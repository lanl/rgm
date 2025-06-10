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


program test

    use libflit
    use librgm

    block

        real, allocatable, dimension(:, :) :: vp
        integer :: n1, n2
        type(rgm2) :: p

        n1 = 400
        n2 = 420

        call make_directory('./model')

        p%n1 = n1
        p%n2 = n2
        p%nf = 17
        p%yn_regular_fault = .true.
        p%nl = 40
        p%unconf = 1
        p%unconf_z = 0.5
        p%disp = [10, 20]
        p%lwv = 0.1
        p%dip = [rand(range=[100.0, 120.0], seed=123), rand(range=[60.0, 80.0], seed=567)]
        p%disp = [10.0, -10.0]
        p%secondary_refl_height_ratio = 0.1
        p%refl_smooth = 10.0
        p%yn_facies = .true.
        p%yn_rgt = .true.
        p%fwidth = 3
        p%seed = 324343
        call p%generate
        vp = rescale(p%facies, [1000.0, 3000.0])
        call output_array(vp, './model/vp2_1.bin')
        call output_array(p%fault, './model/fault2_1.bin')
        call output_array(p%rgt, './model/rgt2_1.bin')
        p%image = p%image/maxval(p%image)
        call output_array(p%image, './model/img2_1.bin')

        p%n1 = n1
        p%n2 = n2
        p%nf = 4
        p%nl = 40
        p%unconf = 0
        p%yn_regular_fault = .false.
        p%disp = [10, 20]
        p%lwv = 0.1
        p%dip = [rand(range=[100.0, 120.0], seed=123), rand(range=[60.0, 80.0], seed=567)]
        p%secondary_refl_height_ratio = 0.05
        p%refl_smooth = 20.0
        p%yn_facies = .true.
        p%yn_rgt = .true.
        p%yn_salt = .true.
        p%nsalt = 3
        p%fwidth = 2
        p%seed = 112
        p%noise_level = 0.2
        call p%generate
        vp = rescale(p%facies, [1000.0, 3000.0])
        call output_array(vp, './model/vp2_2.bin')
        call output_array(p%fault, './model/fault2_2.bin')
        call output_array(p%rgt, './model/rgt2_2.bin')
        p%image = p%image/maxval(p%image)
        call output_array(p%image, './model/img2_2.bin')

    end block

    block

        real, allocatable, dimension(:, :, :) :: vp
        integer :: n1, n2, n3
        type(rgm3) :: p

        n1 = 100
        n2 = 200
        n3 = 300

        call make_directory('./model')

        p%n1 = n1
        p%n2 = n2
        p%n3 = n3
        p%nf = 17
        p%yn_regular_fault = .true.
        p%nl = 40
        p%unconf = 1
        p%unconf_z = 0.5
        p%disp = [10, 20]
        p%lwv = 0.2
        p%dip = [rand(range=[100.0, 120.0], seed=123), rand(range=[60.0, 80.0], seed=567)]
        p%strike = [rand(range=[100.0, 120.0], seed=123), rand(range=[160.0, 180.0], seed=567)]
        p%disp = [10.0, -10.0]
        p%secondary_refl_height_ratio = 0.1
        p%refl_smooth = 10.0
        p%yn_facies = .true.
        p%yn_rgt = .true.
        p%fwidth = 3
        p%seed = 324343
        call p%generate
        vp = rescale(p%facies, [1000.0, 3000.0])
        call output_array(vp, './model/vp3_1.bin')
        call output_array(p%fault, './model/fault3_1.bin')
        call output_array(p%rgt, './model/rgt3_1.bin')
        call output_array(p%image, './model/img3_1.bin')

        p%n1 = n1
        p%n2 = n2
        p%n3 = n3
        p%nf = 4
        p%yn_regular_fault = .false.
        p%nl = 30
        p%unconf = 1
        p%unconf_z = 0.25
        p%disp = [10, 20]
        p%lwv = 0.1
        p%secondary_refl_height_ratio = 0.1
        p%refl_smooth = 10.0
        p%yn_facies = .true.
        p%yn_rgt = .true.
        p%fwidth = 2
        p%seed = 111
        p%noise_level = 0.1
        call p%generate
        vp = rescale(p%facies, [1000.0, 3000.0])
        call output_array(vp, './model/vp3_2.bin')
        call output_array(p%fault, './model/fault3_2.bin')
        call output_array(p%rgt, './model/rgt3_2.bin')
        call output_array(p%image, './model/img3_2.bin')

    end block

end program test
