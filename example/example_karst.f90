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
! Examples of karst cave systems (RGM v2.0): tube-network cave systems
! inserted into the medium parameter models like salt bodies, with
! the karst mask output in %karst and labels zeroed inside the caves.
!

program test

    use libflit
    use librgm

    implicit none

    ! ==================================================================
    ! 3D karst with faults

    block

        type(rgm3_curved) :: p

        p%n1 = 151
        p%n2 = 201
        p%n3 = 251
        p%nl = 30
        p%lwv = 0.4
        p%lwh = 0.2
        p%refl_shape = 'perlin'
        p%refl_shape_top = 'perlin'
        p%refl_smooth = 3
        p%refl_smooth_top = 3
        p%refl_height = [100.0, 150.0]
        p%refl_height_top = [0.0, 4.0]
        p%nf = 3
        p%dip = [60.0, 120.0]
        p%disp = [4.0, 8.0]
        p%yn_karst = .true.
        p%karst_z = [0.45, 0.85]
        p%karst_npassage = 28
        p%karst_connect = 0.3
        p%karst_nctrl = 18
        p%noise_level = 0.01
        p%psf_sigma = [10.0, 1.0, 1.0]
        p%yn_rgt = .true.
        p%yn_facies = .true.
        p%seed = 404040
        call p%generate
        call output_array(p%vp, './example_3d_vp_karst_1_151x201x251.bin')
        call output_array(p%image, './example_3d_image_karst_1_151x201x251.bin')
        call output_array(p%rgt, './example_3d_rgt_karst_1_151x201x251.bin')
        call output_array(p%karst, './example_3d_karst_1_151x201x251.bin')

        print *, '3D karst done'

    end block

    ! ==================================================================
    ! 3D karst + salt + channel unconformity, elastic

    block

        type(rgm3_curved) :: p

        p%n1 = 151
        p%n2 = 201
        p%n3 = 251
        p%nl = 30
        p%lwv = 0.4
        p%lwh = 0.2
        p%refl_shape = 'perlin'
        p%refl_shape_top = 'perlin'
        p%refl_smooth = 3
        p%refl_smooth_top = 3
        p%refl_height = [0.0, 15.0]
        p%refl_height_top = [0.0, 4.0]
        p%nf = 2
        p%disp = [4.0, 8.0]
        p%unconf = 1
        p%unconf_z = [0.15, 0.25]
        p%unconf_shape = 'meander_channel'
        p%yn_salt = .true.
        p%nsalt = 1
        p%salt_radius = [25.0, 35.0]
        p%salt_top_z = [0.5, 0.7]
        p%yn_karst = .true.
        p%karst_z = [0.35, 0.75]
        p%yn_elastic = .true.
        p%noise_level = 0.01
        p%psf_sigma = [10.0, 1.0, 1.0]
        p%seed = 505050
        call p%generate
        call output_array(p%vp, './example_3d_vp_karst_2_151x201x251.bin')
        call output_array(p%vs, './example_3d_vs_karst_2_151x201x251.bin')
        call output_array(p%rho, './example_3d_rho_karst_2_151x201x251.bin')
        call output_array(p%image_pp, './example_3d_image_pp_karst_2_151x201x251.bin')
        call output_array(p%karst, './example_3d_karst_2_151x201x251.bin')
        call output_array(p%salt, './example_3d_salt_karst_2_151x201x251.bin')

        print *, '3D karst+salt+unconf elastic done'

    end block

    ! ==================================================================
    ! 2D karst with faults

    block

        type(rgm2_curved) :: p

        p%n1 = 151
        p%n2 = 301
        p%nl = 30
        p%lwv = 0.4
        p%lwh = 0.2
        p%refl_shape = 'perlin'
        p%refl_shape_top = 'perlin'
        p%refl_smooth = 3
        p%refl_smooth_top = 3
        p%refl_height = [0.0, 12.0]
        p%refl_height_top = [0.0, 4.0]
        p%nf = 3
        p%dip = [60.0, 120.0]
        p%disp = [5.0, 10.0]
        p%yn_karst = .true.
        p%karst_z = [0.4, 0.85]
        p%karst_npassage = 6
        p%noise_level = 0.01
        p%psf_sigma = [10.0, 1.0]
        p%yn_rgt = .true.
        p%seed = 606060
        call p%generate
        call output_array(p%vp, './example_2d_vp_karst_1_151x301.bin')
        call output_array(p%image, './example_2d_image_karst_1_151x301.bin')
        call output_array(p%rgt, './example_2d_rgt_karst_1_151x301.bin')
        call output_array(p%karst, './example_2d_karst_1_151x301.bin')

        print *, '2D karst done'

    end block

    ! ==================================================================
    ! 2D karst below an unconformity, karst_before_unconf = .false.

    block

        type(rgm2_curved) :: p

        p%n1 = 151
        p%n2 = 301
        p%nl = 30
        p%lwv = 0.4
        p%lwh = 0.2
        p%refl_shape = 'perlin'
        p%refl_shape_top = 'perlin'
        p%refl_smooth = 3
        p%refl_smooth_top = 3
        p%refl_height = [100.0, 150.0]
        p%refl_height_top = [0.0, 4.0]
        p%nf = 2
        p%disp = [5.0, 10.0]
        p%unconf = 1
        p%unconf_z = [0.15, 0.25]
        p%yn_karst = .true.
        p%karst_z = [0.4, 0.85]
        p%karst_npassage = 3
        p%karst_before_unconf = .false.
        p%noise_level = 0.01
        p%psf_sigma = [10.0, 1.0]
        p%seed = 707070
        call p%generate
        call output_array(p%vp, './example_2d_vp_karst_2_151x301.bin')
        call output_array(p%image, './example_2d_image_karst_2_151x301.bin')
        call output_array(p%karst, './example_2d_karst_2_151x301.bin')

        print *, '2D karst+unconf done'

    end block

end program test
