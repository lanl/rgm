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
! Examples of channel/canyon/drainage-type unconformities (RGM v2.0):
! the unconformity surface topography is the erosional depth map of a
! geomorphological simulation selected by unconf_shape.
!

program test

    use libflit
    use librgm

    implicit none

    character(len=24), dimension(1:4) :: shapes
    integer :: i

!    block
!
!        type(meandering_channel) :: mc
!
!        mc%nx = 251
!        mc%ny = 241
!        mc%nz = 64
!        ! The length/W ratio sets how many meander bends span the
!        ! domain. The migration width W is capped at its validated
!        ! value so that longer channels develop more bends of the
!        ! same shape (the migration dynamics are width-invariant),
!        ! while shorter channels shrink W proportionally to keep
!        ! >= ~40 centerline nodes; the channel is rendered with
!        ! W_render so that the on-grid width is approximately the
!        ! requested fraction of the lateral extent
!        mc%length = 25000
!        mc%W = 0.02*min(mc%length, 25000.0)
!        mc%W_render = 0.04*mc%length
!        mc%W_shape = 3
!        mc%n_bends = max(5, nint(30.0*mc%length/25000.0))
!        mc%n_iter = nint(1500*1.0)
!        mc%terrain_bg = 0.2*(mc%nz - 1.0)
!        !                mc%seed = seed*61 + isurf
!        mc%orient = irand(range=[0, 3], seed=1)
!        !                mc%yn_depth_only = .true.
!        call mc%generate
!
!        call output_array(mc%depth_map, './c.bin')
!
!    end block
!    stop

    shapes = ['meander_channel ', 'meander_canyon  ', 'drainage_channel', 'drainage_canyon ']

    ! ==================================================================
    ! 3D models with each unconformity shape

    do i = 1, 4

        block

            type(rgm3_curved) :: p

            p%n1 = 151
            p%n2 = 201
            p%n3 = 251
            p%nl = 30
            p%lwv = 0.4
            p%lwh = 0.3
            p%refl_shape = 'perlin'
            p%refl_shape_top = 'perlin'
            p%refl_smooth = 3
            p%refl_smooth_top = 3
            p%refl_height = [0.0, 35.0]
            p%refl_height_top = [0.0, 4.0]
            p%nf = 4
            p%dip = [60.0, 120.0]
            p%disp = [4.0, 8.0]
            p%delta_dip = [0.0, 10.0]
            p%unconf = 1
            p%unconf_z = [0.2, 0.3]
            p%unconf_height = [15.0, 45.0]
            p%unconf_shape = shapes(i)
            p%unconf_channel_width = [0.03, 0.05]
            p%unconf_topo = 0.2
            if (shapes(i) == 'drainage_channel') then
                p%unconf_channel_density = [0.05, 0.15]
            else
                p%unconf_channel_density = [0.02, 0.08]
            end if
            if (shapes(i) == 'meander_channel') then
                p%unconf_channel_length = 25000
            else
                p%unconf_channel_length = 15000
            end if
            p%noise_level = 0.01
            p%psf_sigma = [10.0, 0.5, 0.5]
            p%yn_rgt = .false.
            p%yn_facies = .false.
            p%seed = i
            call p%generate
            call output_array(p%vp, './example_3d_vp_unconf_'//tidy(shapes(i))//'_151x201x251.bin')
            call output_array(p%image, './example_3d_image_unconf_'//tidy(shapes(i))//'_151x201x251.bin')
            ! call output_array(p%rgt, './example_3d_rgt_unconf_'//tidy(shapes(i))//'_151x201x251.bin')
            ! call output_array(p%facies, './example_3d_facies_unconf_'//tidy(shapes(i))//'_151x201x251.bin')

        end block

        print *, '3D ', tidy(shapes(i)), ' done'

    end do
    stop

    ! ==================================================================
    ! 2D models with each unconformity shape

    do i = 1, 4

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
            p%delta_dip = [0.0, 10.0]
            p%unconf = 1
            p%unconf_z = [0.2, 0.3]
            p%unconf_height = [15.0, 25.0]
            p%unconf_shape = shapes(i)
            p%unconf_channel_width = [0.04, 0.08]
            p%unconf_topo = 0.2
            p%noise_level = 0.01
            p%psf_sigma = [10.0, 1.0]
            p%yn_rgt = .true.
            !p%seed = 202020 + i
            call p%generate
            call output_array(p%vp, './example_2d_vp_unconf_'//tidy(shapes(i))//'_151x301.bin')
            call output_array(p%image, './example_2d_image_unconf_'//tidy(shapes(i))//'_151x301.bin')
            call output_array(p%rgt, './example_2d_rgt_unconf_'//tidy(shapes(i))//'_151x301.bin')

        end block

        print *, '2D ', tidy(shapes(i)), ' done'

    end do

    ! ==================================================================
    ! 3D model with two channel-type unconformities + faults + salt

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
        p%nf = 3
        p%disp = [4.0, 8.0]
        p%unconf = 2
        p%unconf_z = [0.1, 0.35]
        p%unconf_height = [10.0, 20.0]
        p%unconf_shape = 'meander_channel'
        p%unconf_topo = 0.2
        p%yn_salt = .true.
        p%nsalt = 1
        p%salt_radius = [25.0, 35.0]
        p%salt_top_z = [0.5, 0.7]
        p%noise_level = 0.01
        p%psf_sigma = [10.0, 1.0, 1.0]
        p%yn_rgt = .true.
        p%seed = 303030
        call p%generate
        call output_array(p%vp, './example_3d_vp_unconf_channel_salt_151x201x251.bin')
        call output_array(p%image, './example_3d_image_unconf_channel_salt_151x201x251.bin')
        call output_array(p%salt, './example_3d_salt_unconf_channel_salt_151x201x251.bin')

        print *, '3D channel+salt done'

    end block

end program test
