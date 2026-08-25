!
! Generate additional data for the manuscript figures:
!  - erosional depth maps of the four geomorphological unconformity
!    types at the 3D mapping settings (for topography renderings)
!  - a high-connectivity karst example (rgm3_curved) with vp + mask
!
program gen_paper_extras

    use libflit
    use librgm

    implicit none

    block

        type(meandering_channel) :: mc

        mc%nx = 251
        mc%ny = 201
        mc%nz = 64
        mc%length = 25000.0
        mc%W = 0.025*25000.0
        mc%W_render = 0.05*25000.0
        mc%W_shape = 3
        mc%n_bends = 30
        mc%n_iter = 1500
        mc%terrain_bg = 0.25*(mc%nz - 1.0)
        mc%seed = 12345
        mc%orient = 0
        mc%yn_depth_only = .true.
        call mc%generate
        call output_array(mc%depth_map, './topo_meander_channel_201x251.bin')
        print *, 'meander channel done'

    end block

    block

        type(meandering_canyon) :: my
        integer :: nlev

        my%nx = 251
        my%ny = 201
        my%nz = 0
        my%length = 15000.0
        my%W = 0.02*15000.0
        my%W_canyon = 0.5*0.05*15000.0
        my%n_bends = 20
        my%n_iter = 4000
        nlev = max(2, nint(my%n_iter/(my%save_every*1.0)))
        my%terrain_bg = 0.25*(nlev - 1.0)
        my%seed = 12345
        my%orient = 0
        my%yn_depth_only = .true.
        call my%generate
        call output_array(my%depth_map, './topo_meander_canyon_201x251.bin')
        print *, 'meander canyon done'

    end block

    block

        type(drainage_channel) :: dc

        dc%nx = 251
        dc%ny = 201
        dc%nz = 64
        dc%W_max = 0.5*0.05*0.5*(201.0 + 251.0)
        dc%D_max = 48.0
        dc%channel_frac = 0.05
        dc%terrain_bg = 0.25*(dc%nz - 1.0)
        dc%seed = 12345
        dc%orient = 0
        dc%yn_depth_only = .true.
        call dc%generate
        call output_array(dc%depth_map, './topo_drainage_channel_201x251.bin')
        print *, 'drainage channel done'

    end block

    block

        type(drainage_canyon) :: dy

        dy%nx = 251
        dy%ny = 201
        dy%nz = 64
        dy%W_max = 0.05*0.5*(201.0 + 251.0)
        dy%channel_frac = 0.25*0.05
        dy%terrain_bg = 0.25*(dy%nz - 1.0)
        dy%seed = 12345
        dy%orient = 0
        dy%yn_depth_only = .true.
        call dy%generate
        call output_array(dy%depth_map, './topo_drainage_canyon_201x251.bin')
        print *, 'drainage canyon done'

    end block

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
        p%dip = [60.0, 120.0]
        p%disp = [4.0, 8.0]
        p%yn_karst = .true.
        p%karst_z = [0.45, 0.85]
        p%karst_npassage = 28
        p%karst_connect = 0.95
        p%noise_level = 0.01
        p%psf_sigma = [10.0, 1.0, 1.0]
        p%seed = 414141
        call p%generate
        call output_array(p%vp, './karst_connected_vp_151x201x251.bin')
        call output_array(p%karst, './karst_connected_mask_151x201x251.bin')
        print *, 'connected karst done'

    end block

end program gen_paper_extras
