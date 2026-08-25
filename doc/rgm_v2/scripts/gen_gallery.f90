!
! Generate a batch of 16 random 3D models (RGM v2.0) covering the new
! feature matrix - strike-varying faults, spatially varying and decaying
! displacement, channel/canyon/drainage unconformities, and karst - for
! the 4x4 gallery figures of the manuscript.
!
program gen_gallery

    use libflit
    use librgm

    implicit none

    integer :: i

    do i = 1, 16

        block

            type(rgm3_curved) :: p

            p%n1 = 128
            p%n2 = 160
            p%n3 = 200
            p%nl = 24
            p%lwv = 0.4
            p%lwh = 0.2
            p%refl_shape = 'perlin'
            p%refl_shape_top = 'perlin'
            p%refl_smooth = 3
            p%refl_smooth_top = 3
            p%refl_height = [0.0, 10.0]
            p%refl_height_top = [0.0, 3.0]
            p%noise_level = 0.01
            p%psf_sigma = [8.0, 1.0, 1.0]
            p%seed = 9000 + i

            select case (i)
                case (1)
                    ! strike-varying, through-going faults
                    p%nf = 4
                    p%dip = [60.0, 120.0]
                    p%disp = [8.0, 15.0]
                    p%delta_strike = [15.0, 25.0]
                    p%refl_shape = 'gaussian'
                    p%refl_smooth = 0
                    p%refl_sigma2 = [35.0, 60.0]
                    p%refl_sigma3 = [35.0, 60.0]
                    p%refl_height = [0.0, 45.0]
                case (2)
                    ! strike-varying, strongly listric faults
                    p%nf = 3
                    p%dip = [55.0, 125.0]
                    p%disp = [8.0, 15.0]
                    p%delta_dip = [10.0, 25.0]
                    p%delta_strike = [20.0, 30.0]
                case (3)
                    ! strike variation + spatially varying displacement
                    p%nf = 5
                    p%dip = [60.0, 120.0]
                    p%disp = [10.0, 18.0]
                    p%delta_strike = [15.0, 25.0]
                    p%yn_vary_disp = .true.
                    p%disp_radius_strike = [0.4, 0.6]
                    p%disp_radius_dip = [0.5, 0.8]
                    p%refl_height = [0.0, 30.0]
                case (4)
                    ! + fault-normal displacement decay (drag/rollover)
                    p%nf = 4
                    p%dip = [60.0, 120.0]
                    p%disp = [10.0, 18.0]
                    p%delta_strike = [15.0, 25.0]
                    p%yn_vary_disp = .true.
                    p%disp_radius_strike = [0.4, 0.6]
                    p%disp_radius_dip = [0.5, 0.8]
                    p%yn_disp_decay = .true.
                    p%disp_decay_width = [0.3, 0.5]
                case (5)
                    ! meander-channel unconformity + strike-varying faults
                    p%nf = 3
                    p%dip = [60.0, 120.0]
                    p%disp = [5.0, 10.0]
                    p%delta_strike = [15.0, 25.0]
                    p%unconf = 1
                    p%unconf_z = [0.2, 0.3]
                    p%unconf_height = [10.0, 18.0]
                    p%unconf_shape = 'meander_channel'
                    p%refl_shape = 'gaussian'
                    p%refl_smooth = 0
                    p%refl_sigma2 = [35.0, 60.0]
                    p%refl_sigma3 = [35.0, 60.0]
                    p%refl_height = [0.0, 100.0]
                case (6)
                    ! meander-canyon unconformity + varying displacement
                    p%nf = 3
                    p%dip = [60.0, 120.0]
                    p%disp = [8.0, 14.0]
                    p%yn_vary_disp = .true.
                    p%disp_radius_strike = [0.4, 0.6]
                    p%disp_radius_dip = [0.5, 0.8]
                    p%unconf = 1
                    p%unconf_z = [0.2, 0.3]
                    p%unconf_height = [12.0, 20.0]
                    p%unconf_shape = 'meander_canyon'
                case (7)
                    ! drainage-channel unconformity
                    p%nf = 4
                    p%dip = [60.0, 120.0]
                    p%disp = [5.0, 10.0]
                    p%unconf = 1
                    p%unconf_z = [0.2, 0.3]
                    p%unconf_height = [10.0, 18.0]
                    p%unconf_shape = 'drainage_channel'
                    p%refl_height = [0.0, 25.0]
                case (8)
                    ! drainage-canyon unconformity + decay
                    p%nf = 3
                    p%dip = [60.0, 120.0]
                    p%disp = [8.0, 14.0]
                    p%yn_vary_disp = .true.
                    p%disp_radius_strike = [0.4, 0.6]
                    p%disp_radius_dip = [0.5, 0.8]
                    p%yn_disp_decay = .true.
                    p%disp_decay_width = [0.3, 0.5]
                    p%unconf = 1
                    p%unconf_z = [0.2, 0.3]
                    p%unconf_height = [12.0, 20.0]
                    p%unconf_shape = 'drainage_canyon'
                    p%refl_height = [0.0, 100.0]
                case (9)
                    ! karst + strike-varying faults
                    p%nf = 3
                    p%dip = [60.0, 120.0]
                    p%disp = [4.0, 8.0]
                    p%delta_strike = [15.0, 25.0]
                    p%yn_karst = .true.
                    p%karst_z = [0.45, 0.85]
                    p%refl_shape = 'gaussian'
                    p%refl_smooth = 0
                    p%refl_sigma2 = [35.0, 60.0]
                    p%refl_sigma3 = [35.0, 60.0]
                    p%refl_height = [0.0, 40.0]
                case (10)
                    ! high-connectivity karst + varying displacement
                    p%nf = 3
                    p%dip = [60.0, 120.0]
                    p%disp = [6.0, 12.0]
                    p%yn_vary_disp = .true.
                    p%disp_radius_strike = [0.4, 0.6]
                    p%disp_radius_dip = [0.5, 0.8]
                    p%yn_karst = .true.
                    p%karst_z = [0.45, 0.85]
                    p%karst_npassage = 16
                    p%karst_connect = 0.95
                case (11)
                    ! karst formed after a meander-channel unconformity
                    p%nf = 2
                    p%dip = [65.0, 115.0]
                    p%disp = [4.0, 8.0]
                    p%unconf = 1
                    p%unconf_z = [0.18, 0.28]
                    p%unconf_height = [10.0, 16.0]
                    p%unconf_shape = 'meander_channel'
                    p%yn_karst = .true.
                    p%karst_z = [0.45, 0.85]
                    p%karst_before_unconf = .false.
                    p%refl_height = [0.0, 25.0]
                case (12)
                    ! salt + karst + faults
                    p%nf = 2
                    p%dip = [60.0, 120.0]
                    p%disp = [4.0, 8.0]
                    p%yn_salt = .true.
                    p%nsalt = 1
                    p%salt_radius = [22.0, 30.0]
                    p%salt_top_z = [0.5, 0.7]
                    p%yn_karst = .true.
                    p%karst_z = [0.4, 0.75]
                case (13)
                    ! two meander-channel unconformities + varying disp
                    p%nf = 3
                    p%dip = [60.0, 120.0]
                    p%disp = [8.0, 14.0]
                    p%yn_vary_disp = .true.
                    p%disp_radius_strike = [0.4, 0.6]
                    p%disp_radius_dip = [0.5, 0.8]
                    p%unconf = 2
                    p%unconf_z = [0.12, 0.35]
                    p%unconf_height = [8.0, 15.0]
                    p%unconf_shape = 'meander_channel'
                    p%refl_height = [0.0, 100.0]
                case (14)
                    ! dense thin faults + drainage-channel unconformity
                    p%nf = 6
                    p%dip = [65.0, 115.0]
                    p%disp = [4.0, 8.0]
                    p%delta_strike = [10.0, 20.0]
                    p%unconf = 1
                    p%unconf_z = [0.2, 0.3]
                    p%unconf_height = [10.0, 16.0]
                    p%unconf_shape = 'drainage_channel'
                case (15)
                    ! meander canyon + strike variation + decay
                    p%nf = 3
                    p%dip = [55.0, 125.0]
                    p%disp = [8.0, 14.0]
                    p%delta_strike = [15.0, 25.0]
                    p%yn_vary_disp = .true.
                    p%disp_radius_strike = [0.4, 0.6]
                    p%disp_radius_dip = [0.5, 0.8]
                    p%yn_disp_decay = .true.
                    p%disp_decay_width = [0.3, 0.5]
                    p%unconf = 1
                    p%unconf_z = [0.2, 0.3]
                    p%unconf_height = [12.0, 20.0]
                    p%unconf_shape = 'meander_canyon'
                    p%refl_shape = 'gaussian'
                    p%refl_smooth = 0
                    p%refl_sigma2 = [35.0, 60.0]
                    p%refl_sigma3 = [35.0, 60.0]
                    p%refl_height = [0.0, 100.0]
                case (16)
                    ! combined: channel unconformity + karst + all fault features
                    p%nf = 3
                    p%dip = [60.0, 120.0]
                    p%disp = [8.0, 14.0]
                    p%delta_strike = [15.0, 25.0]
                    p%yn_vary_disp = .true.
                    p%disp_radius_strike = [0.4, 0.6]
                    p%disp_radius_dip = [0.5, 0.8]
                    p%yn_disp_decay = .true.
                    p%disp_decay_width = [0.3, 0.5]
                    p%unconf = 1
                    p%unconf_z = [0.18, 0.28]
                    p%unconf_height = [10.0, 16.0]
                    p%unconf_shape = 'meander_channel'
                    p%yn_karst = .true.
                    p%karst_z = [0.45, 0.85]
                    p%refl_height = [0.0, 30.0]
            end select

            call p%generate
            call output_array(p%vp, './gallery_vp_'//num2str(i)//'_128x160x200.bin')
            call output_array(p%image, './gallery_img_'//num2str(i)//'_128x160x200.bin')
            call output_array(p%fault_strike, './gallery_fstrike_'//num2str(i)//'_128x160x200.bin')
            call output_array(p%fault_disp, './gallery_fdisp_'//num2str(i)//'_128x160x200.bin')
            print *, 'gallery model ', i, ' done'

        end block

    end do

end program gen_gallery
