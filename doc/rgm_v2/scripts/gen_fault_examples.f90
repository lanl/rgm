!
! Generate the 3D fault example models of the manuscript: strike-varying
! through-going faults, spatially varying displacement, and fault-normal
! displacement decay, sharing one fault framework and a fixed seed.
!
program gen_fault_examples

    use libflit
    use librgm

    implicit none

    type(rgm3_curved) :: p

    p%n1 = 171
    p%n2 = 251
    p%n3 = 251
    p%lwv = 0.4
    p%lwh = 0.2
    p%refl_shape = 'gaussian'
    p%refl_shape_top = 'perlin'
    p%refl_smooth = 0
    p%refl_smooth_top = 0
    p%refl_sigma2 = [50, 70]
    p%refl_sigma3 = [50, 70]
    p%refl_height = [0.0, 70.0]
    p%refl_height_top = [0.0, 4.0]
    p%nl = 30
    p%dip = [60.0, 120.0]
    p%disp = [15.0, 25.0]
    p%delta_dip = [0.0, 15.0]
    p%delta_strike = [15.0, 25.0]
    p%disp_radius_strike = [0.4, 0.6]
    p%disp_radius_dip = [0.5, 0.8]
    p%noise_level = 0.01
    p%psf_sigma = [10.0, 1.0, 1.0]

    ! Strike-varying, through-going faults (constant displacement)
    p%nf = 6
    p%seed = 111
    p%yn_vary_disp = .false.
    p%yn_disp_decay = .false.
    call p%generate
    call output_array(p%vp, './example_3d_vp_vstrike_1.bin')
    call output_array(p%image, './example_3d_image_vstrike_1.bin')
    call output_array(p%fault, './example_3d_fault_vstrike_1.bin')
    call output_array(p%fault_strike, './example_3d_fstrike_vstrike_1.bin')
    print *, 'vstrike done'

    stop

    ! Spatially varying displacement (elliptical slip patch); the decay
    ! case below inherits the same fault framework and seed
    p%nf = 6
    p%seed = 202513
    p%yn_vary_disp = .true.
    p%yn_disp_decay = .false.
    call p%generate
    call output_array(p%vp, './example_3d_vp_vdisp_1.bin')
    call output_array(p%image, './example_3d_image_vdisp_1.bin')
    call output_array(p%fault, './example_3d_fault_vdisp_1.bin')
    call output_array(p%fault_disp, './example_3d_fdisp_vdisp_1.bin')
    print *, 'vdisp done'

    ! Additionally decay the displacement away from the faults
    p%yn_vary_disp = .true.
    p%yn_disp_decay = .true.
    p%disp_decay_width = [0.3, 0.5]
    call p%generate
    call output_array(p%vp, './example_3d_vp_decay_1.bin')
    call output_array(p%image, './example_3d_image_decay_1.bin')
    call output_array(p%fault_disp, './example_3d_fdisp_decay_1.bin')
    print *, 'decay done'

end program gen_fault_examples
