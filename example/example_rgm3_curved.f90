
program test

    use libflit
    use librgm
    use omp_lib

    implicit none

    ! ==================================================================
    ! Scalability test

    block

        type(rgm3_curved) :: p
        integer :: i, j
        real, allocatable, dimension(:, :) :: ctime
        integer, allocatable, dimension(:) :: ns, nf
        double precision :: tbeg, tend

        ns = regspace(50.0, 50.0, 500.0)
        nf = regspace(1, 1, 20)

        ctime = zeros(size(ns), size(nf))

        do i = 1, size(ns)
            do j = 1, size(nf)

                p%n1 = ns(i)
                p%n2 = ns(i)
                p%n3 = ns(i)
                p%lwv = 0.3
                p%lwh = 0.3
                p%refl_shape = 'cauchy'
                p%refl_shape_top = 'perlin'
                p%refl_smooth_top = 3
                p%refl_slope = [rand(range=[-20.0, 20.0], seed=i), rand(range=[-20.0, 20.0], seed=i + 100)]
                p%refl_slope_top = [rand(range=[-5.0, 5.0], seed=i + 200), rand(range=[-5.0, 5.0], seed=i + 300)]
                p%refl_height = [0.0, 0.3*ns(i)]
                p%refl_height_top = [0.0, 0.05*ns(i)]
                p%disp = [5.0, 10.0]
                p%nl = 35
                p%nf = nf(j)
                p%seed = i
                p%delta_dip = [0.0, rand(range=[0.0, 30.0], seed=i + 200)]

                tbeg = omp_get_wtime()
                call p%generate
                tend = omp_get_wtime()

                ctime(i, j) = tend - tbeg

                print *, i, j

            end do
        end do

        call output_array(ctime, './ctime3.bin', transp=.true.)
        call output_array(ctime, './ctime3.txt', transp=.true., ascii=.true.)

    end block

    ! ==================================================================
    ! Types of faults

    block

        type(rgm3_curved) :: p

        p%n1 = 171
        p%n2 = 251
        p%n3 = 251
        p%lwv = 0.4
        p%lwh = 0.3
        p%refl_shape = 'cauchy'
        p%ng = 5
        p%rotate_fold = .true.
        p%refl_sigma2 = [30.0, 60.0]
        p%refl_shape_top = 'perlin'
        p%refl_smooth = 0
        p%refl_smooth_top = 3
        p%refl_slope = [20.0, -10.0]
        p%refl_slope_top = [-2.0, 3.0]
        p%refl_height = [0.0, 50.0]
        p%refl_height_top = [0.0, 4.0]
        p%nl = 35
        p%nf = 2
        p%noise_level = 0.025
        p%disp = [5.0, 12.0]
        p%delta_dip = [0.0, 0.0]
        p%psf_sigma = [10.0, 1.0, 1.0]
        p%seed = 1234567
        call p%generate
        call output_array(p%image, './example_3d_image_type_1.bin')
        call output_array(p%fault_dip, './example_3d_fdip_type_1.bin')

        p%delta_dip = [30.0, 40.0]
        p%seed = 1234567
        call p%generate
        call output_array(p%image, './example_3d_image_type_2.bin')
        call output_array(p%fault_dip, './example_3d_fdip_type_2.bin')

    end block

    ! ==================================================================
    ! Unfaulted models

    block

        type(rgm3_curved) :: p

        p%n1 = 171
        p%n2 = 251
        p%n3 = 251
        p%lwv = 0.4
        p%lwh = 0.3
        p%refl_shape = 'perlin'
        p%refl_shape_top = 'perlin'
        p%refl_smooth = 3
        p%refl_smooth_top = 3
        p%refl_slope = [20.0, -10.0]
        p%refl_slope_top = [-2.0, 3.0]
        p%refl_height = [0.0, 50.0]
        p%refl_height_top = [0.0, 4.0]
        p%nl = 30
        p%nf = 0
        p%noise_level = 0.025
        p%psf_sigma = [10.0, 1.0, 1.0]
        p%seed = 123
        p%yn_rgt = .true.
        call p%generate
        call output_array(p%vp, './example_3d_vp_unfaulted_1.bin')
        call output_array(p%rgt, './example_3d_rgt_unfaulted_1.bin')
        call output_array(p%image, './example_3d_image_unfaulted_1.bin')

    end block

    block

        type(rgm3_curved) :: p

        p%n1 = 171
        p%n2 = 251
        p%n3 = 251
        p%lwv = 0.4
        p%lwh = 0.2
        p%refl_shape = 'cauchy'
        p%ng = 3
        p%refl_sigma2 = [40.0, 60.0]
        p%refl_sigma3 = [30.0, 50.0]
        p%refl_shape_top = 'perlin'
        p%refl_smooth = 0
        p%refl_smooth_top = 4
        p%refl_slope = [20.0, -20.0]
        p%refl_slope_top = [-2.0, -5.0]
        p%refl_height = [0.0, 70.0]
        p%refl_height_top = [0.0, 5.0]
        p%nl = 25
        p%nf = 0
        p%noise_level = 0.02
        p%psf_sigma = [10.0, 1.0, 1.0]
        p%seed = 456
        p%yn_rgt = .true.
        call p%generate
        call output_array(p%vp, './example_3d_vp_unfaulted_2.bin')
        call output_array(p%rgt, './example_3d_rgt_unfaulted_2.bin')
        call output_array(p%image, './example_3d_image_unfaulted_2.bin')

    end block

    block

        type(rgm3_curved) :: p

        p%n1 = 171
        p%n2 = 251
        p%n3 = 251
        p%lwv = 0.5
        p%lwh = 0.3
        p%refl_shape = 'random'
        p%refl_shape_top = 'perlin'
        p%refl_smooth = 40
        p%refl_smooth_top = 2
        p%refl_slope = [-20.0, 10.0]
        p%refl_slope_top = [5.0, 0.0]
        p%refl_height = [0.0, 50.0]
        p%refl_height_top = [0.0, 4.0]
        p%nl = 35
        p%nf = 0
        p%noise_level = 0.05
        p%psf_sigma = [10.0, 1.0, 1.0]
        p%seed = 888
        p%yn_rgt = .true.
        call p%generate
        call output_array(p%vp, './example_3d_vp_unfaulted_3.bin')
        call output_array(p%rgt, './example_3d_rgt_unfaulted_3.bin')
        call output_array(p%image, './example_3d_image_unfaulted_3.bin')

    end block

    block

        type(rgm3_curved) :: p

        p%n1 = 171
        p%n2 = 251
        p%n3 = 251
        p%lwv = 0.4
        p%lwh = 0.4
        p%refl_shape = 'gaussian'
        p%refl_shape_top = 'same'
        p%refl_smooth = 0
        p%refl_smooth_top = 0
        p%ng = 3
        p%refl_sigma2 = [40.0, 60.0]
        p%refl_sigma3 = [40.0, 60.0]
        p%refl_slope = [10.0, -10.0]
        p%refl_slope_top = [20.0, -20.0]
        p%refl_height = [0.0, 50.0]
        p%refl_height_top = [0.0, 0.0]
        p%nl = 40
        p%nf = 0
        p%noise_level = 0.02
        p%psf_sigma = [10.0, 1.0, 1.0]
        p%seed = 789
        p%yn_rgt = .true.
        call p%generate
        call output_array(p%vp, './example_3d_vp_unfaulted_4.bin')
        call output_array(p%rgt, './example_3d_rgt_unfaulted_4.bin')
        call output_array(p%image, './example_3d_image_unfaulted_4.bin')

    end block

    ! ==================================================================
    ! Faulted models

    block

        type(rgm3_curved) :: p

        p%n1 = 171
        p%n2 = 251
        p%n3 = 251
        p%lwv = 0.4
        p%lwh = 0.3
        p%refl_shape = 'perlin'
        p%refl_shape_top = 'perlin'
        p%refl_smooth = 3
        p%refl_smooth_top = 3
        p%refl_slope = [20.0, -10.0]
        p%refl_slope_top = [-2.0, 3.0]
        p%refl_height = [0.0, 50.0]
        p%refl_height_top = [0.0, 4.0]
        p%nl = 30
        p%nf = 7
        p%noise_level = 0.025
        p%disp = [5.0, 12.0]
        p%delta_dip = [0.0, 10.0]
        p%psf_sigma = [10.0, 1.0, 1.0]
        p%seed = 123
        p%yn_rgt = .true.
        call p%generate
        call output_array(p%vp, './example_3d_vp_faulted_1.bin')
        call output_array(p%rgt, './example_3d_rgt_faulted_1.bin')
        call output_array(p%image, './example_3d_image_faulted_1.bin')
        call output_array(p%fault, './example_3d_fault_1.bin')
        call output_array(p%fault_dip, './example_3d_fdip_1.bin')
        call output_array(p%fault_strike, './example_3d_fstrike_1.bin')

    end block

    block

        type(rgm3_curved) :: p

        p%n1 = 171
        p%n2 = 251
        p%n3 = 251
        p%lwv = 0.4
        p%lwh = 0.2
        p%refl_shape = 'cauchy'
        p%ng = 3
        p%refl_sigma2 = [40.0, 60.0]
        p%refl_sigma3 = [30.0, 50.0]
        p%refl_shape_top = 'perlin'
        p%refl_smooth = 0
        p%refl_smooth_top = 4
        p%refl_slope = [20.0, -20.0]
        p%refl_slope_top = [-2.0, -5.0]
        p%refl_height = [0.0, 70.0]
        p%refl_height_top = [0.0, 5.0]
        p%nl = 25
        p%nf = 12
        p%disp = [-6.0, 6.0]
        p%yn_regular_fault = .true.
        p%strike = [45.0, 130.0]
        p%delta_dip = [0.0, 10.0]
        p%noise_level = 0.02
        p%psf_sigma = [10.0, 1.0, 1.0]
        p%seed = 456
        p%yn_rgt = .true.
        call p%generate
        call output_array(p%vp, './example_3d_vp_faulted_2.bin')
        call output_array(p%rgt, './example_3d_rgt_faulted_2.bin')
        call output_array(p%image, './example_3d_image_faulted_2.bin')
        call output_array(p%fault, './example_3d_fault_2.bin')
        call output_array(p%fault_dip, './example_3d_fdip_2.bin')
        call output_array(p%fault_strike, './example_3d_fstrike_2.bin')

    end block

    block

        type(rgm3_curved) :: p

        p%n1 = 171
        p%n2 = 251
        p%n3 = 251
        p%lwv = 0.5
        p%lwh = 0.3
        p%refl_shape = 'random'
        p%refl_shape_top = 'perlin'
        p%refl_smooth = 40
        p%refl_smooth_top = 2
        p%refl_slope = [-20.0, 10.0]
        p%refl_slope_top = [5.0, 0.0]
        p%refl_height = [0.0, 50.0]
        p%refl_height_top = [0.0, 4.0]
        p%nl = 35
        p%nf = 12
        p%disp = [5.0, 10.0]
        p%delta_dip = [20.0, 40.0]
        p%noise_level = 0.05
        p%psf_sigma = [10.0, 1.0, 1.0]
        p%seed = 888
        p%yn_rgt = .true.
        call p%generate
        call output_array(p%vp, './example_3d_vp_faulted_3.bin')
        call output_array(p%rgt, './example_3d_rgt_faulted_3.bin')
        call output_array(p%image, './example_3d_image_faulted_3.bin')
        call output_array(p%fault, './example_3d_fault_3.bin')
        call output_array(p%fault_dip, './example_3d_fdip_3.bin')
        call output_array(p%fault_strike, './example_3d_fstrike_3.bin')

    end block

    block

        type(rgm3_curved) :: p

        p%n1 = 171
        p%n2 = 251
        p%n3 = 251
        p%lwv = 0.4
        p%lwh = 0.4
        p%refl_shape = 'gaussian'
        p%refl_shape_top = 'same'
        p%refl_smooth = 0
        p%refl_smooth_top = 0
        p%ng = 3
        p%refl_sigma2 = [40.0, 60.0]
        p%refl_sigma3 = [40.0, 60.0]
        p%refl_slope = [10.0, -10.0]
        p%refl_slope_top = [20.0, -20.0]
        p%refl_height = [0.0, 50.0]
        p%refl_height_top = [0.0, 0.0]
        p%nl = 40
        p%nf = 7
        p%disp = [10.0, 20.0]
        p%noise_level = 0.02
        p%psf_sigma = [10.0, 1.0, 1.0]
        p%seed = 789
        p%yn_rgt = .true.
        call p%generate
        call output_array(p%vp, './example_3d_vp_faulted_4.bin')
        call output_array(p%rgt, './example_3d_rgt_faulted_4.bin')
        call output_array(p%image, './example_3d_image_faulted_4.bin')
        call output_array(p%fault, './example_3d_fault_4.bin')
        call output_array(p%fault_dip, './example_3d_fdip_4.bin')
        call output_array(p%fault_strike, './example_3d_fstrike_4.bin')

    end block

    ! ==================================================================
    ! Salt models

    block

        type(rgm3_curved) :: p

        p%n1 = 171
        p%n2 = 251
        p%n3 = 251
        p%lwv = 0.5
        p%lwh = 0.3
        p%refl_shape = 'random'
        p%refl_shape_top = 'perlin'
        p%refl_smooth = 40
        p%refl_smooth_top = 2
        p%refl_slope = [-20.0, 10.0]
        p%refl_slope_top = [5.0, 0.0]
        p%refl_height = [0.0, 50.0]
        p%refl_height_top = [0.0, 4.0]
        p%nl = 35
        p%nf = 0
        p%noise_level = 0.05
        p%psf_sigma = [10.0, 1.0, 1.0]
        p%seed = 13579
        p%yn_salt = .true.
        p%nsalt = 3
        p%salt_radius = [15.0, 40.0]
        p%salt_top_z = [0.4, 0.6]
        call p%generate
        call output_array(p%vp, './example_3d_vp_salt_1.bin')
        call output_array(p%image, './example_3d_image_salt_1.bin')
        call output_array(p%salt, './example_3d_salt_1.bin')

    end block

    block

        type(rgm3_curved) :: p

        p%n1 = 171
        p%n2 = 251
        p%n3 = 251
        p%lwv = 0.4
        p%lwh = 0.4
        p%refl_shape = 'gaussian'
        p%refl_shape_top = 'same'
        p%refl_smooth = 0
        p%refl_smooth_top = 0
        p%ng = 3
        p%refl_sigma2 = [40.0, 60.0]
        p%refl_sigma3 = [40.0, 60.0]
        p%refl_slope = [10.0, -10.0]
        p%refl_slope_top = [20.0, -20.0]
        p%refl_height = [0.0, 50.0]
        p%refl_height_top = [0.0, 0.0]
        p%nl = 40
        p%nf = 0
        p%noise_level = 0.02
        p%psf_sigma = [10.0, 1.0, 1.0]
        p%seed = 789789
        p%yn_salt = .true.
        p%nsalt = 4
        p%salt_radius = [25.0, 55.0]
        p%salt_top_z = [0.4, 0.6]
        p%salt_radius_variation = 0.9
        call p%generate
        call output_array(p%vp, './example_3d_vp_salt_2.bin')
        call output_array(p%image, './example_3d_image_salt_2.bin')
        call output_array(p%salt, './example_3d_salt_2.bin')

    end block

    block

        type(rgm3_curved) :: p

        p%n1 = 171
        p%n2 = 251
        p%n3 = 251
        p%lwv = 0.5
        p%lwh = 0.3
        p%refl_shape = 'random'
        p%refl_shape_top = 'perlin'
        p%refl_smooth = 40
        p%refl_smooth_top = 2
        p%refl_slope = [-20.0, 10.0]
        p%refl_slope_top = [5.0, 0.0]
        p%refl_height = [0.0, 50.0]
        p%refl_height_top = [0.0, 4.0]
        p%nl = 35
        p%nf = 12
        p%disp = [5.0, 10.0]
        p%delta_dip = [20.0, 40.0]
        p%noise_level = 0.05
        p%psf_sigma = [10.0, 1.0, 1.0]
        p%seed = 1888
        p%yn_salt = .true.
        p%nsalt = 2
        p%salt_radius = [25.0, 55.0]
        p%salt_top_z = [0.4, 0.6]
        p%salt_top_height = 30
        p%salt_path_variation = 20
        call p%generate
        call output_array(p%vp, './example_3d_vp_salt_3.bin')
        call output_array(p%image, './example_3d_image_salt_3.bin')
        call output_array(p%salt, './example_3d_salt_3.bin')

    end block

    block

        type(rgm3_curved) :: p

        p%n1 = 171
        p%n2 = 251
        p%n3 = 251
        p%lwv = 0.4
        p%lwh = 0.4
        p%refl_shape = 'gaussian'
        p%refl_shape_top = 'same'
        p%refl_smooth = 0
        p%refl_smooth_top = 0
        p%ng = 3
        p%refl_sigma2 = [40.0, 60.0]
        p%refl_sigma3 = [40.0, 60.0]
        p%refl_slope = [10.0, -10.0]
        p%refl_slope_top = [20.0, -20.0]
        p%refl_height = [0.0, 50.0]
        p%refl_height_top = [0.0, 0.0]
        p%nl = 40
        p%nf = 7
        p%disp = [10.0, 20.0]
        p%noise_level = 0.02
        p%psf_sigma = [10.0, 1.0, 1.0]
        p%seed = 11789
        p%yn_salt = .true.
        p%nsalt = 5
        p%salt_radius = [15.0, 45.0]
        p%salt_top_z = [0.4, 0.7]
        p%salt_top_height = 20
        call p%generate
        call output_array(p%vp, './example_3d_vp_salt_4.bin')
        call output_array(p%image, './example_3d_image_salt_4.bin')
        call output_array(p%salt, './example_3d_salt_4.bin')

    end block

    ! ==================================================================
    ! Unconformity models

    block

        type(rgm3_curved) :: p

        p%n1 = 171
        p%n2 = 251
        p%n3 = 251
        p%lwv = 0.4
        p%lwh = 0.3
        p%refl_shape = 'perlin'
        p%refl_shape_top = 'perlin'
        p%refl_smooth = 3
        p%refl_smooth_top = 3
        p%refl_slope = [20.0, -10.0]
        p%refl_slope_top = [-2.0, 3.0]
        p%refl_height = [0.0, 50.0]
        p%refl_height_top = [0.0, 20.0]
        p%nl = 30
        p%nf = 7
        p%noise_level = 0.025
        p%disp = [5.0, 12.0]
        p%delta_dip = [0.0, 10.0]
        p%psf_sigma = [10.0, 1.0, 1.0]
        p%seed = 12357
        p%unconf = 1
        p%unconf_height = [20.0, 30.0]
        p%unconf_z = [0.1, 0.3]
        p%yn_rgt = .true.
        p%yn_salt = .true.
        p%nsalt = 2
        p%salt_radius = [25.0, 55.0]
        p%salt_top_z = [0.4, 0.6]
        p%salt_top_height = 30
        p%salt_path_variation = 20
        p%yn_facies = .true.
        call p%generate
        call output_array(p%vp, './example_3d_vp_unconf_1.bin')
        call output_array(p%rgt, './example_3d_rgt_unconf_1.bin')
        call output_array(p%image, './example_3d_image_unconf_1.bin')
        call output_array(p%fault, './example_3d_fault_unconf_1.bin')
        call output_array(p%fault_dip, './example_3d_fdip_unconf_1.bin')
        call output_array(p%fault_strike, './example_3d_fstrike_unconf_1.bin')
        call output_array(p%facies, './example_3d_fstrike_facies_1.bin')

    end block

    block

        type(rgm3_curved) :: p

        p%n1 = 171
        p%n2 = 251
        p%n3 = 251
        p%lwv = 0.4
        p%lwh = 0.2
        p%refl_shape = 'cauchy'
        p%ng = 3
        p%refl_sigma2 = [40.0, 60.0]
        p%refl_sigma3 = [30.0, 50.0]
        p%refl_shape_top = 'perlin'
        p%refl_smooth = 0
        p%refl_smooth_top = 4
        p%refl_slope = [20.0, -20.0]
        p%refl_slope_top = [-2.0, -5.0]
        p%refl_height = [0.0, 70.0]
        p%refl_height_top = [0.0, 40.0]
        p%nl = 25
        p%nf = 12
        p%disp = [-6.0, 6.0]
        p%yn_regular_fault = .true.
        p%strike = [45.0, 130.0]
        p%delta_dip = [0.0, 10.0]
        p%noise_level = 0.02
        p%psf_sigma = [10.0, 1.0, 1.0]
        p%seed = 4566
        p%unconf = 1
        p%unconf_height = [10.0, 20.0]
        p%unconf_z = [0.1, 0.3]
        p%unconf_smooth = 20
        p%yn_rgt = .true.
        call p%generate
        call output_array(p%vp, './example_3d_vp_unconf_2.bin')
        call output_array(p%rgt, './example_3d_rgt_unconf_2.bin')
        call output_array(p%image, './example_3d_image_unconf_2.bin')
        call output_array(p%fault, './example_3d_fault_unconf_2.bin')
        call output_array(p%fault_dip, './example_3d_fdip_unconf_2.bin')
        call output_array(p%fault_strike, './example_3d_fstrike_unconf_2.bin')

    end block

    block

        type(rgm3_curved) :: p

        p%n1 = 171
        p%n2 = 251
        p%n3 = 251
        p%lwv = 0.5
        p%lwh = 0.3
        p%refl_shape = 'random'
        p%refl_shape_top = 'perlin'
        p%refl_smooth = 40
        p%refl_smooth_top = 2
        p%refl_slope = [-20.0, 10.0]
        p%refl_slope_top = [5.0, 0.0]
        p%refl_height = [0.0, 50.0]
        p%refl_height_top = [0.0, 4.0]
        p%nl = 35
        p%nf = 12
        p%disp = [5.0, 10.0]
        p%delta_dip = [20.0, 40.0]
        p%noise_level = 0.05
        p%psf_sigma = [10.0, 1.0, 1.0]
        p%seed = 888999
        p%unconf = 2
        p%unconf_height = [5.0, 15.0]
        p%unconf_z = [0.1, 0.4]
        p%unconf_smooth = 30
        p%yn_rgt = .true.
        call p%generate
        call output_array(p%vp, './example_3d_vp_unconf_3.bin')
        call output_array(p%rgt, './example_3d_rgt_unconf_3.bin')
        call output_array(p%image, './example_3d_image_unconf_3.bin')
        call output_array(p%fault, './example_3d_fault_unconf_3.bin')
        call output_array(p%fault_dip, './example_3d_fdip_unconf_3.bin')
        call output_array(p%fault_strike, './example_3d_fstrike_unconf_3.bin')

    end block

    block

        type(rgm3_curved) :: p

        p%n1 = 171
        p%n2 = 251
        p%n3 = 251
        p%lwv = 0.4
        p%lwh = 0.4
        p%refl_shape = 'gaussian'
        p%refl_shape_top = 'same'
        p%refl_smooth = 0
        p%refl_smooth_top = 0
        p%ng = 3
        p%refl_sigma2 = [40.0, 60.0]
        p%refl_sigma3 = [40.0, 60.0]
        p%refl_slope = [10.0, -10.0]
        p%refl_slope_top = [20.0, -20.0]
        p%refl_height = [0.0, 50.0]
        p%refl_height_top = [0.0, 0.0]
        p%nl = 40
        p%nf = 7
        p%disp = [10.0, 20.0]
        p%noise_level = 0.02
        p%psf_sigma = [10.0, 1.0, 1.0]
        p%seed = 789789
        p%yn_rgt = .true.
        p%unconf = 2
        p%unconf_height = [5.0, 15.0]
        p%unconf_z = [0.1, 0.5]
        p%unconf_smooth = 0
        call p%generate
        call output_array(p%vp, './example_3d_vp_unconf_4.bin')
        call output_array(p%rgt, './example_3d_rgt_unconf_4.bin')
        call output_array(p%image, './example_3d_image_unconf_4.bin')
        call output_array(p%fault, './example_3d_fault_unconf_4.bin')
        call output_array(p%fault_dip, './example_3d_fdip_unconf_4.bin')
        call output_array(p%fault_strike, './example_3d_fstrike_unconf_4.bin')

    end block

    ! ==================================================================
    ! Batch labeled data

    block

        type(rgm3_curved) :: p

        integer :: i, n

        n = 16

        do i = 1, n

            p%n1 = 128
            p%n2 = 128
            p%n3 = 128
            p%lwv = 0.3
            p%lwh = 0.3
            p%f0 = rand(range=[180.0, 220.0], seed=i + 15*n)
            select case (mod(i, 4))
                case (0)
                    p%refl_shape = 'perlin'
                    p%refl_shape_top = 'perlin'
                    p%refl_smooth = 3
                    p%refl_smooth_top = 3
                case (1)
                    p%refl_shape = 'random'
                    p%refl_shape_top = 'random'
                    p%refl_smooth = 20
                    p%refl_smooth_top = 10
                case (2)
                    p%refl_shape = 'gaussian'
                    p%refl_shape_top = 'random'
                    p%refl_smooth = 0
                    p%refl_smooth_top = 10
                    p%ng = 3
                    p%refl_sigma2 = [20.0, 40.0]
                    p%refl_sigma3 = [20.0, 40.0]
                case (3)
                    p%refl_shape = 'cauchy'
                    p%refl_shape_top = 'random'
                    p%refl_smooth = 0
                    p%refl_smooth_top = 10
                    p%ng = 3
                    p%refl_sigma2 = [20.0, 40.0]
                    p%refl_sigma3 = [20.0, 40.0]
            end select
            p%refl_slope = [rand(range=[5.0, 20.0], seed=i), -rand(range=[5.0, 20.0], seed=i + n)]
            p%refl_slope_top = -[rand(range=[5.0, 10.0], seed=i + 2*n), -rand(range=[5.0, 10.0], seed=i + 3*n)]
            p%refl_height = [0.0, rand(range=[30.0, 50.0], seed=i + 4*n)]
            p%refl_height_top = [0.0, rand(range=[0.0, 10.0], seed=i + 5*n)]
            p%nl = irand(range=[20, 30], seed=i + 6*n)
            p%noise_level = rand(range=[0.05, 0.1], seed=i + 7*n)
            p%yn_conv_noise = .true.
            select case (mod(i, 2))
                case (0)
                    p%yn_regular_fault = .true.
                    p%nf = irand(range=[6, 12], seed=i + 8*n)
                    p%disp = [-2.0, 2.0]*rand(range=[3.0, 5.0], seed=i + 9*n)
                    p%strike = [0.0, rand(range=[80.0, 100.0], seed=i + 13*n)] + rand(range=[10.0, 80.0], seed=i + 14*n)
                case (1)
                    p%yn_regular_fault = .false.
                    p%nf = irand(range=[2, 7], seed=i + 8*n)
                    p%disp = [2.0, 5.0]*rand(range=[3.0, 4.0], seed=i + 9*n)
            end select
            p%delta_dip = [0.0, rand(range=[0.0, 20.0], seed=i + 10*n)]
            p%psf_sigma = [10.0, 1.0, 1.0]
            select case (mod(irand(range=[1, 10*n], seed=i + 16*n), 3))
                case (0, 1)
                    p%unconf = 0
                case (2)
                    p%unconf = 1
                    p%unconf_z = [0.1, 0.4]
            end select
            p%seed = i + 20*n

            call p%generate

            call output_array(p%image, './batch_image_'//num2str(i)//'.bin')
            call output_array(p%fault, './batch_fault_'//num2str(i)//'.bin')

            print *, i

        end do

    end block

    ! ==================================================================
    ! Strike-slip fault

    block

        type(rgm3_curved) :: p

        p%n1 = 171
        p%n2 = 251
        p%n3 = 251
        p%lwv = 0.5
        p%lwh = 0.4
        p%refl_shape = 'cauchy'
        p%ng = 4
        p%refl_sigma2 = [40.0, 60.0]
        p%refl_sigma3 = [30.0, 50.0]
        p%refl_shape_top = 'random'
        p%refl_smooth = 0
        p%refl_smooth_top = 10
        p%refl_slope = [20.0, -20.0]
        p%refl_slope_top = [-2.0, -5.0]
        p%refl_height = [0.0, 70.0]
        p%refl_height_top = [0.0, 5.0]
        p%nl = 35
        p%nf = 2
        p%strike = [20.0, 30.0]
        p%rake = [0.0, 0.0]
        p%disp = [20.0, 30.0]
        p%delta_dip = [0.0, 10.0]
        p%noise_level = 0.02
        p%psf_sigma = [10.0, 1.0, 1.0]
        p%seed = 123456
        p%yn_rgt = .true.
        call p%generate

        call output_array(p%vp, './example_3d_vp_strikeslip_1.bin')
        call output_array(p%rgt, './example_3d_rgt_strikeslip_1.bin')
        call output_array(p%image, './example_3d_image_strikeslip_1.bin')

    end block

end program test

