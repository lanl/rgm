
program test

    use libflit
    use librgm
    use omp_lib

    implicit none

    ! ==================================================================
    ! Types of unconformities

    block

        type(rgm2_curved) :: p

        p%n1 = 128
        p%n2 = 128
        p%lwv = 0.3
        p%lwh = 0.1
        p%refl_shape = 'gaussian'
        p%ng = 1
        p%refl_mu2 = [70, 70]
        p%refl_sigma2 = [40, 40]
        p%refl_shape_top = 'perlin'
        p%refl_smooth_top = 3
        p%nl = 35
        p%seed = 9999

        ! Disconformity
        p%refl_height = [0.0, 1.0]
        p%refl_height_top = [0.0, 1.0]
        p%nf = 0
        p%unconf_refl_height = [0.0, 1.0]
        p%unconf = 1
        p%unconf_z = [0.32, 0.32]
        p%unconf_height = [15, 15]
        p%unconf_smooth = 0
        call p%generate
        call output_array(p%vp, './example_2d_unconf_type_1.bin')

        ! Paraconformity
        p%refl_height = [0.0, 1.0]
        p%refl_height_top = [0.0, 1.0]
        p%nf = 3
        p%disp = [3, 6]
        p%unconf_refl_height = [0.0, 1.0]
        p%unconf = 1
        p%unconf_z = [0.32, 0.32]
        p%unconf_height = [0, 0]
        p%unconf_smooth = 15
        call p%generate
        call output_array(p%vp, './example_2d_unconf_type_2.bin')

        ! Angular unconformity
        p%nf = 4
        p%lwv = 0.4
        p%lwh = 0.3
        p%refl_height = [-60, 0]
        p%refl_height_top = [-10, 0]
        p%disp = [5.0, 10.0]
        p%unconf = 1
        p%unconf_z = [0.3, 0.3]
        p%unconf_height = [5, 5]
        p%unconf_smooth = 15
        call p%generate
        call output_array(p%vp, './example_2d_unconf_type_3.bin')

        ! Angular unconformity
        p%nf = 4
        p%lwv = 0.4
        p%lwh = 0.3
        p%refl_height = [0, 60]
        p%refl_height_top = [0, 10]
        p%disp = [5.0, 10.0]
        p%unconf = 1
        p%unconf_z = [0.3, 0.3]
        p%unconf_height = [20, 20]
        p%unconf_smooth = 0
        call p%generate
        call output_array(p%vp, './example_2d_unconf_type_4.bin')

        ! Nonconformity
        p%nf = 4
        p%lwv = 0.4
        p%lwh = 0.3
        p%refl_height = [0, 60]
        p%refl_height_top = [0, 10]
        p%disp = [5.0, 10.0]
        p%unconf = 1
        p%unconf_z = [0.3, 0.3]
        p%unconf_height = [5, 5]
        p%unconf_smooth = 15
        p%yn_salt = .true.
        p%nsalt = 1
        p%salt_top_z = [0.2, 0.2]
        p%salt_nnode = 5
        p%salt_radius = [36.0, 36.0]
        call p%generate
        call output_array(p%vp, './example_2d_unconf_type_5.bin')

        ! Nonconformity
        p%lwh = 0.1
        p%nf = 0
        p%lwv = 0.4
        p%lwh = 0.3
        p%refl_height = [0, 1]
        p%refl_height_top = [0, 1]
        p%unconf = 0
        p%yn_salt = .true.
        p%nsalt = 1
        p%salt_top_height = 15
        p%salt_top_z = [0.4, 0.4]
        p%salt_radius = [120, 120]
        call p%generate
        call output_array(p%vp, './example_2d_unconf_type_6.bin')

    end block

    ! ==================================================================
    ! Scalability test

    block

        type(rgm2_curved) :: p
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
                p%lwv = 0.3
                p%lwh = 0.3
                p%refl_shape = 'cauchy'
                p%refl_shape_top = 'perlin'
                p%refl_smooth_top = 3
                p%refl_slope = rand(range=[-20.0, 20.0], seed=i)
                p%refl_slope_top = rand(range=[-5.0, 5.0], seed=i + 100)
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

        call output_array(ctime, './ctime2.bin', transp=.true.)
        call output_array(ctime, './ctime2.txt', transp=.true., ascii=.true.)

    end block

    ! ==================================================================
    ! Types of random curves/surfaces

    block

        integer :: n2, n3
        real, allocatable, dimension(:) :: r, rr
        real, allocatable, dimension(:, :) :: t, tt
        type(fractal_noise_1d) :: pn
        type(fractal_noise_2d) :: qn

        n2 = 200
        n3 = 300

        ! Random
        r = random(n2, seed=123)
        rr = gauss_filt(r, 15.0)
        rr = rescale(rr, [0.0, 1.0])
        call output_array(rr, './random_1d.bin')

        t = random(n2, n3, seed=246)
        tt = gauss_filt(t, [15.0, 15.0])
        tt = rescale(tt, [0.0, 1.0])
        call output_array(tt, './random_2d.bin')

        ! Gaussian
        r = rescale(gaussian(linspace(0.0, n2 - 1.0, n2), 50.0, 20.0), [0.0, 10.0]) &
            + rescale(gaussian(linspace(0.0, n2 - 1.0, n2), 150.0, 40.0), [0.0, 20.0])
        r = rescale(r, [0.0, 1.0])
        call output_array(r, 'gaussian_1d.bin')

        t = rescale(gaussian(linspace(0.0, n2 - 1.0, n2), linspace(0.0, n3 - 1.0, n3), &
            [50.0, 60.0], [40.0, 30.0], 30*real(const_deg2rad)), [0.0, 10.0]) &
            + rescale(gaussian(linspace(0.0, n2 - 1.0, n2), linspace(0.0, n3 - 1.0, n3), &
            [150.0, 120.0], [40.0, 60.0] , 45*real(const_deg2rad)), [0.0, 20.0]) &
            + rescale(gaussian(linspace(0.0, n2 - 1.0, n2), linspace(0.0, n3 - 1.0, n3), &
            [50.0, 250.0], [30.0, 40.0], 120*real(const_deg2rad)), [0.0, 15.0])
        t = rescale(t, [0.0, 1.0])
        call output_array(t, 'gaussian_2d.bin')

        ! Cauchy
        r = rescale(cauchy(linspace(0.0, n2 - 1.0, n2), 50.0, 20.0), [0.0, 10.0]) &
            + rescale(cauchy(linspace(0.0, n2 - 1.0, n2), 150.0, 40.0), [0.0, 20.0])
        r = rescale(r, [0.0, 1.0])
        call output_array(r, 'cauchy_1d.bin')

        t = rescale(cauchy(linspace(0.0, n2 - 1.0, n2), linspace(0.0, n3 - 1.0, n3), &
            [50.0, 60.0], [40.0, 30.0], 30*real(const_deg2rad)), [0.0, 10.0]) &
            + rescale(cauchy(linspace(0.0, n2 - 1.0, n2), linspace(0.0, n3 - 1.0, n3), &
            [150.0, 120.0], [40.0, 60.0] , 45*real(const_deg2rad)), [0.0, 20.0]) &
            + rescale(cauchy(linspace(0.0, n2 - 1.0, n2), linspace(0.0, n3 - 1.0, n3), &
            [50.0, 250.0], [30.0, 40.0], 120*real(const_deg2rad)), [0.0, 15.0])
        t = rescale(t, [0.0, 1.0])
        call output_array(t, 'cauchy_2d.bin')

        ! Perlin
        pn%n1 = n2
        pn%octaves = 4
        pn%seed = 222
        r = pn%generate()
        r = gauss_filt(r, 2.0)
        r = rescale(r, [0.0, 1.0])
        call output_array(r, 'perlin_1d.bin')

        qn%n1 = n2
        qn%n2 = n3
        qn%octaves = 4
        qn%seed = 111
        t = qn%generate()
        t = gauss_filt(t, [2.0, 2.0])
        t = rescale(t, [0.0, 1.0])
        call output_array(t, 'perlin_2d.bin')

    end block

    ! ==================================================================
    ! Unfaulted models

    block

        type(rgm2_curved) :: p

        p%n1 = 201
        p%n2 = 301
        p%lwv = 0.4
        p%lwh = 0.2
        p%refl_shape = 'perlin'
        p%refl_shape_top = 'perlin'
        p%refl_smooth = 4
        p%refl_smooth_top = 4
        p%refl_slope = 20.0
        p%refl_slope_top = -2.0
        p%refl_height = [0.0, 50.0]
        p%refl_height_top = [0.0, 4.0]
        p%nl = 30
        p%nf = 0
        p%noise_level = 0.05
        p%psf_sigma = [10.0, 1.0]
        p%seed = 123
        p%yn_rgt = .true.
        call p%generate
        call output_array(p%vp, './example_2d_vp_unfaulted_1.bin')
        call output_array(p%rgt, './example_2d_rgt_unfaulted_1.bin')
        call output_array(p%image, './example_2d_image_unfaulted_1.bin')

    end block

    block

        type(rgm2_curved) :: p

        p%n1 = 201
        p%n2 = 301
        p%lwv = 0.4
        p%lwh = 0.2
        p%refl_shape = 'cauchy'
        p%refl_shape_top = 'perlin'
        p%refl_smooth = 0
        p%refl_smooth_top = 4
        p%refl_slope = 20.0
        p%refl_slope_top = -2.0
        p%refl_height = [0.0, 70.0]
        p%refl_height_top = [0.0, 5.0]
        p%nl = 25
        p%nf = 0
        p%noise_level = 0.1
        p%psf_sigma = [10.0, 1.0]
        p%seed = 456
        p%yn_rgt = .true.
        call p%generate
        call output_array(p%vp, './example_2d_vp_unfaulted_2.bin')
        call output_array(p%rgt, './example_2d_rgt_unfaulted_2.bin')
        call output_array(p%image, './example_2d_image_unfaulted_2.bin')

    end block

    block

        type(rgm2_curved) :: p

        p%n1 = 201
        p%n2 = 301
        p%lwv = 0.2
        p%lwh = 0.2
        p%refl_shape = 'random'
        p%refl_shape_top = 'perlin'
        p%refl_smooth = 20
        p%refl_smooth_top = 2
        p%refl_slope = -20.0
        p%refl_slope_top = 5.0
        p%refl_height = [0.0, 50.0]
        p%refl_height_top = [0.0, 4.0]
        p%nl = 15
        p%nf = 0
        p%noise_level = 0.01
        p%psf_sigma = [10.0, 1.0]
        p%seed = 888
        p%yn_rgt = .true.
        call p%generate
        call output_array(p%vp, './example_2d_vp_unfaulted_3.bin')
        call output_array(p%rgt, './example_2d_rgt_unfaulted_3.bin')
        call output_array(p%image, './example_2d_image_unfaulted_3.bin')

    end block

    block

        type(rgm2_curved) :: p

        p%n1 = 201
        p%n2 = 301
        p%lwv = 0.4
        p%lwh = 0.4
        p%refl_shape = 'gaussian'
        p%refl_shape_top = 'same'
        p%refl_smooth = 0
        p%refl_smooth_top = 0
        p%ng = 3
        p%refl_sigma2 = [50.0, 100.0]
        p%refl_slope = 10.0
        p%refl_slope_top = 20.0
        p%refl_height = [0.0, 50.0]
        p%refl_height_top = [0.0, 0.0]
        p%nl = 40
        p%nf = 0
        p%noise_level = 0.02
        p%psf_sigma = [10.0, 1.0]
        p%seed = 789
        p%yn_rgt = .true.
        call p%generate
        call output_array(p%vp, './example_2d_vp_unfaulted_4.bin')
        call output_array(p%rgt, './example_2d_rgt_unfaulted_4.bin')
        call output_array(p%image, './example_2d_image_unfaulted_4.bin')

    end block

    ! ==================================================================
    ! Fauled models

    block

        type(rgm2_curved) :: p

        p%n1 = 201
        p%n2 = 301
        p%lwv = 0.4
        p%lwh = 0.2
        p%refl_shape = 'perlin'
        p%refl_shape_top = 'perlin'
        p%refl_smooth = 4
        p%refl_smooth_top = 4
        p%refl_slope = 20.0
        p%refl_slope_top = -2.0
        p%refl_height = [0.0, 50.0]
        p%refl_height_top = [0.0, 4.0]
        p%nl = 30
        p%nf = 3
        p%disp = [3.0, 10.0]
        p%delta_dip = [0.0, 10.0]
        p%noise_level = 0.05
        p%psf_sigma = [10.0, 1.0]
        p%seed = 123
        p%yn_rgt = .true.
        call p%generate
        call output_array(p%vp, './example_2d_vp_faulted_1.bin')
        call output_array(p%rgt, './example_2d_rgt_faulted_1.bin')
        call output_array(p%image, './example_2d_image_faulted_1.bin')
        call output_array(p%fault, './example_2d_fault_1.bin')
        call output_array(p%fault_dip, './example_2d_fdip_1.bin')

    end block

    block

        type(rgm2_curved) :: p

        p%n1 = 201
        p%n2 = 301
        p%lwv = 0.4
        p%lwh = 0.2
        p%refl_shape = 'cauchy'
        p%refl_shape_top = 'perlin'
        p%refl_smooth = 0
        p%refl_smooth_top = 4
        p%refl_slope = 20.0
        p%refl_slope_top = -2.0
        p%refl_height = [0.0, 70.0]
        p%refl_height_top = [0.0, 5.0]
        p%nl = 25
        p%nf = 12
        p%yn_regular_fault = .true.
        p%disp = [-5.0, 5.0]
        p%delta_dip = [0.0, 10.0]
        p%noise_level = 0.1
        p%psf_sigma = [10.0, 1.0]
        p%seed = 456
        p%yn_rgt = .true.
        call p%generate
        call output_array(p%vp, './example_2d_vp_faulted_2.bin')
        call output_array(p%rgt, './example_2d_rgt_faulted_2.bin')
        call output_array(p%image, './example_2d_image_faulted_2.bin')
        call output_array(p%fault, './example_2d_fault_2.bin')
        call output_array(p%fault_dip, './example_2d_fdip_2.bin')

    end block

    block

        type(rgm2_curved) :: p

        p%n1 = 201
        p%n2 = 301
        p%lwv = 0.2
        p%lwh = 0.2
        p%refl_shape = 'random'
        p%refl_shape_top = 'perlin'
        p%refl_smooth = 20
        p%refl_smooth_top = 2
        p%refl_slope = -20.0
        p%refl_slope_top = 5.0
        p%refl_height = [0.0, 50.0]
        p%refl_height_top = [0.0, 4.0]
        p%nl = 15
        p%nf = 7
        p%yn_regular_fault = .false.
        p%disp = [5.0, 15.0]
        p%delta_dip = [20.0, 40.0]
        p%noise_level = 0.01
        p%psf_sigma = [10.0, 1.0]
        p%seed = 888
        p%yn_rgt = .true.
        call p%generate
        call output_array(p%vp, './example_2d_vp_faulted_3.bin')
        call output_array(p%rgt, './example_2d_rgt_faulted_3.bin')
        call output_array(p%image, './example_2d_image_faulted_3.bin')
        call output_array(p%fault, './example_2d_fault_3.bin')
        call output_array(p%fault_dip, './example_2d_fdip_3.bin')

    end block

    block

        type(rgm2_curved) :: p

        p%n1 = 201
        p%n2 = 301
        p%lwv = 0.4
        p%lwh = 0.4
        p%refl_shape = 'gaussian'
        p%refl_shape_top = 'same'
        p%refl_smooth = 0
        p%refl_smooth_top = 0
        p%ng = 3
        p%refl_sigma2 = [50.0, 100.0]
        p%refl_slope = 10.0
        p%refl_slope_top = 20.0
        p%refl_height = [0.0, 50.0]
        p%refl_height_top = [0.0, 0.0]
        p%nl = 40
        p%nf = 5
        p%yn_regular_fault = .false.
        p%disp = [5.0, 15.0]
        p%delta_dip = [0.0, 40.0]
        p%noise_level = 0.02
        p%psf_sigma = [10.0, 1.0]
        p%seed = 789
        p%yn_rgt = .true.
        call p%generate
        call output_array(p%vp, './example_2d_vp_faulted_4.bin')
        call output_array(p%rgt, './example_2d_rgt_faulted_4.bin')
        call output_array(p%image, './example_2d_image_faulted_4.bin')
        call output_array(p%fault, './example_2d_fault_4.bin')
        call output_array(p%fault_dip, './example_2d_fdip_4.bin')

    end block

    ! ==================================================================
    ! Salt models

    block

        type(rgm2_curved) :: p

        p%n1 = 201
        p%n2 = 301
        p%lwv = 0.4
        p%lwh = 0.2
        p%refl_shape = 'perlin'
        p%refl_shape_top = 'perlin'
        p%refl_smooth = 5
        p%refl_smooth_top = 15
        p%refl_slope = -20.0
        p%refl_slope_top = 4.0
        p%refl_height = [0.0, 65.0]
        p%refl_height_top = [0.0, 4.0]
        p%nl = 40
        p%nf = 6
        p%disp = [5.0, 10.0]
        p%noise_level = 0
        p%psf_sigma = [10.0, 1.0]
        p%seed = 1123
        p%yn_rgt = .true.
        p%yn_salt = .true.
        p%nsalt = 2
        call p%generate
        call output_array(p%vp, './example_2d_vp_salt_1.bin')
        call output_array(p%rgt, './example_2d_rgt_salt_1.bin')
        call output_array(p%image, './example_2d_image_salt_1.bin')

    end block

    block

        type(rgm2_curved) :: p

        p%n1 = 201
        p%n2 = 301
        p%lwv = 0.4
        p%lwh = 0.2
        p%refl_shape = 'perlin'
        p%refl_shape_top = 'perlin'
        p%refl_smooth = 5
        p%refl_smooth_top = 15
        p%refl_slope = -20.0
        p%refl_slope_top = 4.0
        p%refl_height = [0.0, 65.0]
        p%refl_height_top = [0.0, 4.0]
        p%nl = 40
        p%nf = 11
        p%yn_regular_fault = .true.
        p%disp = [-2.0, 4.0]
        p%noise_level = 0
        p%psf_sigma = [10.0, 1.0]
        p%seed = 778899
        p%yn_rgt = .true.
        p%yn_salt = .true.
        p%nsalt = 3
        p%salt_radius = [25.0, 45.0]
        p%salt_top_z = [0.5, 0.8]
        call p%generate
        call output_array(p%vp, './example_2d_vp_salt_2.bin')
        call output_array(p%rgt, './example_2d_rgt_salt_2.bin')
        call output_array(p%image, './example_2d_image_salt_2.bin')

    end block

    block

        type(rgm2_curved) :: p

        p%n1 = 201
        p%n2 = 301
        p%lwv = 0.2
        p%lwh = 0.3
        p%refl_shape = 'gaussian'
        p%refl_shape_top = 'perlin'
        p%refl_smooth = 0
        p%refl_smooth_top = 15
        p%refl_slope = -20.0
        p%refl_slope_top = 4.0
        p%refl_height = [0.0, 65.0]
        p%refl_height_top = [0.0, 4.0]
        p%nl = 40
        p%nf = 11
        p%yn_regular_fault = .true.
        p%disp = [-2.0, 4.0]
        p%noise_level = 0.02
        p%yn_conv_noise = .false.
        p%psf_sigma = [10.0, 1.0]
        p%seed = 135
        p%yn_rgt = .true.
        p%yn_salt = .true.
        p%nsalt = 3
        p%salt_radius = [15.0, 45.0]
        p%salt_top_z = [0.4, 0.7]
        call p%generate
        call output_array(p%vp, './example_2d_vp_salt_3.bin')
        call output_array(p%rgt, './example_2d_rgt_salt_3.bin')
        call output_array(p%image, './example_2d_image_salt_3.bin')

    end block

    block

        type(rgm2_curved) :: p

        p%n1 = 201
        p%n2 = 301
        p%lwv = 0.2
        p%lwh = 0.3
        p%refl_shape = 'random'
        p%refl_shape_top = 'random'
        p%refl_smooth = 15
        p%refl_smooth_top = 20
        p%refl_slope = -20.0
        p%refl_slope_top = 7.0
        p%refl_height = [0.0, 65.0]
        p%refl_height_top = [0.0, 4.0]
        p%nl = 35
        p%nf = 13
        p%disp = [3.0, 10.0]
        p%noise_level = 0.02
        p%yn_conv_noise = .true.
        p%psf_sigma = [10.0, 1.0]
        p%seed = 1357
        p%yn_rgt = .true.
        p%yn_salt = .true.
        p%nsalt = 4
        p%salt_top_height = 30
        p%salt_radius = [15.0, 35.0]
        p%salt_radius_variation = 0.9
        p%salt_top_z = [0.4, 0.9]
        call p%generate
        call output_array(p%vp, './example_2d_vp_salt_4.bin')
        call output_array(p%rgt, './example_2d_rgt_salt_4.bin')
        call output_array(p%image, './example_2d_image_salt_4.bin')

    end block

    ! ==================================================================
    ! Unconformity models

    block

        type(rgm2_curved) :: p

        p%n1 = 201
        p%n2 = 301
        p%lwv = 0.4
        p%lwh = 0.2
        p%refl_shape = 'perlin'
        p%refl_shape_top = 'perlin'
        p%refl_smooth = 5
        p%refl_smooth_top = 15
        p%refl_slope = -20.0
        p%refl_slope_top = 4.0
        p%refl_height = [0.0, 60.0]
        p%refl_height_top = [0.0, 2.0]
        p%nl = 40
        p%nf = 6
        p%disp = [5.0, 10.0]
        p%noise_level = 0
        p%psf_sigma = [10.0, 1.0]
        p%seed = 7774
        p%yn_rgt = .true.
        p%unconf = 1
        p%unconf_z = [0.1, 0.3]
        p%unconf_height = [20.0, 30.0]
        call p%generate
        call output_array(p%vp, './example_2d_vp_unconf_1.bin')
        call output_array(p%rgt, './example_2d_rgt_unconf_1.bin')
        call output_array(p%image, './example_2d_image_unconf_1.bin')
        call output_array(p%fault, './example_2d_fault_unconf_1.bin')
        call output_array(p%fault_dip, './example_2d_fdip_unconf_1.bin')

    end block

    block

        type(rgm2_curved) :: p

        p%n1 = 201
        p%n2 = 301
        p%lwv = 0.4
        p%lwh = 0.2
        p%refl_shape = 'random'
        p%refl_shape_top = 'random'
        p%refl_smooth = 10
        p%refl_smooth_top = 20
        p%refl_slope = -20.0
        p%refl_slope_top = 4.0
        p%refl_height = [0.0, 14.0]
        p%refl_height_top = [0.0, 10.0]
        p%nl = 40
        p%nf = 6
        p%disp = [5.0, 10.0]
        p%delta_dip = [0.0, 0.0]
        p%noise_level = 0
        p%psf_sigma = [10.0, 1.0]
        p%seed = 8881
        p%yn_rgt = .true.
        p%unconf = 2
        p%unconf_z = [0.05, 0.4]
        p%unconf_smooth = 3
        p%unconf_height = [5.0, 10.0]
        call p%generate
        call output_array(p%vp, './example_2d_vp_unconf_2.bin')
        call output_array(p%rgt, './example_2d_rgt_unconf_2.bin')
        call output_array(p%image, './example_2d_image_unconf_2.bin')
        call output_array(p%fault, './example_2d_fault_unconf_2.bin')
        call output_array(p%fault_dip, './example_2d_fdip_unconf_2.bin')

    end block

    block

        type(rgm2_curved) :: p

        p%n1 = 201
        p%n2 = 301
        p%lwv = 0.3
        p%lwh = 0.3
        p%refl_shape = 'cauchy'
        p%refl_shape_top = 'same'
        p%refl_smooth = 0
        p%refl_smooth_top = 0
        p%refl_slope = -20.0
        p%refl_slope_top = 4.0
        p%refl_height = [0.0, 50.0]
        p%refl_height_top = [0.0, 50.0]
        p%nl = 45
        p%nf = 13
        p%yn_regular_fault = .true.
        p%disp = [-5.0, 5.0]
        p%noise_level = 0.025
        p%yn_conv_noise = .true.
        p%psf_sigma = [10.0, 1.0]
        p%seed = 9991
        p%f0 = 220
        p%delta_dip = [0.0, 0.0]
        p%noise_smooth = [0, 0]
        p%yn_rgt = .true.
        p%yn_salt = .true.
        p%nsalt = 3
        p%salt_top_height = 30
        p%salt_radius = [15.0, 35.0]
        p%salt_radius_variation = 0.9
        p%salt_top_z = [0.1, 0.7]
        p%unconf = 1
        p%unconf_z = [0.3, 0.5]
        p%unconf_smooth = 2
        p%unconf_height = [5.0, 15.0]
        p%yn_facies = .true.
        call p%generate
        call output_array(p%vp, './example_2d_vp_unconf_3.bin')
        call output_array(p%rgt, './example_2d_rgt_unconf_3.bin')
        call output_array(p%image, './example_2d_image_unconf_3.bin')
        call output_array(p%fault, './example_2d_fault_unconf_3.bin')
        call output_array(p%fault_dip, './example_2d_fdip_unconf_3.bin')
        call output_array(p%facies, './example_2d_facies_unconf_3.bin')

    end block

    block

        type(rgm2_curved) :: p

        p%n1 = 201
        p%n2 = 301
        p%lwv = 0.5
        p%lwh = 0.2
        p%refl_shape = 'cauchy'
        p%refl_shape_top = 'same'
        p%refl_smooth = 0
        p%refl_smooth_top = 0
        p%refl_sigma2 = [50.0, 80.0]
        p%refl_slope = -10.0
        p%refl_slope_top = 10.0
        p%refl_height = [0.0, 50.0]
        p%refl_height_top = [0.0, 50.0]
        p%nl = 35
        p%nf = 17
        p%yn_regular_fault = .true.
        p%yn_conv_noise = .true.
        p%disp = [-5.0, 5.0]
        p%noise_level = 0.05
        p%psf_sigma = [10.0, 1.0]
        p%seed = 1357
        p%yn_rgt = .true.
        p%unconf = 2
        p%unconf_z = [0.05, 0.4]
        p%unconf_refl_shape = 'random'
        p%unconf_refl_height = [0.0, 3.0]
        p%unconf_smooth = 20
        p%unconf_height = [10.0, 20.0]
        call p%generate
        call output_array(p%vp, './example_2d_vp_unconf_4.bin')
        call output_array(p%rgt, './example_2d_rgt_unconf_4.bin')
        call output_array(p%image, './example_2d_image_unconf_4.bin')
        call output_array(p%fault, './example_2d_fault_unconf_4.bin')
        call output_array(p%fault_dip, './example_2d_fdip_unconf_4.bin')

    end block

    ! ==================================================================
    ! Elastic models

    block

        type(rgm2_curved) :: p

        p%n1 = 201
        p%n2 = 301
        p%lwv = 0.3
        p%lwh = 0.4
        p%refl_shape = 'cauchy'
        p%refl_shape_top = 'perlin'
        p%refl_smooth = 0
        p%refl_smooth_top = 0
        p%refl_sigma2 = [50.0, 80.0]
        p%refl_slope = -10.0
        p%refl_slope_top = 10.0
        p%refl_height = [0.0, 50.0]
        p%refl_height_top = [0.0, 10.0]
        p%nl = 35
        p%nf = 11
        p%yn_conv_noise = .true.
        p%disp = [10.0, 20.0]
        p%noise_level = 0.025
        p%psf_sigma = [10.0, 1.0]
        p%seed = 135799
        p%unconf = 1
        p%unconf_z = [0.05, 0.4]
        p%unconf_refl_height = [0.0, 3.0]
        p%unconf_height = [10.0, 20.0]
        p%yn_elastic = .true.
        call p%generate
        call output_array(p%vp, './example_2d_vp_elastic_1.bin')
        call output_array(p%vs, './example_2d_vs_elastic_1.bin')
        call output_array(p%rho, './example_2d_rho_elastic_1.bin')
        call output_array(p%fault_dip, './example_2d_fdip_elastic_1.bin')
        call output_array(p%image_pp, './example_2d_image_pp_elastic_1.bin')
        call output_array(p%image_ps, './example_2d_image_ps_elastic_1.bin')
        call output_array(p%image_sp, './example_2d_image_sp_elastic_1.bin')
        call output_array(p%image_ss, './example_2d_image_ss_elastic_1.bin')

        call output_array(fftshift(abs(fft(p%image_pp))), './spec_pp.bin')
        call output_array(fftshift(abs(fft(p%image_ps))), './spec_ps.bin')
        call output_array(fftshift(abs(fft(p%image_sp))), './spec_sp.bin')
        call output_array(fftshift(abs(fft(p%image_ss))), './spec_ss.bin')

    end block

end program test

