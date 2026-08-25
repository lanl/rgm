
program test

    use libflit
    use librgm

    implicit none

    block

        type(rgm3_curved) :: p
        integer :: iter

        do iter = 1, 10

            p%n1 = 121
            p%n2 = 181
            p%n3 = 181
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
            p%nf = 5
            p%dip = [60.0, 120.0]
            p%disp = [15.0, 25.0]
            p%delta_dip = [0.0, 15.0]
            p%delta_strike = [15.0, 25.0]
            p%yn_vary_disp = .true.
            p%psf_sigma = [10.0, 1.0, 1.0]

            p%seed = 1
            call p%generate
            call plot_histogram(p%vp)
            call output_array(p%vp, './random_repeatability_'//num2str(iter)//'_seed='//num2str(p%seed)//'.bin')
            print *, iter, p%seed
            print *

            p%seed = 123
            call p%generate
            call plot_histogram(p%vp)
            call output_array(p%vp, './random_repeatability_'//num2str(iter)//'_seed='//num2str(p%seed)//'.bin')
            print *, iter, p%seed
            print *

            p%seed = 12345
            call p%generate
            call plot_histogram(p%vp)
            call output_array(p%vp, './random_repeatability_'//num2str(iter)//'_seed='//num2str(p%seed)//'.bin')
            print *, iter, p%seed
            print *

            p%seed = 123456
            call p%generate
            call plot_histogram(p%vp)
            call output_array(p%vp, './random_repeatability_'//num2str(iter)//'_seed='//num2str(p%seed)//'.bin')
            print *, iter, p%seed
            print *

            p%seed = 1234567
            call p%generate
            call plot_histogram(p%vp)
            call output_array(p%vp, './random_repeatability_'//num2str(iter)//'_seed='//num2str(p%seed)//'.bin')
            print *, iter, p%seed
            print *

        end do

    end block

end program

