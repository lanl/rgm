!
! Generate the meandering-channel erosional depth maps for the three
! channel-length values of the calibrated length-mapping figure, using
! the same mapping formulas as the generator and a fixed seed.
!
program gen_meander_length

    use libflit
    use librgm

    implicit none

    type(meandering_channel) :: mc
    real, dimension(1:3) :: lens
    real :: len_eff, wf
    integer :: i

    lens = [12000.0, 25000.0, 100000.0]
    wf = 0.05

    do i = 1, 3

        len_eff = max(lens(i), 10000.0)
        mc%nx = 251
        mc%ny = 201
        mc%nz = 64
        mc%length = len_eff
        mc%W = 0.025*min(len_eff, 25000.0)
        mc%W_render = wf*len_eff
        mc%W_shape = 3
        mc%n_bends = max(5, nint(30.0*len_eff/25000.0))
        mc%n_iter = 1500
        mc%terrain_bg = 0.25*(mc%nz - 1.0)
        mc%seed = 12345
        mc%orient = 0
        mc%yn_depth_only = .true.
        call mc%generate
        call output_array(mc%depth_map, './mlen_channel_'//num2str(nint(lens(i)))//'.bin')

        print *, 'length = ', lens(i), ' done'

    end do

end program gen_meander_length
