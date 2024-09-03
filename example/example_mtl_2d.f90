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


program main

    use libflit
    use librgm

    implicit none

    type(rgm2) :: p, q
    integer :: i, ibeg, iend, nt, nv, nm
    real, allocatable, dimension(:) :: psfm, zr1, zr2, sr, rx, fs, noise, height, slope
    real, allocatable, dimension(:) :: lwv, nsmooth1, nsmooth2, height2, lwv2
    real, allocatable, dimension(:) :: nf, nff, nl
    real, allocatable, dimension(:, :) :: r, iter_rgt, iter_dhr, iter_fsem, iter_fdip
    character(len=1024) :: dir_output
    integer, allocatable, dimension(:) :: sd

    dir_output = './dataset2_mtl'
    nt = 6000
    nv = 600
    nm = nt + nv

    call make_directory(tidy(dir_output)//'/data_train')
    call make_directory(tidy(dir_output)//'/data_valid')
    call make_directory(tidy(dir_output)//'/target_train')
    call make_directory(tidy(dir_output)//'/target_valid')

    fs = random(nm, range=[100.0, 180.0], seed=101)
    nf = irandom(nm, range=[1, 12], seed=1111)
    nff = irandom(nm, range=[1, 12], seed=2222)
    nl = nint(rescale(fs, [20.0, 60.0]))
    noise = random(nm, range=[0.0, 0.6], seed=3333)
    height = random(nm, range=[2.0, 12.0], seed=4444)
    height2 = random(nm, range=[10.0, 20.0], seed=5555)
    lwv = random(nm, range=[0.0, 0.1], seed=6666)
    lwv2 = random(nm, range=[0.0, 0.1], seed=7777)
    slope = random(nm, range=[-50.0, 50.0], seed=8888)
    nsmooth1 = random(nm, range=[0.0, 3.5], seed=1234)
    nsmooth2 = random(nm, range=[0.0, 3.5], seed=4321)
    where (nsmooth1 < 0.5)
        nsmooth1 = 0.0
    end where
    where (nsmooth2 < 0.5)
        nsmooth2 = 0.0
    end where
    rx = random(nm, range=[0.0, 0.1], seed=9999)
    sr = random(nm, range=[0.0, 0.1], seed=111)
    sd = irandom(nm, range=[1, 10*nm], seed=222)
    zr1 = random(nm, range=[0.05, 0.45], seed=333)
    zr2 = random(nm, range=[0.85, 0.95], seed=444)
    psfm = random(nm, range=[1.0, 3.0], seed=555)

    call getpar_int('ibeg', ibeg, 1)
    call getpar_int('iend', iend, nm)

    do i = ibeg, iend

        ! ===============================================================
        ! Generate labels for for inference and refinement NN

        p%seed = sd(i)

        p%n1 = 360
        p%n2 = 256
        p%nf = nf(i)
        p%refl_slope = slope(i)
        p%nl = nl(i)
        p%refl_amp = [0.1, 1.0]
        p%secondary_refl_amp = sr(i)
        p%fwidth = 3.0

        if (mod(i, 4) == 0) then
            p%refl_shape = 'gaussian'
            p%refl_mu2 = [0.0, p%n2 - 1.0]
            p%refl_sigma2 = [20.0, 50.0]
            p%ng = irand(range=[1, 3], seed=sd(i)*8)
            p%refl_height = [0.0, height2(i)]
            p%lwv = lwv2(i)
        else
            p%refl_shape = 'random'
            p%refl_smooth = 25
            p%refl_height = [0.0, height(i)]
            p%lwv = lwv(i)
        end if

        select case (mod(i, 3))
            case (0)
                p%wave = 'ricker'
            case (1)
                p%wave = 'ricker_deriv'
            case (2)
                p%wave = 'gaussian_deriv'
        end select

        if (mod(i, 3) == 0) then
            p%refl_amp_dist = 'normal'
        else
            p%refl_amp_dist = 'uniform'
        end if

        if (mod(i + 1, 10) == 0) then
            p%yn_regular_fault = .true.
            p%nf = irand(range=[10, 25], seed=sd(i)*2)
            if (mod(irand(range=[1, 10], seed=sd(i)*3), 2) == 0) then
                p%dip = [rand(range=[100.0, 120.0], seed=sd(i) + 2), rand(range=[60.0, 80.0], seed=sd(i) + 3)]
                p%disp = [5.0, -5.0]
            else
                p%dip = [rand(range=[60.0, 80.0], seed=sd(i) + 4), rand(range=[100.0, 120.0], seed=sd(i) + 5)]
                p%disp = [-5.0, 5.0]
            end if
        else
            p%yn_regular_fault = .false.
            p%nf = nf(i)
            p%disp = [5.0, 30.0]
            p%dip = [50.0, 130.0]
        end if

        if (mod(i + 3, 5) == 0) then
            p%unconf = irand(range=[1, 2], seed=sd(i) + 11)
            p%unconf_amp = [0.05, 0.1]
            p%unconf_z = [0.05, 0.35]
        else
            p%unconf = 0
        end if

        p%f0 = fs(i)*1.25
        p%noise_level = 0
        p%psf_sigma = [12.0, 0.0]
        p%yn_rgt = .true.
        p%yn_fault = .true.
        call p%generate

        p%image = p%image/norm2(p%image)*100

        where (p%fault /= 0)
            p%fault = 1.0
        end where

        if (i <= nt) then
            call output_array(p%image, tidy(dir_output)//'/target_train/'//num2str(i - 1)//'_dhr.bin')
            call output_array(p%rgt, tidy(dir_output)//'/target_train/'//num2str(i - 1)//'_rgt.bin')
            call output_array(p%fault, tidy(dir_output)//'/target_train/'//num2str(i - 1)//'_fsem.bin')
            call output_array(p%fault_dip/180.0, tidy(dir_output)//'/target_train/'//num2str(i - 1)//'_fdip.bin')
        else
            call output_array(p%image, tidy(dir_output)//'/target_valid/'//num2str(i - nt - 1)//'_dhr.bin')
            call output_array(p%rgt, tidy(dir_output)//'/target_valid/'//num2str(i - nt - 1)//'_rgt.bin')
            call output_array(p%fault, tidy(dir_output)//'/target_valid/'//num2str(i - nt - 1)//'_fsem.bin')
            call output_array(p%fault_dip/180.0, tidy(dir_output)//'/target_valid/'//num2str(i - nt - 1)//'_fdip.bin')
        end if

        ! ===============================================================
        ! Generate data for refinement NN

        ! Create noisy Fault for refinement NN
        r = random_mask_smooth(p%n1, p%n2, gs=[4.0, 4.0], mask_out=zr1(i), seed=sd(i)*5)
        iter_fsem = p%fault*r
        iter_fdip = p%fault_dip/180.0*r

        q%seed = p%seed + nm
        q%n1 = p%n1
        q%n2 = p%n2
        q%nf = nff(i)
        q%fwidth = 3.0
        call q%generate
        where (q%fault /= 0)
            q%fault = 1.0
        end where
        r = random_mask_smooth(p%n1, p%n2, gs=[4.0, 4.0], mask_out=zr2(i), seed=sd(i)*3)
        q%fault = q%fault*r
        q%fault_dip = q%fault_dip/180*r

        where (iter_fsem /= 1)
            iter_fsem = iter_fsem + q%fault
            iter_fdip = iter_fdip + q%fault_dip
        end where

        ! Create noisy RGT for refinement NN
        r = random(p%n1, p%n2, dist='normal', seed=(i + nint(i/2.0))*i)
        r = gauss_filt(r, [5.0, 5.0])
        r = r/maxval(abs(r))*rx(i)
        r = r - mean(r) + 1.0
        iter_rgt = p%rgt*r
        iter_rgt = clip(iter_rgt, 0.0, 1.0)

        ! Create noisy DHR for refinement NN
        r = random(p%n1, p%n2, dist='normal', seed=(i + nint(i/3.0))*i)
        r = gauss_filt(r, [5.0, 5.0])
        r = rescale(r, [0.0, 1.0])
        where (r < 0.2)
            r = 0.0
        end where
        where (r > 0.5)
            r = 1.0
        end where
        iter_dhr = r
        r = random(p%n1, p%n2, dist='normal', seed=(i + nint(i/2.0))**2)
        r = gauss_filt(r, [2.0, 2.0])
        r = r - mean(r)
        r = r/std(r)
        r = r/maxval(abs(r))*0.05*maxval(p%image)
        iter_dhr = iter_dhr*(p%image + r)

        ! ===============================================================
        ! Generate noisy image for inference NN

        p%f0 = fs(i)
        p%yn_rgt = .false.
        p%yn_fault = .false.
        select case (mod(i + 1, 2))
            case (0)
                p%noise_type = 'uniform'
            case (1)
                p%noise_type = 'normal'
            case (2)
                p%noise_type = 'exp'
        end select
        p%noise_smooth = [nsmooth1(i), nsmooth2(i)]
        p%psf_sigma = [12.0, psfm(i)]
        p%noise_level = noise(i)
        !        if (mod(i + 2, 3) == 0) then
        p%yn_conv_noise = .false.
        !        else
        !            p%yn_conv_noise = .true.
        !        end if
        call p%generate

        p%image = p%image/norm2(p%image)*100

        if (i <= nt) then
            call output_array(p%image, tidy(dir_output)//'/data_train/'//num2str(i - 1)//'.bin')
            call output_array(iter_dhr, tidy(dir_output)//'/data_train/'//num2str(i - 1)//'_dhr.bin')
            call output_array(iter_rgt, tidy(dir_output)//'/data_train/'//num2str(i - 1)//'_rgt.bin')
            call output_array(iter_fsem, tidy(dir_output)//'/data_train/'//num2str(i - 1)//'_fsem.bin')
            call output_array(iter_fdip, tidy(dir_output)//'/data_train/'//num2str(i - 1)//'_fdip.bin')
        else
            call output_array(p%image, tidy(dir_output)//'/data_valid/'//num2str(i - nt - 1)//'.bin')
            call output_array(iter_dhr, tidy(dir_output)//'/data_valid/'//num2str(i - nt - 1)//'_dhr.bin')
            call output_array(iter_rgt, tidy(dir_output)//'/data_valid/'//num2str(i - nt - 1)//'_rgt.bin')
            call output_array(iter_fsem, tidy(dir_output)//'/data_valid/'//num2str(i - nt - 1)//'_fsem.bin')
            call output_array(iter_fdip, tidy(dir_output)//'/data_valid/'//num2str(i - nt - 1)//'_fdip.bin')
        end if

        print *, date_time_compact(), ' >> ', num2str(i), ' of ', num2str(nm)

    end do

end program main
