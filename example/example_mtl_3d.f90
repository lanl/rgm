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
    use rgm

    implicit none

    type(rgm3) :: p, q
    integer :: i, ibeg, iend, nt, nv, nm
    real, allocatable, dimension(:) :: psfm, zr1, zr2, sr, rx, fs, noise, smooth, height, height2
    real, allocatable, dimension(:) :: slope1, slope2, lwv, lwv2, nsmooth1, nsmooth2, nsmooth3
    integer, allocatable, dimension(:) :: nf, nl, nff
    real, allocatable, dimension(:, :, :) :: r, iter_rgt, iter_dhr, iter_fsem, iter_fdip, iter_fstrike
    character(len=1024) :: dir_output
    integer, allocatable, dimension(:) :: sd

    dir_output = './dataset3_mtl'
    nt = 4000
    nv = 400
    nm = nt + nv

    call make_directory(tidy(dir_output)//'/data_train')
    call make_directory(tidy(dir_output)//'/data_valid')
    call make_directory(tidy(dir_output)//'/target_train')
    call make_directory(tidy(dir_output)//'/target_valid')

    fs = random(nm, range=[100.0, 180.0], seed=101)
    nf = irandom(nm, range=[1, 12], seed=1111)
    nff = irandom(nm, range=[6, 16], seed=2222)
    nl = nint(rescale(fs, [15.0, 50.0]))
    noise = random(nm, range=[0.0, 0.5], seed=3333)
    smooth = random(nm, range=[15.0, 25.0], seed=4444)
    height = random(nm, range=[2.0, 20.0], seed=5555)
    height2 = random(nm, range=[15.0, 35.0], seed=6666)
    lwv = random(nm, range=[0.0, 0.1], seed=7777)
    lwv2 = random(nm, range=[0.1, 0.2], seed=8888)
    slope1 = random(nm, range=[-25.0, 25.0], seed=9999)
    slope2 = random(nm, range=[-25.0, 25.0], seed=10999)
    nsmooth1 = random(nm, range=[0.0, 3.5], seed=111)
    nsmooth2 = random(nm, range=[0.0, 3.5], seed=222)
    nsmooth3 = random(nm, range=[0.0, 3.5], seed=333)
    where (nsmooth1 < 0.5)
        nsmooth1 = 0.0
    end where
    where (nsmooth2 < 0.5)
        nsmooth2 = 0.0
    end where
    where (nsmooth3 < 0.5)
        nsmooth3 = 0.0
    end where
    rx = random(nm, range=[0.0, 0.1], seed=444)
    sr = random(nm, range=[0.0, 0.1], seed=555)
    sd = irandom(nm, range=[1, 10*nm], seed=666)
    zr1 = random(nm, range=[0.05, 0.45], seed=777)
    zr2 = random(nm, range=[0.85, 0.95], seed=888)
    psfm = random(nm, range=[1.5, 3.0], seed=999)

    call getpar_int('ibeg', ibeg, 1)
    call getpar_int('iend', iend, nm)

    do i = ibeg, iend

        ! ===============================================================
        ! Generate labels for for inference and refinement NN

        p%seed = sd(i)

        p%n1 = 360
        p%n2 = 128
        p%n3 = 128
        p%nf = nf(i)
        p%refl_slope = [slope1(i), slope2(i)]
        p%nl = nl(i)
        p%refl_amp = [0.1, 1.0]
        p%secondary_refl_amp = sr(i)
        p%fwidth = 3.0

        if (mod(i, 4) == 0) then
            p%refl_shape = 'gaussian'
            p%refl_mu2 = [0.0, p%n2 - 1.0]
            p%refl_mu3 = [0.0, p%n3 - 1.0]
            p%refl_sigma2 = [25.0, 50.0]
            p%refl_sigma3 = [25.0, 50.0]
            p%ng = irand(range=[1, 3], seed=sd(i)*8)
            p%refl_height = [0.0, height2(i)]
            p%lwv = lwv2(i)
        else
            p%refl_shape = 'random'
            p%refl_smooth = smooth(i)
            p%refl_height = [0.0, height(i)]
            p%lwv = lwv(i)
        end if
        select case (mod(i, 3))
            case (0)
                p%wave = 'ricker'
            case (1)
                p%wave = 'gaussian_deriv'
            case (2)
                p%wave = 'ricker_deriv'
        end select

        if (mod(i, 2) == 0) then
            p%refl_amp_dist = 'uniform'
        else
            p%refl_amp_dist = 'normal'
        end if

        !        if (mod(i + 1, 10) <= 2) then
        !            p%yn_regular_fault = .true.
        !            p%nf = irand(range=[10, 25], seed=sd(i)*2)
        !            if (mod(irand(range=[1, 10], seed=sd(i)*3), 2) == 0) then
        !                p%dip = [rand(range=[100.0, 120.0], seed=sd(i) + 1), rand(range=[60.0, 80.0], seed=sd(i) + 2)]
        !                p%strike = [rand(range=[0.0, 90.0], seed=sd(i) + 3), rand(range=[90.0, 180.0], seed=sd(i) + 4)]
        !                p%rake = [rand(range=[0.0, 90.0], seed=sd(i) + 5), rand(range=[90.0, 180.0], seed=sd(i) + 6)]
        !                p%disp = [5.0, -5.0]
        !            else
        !                p%dip = [rand(range=[60.0, 80.0], seed=sd(i) + 7), rand(range=[100.0, 120.0], seed=sd(i) + 8)]
        !                p%strike = [rand(range=[90.0, 180.0], seed=sd(i) + 9), rand(range=[0.0, 90.0], seed=sd(i) + 10)]
        !                p%rake = [rand(range=[90.0, 180.0], seed=sd(i) + 11), rand(range=[0.0, 90.0], seed=sd(i) + 12)]
        !                p%disp = [-5.0, 5.0]
        !            end if
        !        else
        p%yn_regular_fault = .false.
        p%nf = nf(i)
        p%disp = [2.0, 20.0]
        p%dip = [50.0, 130.0]
        p%strike = [0.0, 180.0]
        !        end if

        !                if (mod(i + 3, 5) == 0) then
        !                    p%unconf = irand(range=[1, 2], seed=sd(i)+11)
        !                    p%unconf_amp = [0.05, 0.1]
        !                    p%unconf_z = [0.05, 0.35]
        !                else
        p%unconf = 0
        !                end if

        p%f0 = fs(i)*1.25
        p%noise_level = 0
        p%psf_sigma = [12.0, 0.0, 0.0]
        p%yn_rgt = .true.
        p%yn_fault = .true.
        p%image_smooth = [0.0, 0.0, 0.0]
        call p%generate

        p%image = p%image/norm2(p%image)*1000

        where (p%fault /= 0)
            p%fault = 1.0
        end where

        if (i <= nt) then
            call output_array(p%image, tidy(dir_output)//'/target_train/'//num2str(i - 1)//'_dhr.bin')
            call output_array(p%rgt, tidy(dir_output)//'/target_train/'//num2str(i - 1)//'_rgt.bin')
            call output_array(p%fault, tidy(dir_output)//'/target_train/'//num2str(i - 1)//'_fsem.bin')
            call output_array(p%fault_dip/180.0, tidy(dir_output)//'/target_train/'//num2str(i - 1)//'_fdip.bin')
            call output_array(p%fault_strike/180.0, tidy(dir_output)//'/target_train/'//num2str(i - 1)//'_fstrike.bin')
        else
            call output_array(p%image, tidy(dir_output)//'/target_valid/'//num2str(i - nt - 1)//'_dhr.bin')
            call output_array(p%rgt, tidy(dir_output)//'/target_valid/'//num2str(i - nt - 1)//'_rgt.bin')
            call output_array(p%fault, tidy(dir_output)//'/target_valid/'//num2str(i - nt - 1)//'_fsem.bin')
            call output_array(p%fault_dip/180.0, tidy(dir_output)//'/target_valid/'//num2str(i - nt - 1)//'_fdip.bin')
            call output_array(p%fault_strike/180.0, tidy(dir_output)//'/target_valid/'//num2str(i - nt - 1)//'_fstrike.bin')
        end if

        ! ===============================================================
        ! Generate data for refinement NN

        ! Create noisy Fault for refinement NN
        r = random_mask_smooth(p%n1, p%n2, p%n3, gs=[4.0, 4.0, 4.0], mask_out=zr1(i), seed=sd(i)*5)
        iter_fsem = p%fault*r
        iter_fdip = p%fault_dip/180.0*r
        iter_fstrike = p%fault_strike/180.0*r

        q%seed = p%seed + nm
        q%n1 = p%n1
        q%n2 = p%n2
        q%n3 = p%n3
        q%nf = nff(i)
        q%dip = [40.0, 130.0]
        q%strike = [0.0, 180.0]
        q%rake = [0.0, 180.0]
        q%disp = p%disp
        q%fwidth = p%fwidth
        call q%generate
        where (q%fault /= 0)
            q%fault = 1.0
        end where
        r = random_mask_smooth(p%n1, p%n2, p%n3, gs=[4.0, 4.0, 4.0], mask_out=zr2(i), seed=sd(i)*3)
        q%fault = q%fault*r
        q%fault_dip = q%fault_dip/180*r
        q%fault_strike = q%fault_strike/180*r

        where (iter_fsem /= 1)
            ! Must only add extra fault pixels to zero pixels
            iter_fsem = iter_fsem + q%fault
            iter_fdip = iter_fdip + q%fault_dip
            iter_fstrike = iter_fstrike + q%fault_strike
        end where

        ! Create noisy RGT for refinement NN
        r = random(p%n1, p%n2, p%n3, dist='normal', seed=(i + nint(i/2.0))*i)
        r = gauss_filt(r, [5.0, 5.0, 5.0])
        r = r/maxval(abs(r))*rx(i)
        r = r - mean(r) + 1.0
        iter_rgt = p%rgt*r
        iter_rgt = clip(iter_rgt, 0.0, 1.0)

        ! Create noisy DHR for refinement NN
        r = random(p%n1, p%n2, p%n3, dist='normal', seed=(i + nint(i/3.0))*i)
        r = gauss_filt(r, [5.0, 5.0, 5.0])
        r = rescale(r, [0.0, 1.0])
        where (r < 0.2)
            r = 0.0
        end where
        where (r > 0.5)
            r = 1.0
        end where
        iter_dhr = r
        r = random(p%n1, p%n2, p%n3, dist='normal', seed=(i + nint(i/2.0))**2)
        r = gauss_filt(r, [2.0, 2.0, 2.0])
        r = r - mean(r)
        r = r/std(r)
        r = r/maxval(abs(r))*0.05*maxval(p%image)
        iter_dhr = iter_dhr*(p%image + r)

        ! ===============================================================
        ! Generate noisy image for inference NN

        p%f0 = fs(i)
        p%yn_rgt = .false.
        p%yn_fault = .false.
        select case (mod(i + 2, 3))
            case (0)
                p%noise_type = 'uniform'
            case (1)
                p%noise_type = 'normal'
            case (2)
                p%noise_type = 'exp'
        end select
        p%noise_smooth = [nsmooth1(i), nsmooth2(i), nsmooth3(i)]
        p%psf_sigma = [12.0, psfm(i), psfm(i)]
        p%noise_level = noise(i)
        !        if (mod(i, 5) == 0) then
        p%yn_conv_noise = .false.
        !        else
        !            p%yn_conv_noise = .true.
        !        end if
        call p%generate

        p%image = p%image/norm2(p%image)*1000

        if (i <= nt) then
            call output_array(p%image, tidy(dir_output)//'/data_train/'//num2str(i - 1)//'.bin')
            call output_array(iter_dhr, tidy(dir_output)//'/data_train/'//num2str(i - 1)//'_dhr.bin')
            call output_array(iter_rgt, tidy(dir_output)//'/data_train/'//num2str(i - 1)//'_rgt.bin')
            call output_array(iter_fsem, tidy(dir_output)//'/data_train/'//num2str(i - 1)//'_fsem.bin')
            call output_array(iter_fdip, tidy(dir_output)//'/data_train/'//num2str(i - 1)//'_fdip.bin')
            call output_array(iter_fstrike, tidy(dir_output)//'/data_train/'//num2str(i - 1)//'_fstrike.bin')
        else
            call output_array(p%image, tidy(dir_output)//'/data_valid/'//num2str(i - nt - 1)//'.bin')
            call output_array(iter_dhr, tidy(dir_output)//'/data_valid/'//num2str(i - nt - 1)//'_dhr.bin')
            call output_array(iter_rgt, tidy(dir_output)//'/data_valid/'//num2str(i - nt - 1)//'_rgt.bin')
            call output_array(iter_fsem, tidy(dir_output)//'/data_valid/'//num2str(i - nt - 1)//'_fsem.bin')
            call output_array(iter_fdip, tidy(dir_output)//'/data_valid/'//num2str(i - nt - 1)//'_fdip.bin')
            call output_array(iter_fstrike, tidy(dir_output)//'/data_valid/'//num2str(i - nt - 1)//'_fstrike.bin')
        end if
        print *, date_time_compact(), ' >> ', num2str(i), ' of ', num2str(nm)

    end do

end program main
