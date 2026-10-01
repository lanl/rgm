
import rgm
from pymplot import showmatrix, showslice

# ==================================================================
# Reflector shapes (2D); all models share the same faults. The images are
# dip-independent, so the steep limbs are imaged. Arrays are written with .T
# in Fortran (column-major) order.

# 1 - fold train with sharp anticlines and broad synclines (crest > 0)
p = rgm.rgm2(n1=201, n2=401, nl=60, nf=3, seed=1357, dip=[60.0, 120.0], disp=[5.0, 12.0],
             lwv=0.3, lwh=0.2, refl_shape='fold', refl_shape_top='same',
             refl_height=[0.0, 40.0], refl_height_top=[0.0, 12.0],
             refl_fold_lambda=[110.0, 130.0], refl_fold_crest=[0.2, 0.25], refl_fold_vergence=[0.3, 0.4],
             yn_dip_independent=True, noise_level=0.0)
p.generate()
p.vp.T.tofile('./example_2d_vp_shape_1.bin')
p.image.T.tofile('./example_2d_image_shape_1.bin')
showmatrix(p.vp, colormap='jet', legend=True, unit='Vp (m/s)', label1='Z', label2='X',
           outfile='example_2d_vp_shape_1.png')
showmatrix(p.image, colormap='binary', cperc=99, label1='Z', label2='X',
           outfile='example_2d_image_shape_1.png')

# 2 - fold train with broad (box) anticlines and sharp synclines (crest < 0)
p = rgm.rgm2(n1=201, n2=401, nl=60, nf=3, seed=1357, dip=[60.0, 120.0], disp=[5.0, 12.0],
             lwv=0.3, lwh=0.2, refl_shape='fold', refl_shape_top='same',
             refl_height=[0.0, 40.0], refl_height_top=[0.0, 12.0],
             refl_fold_lambda=[110.0, 130.0], refl_fold_crest=[-0.25, -0.2], refl_fold_vergence=[0.3, 0.4],
             yn_dip_independent=True, noise_level=0.0)
p.generate()
p.vp.T.tofile('./example_2d_vp_shape_2.bin')
p.image.T.tofile('./example_2d_image_shape_2.bin')
showmatrix(p.vp, colormap='jet', legend=True, unit='Vp (m/s)', label1='Z', label2='X',
           outfile='example_2d_vp_shape_2.png')
showmatrix(p.image, colormap='binary', cperc=99, label1='Z', label2='X',
           outfile='example_2d_image_shape_2.png')

# 3 - strongly verging folds with planar limbs
p = rgm.rgm2(n1=201, n2=401, nl=60, nf=3, seed=1357, dip=[60.0, 120.0], disp=[5.0, 12.0],
             lwv=0.3, lwh=0.2, refl_shape='fold', refl_shape_top='same',
             refl_height=[0.0, 40.0], refl_height_top=[0.0, 12.0],
             refl_fold_lambda=[110.0, 130.0], refl_fold_crest=[0.0, 0.0], refl_fold_vergence=[0.5, 0.5],
             refl_fold_mode='limb', yn_dip_independent=True, noise_level=0.0)
p.generate()
p.vp.T.tofile('./example_2d_vp_shape_3.bin')
p.image.T.tofile('./example_2d_image_shape_3.bin')
showmatrix(p.vp, colormap='jet', legend=True, unit='Vp (m/s)', label1='Z', label2='X',
           outfile='example_2d_vp_shape_3.png')
showmatrix(p.image, colormap='binary', cperc=99, label1='Z', label2='X',
           outfile='example_2d_image_shape_3.png')

# 4 - skewed Gaussian bumps steep on the same side, on a smooth background
p = rgm.rgm2(n1=201, n2=401, nl=60, nf=3, seed=1357, dip=[60.0, 120.0], disp=[5.0, 12.0],
             lwv=0.3, lwh=0.2, refl_shape='gaussian', refl_shape_top='same',
             refl_height=[0.0, 40.0], refl_height_top=[0.0, 12.0],
             ng=3, refl_skew=[0.4, 0.7], refl_skew_common=True, refl_background=0.2,
             yn_dip_independent=True, noise_level=0.0)
p.generate()
p.vp.T.tofile('./example_2d_vp_shape_4.bin')
p.image.T.tofile('./example_2d_image_shape_4.bin')
showmatrix(p.vp, colormap='jet', legend=True, unit='Vp (m/s)', label1='Z', label2='X',
           outfile='example_2d_vp_shape_4.png')
showmatrix(p.image, colormap='binary', cperc=99, label1='Z', label2='X',
           outfile='example_2d_image_shape_4.png')

# Model 1 imaged with the default separable PSF, with which the steep limbs
# fade and resemble faults; compare with example_2d_image_shape_1
p = rgm.rgm2(n1=201, n2=401, nl=60, nf=3, seed=1357, dip=[60.0, 120.0], disp=[5.0, 12.0],
             lwv=0.3, lwh=0.2, refl_shape='fold', refl_shape_top='same',
             refl_height=[0.0, 40.0], refl_height_top=[0.0, 12.0],
             refl_fold_lambda=[110.0, 130.0], refl_fold_crest=[0.2, 0.25], refl_fold_vergence=[0.3, 0.4],
             noise_level=0.0)
p.generate()
p.image.T.tofile('./example_2d_image_shape_1_separable.bin')
showmatrix(p.image, colormap='binary', cperc=99, label1='Z', label2='X',
           outfile='example_2d_image_shape_1_separable.png')
print('2D reflector shapes done')

# ==================================================================
# Image noise (2D) on fold model 2

# 1 - no noise
p = rgm.rgm2(n1=201, n2=401, nl=60, nf=3, seed=1357, dip=[60.0, 120.0], disp=[5.0, 12.0],
             lwv=0.3, lwh=0.2, refl_shape='fold', refl_shape_top='same',
             refl_height=[0.0, 40.0], refl_height_top=[0.0, 12.0],
             refl_fold_lambda=[110.0, 130.0], refl_fold_crest=[-0.25, -0.2], refl_fold_vergence=[0.3, 0.4],
             yn_dip_independent=True, noise_level=0.0)
p.generate()
p.image.T.tofile('./example_2d_image_noise_1.bin')
showmatrix(p.image, colormap='binary', cperc=99, label1='Z', label2='X',
           outfile='example_2d_image_noise_1.png')

# 2 - random noise
p = rgm.rgm2(n1=201, n2=401, nl=60, nf=3, seed=1357, dip=[60.0, 120.0], disp=[5.0, 12.0],
             lwv=0.3, lwh=0.2, refl_shape='fold', refl_shape_top='same',
             refl_height=[0.0, 40.0], refl_height_top=[0.0, 12.0],
             refl_fold_lambda=[110.0, 130.0], refl_fold_crest=[-0.25, -0.2], refl_fold_vergence=[0.3, 0.4],
             yn_dip_independent=True, noise_type='normal', noise_level=0.5)
p.generate()
p.image.T.tofile('./example_2d_image_noise_2.bin')
showmatrix(p.image, colormap='binary', cperc=99, label1='Z', label2='X',
           outfile='example_2d_image_noise_2.png')

# 3 - migration-like noise: reflector-parallel worms, steeply dipping swings
#     (crossing in both directions by default) and band-limited background noise
p = rgm.rgm2(n1=201, n2=401, nl=60, nf=3, seed=1357, dip=[60.0, 120.0], disp=[5.0, 12.0],
             lwv=0.3, lwh=0.2, refl_shape='fold', refl_shape_top='same',
             refl_height=[0.0, 40.0], refl_height_top=[0.0, 12.0],
             refl_fold_lambda=[110.0, 130.0], refl_fold_crest=[-0.25, -0.2], refl_fold_vergence=[0.3, 0.4],
             yn_dip_independent=True, noise_type='migration', noise_level=0.5)
p.generate()
p.image.T.tofile('./example_2d_image_noise_3.bin')
showmatrix(p.image, colormap='binary', cperc=99, label1='Z', label2='X',
           outfile='example_2d_image_noise_3.png')

# 4 - migration-like noise, plus a wavelet frequency decreasing with depth,
#     uneven illumination and trace jitter
p = rgm.rgm2(n1=201, n2=401, nl=60, nf=3, seed=1357, dip=[60.0, 120.0], disp=[5.0, 12.0],
             lwv=0.3, lwh=0.2, refl_shape='fold', refl_shape_top='same',
             refl_height=[0.0, 40.0], refl_height_top=[0.0, 12.0],
             refl_fold_lambda=[110.0, 130.0], refl_fold_crest=[-0.25, -0.2], refl_fold_vergence=[0.3, 0.4],
             yn_dip_independent=True, noise_type='migration', noise_level=0.5,
             f0_bottom=100.0, illum_level=0.35, jitter_shift=[0.8, 0.25], jitter_gain=0.05)
p.generate()
p.image.T.tofile('./example_2d_image_noise_4.bin')
showmatrix(p.image, colormap='binary', cperc=99, label1='Z', label2='X',
           outfile='example_2d_image_noise_4.png')

# 5 - migration-like noise with swings dipping in a single direction
p = rgm.rgm2(n1=201, n2=401, nl=60, nf=3, seed=1357, dip=[60.0, 120.0], disp=[5.0, 12.0],
             lwv=0.3, lwh=0.2, refl_shape='fold', refl_shape_top='same',
             refl_height=[0.0, 40.0], refl_height_top=[0.0, 12.0],
             refl_fold_lambda=[110.0, 130.0], refl_fold_crest=[-0.25, -0.2], refl_fold_vergence=[0.3, 0.4],
             yn_dip_independent=True, noise_type='migration', noise_level=0.5, noise_swing_direction='single')
p.generate()
p.image.T.tofile('./example_2d_image_noise_5.bin')
showmatrix(p.image, colormap='binary', cperc=99, label1='Z', label2='X',
           outfile='example_2d_image_noise_5.png')
print('2D image noise done')

# ==================================================================
# 3D, with strike-varying faults; dip-independent images with migration-like noise

# Fold train: the anticlines plunge out along strike and the fold axes bend
q = rgm.rgm3(n1=151, n2=251, n3=251, nl=45, nf=3, seed=2468, dip=[60.0, 120.0], disp=[4.0, 10.0],
             delta_strike=[15.0, 25.0], lwv=0.3, refl_shape='fold', refl_shape_top='same',
             refl_height=[0.0, 45.0], refl_height_top=[0.0, 15.0],
             refl_fold_lambda=[70.0, 90.0], refl_fold_crest=[0.15, 0.25], refl_fold_vergence=[0.3, 0.3],
             refl_fold_strike=[65.0, 75.0], refl_fold_plunge=0.8, refl_fold_wobble=0.2,
             yn_dip_independent=True, noise_type='migration', noise_level=0.4,
             f0_bottom=110.0, illum_level=0.3, jitter_shift=[0.6, 0.2], jitter_gain=0.04)
q.generate()
q.vp.T.tofile('./example_3d_vp_fold.bin')
q.image.T.tofile('./example_3d_image_fold.bin')
showslice(q.vp, colormap='jet', legend=True, unit='Vp (m/s)', label1='Z', label2='X', label3='Y',
          slice1=120, outfile='example_3d_vp_fold.png')
showslice(q.image, colormap='binary', cperc=99, label1='Z', label2='X', label3='Y',
          slice1=120, outfile='example_3d_image_fold.png')
print('3D fold train done')

# Skewed Gaussian ridges, elongated along rotated axes and steep on the same side
q = rgm.rgm3(n1=151, n2=251, n3=251, nl=45, nf=3, seed=2468, dip=[60.0, 120.0], disp=[4.0, 10.0],
             delta_strike=[15.0, 25.0], lwv=0.3, refl_shape='gaussian', refl_shape_top='same',
             refl_height=[0.0, 45.0], refl_height_top=[0.0, 15.0],
             rotate_fold=True, ng=4, refl_sigma2=[15.0, 25.0], refl_sigma3=[50.0, 90.0],
             refl_skew=[0.4, 0.7], refl_skew_common=True, refl_background=0.2,
             yn_dip_independent=True, noise_type='migration', noise_level=0.4,
             f0_bottom=110.0, illum_level=0.3, jitter_shift=[0.6, 0.2], jitter_gain=0.04)
q.generate()
q.vp.T.tofile('./example_3d_vp_skew.bin')
q.image.T.tofile('./example_3d_image_skew.bin')
showslice(q.vp, colormap='jet', legend=True, unit='Vp (m/s)', label1='Z', label2='X', label3='Y',
          slice1=120, outfile='example_3d_vp_skew.png')
showslice(q.image, colormap='binary', cperc=99, label1='Z', label2='X', label3='Y',
          slice1=120, outfile='example_3d_image_skew.png')
print('3D skewed ridges done')

# Box folds, 128 x 128 x 128: broad anticlines and sharp synclines (crest < 0),
# with swings dipping in a single direction
q = rgm.rgm3(n1=128, n2=128, n3=128, nl=40, nf=3, seed=1357, dip=[60.0, 120.0], disp=[3.0, 8.0],
             delta_strike=[15.0, 25.0], lwv=0.3, refl_shape='fold', refl_shape_top='same',
             refl_height=[0.0, 30.0], refl_height_top=[0.0, 10.0],
             refl_fold_lambda=[45.0, 60.0], refl_fold_crest=[-0.25, -0.2], refl_fold_vergence=[0.4, 0.4],
             refl_fold_strike=[60.0, 80.0],
             yn_dip_independent=True, noise_type='migration', noise_level=0.4, noise_swing_direction='single',
             f0_bottom=110.0, illum_level=0.3, jitter_shift=[0.6, 0.2], jitter_gain=0.04)
q.generate()
q.vp.T.tofile('./example_3d_vp_box.bin')
q.image.T.tofile('./example_3d_image_box.bin')
showslice(q.vp, colormap='jet', legend=True, unit='Vp (m/s)', label1='Z', label2='X', label3='Y',
          slice1=90, outfile='example_3d_vp_box.png')
showslice(q.image, colormap='binary', cperc=99, label1='Z', label2='X', label3='Y',
          slice1=90, outfile='example_3d_image_box.png')
print('3D box folds done')

# Skewed Cauchy domes, 128 x 128 x 128, each steep on its own random side, on a
# smooth background
q = rgm.rgm3(n1=128, n2=128, n3=128, nl=40, nf=3, seed=8642, dip=[60.0, 120.0], disp=[3.0, 8.0],
             delta_strike=[15.0, 25.0], lwv=0.3, refl_shape='cauchy', refl_shape_top='same',
             refl_height=[0.0, 30.0], refl_height_top=[0.0, 10.0],
             rotate_fold=True, ng=3, refl_sigma2=[12.0, 20.0], refl_sigma3=[25.0, 40.0],
             refl_skew=[0.3, 0.6], refl_background=0.25,
             yn_dip_independent=True, noise_type='migration', noise_level=0.4,
             f0_bottom=110.0, illum_level=0.3, jitter_shift=[0.6, 0.2], jitter_gain=0.04)
q.generate()
q.vp.T.tofile('./example_3d_vp_cauchy.bin')
q.image.T.tofile('./example_3d_image_cauchy.bin')
showslice(q.vp, colormap='jet', legend=True, unit='Vp (m/s)', label1='Z', label2='X', label3='Y',
          slice1=90, outfile='example_3d_vp_cauchy.png')
showslice(q.image, colormap='binary', cperc=99, label1='Z', label2='X', label3='Y',
          slice1=90, outfile='example_3d_image_cauchy.png')
print('3D skewed Cauchy domes done')
