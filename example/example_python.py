
#
# Examples of the RGM v2.0 features through the Python interface: folded
# layers, faults, unconformities, salt bodies, karst cave systems, and elastic
# models. The figures are plotted with pymplot and written to the current
# directory.
#

import rgm
from pymplot import showmatrix, showslice

# ==================================================================
# Layers and folds

# Folded layers without faults (nf = 0), with RGT and facies
p = rgm.rgm2(n1=151, n2=301, nl=25, nf=0, seed=101, lwv=0.4, lwh=0.2,
             refl_shape='perlin', refl_shape_top='perlin', refl_smooth=3, refl_smooth_top=3,
             refl_height=[0.0, 25.0], refl_height_top=[0.0, 4.0],
             noise_level=0.01, psf_sigma=[10.0, 1.0], yn_rgt=True, yn_facies=True)
p.generate()
showmatrix(p.vp, colormap='jet', legend=True, unit='Vp (m/s)', label1='Z', label2='X',
           outfile='layers_2d_vp.png')
showmatrix(p.image, colormap='binary', cperc=99, label1='Z', label2='X',
           outfile='layers_2d_image.png')
showmatrix(p.rgt, colormap='jet', legend=True, unit='RGT', label1='Z', label2='X',
           outfile='layers_2d_rgt.png')
showmatrix(p.facies, colormap='jet', legend=True, unit='Facies', label1='Z', label2='X',
           outfile='layers_2d_facies.png')

# Gaussian anticlines and synclines along rotated axes
q = rgm.rgm3(n1=128, n2=160, n3=160, nl=25, nf=0, seed=102, lwv=0.4, lwh=0.2,
             refl_shape='gaussian', refl_shape_top='perlin', refl_smooth=0, refl_smooth_top=3,
             refl_height=[0.0, 50.0], refl_height_top=[0.0, 4.0],
             ng=3, rotate_fold=True, refl_sigma2=[30.0, 60.0], refl_sigma3=[30.0, 60.0],
             noise_level=0.01, psf_sigma=[10.0, 1.0, 1.0], yn_rgt=True)
q.generate()
showslice(q.vp, colormap='jet', legend=True, unit='Vp (m/s)', label1='Z', label2='X', label3='Y',
          outfile='fold_3d_vp.png')
showslice(q.image, colormap='binary', cperc=99, label1='Z', label2='X', label3='Y',
          outfile='fold_3d_image.png')
showslice(q.rgt, colormap='jet', legend=True, unit='RGT', label1='Z', label2='X', label3='Y',
          outfile='fold_3d_rgt.png')

# ==================================================================
# Faults; the fault labels are displayed on top of the images

# Listric normal and reverse faults in 2D
p = rgm.rgm2(n1=151, n2=301, nl=25, nf=4, seed=201, lwv=0.4, lwh=0.2,
             refl_shape='perlin', refl_shape_top='perlin', refl_smooth=3, refl_smooth_top=3,
             refl_height=[0.0, 15.0], refl_height_top=[0.0, 4.0],
             dip=[60.0, 120.0], disp=[5.0, 12.0], delta_dip=[0.0, 20.0],
             noise_level=0.01, psf_sigma=[10.0, 1.0])
p.generate()
showmatrix(p.vp, colormap='jet', legend=True, unit='Vp (m/s)', label1='Z', label2='X',
           outfile='fault_2d_vp.png')
showmatrix(p.image, colormap='binary', cperc=99, label1='Z', label2='X',
           outfile='fault_2d_image.png')
showmatrix(p.fault, background=p.image, backcolormap='binary', backcperc=99, colormap='jet', alphas='0:0,0.5:1',
           label1='Z', label2='X', outfile='fault_2d_fault.png')

# Faults in 3D with prescribed dip, strike and rake
q = rgm.rgm3(n1=128, n2=160, n3=160, nl=25, nf=5, seed=202, lwv=0.4, lwh=0.2,
             refl_shape='perlin', refl_shape_top='perlin', refl_smooth=3, refl_smooth_top=3,
             refl_height=[0.0, 15.0], refl_height_top=[0.0, 4.0],
             dip=[60.0, 120.0], strike=[20.0, 60.0], rake=[0.0, 30.0], disp=[8.0, 16.0], delta_dip=[0.0, 15.0],
             noise_level=0.01, psf_sigma=[10.0, 1.0, 1.0])
q.generate()
showslice(q.vp, colormap='jet', legend=True, unit='Vp (m/s)', label1='Z', label2='X', label3='Y',
          outfile='fault_3d_vp.png')
showslice(q.image, colormap='binary', cperc=99, label1='Z', label2='X', label3='Y',
          outfile='fault_3d_image.png')
showslice(q.fault, background=q.image, backcolormap='binary', backcperc=99, colormap='jet', alphas='0:0,0.5:1',
          label1='Z', label2='X', label3='Y', outfile='fault_3d_fault.png')

# Faults with strike varying along the fault (curved in map view)
q = rgm.rgm3(n1=128, n2=160, n3=160, nl=25, nf=4, seed=203, lwv=0.4, lwh=0.2,
             refl_shape='perlin', refl_shape_top='perlin', refl_smooth=3, refl_smooth_top=3,
             refl_height=[0.0, 15.0], refl_height_top=[0.0, 4.0],
             dip=[60.0, 120.0], disp=[10.0, 20.0], delta_dip=[0.0, 15.0],
             delta_strike=[15.0, 25.0], strike_nperiod=2,
             noise_level=0.01, psf_sigma=[10.0, 1.0, 1.0])
q.generate()
showslice(q.vp, colormap='jet', legend=True, unit='Vp (m/s)', label1='Z', label2='X', label3='Y',
          outfile='fault_strike_3d_vp.png')
showslice(q.image, colormap='binary', cperc=99, label1='Z', label2='X', label3='Y',
          outfile='fault_strike_3d_image.png')
showslice(q.fault, background=q.image, backcolormap='binary', backcperc=99, colormap='jet', alphas='0:0,0.5:1',
          label1='Z', label2='X', label3='Y', outfile='fault_strike_3d_fault.png')

# Faults with an elliptical slip patch: the displacement dies out toward the
# fault tips
q = rgm.rgm3(n1=128, n2=160, n3=160, nl=25, nf=4, seed=204, lwv=0.4, lwh=0.2,
             refl_shape='perlin', refl_shape_top='perlin', refl_smooth=3, refl_smooth_top=3,
             refl_height=[0.0, 15.0], refl_height_top=[0.0, 4.0],
             dip=[60.0, 120.0], disp=[10.0, 20.0], delta_strike=[10.0, 20.0],
             yn_vary_disp=True, disp_radius_strike=[0.4, 0.6], disp_radius_dip=[0.5, 0.8],
             disp_center_dip=[0.3, 0.6],
             noise_level=0.01, psf_sigma=[10.0, 1.0, 1.0])
q.generate()
showslice(q.vp, colormap='jet', legend=True, unit='Vp (m/s)', label1='Z', label2='X', label3='Y',
          outfile='fault_vary_disp_3d_vp.png')
showslice(q.image, colormap='binary', cperc=99, label1='Z', label2='X', label3='Y',
          outfile='fault_vary_disp_3d_image.png')
showslice(q.fault_disp, background=q.image, backcolormap='binary', backcperc=99, colormap='jet', alphas='0:0,0.5:1',
          legend=True, unit='Displacement', label1='Z', label2='X', label3='Y',
          outfile='fault_vary_disp_3d_disp.png')

# Displacement decaying away from the faults, which creates drag folds and
# rollovers
q = rgm.rgm3(n1=128, n2=160, n3=160, nl=25, nf=3, seed=205, lwv=0.4, lwh=0.2,
             refl_shape='perlin', refl_shape_top='perlin', refl_smooth=3, refl_smooth_top=3,
             refl_height=[0.0, 15.0], refl_height_top=[0.0, 4.0],
             dip=[60.0, 120.0], disp=[12.0, 22.0], delta_strike=[10.0, 20.0],
             yn_vary_disp=True, disp_radius_strike=[0.4, 0.6], disp_radius_dip=[0.5, 0.8],
             yn_disp_decay=True, disp_decay_width=[0.3, 0.5],
             noise_level=0.01, psf_sigma=[10.0, 1.0, 1.0])
q.generate()
showslice(q.vp, colormap='jet', legend=True, unit='Vp (m/s)', label1='Z', label2='X', label3='Y',
          outfile='fault_decay_3d_vp.png')
showslice(q.image, colormap='binary', cperc=99, label1='Z', label2='X', label3='Y',
          outfile='fault_decay_3d_image.png')
showslice(q.fault_disp, background=q.image, backcolormap='binary', backcperc=99, colormap='jet', alphas='0:0,0.5:1',
          legend=True, unit='Displacement', label1='Z', label2='X', label3='Y',
          outfile='fault_decay_3d_disp.png')

# ==================================================================
# Unconformities

# Two unconformities with random erosional topography
p = rgm.rgm2(n1=151, n2=301, nl=25, nf=3, seed=301, lwv=0.4, lwh=0.2,
             refl_shape='perlin', refl_shape_top='perlin', refl_smooth=3, refl_smooth_top=3,
             refl_height=[0.0, 15.0], refl_height_top=[0.0, 4.0], disp=[5.0, 10.0],
             unconf=2, unconf_z=[0.15, 0.4], unconf_height=[5.0, 15.0], unconf_nl=10,
             noise_level=0.01, psf_sigma=[10.0, 1.0], yn_rgt=True)
p.generate()
showmatrix(p.vp, colormap='jet', legend=True, unit='Vp (m/s)', label1='Z', label2='X',
           outfile='unconf_2d_vp.png')
showmatrix(p.image, colormap='binary', cperc=99, label1='Z', label2='X',
           outfile='unconf_2d_image.png')
showmatrix(p.rgt, colormap='jet', legend=True, unit='RGT', label1='Z', label2='X',
           outfile='unconf_2d_rgt.png')

# Unconformity carved by meandering river channels; the depth slice cuts
# through the channels
q = rgm.rgm3(n1=128, n2=160, n3=160, nl=25, nf=2, seed=302, lwv=0.4, lwh=0.2,
             refl_shape='perlin', refl_shape_top='perlin', refl_smooth=3, refl_smooth_top=3,
             refl_height=[0.0, 15.0], refl_height_top=[0.0, 4.0], disp=[5.0, 10.0],
             unconf=1, unconf_z=[0.2, 0.3], unconf_shape='meander_channel',
             unconf_channel_width=[0.04, 0.08], unconf_channel_sinuosity=1.2, unconf_topo=0.25,
             noise_level=0.01, psf_sigma=[10.0, 1.0, 1.0])
q.generate()
showslice(q.vp, colormap='jet', legend=True, unit='Vp (m/s)', label1='Z', label2='X', label3='Y',
          slice1=38, outfile='meander_channel_3d_vp.png')
showslice(q.image, colormap='binary', cperc=99, label1='Z', label2='X', label3='Y',
          slice1=38, outfile='meander_channel_3d_image.png')

# Unconformity carved by a meandering incised canyon
q = rgm.rgm3(n1=128, n2=160, n3=160, nl=25, nf=2, seed=303, lwv=0.4, lwh=0.2,
             refl_shape='perlin', refl_shape_top='perlin', refl_smooth=3, refl_smooth_top=3,
             refl_height=[0.0, 15.0], refl_height_top=[0.0, 4.0], disp=[5.0, 10.0],
             unconf=1, unconf_z=[0.2, 0.3], unconf_shape='meander_canyon',
             unconf_channel_width=[0.05, 0.10], unconf_height=[10.0, 20.0],
             noise_level=0.01, psf_sigma=[10.0, 1.0, 1.0])
q.generate()
showslice(q.vp, colormap='jet', legend=True, unit='Vp (m/s)', label1='Z', label2='X', label3='Y',
          slice1=38, outfile='meander_canyon_3d_vp.png')
showslice(q.image, colormap='binary', cperc=99, label1='Z', label2='X', label3='Y',
          slice1=38, outfile='meander_canyon_3d_image.png')

# Unconformity carved by a dendritic drainage network
q = rgm.rgm3(n1=128, n2=160, n3=160, nl=25, nf=2, seed=304, lwv=0.4, lwh=0.2,
             refl_shape='perlin', refl_shape_top='perlin', refl_smooth=3, refl_smooth_top=3,
             refl_height=[0.0, 15.0], refl_height_top=[0.0, 4.0], disp=[5.0, 10.0],
             unconf=1, unconf_z=[0.2, 0.3], unconf_shape='drainage_channel',
             unconf_channel_density=[0.03, 0.08],
             noise_level=0.01, psf_sigma=[10.0, 1.0, 1.0])
q.generate()
showslice(q.vp, colormap='jet', legend=True, unit='Vp (m/s)', label1='Z', label2='X', label3='Y',
          slice1=30, outfile='drainage_channel_3d_vp.png')
showslice(q.image, colormap='binary', cperc=99, label1='Z', label2='X', label3='Y',
          slice1=30, outfile='drainage_channel_3d_image.png')

# Unconformity carved by a dendritic drainage canyon system
q = rgm.rgm3(n1=128, n2=160, n3=160, nl=25, nf=2, seed=305, lwv=0.4, lwh=0.2,
             refl_shape='perlin', refl_shape_top='perlin', refl_smooth=3, refl_smooth_top=3,
             refl_height=[0.0, 15.0], refl_height_top=[0.0, 4.0], disp=[5.0, 10.0],
             unconf=1, unconf_z=[0.2, 0.3], unconf_shape='drainage_canyon',
             unconf_channel_density=[0.03, 0.08], unconf_height=[10.0, 20.0],
             noise_level=0.01, psf_sigma=[10.0, 1.0, 1.0])
q.generate()
showslice(q.vp, colormap='jet', legend=True, unit='Vp (m/s)', label1='Z', label2='X', label3='Y',
          slice1=46, outfile='drainage_canyon_3d_vp.png')
showslice(q.image, colormap='binary', cperc=99, label1='Z', label2='X', label3='Y',
          slice1=46, outfile='drainage_canyon_3d_image.png')

# ==================================================================
# Salt bodies; the salt mask is displayed on top of the image

q = rgm.rgm3(n1=128, n2=160, n3=160, nl=25, nf=2, seed=401, lwv=0.4, lwh=0.2,
             refl_shape='perlin', refl_shape_top='perlin', refl_smooth=3, refl_smooth_top=3,
             refl_height=[0.0, 15.0], refl_height_top=[0.0, 4.0], disp=[5.0, 10.0],
             yn_salt=True, nsalt=2, salt_radius=[20.0, 30.0], salt_top_z=[0.4, 0.6],
             salt_nnode=8, salt_path_variation=6.0,
             noise_level=0.01, psf_sigma=[10.0, 1.0, 1.0])
q.generate()
showslice(q.vp, colormap='jet', legend=True, unit='Vp (m/s)', label1='Z', label2='X', label3='Y',
          outfile='salt_3d_vp.png')
showslice(q.image, colormap='binary', cperc=99, label1='Z', label2='X', label3='Y',
          outfile='salt_3d_image.png')
showslice(q.salt, background=q.image, backcolormap='binary', backcperc=99, colormap='jet', alphas='0:0,0.5:1',
          label1='Z', label2='X', label3='Y', outfile='salt_3d_salt.png')

# ==================================================================
# Karst cave system: a connected network of tubes; the karst mask is
# displayed on top of the image

q = rgm.rgm3(n1=128, n2=160, n3=160, nl=25, nf=3, seed=501, lwv=0.4, lwh=0.2,
             refl_shape='perlin', refl_shape_top='perlin', refl_smooth=3, refl_smooth_top=3,
             refl_height=[0.0, 15.0], refl_height_top=[0.0, 4.0], dip=[60.0, 120.0], disp=[5.0, 10.0],
             yn_karst=True, karst_z=[0.45, 0.85], karst_npassage=25, karst_nctrl=18,
             karst_connect=0.4, karst_tortuosity=0.6,
             noise_level=0.01, psf_sigma=[10.0, 1.0, 1.0])
q.generate()
showslice(q.vp, colormap='jet', legend=True, unit='Vp (m/s)', label1='Z', label2='X', label3='Y',
          outfile='karst_3d_vp.png')
showslice(q.image, colormap='binary', cperc=99, label1='Z', label2='X', label3='Y',
          outfile='karst_3d_image.png')
showslice(q.karst, background=q.image, backcolormap='binary', backcperc=99, colormap='jet', alphas='0:0,0.5:1',
          label1='Z', label2='X', label3='Y', outfile='karst_3d_karst.png')

# ==================================================================
# Elastic model: Vp, Vs, density, and the PP, PS, SP and SS images

q = rgm.rgm3(n1=128, n2=160, n3=160, nl=25, nf=3, seed=601, lwv=0.4, lwh=0.2,
             refl_shape='perlin', refl_shape_top='perlin', refl_smooth=3, refl_smooth_top=3,
             refl_height=[0.0, 15.0], refl_height_top=[0.0, 4.0],
             dip=[60.0, 120.0], disp=[8.0, 16.0], delta_strike=[10.0, 20.0],
             yn_elastic=True, vpvsratio=[1.6, 1.9],
             noise_level=0.01, psf_sigma=[10.0, 1.0, 1.0])
q.generate()
showslice(q.vp, colormap='jet', legend=True, unit='Vp (m/s)', label1='Z', label2='X', label3='Y',
          outfile='elastic_3d_vp.png')
showslice(q.vs, colormap='jet', legend=True, unit='Vs (m/s)', label1='Z', label2='X', label3='Y',
          outfile='elastic_3d_vs.png')
showslice(q.rho, colormap='jet', legend=True, unit='Density (kg/m$^3$)', label1='Z', label2='X', label3='Y',
          outfile='elastic_3d_rho.png')
showslice(q.image_pp, colormap='binary', cperc=99, label1='Z', label2='X', label3='Y',
          outfile='elastic_3d_image_pp.png')
showslice(q.image_ps, colormap='binary', cperc=99, label1='Z', label2='X', label3='Y',
          outfile='elastic_3d_image_ps.png')
showslice(q.image_sp, colormap='binary', cperc=99, label1='Z', label2='X', label3='Y',
          outfile='elastic_3d_image_sp.png')
showslice(q.image_ss, colormap='binary', cperc=99, label1='Z', label2='X', label3='Y',
          outfile='elastic_3d_image_ss.png')
