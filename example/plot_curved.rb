
# =============================================================================
# Types of unconformities

opts = " -n1=128 -label1='Z' -label2='X' -color=jet -tick1d=50 -tick2d=50 -mtick1=4 -mtick2=4 -cmin=2000 -cmax=4000 -ctruncend=0.9 "

for i  in 1..6

    system "x_showmatrix -in=example_2d_unconf_type_#{i}.bin #{opts} -out=example_2d_unconf_type_#{i}.pdf &"

end


# =============================================================================
# Types of random interfaces/surfaces

opts = " -n1=200 -label1='Y' -label2='X' -legend=y -lloc=bottom -color=jet -tick1d=50 -tick2d=50 -mtick1=4 -mtick2=4 "

for c in ['random', 'gaussian', 'cauchy', 'perlin']

    system "x_showcontour -in=#{c}_2d.bin -background=#{c}_2d.bin -backcolor=jet -backlegend=y -clabelsize=0 -contourlevel=0.1 #{opts} -out=surface_#{c}.pdf -unit='Normalized Amplitude' -lloc=right &"

    system "x_showgraph -in=#{c}_1d.bin -size2=#{5*2.0/3.0} -size1=5 -label1='X' -linewidth=2 -label2='Normalized Amplitude' -x2beg=-0.5 -x2end=1.5 -tick2d=0.5 -mtick2=4 -tick1d=50 -mtick1=4 -label1loc=top -label1pad=10 -ticktop=y -tickbottom=n -out=interface_#{c}.pdf & "

end


# =============================================================================
# Scalability test

system "x_showgraph -in=ctime2.bin -n1=20,20,20,20,20,20,20,20,20,20 -label1='Number of Faults' -tick1d=5 -o1=1 -mtick1=4 -linewidth=2,2,2,2,2,2,2,2,2,2,2,2 -label2='CPU Wallclock Time (s)' -size1=4 -size2=5 -x1end=30 -x2beg=-0.01 -tick2beg=-0.05 -tick2d=0.05 -x2end=0.2 -mtick2=4 -plotlabel='$N = 50$':'$N = 100$':'$N = 150$':'$N = 200$':'$N = 250$':'$N = 300$':'$N = 350$':'$N = 400$':'$N = 450$':'$N = 500$' -out=ctime2.pdf "

system "x_showgraph -in=ctime3.bin -n1=20,20,20,20,20,20,20,20,20,20 -label1='Number of Faults' -tick1d=5 -o1=1 -mtick1=4 -linewidth=2,2,2,2,2,2,2,2,2,2,2,2 -label2='CPU Wallclock Time (s)' -size1=4 -size2=5 -x1end=30 -x2beg=-4 -tick2beg=-25 -tick2d=25 -x2end=150 -mtick2=4 -plotlabel='$N = 50$':'$N = 100$':'$N = 150$':'$N = 200$':'$N = 250$':'$N = 300$':'$N = 350$':'$N = 400$':'$N = 450$':'$N = 500$' -out=ctime3.pdf "

system "x_showmatrix -o2=50 -d2=50 -tick2d=50 -o1=1 -tick1d=1 -in=ctime2.bin -n1=20 -size1=4 -size2=5 -label1='Number of Faults' -label2='Number of Grid Points Per Dimension' -ctruncend=0.9 -legend=y -unit='CPU Wallclock Time (s)' -ld=0.03 -lmtick=2 -lloc=bottom -out=ctime2_image.pdf &"

system "x_showmatrix -o2=50 -d2=50 -tick2d=50 -o1=1 -tick1d=1 -in=ctime3.bin -n1=20 -size1=4 -size2=5 -label1='Number of Faults' -label2='Number of Grid Points Per Dimension' -ctruncend=0.9 -legend=y -unit='CPU Wallclock Time (s)' -lloc=bottom -ld=25 -lmtick=4  -out=ctime3_image.pdf &"


# =============================================================================
# Medium property models

opts = " -n1=201 -label1='Z' -label2='X' -color=jet -tick1d=50 -tick2d=50 -mtick1=4 -mtick2=4 -cmin=2000 -cmax=4000 -legend=y -lloc=bottom -lwidth=4 -unit='P-wave Velocity (m/s)' -ld=1000 -lmtick=9 "

for i  in 1..4

    system "x_showmatrix -in=example_2d_vp_unfaulted_#{i}.bin #{opts} -out=example_2d_vp_unfaulted_#{i}.pdf &"
    system "x_showmatrix -in=example_2d_vp_faulted_#{i}.bin #{opts} -out=example_2d_vp_faulted_#{i}.pdf &"

    system "x_showmatrix -in=example_2d_vp_salt_#{i}.bin #{opts} -cmax=5000 -out=example_2d_vp_salt_#{i}.pdf &"
    system "x_showmatrix -in=example_2d_vp_unconf_#{i}.bin #{opts} -out=example_2d_vp_unconf_#{i}.pdf &"

end


# =============================================================================
# RGT models

opts = " -n1=201 -label1='Z' -label2='X' -color=jet -tick1d=50 -tick2d=50 -mtick1=4 -mtick2=4 -cmin=0 -cmax=1 -legend=y -lloc=bottom -lwidth=4 -unit='Relative Geological Time' -ld=0.2 -lmtick=1 "

for i  in 1..4

    system "x_showmatrix -in=example_2d_rgt_unfaulted_#{i}.bin #{opts} -out=example_2d_rgt_unfaulted_#{i}.pdf &"
    system "x_showmatrix -in=example_2d_rgt_faulted_#{i}.bin #{opts} -out=example_2d_rgt_faulted_#{i}.pdf &"

    system "x_showmatrix -in=example_2d_rgt_unconf_#{i}.bin #{opts} -out=example_2d_rgt_unconf_#{i}.pdf &"

end


# =============================================================================
# Seismic images

opts = " -n1=201 -label1='Z' -label2='X' -color=binary -tick1d=50 -tick2d=50 -mtick1=4 -mtick2=4 -clip=0.03 -legend=y -lloc=bottom -lwidth=4 -unit='Reflectivity' -ld=0.015 -lmtick=2 "

for i  in 1..4

    system "x_showmatrix -in=example_2d_image_unfaulted_#{i}.bin #{opts} -out=example_2d_image_unfaulted_#{i}.pdf &"
    system "x_showmatrix -in=example_2d_image_faulted_#{i}.bin #{opts} -out=example_2d_image_faulted_#{i}.pdf &"

    system "x_showmatrix -in=example_2d_image_salt_#{i}.bin #{opts} -out=example_2d_image_salt_#{i}.pdf &"
    system "x_showmatrix -in=example_2d_image_unconf_#{i}.bin #{opts}  -out=example_2d_image_unconf_#{i}.pdf &"

end


# =============================================================================
# Fault attributes

opts = " -n1=201 -label1='Z' -label2='X' -legend=y -lloc=bottom -color=jet -tick1d=50 -tick2d=50 -mtick1=4 -mtick2=4 -lwidth=4 "

for i  in 1..4

    system "x_showmatrix -unit='Fault Index' -in=example_2d_fault_#{i}.bin #{opts} -cmin=0 -ld=2 -lmtick=1 -alpha=0:0,0.1:1 -out=example_2d_fault_#{i}.pdf &"
    system "x_showmatrix -unit='Dip (degree)' -in=example_2d_fdip_#{i}.bin #{opts} -ld=30 -lmtick=5 -alpha=30:0,31:1 -cmin=30 -cmax=150 -out=example_2d_fdip_#{i}.pdf &"

    system "x_showmatrix -unit='Dip (degree)' -in=example_2d_fdip_unconf_#{i}.bin #{opts} -ld=30 -lmtick=5 -alpha=30:0,31:1 -cmin=30 -cmax=150 -out=example_2d_fdip_unconf_#{i}.pdf &"

end


# =============================================================================
# Batch labeled data

for i in 1..16

    system "x_showslice -in=batch_fault_#{i}.bin -background=batch_image_#{i}.bin -n1=128 -n2=128 -render=3d -slice1=100 -slice2=100 -slice3=100  -color=jet -backcolor=binary -backclip=0.075 -alpha=0:0,0.1:1 -ctruncend=0.9 -out=batch_#{i}.png; \
        x_trim batch_#{i}.png &"

end


# =============================================================================
# Elastic models

opts = " -n1=201 -label1='Z' -label2='X' -legend=y -lloc=bottom -color=jet -tick1d=50 -tick2d=50 -mtick1=4 -mtick2=4 -lwidth=4 "

for i  in 1..1


    system "x_showmatrix -in=example_2d_vp_elastic_#{1}.bin #{opts} -unit='P-wave Velocity (m/s)' -out=example_2d_vp_elastic_#{i}.pdf &"
    system "x_showmatrix -in=example_2d_vs_elastic_#{1}.bin #{opts} -unit='S-wave Velocity (m/s)' -out=example_2d_vs_elastic_#{i}.pdf &"
    system "x_showmatrix -in=example_2d_rho_elastic_#{1}.bin #{opts} -unit='Density (kg/m$^3$)' -out=example_2d_rho_elastic_#{i}.pdf &"


    for c in ['pp', 'ps', 'sp', 'ss']

        system "x_showmatrix -unit='Dip (degree)' -in=example_2d_fdip_elastic_#{i}.bin #{opts} -ld=30 -lmtick=5 -alpha=30:0,31:1 -cmin=30 -cmax=150 -background=example_2d_image_#{c}_elastic_#{i}.bin -backcolor=binary -backclip=0.03 -out=example_2d_fdip_elastic_#{i}_#{c}.pdf &"

        if c == 'ps' or c == 'sp'
            clip = 3
        else
            clip = 10
        end
        system "x_showmatrix -in=spec_#{c}.bin #{opts} -label1='Normalized $k_z$' -label2='Normalized $k_x$' -d1=#{1.0/201} -d2=#{1.0/301} -o1=-0.5 -o2=-0.5 -tick1d=0.25 -tick2d=0.25 -mtick1=4 -mtick2=4 -unit='Spectrum Magnitude' -cmin=0 -cmax=#{clip} -color=magma -out=spec_#{c}.pdf &  "

    end

end





