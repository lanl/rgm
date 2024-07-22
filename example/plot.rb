
opts = "-label1='' -label2='' -ticktop=n -color=rainbowcmyk -tickbottom=n -tickleft=n -tickright=n "

system "x_showmatrix -in=model/vp2_1.bin -n1=400 #{opts} -out=figures/9.png &"
system "x_showmatrix -in=model/fault2_1.bin -n1=400 #{opts} -out=figures/10.png &"
system "x_showmatrix -in=model/vp2_2.bin -n1=400 #{opts} -out=figures/11.png &"
system "x_showmatrix -in=model/fault2_2.bin -n1=400 #{opts} -out=figures/12.png &"

system "x_showmatrix -in=model/img2_1.bin -n1=400 #{opts} -color=binary -clip=0.2 -out=figures/13.png &"
system "x_showmatrix -in=model/img2_2.bin -n1=400 #{opts} -color=binary -clip=0.2 -out=figures/14.png &"
system "x_showmatrix -in=model/rgt2_1.bin -n1=400 #{opts} -color=rainbowcmyk -ncolor=30 -cmin=0 -cmax=1 -out=figures/15.png "
system "x_showmatrix -in=model/rgt2_2.bin -n1=400 #{opts} -color=rainbowcmyk -ncolor=30 -cmin=0 -cmax=1 -out=figures/16.png "
