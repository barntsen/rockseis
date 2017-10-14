cd Results
rsrss2rsf <Pres.rss out=stdout > Pres.rsf 
rsrss2rsf <Pmod.rss out=stdout > Pmod.rsf 
rsrss2rsf <Pshot.rss out=stdout > Prec.rsf 
rsrss2rsf <Vpgrad3d.rss out=stdout > Vpgrad3d.rsf 
rsrss2rsf <Rhograd3d.rss out=stdout > Rhograd3d.rsf 

sfgraph < Pres.rsf title=residual | sfpen &
sfgraph < Pmod.rsf title=modelled | sfpen &
sfgraph < Prec.rsf title=recorded | sfpen &
sfput < Vpgrad3d.rsf label1=x label2=y label3=z |  sfbyte mean=y | sfgrey3 frame1=22 frame2=22 frame3=22 color=j title=VP | sfpen &
sfput < Rhograd3d.rsf label1=x label2=y label3=z | sfbyte mean=y | sfgrey3  frame1=22 frame2=22 frame3=22 color=j title=RHO | sfpen &
cd ..
