cd Results
rsrss2rsf <Simage3d.rss out=stdout > Simage3d.rsf
sfwindow <Simage3d.rsf f3=22 n3=1 | sftransp | sfimage perc=99 &
sfwindow <Simage3d.rsf f3=22 n3=1 | sftransp | sfgrey  pclip=99 | sfpen &
sfbyte <Simage3d.rsf pclip=99 | sfgrey3 frame1=22 frame2=22 frame3=22 | sfpen &

cd ..
