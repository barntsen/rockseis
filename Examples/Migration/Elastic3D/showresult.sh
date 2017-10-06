cd Results
rsrss2rsf <Pimage3d.rss out=stdout > Pimage3d.rsf
#sfwindow <Pimage3d.rsf f3=22 n3=1 | sftransp | sfgrey | sfpen &
sfwindow <Pimage3d.rsf f2=22 n2=1 | sftransp | sfimage perc=99 &

rsrss2rsf <Simage3d.rss out=stdout > Simage3d.rsf
#sfwindow <Simage3d.rsf f3=22 n3=1 | sftransp | sfgrey | sfpen &
sfwindow <Simage3d.rsf f2=22 n2=1 | sftransp | sfimage perc=99 &
cd ..
