cd Results
rsrss2rsf <Pimage2d.rss out=stdout > Pimage2d.rsf
sfwindow <Pimage2d.rsf | sftransp | sfgrey | sfpen &

rsrss2rsf <Simage2d.rss out=stdout > Simage2d.rsf
sfwindow <Simage2d.rsf | sftransp | sfgrey | sfpen &
sfattr < Pimage2d.rsf
sfattr < Simage2d.rsf
cd ..
