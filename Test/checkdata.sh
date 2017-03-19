rm *.rss *.rsf
../bin/rshello_models
cp Shot2d.rss Pshot2d.rss
cp Shot2d.rss Axshot2d.rss
cp Shot2d.rss Azshot2d.rss
cp Shot3d.rss Pshot3d.rss
#../bin/rshello_acoustic
#../bin/rsrss2rsf <Psnaps2d.rss out=stdout > Psnap.rsf
#../bin/rsrss2rsf <Axsnaps2d.rss out=stdout > Axsnaps.rsf
#../bin/rsrss2rsf <Azsnaps2d.rss out=stdout > Azsnaps.rsf
#../bin/rsrss2rsf <Pshot2d.rss out=stdout > Pshot.rsf
#../bin/rsrss2rsf <Axshot2d.rss out=stdout > Axshot.rsf
#../bin/rsrss2rsf <Azshot2d.rss out=stdout > Azshot.rsf
#
#../bin/rshello_elastic
#../bin/rsrss2rsf <snaps.rss out=stdout > data.rsf

valgrind ../bin/rshello3d_acoustic
#../bin/rsrss2rsf <Pshot2d.rss out=stdout > Pshot.rsf
#../bin/rsrss2rsf <Psnaps3d.rss out=stdout > Psnap.rsf

#../bin/rshello3d_elastic
#echo data_format=float in=data.bin | sfput d1=10 d2=10 d3=10 n1=65 n2=31 n3=61 n4=501 d4=1e-3 o1=0 o2=0 o3=0 o4=0 out=stdout > data.rsf
#../bin/rsrss2rsf <snaps.rss out=stdout > data.rsf
