rm *.rss *.rsf
../bin/rshello_models

#cp Shot2d.rss Pshot2d.rss
#cp Shot2d.rss Axshot2d.rss
#cp Shot2d.rss Azshot2d.rss
#../bin/rshello_acoustic
#../bin/rsrss2rsf <Psnaps2d.rss out=stdout > Psnap.rsf
#../bin/rsrss2rsf <Axsnaps2d.rss out=stdout > Axsnaps.rsf
#../bin/rsrss2rsf <Azsnaps2d.rss out=stdout > Azsnaps.rsf
#../bin/rsrss2rsf <Pshot2d.rss out=stdout > Pshot.rsf
#../bin/rsrss2rsf <Axshot2d.rss out=stdout > Axshot.rsf
#../bin/rsrss2rsf <Azshot2d.rss out=stdout > Azshot.rsf
##
##### 
#cp Shot2d.rss Pshot2d.rss
#cp Shot2d.rss Vxshot2d.rss
#cp Shot2d.rss Vzshot2d.rss
#../bin/rshello_elastic
#../bin/rsrss2rsf <Psnaps2d.rss out=stdout > Psnap.rsf
#../bin/rsrss2rsf <Vxsnaps2d.rss out=stdout > Vxsnaps.rsf
#../bin/rsrss2rsf <Vzsnaps2d.rss out=stdout > Vzsnaps.rsf
#../bin/rsrss2rsf <Pshot2d.rss out=stdout > Pshot.rsf
#../bin/rsrss2rsf <Vxshot2d.rss out=stdout > Vxshot.rsf
#../bin/rsrss2rsf <Vzshot2d.rss out=stdout > Vzshot.rsf
##
##### 
cp Shot3d.rss Pshot3d.rss
cp Shot3d.rss Ayshot3d.rss
../bin/rshello3d_acoustic
../bin/rsrss2rsf <Pshot3d.rss out=stdout > Pshot.rsf
../bin/rsrss2rsf <Psnaps3d.rss out=stdout > Psnap.rsf
#
##### 
#cp Shot3d.rss Pshot3d.rss
#cp Shot3d.rss Vxshot3d.rss
#cp Shot3d.rss Vyshot3d.rss
#cp Shot3d.rss Vzshot3d.rss
#../bin/rshello3d_elastic
#../bin/rsrss2rsf <Psnaps3d.rss out=stdout > Psnap.rsf
#../bin/rsrss2rsf <Vxsnaps3d.rss out=stdout > Vxsnaps.rsf
#../bin/rsrss2rsf <Vysnaps3d.rss out=stdout > Vysnaps.rsf
#../bin/rsrss2rsf <Vzsnaps3d.rss out=stdout > Vzsnaps.rsf
#../bin/rsrss2rsf <Pshot3d.rss out=stdout > Pshot.rsf
#../bin/rsrss2rsf <Vxshot3d.rss out=stdout > Vxshot.rsf
#../bin/rsrss2rsf <Vyshot3d.rss out=stdout > Vyshot.rsf
#../bin/rsrss2rsf <Vzshot3d.rss out=stdout > Vzshot.rsf
