../bin/rshello_models
valgrind ../bin/rshello_acoustic
#echo data_format=float in=data.bin | sfput d1=10 d2=10 n1=77 n2=73 n3=501 d3=1e-3 o1=0 o2=0 o3=0 out=stdout > data.rsf
valgrind ../bin/rsrss2rsf <snaps.rss out=stdout > data.rsf

#../bin/rshello_elastic
#echo data_format=float in=data.bin | sfput d1=10 d2=10 n1=77 n2=73 n3=501 d3=1e-3 o1=0 o2=0 o3=0 out=stdout > data.rsf

#../bin/rshello3d_acoustic
#echo data_format=float in=data.bin | sfput d1=10 d2=10 d3=10 n1=65 n2=31 n3=61 n4=501 d4=1e-3 o1=0 o2=0 o3=0 o4=0 out=stdout > data.rsf

#../bin/rshello3d_elastic
#echo data_format=float in=data.bin | sfput d1=10 d2=10 d3=10 n1=65 n2=31 n3=61 n4=501 d4=1e-3 o1=0 o2=0 o3=0 o4=0 out=stdout > data.rsf
