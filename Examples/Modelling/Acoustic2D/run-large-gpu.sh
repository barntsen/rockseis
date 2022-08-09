#!/bin/sh

./clean.sh

echo "** Creating input files for large model **"
# Creating the input rss files
python3 large.py


echo "** Start gpu modelling for large model **"
# Running single node acoustic cpu modeling
/usr/bin/time --format='Wall time: %e sec(s)' rsAcousticmod2dgpu mod-large.cfg
#nvprof --unified-memory-profiling off rsAcousticmod2dgpu mod-large.cfg


# Plotting seismograms
python3 Ploting_seismo.py

