sh makesurvey.sh 
if [ ! -d "Results" ]; then
    mkdir Results
fi
rsAcousticmod2d mod.cfg 
mpirun rsmpiAcousticrtm2d rtm.cfg 
