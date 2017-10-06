sh makesurvey.sh 
if [ ! -d "Results" ]; then
    mkdir Results
fi
mpirun rsmpiAcousticmod2d mod.cfg 
mpirun rsmpiAcousticrtm2d rtm.cfg 
