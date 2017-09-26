sh makesurvey.sh 
if [ ! -d "Results" ]; then
    mkdir Results
fi
mpirun rsmpiAcousticmod3d mod.cfg 
mpirun rsmpiAcousticrtm3d rtm.cfg 

