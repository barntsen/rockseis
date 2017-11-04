sh makesurvey.sh 
if [ ! -d "Results" ]; then
    mkdir Results
fi
mpirun rsmpiAcousticmod2d mod.cfg 
