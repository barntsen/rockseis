sh makesurvey.sh 
if [ ! -d "Results" ]; then
    mkdir Results
fi
mpirun -np 2 rsmpiAcousticmod3d mod.cfg 
mpirun -np 2 rsmpiAcousticrtm3d rtm.cfg 

