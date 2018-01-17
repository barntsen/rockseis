sh makesurvey.sh 
if [ ! -d "Results" ]; then
    mkdir Results
fi
mpirun -np 4 rsmpiElasticmod3d mod.cfg 
