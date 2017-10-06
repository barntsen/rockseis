sh makesurvey.sh 
if [ ! -d "Results" ]; then
    mkdir Results
fi
mpirun -np 2 rsmpiElasticmod2d mod.cfg 
mpirun -np 2 rsmpiElasticrtm2d rtm.cfg 

