sh makesurvey.sh 
if [ ! -d "Results" ]; then
    mkdir Results
fi
mpirun -np 2 rsmpiElasticmod3d mod.cfg 
mpirun -np 2 rsmpiElasticrtm3d rtm.cfg 

