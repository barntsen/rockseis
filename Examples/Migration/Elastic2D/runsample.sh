sh makesurvey.sh 
if [ ! -d "Results" ]; then
    mkdir Results
fi
mpirun rsmpiElasticmod2d mod.cfg 
mpirun rsmpiElasticrtm2d rtm.cfg 

