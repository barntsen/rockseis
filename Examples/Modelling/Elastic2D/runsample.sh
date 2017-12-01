sh makesurvey.sh 
if [ ! -d "Results" ]; then
    mkdir Results
fi
rsElasticmod2d mod.cfg 
