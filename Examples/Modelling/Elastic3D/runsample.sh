sh makesurvey.sh 
if [ ! -d "Results" ]; then
    mkdir Results
fi
rsElasticmod3d mod.cfg 
