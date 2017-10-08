sh makesurvey.sh 
if [ ! -d "Results" ]; then
    mkdir Results
fi
valgrind rsElasticmod2d mod.cfg 
