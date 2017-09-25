sh makesurvey.sh 
if [ ! -d "Results" ]; then
    mkdir Results
fi
rsAcousticmod2d mod.cfg 
