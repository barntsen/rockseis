sh makesurvey.sh 
if [ ! -d "Results" ]; then
    mkdir Results
fi
rsAcousticmod3d mod.cfg 
