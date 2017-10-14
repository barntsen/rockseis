sh makesurvey.sh 
if [ ! -d "Results" ]; then
    mkdir Results
fi
# Make perturbed model 
rsrss2rsf <../../Models/Vp3d.rss out=stdout > Vp3d.rsf
sfspike <Vp3d.rsf k3=25 l3=26 k1=30 k2=22 l2=30 l1=35 mag=200. | sfsmooth rect1=5 rect2=5 rect3=5 out=stdout > pert.rsf
sfmath <Vp3d.rsf x=pert.rsf output="in+x" out=stdout > Vp3d_pert.rsf 
rsrsf2rss <Vp3d_pert.rsf > Vp3d_pert.rss

rsrss2rsf <../../Models/Rho3d.rss out=stdout > Rho3d.rsf
sfspike <Rho3d.rsf k3=25 l3=26 mag=0. | sfsmooth rect3=5 out=stdout > pert.rsf
sfmath <Rho3d.rsf x=pert.rsf output="in+x" out=stdout > Rho3d_pert.rsf 
rsrsf2rss <Rho3d_pert.rsf > Rho3d_pert.rss

mpirun -np 2 rsmpiAcousticmod3d mod.cfg 
mpirun -np 2 rsmpiAcousticfwigrad3d fwi.cfg 
