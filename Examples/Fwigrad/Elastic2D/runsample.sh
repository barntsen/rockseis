sh makesurvey.sh 
if [ ! -d "Results" ]; then
    mkdir Results
fi
# Make perturbed model 
rsrss2rsf <../../Models/Vp2d.rss out=stdout > Vp2d.rsf
sfspike <Vp2d.rsf k3=25 l3=26 k1=30 l1=35 mag=-100. | sfsmooth rect1=5 rect3=5 out=stdout > pert.rsf
sfmath <Vp2d.rsf x=pert.rsf output="in+x" out=stdout > Vp2d_pert.rsf 
rsrsf2rss <Vp2d_pert.rsf > Vp2d_pert.rss

rsrss2rsf <../../Models/Vs2d.rss out=stdout > Vs2d.rsf
sfspike <Vs2d.rsf k3=25 l3=26 k1=30 l1=35 mag=+100. | sfsmooth rect1=5 rect3=5 out=stdout > pert.rsf
sfmath <Vs2d.rsf x=pert.rsf output="in+x" out=stdout > Vs2d_pert.rsf 
rsrsf2rss <Vs2d_pert.rsf > Vs2d_pert.rss

rsrss2rsf <../../Models/Rho2d.rss out=stdout > Rho2d.rsf
sfspike <Rho2d.rsf k3=25 l3=26 k1=30 l1=35 mag=0. | sfsmooth rect1=5 rect3=5 out=stdout > pert.rsf
sfmath <Rho2d.rsf x=pert.rsf output="in+x" out=stdout > Rho2d_pert.rsf 
rsrsf2rss <Rho2d_pert.rsf > Rho2d_pert.rss

mpirun -np 4 rsmpiElasticmod2d mod.cfg 
mpirun -np 4 rsmpiElasticfwigrad2d fwi.cfg 
