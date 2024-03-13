sh makesurvey.sh 
# Make perturbed model 
rsrss2rsf <../../Models/Vp2d.rss out=stdout > Vp2d.rsf
sfspike <Vp2d.rsf k3=25 l3=26 k1=30 l1=35 mag=600. | sfsmooth rect1=5 rect3=5 out=stdout > pert.rsf
sfmath <Vp2d.rsf x=pert.rsf output="in+x" out=stdout > Vp2d_pert.rsf 
rsrsf2rss <Vp2d_pert.rsf > Vp2d_pert.rss

rsrss2rsf <../../Models/Vs2d.rss out=stdout > Vs2d.rsf
sfspike <Vs2d.rsf k3=25 l3=26 k1=10 l1=15 mag=-400. | sfsmooth rect1=5 rect3=5 out=stdout > pert.rsf
sfmath <Vs2d.rsf x=pert.rsf output="in+x" out=stdout > Vs2d_pert.rsf 
rsrsf2rss <Vs2d_pert.rsf > Vs2d_pert.rss

rsrss2rsf <../../Models/Rho2d.rss out=stdout > Rho2d.rsf
sfspike <Rho2d.rsf k3=25 l3=26 k1=20 l1=25 mag=00. | sfsmooth rect1=5 rect3=5 out=stdout > pert.rsf
sfmath <Rho2d.rsf x=pert.rsf output="in+x" out=stdout > Rho2d_pert.rsf 
rsrsf2rss <Rho2d_pert.rsf > Rho2d_pert.rss

sfspike <Vp2d.rsf nsp=3 k3=1,6,36 l3=5,35,45 mag=0,1,0 out=stdout > mute.rsf
rsrsf2rss <mute.rsf > mute.rss
rm *.rsf

cp ../../Models/Vp2d.rss .
cp ../../Models/Vs2d.rss .
cp ../../Models/Rho2d.rss .
cp ../../Models/Wav2d.rss .

mpirun -np 4 rsmpiElasticmod2d mod.cfg 

#Create a weight file
rsrss2rsf < Vxshot.rss out=stdout > temp.rsf
sfmath <temp.rsf output=1.0  | sfpow1d tpow=2.0 out=stdout > temp2.rsf 
sfheadermutter <temp2.rsf head=tfile.rsf v0=100 delay=0.0 type=0 out=stdout > weight.rsf
sfsegywrite < weight.rsf tfile=tfile.rsf tape=temp.sgy
rssegy2rss <temp.sgy > weight.rss segy.cfg
rm weight.rsf temp.rsf temp2.rsf temp.sgy tfile.rsf

# Add noise
#rsrss2rsf < Vxshot.rss out=stdout > temp.rsf
#sfnoise <temp.rsf range=900 out=stdout > noisy.rsf
#sfsegywrite < noisy.rsf tfile=tfile.rsf tape=temp.sgy
#rssegy2rss <temp.sgy > Vxshot.rss segy.cfg

#rsrss2rsf < Vzshot.rss out=stdout > temp.rsf
#sfnoise <temp.rsf range=900 out=stdout > noisy.rsf
#sfsegywrite < noisy.rsf tfile=tfile.rsf tape=temp.sgy
#rssegy2rss <temp.sgy > Vzshot.rss segy.cfg
#rm temp.rsf temp.sgy noisy.rsf tfile.rsf
