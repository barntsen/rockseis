rsrss2rsf <Vpgrad2d.rss out=stdout > Vpgrad2d.rsf
sfwindow <Vpgrad2d.rsf | sftransp | sfgrey title=Vp color=j mean=y | sfpen &

rsrss2rsf <Vsgrad2d.rss out=stdout > Vsgrad2d.rsf
sfwindow <Vsgrad2d.rsf | sftransp | sfgrey title=Vs color=j mean=y | sfpen &

rsrss2rsf <Rhograd2d.rss out=stdout > Rhograd2d.rsf
sfwindow <Rhograd2d.rsf | sftransp | sfgrey title=Rho color=j mean=y | sfpen &

#rsrss2rsf <Vxres2d.rss out=stdout > Vxres2d.rsf
#rsrss2rsf <Vxmod2d.rss out=stdout > Vxmod2d.rsf
#rsrss2rsf <Results/Vxshot.rss out=stdout > Vxshot2d.rsf

#rsrss2rsf <Vzres2d.rss out=stdout > Vzres2d.rsf
#rsrss2rsf <Vzmod2d.rss out=stdout > Vzmod2d.rsf
#rsrss2rsf <Results/Vzshot.rss out=stdout > Vzshot2d.rsf

#sfgrey <Vxshot2d.rsf | sfpen 
#sfgrey <Vxmod2d.rsf | sfpen 

#sfattr < Vxshot2d.rsf
#sfattr < Vxmod2d.rsf
