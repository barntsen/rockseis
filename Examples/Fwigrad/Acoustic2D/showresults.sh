cd Results
rsrss2rsf <Pres.rss out=stdout > Pres.rsf 
rsrss2rsf <Pmod.rss out=stdout > Pmod.rsf 
rsrss2rsf <Pshot.rss out=stdout > Prec.rsf 
rsrss2rsf <Vpgrad2d.rss out=stdout > Vpgrad2d.rsf 
rsrss2rsf <Rhograd2d.rss out=stdout > Rhograd2d.rsf 

sfgraph < Pres.rsf title=residual | sfpen &
sfgraph < Pmod.rsf title=modelled | sfpen &
sfgraph < Prec.rsf title=recorded | sfpen &
sfwindow < Vpgrad2d.rsf | sftransp | sfgrey color=j mean=y | sfpen &
sfwindow < Rhograd2d.rsf | sftransp | sfgrey color=j mean=y | sfpen &
cd ..
