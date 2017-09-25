#sfmakehdr3d dsx=10 nsx=5 dgx=10 ngx=10 scalco=-100 out=stdout > hdr.rsf 
#sfspike n1=501 n2=50 d2=1 d1=1e-3 o2=0 o1=0 out=stdout > traces.rsf 
#sfsegywrite <traces.rsf tfile=hdr.rsf tape=Survey.sgy 
#rm traces.rsf hdr.rsf
#../../bin/rssegy2rss <Survey.sgy >2DSurvey.rss
#../../bin/rsinfo <2DSurvey.rss
#rm Survey.sgy
#
rsmakesurvey survey.cfg >2DSurvey.rss
