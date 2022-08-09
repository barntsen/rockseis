int Ac2dFwstepvx(float [*,*] Pleft, float [*,*] Pright,   
              float [*,*] Vx,float [*,*] Rx, float [*,*] df,    
              float [*] Alftstag,float [*] Blftstag,float [*] Clftstag, 
              float [*] Arbbstag,float [*] Brbbstag,float [*] Crbbstag, 
              int [*] Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt){} 
              
int Ac2dFwstepvz(float [*,*] Ptop, float [*,*] Pbot,   
              float [*,*] Vz,float [*,*] Rz, float [*,*] df,    
              float [*] Alftstag,float [*] Blftstag,float [*] Clftstag, 
              float [*] Arbbstag,float [*] Brbbstag,float [*] Crbbstag, 
              int [*] Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt){} 


int Ac2dFwstepstressx(float [*,*] Vxxleft, float [*,*] Vxxright,   
              float [*,*] P,float [*,*] L, float [*,*] df,    
              float [*] Alft,float [*] Blft,float [*] Clft, 
              float [*] Arbb,float [*] Brbb,float [*] Crbb, 
              int [*] Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt){} 
              
int Ac2dFwstepstressz(float [*,*] Vzztop, float [*,*] Vzzbot,   
              float [*,*] P,float [*,*] L, float [*,*] df,    
              float [*] Alft,float [*] Blft,float [*] Clft, 
              float [*] Arbb,float [*] Brbb,float [*] Crbb, 
              int [*] Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt){} 
              
