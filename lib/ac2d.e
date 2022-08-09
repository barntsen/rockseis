include "ac2d.i"



int Ac2dFwstepvx(float [*,*] Pleft, float [*,*] Pright,   
              float [*,*] Vx,float [*,*] Rx, float [*,*] df,    
              float [*] Alftstag,float [*] Blftstag,float [*] Clftstag, 
              float [*] Arbbstag,float [*] Brbbstag,float [*] Crbbstag, 
              int [*] Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt) 
              
{
  int nx, nz;
  int ix,iz,i;
  int lpml;

  nx   = len(Vx,0);
  nz   = len(Vx,1);
  lpml = len(Alftstag,0);

  // Compute Vx
  parallel(ix=0:nx,iz=0:nz){
    Vx[ix,iz] = Vx[ix,iz] + dt*Rx[ix,iz]*df[ix,iz];
  }

  // Attenuate left and right using staggered variables
  if(Getapplypml[0] || Getapplypml[1]){

    parallel(ix=0:lpml,iz=0:nz){
      if(Getapplypml[0]){
         if(ix >= ix0 && ix < (ix0 + nx)){
          // Left
          Pleft[ix,iz]  = Blftstag[ix]*Pleft[ix,iz] 
                        + Alftstag[ix]*df[ix-ix0,iz];
          Vx[ix-ix0,iz] = Vx[ix-iz0,iz] 
                        - dt*Rx[ix-ix0,iz]*(Pleft[ix,iz] 
                        + Clftstag[ix]*df[ix-ix0,iz]);
        }
      }
      if(Getapplypml[1]){
        i = ix + nxo - lpml;
        if(i >= ix0 && i < (ix0 + nx)){
          // Right
          Pright[ix,iz]  = Brbbstag[ix]*Pright[ix,iz] 
                         + Arbbstag[ix]*df[i-ix0,iz];
          Vx[i-ix0,iz]   = Vx[ix-ix0,iz] 
                         - dt*Rx[i-ix0,iz]*(Pright[ix,iz] 
                         + Crbbstag[ix]*df[i-ix0,iz]);
        }
      }
    }
  }
}

int Ac2dFwstepvz(float [*,*] Ptop, float [*,*] Pbot,   
              float [*,*] Vz,float [*,*] Rz, float [*,*] df,    
              float [*] Alftstag,float [*] Blftstag,float [*] Clftstag, 
              float [*] Arbbstag,float [*] Brbbstag,float [*] Crbbstag, 
              int [*] Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt) 
{
  int nx, nz;
  int ix,iz,i;
  int lpml;

  nx = len(Vz,0);
  nz = len(Vz,1);
  lpml = len(Alftstag,0);
  // Compute Vz
  parallel(ix=0:nx,iz=0:nz){
    Vz[ix,iz] = Vz[ix,iz] + dt*Rz[ix,iz]*df[ix,iz];
  }

  // Attenuate bottom and top using staggered variables
  if(Getapplypml[4] || Getapplypml[5]){
    parallel(iz=0:lpml,ix=0:nx){
      if(Getapplypml[4]){
        if((iz >= iz0) && (iz < (iz0 + nz))){
          // Top
          Ptop[ix,iz] = Blftstag[iz]*Ptop[ix,iz] 
                      + Alftstag[iz]*df[ix,iz-iz0];
          Vz[ix,iz-iz0] = Vz[ix,iz-iz0]
                        -dt*Rz[ix,iz-iz0]*(Ptop[ix,iz] 
                        + Clftstag[iz]*df[ix,iz-iz0]);
                        
        }
      }
      if(Getapplypml[5]){
        i = iz + nzo - lpml;
        if(i >= iz0 && i < (iz0 + nz)){
          //Bottom
          Pbot[ix,iz] = Brbbstag[iz]*Pbot[ix,iz] 
                         + Arbbstag[iz]*df[ix,i-iz0];
          Vz[ix,i-iz0]   = Vz[ix,i-iz0]-dt*Rz[ix,i-iz0]*(Pbot[ix,iz] 
                         + Crbbstag[iz]*df[ix,i-iz0]);
        }
      }
    }
  }
}

int Ac2dFwstepstressx(float [*,*] Vxxleft, float [*,*] Vxxright,   
              float [*,*] P,float [*,*] L, float [*,*] df,    
              float [*] Alft,float [*] Blft,float [*] Clft, 
              float [*] Arbb,float [*] Brbb,float [*] Crbb, 
              int [*] Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt) 
              
{
  int nx, nz;
  int ix,iz,i;
  int lpml;

  nx   = len(P,0);
  nz   = len(P,1);
  lpml = len(Alft,0);

  // Compute Vx
  parallel(ix=0:nx,iz=0:nz){
    P[ix,iz] = P[ix,iz] + dt*L[ix,iz]*df[ix,iz];
  }

  // Attenuate left and right using staggered variables
  if(Getapplypml[0] || Getapplypml[1]){

    parallel(ix=0:lpml,iz=0:nz){
      if(Getapplypml[0]){
         if(ix >= ix0 && ix < (ix0 + nx)){
          // Left
          Vxxleft[ix,iz]  = Blft[ix]*Vxxleft[ix,iz] 
                        + Alft[ix]*df[ix-ix0,iz];
          P[ix-ix0,iz] = P[ix-iz0,iz] 
                        - dt*L[ix-ix0,iz]*(Vxxleft[ix,iz] 
                        + Clft[ix]*df[ix-ix0,iz]);
        }
      }
      if(Getapplypml[1]){
        i = ix + nxo - lpml;
        if(i >= ix0 && i < (ix0 + nx)){
          // Right
          Vxxright[ix,iz]  = Brbb[ix]*Vxxright[ix,iz] 
                         + Arbb[ix]*df[i-ix0,iz];
          P[i-ix0,iz]   = P[ix-ix0,iz] 
                         - dt*L[i-ix0,iz]*(Vxxright[ix,iz] 
                         + Crbb[ix]*df[i-ix0,iz]);
        }
      }
    }
  }
}


int Ac2dFwstepstressz(float [*,*] Vzztop, float [*,*] Vzzbot,   
              float [*,*] P,float [*,*] L, float [*,*] df,    
              float [*] Alft,float [*] Blft,float [*] Clft, 
              float [*] Arbb,float [*] Brbb,float [*] Crbb, 
              int [*] Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt) 
              
{
  int nx, nz;
  int ix,iz,i;
  int lpml;

  nx   = len(P,0);
  nz   = len(P,1);
  lpml = len(Alft,0);

  // Compute Vx
  parallel(ix=0:nx,iz=0:nz){
    P[ix,iz] = P[ix,iz] + dt*L[ix,iz]*df[ix,iz];
  }

  // Attenuate left and right using staggered variables
  if(Getapplypml[4] || Getapplypml[5]){

    parallel(ix=0:nx,iz=0:lpml){
      if(Getapplypml[4]){
         if(iz >= iz0 && iz < (iz0 + nz)){
          // Top
          Vzztop[ix,iz]  = Blft[iz]*Vzztop[ix,iz] 
                        + Alft[iz]*df[ix,iz-iz0];
          P[ix,iz-iz0] = P[ix,iz-iz0] 
                        - dt*L[ix,i-iz0]*(Vzztop[ix,iz] 
                        + Clft[iz]*df[ix,iz-iz0]);
        }
      }
      if(Getapplypml[5]){
        i = iz + nzo - lpml;
        if(i >= iz0 && i < (iz0 + nz)){
          // Bottom
          Vzzbot[ix,iz]  = Brbb[iz]*Vzzbot[ix,iz] 
                         + Arbb[iz]*df[ix,i-iz0];
          P[ix,i-iz0]   = P[ix,i-iz0] 
                         - dt*L[ix,i-iz0]*(Vzzbot[ix,iz] 
                         + Crbb[iz]*df[ix,i-iz0]);
        }
      }
    }
  }
}
