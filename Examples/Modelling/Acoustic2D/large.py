#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import RSSPython as rs


# In[19]:
## Some function definitions
def Ricker(f0, t0, Nt, dt):
    t=np.linspace(0,(Nt-1)*dt,Nt) - t0;
    s = -np.pi*np.pi*f0*f0
    S = np.zeros([Nt,1])
    S[:,0]=(1+2*s*t*t)*np.exp(s*t*t)
    return S,t
    

def fd2dstability(dx, dz, Cmax):
    dt=2/(np.pi*np.sqrt((1./(dx**2))+(1./(dz**2)))*Cmax)
    return dt

# In[24]:
## Creating a model 
nx = 1001;
nz = 1001;
dx = 10;
dz = 10;
z = np.zeros([nz,1])
z[:,0] = np.linspace(0,nz-1*dz, nz);
G=1.5;
G=0.0
vp = 2000 + G*z;
rho = 1000*np.ones(vp.shape);


### Expanding these models into 2D models in rss shape

vp = np.reshape(np.dot(vp,np.ones([1,nx])).T, [nx,1,nz]);
rho = np.reshape(np.dot(rho,np.ones([1,nx])).T, [nx,1,nz]);

# In[26]:

### Creating a rockseis model header for writting out the model
mod = rs.RSSdata(vp);
# Setting origin
mod.geomO[0] = 0.0;
mod.geomO[2] = 0.0;

# Setting sampling
mod.geomD[0] = dx;
mod.geomD[2] = dz;

### Writting out the models to disk
mod.data = vp
mod.write('vp.rss');
mod.data = rho
mod.write('rho.rss');

## Wavelet (source signature)
f0 = 15;
t0 = 2/f0;
dt = 0.9*fd2dstability(dx, dz, 4000);
rec_time = 2.0 
nt = int(np.floor(rec_time/dt)); 
wav,t = Ricker(f0, t0, nt, dt); 

wavfile = rs.RSSdata(wav,2); # 2 for 2D or 3 for 3D
wavfile.geomD[0] = dt;
wavfile.write('wav2d.rss')

## Survey (source and receiver positions)
nsources = 1;
nreceivers = 801;
data = np.ones([1, nsources*nreceivers]);
surveyfile = rs.RSSdata(data,2); # 2 for 2D or 3 for 3D
count=0;
for i in range(0,nsources):
    sx = dx*(nx-1)/2;
    sz = dz*2.0
    minoff = (nreceivers*dx)/2
    for j in range(0,nreceivers):
        rx = sx+(j)*dx-minoff;
        rz = 30.0;        
        surveyfile.srcX[count] = sx;
        surveyfile.srcZ[count] = sz;
        surveyfile.GroupX[count] = rx;
        surveyfile.GroupZ[count] = rz;
        count = count+1;

surveyfile.write('2DSurvey.rss');


