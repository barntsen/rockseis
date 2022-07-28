import numpy as np
import matplotlib.pyplot as plt
import RSSPython as rs

# Reading the data
rsfile = rs.RSSdata()
rsfile.read("Psnaps.rss")

# Getting the information about the snap
nx = int(rsfile.geomN[0])
dx = np.float64(rsfile.geomD[0])
ox = np.float64(rsfile.geomO[0])

nz = int(rsfile.geomN[2])
dz = np.float64(rsfile.geomD[2])
oz = np.float64(rsfile.geomO[2])

nt = int(rsfile.geomN[3])
dt = np.float64(rsfile.geomD[3])
ot = np.float64(rsfile.geomO[3])

x = np.linspace(0, (nx-1)*dx, nx) + ox 
z = np.linspace(0, (nz-1)*dz, nz) + oz 
t = np.linspace(0, (nt-1)*dt, nt) + ot 

# Displaying the snaps using matplotlib
fig,axs = plt.subplots(3,1,figsize=(12, 6))

snapnum = 5 

# plotting Psnap
snap = rsfile.data[:,:,:,snapnum].squeeze().T
ax = axs[0]
vm = np.percentile(snap, 99)
extent = [x[0], x[-1], z[-1], z[0]]  # define extent",
ax.imshow(snap, cmap="RdBu", vmin=-vm, vmax=vm, aspect=1, extent=extent)
ax.set_xlabel('X [m]')
ax.set_ylabel('Z [m]')
ax.set_title('Psnap at ' + str(t[snapnum]) + ' s')

# plotting Vxsnap
rsfile.read("Vxsnaps.rss")
snap = rsfile.data[:,:,:,snapnum].squeeze().T
ax = axs[1]
vm = np.percentile(snap, 99)
extent = [x[0], x[-1], z[-1], z[0]]  # define extent",
ax.imshow(snap, cmap="RdBu", vmin=-vm, vmax=vm, aspect=1, extent=extent)
ax.set_xlabel('X [m]')
ax.set_ylabel('Z [m]')
ax.set_title('Vxsnap')
ax.set_title('Vxsnap at ' + str(t[snapnum]) + ' s')

# plotting Vzsnap
rsfile.read("Vzsnaps.rss")
snap = rsfile.data[:,:,:,snapnum].squeeze().T
ax = axs[2]
vm = np.percentile(snap, 99)
extent = [x[0], x[-1], z[-1], z[0]]  # define extent",
ax.imshow(snap, cmap="RdBu", vmin=-vm, vmax=vm, aspect=1, extent=extent)
ax.set_xlabel('X [m]')
ax.set_ylabel('Z [m]')
ax.set_title('Vzsnap at ' + str(t[snapnum]) + ' s')

plt.show()

