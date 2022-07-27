import numpy as np
import matplotlib.pyplot as plt
import RSSPython as rs

# Reading the data
rsfile = rs.RSSdata()
rsfile.read("Pshot.rss")

# Getting the information about the shot
nt = int(rsfile.geomN[0])
dt = np.float64(rsfile.geomD[0])
ot = np.float64(rsfile.geomO[0])
ntrc = int(rsfile.geomN[1])
twt = np.linspace(0, (nt-1)*dt, nt) + ot # two-way time

# You can also access source and receiver coordinates of the shot in the ndarrays on the fields:
# rsfile.srcX[]
# rsfile.srcZ[]
# rsfile.GroupX[]
# rsfile.GroupZ[]

# Displaying the shots using matplotlib
fig,axs = plt.subplots(1,3,figsize=(24, 8))

# plotting Pshot
ax = axs[0]
vm = np.percentile(rsfile.data, 99)
extent = [0, ntrc, twt[-1], twt[0]]  # define extent",
ax.imshow(rsfile.data, cmap="RdBu", vmin=-vm, vmax=vm, aspect='auto', extent=extent)
ax.set_xlabel('TRC number')
ax.set_ylabel('TWT [s]')
ax.set_title('Pshot')

# plotting Vxshot
rsfile.read("Vxshot.rss") # Read the data
ax = axs[1]
vm = np.percentile(rsfile.data, 99)
extent = [0, ntrc, twt[-1], twt[0]]  # define extent",
ax.imshow(rsfile.data, cmap="RdBu", vmin=-vm, vmax=vm, aspect='auto', extent=extent)
ax.set_xlabel('TRC number')
ax.set_ylabel('TWT [s]')
ax.set_title('Vxshot')

# plotting Vzshot
rsfile.read("Vzshot.rss") # Read the data
ax = axs[2]
vm = np.percentile(rsfile.data, 99)
extent = [0, ntrc, twt[-1], twt[0]]  # define extent",
ax.imshow(rsfile.data, cmap="RdBu", vmin=-vm, vmax=vm, aspect='auto', extent=extent)
ax.set_xlabel('TRC number')
ax.set_ylabel('TWT [s]')
ax.set_title('Vzshot')

plt.show()

