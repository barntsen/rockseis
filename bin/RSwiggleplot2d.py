#!/usr/bin/env python
import RSSPython as rs
import argparse
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Plot a 2d rss model file .')
parser.add_argument("-i", dest="inputfile", required=True, help="input RSS file", metavar="FILE")
parser.add_argument("-o", dest="figfile", required=False, help="output image file", metavar="FILE")
parser.add_argument("--format", dest="format", required=False, help="output image file format: pdf (default), png, eps.", default="pdf")
parser.add_argument("--showclip", dest="showclip", required=False, help="show current limits of the color scale.", default=False)
parser.add_argument("--vmin", dest="minclip", type=float, required=False, help="set minimum clip value for the color scale.", default=argparse.SUPPRESS)
parser.add_argument("--vmax", dest="maxclip", type=float, required=False, help="set maximum clip value for the color scale.", default=argparse.SUPPRESS)

args = parser.parse_args()
model = rs.RSSdata()
model.read(args.inputfile)

nt = model.geomN[0];
dt = model.geomD[0];
ot = model.geomO[0];

nx = model.geomN[1];
dx = model.geomD[1];
ox = model.geomO[1];

min_x = float(ox)
max_x = float(min_x + (nx-1.)*dx)
min_t = float(ot)
max_t = float(min_t + (nt-1.)*dt)
t = np.linspace(min_t, max_t, nt)
fig =  plt.figure() 
ax = fig.add_subplot(111)

# Normalize data
maxval = np.amax(abs(model.data))
if(not maxval == 0):
    model.data = model.data/maxval
for i in range(0,nx):
    trace = i + model.data[:,i]
    ax.plot(trace,t,'k-')
    ax.fill_betweenx(t,i,trace,where=(trace>i),color='k')

ax.set_xlim(-1.1,nx+0.1)
ax.set_ylim(min_t,max_t)
ax.invert_yaxis()
plt.xlabel('Trace')
plt.ylabel('Time (s)')
plt.rc('font', size=15)

if(args.figfile is not None):
    plt.savefig(args.figfile,bbox_inches='tight', format=args.format)
plt.show()

