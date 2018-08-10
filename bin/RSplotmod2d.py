#!/usr/bin/env python
import RSSPython as rs
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Plot a 2d rss model file .')
parser.add_argument("-i", dest="inputfile", required=True, help="input RSS file", metavar="FILE")
parser.add_argument("-cbar", type=int, dest="cbar", required=False, help="Add a colorbar", default=False)
parser.add_argument("-o", dest="figfile", required=False, help="output image file", metavar="FILE")
parser.add_argument("--format", dest="format", required=False, help="output image file format: pdf (default), png, eps.", default="pdf")
parser.add_argument("--interp", dest="interp", required=False, help="interpolation method: neareast (default), bicubic, lanczos.", default="nearest")
parser.add_argument("--showclip", dest="showclip", required=False, help="show current limits of the color scale.", default=False)
parser.add_argument("--vmin", dest="minclip", type=float, required=False, help="set minimum clip value for the color scale.", default=argparse.SUPPRESS)
parser.add_argument("--vmax", dest="maxclip", type=float, required=False, help="set maximum clip value for the color scale.", default=argparse.SUPPRESS)
parser.add_argument("--aspect", dest="aspect", type=float, required=False, help="set aspect ratio (default = 1).", default=1)
parser.add_argument("--cmap", dest="cmap", required=False, help="set colormap (default = 'gray').", default='jet')

args = parser.parse_args()
model = rs.RSSdata()
model.read(args.inputfile)

nx = model.geomN[0];
dx = model.geomD[0];
ox = model.geomO[0];

nz = model.geomN[2];
dz = model.geomD[2];
oz = model.geomO[2];

min_x = float(ox)
max_x = float(min_x + (nx-1.)*dx)
min_z = float(oz)
max_z = float(min_z + (nz-1.)*dz)

fig = plt.imshow((model.data.squeeze()).transpose(),interpolation=args.interp, extent=[min_x,max_x,max_z,min_z])
vmin, vmax = fig.get_clim()

fig.set_cmap(args.cmap)

if(hasattr(args, 'minclip')):
    vmin = args.minclip
if(hasattr(args, 'maxclip')):
    vmax = args.maxclip
fig.set_clim(vmin,vmax)

if(args.showclip):
    print "--vmin",vmin," --vmax",vmax

plt.axes().set_aspect(args.aspect)

if(args.cbar == True and not vmin == vmax):
    plt.colorbar(orientation='horizontal')
plt.xlabel('x (m)')
plt.ylabel('z (m)')
plt.rc('font', size=15)
if(args.figfile is not None):
    plt.savefig(args.figfile,bbox_inches='tight', format=args.format)
plt.show()

