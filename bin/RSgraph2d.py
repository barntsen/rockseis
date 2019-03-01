#!/usr/bin/env python
import RSSPython as rs
import argparse
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Plot a 2d rss model file .')
parser.add_argument("-i", dest="inputfile", required=True, help="input RSS file", metavar="FILE")
parser.add_argument("-o", dest="figfile", required=False, help="output image file", metavar="FILE")
parser.add_argument("--format", dest="format", required=False, help="output image file format: pdf (default), png, eps.", default="pdf")
parser.add_argument("--aspect", dest="aspect", type=float, required=False, help="set aspect ratio (default = 1).", default=1)

args = parser.parse_args()
model = rs.RSSdata()
model.read(args.inputfile)

n = model.geomN[0];
d = model.geomD[0];
o = model.geomO[0];

if(model.geomN[1] > n):
    n = model.geomN[1]
    d = model.geomD[1];
    o = model.geomO[1];


min_x = float(o)
max_x = float(min_x + (n-1.)*d)
min_t = np.amin(model.data)
max_t = np.amax(model.data)
if(min_t == max_t):
    min_t = min_t*0.9;
    max_t = max_t*1.1;
    
t = np.arange(o, (n*d + o), d);
#fig = plt.plot(t, (model.data).squeeze(), '+')
fig = plt.plot(t, (model.data).squeeze())

delta = max_x - min_x
if(delta == 0):
    delta = 1

data_r = (max_t-min_t)/(delta)
disp_r = args.aspect;
plt.axes().set_aspect(disp_r/data_r)

plt.xlabel('Trace')
plt.ylabel('Time (s)')
plt.rc('font', size=15)
if(args.figfile is not None):
    plt.savefig(args.figfile,bbox_inches='tight', format=args.format)
else:
    plt.show()

