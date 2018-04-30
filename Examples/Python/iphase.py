import RSSPython as rs
from scipy.signal import hilbert
import numpy as np

indata  = rs.RSSdata();
indata.read('Vxshot.rss')
outdata = indata;
ntraces = indata.geomN[1] # Number of traces
fs = 1.0/indata.geomD[0] # Sampling frequency
for i in range(0, ntraces):
    trc=indata.data[:,i];
    z= hilbert(trc) #form the analytical signal
    inst_amplitude = np.abs(z) #envelope extraction
    inst_phase = np.unwrap(np.angle(z))#inst phase
    inst_freq = np.diff(inst_phase)/(2*np.pi)*fs #inst frequency
    outdata.data[:,i] = inst_phase

outdata.write('Vxphase.rss')
