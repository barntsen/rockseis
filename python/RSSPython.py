import numpy as np
import struct as st
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

MAGICNUMBER='R0CKS' #rss identifier
MAGICNUMBERLENGTH=5
MAXDIMS=9

class RSSdata:
    def __init__(self):
        self.data_format=0;
        self.header_format=0;
        self.type=0;
        self.Nheader=0;
        self.Ndims=0;
        self.fullsize=0
        self.geomN = np.zeros([MAXDIMS,1], dtype='uint64')
        self.geomD = np.zeros([MAXDIMS,1], dtype='float64')
        self.geomO = np.zeros([MAXDIMS,1], dtype='float64')

    def read(self,filename):
        f = open(filename, 'rb')
        buff= f.read(MAGICNUMBERLENGTH)
    
        # Check if file has rss identifier
        if(buff != MAGICNUMBER):
            print 'Not a valid rss file!'
            return -1
    
        # Reading header
        head = {}
        self.data_format = st.unpack('i', f.read(4))[0]
        self.header_format = st.unpack('i', f.read(4))[0]
        self.type = st.unpack('i', f.read(4))[0]
        self.Nheader = st.unpack('Q', f.read(8))[0]
        self.Ndims = st.unpack('Q', f.read(8))[0]
        for i in range(0,MAXDIMS):
            self.geomN[i] = st.unpack('Q', f.read(8))[0]
        for i in range(0,MAXDIMS):
            self.geomD[i] = st.unpack('d', f.read(8))[0]
        for i in range(0,MAXDIMS):
            self.geomO[i] = st.unpack('d', f.read(8))[0]

        #Reading coordinates and data
        self.fullsize = 1;
        self.Ndims = 0
        for i in range(0, MAXDIMS):
            if(self.geomN[i] > 0):    
                self.fullsize = self.fullsize*self.geomN[i];
                self.Ndims = self.Ndims + 1
        resize = np.zeros([self.Ndims,1], dtype='uint64')
        for i in range(0,self.Ndims):
            resize[i] = self.geomN[i];

        # Allocate 1d array for data
        self.data = np.zeros([self.fullsize,1], dtype='float32')
        if(self.Nheader):
            self.srcX = np.zeros([self.geomN[1],1], dtype='float32')
            self.srcY = np.zeros([self.geomN[1],1], dtype='float32')
            self.srcZ = np.zeros([self.geomN[1],1], dtype='float32')
            self.GroupX = np.zeros([self.geomN[1],1], dtype='float32')
            self.GroupY = np.zeros([self.geomN[1],1], dtype='float32')
            self.GroupZ = np.zeros([self.geomN[1],1], dtype='float32')
            for i in range(0,self.geomN[1]):
                if(self.Nheader == 4):
                    self.srcX[i] = st.unpack('f', f.read(4))[0];
                    self.srcZ[i] = st.unpack('f', f.read(4))[0];
                    self.GroupX[i] = st.unpack('f', f.read(4))[0];
                    self.GroupZ[i] = st.unpack('f', f.read(4))[0];
        
                if(self.Nheader == 6):
                    self.srcX[i] = st.unpack('f', f.read(4))[0];
                    self.srcY[i] = st.unpack('f', f.read(4))[0];
                    self.srcZ[i] = st.unpack('f', f.read(4))[0];
                    self.GroupX[i] = st.unpack('f', f.read(4))[0];
                    self.GroupY[i] = st.unpack('f', f.read(4))[0];
                    self.GroupZ[i] = st.unpack('f', f.read(4))[0];
                            
                for j in range(0,self.geomN[0]):
                    self.data[i*self.geomN[0] + j] = st.unpack('f', f.read(4))[0];
        else:
            for i in range(0,self.fullsize):
               self.data[i] = st.unpack('f', f.read(4))[0];

        # Reshape data to correct dimensions
        self.data=self.data.reshape(resize, order='F')

    def write(self,filename):
        f = open(filename, 'wb')
        f.write(MAGICNUMBER)
    
        # Writting header
        f.write(st.pack('i', self.data_format))
        f.write(st.pack('i', self.header_format))
        f.write(st.pack('i', self.type))
        f.write(st.pack('Q', self.Nheader))
        f.write(st.pack('Q', self.Ndims))
        for i in range(0,MAXDIMS):
            f.write(st.pack('Q', self.geomN[i]))
        for i in range(0,MAXDIMS):
            f.write(st.pack('d', self.geomD[i]))
        for i in range(0,MAXDIMS):
            f.write(st.pack('d', self.geomO[i]))

        #Estimating sizes for output
        self.fullsize = 1;
        self.Ndims = 0
        for i in range(0, MAXDIMS):
            if(self.geomN[i] > 0):    
                self.fullsize = self.fullsize*self.geomN[i];
                self.Ndims = self.Ndims + 1

        #Writting Coordinates and data
        dataout=self.data.reshape([self.fullsize,1], order='F')
        if(self.Nheader):
            for i in range(0,self.geomN[1]):
                if(self.Nheader == 4):
                    f.write(st.pack('f', self.srcX[i]))
                    f.write(st.pack('f', self.srcZ[i]))
                    f.write(st.pack('f', self.GroupX[i]))
                    f.write(st.pack('f', self.GroupZ[i]))
        
                if(self.Nheader == 6):
                    f.write(st.pack('f', self.srcX[i]))
                    f.write(st.pack('f', self.srcY[i]))
                    f.write(st.pack('f', self.srcZ[i]))
                    f.write(st.pack('f', self.GroupX[i]))
                    f.write(st.pack('f', self.GroupY[i]))
                    f.write(st.pack('f', self.GroupZ[i]))
                            
                for j in range(0,self.geomN[0]):
                    f.write(st.pack('f', dataout[j*self.geomN[1] + i]))
        else:
            for i in range(0,self.fullsize):
               f.write(st.pack('f',  dataout[i]))
        
    def imageshow(self):
        nx = self.geomN[0];
        dx = self.geomD[0];
        ox = self.geomO[0];

        nz = self.geomN[2];
        dz = self.geomD[2];
        oz = self.geomO[2];

        min_x = float(ox)
        max_x = float(min_x + (nx-1.)*dx)
        min_z = float(oz)
        max_z = float(min_z + (nz-1.)*dz)

        #plt.imshow(self.data.squeeze(), extent=[min_x,max_x,max_z,min_z])
        plt.imshow(self.data.squeeze())
        plt.colorbar()
        plt.xlabel('x (m)')
        plt.ylabel('z (m)')
        plt.rc('font', size=15)
        plt.show()

def uiscolor():
    colors = [(0.99609375, 0.5390625, 0.99609375), (0.99609375, 0.52734375, 0.99609375), (0.99609375, 0.51953125, 0.99609375), (0.99609375, 0.51171875, 0.99609375), (0.99609375, 0.50390625, 0.99609375), (0.99609375, 0.49609375, 0.99609375), (0.99609375, 0.48828125, 0.99609375), (0.99609375, 0.48046875, 0.99609375), (0.99609375, 0.47265625, 0.99609375), (0.99609375, 0.4609375, 0.99609375), (0.99609375, 0.453125, 0.99609375), (0.99609375, 0.4453125, 0.99609375), (0.99609375, 0.4375, 0.99609375), (0.99609375, 0.4296875, 0.99609375), (0.98828125, 0.421875, 0.99609375), (0.98046875, 0.41015625, 0.99609375), (0.97265625, 0.40234375, 0.99609375), (0.95703125, 0.39453125, 0.99609375), (0.93359375, 0.38671875, 0.99609375), (0.90625, 0.375, 0.99609375), (0.87890625, 0.3671875, 0.99609375), (0.85546875, 0.359375, 0.99609375), (0.828125, 0.34765625, 0.99609375), (0.8046875, 0.33984375, 0.99609375), (0.78125, 0.33203125, 0.99609375), (0.7578125, 0.3203125, 0.99609375), (0.734375, 0.3125, 0.99609375), (0.70703125, 0.30078125, 0.99609375), (0.68359375, 0.29296875, 0.99609375), (0.66015625, 0.28125, 0.99609375), (0.63671875, 0.2734375, 0.99609375), (0.61328125, 0.26171875, 0.99609375), (0.5859375, 0.25390625, 0.99609375), (0.5625, 0.2421875, 0.99609375), (0.5390625, 0.23046875, 0.99609375), (0.515625, 0.22265625, 0.99609375), (0.48828125, 0.2109375, 0.99609375), (0.46484375, 0.19921875, 0.99609375), (0.4375, 0.1875, 0.99609375), (0.4140625, 0.17578125, 0.99609375), (0.38671875, 0.1640625, 0.99609375), (0.36328125, 0.15234375, 0.99609375), (0.3359375, 0.140625, 0.99609375), (0.30859375, 0.125, 0.99609375), (0.28125, 0.109375, 0.99609375), (0.24609375, 0.09375, 0.99609375), (0.20703125, 0.0703125, 0.99609375), (0.16796875, 0.04296875, 0.99609375), (0.12890625, 0.0234375, 0.99609375), (0.09375, 0.0234375, 0.99609375), (0.05859375, 0.03515625, 0.9921875), (0.02734375, 0.05859375, 0.98828125), (0.00390625, 0.078125, 0.984375), (0.00390625, 0.08984375, 0.984375), (0.02734375, 0.09765625, 0.98828125), (0.06640625, 0.10546875, 0.9921875), (0.10546875, 0.1171875, 0.99609375), (0.12890625, 0.12890625, 0.99609375), (0.15234375, 0.1484375, 0.99609375), (0.16796875, 0.17578125, 0.99609375), (0.18359375, 0.203125, 0.99609375), (0.1953125, 0.2265625, 0.99609375), (0.203125, 0.25, 0.99609375), (0.2109375, 0.26953125, 0.99609375), (0.21484375, 0.29296875, 0.99609375), (0.21875, 0.3125, 0.99609375), (0.21875, 0.33203125, 0.99609375), (0.22265625, 0.35546875, 0.99609375), (0.22265625, 0.375, 0.99609375), (0.22265625, 0.39453125, 0.99609375), (0.22265625, 0.4140625, 0.99609375), (0.22265625, 0.43359375, 0.99609375), (0.21875, 0.453125, 0.99609375), (0.21484375, 0.46875, 0.99609375), (0.21484375, 0.48828125, 0.99609375), (0.2109375, 0.5078125, 0.99609375), (0.20703125, 0.52734375, 0.99609375), (0.203125, 0.54296875, 0.99609375), (0.19921875, 0.5625, 0.99609375), (0.1953125, 0.578125, 0.99609375), (0.1875, 0.59765625, 0.99609375), (0.18359375, 0.61328125, 0.99609375), (0.1796875, 0.6328125, 0.99609375), (0.171875, 0.6484375, 0.99609375), (0.16796875, 0.66796875, 0.99609375), (0.16015625, 0.68359375, 0.99609375), (0.15625, 0.69921875, 0.99609375), (0.1484375, 0.71875, 0.99609375), (0.140625, 0.734375, 0.99609375), (0.1328125, 0.75, 0.99609375), (0.125, 0.765625, 0.99609375), (0.1171875, 0.78515625, 0.99609375), (0.109375, 0.80078125, 0.99609375), (0.09765625, 0.81640625, 0.99609375), (0.0859375, 0.83203125, 0.99609375), (0.07421875, 0.8515625, 0.99609375), (0.05859375, 0.8671875, 0.99609375), (0.04296875, 0.8828125, 0.99609375), (0.03125, 0.90234375, 0.99609375), (0.015625, 0.921875, 0.99609375), (0.00390625, 0.9375, 0.9921875), (0.0, 0.94921875, 0.9921875), (0.0, 0.95703125, 0.984375), (0.0, 0.9609375, 0.97265625), (0.0, 0.96484375, 0.95703125), (0.0, 0.96875, 0.93359375), (0.0, 0.96484375, 0.9140625), (0.0, 0.95703125, 0.88671875), (0.0, 0.9453125, 0.859375), (0.0, 0.93359375, 0.828125), (0.0, 0.92578125, 0.80078125), (0.0, 0.9140625, 0.77734375), (0.0, 0.90625, 0.75), (0.0, 0.89453125, 0.72265625), (0.0, 0.88671875, 0.69921875), (0.0, 0.875, 0.671875), (0.0, 0.8671875, 0.6484375), (0.0, 0.85546875, 0.625), (0.0, 0.84375, 0.6015625), (0.0, 0.8359375, 0.578125), (0.0, 0.82421875, 0.5546875), (0.0, 0.81640625, 0.53125), (0.0, 0.8046875, 0.5078125), (0.0, 0.79296875, 0.48828125), (0.0, 0.78515625, 0.46875), (0.0, 0.7734375, 0.44921875), (0.0, 0.76171875, 0.4296875), (0.0, 0.75390625, 0.41015625), (0.0, 0.7421875, 0.39453125), (0.0, 0.73046875, 0.37890625), (0.0, 0.72265625, 0.36328125), (0.0, 0.7109375, 0.3515625), (0.0, 0.69921875, 0.33984375), (0.0, 0.69140625, 0.328125), (0.0, 0.6796875, 0.31640625), (0.0, 0.66796875, 0.30859375), (0.0, 0.66015625, 0.30078125), (0.0, 0.6484375, 0.296875), (0.0, 0.63671875, 0.2890625), (0.0, 0.625, 0.28515625), (0.0, 0.6171875, 0.28125), (0.0, 0.60546875, 0.27734375), (0.0, 0.59375, 0.2734375), (0.0, 0.58203125, 0.265625), (0.0, 0.5703125, 0.2578125), (0.0, 0.55859375, 0.25390625), (0.0, 0.546875, 0.2421875), (0.0, 0.53515625, 0.234375), (0.0, 0.51953125, 0.2265625), (0.0, 0.5078125, 0.21484375), (0.0, 0.49609375, 0.203125), (0.0, 0.4765625, 0.18359375), (0.0, 0.4609375, 0.16796875), (0.0, 0.44921875, 0.15234375), (0.0, 0.4453125, 0.1484375), (0.0, 0.4453125, 0.15234375), (0.0, 0.453125, 0.15625), (0.0, 0.4609375, 0.16796875), (0.0, 0.47265625, 0.17578125), (0.0, 0.484375, 0.18359375), (0.0, 0.49609375, 0.1953125), (0.0, 0.5078125, 0.203125), (0.0, 0.5234375, 0.2109375), (0.0, 0.53515625, 0.21875), (0.0, 0.546875, 0.22265625), (0.0, 0.5546875, 0.2265625), (0.0, 0.56640625, 0.2265625), (0.0, 0.578125, 0.2265625), (0.0, 0.5859375, 0.2265625), (0.0, 0.59765625, 0.22265625), (0.0, 0.60546875, 0.21875), (0.0, 0.61328125, 0.21484375), (0.0, 0.62109375, 0.203125), (0.0, 0.62890625, 0.19140625), (0.0, 0.63671875, 0.17578125), (0.00390625, 0.640625, 0.140625), (0.01171875, 0.6484375, 0.08203125), (0.01953125, 0.65234375, 0.02734375), (0.03515625, 0.66015625, 0.0), (0.06640625, 0.66796875, 0.0), (0.12109375, 0.6796875, 0.0), (0.18359375, 0.69140625, 0.0), (0.23828125, 0.703125, 0.0), (0.2734375, 0.71484375, 0.0), (0.3046875, 0.72265625, 0.0), (0.3359375, 0.734375, 0.0), (0.36328125, 0.74609375, 0.0), (0.39453125, 0.75390625, 0.0), (0.421875, 0.765625, 0.0), (0.44921875, 0.7734375, 0.0), (0.4765625, 0.78515625, 0.0), (0.50390625, 0.79296875, 0.0), (0.53125, 0.8046875, 0.0), (0.55859375, 0.8125, 0.0), (0.5859375, 0.8203125, 0.0), (0.61328125, 0.83203125, 0.0), (0.640625, 0.83984375, 0.0), (0.6640625, 0.84765625, 0.0), (0.69140625, 0.859375, 0.0), (0.71875, 0.8671875, 0.0), (0.74609375, 0.875, 0.0), (0.7734375, 0.8828125, 0.0), (0.80078125, 0.890625, 0.0), (0.828125, 0.8984375, 0.0), (0.86328125, 0.91015625, 0.0), (0.8984375, 0.921875, 0.0), (0.91796875, 0.92578125, 0.0), (0.921875, 0.92578125, 0.0), (0.921875, 0.91015625, 0.0), (0.921875, 0.88671875, 0.0), (0.921875, 0.86328125, 0.0), (0.92578125, 0.84375, 0.0), (0.92578125, 0.828125, 0.0), (0.92578125, 0.80859375, 0.0), (0.92578125, 0.7890625, 0.0), (0.92578125, 0.7734375, 0.0), (0.92578125, 0.75390625, 0.0), (0.92578125, 0.734375, 0.0), (0.92578125, 0.71484375, 0.0), (0.92578125, 0.69921875, 0.0), (0.92578125, 0.6796875, 0.0), (0.92578125, 0.66015625, 0.0), (0.92578125, 0.640625, 0.0), (0.92578125, 0.62109375, 0.0), (0.92578125, 0.6015625, 0.0), (0.92578125, 0.58203125, 0.0), (0.92578125, 0.5625, 0.0), (0.921875, 0.54296875, 0.0), (0.921875, 0.5234375, 0.0), (0.91796875, 0.50390625, 0.0), (0.91796875, 0.484375, 0.0), (0.9140625, 0.46484375, 0.0), (0.9140625, 0.44140625, 0.0), (0.91015625, 0.421875, 0.0), (0.90625, 0.40234375, 0.0), (0.90234375, 0.37890625, 0.0), (0.8984375, 0.35546875, 0.0), (0.89453125, 0.33203125, 0.0), (0.890625, 0.3125, 0.0), (0.88671875, 0.28515625, 0.0), (0.87890625, 0.26171875, 0.0), (0.875, 0.23828125, 0.0), (0.8671875, 0.20703125, 0.0), (0.859375, 0.17578125, 0.0), (0.8515625, 0.125, 0.0), (0.84375, 0.06640625, 0.00390625), (0.8359375, 0.01953125, 0.01171875), (0.828125, 0.0, 0.015625), (0.8203125, 0.0, 0.0390625), (0.8125, 0.0, 0.08203125), (0.8046875, 0.0, 0.12109375), (0.796875, 0.0, 0.14453125), (0.7890625, 0.0, 0.15625), (0.78125, 0.0, 0.16796875), (0.76953125, 0.0, 0.171875), (0.7578125, 0.0, 0.17578125)]
    cm = LinearSegmentedColormap.from_list('uiscolor', colors, N=256)
    return cm
