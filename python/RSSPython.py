import numpy as np
import struct as st
import matplotlib.pyplot as plt

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
                    self.data[j*self.geomN[1] + i] = st.unpack('f', f.read(4))[0];
        else:
            for i in range(0,self.fullsize):
                self.data[i] = st.unpack('f', f.read(4))[0];

        # Reshape data to correct dimensions
        self.data=self.data.reshape(resize)

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
        dataout=self.data.reshape([self.fullsize,1])
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

