import numpy as np
import struct as st

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
        fullsize = 1;
        self.Ndims = 0
        for i in range(0, MAXDIMS):
            if(self.geomN[i] > 0):    
                fullsize = fullsize*self.geomN[i];
                self.Ndims = self.Ndims + 1
        resize = np.zeros([self.Ndims,1], dtype='uint64')
        for i in range(0,self.Ndims):
            resize[i] = self.geomN[i];

        # Allocate 1d array for data
        self.data = np.zeros([fullsize,1], dtype='float32')
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
            for i in range(1,fullsize):
                self.data[i] = st.unpack('f', f.read(4))[0];

        # Reshape data to correct dimensions
        self.data=self.data.reshape(resize)
        print self.data.shape

