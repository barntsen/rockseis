import numpy as np
import struct as st
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

MAGICNUMBER='R0CKS' #rss identifier
MAGICNUMBERLENGTH=5
MAXDIMS=9


class RSSdata:
    def __init__(self, data=None, datatype=None):
        ndim = 0
        dims = []
        if(data is not None):
            ndim = data.ndim
            dims = data.shape
            self.data = data
            if(datatype is not None):
                if (datatype < 2 or datatype > 3):
                    raise TypeError('Invalid datatype: use 2 for DATA2D and 3 for DATA3D formats')
                if(ndim > 2):
                    raise TypeError('For DATA2D and DATA3D formats ndim = 2 is expeted')
        self.data_format=4
        self.header_format=4
        self.Nheader=0
        self.type=0
        if (datatype == 2):
            self.type=2
            self.Nheader=4
        if (datatype == 3):
            self.type=3
            self.Nheader=6
        self.Ndims=ndim;
        self.geomN = np.zeros([MAXDIMS,1], dtype='uint64')
        self.geomD = np.zeros([MAXDIMS,1], dtype='float64')
        self.geomO = np.zeros([MAXDIMS,1], dtype='float64')
        self.fullsize = 1
        for i in range(0,ndim):
            self.geomN[i] = np.uint64(dims[i])
            self.geomD[i] = 1.0
            self.fullsize = self.fullsize * self.geomN[i]

        # Creating coordinates arrays
        self.srcX = np.zeros([int(self.geomN[1]),1], dtype='float32')
        self.srcZ = np.zeros([int(self.geomN[1]),1], dtype='float32')
        self.GroupX = np.zeros([int(self.geomN[1]),1], dtype='float32')
        self.GroupZ = np.zeros([int(self.geomN[1]),1], dtype='float32')
        if(self.Nheader == 6):
            self.srcY = np.zeros([int(self.geomN[1]),1], dtype='float32')
            self.GroupY = np.zeros([int(self.geomN[1]),1], dtype='float32')

    def read(self,filename):
        f = open(filename, 'rb')
        buff= f.read(MAGICNUMBERLENGTH)
    
        # Check if file has rss identifier
        if(buff.decode("utf-8") != MAGICNUMBER):
            print('Not a valid rss file!')
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
#resize = np.zeros([self.Ndims,1], dtype='uint64')
        resize = []
        for i in range(0,self.Ndims):
            resize.append(int(self.geomN[i]))

        # Making size arrays integers
        self.fullsize = int(self.fullsize)

        # Allocate 1d array for data
        self.data = np.zeros([self.fullsize,1], dtype='float32')
        if(self.Nheader):
            self.srcX = np.zeros([int(self.geomN[1]),1], dtype='float32')
            self.srcY = np.zeros([int(self.geomN[1]),1], dtype='float32')
            self.srcZ = np.zeros([int(self.geomN[1]),1], dtype='float32')
            self.GroupX = np.zeros([int(self.geomN[1]),1], dtype='float32')
            self.GroupY = np.zeros([int(self.geomN[1]),1], dtype='float32')
            self.GroupZ = np.zeros([int(self.geomN[1]),1], dtype='float32')
            for i in range(0,int(self.geomN[1])):
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
                            
                self.data[i*int(self.geomN[0]):(i+1)*int(self.geomN[0]),0] = st.unpack('f'*int(self.geomN[0]), f.read(4*int(self.geomN[0])));
        else:
            self.data[0:self.fullsize,0] = st.unpack('f'*self.fullsize, f.read(4*self.fullsize));


        # Reshape data to correct dimensions
        self.data=np.reshape(self.data, resize, order='F')

    def write(self,filename):
        f = open(filename, 'wb')
        f.write(MAGICNUMBER.encode())
    
        # Writting header
        f.write(st.pack('i', self.data_format))
        f.write(st.pack('i', self.header_format))
        f.write(st.pack('i', self.type))
        f.write(st.pack('Q', self.Nheader))
        f.write(st.pack('Q', self.Ndims))
        for i in range(0,MAXDIMS):
            f.write(st.pack('Q', int(self.geomN[i])))
        for i in range(0,MAXDIMS):
            f.write(st.pack('d', (self.geomD[i])))
        for i in range(0,MAXDIMS):
            f.write(st.pack('d', (self.geomO[i])))

        #Estimating sizes for output
        self.fullsize = 1;
        self.Ndims = 0
        for i in range(0, MAXDIMS):
            if(self.geomN[i] > 0):    
                self.fullsize = self.fullsize*int(self.geomN[i]);
                self.Ndims = self.Ndims + 1

        #Writting Coordinates and data
        dataout=self.data.reshape([self.fullsize,1], order='F')
        if(self.Nheader):
            for i in range(0,int(self.geomN[1])):
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
                            
                for j in range(0,int(self.geomN[0])):
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
def oasiscolor():
   colors = [(0.00000000,0.12941176,0.95686275),(0.01569443,0.16758723,0.95686275),(0.03121853,0.20429106,0.95686275),(0.04656920,0.23937529,0.95686275),(0.06174333,0.27269199,0.95686275),(0.07673779,0.30409322,0.95686275),(0.09154947,0.33343103,0.95686275),(0.10617674,0.36053653,0.95689130),(0.12064031,0.38514503,0.95738630),(0.13494283,0.40776276,0.95837717),(0.14907780,0.42909385,0.95970819),(0.16303873,0.44984243,0.96122362),(0.17681911,0.47071263,0.96276774),(0.19041245,0.49240860,0.96418482),(0.20379191,0.51539258,0.96538255),(0.21689444,0.53895227,0.96655133),(0.22978137,0.56286291,0.96772011),(0.24252718,0.58700325,0.96888889),(0.25520635,0.61125206,0.97005767),(0.26789337,0.63548811,0.97122645),(0.28066270,0.65959016,0.97239523),(0.29398959,0.68435952,0.97356401),(0.30818426,0.71060882,0.97473280),(0.32261560,0.73726304,0.97590158),(0.33664122,0.76322313,0.97707036),(0.34961875,0.78739004,0.97823914),(0.36090578,0.80866473,0.97940792),(0.36987909,0.82597959,0.98057811),(0.37708543,0.84020544,0.98180989),(0.38328357,0.85261264,0.98309767),(0.38880259,0.86377531,0.98438955),(0.39397156,0.87426754,0.98563362),(0.39911955,0.88466344,0.98677797),(0.40457563,0.89553711,0.98777069),(0.41074878,0.90767594,0.98858891),(0.41818632,0.92246147,0.98939764),(0.42631079,0.93865861,0.99018930),(0.43436765,0.95459232,0.99091197),(0.44160233,0.96858756,0.99151375),(0.44726031,0.97896927,0.99194271),(0.45058702,0.98406243,0.99214697),(0.45173745,0.98431373,0.98777763),(0.45265475,0.98431373,0.96974354),(0.45343339,0.98431373,0.94126738),(0.45406041,0.98431373,0.90616450),(0.45452281,0.98431373,0.86825030),(0.45480763,0.98431373,0.83134014),(0.45490196,0.98431336,0.79924672),(0.45490196,0.98381947,0.77096383),(0.45490196,0.98260019,0.74310860),(0.45490196,0.98096701,0.71532774),(0.45490196,0.97923137,0.68726797),(0.45490196,0.97770474,0.65857603),(0.45490196,0.97669857,0.62889864),(0.45492720,0.97649582,0.59748583),(0.45529515,0.97686378,0.55955646),(0.45598071,0.97754934,0.51695002),(0.45682814,0.97839677,0.47371144),(0.45768173,0.97925036,0.43388560),(0.45838573,0.97995436,0.40151742),(0.45878442,0.98035304,0.38065180),(0.45882353,0.98039216,0.37234602),(0.45882353,0.98039216,0.36671016),(0.45882353,0.98039216,0.36224873),(0.45882353,0.98039216,0.35862314),(0.45882353,0.98039216,0.35549479),(0.45882353,0.98039216,0.35252509),(0.45882353,0.98039216,0.34937544),(0.45882353,0.98039216,0.34614766),(0.45882353,0.98039216,0.34326694),(0.45882353,0.98039216,0.34058429),(0.45882353,0.98039216,0.33794499),(0.45882353,0.98039216,0.33519431),(0.45882353,0.98039216,0.33217751),(0.45882353,0.98039216,0.32871554),(0.45882353,0.98039216,0.32393732),(0.45882353,0.98039216,0.31803054),(0.45882353,0.98039216,0.31176185),(0.45882353,0.98039216,0.30589794),(0.45882353,0.98039216,0.30120546),(0.45882353,0.98039216,0.29845110),(0.45941931,0.98039216,0.29809508),(0.46412035,0.98039216,0.29852353),(0.47242434,0.98039216,0.29924500),(0.48303354,0.98039216,0.30010376),(0.49465021,0.98039216,0.30094408),(0.50597658,0.98039216,0.30161023),(0.51571493,0.98039216,0.30194647),(0.52364573,0.98025153,0.30196078),(0.53155665,0.97972224,0.30196078),(0.53933369,0.97894502,0.30196078),(0.54676921,0.97807558,0.30196078),(0.55365559,0.97726966,0.30196078),(0.55978517,0.97668300,0.30196078),(0.56495016,0.97647107,0.30196126),(0.56899242,0.97664844,0.30213864),(0.57218579,0.97709219,0.30258239),(0.57495294,0.97772446,0.30321465),(0.57771658,0.97846737,0.30395757),(0.58089940,0.97924307,0.30473327),(0.58492411,0.97997369,0.30546389),(0.59035228,0.98060175,0.30608228),(0.59894977,0.98130196,0.30666667),(0.61048494,0.98207184,0.30725106),(0.62412169,0.98283352,0.30783545),(0.63902389,0.98350914,0.30841984),(0.65435542,0.98402083,0.30900423),(0.66928015,0.98429073,0.30958862),(0.68341956,0.98431373,0.31016252),(0.69794789,0.98431373,0.31069891),(0.71291032,0.98431373,0.31121960),(0.72821312,0.98431373,0.31175053),(0.74376256,0.98431373,0.31231766),(0.75946493,0.98431373,0.31294694),(0.77522650,0.98431373,0.31366434),(0.79109493,0.98445134,0.31450941),(0.80719891,0.98486449,0.31549777),(0.82349694,0.98547640,0.31660881),(0.83994651,0.98620920,0.31782175),(0.85650513,0.98698504,0.31911584),(0.87313029,0.98772604,0.32047030),(0.88989206,0.98836206,0.32187967),(0.90958111,0.98904731,0.32371954),(0.93171693,0.98981247,0.32595159),(0.95400183,0.99057967,0.32829553),(0.97413811,0.99127106,0.33047103),(0.98982808,0.99180878,0.33219779),(0.99877405,0.99211494,0.33319548),(0.99987424,0.99126621,0.33302904),(0.99895574,0.98476912,0.33081480),(0.99734261,0.97338019,0.32694886),(0.99526843,0.95877446,0.32201818),(0.99296681,0.94262698,0.31660975),(0.99067133,0.92661279,0.31131056),(0.98861560,0.91240693,0.30670758),(0.98685121,0.90042485,0.30294841),(0.98509804,0.88871431,0.29934828),(0.98334487,0.87712273,0.29584462),(0.98159170,0.86557412,0.29240149),(0.97983852,0.85399252,0.28898294),(0.97808535,0.84230198,0.28555305),(0.97633138,0.83042226,0.28207505),(0.97445139,0.81764428,0.27839962),(0.97243617,0.80397356,0.27454388),(0.97040250,0.79010407,0.27066354),(0.96846719,0.77672974,0.26691435),(0.96674703,0.76454453,0.26345202),(0.96535882,0.75424238,0.26043230),(0.96439815,0.74638958,0.25797861),(0.96369606,0.73999790,0.25583976),(0.96313741,0.73447668,0.25390213),(0.96265731,0.72954559,0.25213458),(0.96219087,0.72492431,0.25050595),(0.96167320,0.72033255,0.24898510),(0.96103942,0.71548998,0.24754089),(0.96022390,0.71038767,0.24619878),(0.95921277,0.70565459,0.24508163),(0.95805548,0.70108379,0.24409726),(0.95680394,0.69639496,0.24313667),(0.95551004,0.69130779,0.24209083),(0.95422571,0.68554197,0.24085074),(0.95300285,0.67881719,0.23930738),(0.95178916,0.66973500,0.23691755),(0.95050885,0.65760235,0.23340439),(0.94921352,0.64364937,0.22923347),(0.94795507,0.62911034,0.22487199),(0.94678542,0.61521955,0.22078713),(0.94575646,0.60321127,0.21744609),(0.94491669,0.59423480,0.21528253),(0.94421494,0.58708679,0.21380323),(0.94359719,0.58079213,0.21260712),(0.94303750,0.57511549,0.21158517),(0.94250992,0.56982157,0.21062839),(0.94198848,0.56467505,0.20962775),(0.94144723,0.55944062,0.20847426),(0.94086890,0.55403754,0.20707715),(0.94028451,0.54905310,0.20548666),(0.93970012,0.54434334,0.20375847),(0.93911572,0.53963779,0.20193930),(0.93853133,0.53466596,0.20007587),(0.93794694,0.52915739,0.19821490),(0.93736255,0.52284159,0.19640311),(0.93677816,0.51519455,0.19467003),(0.93619377,0.50585567,0.19298480),(0.93560938,0.49527506,0.19131362),(0.93502499,0.48391990,0.18962311),(0.93444060,0.47225738,0.18787992),(0.93385621,0.46075468,0.18605067),(0.93327214,0.44987557,0.18410215),(0.93271705,0.43958260,0.18202712),(0.93219133,0.42951620,0.17984137),(0.93166903,0.41947614,0.17755480),(0.93112419,0.40926221,0.17517729),(0.93053085,0.39867417,0.17271873),(0.92986307,0.38751180,0.17018900),(0.92906109,0.37550409,0.16756020),(0.92785421,0.36203417,0.16454868),(0.92635636,0.34747361,0.16126346),(0.92477518,0.33240841,0.15791465),(0.92331831,0.31742461,0.15471236),(0.92219339,0.30310822,0.15186672),(0.92160807,0.29004526,0.14958782),(0.92156863,0.27824352,0.14791911),(0.92156863,0.26643545,0.14645448),(0.92156863,0.25490795,0.14514734),(0.92156863,0.24406808,0.14398285),(0.92156863,0.23432292,0.14294619),(0.92156863,0.22607955,0.14202252),(0.92156863,0.21974505,0.14119700),(0.92134490,0.21457692,0.14039825),(0.92075047,0.20960151,0.13959204),(0.91994089,0.20506175,0.13883022),(0.91907188,0.20120139,0.13816471),(0.91829918,0.19826419,0.13764741),(0.91777852,0.19649388,0.13733024),(0.91766468,0.19609026,0.13757306),(0.91801160,0.19633018,0.14409491),(0.91868433,0.19682009,0.15763753),(0.91952713,0.19748212,0.17635676),(0.92038427,0.19823841,0.19840843),(0.92110002,0.19901109,0.22194838),(0.92151866,0.19972229,0.24513244),(0.92156863,0.20030831,0.26691552),(0.92156863,0.20081053,0.29006913),(0.92156863,0.20128078,0.31496995),(0.92156863,0.20176578,0.34152956),(0.92156863,0.20231226,0.36965954),(0.92156863,0.20296693,0.39927148),(0.92156863,0.20377651,0.43027694),(0.92139886,0.20496278,0.46476217),(0.92084433,0.20677247,0.50505849),(0.92005416,0.20906027,0.54857902),(0.91918408,0.21167391,0.59264971),(0.91838981,0.21446114,0.63459648),(0.91782709,0.21726968,0.67174526),(0.91764706,0.21994483,0.70147319),(0.91764706,0.22225779,0.72551033),(0.91764706,0.22434164,0.74657297),(0.91764706,0.22649122,0.76560173),(0.91764706,0.22900136,0.78353727),(0.91764706,0.23216691,0.80132020),(0.91764706,0.23628269,0.81989117),(0.91777496,0.24373393,0.84047445),(0.91898419,0.26996005,0.86520302),(0.92118047,0.31248967,0.89222213),(0.92400788,0.36450184,0.91894362),(0.92711044,0.41917558,0.94277934),(0.93013222,0.46968991,0.96114113),(0.93271725,0.50922386,0.97144083),(0.93474093,0.53609892,0.97405211),(0.93666495,0.56093310,0.97592104),(0.93851725,0.58395169,0.97748693),(0.94029227,0.60461412,0.97873338),(0.94198445,0.62237982,0.97964399),(0.94358823,0.63670825,0.98020238),(0.94509804,0.64705882,0.98039216)]
   cm = LinearSegmentedColormap.from_list('oasiscolor', colors, N=256)
   return cm
