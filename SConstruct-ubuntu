# SConstruct initializing
import os
import sys

#--- Paths
root = os.getcwd()
libdir = 'lib'
blddir = os.path.join(root,'build')
incdir = 'include'
bindir = os.path.join(root,'bin')
srcdir = os.path.join(root,'src')
devdir = os.path.join(root,'dev')

#--- Create standard environment
env = Environment(ENV=os.environ,CPPFLAGS='-std=c++11')

# Set program prefix for all programs
env['program_prefix'] = "rs" 

# Configure libraries
# Build GUI programs 
env['WX'] = False;
# Include FFTW programs
env['FFTW'] = False;

# Include paths
env.Append(CPPPATH=['../include'])
env.Append(CPPPATH=['../config4cpp/include'])
env.Append(CPPPATH=['../config4cpp/src'])
env.Append(CPPPATH=['../madagascar/include'])
env.Append(LIBPATH=['../build'])
env.Append(LIBPATH=['../config4cpp/lib'])
env.Append(LIBPATH=['../madagascar/build'])
if env['FFTW'] == True :
    env.Append(CPPPATH=['../../fftw-3.3.6-pl2/fftw-build/include/'])

#--- Create mpi environment
mpi = env.Clone(CC='mpicc', CXX='mpicxx')
mpi['ON'] = True

#--- Create GPU environment
gpu = env.Clone(CXX='nvcc', CPPFLAGS='-x cu -arch=sm_60')
gpu['ON'] = True

#--- Create clean environment
envc = Environment(ENV=os.environ, CC ='gcc', CXX='g++', CPPFLAGS='-g')

# Include paths 
gpu.Append(CPPPATH=['../include','/usr/lib/x86_64-linux-gnu/openmpi/include'])
gpu.Append(LIBPATH=['../build','/usr/lib/x86_64-linux-gnu/openmpi/lib'])
gpu.Append(CPPPATH=['../config4cpp/include'])
gpu.Append(CPPPATH=['../config4cpp/src'])
gpu.Append(LIBPATH=['../config4cpp/lib'])

#--- Building

# Build libraries
SConscript('lib/SConscript', ['mpi', 'gpu'])

# Build programs
SConscript('src/SConscript',['env', 'mpi', 'gpu'])

# Build gui
SConscript('gui/SConscript','env')

#Build parser
SConscript('config4cpp/src/SConscript','env')

#Build madagscar libraries
#Do not build on Idun
SConscript('madagascar/lib/SConscript',['env','envc'])
