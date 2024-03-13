# SConstruct initializing
import os
import sys

# Paths
root = os.getcwd()
libdir = 'lib'
blddir = os.path.join(root,'build')
incdir = 'include'
bindir = os.path.join(root,'bin')
srcdir = os.path.join(root,'src')
devdir = os.path.join(root,'dev')

# System dependent flags
if sys.platform == 'darwin' :    # Mac
    env = Environment(ENV=os.environ,CPPFLAGS='-std=c++11 -stdlib=libc++') 
else:
    env = Environment(ENV=os.environ,CPPFLAGS='-std=c++11')


# Setup
env['program_prefix'] = "rs" # Prefix for all programs

# Configure (checking libraries for environment)
## Build GUI programs 
env['WX'] = False;
## Include FFTW programs
env['FFTW'] = False;
## Build test codes
env['TEST'] = True;
## MPI
env['MPI'] = True;

# Compiler flags
env.Append(CPPFLAGS=' -O3 -Wall' ) # Optimized 
#env.Append(CPPFLAGS=' -g -Wall' ) # For debugging

# Include path
env.Append(CPPPATH=['../include'])
env.Append(CPPPATH=['../config4cpp/include'])
env.Append(CPPPATH=['../config4cpp/src'])
env.Append(CPPPATH=['../madagascar/include'])
#env.Append(CPPPATH=['../../fftw-3.3.6-pl2/fftw-build/include/'])
env.Append(LIBPATH=['../build'])
env.Append(LIBPATH=['../config4cpp/lib'])
env.Append(LIBPATH=['../madagascar/build'])

#Create other environments
if sys.platform == 'darwin' :    # Mac
    mpi = env.Clone(CC='mpicc', CXX='mpic++')
else:
    mpi = env.Clone(CC='mpicc', CXX='mpicxx')
    
envc = Environment(ENV=os.environ, CC = 'gcc', CXX='g++') 
envc.Append(CPPFLAGS=' -O3 -Wall' ) # Optimized 
#envc.Append(CPPFLAGS=' -g -Wall' ) # For debugging

# Programs
SConscript('lib/SConscript', 'env mpi')
SConscript('src/SConscript','env mpi')
SConscript('dev/SConscript','env')
SConscript('gui/SConscript','env')
SConscript('config4cpp/src/SConscript','env')
SConscript('madagascar/lib/SConscript','env envc')
