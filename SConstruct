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

#conf = Configure(env)
#if conf.CheckLibWithHeader('wx','wx.h','cpp'):
#	print "Yes"
#
#
## CUDA
#if "CUDA_HOME" in os.environ:
#	env['CXX'] = 'nvcc'
#	env['BUILD_CUDA'] = True
#	env['CUDA_TOOLKIT_PATH'] = os.environ['CUDA_HOME']
#	env['CUDA_SDK_PATH'] = os.environ['CUDA_HOME']
#	env.Tool('cuda', toolpath='./')
#	
#	#env['NVCCFLAGS'] = '-arch=sm_52 '
#
#	# OpenMP
#	if conf.CheckLibWithHeader('pthread','pthread.h','cpp'):
#		env.Append(CPPFLAGS = ' -Xcompiler -fopenmp')
#		env.Append(LINKFLAGS = '-lgomp')
#
#else:
#	env['BUILD_CUDA'] = False
#
#	# OpenMP
#	if conf.CheckLibWithHeader('pthread','pthread.h','cpp'):
#		env.Append(CPPFLAGS = ' -fopenmp')
#		env.Append(LINKFLAGS = ' -lpthread -fopenmp')
#
### TODO: Add autoconfig of MPI
#env['MPI'] = True;
#
## Finishing configure
#env = conf.Finish()



# Compiler flags
env.Append(CPPFLAGS=' -O3 -Wall' ) # Optimized 
#env.Append(CPPFLAGS=' -g -Wall' ) # For debugging

# Include path
env.Append(CPPPATH=['../include'])
env.Append(CPPPATH=['../config4cpp/include'])
env.Append(CPPPATH=['../config4cpp/src'])
env.Append(CPPPATH=['../madagascar/include'])
env.Append(CPPPATH=['../lbfgs'])
env.Append(CPPPATH=['../../fftw-3.3.6-pl2/fftw-build/include/'])
env.Append(LIBPATH=['../build'])

# Programs
SConscript('lib/SConscript','env')
SConscript('src/SConscript','env')
SConscript('dev/SConscript','env')
SConscript('gui/SConscript','env')
SConscript('config4cpp/src/SConscript','env')
SConscript('madagascar/lib/SConscript','env')
SConscript('lbfgs/SConscript','env')
#SConscript('python/SConscript','env')

