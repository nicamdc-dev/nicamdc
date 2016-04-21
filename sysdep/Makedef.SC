#
# ------ FOR JAMSTEC SC system -----
#

FFLAGS_FAST  = -fpp3 -O3 -m64 -axAVX -ip -traceback           \
               -assume byterecl -convert big_endian \
               -ftz -fp-model precise -pc 80        \
               -mcmodel=medium -shared-intel        \
               -heap-arrays -fno-alias -lmpi

FFLAGS_DEBUG = -fpp3 -O0 -m64                       \
               -assume byterecl -convert big_endian \
               -fp-model precise -pc 80             \
               -check all -warn all -g -traceback   \
               -mcmodel=medium -shared-intel        \
               -heap-arrays -fno-alias -lmpi

FC     = ifort
FFLAGS = $(FFLAGS_FAST)
#FFLAGS = $(FFLAGS_DEBUG)

CC     = icc
CFLAGS = -O3 -axAVX -ip -ftz -mcmodel=medium -shared-intel
NETCDF4 = /home/G10202/masuo/cc-env
HDF5 = /home/G10202/masuo/cc-env


##### Special Options #####

ifeq ($(ENABLE_OPENMP),T)
	FFLAGS += -openmp -openmp-report2
endif

ifeq ($(ENABLE_NETCDF),T)
	NETCDF_INCLUDE ?= -I$(NETCDF4)/include
	NETCDF_LIBS    ?= -L$(NETCDF4)/lib -L$(HDF5)/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lm -lz
	FFLAGS         += $(NETCDF_INCLUDE)
	CFLAGS         += $(NETCDF_INCLUDE)
endif

###########################

MODDIROPT ?= -module $(BUILDDIR)
OBJDIR    ?= $(BUILDDIR)

LD     = $(FC)
LFLAGS = $(FFLAGS)

MAKE    = make
INSTALL = cp -u
AR      = ar
ARFLAGS = r
RANLIB  = ranlib
