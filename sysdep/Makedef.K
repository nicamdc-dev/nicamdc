#
# ------ FOR K computer, FX10 -----
#

FFLAGS_FAST  = -Kfast,parallel,auto,ocl,preex,array_private,mfunc=2 \
               -Qi -Qt -X03 -Ncompdisp -Koptmsg=1 -Cpp              \
               -Ksimd

FFLAGS_DEBUG = -O0 \
               -Qi -Qt -X03 -v03s -v03d -v03o -Ncompdisp -Koptmsg=1 -Cpp \
               -Ec -Eg -Ha -He -Hf -Ho -Hs -Hu -Hx -Ncheck_global

FC     = mpifrtpx
FFLAGS = $(FFLAGS_FAST)
#FFLAGS = $(FFLAGS_DEBUG)

CC     = mpifccpx
CFLAGS = -Kfast,parallel,ocl,preex,array_private,region_extension,restp=all -Qt -Ksimd $(PERF_MONIT)

##### Special Options #####

# Performance monitor
# disable
PERF_MONIT = -Ntl_notrt -U_FIPP_ -U_FAPP_
# fipp
#PERF_MONIT = -Ntl_trt -D_FIPP_ -U_FAPP_
# fapp
#PERF_MONIT = -Ntl_trt -U_FIPP_ -D_FAPP_

FFLAGS += $(PERF_MONIT)
CFLAGS += $(PERF_MONIT)

ifeq ($(ENABLE_OPENMP),T)
	FFLAGS += -Kopenmp
endif

ifeq ($(ENABLE_NETCDF),T)
	NETCDF_INCLUDE ?= -I/opt/aics/netcdf/k-serial-noszip/include
	NETCDF_LIBS    ?= -L/opt/aics/netcdf/k-serial-noszip/lib-static -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -lm
	FFLAGS += $(NETCDF_INCLUDE)
endif

###########################

MODDIROPT ?= -M $(BUILD_DIR)
OBJDIR    ?= $(BUILD_DIR)

LD     = $(FC)
LFLAGS = $(FFLAGS)

INSTALL = install
AR      = ar
ARFLAGS = r
RANLIB  = ranlib
JOBSUB  = pjsub
