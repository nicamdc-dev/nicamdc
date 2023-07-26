################################################################################
#
# Makefile for each test program
#
################################################################################
PWD       = $(shell pwd)
TOPDIR    = $(abspath ../../../..)
TESTDIR   = ../../..


# parameters for run
glevel     = 9
rlevel_in  = 2
nmpi_in    = 32
rlevel_out = 6
nmpi_out   = 8192

BINNAME   = fio_pe2pe

# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common

jobshell:
	sh ../Mkjobshell.$(NICAM_SYS).sh $(glevel) $(rlevel_in) $(nmpi_in) $(rlevel_out) $(nmpi_out) $(TOPDIR) $(BINNAME); \
