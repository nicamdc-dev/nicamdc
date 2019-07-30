################################################################################
#
# Makefile for each test program
#
################################################################################
PWD       = $(shell pwd)
TOPDIR    = $(abspath ../../../..)
TESTDIR   = ../../..


# parameters for run
glevel     = 6
rlevel_in  = 0
nmpi_in    = 1
rlevel_out = 1
nmpi_out   = 8

BINNAME   = fio_pe2pe

# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common

jobshell:
	sh ../Mkjobshell.$(NICAM_SYS).sh $(glevel) $(rlevel_in) $(nmpi_in) $(rlevel_out) $(nmpi_out) $(TOPDIR) $(BINNAME); \
