################################################################################
#
# Makefile for testcase
#
################################################################################

TOPDIR = ../../../..

include $(TOPDIR)/sysdep/Makedef.$(NICAM_SYS)
include $(TOPDIR)/src/Mkinclude

# for run
exename = nhm_driver
glevel  = 5
rlevel  = 0
nmpi    = 10
zlayer  = 40
vgrid   = vgrid40_24000-600m.dat



all: build

# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
