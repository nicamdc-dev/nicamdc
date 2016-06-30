################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../../..)
TESTDIR     = ../../..
SYSDEP_DIR  = $(abspath ..)

# user-defined source files
CODE_DIR    = .
ORG_SRCS    =

# parameters for run
glevel      = 7
rlevel      = 1
nmpi        = 8
zlayer      = 78
vgrid       = vgrid78.dat



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
