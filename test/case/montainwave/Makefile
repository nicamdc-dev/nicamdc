################################################################################
#
# Makefile for testcase
#
################################################################################

TOPDIR = ../../..

include $(TOPDIR)/sysdep/Makedef.$(NICAM_SYS)
include $(TOPDIR)/src/Mkinclude

include $(TESTDIR)/test.conf
RUNDIR := $(shell sh $(TESTDIR)/util/dirname.sh $(glevel) $(rlevel) $(nmpi) $(zlayer))
testcase = montainwave

all: build

makedir:
	mkdir -p $(RUNDIR)

build: makedir
	$(MAKE) -C $(TOPDIR)/src || exit 1

jobshell: build checkconfig
	sh $(TESTDIR)/sysdep/Mkjobshell.$(NICAM_SYS).sh $(glevel) $(rlevel) $(nmpi) $(zlayer) $(vgrid) $(LSMAX) $(DTL) $(DIFCF) $(NHIST)
	sh Mkconf.$(testcase).sh $(glevel) $(rlevel) $(nmpi) $(zlayer) $(vgrid) $(LSMAX) $(DTL) $(DIFCF) $(NHIST)

run: jobshell
	cd $(RUNDIR); $(JOBSUB) run.sh

checkconfig:
	@echo "##### Run configuration #####"
	@echo "+system           :" $(NICAM_SYS)
	@echo "+testcase         :" $(testcase)
	@echo "+glevel           :" $(glevel)
	@echo "+rlevel           :" $(rlevel)
	@echo "+MPI process      :" $(nmpi)
	@echo "+zlayer           :" $(zlayer)
	@echo "+vgrid            :" $(vgrid)
	@echo "+large step max   :" $(LSMAX)
	@echo "+DTL              :" $(DTL)
	@echo "+diffusion coef.  :" $(DIFCF)
	@echo "+history steps    :" $(NHIST)
	@echo
