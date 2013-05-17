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
testcase = nhm_driver

all: build

getrundir:
	@echo $(RUNDIR)

makedir:
	mkdir -p $(RUNDIR)

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

jobshell: makedir checkconfig
	sh $(TOPDIR)/sysdep/Mkjobshell.$(NICAM_SYS).sh $(glevel) $(rlevel) $(nmpi) $(zlayer) $(vgrid) ../${TOPDIR} ${testcase}
	sh Mkconf.$(testcase).sh                       $(glevel) $(rlevel) $(nmpi) $(zlayer) $(vgrid) ../${TOPDIR} ${testcase} $(LSMAX) $(DTL) $(DIFCF) $(NHIST)

run: jobshell
	cd $(RUNDIR); $(JOBSUB) run.sh

build:
	$(MAKE) -C $(TOPDIR)/src || exit 1
