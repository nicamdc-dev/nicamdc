################################################################################
#
# Makefile for framework test
#
################################################################################

TOPDIR = ../../..

include $(TOPDIR)/sysdep/Makedef.$(NICAM_SYS)
include $(TOPDIR)/src/Mkinclude

include $(TESTDIR)/test.conf
RUNDIR := $(shell sh $(TESTDIR)/util/dirname.sh $(glevel) $(rlevel) $(nmpi) $(zlayer))
testcase = mkmnginfo

all: depends buildlocal

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
	sh Mkjobshell.$(NICAM_SYS).sh $(glevel) $(rlevel) $(nmpi) $(zlayer) $(vgrid) ../${TOPDIR} ${testcase}
	sh Mkconf.$(testcase).sh      $(glevel) $(rlevel) $(nmpi) $(zlayer) $(vgrid) ../${TOPDIR} ${testcase}

run: jobshell
	cd $(RUNDIR); $(JOBSUB) run.sh

################################################################################
# Make in local
################################################################################
VPATH = $(SRCDIR)/share:$(SRCDIR)/nhm:$(SRCDIR)/prep/mkmnginfo:$(SRCDIR)/prep/mkgrid:$(SRCDIR)/prep/mkvlayer

PRJ = $(testcase)

buildlocal: $(BINDIR)/$(PRJ)

$(BINDIR)/$(PRJ): $(PRJ)
	mkdir -p $(BINDIR)
	install $< $@
	@echo; echo "##### $< is Installed Successfully #####"; echo

$(PRJ) : $(PRJ).o
	$(LD) $(LDFLAGS) -o $@ $^ $(MODS)

.SUFFIXES: .o .f90 .c
.f90.o:
	$(FC)  $(FFLAGS) -o $@ -c $<
.c.o:
	$(CC)  $(CFLAGS) -o $@ -c $<

%.o: %.mod

depends.info:
	@$(TESTDIR)/util/mkdepends.pl $(PRJ) $(VPATH)

################################################################################

.PHONY: depends freeze clean copylocal mktar

clean:
	rm -f $(PRJ) depends.info
	$(CLEAN)

freeze: clean depends copylocal mktar

copylocal:
	@for i in $(MODSRC) ; do \
		cp $$i . ; \
	done

mktar:
	tar -czf ../$(PRJ).tgz ./

depends: depends.info

-include depends.info
