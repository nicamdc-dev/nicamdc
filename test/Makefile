-include ../src/Mkinclude
-include ../src/sysdep/Mkinclude.$(NICAM_SYS)
include test.conf

ifeq ($(origin NICAMRUNDIR), undefined)
	export NICAMRUNDIR = $(NICAMDIR)/test
endif

all: heldsuarez jablonowski tracer

ifeq ($(runnow), yes)
rundir := $(shell sh ./util/dirname.sh $(glevel) $(rlevel) $(nmpi) $(zlayer))
else
runnow="no"
endif


heldsuarez: makedir checkconfig
	@echo "prepare heldsuarez test"
	@export NICAMDIR="$(NICAMDIR)";       \
	 export NICAMBINDIR="$(NICAMBINDIR)"; \
	 export NICAMRUNDIR="$(NICAMRUNDIR)"; \
	 export MPIEXEC="$(MPIEXEC)";	        \
	 sh sysdep/Mkshell_$(NICAM_SYS).sh $(glevel) $(rlevel) $(nmpi) $(zlayer) $(vgrid) heldsuarez
	@export NICAMRUNDIR="$(NICAMRUNDIR)"; \
	 sh case/Mkconf_heldsuarez.sh $(glevel) $(rlevel) $(nmpi) $(zlayer) $(vgrid) heldsuarez
ifeq ($(runnow), yes)
	@echo "run immediately"
	cd $(NICAMRUNDIR)/jablonowski/$(rundir); $(JOBSUB) run.sh
endif

jablonowski: makedir checkconfig
	@echo "prepare jablonowski test"
	@export NICAMDIR="$(NICAMDIR)";       \
	 export NICAMBINDIR="$(NICAMBINDIR)"; \
	 export NICAMRUNDIR="$(NICAMRUNDIR)"; \
	 export MPIEXEC="$(MPIEXEC)";	        \
	 sh sysdep/Mkshell_$(NICAM_SYS).sh $(glevel) $(rlevel) $(nmpi) $(zlayer) $(vgrid) $(LSMAX) $(DTL) $(DIFCF) $(NHIST) jablonowski
	@export NICAMRUNDIR="$(NICAMRUNDIR)"; \
	 sh case/Mkconf_jablonowski.sh $(glevel) $(rlevel) $(nmpi) $(zlayer) $(vgrid) $(LSMAX) $(DTL) $(DIFCF) $(NHIST) jablonowski
ifeq ($(runnow), yes)
	@echo "run immediately"
	cd $(NICAMRUNDIR)/jablonowski/$(rundir); $(JOBSUB) run.sh
endif

tracer: makedir checkconfig
	@echo "prepare tracer test"
	@export NICAMDIR="$(NICAMDIR)";       \
	 export NICAMBINDIR="$(NICAMBINDIR)"; \
	 export NICAMRUNDIR="$(NICAMRUNDIR)"; \
	 export MPIEXEC="$(MPIEXEC)";	        \
	 sh sysdep/Mkshell_$(NICAM_SYS).sh $(glevel) $(rlevel) $(nmpi) $(zlayer) $(vgrid) tracer
	@export NICAMRUNDIR="$(NICAMRUNDIR)"; \
	 sh case/Mkconf_tracer.sh $(glevel) $(rlevel) $(nmpi) $(zlayer) $(vgrid) tracer
ifeq ($(runnow), yes)
	@echo "run immediately"
	cd $(NICAMRUNDIR)/tracer/$(rundir); $(JOBSUB) run.sh
endif

makedir:
	mkdir -p $(NICAMRUNDIR)

checkconfig:
	@echo "+glevel           :" $(glevel)
	@echo "+rlevel           :" $(rlevel)
	@echo "+zlayer           :" $(zlayer)
	@echo "+vgrid            :" $(vgrid)
	@echo "+MPI process      :" $(nmpi)
	@echo "+system           :" $(NICAM_SYS)
	@echo "+run immediately? :" $(runnow)
	@echo "+large step max   :" $(LSMAX)
	@echo "+DTL              :" $(DTL)
	@echo "+diffusion coef.  :" $(DIFCF)
	@echo "+history steps    :" $(NHIST)
	@echo
