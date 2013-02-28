-include ../Mkinclude
-include ../../sysdep/Mkinclude.$(NICAM_SYS)

VPATH = ../share

PRJ12 = mkllmap
PRJ13 = ico2ll
PRJ26 = fio_dump
PRJ30 = fio_ico2ll
PRJ37 = fio_ico2ll_mpi

PRJS = $(PRJ12) $(PRJ13) $(PRJ26) $(PRJ30) $(PRJ37)

MODS_TOOL =	\
	mod_mnginfo_light.o	\
	mod_option.o

all: $(PRJS)
	$(INSTALL) $(PRJS) $(NICAMBINDIR)
	@echo; echo "##### Tools are Installed Successfully #####"; echo

$(PRJ12): prg_$(PRJ12).o $(MODS_TOOL) $(NICAMLIBDIR)/$(LIBNAME)
	$(LD) $(LFLAGS) -o $@ $^
$(PRJ13): prg_$(PRJ13).o $(MODS_TOOL) $(NICAMLIBDIR)/$(LIBNAME)
	$(LD) $(LFLAGS) -o $@ $^
$(PRJ26): prg_$(PRJ26).o $(MODS_TOOL) $(NICAMLIBDIR)/$(LIBNAME)
	$(LD) $(LFLAGS) -o $@ $^
$(PRJ30): prg_$(PRJ30).o $(MODS_TOOL) $(NICAMLIBDIR)/$(LIBNAME)
	$(LD) $(LFLAGS) -o $@ $^
$(PRJ37): prg_$(PRJ37).o $(MODS_TOOL) $(NICAMLIBDIR)/$(LIBNAME)
	$(LD) $(LFLAGS) -o $@ $^

prg_$(PRJ12).o	: prg_$(PRJ12).f90 $(NICAMLIBDIR)/$(LIBNAME)
prg_$(PRJ13).o	: prg_$(PRJ13).f90 $(NICAMLIBDIR)/$(LIBNAME)
prg_$(PRJ26).o	: prg_$(PRJ26).f90 $(NICAMLIBDIR)/$(LIBNAME)
prg_$(PRJ30).o	: prg_$(PRJ30).f90 $(NICAMLIBDIR)/$(LIBNAME)
prg_$(PRJ37).o	: prg_$(PRJ37).f90 $(NICAMLIBDIR)/$(LIBNAME)

mod_mnginfo_light.o	: mod_mnginfo_light.f90 $(NICAMLIBDIR)/$(LIBNAME)
mod_option.o	: mod_option.f90



.SUFFIXES: .o .f90 .c
.f90.o:
	$(FC) $(FFLAGS) -I$(NICAMDIR)/src -o $@ -c $<
.c.o:
	$(CC) $(CFLAGS) -I$(NICAMDIR)/src -o $@ -c $<

%.o: %.mod

.PHONY: clean

clean:
	rm -f $(PRJS)
	$(CLEAN)
