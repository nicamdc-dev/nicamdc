################################################################################
#
# Makefile for tools
#
################################################################################

TOPDIR = ../..

include $(TOPDIR)/sysdep/Makedef.$(NICAM_SYS)
include $(TOPDIR)/src/Mkinclude

PRJ13 = ico2ll
PRJ26 = fio_dump
PRJ30 = fio_ico2ll
PRJ37 = fio_ico2ll_mpi

PRJS = $(PRJ13) $(PRJ26) $(PRJ30) $(PRJ37)

MODS_TOOL =	\
	mod_mnginfo_light.o	\
	mod_option.o

all: $(PRJS)
	$(INSTALL) $(PRJS) $(BINDIR)
	@echo; echo "##### Tools are Installed Successfully #####"; echo

$(LIBDIR)/$(LIBNAME):
	$(MAKE) -C $(TOPDIR)/src || exit 1

$(PRJ13): prg_$(PRJ13).o $(MODS_TOOL) $(LIBDIR)/$(LIBNAME)
	$(LD) $(LFLAGS) -o $@ $^
$(PRJ26): prg_$(PRJ26).o $(MODS_TOOL) $(LIBDIR)/$(LIBNAME)
	$(LD) $(LFLAGS) -o $@ $^
$(PRJ30): prg_$(PRJ30).o $(MODS_TOOL) $(LIBDIR)/$(LIBNAME)
	$(LD) $(LFLAGS) -o $@ $^
$(PRJ37): prg_$(PRJ37).o $(MODS_TOOL) $(LIBDIR)/$(LIBNAME)
	$(LD) $(LFLAGS) -o $@ $^

prg_$(PRJ13).o	: prg_$(PRJ13).f90 $(LIBDIR)/$(LIBNAME)
prg_$(PRJ26).o	: prg_$(PRJ26).f90 $(LIBDIR)/$(LIBNAME)
prg_$(PRJ30).o	: prg_$(PRJ30).f90 $(LIBDIR)/$(LIBNAME)
prg_$(PRJ37).o	: prg_$(PRJ37).f90 $(LIBDIR)/$(LIBNAME)

mod_mnginfo_light.o	: mod_mnginfo_light.f90 $(LIBDIR)/$(LIBNAME)
mod_option.o	: mod_option.f90



.SUFFIXES: .o .f90 .c
.f90.o:
	$(FC) $(FFLAGS) -I$(TOPDIR)/src -o $@ -c $<
.c.o:
	$(CC) $(CFLAGS) -I$(TOPDIR)/src -o $@ -c $<

%.o: %.mod

.PHONY: clean

clean:
	rm -f $(PRJS) *~ *.h *.mod *.o *__genmod.f90 *.lst *.L
