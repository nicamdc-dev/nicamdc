################################################################################
#
# Makefile grand central
#
################################################################################

TOPDIR      = $(abspath ..)
BUILD_DIR   = ./.libs
SYSDEP_DIR  = ../sysdep

VPATH = $(BUILD_DIR)

include $(SYSDEP_DIR)/Makedef.$(NICAM_SYS)
include $(TOPDIR)/Mkinclude

all:
	$(MAKE) makedir
	$(MAKE) makelib
	$(MAKE) makebin

allclean: clean cleandir

# share libraries
include $(SRCDIR)/share/dependency
VPATH +=:$(VPATH_SHARE)
PRJS  += $(PRJ_SHARE)
LIBS  += $(LIB_SHARE)

# prep libraries & programs
include $(SRCDIR)/prep/dependency
VPATH +=:$(VPATH_PREP)
PRJS  += $(PRJ_PREP)
LIBS  += $(LIB_PREP)

# nhm libraries & programs
include $(SRCDIR)/nhm/dependency
VPATH +=:$(VPATH_NHM)
PRJS  += $(PRJ_NHM)
LIBS  += $(LIB_NHM)

# tool libraries & programs
include $(SRCDIR)/tool/dependency
VPATH +=:$(VPATH_TOOL)
PRJS  += $(PRJ_TOOL)
LIBS  += $(LIB_TOOL)

makedir:
	mkdir -p $(BUILD_DIR)
	mkdir -p $(NICAM_BINDIR)
	mkdir -p $(NICAM_LIBDIR)

makelib: makedir $(addprefix $(NICAM_LIBDIR)/,$(LIBS))

makebin: makedir $(addprefix $(NICAM_BINDIR)/,$(PRJS))

nhm    : makedir $(addprefix $(NICAM_BINDIR)/,$(PRJ_NHM))
prep   : makedir $(addprefix $(NICAM_BINDIR)/,$(PRJ_PREP))
tool   : makedir $(addprefix $(NICAM_BINDIR)/,$(PRJ_TOOL))

cleandir:
	-rm $(NICAM_LIBDIR)/*
	-rm $(NICAM_BINDIR)/*

clean:
	rm -f $(LIBS) $(PRJS)
	rm -f *~ *.mod *.o *.lst *.L
	cd $(BUILD_DIR); rm -f $(LIBS) $(PRJS); rm -f *~ *.mod *.o *.lst *.L

.SUFFIXES: .o .f90 .c .mod
%.o : %.f90
	$(FC) $(FFLAGS) -I$(BUILD_DIR) $(MODDIROPT) -o $(OBJDIR)/$@ -c $<
%.o : %.c
	$(CC) $(CFLAGS) -o $(OBJDIR)/$@ -c $<

%.o : %.mod

%.mod : %.f90
	$(MAKE) $(patsubst %.f90,%.o,$<)

.PHONY : clean cleandir allclean share prep nhm tool
