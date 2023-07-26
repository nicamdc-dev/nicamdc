################################################################################
#
# Makefile grand central
#
################################################################################

TOPDIR    = $(abspath ..)
BUILDDIR  = ./.libs
SYSDEPDIR = ../sysdep

VPATH = $(BUILDDIR)

include $(SYSDEPDIR)/Makedef.$(NICAM_SYS)
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
	mkdir -p $(BUILDDIR)
	mkdir -p $(NICAM_BINDIR)
	mkdir -p $(NICAM_LIBDIR)

makelib: makedir $(addprefix $(NICAM_LIBDIR)/,$(LIBS))

makebin: makedir $(addprefix $(NICAM_BINDIR)/,$(PRJS))

nhm    : makedir $(addprefix $(NICAM_BINDIR)/,$(PRJ_NHM))
prep   : makedir $(addprefix $(NICAM_BINDIR)/,$(PRJ_PREP))
tool   : makedir $(addprefix $(NICAM_BINDIR)/,$(PRJ_TOOL))

# parameters for fixed index
glevel  =  5
rlevel  =  0
nmpi    = 10
zlayer  = 38
diamond = 10

fixedindex:
	@cd ./include; ruby ./setindex.rb $(glevel) $(rlevel) $(nmpi) $(zlayer) $(diamond)

cleandir:
	-rm -f $(addprefix $(NICAM_LIBDIR)/,$(LIBS))
	-rm -f $(addprefix $(NICAM_BINDIR)/,$(PRJS))

clean:
	rm -f $(LIBS) $(PRJS)
	rm -f *~ *.mod *.o *.lst *.L
	cd $(BUILDDIR); rm -f $(LIBS) $(PRJS); rm -f *~ *.mod *.o *.lst *.L

.SUFFIXES: .o .f90 .c .mod
%.o : %.f90
	$(FC) $(FFLAGS) -I$(BUILDDIR) $(MODDIROPT) -o $(OBJDIR)/$@ -c $< -I$(TOPDIR)/src/share/io
%.o : %.c
	$(CC) $(CFLAGS) -o $(OBJDIR)/$@ -c $< -I$(TOPDIR)/NICAM/share/io

%.o : %.mod

%.mod : %.f90
	$(MAKE) $(patsubst %.f90,%.o,$<)

.PHONY : clean cleandir allclean share prep nhm tool
