#*********************************************************************#
#*  Makefile for DG 3D program                                       *#
#*                                                                   *#
#*  Author: Ruichao Ye     Email: rye@purdue.edu                     *#
#*********************************************************************#

COMPILER := GNU
STATIC := ON
WithOMP := OFF
DEBUG := ON

#######################################################################
#     list source files and the target names                          #
#######################################################################

SRC_MOD = mod_datatype.f90 mod_jacobi.f90 mod_string.f90 mod_para.f90 \
	mod_meshfile.f90 mod_geometry.f90 mod_ssort.f90 mod_bitree.f90 \
	mod_domain_decomp.f90 mod_parallel.f90 mod_conv_mpi.f90 \
	mod_source.f90 mod_PML.f90 mod_friction.f90 mod_ruptsolve.f90 \
	mod_surface.f90 mod_checkpoint.f90 mod_prem.f90 mod_grav.f90 \
	mod_solve.f90 mod_ODEsolve.f90

SRC_WAVE = DGseis3D.f90
EXE_WAVE := DGseis3Dmpirupt

OBJDIR := ../OBJ
BINDIR := ../bin
SRCDIR := ./
OBJ_MOD :=  $(foreach file,$(SRC_MOD),$(OBJDIR)/$(file:.f90=.o))
OBJ_WAVE :=  $(foreach file,$(SRC_WAVE),$(OBJDIR)/$(file:.f90=.o))
MKLPATH := /opt/intel/mkl/lib/intel64
LIBPATH := /opt/intel/lib/intel64/
LINKS := -lmkl_sequential -lmkl_core -lmkl_gf_lp64

#######################################################################
#                        GNU compiler                               #
#######################################################################
ifeq ($(COMPILER),GNU)
FC := mpifort
ifneq ($(DEBUG),ON)
  FFLAGS := $(OBJDIR) -I$(OBJDIR) -L${MKLPATH} ${LINKS} 
  FFLAGS := ${FFLAGS} -O3
else
  FFLAGS := -I$(OBJDIR) -L${MKLPATH} ${LINKS} 
  FFLAGS := ${FFLAGS} -O2 -fimplicit-none -Wall  -Wline-truncation  
  FFLAGS := ${FFLAGS} -Wcharacter-truncation -Wsurprising -Waliasing 
#  FFLAGS := ${FFLAGS} -Wimplicit-interface -std=f2008 
  FFLAGS := ${FFLAGS} -Wunused-parameter -fwhole-file -fcheck=all  
  FFLAGS := ${FFLAGS} -pedantic -fbacktrace
  FFLAGS := ${FFLAGS} -fallow-argument-mismatch 
endif
FLAGS := -c ${FFLAGS}

########################################################################
##                        Intel compiler                               #
########################################################################
#ifeq ($(COMPILER),Intel)
#FC := mpif90
#CPP := /usr/bin/cpp
#ifneq ($(DEBUG),ON)
#  FFLAGS := -f90=ifort -c -module $(OBJDIR) -I$(OBJDIR) -L${MKLPATH} ${LINKS} 
#  FFLAGS := -O3 -ftz -align records -no-prec-div -unroll -ipo $(FFLAGS) 
#else
#  FFLAGS := -BC -f90=ifort -g -c -module $(OBJDIR) -I$(OBJDIR) -L${MKLPATH} ${LINKS} 
#  FFLAGS := -O2 -assume realloc_lhs -check noarg_temp_created -traceback -warn all -fstack-protector -assume protect_parens -implicitnone ${FFLAGS} 
##  FFLAGS := -O2 -assume realloc_lhs -check all -traceback -warn all -fstack-protector -assume protect_parens -implicitnone ${FFLAGS} 
##  FFLAGS := -c -module $(OBJDIR) -I$(OBJDIR) -L${MKLPATH} ${LINKS} 
##  LDFLAGS := -L${LIBPATH} -lifcore -limf -lintlc
#endif
#ifeq ($(STATIC),ON)
#  FFLAGS := -static $(FFLAGS)
#endif
#  LDFLAGS := -f90=ifort $(LDFLAGS) -I$(OBJDIR) -L${MKLPATH} ${LINKS} 
#ifeq ($(VCTRIZE),ON)
#  FFLAGS := $(FFLAGS) -vec-report2
#  LDFLAGS := $(LDFLAGS) -vec-report2
#endif
#ifeq ($(WithOMP),ON)
#  FFLAGS := $(FFLAGS) -openmp -openmp_report2
#  LDFLAGS := $(LDFLAGS) -openmp -openmp_report2  
#endif

endif

#######################################################################
#                        subffixes rules                              #
#######################################################################

.SUFFIXES:
.SUFFIXES: .F90 .f90 .o

#%.f90 : %.F90 
#	$(CPP) $(CPPINCS) $(FPPFLAGS) $< > $(SRCDIR)$(@F)

$(OBJDIR)/%.o : %.f90
	$(FC) $(FLAGS) -o $(OBJDIR)/$(@F) $(SRCDIR)$(<F)

#######################################################################
#                        the target                                   #
#######################################################################
simu: $(BINDIR)/$(EXE_WAVE)

$(BINDIR)/$(EXE_WAVE): $(OBJ_MOD)  $(OBJ_WAVE) 
	$(FC) -o $@ $(OBJ_MOD) $(OBJ_WAVE) $(FFLAGS) 

rtm: $(BINDIR)/$(EXE_RTM)

$(BINDIR)/$(EXE_RTM): $(OBJ_MOD)  $(OBJ_RTM)
	$(FC) -o $@ $(OBJ_MOD) $(OBJ_RTM) $(LDFLAGS)

postproc:
	make -f ../postproc/makefile

RM := rm
cleanobj:
	$(RM) -f $(OBJDIR)/*
cleanout:
	$(RM) -f ../records/* ../snapshot/*
clean:
	$(RM) -f $(OBJDIR)/*  

