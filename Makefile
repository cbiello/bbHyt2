#-*- Makefile -*-
## Choose compiler: gfortran,ifort (g77 not supported, F90 constructs in use!)
COMPILER=gfortran
FC=$(COMPILER)
## Choose PDF: native,lhapdf
## LHAPDF package has to be installed separately
PDF=lhapdf
#Choose Analysis: dummy, process specific
## default analysis may require FASTJET package, that has to be installed separately (see below)
ANALYSIS=MiNNLO #basic analysis
# NOTE: the inline comment above makes $(ANALYSIS) end in a space, which
# makes every ifeq ("$(ANALYSIS)","...") below fail silently and leaves
# PWHGANAL empty (no analysis linked at all). Strip it.
ANALYSIS:=$(strip $(ANALYSIS))
#ANALYSIS=PCI #analysis for the isolation of power corrections with Rhorry
#ANALYSIS=flavantikt-bbgamgam #MiNNLO 4FS vs 5FS proj
#ANALYSIS=bbgamgam #HH background
## For static linking uncomment the following
#STATIC= -static
#

#######################################################################
#  Recola2 -- amplitude provider.  SET THESE TWO PATHS.
#
#  Recola2 must be built with the HEFT model file; see the README for
#  where to download it and the exact build commands.
#
#  RECOLALOCATION : the install prefix (contains lib/ and include/)
#  RECOLASRC      : the unpacked recola2-X.Y.Z source directory.
#                   Its include/ holds the internal modules (globals_rcl
#                   in particular) that the install prefix does not ship
#                   but that recola.f needs.
#
RECOLALOCATION = $(HOME)/Packages/recola-heft-install
RECOLASRC      = $(HOME)/Packages/recola2-collier-2.2.4/recola2-2.2.4
#
#######################################################################
RECOLAMODULEDIR=$(RECOLALOCATION)/include
RECOLALIBDIR=$(RECOLALOCATION)/lib
# static build -> librecola needs libmodelfile and libcollier after it
RECOLALIB= -L$(RECOLALIBDIR) -lrecola -lmodelfile -lcollier -lgfortran -lquadmath
RECFLAGS=-I$(RECOLAMODULEDIR) -I$(RECOLASRC)/include

STDCLIB=-lstdc++

OBJ=obj-$(COMPILER)
OBJDIR:=$(OBJ)

ifeq ("$(COMPILER)","gfortran")
F77=gfortran -ffixed-line-length-none -ffree-line-length-none -J$(OBJ) -I$(OBJ) -fbounds-check -fno-align-commons -fdec-static -fno-automatic 
#F77=gfortran -mcmodel=medium -ffixed-line-length-none -ffree-line-length-none -J$(OBJ) -I$(OBJ) -fbounds-check -fno-align-commons -fdec-static -fno-automatic
## -fbounds-check sometimes causes a weird error due to non-lazy evaluation
## of boolean in gfortran.
#FFLAGS= -Wall -Wimplicit-interface -fbounds-check
## For floating point exception trapping  uncomment the following
#FPE=-ffpe-trap=invalid,zero,overflow,underflow
## gfortran 4.4.1 optimized with -O3 yields erroneous results
## Use -O2 to be on the safe side
OPT=-O2
## For debugging uncomment the following
#DEBUG= -ggdb -ffpe-trap=invalid,zero,overflow
ifdef DEBUG
OPT=-O0
#FPE=-ffpe-trap=invalid,zero,overflow
#,underflow
endif
endif

CC = gcc
CXX = g++

ifeq ("$(COMPILER)","ifort")
F77 = ifort -save  -extend_source  -module $(OBJ)
#CXX = g++
#LIBS = -limf
#FFLAGS =  -checkm
## For floating point exception trapping  uncomment the following
#FPE = -fpe0
OPT = -O3 #-fast
## For debugging uncomment the following
#DEBUG= -debug -g
ifdef DEBUG
OPT=-O0
FPE = -fpe0
endif
endif

PWD=$(shell pwd)
WDNAME=$(shell basename $(PWD))
RESDIR=../
VPATH=./:./analysis-yb2:$(POWHEG-BOX_tmp):$(RESDIR):$(OBJDIR)

INCLUDE0=$(PWD)
INCLUDE1=$(RESDIR)/include
FF=$(F77) $(FFLAGS) $(FPE) $(OPT) $(DEBUG) -I$(INCLUDE0) -I$(INCLUDE1)

INCLUDE =$(wildcard $(RESDIR)/include/*.h *.h include/*.h)

ifeq ("$(PDF)","lhapdf")
# lhapdf-config is not always on PATH; fall back to the local install.
# Override from the command line with e.g.
#     make LHAPDF_CONFIG=/path/to/bin/lhapdf-config
LHAPDF_CONFIG?=$(shell which lhapdf-config 2>/dev/null)
ifeq ($(LHAPDF_CONFIG),)
LHAPDF_CONFIG=$(HOME)/Packages/lhapdf-install/bin/lhapdf-config
endif
PDFPACK=lhapdf6if.o lhapdf6ifcc.o
# NOTE: "lhapdf-config --cxxflags" returns only optimisation flags for
# LHAPDF >= 6.5 (it gives "-O3", no -I), so the include directory has to
# be taken from --incdir explicitly, otherwise lhapdf6ifcc.cc fails with
#     fatal error: LHAPDF/LHAPDF.h: No such file or directory
LHAPDFCXXFLAGS+= $(shell $(LHAPDF_CONFIG) --cxxflags) -I$(shell $(LHAPDF_CONFIG) --incdir)
LIBSLHAPDF= -Wl,-rpath,$(shell $(LHAPDF_CONFIG) --libdir)  -L$(shell $(LHAPDF_CONFIG) --libdir) -lLHAPDF
ifeq  ("$(STATIC)","-static")
## If LHAPDF has been compiled with gfortran and you want to link it statically, you have to include
## libgfortran as well. The same holds for libstdc++.
## One possible solution is to use fastjet, since $(shell $(FASTJET_CONFIG) --libs --plugins ) $(STDCLIB)
## does perform this inclusion. The path has to be set by the user.
# LIBGFORTRANPATH= #/usr/lib/gcc/x86_64-redhat-linux/4.1.2
# LIBSTDCPP=/lib64
LIBSLHAPDF+=  -L$(LIBGFORTRANPATH)  -lgfortranbegin -lgfortran -L$(LIBSTDCPP) $(STDCLIB)
endif
LIBS+=$(LIBSLHAPDF)
else
PDFPACK=mlmpdfif.o hvqpdfpho.o
endif

# fastjet-config is not always on PATH; fall back to the local install.
# Override with e.g.  make FASTJET_CONFIG=/path/to/bin/fastjet-config
FASTJET_CONFIG?=$(shell which fastjet-config 2>/dev/null)
ifeq ($(FASTJET_CONFIG),)
FASTJET_CONFIG=$(HOME)/Packages/fastjet-install/bin/fastjet-config
endif

ifeq ("$(ANALYSIS)","MiNNLO")
##To include Fastjet configuration uncomment the following lines.                                                                                                                                                    
#FASTJET_CONFIG=$(shell which fastjet-config)  # set globally above
#FASTJET_CONFIG=/u/asankar/Pkgs/fastjet/bin/fastjet-config
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) $(STDCLIB)
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=pwhg_bookhist-multi.o pwhg_analysis-minnlo.o
## Also add required Fastjet drivers to PWHGANAL (examples are reported)                                                                                                                                             
#PWHGANAL+= fastjetsisconewrap.o fastjetktwrap.o fastjetCDFMidPointwrap.o fastjetD0RunIIConewrap.o fastjetfortran.o                                                                                                  
#PWHGANAL+= fastjetsisconewrap.o fastjetktwrap.o fastjetfortran.o                                                                                                                                                    
PWHGANAL+=  fastjetfortran.o
else
ifeq ("$(ANALYSIS)","bjets")
##To include Fastjet configuration uncomment the following lines.                                                                                                                                                    
#FASTJET_CONFIG=$(shell which fastjet-config)  # set globally above
#FASTJET_CONFIG=/u/asankar/Pkgs/fastjet/bin/fastjet-config
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) $(STDCLIB)
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=pwhg_bookhist-multi.o pwhg_analysis-bjets.o
## Also add required Fastjet drivers to PWHGANAL (examples are reported)                                                                                                                                             
#PWHGANAL+= fastjetsisconewrap.o fastjetktwrap.o fastjetCDFMidPointwrap.o fastjetD0RunIIConewrap.o fastjetfortran.o                                                                                                  
#PWHGANAL+= fastjetsisconewrap.o fastjetktwrap.o fastjetfortran.o                                                                                                                                                    
PWHGANAL+=  fastjetfortran.o
else
ifeq ("$(ANALYSIS)","bbgamgam")
##To include Fastjet configuration uncomment the following lines.                                                                                                                                                    
FASTJET_CONFIG=/u/asankar/Pkgs/fastjet/bin/fastjet-config
##FASTJET_CONFIG=$(shell which fastjet-config)  # set globally above
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) $(STDCLIB)
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=pwhg_bookhist-multi.o pwhg_analysis-bbgamgam.o
## Also add required Fastjet drivers to PWHGANAL (examples are reported)                                                                                                                                             
#PWHGANAL+= fastjetsisconewrap.o fastjetktwrap.o fastjetCDFMidPointwrap.o fastjetD0RunIIConewrap.o fastjetfortran.o                                                                                                  
#PWHGANAL+= fastjetsisconewrap.o fastjetktwrap.o fastjetfortran.o                                                                                                                                                    
PWHGANAL+=  fastjetfortran.o
else
ifeq ("$(ANALYSIS)","flavantikt")
##To include Fastjet configuration uncomment the following lines.                                                                                                                                                    #FASTJET_CONFIG=$(shell which fastjet-config)  # set globally above
#FASTJET_CONFIG=$(shell which fastjet-config)  # set globally above
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) $(STDCLIB)
LIBSFASTJET += -L/u/asankar/POWHEG-BOX-RES/MiNNLOPS_res/bbH4FS_local/IFNPlugin -lIFNPlugin
FJCXXFLAGS += $(shell $(FASTJET_CONFIG) --cxxflags)
FJCXXFLAGS += -I/u/asankar/POWHEG-BOX-RES/MiNNLOPS_res/bbH4FS_local/IFNPlugin
#FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=pwhg_bookhist-multi.o pwhg_analysis-flavantikt.o
## Also add required Fastjet drivers to PWHGANAL (examples are reported)                                                                                                                                             
#PWHGANAL+= fastjetsisconewrap.o fastjetktwrap.o fastjetCDFMidPointwrap.o fastjetD0RunIIConewrap.o fastjetfortran.o                                                                                                  
#PWHGANAL+= fastjetsisconewrap.o fastjetktwrap.o fastjetfortran.o                                                                                                                                                    
PWHGANAL+=  fastjetfortranifn.o
else
ifeq ("$(ANALYSIS)","flavantikt-bbgamgam")
##To include Fastjet configuration uncomment the following lines.                                                                                                                                                    #FASTJET_CONFIG=$(shell which fastjet-config)  # set globally above
#FASTJET_CONFIG=$(shell which fastjet-config)  # set globally above
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) $(STDCLIB)
LIBSFASTJET += -L/u/asankar/POWHEG-BOX-RES/MiNNLOPS_res/bbH4FS_local/IFNPlugin -lIFNPlugin
FJCXXFLAGS += $(shell $(FASTJET_CONFIG) --cxxflags)
FJCXXFLAGS += -I/u/asankar/POWHEG-BOX-RES/MiNNLOPS_res/bbH4FS_local/IFNPlugin
#FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=pwhg_bookhist-multi.o pwhg_analysis-flavantikt-bbgamgam.o
## Also add required Fastjet drivers to PWHGANAL (examples are reported)                                                                                                                                             
#PWHGANAL+= fastjetsisconewrap.o fastjetktwrap.o fastjetCDFMidPointwrap.o fastjetD0RunIIConewrap.o fastjetfortran.o                                                                                                  
#PWHGANAL+= fastjetsisconewrap.o fastjetktwrap.o fastjetfortran.o                                                                                                                                                    
PWHGANAL+=  fastjetfortranifn.o

else
ifeq ("$(ANALYSIS)","PCI")

#FASTJET_CONFIG=$(shell which fastjet-config)  # set globally above
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) $(STDCLIB)
FJCXXFLAGS += $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=pwhg_bookhist-multi.o pwhg_analysis-PCI.o
PWHGANAL+=  fastjetfortran.o

endif
endif
endif
endif
endif
endif

LIBS+=-lz

%.o: %.f $(INCLUDE) | $(OBJDIR)
	$(FF) -c -o $(OBJ)/$@ $<

%.o: %.f90 $(INCLUDE) | $(OBJDIR)
	$(FF) -c -o $(OBJ)/$@ $<

%.o: %.F90 $(INCLUDE) | $(OBJDIR)
	$(FF) -c -o $(OBJ)/$@ $<

%.o: %.c | $(OBJDIR)
	$(CC) $(DEBUG) -c -o $(OBJ)/$@ $^

%.o: %.cc | $(OBJDIR)
	$(CXX) $(DEBUG) -c -o $(OBJ)/$@ $^ $(FJCXXFLAGS) $(LHAPDFCXXFLAGS) 

%.o: %.cpp | $(OBJDIR)
	$(CPP) $(DEBUG) -c -o $(OBJ)/$@ $^ 

USER=init_couplings.o init_processes.o Born_phsp.o Born.o virtual.o	\
     real.o Check_LesHouches.o $(PWHGANAL)


ifdef RECOLALOCATION
FF+= $(RECFLAGS)
# recola.f is kept LOCAL to this process directory (as in
# bblnulnu_recola): it carries the generalisation of the coupling-power
# bookkeeping to res_powst > 1, which the copy in RecolaStuff does not
# have. VPATH puts ./ first, so the local one is picked up.
USER+=recola.o
endif



# PYTHIA 8

FJCXXFLAGS+=$(shell  pythia8-config --cxxflags)
#LIBPYTHIA8=$(shell pythia8-config --ldflags) -ldl $(STDCLIB) #-llhapdfdummy
LIBPYTHIA8= -L$(shell pythia8-config --libdir)  -L$(shell pythia8-config --libdir)/archive -lpythia8 -ldl $(STDCLIB) #-llhapdfdummy


PWHG=pwhg_main.o pwhg_init.o bbinit.o btilde.o lhefwrite.o		\
	LesHouches.o LesHouchesreg.o gen_Born_phsp.o find_regions.o	\
	fill_res_histories.o                                        	\
	test_Sudakov.o pt2maxreg.o sigborn.o gen_real_phsp.o maxrat.o	\
	gen_index.o gen_radiation.o Bornzerodamp.o sigremnants.o	\
	sigregular.o build_resonance_hists.o resize_arrays.o		\
	random.o boostrot.o bra_ket_subroutines.o cernroutines.o	\
	init_phys.o powheginput.o pdfcalls.o sigreal.o sigcollremn.o	\
	pwhg_analysis_driver.o checkmomzero.o		                \
	setstrongcoupl.o integrator.o mintwrapper.o newunit.o mwarn.o  	\
	sigsoftvirt.o reshufflemoms.o                                 	\
	sigcollsoft.o sigvirtual.o  ubprojections-new.o	            	\
        resweights.o locks.o genericphsp.o PhaseSpaceUtils.o boostrot4.o\
	setlocalscales.o mint_upb.o opencount.o         	        \
        lhefread.o pwhg_io_interface.o rwl_weightlists.o rwl_setup_param_weights.o \
	fullrwgt.o rwl_setup_param_weights_user.o sigequiv_hook.o	\
	validflav.o cache_similar.o utils.o  $(PDFPACK) $(USER) $(FPEOBJ)

LIBDIRMG=$(OBJ)

# target to generate LHEF output
pwhg_main:check_recola $(PWHG)
	$(FF) $(patsubst %,$(OBJ)/%,$(PWHG)) $(LIBS) $(RECOLALIB) $(LIBSFASTJET) $(STATIC) -o $@ $(STDCLIB)

LHEF=lhef_analysis.o boostrot.o random.o locks.o cernroutines.o utils.o \
	opencount.o powheginput.o $(PWHGANAL)   \
        lhefread.o pwhg_io_interface.o newunit.o pwhg_analysis_driver.o $(FPEOBJ) \
	rwl_weightlists.o

# target to analyze LHEF output
lhef_analysis:$(LHEF)
	$(FF) $(patsubst %,$(OBJ)/%,$(LHEF)) $(LIBS) $(LIBSFASTJET) $(STATIC)  -o $@



# target to read event file, shower events with HERWIG + analysis
HERWIG=main-HERWIG.o setup-HERWIG-lhef.o herwig.o boostrot.o	\
	powheginput.o $(PWHGANAL) lhefread.o	\
	pdfdummies.o opencount.o $(FPEOBJ)

main-HERWIG-lhef: $(HERWIG)
	$(FF) $(patsubst %,$(OBJ)/%,$(HERWIG))  $(LIBSFASTJET)  $(STATIC) -o $@

# target to read event file, shower events with PYTHIA + analysis
PYTHIA=main-PYTHIA.o setup-PYTHIA-lhef.o pythia.o boostrot.o powheginput.o \
	$(PWHGANAL) lhefread.o newunit.o 	\
	pwhg_analysis_driver.o random.o cernroutines.o opencount.o	\
	$(FPEOBJ)

main-PYTHIA-lhef: $(PYTHIA)
	$(FF) $(patsubst %,$(OBJ)/%,$(PYTHIA)) $(LIBS) $(LIBSFASTJET)  $(STATIC) -o $@



# target to read event file, shower events with PYTHIA8.1 + analysis
PYTHIA8=main-PYTHIA8.o boostrot.o powheginput.o locks.o \
	$(PWHGANAL) opencount.o lhefread.o pwhg_io_interface.o newunit.o pdfdummies.o \
	random.o cernroutines.o bra_ket_subroutines.o boostrot4.o utils.o\
	$(FPEOBJ) $(LIBZDUMMY) \
	rwl_weightlists.o

main-PYTHIA8-lhef: $(PYTHIA8) pythia8F77.o
	$(FF) $(patsubst %,$(OBJ)/%,$(PYTHIA8) pythia8F77.o ) $(LIBSFASTJET) $(LIBPYTHIA8) $(STATIC) $(LIBS) -o $@

main-PYTHIA8-lhef-gamgam: $(PYTHIA8) pythia8F77-gamgam.o
	$(FF) $(patsubst %,$(OBJ)/%,$(PYTHIA8) pythia8F77-gamgam.o ) $(LIBSFASTJET) $(LIBPYTHIA8) $(STATIC) $(LIBS) -o $@


# target to read event file, shower events with PYTHIA8.2 + analysis
main-PYTHIA82-lhef: $(PYTHIA8) pythia82F77.o
	$(FF) $(patsubst %,$(OBJ)/%,$(PYTHIA8) pythia82F77.o ) $(LIBSFASTJET) $(LIBPYTHIA8) $(STATIC) $(LIBS) -o $@

# Recola is not built by this Makefile: check that it is there, and that
# it was configured with the HEFT model file.
check_recola:
	@if [ ! -f "$(RECOLALIBDIR)/librecola.a" ] && [ ! -f "$(RECOLALIBDIR)/librecola.so" ]; then \
	  echo "ERROR: Recola library not found in $(RECOLALIBDIR)."; \
	  echo "       Build it with the HEFT model file:"; \
	  echo "       See the README; set RECOLALOCATION in this Makefile."; \
	  exit 1; \
	fi
	@if [ ! -f "$(RECOLAMODULEDIR)/recola.mod" ]; then \
	  echo "ERROR: recola.mod not found in $(RECOLAMODULEDIR)."; \
	  echo "       Check RECOLAMODULEDIR at the top of this Makefile."; \
	  exit 1; \
	fi

ifeq ("$(COMPILER)","gfortran")
XFFLAGS +=-ffixed-line-length-132
else
XFFLAGS +=-extend-source
endif

clean:
	rm -f $(patsubst %,$(OBJ)/%,$(USER) $(PWHG) $(LHEF) $(HERWIG) $(PYTHIA) $(PYTHIA8) pythia8?F77.o) \
        pwhg_main lhef_analysis main-HERWIG-lhef main-PYTHIA*-lhef


veryclean: clean
	rm -f $(OBJ)/*.o $(OBJ)/*.mod $(OBJ)/*.a $(OBJ)/*.so pwhg_main lhef_analysis main-HERWIG-lhef	\
	main-PYTHIA-lhef *.a MODEL/*.o


# target to generate object directory if it does not exist
$(OBJDIR):
	mkdir -p $(OBJDIR)

##########################################################################


# --- Recola module dependencies --------------------------------------
# every file that does "use recola_powheg" must be compiled after recola.o
init_processes.o: recola.o
Born.o: recola.o
virtual.o: recola.o
real.o: recola.o
