#Last-modified: 13 Jun 2011 11:01:28 PM
SHELL = /bin/sh

#****************************************************
# @file  Makefile
#****************************************************
# Compilation Options
F90     = gfortran
LIBDIR  = /Users/danielpe/Dropbox/sdssidl_v170/src/atlas/ # /usr/lib/atlas-base/  # for my Debian laptop
#LIBDIR  = /usr/local/atlas/lib/ # for arjuna
#LIBDIR = 
OS     :=$(shell uname -s)
NPROCS :=1

# Directories
BINDIR = $(HOME)/bin
#OBJDIR = $(CURDIR)/../trash-dev
OBJDIR = $(CURDIR)
VPATH  = $(OBJDIR)

#****************************************************

#****************************************************
# System specific settings
ifeq ($(OS),Linux)
	NPROCS:=$(shell grep -c ^processor /proc/cpuinfo)
	LIB_FLAGS = 
endif

ifeq ($(OS),Darwin) # Assume Mac Os X
	NPROCS:=$(shell system_profiler | awk '/Number Of CPUs/{print $4}{next;}')
	LIB_FLAGS = -framework vecLib
endif

ifeq ($(NPROCS), 1)
	FLAG = -L$(LIBDIR) -llapack -lf77blas -lcblas -latlas $(LIB_FLAGS)
#    FLAG  = -llapack -lblas 
else
	FLAG = -L$(LIBDIR) -llapack -lpthread -lptf77blas -lptcblas -latlas $(LIB_FLAGS)
endif

#****************************************************


#****************************************************
# Mail Target
$(BINDIR)/spear : main.o description.o readfiles.o sort_data.o func.o fastlikelin.o getcmatrix.o jd2lag.o copy_modelpara.o run_amoeba.o amoeba.o cholu.o highres.o mcmc.o 
	$(F90) $(OBJDIR)/*.o  $(FLAG) -o $@
	cp -Rf  $(CURDIR)/scripts/SinWrapper.pl $(BINDIR)/SinWrapper.pl
	cp -Rf  $(CURDIR)/scripts/TopWrapper.pl $(BINDIR)/TopWrapper.pl
	cp -Rf  $(CURDIR)/scripts/GenLC.dat  $(BINDIR)/GenLC.pl
	cp -Rf  $(CURDIR)/scripts/FindCL.pl  $(BINDIR)/FindCL.pl
#****************************************************

#****************************************************
# OBJECT Compilation
$(OBJDIR)/amoeba.o : amoeba.f90
	$(F90) $(FLAG) -c $< -o $@

$(OBJDIR)/cholu.o : cholu.f90
	$(F90) $(FLAG) -c $< -o $@

$(OBJDIR)/sort_data.o : sort_data.f
	$(F90) $(FLAG) -c $< -o $@

$(OBJDIR)/description.o : description.f90
	$(F90) $(FLAG) -c $< -o $@

$(OBJDIR)/readfiles.o : readfiles.f90 description.o
	$(F90) $(FLAG) -c readfiles.f90  -o $@

$(OBJDIR)/copy_modelpara.o : copy_modelpara.f90 description.o
	$(F90) $(FLAG) -c copy_modelpara.f90 -o $@

$(OBJDIR)/jd2lag.o : jd2lag.f90 description.o
	$(F90) $(FLAG) -c jd2lag.f90 -o $@

$(OBJDIR)/getcmatrix.o : getcmatrix.f90 description.o jd2lag.o
	$(F90) $(FLAG) -c getcmatrix.f90 -o $@

$(OBJDIR)/highres.o :: highres.f90  description.o getcmatrix.o cholu.o jd2lag.o
	$(F90) $(FLAG) -c highres.f90 -o $@

$(OBJDIR)/fastlikelin.o : fastlikelin.f90 description.o getcmatrix.o cholu.o
	$(F90) $(FLAG) -c fastlikelin.f90 -o $@

$(OBJDIR)/func.o : func.f90 description.o copy_modelpara.o fastlikelin.o
	$(F90) $(FLAG) -c func.f90 -o $@

$(OBJDIR)/run_amoeba.o : run_amoeba.f90 description.o func.o amoeba.o
	$(F90) $(FLAG) -c run_amoeba.f90 -o $@

$(OBJDIR)/mcmc.o : mcmc.f90 description.o func.o
	$(F90) $(FLAG) -c mcmc.f90 -o $@

$(OBJDIR)/main.o : main.f90 description.o readfiles.o run_amoeba.o copy_modelpara.o highres.o mcmc.o
	$(F90) $(FLAG) -c main.f90 -o $@
#****************************************************

#****************************************************
# CLEAN
clean:
	rm -rf $(OBJDIR)/*.o $(CURDIR)/*.mod $(CURDIR)/fort.* $(CURDIR)/best3.new  $(BINDIR)/spear $(BINDIR)/SinWrapper.pl  $(BINDIR)/TopWrapper.pl $(BINDIR)/GenLC.pl $(BINDIR)/FindCL.pl

devclean:
	rm -rf $(OBJDIR)/*.o $(CURDIR)/*.mod $(CURDIR)/fort.* $(CURDIR)/best3.new 

#****************************************************



