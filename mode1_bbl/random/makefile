SHELL=bash

# Kris R. -- Add option for output verbosity control at some point
DIM         = 3      #default is 3
DEBUG       = TRUE  #default is false
OPT         = FALSE   #default is false, can be TRUE, or HIGH
MPI         = TRUE   #default is false 
OPENMP      = FALSE #default is false 
FFTW        = FALSE  #default is false
NETCDF      = FALSE  #default is false
ADJUSTMENT  = FALSE #default is false

VPATH = ../src
SRCDIR = ../src
EXECDIR = ../exec

SRCS=   physicalconstants.f random.f dealiasing.f \
        create_MPI_data_type.f fft.f setup.f \
        filtering.f init.f\
        horfft.f dataflip.f various.f multistep.f nlterm.f derv.f\
        send_receive_nbh.f penalty.f testhelm.f mpi_setup.f\
        calcuvwt.f output.f restart.f postp.f lhs_gmres.f\
        solver_gmres.f band_lu.f bandludcmp.f ludcmp.f curvefit.f\
        regrid.f ratint.f vertgrad_smooth.f output_Dx.f vertfft.f\
        main.f

#Add special bridge to C++ routine object file to call SuperLU
#This is accessed through a symbolic link to
#the EXAMPLES subdirectory of SuperLU

OBJS=   physicalconstants.o random.o dealiasing.o \
        create_MPI_data_type.o fft.o setup.o \
        filtering.o init.o\
        horfft.o dataflip.o various.o multistep.o nlterm.o derv.o\
        send_receive_nbh.o penalty.o testhelm.o mpi_setup.o\
        calcuvwt.o output.o restart.o postp.o\
        lhs_gmres.o solver_gmres.o band_lu.o bandludcmp.o ludcmp.o curvefit.o\
        regrid.o ratint.o vertgrad_smooth.o output_Dx.o vertfft.o\
         main.o

MODS = physicalconstants.mod rng.mod dealias.mod 


NETCDF_FORTRAN_HOME = /opt/netcdf-fortran-4.4.3.mpich/netcdf-fortran
NETCDF_C_HOME  = /opt/netcdf-4.4.0.mpich/netcdf-4.4.0
HDF5_HOME = /opt/hdf5-1.8.13.mpich

# Tunable parameters
#
# CF		Name of the fortran compiling system to use
# LDFLAGS	Flags to the loader
# LIBS		List of libraries
# CMD		Name of the executable
# PROFLIB	Library needed for profiling
#
# On T90: Use CF=f90
# On ferret (USC machine): Use CF=f77
# FC = mpif90 #$(MPI_HOME)/bin/mpif90 		 # -L/usc/usr/lib
# CF = mpif90 #$(MPI_HOME)/bin/mpif90
FC = mpif90 #-m64 #SMJ
CF = mpif90 #-m64 #SMJ
LDFLAGS =    
PROF =

export $(LDFLAGS)

#export MPI

#LIBS = -L/usr/usc/absoft/default/lib -lV77      
#
#(PD: 7/25/03)-ABOUT THE LIBRARIES:
# 1) Use SuperLU and relevant CBLAS libraries to compile
#    Super-LU fast linear system solver. See SuperLU directory.
#    I copied the specific libraries into this directory for now.
# 2) The Atlas Lib libraries (may actually need libcblas
#    eventually) are used for matrix-matrix multiplication
#    in spectral differentiation & filtering. See ATLAS directory
# 3) FFTW libraries are used for implementation of FFTW. See
#    FFTW directory.
#
LIBS = # -L/usr/usc/absoft/default/lib -lV77
CMD = smpm.x     
PROFLIB = #	-lprof

# To perform the default compilation, use the first line
# To compile with flowtracing turned on, use the second line
# To compile giving profile additonal information, use the third line
# WARNING:  SIMULTANEOUSLY PROFILING AND FLOWTRACING IS NOT RECOMMENDED 
#All these commands must be commented out when compiling and running
#on the T90
#When compiling of ferret use first line
#
# Listed below are the compilation options for Linux-Absoft Fortran
# -f: fold all symbolic names to lower case.
# -N15: Appends underscore to names
# -N113: Set double precision as default
# -03: Advanced optimization level
# -C: Check array bounds
# -s: Force all program storage to be treated as static and initialized



# -w option quiets all warnings for now, but should double check these at
# some point -- Kris R.       
FFLAGS = -w -cpp -ffixed-line-length-none -fdefault-real-8 \
         -I$(SRCDIR) -I$(EXECDIR)\
         $(subst FALSE,-O0 ,$(subst TRUE, -O2,$(subst HIGH, -O3 -flto, $(OPT)) )) \
         $(subst FALSE,,$(subst TRUE, -g3 -ggdb -fbounds-check,$(DEBUG) )) \
         $(subst FALSE,,$(subst TRUE, -DUSE_MPI,$(MPI) )) \
         $(subst FALSE,,$(subst TRUE, -DADJUSTMENT,$(ADJUSTMENT) )) \
         $(subst FALSE,,$(subst TRUE, -DUSE_NETCDF -I$(NETCDF_FORTRAN_HOME)/include -I$(HDF5_HOME)/include,$(NETCDF) )) \
         $(subst FALSE,,$(subst TRUE, -fopenmp,$(OPENMP) )) \

LIBS = $(subst FALSE, ,$(subst TRUE, -L$(NETCDF_FORTRAN_HOME)/lib -L$(NETCDF_C_HOME)/lib -L$(HDF5_HOME)/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -ldl -lm, $(NETCDF) ))
 
#FFLAGS =    -f  -s -N15 -N113 -O3  #-r8  -i4  -fast   
#FFLAGS =	 -F
#FFLAGS =	 -Wf"-ez"

# Lines from here on down should not need to be changed.  They are the
# actual rules which make uses to build a.out.
# Note (PD: 7/03/03) that SuperLU library and Blas are needed to
# run Super LU. These were compiled in a separate directory
#
all:		$(CMD)

$(CMD):		$(OBJS) 
#$(SUPERLULIB) $(BLASLIB)
	$(CF) $(LDFLAGS) -o $(@) $(LIBS) $(PROF) $(OBJS)  
#$(SUPERLULIB) $(BLASLIB)

# Make the profiled version of the command and call it a.out.prof
#
$(CMD).prof:	$(OBJS)
	$(CF) $(LDFLAGS) -o $(@) $(OBJS) $(PROFLIB) $(LIBS)

clean:
	-rm -f $(OBJS)

realclean:	clean
	-rm -f $(CMD) $(CMD).prof $(MODS) *.nc fort.* *.x

run: ; mpirun -n 1 ./smpm.x

#void:	realclean
#	-rm -f $(SRCS) makefile
