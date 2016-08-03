SHELL=/bin/sh


####################################################################################
# INTEL COMPILER
F90_SRC = 

F77_SRC = main.f setup.f fft.f testhelm.f various.f solver.f horfft.f\
        derv.f filtering.f init.f calcuvwt.f\
        output.f restart.f postp.f erfc.f sgscalc.f curvefit.f penalty.f\
        testderv.f band_lu.f bandludcmp.f ludcmp.f multistep.f bcreset.f\
        forcing.f

F90_OBJS = $(patsubst %.f90, %.o, $(F90_SRC))
F90_MODS = $(patsubst %.f90, %.mod, $(F90_SRC))

F77_OBJS =  $(patsubst %.f, %.o, $(F77_SRC))


OBJS = $(F90_OBJS) $(F77_OBJS)


# Intel Compiler
FC90 = ifort
FC77 = ifort

# Intel Compiler Flags
FFLAGS90 = 
FFLAGS77 = -cpp -c -nogen-interfaces -autodouble -wrap-margin-  -O3 

LIBS = -L/usr/local/lib -lm -mkl
CMD = mode1.exe

all: $(CMD)


# Link the executable. 
$(CMD):	$(OBJS) 
	$(FC90) $(LDFLAGS) -o $(@) $(PROF) $(OBJS) $(LIBS) 

%.o: %.f90 
	$(FC90) -o $@ $(FFLAGS90) $< 

%.o: %.f 
	$(FC77) -o $@ $(FFLAGS77) $< 
	

clean: 
	-rm -f $(OBJS) $(CMD)  *.mod




