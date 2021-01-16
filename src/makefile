target  = ../main
object  = \
  variablesMod.o \
  PBiCGStabMod.o \
  xyzPois3DMod.o \
  ioUtilityMod.o \
  main.o

# FC      = gfortran
FC      = mpif90

FFLAGS  = -O3 -llapack -lopenblas -lpthread
LDFLAGS = -L/usr/local/Cellar/openblas/0.3.10_1/lib/


.SUFFIXES: .o .f90

%.o: %.f90
	$(FC) -c $< $(FFLAGS) $(LDFLAGS)

%.mod: %.f90 %.o
	@:

$(target): $(object)
	$(FC) -o $@ $(object) $(FFLAGS) $(LDFLAGS)


.PHONY:clean all

clean:
	$(RM) $(target) $(object) *.mod

all: clean $(target)
