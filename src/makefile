target = ../pylib.so
object = \
  linearInterp.f90 \
  barycentric__interpolator.f90 \
  interpolate__fromTriMesh.f90 \
  cubicinterpolation_3d.f90 \
  interpolate__trilinear.f90

FC = gfortran
FFLAGS  = -O3 -llapack -lblas
LDFLAGS =

#FC = ifort
# FFLAGS  = -O3 -mkl
# LDFLAGS =

PYFLAGS = -fPIC -shared 

.SUFFIXES: .so .f90
$(target): $(object)
	$(FC) -o $(target) $(LDFLAGS) $(PYFLAGS) $(FFLAGS) $(object)

.PHONY:all clean
clean:
	$(RM) $(target) 
all: clean $(target)
