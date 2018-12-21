#=======================================================================
#                   define the compiler names
#=======================================================================

CC       = gcc
F90      = gfortran
#F90      = ifort
PYTHON   = python

#=======================================================================
#                     additional flags
#=======================================================================

ifeq ($(F90),gfortran)
	FPP      = gfortran -E -DPYEXT
	FPP_F90FLAGS = -x f95-cpp-input -fPIC -DPYEXT
	F90FLAGS = -fPIC -cpp -ffree-line-length-none -DPYEXT
    FCOMP    = gfortran
    LIBS     =
endif

ifeq ($(F90),ifort)

	FPP      = gfortran -E # gfortran f90wrap temp files only. not compilation
	FPP_F90FLAGS = -x f95-cpp-input -fPIC
	F90FLAGS = -fpscomp logicals -fPIC # use 1 and 0 for True and False
    FCOMP    = intelem # for f2py
    LIBS =
endif

CFLAGS = -fPIC  #     ==> universal for ifort, gfortran, pgi

#=======================================================================
#=======================================================================

UNAME = $(shell uname)

ifeq (${UNAME}, Darwin)
  LIBTOOL = libtool -static -o
else
  LIBTOOL = ar src
endif

# ======================================================================
# PROJECT CONFIG, do not put spaced behind the variables
# ======================================================================
# Python module name
PYTHON_MODN = CompDam_DGD
# mapping between Fortran and C types
KIND_MAP = kind_map

#=======================================================================
#       Relative path to the .for source files (needs trailing '/')
#=======================================================================

PATH_TO_FOR = ../for/

#=======================================================================
#       List all source files required for the project
#=======================================================================

# names (without suffix)
LIBSRC_SOURCES = forlog matrixUtil matProp schapery stress stateVar parameters strain friction cohesive fiberDamage schaefer plasticity DGD

# file names
LIBSRC_FILES = $(addsuffix .f90,${LIBSRC_SOURCES})

# object files
LIBSRC_OBJECTS = $(addsuffix .o,${LIBSRC_SOURCES})

# only used when cleaning up
LIBSRC_FPP_FILES = $(addsuffix .fpp,${LIBSRC_SOURCES})

#=======================================================================
#       List all source files that require a Python interface
#=======================================================================

# names (without suffix)
LIBSRC_WRAP_SOURCES = matrixUtil matProp stress stateVar parameters strain plasticity DGD

# file names
LIBSRC_WRAP_FILES = $(addsuffix .f90,${LIBSRC_WRAP_SOURCES})

# object files
LIBSRC_WRAP_OBJECTS = $(addsuffix .o,${LIBSRC_WRAP_SOURCES})

# fpp files
LIBSRC_WRAP_FPP_FILES = $(addsuffix .fpp,${LIBSRC_WRAP_SOURCES})

#=======================================================================
#                 Relevant suffixes
#=======================================================================

.SUFFIXES: .f90 .fpp

#=======================================================================
#
#=======================================================================

.PHONY: all clean


all: _${PYTHON_MODN}.so


clean:
	-rm ${LIBSRC_OBJECTS} \
		${LIBSRC_FPP_FILES} \
		libsrc.a \
		_${PYTHON_MODN}.so \
		_${PYTHON_MODN}_pkg.so \
		*.mod *.fpp *.f90 f90wrap*.o *.o \
		.f2py_f2cmap ${PYTHON_MODN}.py ${PYTHON_MODN}.pyc \
		-rf src.*


.f90.o: $(LIBSRC_FILES)
	${F90} ${F90FLAGS} -c $< -o $@


.c.o:
	${CC} ${CFLAGS} -c $< -o $@


.f90.fpp:
	${FPP} ${FPP_F90FLAGS} $<  -o $@


$(LIBSRC_FILES):
	-ln -s $(PATH_TO_FOR)$(basename $@).for $@


libsrc.a: ${LIBSRC_OBJECTS}
	${LIBTOOL} $@ $?


_${PYTHON_MODN}.so: libsrc.a ${LIBSRC_FPP_FILES}
	f90wrap -m ${PYTHON_MODN} ${LIBSRC_WRAP_FPP_FILES} -k ${KIND_MAP}
	f2py-f90wrap --fcompiler=$(FCOMP) --build-dir . -c -m _${PYTHON_MODN} -L. -lsrc f90wrap*.f90


_${PYTHON_MODN}_pkg.so: libsrc.a ${LIBSRC_FPP_FILES}
	f90wrap -m ${PYTHON_MODN}_pkg ${LIBSRC_WRAP_FPP_FILES} -k ${KIND_MAP} -P
	f2py-f90wrap --fcompiler=$(FCOMP) --build-dir . -c -m _${PYTHON_MODN}_pkg -L. -lsrc f90wrap*.f90


test: _${PYTHON_MODN}.so _${PYTHON_MODN}_pkg.so
	${PYTHON} tests.py
