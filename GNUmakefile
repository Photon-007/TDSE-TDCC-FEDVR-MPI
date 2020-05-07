MACHINE=$(shell uname -n)

#MODLOCFLAG = -module # the space at the end here is required and important!
#MODLOCFLAG = -J #use this one when compiling wit gfortran

ifeq ($(MACHINE),ASUS-UX550VE)
	GSL = /usr
#	FC  = mpif90
#	F77 = mpif90
#	CC  = mpicc

	FC  = h5pfc
	F77 = h5pfc
	CC  = h5pcc

	FCFLAGS = $(F77FLAGS) -ffree-line-length-none -fopenmp 
	F77FLAGS = -C -O3
	FPP = gfortran -cpp -E -P
#	FPPFLAGS +=-D__LANCZOS
	CXXFLAGS += -O3
	MODLOCFLAG = -J 
	LDFLAGS+= -llapack -lblas -lgsl -lgslcblas 
	ifeq ($(DEBUG),yes)
#		FCFLAGS += -C -traceback -g -debug all -ftrapuv
		FCFLAGS += -O0 -g -Wall -Wno-tabs -Wextra -fcheck=all -fbacktrace #-fimplicit-none
#		FCFLAGS += -traceback -g -debug all -ftrapuv
	endif
endif


ifeq ($(MACHINE),ws4315-z640-l)
	GSL = /usr
#	LIBS = /usr/lib
#	LIBLAPACK = /usr/lib/lapack
#	LIBBLAS = /usr/lib/libblas
	FC = mpif90
	F77 = mpif90
	CC = mpicc

	FCFLAGS = $(F77FLAGS) -ffree-line-length-none -fopenmp #-I$(GSL)/include
	#F77FLAGS = -Ofast -C
	F77FLAGS = -C -O3
	FPP = gfortran -cpp -E -P
	FPPFLAGS +=-D__LANCZOS
	CXXFLAGS += -O3
	MODLOCFLAG = -J 
#	LDFLAGS+= -llapack -lblas -lgsl -lgslcblas -L$(GSL)/lib64
#	LDFLAGS+= -L/usr -L$(LIBS) -L$(LIBLAPACK) -L$(LIBBLAS) -llapack -lblas -lgsl -lgslcblas #-L$(GSL)/lib64
	LDFLAGS+= -llapack -lblas -lgsl -lgslcblas #-L$(GSL)/lib64
endif


ifeq ($(MACHINE),photon00)
    MKL = /opt/intel/Compiler/11.1/072/mkl
    MKLIB = $(MKL)/lib/intel64
    MKINC = $(MKL)/include
    FC  = h5pfc
    F77 = h5pfc
    FPP = h5pfc
    CC  = h5pcc
	MODLOCFLAG = -module # the space at the end here is required and important!
    FCFLAGS = -O3 -ip -vec_report0 -openmp# -p -g
    F77FLAGS = $(FCFLAGS)
    LDFLAGS += -L$(MKLIB) -I$(MKINC) -lmkl_lapack -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lguide -lpthread -lgsl -lgslcblas
    LDFLAGS += -L/usr/local/lib
#   LDFLAGS += -llapack
endif


#FCFLAGS += -Iutils -Icore -Iio -Ireprezentation -Ioper -Iprop -Impi
FCFLAGS += -Iutils -Icore  -Ireprezentation -Ioper -Iprop 

DEFAULT : TDSE_MPI.exec
all: TDSE_MPI.exec 

echo :
	echo $(F77)

rmod:
	rm reprezentation/*.o reprezentation/*.mod
	rm core/*.o core/*.mod
	rm oper/*.o oper/*.mod
	rm prop/*.o prop/*.mod
	rm utils/*.o utils/*.mod
	rm *.o *.mod

clean :
	rm reprezentation/*.o reprezentation/*.mod
	rm core/*.o core/*.mod
	rm oper/*.o oper/*.mod
	rm prop/*.o prop/*.mod
	rm utils/*.o utils/*.mod
	rm *.exec *.o #*.mod

srctar : $(shell find . -name '*.f90' -or -name '*.F90' -or -name '*.c' -or -name '*.f' -or -name '*.py' -or -name '*.input' -or -name '*.bz2') GNUmakefile
	tar -cvf src_MPI.tar $^
	bzip2 -f src_MPI.tar

# .f files are Fortran77 and do not contain modules
%.o: %.f
	$(F77) $(F77FLAGS) -c $< $(MODLOCFLAG)`dirname $<` -o $@

%.o: %.f90
	$(FC) $(FCFLAGS) -c $< $(MODLOCFLAG)`dirname $<` -o $@

%.o: %.c
	$(CC) $(CCFLAGS) -c $< -o $@

# this tells make that it can make a .mod file by making the associated .o file
# $(NOOP) is not defined and so does nothing, but it's not an empty rule, which
# would mislead make
%.mod: %.o
	$(NOOP)

# this removes the implicit rule that wants to make a .o file from a .mod file with m2c
# (this is for modula files, which also use .mod as the suffix)
%.o: %.mod

%.f90: %.F90
	$(FPP) $(FPPFLAGS) $< -o $@

%.exec : %.o
	$(FC) $(FCFLAGS) $^ $(LDFLAGS) -o $@

#set the dependencies:

TDSE_MPI.exec : utils/nrtype.o   core/gvars.o   core/BoundStates.o   core/MPI_pars.o   core/Setup.o   reprezentation/FEDVR_module.o   reprezentation/TDCC_module.o   reprezentation/TDCC_FEDVR.o   oper/Operator_typedef.o   oper/Operator_math.o   oper/Operator_diag.o   oper/Operator_expect.o   oper/FEDVR_derivatives.o   oper/FEDVR_potentials.o   oper/LASER_potentials.o   utils/gridsub.o   utils/eigen_calc.o   utils/unit_converions.o   utils/specfun.o   utils/utils.o   prop/tProp.o   prop/Lanczos.o

TDSE_MPI.o :    utils/nrtype.mod core/gvars.mod core/BoundStates.mod core/MPI_pars.mod core/Setup.mod reprezentation/FEDVR_module.mod reprezentation/TDCC_module.mod reprezentation/TDCC_FEDVR.mod oper/Operator_typedef.mod oper/Operator_math.mod oper/Operator_diag.mod oper/Operator_expect.mod oper/FEDVR_derivatives.mod oper/FEDVR_potentials.mod oper/LASER_potentials.mod utils/gridsub.mod utils/eigen_calc.mod utils/unit_converions.mod utils/specfun.mod utils/utils.mod prop/tProp.mod prop/Lanczos.mod

core/gvars.o :    utils/nrtype.mod utils/unit_converions.mod core/MPI_pars.mod

core/MPI_pars.o : utils/nrtype.mod

core/BoundStates.o : utils/nrtype.mod core/gvars.mod reprezentation/FEDVR_module.mod

core/Setup.o :    utils/nrtype.mod core/gvars.mod core/BoundStates.mod core/MPI_pars.mod reprezentation/FEDVR_module.mod reprezentation/TDCC_module.mod reprezentation/TDCC_FEDVR.mod oper/Operator_typedef.mod oper/Operator_math.mod oper/FEDVR_derivatives.mod oper/FEDVR_potentials.mod oper/LASER_potentials.mod utils/gridsub.mod utils/eigen_calc.mod prop/tProp.mod

reprezentation/FEDVR_module.o : utils/nrtype.mod core/gvars.mod utils/gridsub.mod core/MPI_pars.mod

reprezentation/TDCC_module.o :  utils/nrtype.mod core/gvars.mod

reprezentation/TDCC_FEDVR.o :   utils/nrtype.mod core/gvars.mod reprezentation/FEDVR_module.mod reprezentation/TDCC_module.mod core/BoundStates.mod core/MPI_pars.mod

oper/Operator_typedef.o :  utils/nrtype.mod core/gvars.mod reprezentation/FEDVR_module.mod reprezentation/TDCC_module.mod

oper/Operator_math.o :     utils/nrtype.mod core/gvars.mod reprezentation/FEDVR_module.mod reprezentation/TDCC_module.mod oper/Operator_typedef.mod

oper/Operator_diag.o :     utils/nrtype.mod core/gvars.mod core/MPI_pars.mod reprezentation/FEDVR_module.mod reprezentation/TDCC_module.mod reprezentation/TDCC_FEDVR.mod oper/Operator_typedef.mod oper/Operator_math.mod oper/FEDVR_derivatives.mod oper/FEDVR_potentials.mod oper/LASER_potentials.mod utils/eigen_calc.mod

oper/Operator_expect.o :   utils/nrtype.mod core/gvars.mod core/MPI_pars.mod reprezentation/FEDVR_module.mod reprezentation/TDCC_module.mod reprezentation/TDCC_FEDVR.mod oper/Operator_typedef.mod oper/Operator_math.mod

oper/FEDVR_derivatives.o : utils/nrtype.mod core/gvars.mod reprezentation/FEDVR_module.mod oper/Operator_typedef.mod

oper/FEDVR_potentials.o :  utils/nrtype.mod core/gvars.mod reprezentation/FEDVR_module.mod oper/Operator_typedef.mod

oper/LASER_potentials.o :  utils/nrtype.mod core/gvars.mod reprezentation/TDCC_module.mod oper/Operator_typedef.mod 

prop/tProp.o   : utils/nrtype.mod core/gvars.mod prop/Lanczos.mod oper/LASER_potentials.mod oper/FEDVR_potentials.mod core/MPI_pars.mod

prop/Lanczos.o : utils/nrtype.mod core/gvars.mod reprezentation/FEDVR_module.mod reprezentation/TDCC_module.mod reprezentation/TDCC_FEDVR.mod oper/Operator_typedef.mod oper/Operator_math.mod utils/eigen_calc.mod core/MPI_pars.mod

utils/gridsub.o : utils/nrtype.mod

utils/eigen_calc.o : utils/nrtype.mod

utils/unit_converions.o : utils/nrtype.mod

utils/utils.o : utils/nrtype.mod

utils/specfun.o : utils/nrtype.mod utils/utils.mod



















