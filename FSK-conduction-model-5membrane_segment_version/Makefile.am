#INCLUDES = syspara.h

bin_PROGRAMS = 1D_FSK5S_LP
#1D_FSK5S_LP_CFLAGS = -parallel -openmp -par-report2 -xsse4.2 -pg -inline-level=2
1D_FSK5S_LP_CFLAGS = -O2 -axsse4.2 -mkl -parallel -par-report2 -Winline -inline-level=2 -ip -vec-report2 -opt-mem-layout-trans=3 -pg

#1D_FSK5S_LP_LDADD = -L/opt/intel/mkl/lib/intel64 -liomp5 -lmkl_intel_lp64 -lmkl_solver_lp64 -lmkl_intel_thread -lmkl_core
1D_FSK5S_LP_LDADD = -L/opt/intel/mkl/lib/intel64 -liomp5 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core

1D_FSK5S_LP_SOURCES = ./src/fsk.c \
	./src/exptable.c ./src/linear_mat.c ./src/pardiso_mat.c\
	./src/eular.c ./src/input_para.c ./src/system.c ./src/data_out.c\
	./src/mem.c ./src/var_consts.c ./src/syspara.h

