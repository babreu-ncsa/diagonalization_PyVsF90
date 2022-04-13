FC=ifort
FFLAGS=

MKLROOT=/sw/spack/delta-2022-03/apps/intel-oneapi-mkl/2022.0.2-gcc-8.4.1-pes6zb6/mkl/2022.0.2
MKLPATH=${MKLROOT}/lib/intel64
MKLINCUDE=${MKLROOT}/include
MKLLDFLAGS=-qmkl=parallel -L${MKLPATH} -lmkl_lapack95_lp64 -Wl,--start-group ${MKLPATH}/libmkl_intel_lp64.a ${MKLPATH}/libmkl_intel_thread.a ${MKLPATH}/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm

exec = timing_diag.exe

$(exec): timing_diag.o
	$(FC) -o $@ $^ ${MKLLDFLAGS}

timing_diag.o: timing_diag.f90
	$(FC) $(FFLAGS) -I${MKLINCLUDE} -I${MKLINCLUDE}/intel64/lp64 -c $^

clean:
	rm *.o $(exec)
