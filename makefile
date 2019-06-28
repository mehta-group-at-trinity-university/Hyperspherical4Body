#LINUX
#
#MKL := /opt/intel.190611/mkl/lib/intel64
#1DSolver.x:	1DSolver.o
#	ifort 1DSolver.o  -Wl,--start-group $(MKL)/libmkl_intel_lp64.a $(MKL)/libmkl_sequential.a $(MKL)/libmkl_core.a $(MKL)/libmkl_blas95_ilp64.a $(MKL)/libmkl_lapack95_ilp64.a -Wl,--end-group -lpthread -lm -parallel -o Solver.x

1DSolver.x:	1DSolver.f90
	ifort -O4 -o Solver.x 1DSolver.f90 -mkl 
#-check arg_temp_created -g -traceback -warn all -check all -heap-arrays -debug extended
1DSolver.o:	1DSolver.f90 
	ifort 1DSolver.f90 -mkl
#Solver.x:     1DSolver.f90
#	ifort -openMP -parallel -O2 -C -check -check arg_temp_created -gen-interfaces -heap-arrays -debug inline-debug-info -g -traceback -check bounds -no-wrap-margin -i8 -I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include -O4 1DSolver.f90 -mkl ${MKLROOT}/lib/intel64/libmkl_blas95_ilp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_ilp64.a -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core  -lpthread -lm -ldl 

#Solver.x:     1DSolver.f90
#	ifort -mkl 1DSolver.f90
#ifort -O0 -qopenmp -g -traceback -debug full -gen-interface -heap-arrays -check bounds -mkl 1DSolver.f90

#WINDOWS
#
#1DSolver.exe:	1DSolver.f90 
#	ifort /O3 /fp:precise /fpscomp:ioformat /fpscomp:logicals  /debug:full /traceback /Qsave /Qzero /gen-interfaces /warn:interfaces /check /fpe0 1DSolver.f90 -Qmkl -link -libpath:c:/FEAST/lib/x64 libfeast_sparse.a

#1DSolver.exe:	1DSolver.f90 
#	ifort 1DSolver.f90 -Qmkl 
#-traceback -debug:full -check:bounds -gen-interfaces 
clean:
	rm *.x *.exe *.out *.obj *.un~ *.o *__genmod.mod *__genmod.f90 *.pdb
