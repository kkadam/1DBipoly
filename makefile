#compiler
FC=/opt/intel/composerxe-2011.5.209/bin/intel64/ifort

F90FILES= main.F90 poisson_solve.F90 findpar.F90

OFILES= $(F90FILES:.F90=.o) 

bipoly:$(OFILES)
	$(FC) $(OFILES) -o bipoly

$(OFILES):$(F90FILES)
	$(FC) -g -c -r8 -O3 $(F90FILES)
 
cl:
	rm -f *.o scf	
