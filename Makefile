# For intel compiler:
FC = ifort   
# With intel byterecl should match VASP compilation
FLAGS =  -O -assume byterecl

# For gfortran:
#FC = gfortran
#FLAGS = -O

# It must be linked with a fftw implementation, crutial for performance
#LINK = /usr/lib/x86_64-linux-gnu/libfftw3.so.3

LINK = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl

OBJ = declaration.o declfft.o determinequantities.o  \
	volumen.o Fourier.o currentBardeen.o  \
	norm_name.o mappingscar_gen.o  \
	STMpw.o  

STMpw.out: $(OBJ)
	$(FC) $(FLAGS) $(OBJ) $(LINK) -o STMpw.out

utils: Imagen_Bardeen_gnu.out Cond_gnu.out

all: utils STMpw.out

Imagen_Bardeen_gnu.out: 
	$(FC) $(FLAGS) Utils/Imagen_Bardeen_gnu.f90 -o Utils/Imagen_Bardeen_gnu.out

Cond_gnu.out:
	$(FC) $(FLAGS) Utils/Cond_gnu.f90 -o Utils/Cond_gnu.out

Fourier.o: Fourier.f90
	$(FC) Fourier.f90 $(FLAGS) -c

declaration.o: declaration.f90
	$(FC) declaration.f90 $(FLAGS) -c

declfft.o: declfft.f90
	$(FC) declfft.f90 $(FLAGS) -c

determinequantities.o: determinequantities.f90
	$(FC) determinequantities.f90 $(FLAGS) -c

currentBardeen.o: currentBardeen.f90
	$(FC) currentBardeen.f90 $(FLAGS) -c

volumen.o: volumen.f90
	$(FC) volumen.f90 $(FLAGS) -c

STMpw.o: STMpw.f90
	$(FC) STMpw.f90 $(FLAGS) -c	

mappingscar_gen.o: mappingscar_gen.f90
	$(FC) mappingscar_gen.f90 $(FLAGS) -c	

norm_name.o: norm_name.f90
	$(FC) norm_name.f90 $(FLAGS) -c	

clean: 
	@echo "Cleaning objects..."
	rm -f *.o *.mod

cleanall: 
	@echo "Cleaning all..."
	rm -f *.o *.mod *.out */*.out
