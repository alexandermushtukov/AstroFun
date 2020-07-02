
objects_astro =	./obj/af_ref_frame.o \
		./obj/af_Rm.o \
		./obj/af_RT_layer_an.o



./obj/af_ref_frame.o : ../AstroFun/af_ref_frame.f90	
	gfortran -c -o ./obj/af_ref_frame.o ../AstroFun/af_ref_frame.f90
./obj/af_Rm.o : ../AstroFun/af_Rm.f90	
	gfortran -c -o ./obj/af_Rm.o ../AstroFun/af_Rm.f90
./obj/af_RT_layer_an.o : ../AstroFun/af_RT_layer_an.f90	
	gfortran -c -o ./obj/af_RT_layer_an.o ../AstroFun/af_RT_layer_an.f90
