
objects_astro =	./obj/af_mag_accretion.o \
		./obj/af_NS_atm_structure.o \
		./obj/af_orbit_bin.o \
		./obj/af_ref_frame.o \
		./obj/af_RT_layer_an.o

./obj/af_mag_accretion.o : ../AstroFun/af_mag_accretion.f90
	gfortran -c -o ./obj/af_mag_accretion.o ../AstroFun/af_mag_accretion.f90
./obj/af_NS_atm_structure.o : ../AstroFun/af_NS_atm_structure.f90
	gfortran -c -o ./obj/af_NS_atm_structure.o ../AstroFun/af_NS_atm_structure.f90
./obj/af_orbit_bin.o : ../AstroFun/af_orbit_bin.f90
	gfortran -c -o ./obj/af_orbit_bin.o ../AstroFun/af_orbit_bin.f90
./obj/af_ref_frame.o : ../AstroFun/af_ref_frame.f90
	gfortran -c -o ./obj/af_ref_frame.o ../AstroFun/af_ref_frame.f90
./obj/af_RT_layer_an.o : ../AstroFun/af_RT_layer_an.f90
	gfortran -c -o ./obj/af_RT_layer_an.o ../AstroFun/af_RT_layer_an.f90
