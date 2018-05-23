lib = lib/libde.a
OBJS = de_settings.o de_tools.o de_types.o de_hde.o de_wcdm3.o de_model_init.o de_chisqs_JLA.o de_chisqs.o de_srom.o de_ICG.o de_qz.o de_mauricehde.o #de_chisqs_AP.o
#OBJS = de_settings.o de_tools.o de_types.o de_hde.o de_wcdm3.o de_model_init.o de_srom.o de_ICG.o de_qz.o de_mauricehde.o

F90C = ifort #gfortran #ifort
F90FLAGS = -heap-arrys #-mkl_lp64th

default: $(lib)

de_tools.o : de_settings.o
de_types.o: de_tools.o
de_qz.o de_hde.o de_wcdm3.o de_srom.o de_ICG.o de_mauricehde.o: de_types.o
de_model_init.o : de_hde.o  de_wcdm3.o de_srom.o de_ICG.o de_qz.o de_mauricehde.o
de_chisqs_JLA.o : de_model_init.o
de_chisqs.o: de_model_init.o de_chisqs_JLA.o
#de_chisqs_AP.o: de_model_init.o
#chisqs.o : DEModel.o
#find_bf.o : chisqs.o
#main.o : DEModel.o chisqs.o find_bf.o

%.o: %.f90
	$(F90C) $(F90FLAGS) -c $*.f90 

$(lib): $(OBJS)
	ar cr $(lib) $(OBJS)
	cp *.mod mods/
#	echo $(str)
		
clean :
	rm -rf *.o *.mod *.a *.so mods/*.mod lib/*.a
	
