gfortran -c mesh.f90 mod_func.f90
gfortran -c ./elem_proc/extend_mesh.f90
gfortran -c ./elem_proc/tripole_mod.f90
gfortran -c ./hi_intgl/hi_funcs.f90
gfortran -c ./hi_intgl/hi_integral.f90
gfortran -c ./elem_proc/ElemIntgl0.f90
