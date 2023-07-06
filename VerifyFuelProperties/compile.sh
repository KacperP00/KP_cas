rm *.mod
ifort -g ../src/precision.f90 ../src/pc_defs.f90 ../src/pc_database.f90 ../src/pc_func.f90 fuel_properties.f90 -o fuel_properties
