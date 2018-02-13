#!/bin/bash
rm *.mod
ifort -g precision.f90 math.f90 pc_defs.f90 pc_database.f90 pc_func.f90 rk_defs.f90 solver_defs.f90 spray_defs.f90 rk_func.f90 solver_func.f90 spray_func.f90  parser.f90 cas.f90 -o cas

rm *.mod
ifort -O2 precision.f90 math.f90 pc_defs.f90 pc_database.f90 pc_func.f90 rk_defs.f90 solver_defs.f90 spray_defs.f90 rk_func.f90 solver_func.f90 spray_func.f90  parser.f90 cas.f90 -o cas_opt
