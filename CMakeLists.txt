cmake_minimum_required(VERSION 3.21)
project(finite_volumes_f90 Fortran)


enable_language(Fortran)
#add_executable(finite_volumes_f90 main.f90)
include_directories(${PROJECT_SOURCE_DIR})
add_executable(finite_volumes_f90 ch2_main.f90 ch2_general_vars.f90 ch2_subroutines.f90)