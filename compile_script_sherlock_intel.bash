#! /bin/bash
# Need to build SUNDIALS with -fPIC and set(CMAKE_MACOSX_RPATH 1)
export CANTERA_DIR=~/local/cantera/2.4.0/
export SUN_INCLUDE=~/local/sundials/2.7.0/include
export SUN_LIB=~/local/sundials/2.7.0/lib
export BOOST_INCLUDE=~/tools/boost_1_55_0/installationDir/include
export EIGEN_INCLUDE=~/local/eigen/3.2.9/eigen/

ml py-numpy/1.14.3_py27

scons -j4 build prefix=$CANTERA_DIR \
  CXX=mpiicpc CC=mpiicc FORTRAN=ifort python_package=none python3_package='n' \
  optimize_flags='-O3' \
  env_vars='all' \
  sundials_include=$SUN_INCLUDE sundials_libdir=$SUN_LIB \
  boost_inc_dir=$BOOST_INCLUDE f90_interface=y system_eigen=n \
  extra_inc_dirs=$EIGEN_INCLUDE
if [ -z "SCONS_TEST" ]; then
    scons -j4 test
fi
scons -j4 install
