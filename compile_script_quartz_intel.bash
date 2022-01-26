#! /bin/bash
# Need to build SUNDIALS with -fPIC and set(CMAKE_MACOSX_RPATH 1)
export CANTERA_DIR=~/local/cantera/2.4.0/
export SUN_INCLUDE=~/local/sundials/2.7.0/include
export SUN_LIB=~/local/sundials/2.7.0/lib
export BOOST_INCLUDE=~/local/boost/1.77.0/build/include
export EIGEN_INCLUDE=~/local/eigen/3.2.9

source ~/envs/charlesx_py27/bin/activate
scons -j8 build prefix=$CANTERA_DIR \
  CXX=icpc CC=icc FORTRAN=ifort python_package=none \
  optimize_flags='-O3 -ip' \
  env_vars='all' \
  sundials_include=$SUN_INCLUDE sundials_libdir=$SUN_LIB \
  boost_inc_dir=$BOOST_INCLUDE f90_interface=y system_eigen=y \
  extra_inc_dirs=$EIGEN_INCLUDE
if [ -z "SCONS_TEST" ]; then
    scons -j8 test
fi
scons -j8 install
