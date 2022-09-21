#! /bin/bash
# Need to build SUNDIALS with -fPIC and set(CMAKE_MACOSX_RPATH 1)
export CANTERA_DIR=/p/home/mbonanni/local/cantera/2.4.0
export SUN_INCLUDE=/p/home/mbonanni/local/sundials/2.7.0/include
export SUN_LIB=/p/home/mbonanni/local/sundials/2.7.0/lib
export BOOST_DIR=/p/home/mbonanni/local/boost/1.55.0/include
export EIGEN_INCLUDE=/p/home/mbonanni/local/eigen/3.2.9

export NUMPY_INCLUDE=/p/home/mbonanni/local/miniconda3/envs/charlesx/lib/python2.7/site-packages/numpy/core/include

scons -j44 build prefix=$CANTERA_DIR \
  CXX=g++ CC=gcc FORTRAN=gfortran python_package=full \
  optimize_flags='-O3 -march=broadwell' \
  env_vars='all' \
  sundials_include=$SUN_INCLUDE sundials_libdir=$SUN_LIB \
  boost_inc_dir=$BOOST_DIR f90_interface=y system_eigen=n extra_inc_dirs=$EIGEN_INCLUDE:$NUMPY_INCLUDE
if [ -z "SCONS_TEST" ]; then
    scons -j44 test
fi
scons -j44 install
