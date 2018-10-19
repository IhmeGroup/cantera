#! /bin/bash
# Need to build SUNDIALS with -fPIC and set(CMAKE_MACOSX_RPATH 1)
export CANTERA_DIR=$HOME/numLibs/cantera/2.4.0
export SUN_INCLUDE=$HOME/numLibs/sundials/2.7.0/include
export SUN_LIB=$HOME/numLibs/sundials/2.7.0/lib/
export BOOST_DIR=$HOME/numLibs/boost_1_60_0

scons -j20 build prefix=$CANTERA_DIR \
  CXX=CC CC=cc FORTRAN=ftn python_package=full \
  optimize_flags='-O3' \
  thread_flags='' \
  warning_flags='-Wall' \
  env_vars='all' \
  sundials_include=$SUN_INCLUDE sundials_libdir=$SUN_LIB \
  boost_inc_dir=$BOOST_DIR f90_interface=n system_eigen=n
if [ -z "SCONS_TEST" ]; then
    scons -j20 test
fi
scons -j20 install
