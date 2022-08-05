#! /bin/bash
# Need to build SUNDIALS with -fPIC and set(CMAKE_MACOSX_RPATH 1)

export CANTERA_DIR=$HOME/local/cantera/2.4.0
export SUN_INCLUDE=$HOME/local/sundials/2.7.0/include
export SUN_LIB=$HOME/local/sundials/2.7.0/lib
export BOOST_DIR=${HOME}/local/boost_1_46_1/include
export EIGEN_INCLUDE=$HOME/local/eigen/3.2.9/eigen  

scons -j4 build prefix=$CANTERA_DIR CXX=c++ CC=c python_package=full \
    sundials_include=$SUN_INCLUDE sundials_libdir=$SUN_LIB \
    boost_inc_dir=$BOOST_DIR f90_interface=n system_eigen=y extra_inc_dirs=$EIGEN_INCLUDE

# scons test
scons -j4 install
