#!/bin/bash
msis_cmake_dir=msis2.1/src

pushd $msis_cmake_dir
cmake -B build && cmake --build build
popd

ln -fs ${msis_cmake_dir}/int_driver .
ln -fs ${msis_cmake_dir}/msis21.parm .
