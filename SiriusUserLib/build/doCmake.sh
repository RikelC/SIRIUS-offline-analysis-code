#!/bin/csh
rm -rf CMakeFiles
rm -f CMakeCache.txt cmake_install.cmake libUserLIB.a  libUserLIB.so Makefile
sleep 1
cmake -DCMAKE_INSTALL_PREFIX=../ ../
sleep 1
make
make install
