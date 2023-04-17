#!/bin/csh
rm -rf CMakeFiles
rm -f CMakeCache.txt cmake_install.cmake Makefile siriusTest
sleep 1
cmake -DNO_MFMNXML=YES -DCMAKE_INSTALL_PREFIX=../ ../program
sleep 1
make
