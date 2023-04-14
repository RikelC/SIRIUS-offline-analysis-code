-------------------------------------
Build and Installation of MFMlibrary
-------------------------------------
1. Download MFMLibrary
2. cd MFMlib-master
3. mkdir build 
4. cd build
     If tinyxml is installed
5. cmake -DCMAKE_INSTALL_PREFIX=../ ../sources/; make; make install
     else
5. cmake -DNO_MFMNXML=YES -DCMAKE_INSTALL_PREFIX=../ ../sources/

-----------------------
Build of UserLibrary
-----------------------
1. Download UserLib from Github: GUser-SIRIUS/GUser-sources/UserLib
2. cd UserLib
3. cd build
4. bash doCmake.sh

----------------------------------------------
Build of sirius_offline_data_analysis_program
----------------------------------------------
1. cd sirius_offline_data_analysis_program/program
2. open CMakeLists.txt
3. goto line 40. set(UserLib_Path path to /UserLib/build )
4. goto line 60. set(MFM_DIR path to /MFMlib-master/)
5. cd sirius_offline_data_analysis_program/build
6. bash doCmake.sh
