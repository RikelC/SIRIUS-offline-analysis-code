#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(SIRIUS)
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#----------------------------------------------------------------------------
# find root package
# Setup include directory for this project
#find_package(ROOT REQUIRED)
#
#SET( ROOT_USE_COMPONENTS MathMore Gdml Geo)
# or
find_package(ROOT REQUIRED COMPONENTS  MathMore MathCore)

include_directories(${ROOT_INCLUDE_DIRS})

include_directories(${PROJECT_SOURCE_DIR}/includes)
#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/sources/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/includes/*.hh)

#-----------------------
# User Library
#--------------------------
set(UserLib_Path /home/sirius/Chakma/Gru_dev/GRU_SIRIUS/GUser-sources/UserLib/build)
find_library(User_LIBRARIES NAMES UserLIB  HINTS ${UserLib_Path})

if (User_LIBRARIES)
        message(STATUS "-- found  ${User_LIBRARIES}")
        else()
        message(FATAL_ERROR "-- User_LIBRARIES not found in  ${User_LIBRARIES}")
endif()

set (UserLibIncludes ${UserLib_Path}/../includes)
set (UserLibSources ${UserLib_Path}/../sources)

include_directories(${UserLibIncludes} ${UserLibSources})

#add_subdirectory(UserLib)


#add MFM library
#find the library
#include(GNUInstallDirs)
#set(MFM_DIR MFMlib-master/sources/)
#include_directories(${CMAKE_BINARY_DIR})
include_directories(${CMAKE_SOURCE_DIR}/../MFMlib-master/sources/)
include_directories(${CMAKE_SOURCE_DIR}/../MFMlib-master/sources/GanTape/)

find_library(MFM_LIBRARY NAMES MFM HINTS ${PROJECT_BINARY_DIR}/../ PATH_SUFFIXES MFMlib-master/build/ )
#find the headers
find_path(MFM_INCLUDE_DIR NAMES MFMAllFrames.h HINTS ${PROJECT_BINARY_DIR}/../ ${CMAKE_INSTALL_INCLUDEDIR} PATH_SUFFIXES MFMlib-master/sources/)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(libMFM DEFAULT_MSG
                                  MFM_LIBRARY
                                  MFM_INCLUDE_DIR)
#make some checks
mark_as_advanced(MFM_LIBRARY MFM_INCLUDE_DIR)
if(LIBMFM_FOUND AND NOT TARGET libMFM::libMFM)#1st arg is set if find_package_handle_standard_args() succeeds and 2nd part NOT TARGET ensures that the library is not already declared as a target withing the current project
  add_library(libMFM::libMFM SHARED IMPORTED)#SHARED= dynamically linked library and IMPORTED= library is not built by the current project
set_target_properties(
    libMFM::libMFM
    PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${MFM_INCLUDE_DIR}"
      IMPORTED_LOCATION ${MFM_LIBRARY})
endif()
#list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMake/")
#find_package(libMFM REQUIRED)
#--------------
#or
#add_library(MFM_LIBRARIES STATIC IMPORTED)
#set_target_properties(MFM_LIBRARIES PROPERTIES IMPORTED_LOCATION /home/sirius/Chakma/sirius_offline_data_analysis_program/MFMlib-master/build/libMFM.a)
#set_target_properties(MFM_LIBRARIES PROPERTIES INTERFACE_INCLUDE_DIRECTORIES /home/sirius/Chakma/sirius_offline_data_analysis_program/MFMlib-master/sources/)
#----------------------------------------------------------------------------
# Add the executable, and link it to the ROOT libraries
#
add_executable(SiriusTest siriusTest.cc ${sources} ${headers})
target_link_libraries(SiriusTest ${User_LIBRARIES} ${ROOT_LIBRARIES} Spectrum)
target_link_libraries(SiriusTest libMFM::libMFM)

add_custom_target(SIRIUS DEPENDS SiriusTest)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS siriusTest DESTINATION bin)



