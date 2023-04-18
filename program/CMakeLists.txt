#----------------------------------------------------------------------------
# Setup the project
#----------------------------------------------------------------------------


cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(SIRIUS)


include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)


if(COMPILER_SUPPORTS_CXX14)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++14 support. Please use a different C++ compiler.")
endif()

#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)


set(CMAKE_POSITION_INDEPENDENT_CODE ON) 


#----------------------------------------------------------------------------
# Tinyxml
#----------------------------------------------------------------------------

if(NO_MFMNXML)
message(STATUS "compilation without tinyxml")

add_definitions( -DNO_MFMXML )
else()
find_library(TinyXML_LIBRARIES tinyxml HINTS ${MYTINYXMLDIR} )
	if(TinyXML_LIBRARIES)
   	message(STATUS "-- found local ${TinyXML_LIBRARIES}")
	else()
	find_library(TinyXML_LIBRARIES tinyxml  )
		if (TinyXML_LIBRARIES)
		message(STATUS "-- found  ${TinyXML_LIBRARIES}")
		else()
   		message(WARNING "-- TinyXML_LIBRARIES not found")
   		endif()
	endif()
endif()


#----------------------------------------------------------------------------
#                                      ROOT
#----------------------------------------------------------------------------

find_package(ROOT REQUIRED COMPONENTS  MathMore MathCore)

include_directories(${ROOT_INCLUDE_DIRS})


#----------------------------------------------------------------------------
#         Locate sources and headers for this project
#----------------------------------------------------------------------------
include_directories(${PROJECT_SOURCE_DIR}/includes)
#include_directories(${CMAKE_BINARY_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/sources/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/includes/*.hh)

#----------------------------------------------------------------------------
#                      User Library
#----------------------------------------------------------------------------

set(UserLib_Path /Users/chakma/Desktop/SiriusTest/SIRIUS-offline-analysis-code-main/SiriusUserLib/lib )

find_library(User_LIBRARIES NAMES SiriusUserLIB  HINTS ${UserLib_Path})

if (User_LIBRARIES)
        message(STATUS "-- found  ${User_LIBRARIES}")
        else()
        message(FATAL_ERROR "-- User_LIBRARIES not found in  ${User_LIBRARIES}")
endif()

set (UserLibIncludes ${UserLib_Path}/../includes)
set (UserLibSources ${UserLib_Path}/../sources)

include_directories(${UserLibIncludes} ${UserLibSources})

#add_subdirectory(UserLib)

#----------------------------------------------------------------------------
#                          MFM library
#----------------------------------------------------------------------------

set(MFM_DIR /Users/chakma/Desktop/SiriusTest/MFMlib-master)
include_directories(${MFM_DIR}/include/)
include_directories(${MFM_DIR}/include/GanTape/)

find_library(MFM_LIBRARY NAMES MFM HINTS ${MFM_DIR}/lib )

#find the headers

find_path(MFM_INCLUDE_DIR NAMES MFMAllFrames.h HINTS ${MFM_DIR}/include/)

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

#----------------------------------------------------------------------------
# Add the executable, and link it to the ROOT, User, MFM and Tinyxml libraries
#----------------------------------------------------------------------------

add_executable(SiriusTest siriusTest.cc ${sources} ${headers})
target_link_libraries(SiriusTest ${User_LIBRARIES} ${ROOT_LIBRARIES} ${TinyXML_LIBRARIES})
target_link_libraries(SiriusTest libMFM::libMFM)

add_custom_target(SIRIUS DEPENDS SiriusTest)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#----------------------------------------------------------------------------

install(TARGETS SiriusTest DESTINATION bin)



