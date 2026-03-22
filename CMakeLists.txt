set(CMAKE_C_COMPILER "$ENV{GCC_FQ_DIR}/bin/gcc")
set(CMAKE_CXX_COMPILER "$ENV{GCC_FQ_DIR}/bin/g++")
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++1y -g -rdynamic " )
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} $ENV{G4LIB} $ENV{ROOTSYS}/etc/cmake)

cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()

project(NuBeam)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Setup root
#
find_package(ROOT REQUIRED)   

add_subdirectory(src)
add_subdirectory(scripts)
