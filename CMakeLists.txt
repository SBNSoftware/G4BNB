set(CMAKE_C_COMPILER "$ENV{GCC_FQ_DIR}/bin/gcc")
set(CMAKE_CXX_COMPILER "$ENV{GCC_FQ_DIR}/bin/g++")
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++1y -g -rdynamic " )
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} $ENV{G4LIB} $ENV{ROOTSYS}/etc/cmake)

cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(NuBeam)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#

# Check G4 version - need to do some patching for G4.11+
set(G4_VERSION $ENV{GEANT4_VERSION})
string(REPLACE "_" ";" G4_VERSION_LIST ${G4_VERSION})
list(GET G4_VERSION_LIST 1 G4_VERSION_MINOR)

if( G4_VERSION_MINOR GREATER 10 )
  # Configure PTL as G4.11.2 needs it
  set(PTL_INCLUDE_DIR $ENV{GEANT4_FQ_DIR}"/include/Geant4")
  if( NOT TARGET PTL::ptl )
    add_library(PTL::ptl INTERFACE IMPORTED)
    set_target_properties(PTL::ptl PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${PTL_INCLUDE_DIR}"
    )
  endif()

  if( TARGET PTL::ptl )
    message(STATUS "PTL::ptl target OK")
  else()
    message(FATAL_ERROR "PTL::ptl target not found!")
  endif()

endif() # for G4.11+

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
