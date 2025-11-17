# Install script for directory: /exp/sbnd/app/users/fnicolas/BNBFluxUncertainties/G4BNB/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/exp/sbnd/app/users/fnicolas/BNBFluxUncertainties/G4BNB/bin/v4_10_4_p02d/NuBeam" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/exp/sbnd/app/users/fnicolas/BNBFluxUncertainties/G4BNB/bin/v4_10_4_p02d/NuBeam")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/exp/sbnd/app/users/fnicolas/BNBFluxUncertainties/G4BNB/bin/v4_10_4_p02d/NuBeam"
         RPATH "")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/exp/sbnd/app/users/fnicolas/BNBFluxUncertainties/G4BNB/bin/v4_10_4_p02d/NuBeam")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/exp/sbnd/app/users/fnicolas/BNBFluxUncertainties/G4BNB/bin/v4_10_4_p02d" TYPE EXECUTABLE FILES "/exp/sbnd/app/users/fnicolas/BNBFluxUncertainties/G4BNB/src/NuBeam")
  if(EXISTS "$ENV{DESTDIR}/exp/sbnd/app/users/fnicolas/BNBFluxUncertainties/G4BNB/bin/v4_10_4_p02d/NuBeam" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/exp/sbnd/app/users/fnicolas/BNBFluxUncertainties/G4BNB/bin/v4_10_4_p02d/NuBeam")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/exp/sbnd/app/users/fnicolas/BNBFluxUncertainties/G4BNB/bin/v4_10_4_p02d/NuBeam"
         OLD_RPATH "/cvmfs/nova.opensciencegrid.org/externals/geant4/v4_10_4_p02d/Linux64bit+3.10-2.17-e19-prof/lib64:/cvmfs/larsoft.opensciencegrid.org/products/root/v6_28_10a/Linux64bit+3.10-2.17-e26-p3915-prof/lib:/cvmfs/larsoft.opensciencegrid.org/products/dk2nudata/v01_10_01f/Linux64bit+3.10-2.17-e26-prof/lib:/cvmfs/larsoft.opensciencegrid.org/products/xerces_c/v3_2_2/Linux64bit+3.10-2.17-e19-prof/lib:/cvmfs/larsoft.opensciencegrid.org/products/clhep/v2_4_1_2/Linux64bit+3.10-2.17-e19-prof/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/exp/sbnd/app/users/fnicolas/BNBFluxUncertainties/G4BNB/bin/v4_10_4_p02d/NuBeam")
    endif()
  endif()
endif()

