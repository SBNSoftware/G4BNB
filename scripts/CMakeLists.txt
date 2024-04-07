set(BOOSTLIB "-L$ENV{BOOST_LIB} -lboost_program_options")

include_directories($ENV{DK2NUDATA_INC}
		    $ENV{ROOT_INC})
#include BOOST as system (-isystem instead of -I) to suppress compiler warnings
include_directories(SYSTEM $ENV{BOOST_INC})

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O3")

add_executable(beamHist beamHist.cc)
target_link_libraries(beamHist ${ROOT_LIBRARIES} $ENV{DK2NUDATA_LIB}/libdk2nuTree.so ${BOOSTLIB})

install(TARGETS beamHist  DESTINATION ${PROJECT_SOURCE_DIR}/bin/$ENV{GEANT4_VERSION})
install(PROGRAMS submitBeam.py compare.py mergeHist.py DESTINATION ${PROJECT_SOURCE_DIR}/bin/$ENV{GEANT4_VERSION})
