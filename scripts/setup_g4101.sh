source /cvmfs/larsoft.opensciencegrid.org/products/setups
source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setup
setup gcc v4_9_3
setup geant4 v4_10_1_p03 -q e9:prof
setup dk2nu v01_03_00b -q e9:prof:r5
setup cmake v3_3_2
setup jobsub_client
setup ifdhc

#use boost for beamHist.c
setup boost v1_57_0a -q e9:prof

# Set the G4BNB_DIR env variable
export G4BNB_DIR=$(dirname $( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ))
