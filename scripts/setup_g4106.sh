source /cvmfs/larsoft.opensciencegrid.org/products/setups
source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setup

setup gcc v4_9_3
setup geant4 v4_10_6_p01g -q e26:prof
setup dk2nudata v01_10_01f -q e26:prof
setup cmake v3_9_5
setup jobsub_client
setup ifdhc
# For beamHist
setup boost v1_82_0 -q e26:prof
#For debug
setup gdb v10_1
