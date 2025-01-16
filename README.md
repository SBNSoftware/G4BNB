# Booster Neutrino Beam Simulation project

This repository has been migrated from [Redmine](https://cdcvs.fnal.gov/redmine/projects/booster-neutrino-beamline/wiki/Git_Help).

The current version is based on the MiniBooNE BooNEG4Beam.
The original code can be found under "BooNE-BackTrack subproject":https://cdcvs.fnal.gov/redmine/projects/boone-backtrack/repository.

The old code has been cleaned up, partly rewritten and rearranged to get it up to date and working against current versions of geant4.
[Wiki_old](https://cdcvs.fnal.gov/redmine/projects/booster-neutrino-beamline/wiki/Wiki_old) covers some of the validation work by Paul Lebrun while transitioning.
There were two repositories within the main Redmine project page, one used during the upgrade (booster-neutrino-beamline) and one with the final working version where code has been further restructured (booster-neutrino-beamline-g4bnb).
The latter one has been migrated here to GitHub.

The current working version of G4BNB v1.0 uses GEANT4.10.4, and can be cloned using the geant4_10_4 branch.
A branch is also available that uses GEANT4.10.6, named g4bnb-geant4-10-6-build. The master branch can be built with GEANT4.10.1.
Relevant set up scripts for each version are provided in the branch.

## Downloading the G4BNB code

The G4BNB software is maintained in this git repository.

Authenticated clone (i.e., allows to modify the software and upload your modification on this site, via git push):

```
git clone https://github.com/SBNSoftware/G4BNB.git
``` 

To check out a preexisting branch from the repository (recommended to run the current working version or development versions), you can clone the repository as above, then do:
```
$ cd <repository directory name>
$ git checkout origin/<branch name>
$ git checkout -b <branch name>
```
This will create a local branch from the remote branch. The current working version of G4BNB is on the branch "geant4_10_4". 

To check out a fixed release of G4BNB (which is recommended for large production runs) into a directory with the name of the tag, you can do:
```
$ git clone https://github.com/SBNSoftware/G4BNB.git <tag_name>
$ cd <tag_name>
$ git checkout <tag_name>
```
You can see a list of available tags by doing:
```
git tag
```

## Building G4BNB at Fermilab

If building G4BNB on the Fermilab gpvms, an SL7 container will be needed. An apptainer can be activated by:
```
sh /exp/$(id -ng)/data/users/vito/podman/start_SL7dev.sh
```
This apptainer does not support jobsub.

To setup the correct software dependancies, use the setup scripts provided in the `/scripts` directory.
The setup scripts have the relevant GEANT4 version number in their name.
Assuming you are using the current `v1.0` working version (branch `geant4_10_4`), setup by running:
```
source scripts/setup_g4104.sh
```

cd into build directory and run:
```
cmake ../
make -j4
make install
```

## Submit grid jobs

It is recommended to submit jobs directly from the AL9 gpvm node without using an SL7 container.
The `submitBeam.py` script sets up a singularity image of SL7 on the grid.
To submit grid jobs use submitBeam.py script:

```
usage: submitBeam.py [-h] -n [1-10000] -b BIN -i INPUT -g GEOMETRY
                     [-c HORNCURRENT] [-j JOBID] [-o OUTPUTPATH] [-p POT]
                     [-r RANDOMSEED] [-d]

Submit beam MC jobs.

optional arguments:
  -h, --help            show this help message and exit
  -n [1-10000]          Number of jobs to submit.
  -b BIN, --bin BIN     Path to beam g4 executable. Bin dir where exe files
                        are installed.
  -i INPUT, --input INPUT
                        Job options input file.
  -g GEOMETRY, --geometry GEOMETRY
                        GDML geometry file
  -c HORNCURRENT, --horn-current HORNCURRENT
                        Horn current in A. If specified overwrites
                        /boone/field/horncurrent and
                        /boone/field/skin/SkinDepthHornCurrent commands.
  -j JOBID, --jobidoffset JOBID
                        Id for running job. If submitting multiple jobs each
                        one is offset by its process number (0..n)
  -o OUTPUTPATH, --output-path OUTPUTPATH
                        Path where to copy final output.
                        (default=/pnfs/uboone/scratch/users/${USER}/beammc)
  -p POT, --pot POT     Protons on target per job. If specified overwrites
                        /run/beamOn command.
  -r RANDOMSEED, --random-seed RANDOMSEED
                        Random seed RS. Each process gets seed=$((PROCESS+RS))
  -d, --debug           Will not delete submission files in the end. Useful
                        for debugging and will only print the submission
                        command on screen.
```

Example to submit 100 jobs with 50000 POT each using `production.in` input file as template (note that `submitBeam.py` overrides POT, horn current and geometry):
```
bin/v4_10_4_p02d/submitBeam.py -n 100 -b bin/v4_10_4_p02d/ -i input/production.in -g geometry/BooNE_50m.gdml  -p 50000
```

Different production scripts have been provided to enable decay at rest (DAR) for only muons and both muons and pions. These are called `production_muDAR.in` and `production_muDAR_piDAR.in` respectively.
The production ntuple can be saved in the output files by uncommenting the line that reads
```
#/boone/output/saveProductionNtuple true
```

## Analyzing beam MC

The output of the beam MC are [dk2nu beam ntuples](https://cdcvs.fnal.gov/redmine/projects/dk2nu/wiki/Wiki).
Code to produce standard set of histograms is in `scripts/beamHist.cc` and gets compiled with beam MC.
```
Options:
  -h [ --help ]                     Print help message
  -i [ --input ] arg                Path pattern for input files. Put it in 
                                    quotes or escape *.
  -o [ --output ] arg (=hist.root)  Output file name.
  -p [ --pot ] arg                  POT used for normalization (overides 
                                    counting using info in meta tree and speeds
                                    up process). 
                                    Total POT should be given (number of files 
                                    X POT per file).
  -n [ --nredecays ] arg (=1)       Number of redecays.
  -r [ --detector-radius ] arg (=0) Detector radius (in cm).
  -x [ --detector-position ] arg    Detector position (in cm).
  -t [ --thread ] arg (=1)          Number of threads to run. (max set to 8)
```

Example to produce histograms from previously generated jobs:
```
bin/v4_10_4_p02d/beamHist -i /pnfs/uboone/scratch/users/${USER}/beammc/production_BooNE_50m_I174000A/\*/\*root -r 610 -x 0 189.614 54134 -t4 -o hist_miniboone.root
```

This might take a while if analyzing many jobs. submitBeam.py actually runs the beamHist command on the grid and produces histogram files for sbnd, uboone, miniboone, icarus locations:

 |                    |        x/cm   |         y/cm   |            z/cm | r/cm |
 | ---|---|---|---|---|
 |  sbnd             |         0          |         0            |    11000   | 200   |
 | uboone          |         0          |         0            |    47000 |  200   |
 | miniboone     |         0          | 189.614          |    54134  |  610   |
 | icarus             |        0           |        0             |   60000  |  200   |


One can merge those histograms to save some time. For this use `mergeHist.py` script:
```
usage: mergeHist.py [-h] -i INPUT [-o OUTPUT] [-l LOCATION]

Merge beam MC histograms.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to top directory containing histograms to merge.
  -o OUTPUT, --output OUTPUT
                        Output root file.
  -l LOCATION, --location LOCATION
                        Used to pick the correct hist file. (grid jobs create
                        sbnd, uboone, icarus, miniboone by default.)
```

For example to merge histograms from previously ran jobs:
```
bin/v4_10_4_p02d/mergeHist.py -i /pnfs/uboone/scratch/users/${USER}/beammc/production_BooNE_50m_I174000A/  -l miniboone
```
This creates hist_miniboone.root. Merge script scales histograms dividing out the number of files so normalization is always the same, assuming the files each have the same POT in them.

## Making plots

Histograms can be plotted in the usual way using root. 
`compare.py` script is provided that can overlay multiple histograms and produce standard set of plots.
```
usage: compare.py [-h] -d DATA -c COMPARE [COMPARE ...] [-o OUTPUT]
compare.py: error: argument -d/--data is required
```

Example to compare to the original MiniBooNE flux histogram:
```
bin/v4_10_4_p02d/compare.py -d validation/april07_baseline_rgen610.6_fixrnd_20171212.root -c hist_miniboone.root -o miniboone_plot.pdf
```

Multiple histogram files can be listed after -c option and those will all be compared to file under `-d` option.
