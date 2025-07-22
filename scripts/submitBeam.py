#!/usr/bin/env python

import argparse
import os.path
import sys
import tarfile
import re
import subprocess
from  datetime import datetime

parser = argparse.ArgumentParser(description='Submit beam MC jobs.')

parser.add_argument("-n", type=int,
                    choices=range(1,10000),
                    metavar="[1-10000]",
                    required=True,
                    help="Number of jobs to submit.")
parser.add_argument("-b", "--bin",
                    required=True,
                    help="Path to beam g4 executable. Bin dir where exe files are installed.")
parser.add_argument("-i", "--input",
                    required=True,
                    help="Job options input file.")
parser.add_argument("-g", "--geometry",
                    required=True,
                    help="GDML geometry file")
parser.add_argument("-c", "--horn-current", type=float, dest="horncurrent",
                    help="Horn current in A. If specified overwrites /boone/field/horncurrent and /boone/field/skin/SkinDepthHornCurrent commands.")
parser.add_argument("-j", "--jobidoffset", type=int, dest="jobid",
                    default=0,
                    help="Id for running job. If submitting multiple jobs each one is offset by its process number (0..n)")
parser.add_argument("-o", "--output-path", dest="outputpath",
                    default="/pnfs/uboone/scratch/users/%s/beammc"%os.environ['USER'],
                    help="Path where to copy final output. (default=/pnfs/uboone/scratch/users/%s/beammc)"%os.environ['USER'])
parser.add_argument("-gg", "--grid-group", dest="gridgroup",
                    default="sbnd",
                    help="The group to pass to jobsub_submit. Default: sbnd.")
parser.add_argument("-p", "--pot", type=int,
                    help="Protons on target per job. If specified overwrites /run/beamOn command.")
parser.add_argument("-r", "--random-seed", type=int, dest="randomseed",
                    default=1000000,
                    help="Random seed RS. Each process gets seed=$((PROCESS+RS))")
parser.add_argument("-d", "--debug",action='store_true',
                    help="Will not delete submission files in the end. Useful for debugging and will only print the submission command on screen.")

args = parser.parse_args()

flist=[ args.bin+"/NuBeam", args.bin+"beamHist", args.geometry, args.input ]
for ff in  flist:
    if (not os.path.isfile(ff)):
        print("%s does not exist."%ff)
        sys.exit(1)

#prepare the input based on provided file
#this may override some of the options (geometry file, horn current, POT)
inpfilename="production_%i.in"%os.getpid()
jobinpfile=open(inpfilename,'w')
inptemplate=open(args.input)
horncurrent=args.horncurrent
for line in inptemplate:
    if (line.startswith("#")):
        jobinpfile.write(line)
    elif ("/boone/rndm/setRndmSeed" in line):
        jobinpfile.write("/boone/rndm/setRndmSeed RNDSEED\n")
    elif ("/boone/geometry/geometryfile" in line):
        jobinpfile.write("/boone/geometry/geometryfile %s\n"%os.path.basename(args.geometry))
    elif (args.horncurrent is not None and "/boone/field/horncurrent" in line):
        jobinpfile.write("/boone/field/horncurrent %f A\n"%args.horncurrent)
    elif (args.horncurrent is not None and "/boone/field/skin/SkinDepthHornCurrent" in line):
        jobinpfile.write("/boone/field/skin/SkinDepthHornCurrent %f A\n"%args.horncurrent)
    elif (args.horncurrent is None and "/boone/field/horncurrent" in line):
        horncurrent=float(line.split(" ")[1])
        jobinpfile.write(line)
    elif ("/boone/output/filename" in line):
        jobinpfile.write("/boone/output/filename OUTPUTFILE\n")
    elif (args.pot is not None and "/run/beamOn" in line):
        jobinpfile.write("/run/beamOn %i\n"%args.pot)
    else:
        jobinpfile.write(line)

jobinpfile.close()

#now create jobfiles_*.tar that is shipped with the job
#this includes the executable, input file and geometry gdml file
#tar -cf jobfiles.tar --transform '!^[^/]*/!!' file1 file2 file3
if not os.path.exists(args.outputpath):
    os.makedirs(args.outputpath)

tarfilename="%s/jobfiles_%i.tar"%(args.outputpath,os.getpid())
outtar = tarfile.open(tarfilename, mode='w')
outtar.add(args.bin+"/NuBeam", arcname=os.path.basename("NuBeam"))
outtar.add(args.bin+"/beamHist", arcname=os.path.basename("beamHist"))
outtar.add(args.geometry, arcname=os.path.basename(args.geometry))
outtar.add(inpfilename, arcname=os.path.basename(inpfilename))
outtar.close()

#create the sh script that gets submitted
jobinpfn=re.sub('\.in','',os.path.basename(args.input))
geofn=re.sub('\.gdml','',os.path.basename(args.geometry))
jobname="%s_%s_I%iA"%(jobinpfn,geofn,int(horncurrent))

#check setup
upsp=subprocess.Popen(['objcopy', '--dump-section','ups_setup=/dev/stdout',args.bin+"/NuBeam"], 
                      stdout=subprocess.PIPE,
                      stderr=subprocess.PIPE)
upsout,upserr=upsp.communicate()
g4setup, gccsetup, dk2nusetup, ifdhcsetup, boostsetup=None,None,None,None,None
for line in upsout.decode().split("\n"):
    if "geant4" in line.split(" ")[0]:
        g4setup=line
    elif "gcc" in line.split(" ")[0]:
        gccsetup=line
    elif "dk2nu" in line.split(" ")[0]:
        dk2nusetup=line
    elif "boost" in line.split(" ")[0]:
        boostsetup=line
    elif "ifdhc" in line.split(" ")[0] and "ifdhc_config" not in line.split(" ")[0]:
        ifdhcsetup=line

if (g4setup==None or gccsetup==None or dk2nusetup==None or ifdhcsetup==None or boostsetup==None):
    print("ups setup was not embedded in executable.")
    sys.exit(1)

ofstr='''
#!/bin/bash
echo "Running $0 on "$HOSTNAME
echo "Cluster: " ${CLUSTER}
echo "Process: " ${PROCESS}

source /cvmfs/larsoft.opensciencegrid.org/products/setups
setup %(gccsetup)s
setup %(g4setup)s
setup %(dk2nusetup)s
setup %(ifdhcsetup)s
setup %(boostsetup)s

export JOBID=$((PROCESS+%(jobidoffset)s))
let "RNDNO=%(randomseed)s+${JOBID}"
echo "Random number: " ${RNDNO}
cd ${_CONDOR_SCRATCH_DIR}
mkdir ${CLUSTER}
cd ${CLUSTER}
mkdir ${JOBID}
cd ${JOBID}

cp ${INPUT_TAR_DIR_LOCAL}/* . -r

sed -i "s/RNDSEED/${RNDNO}/g" %(inputfile)s
sed -i "s/OUTPUTFILE/NuBeam_%(jobname)s_${JOBID}.dk2nu.root/g" %(inputfile)s

./NuBeam %(inputfile)s
rm NuBeam
./beamHist --input ./\*dk2nu.root --detector-radius 200 --detector-position 73.78 0 11000 --thread 1 --output hist_sbnd_NuBeam_%(jobname)s_${JOBID}.root
./beamHist --input ./\*dk2nu.root --detector-radius 200 --detector-position 0 0 47000 --thread 1 --output hist_uboone_NuBeam_%(jobname)s_${JOBID}.root
./beamHist --input ./\*dk2nu.root --detector-radius 610 --detector-position 0 189.614 54134 --thread 1 --output hist_miniboone_NuBeam_%(jobname)s_${JOBID}.root
./beamHist --input ./\*dk2nu.root --detector-radius 200 --detector-position 0 0 60000 --thread 1 --output hist_icarus_NuBeam_%(jobname)s_${JOBID}.root

rm beamHist

cd ..
ifdh mkdir %(outputdir)s
ifdh mkdir %(outputdir)s/%(jobname)s
ifdh mkdir %(outputdir)s/%(jobname)s/%(now)s
ifdh cp -r ${JOBID} %(outputdir)s/%(jobname)s/%(now)s/${JOBID}
'''%{'g4setup':g4setup,
     'gccsetup':gccsetup,
     'dk2nusetup':dk2nusetup,
     'ifdhcsetup':ifdhcsetup,
     'boostsetup':boostsetup,
     'jobidoffset':str(args.jobid),
     'randomseed':str(args.randomseed),
     'inputfile':inpfilename,
     'jobname':jobname,
     'outputdir':args.outputpath,
     'now':datetime.now().strftime('%Y%m%d-%H:%M')}

runjobfname="runjob_%i.sh"%os.getpid()
of=open(runjobfname,'w')
of.write(ofstr)
of.close()

#Create submit command
jobsub_options=" --group=%s --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --role=Analysis --memory 2000MB --expected-lifetime=8h --singularity-image /cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest "%(args.gridgroup)
cmd="jobsub_submit" + jobsub_options + "-N %i --tar_file_name=dropbox://%s file://%s"%(args.n,os.path.abspath(tarfilename),os.path.abspath(runjobfname))

if (not args.debug):
    print("Running submit cmd:")
    print(cmd)
    os.system(cmd)
else:
    print("Would have ran:")
    print(cmd)

#Delete temp files unless debugging
if (not args.debug):
   os.remove(inpfilename)
   os.remove(runjobfname)
   os.remove(tarfilename)
