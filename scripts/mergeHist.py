#!/usr/bin/env python

import argparse
import os.path
import sys
import tarfile
import re
import subprocess
tmpargv=sys.argv[:]
#add -b argument for root to see and then remove it so argparse doesn't complain 
#also remove -h because stupid root hijacks it
sys.argv=[tmpargv[0],'-b']
from ROOT import TFile, TH1D, TCanvas, TPad, TLine, TLegend, gStyle, gROOT, kRed, kBlue, kMagenta, kBlack, kGreen
sys.argv=tmpargv

parser = argparse.ArgumentParser(description='Submit beam MC jobs.')

parser.add_argument("-i", "--input",
                    required=True,
                    help="Path to top directory containing histograms to merge.")
parser.add_argument("-o", "--output",
                    help="Output root file.")
parser.add_argument("-l", "--location",
                    default="uboone",
                    help="Used to pick the correct hist file. (grid jobs create sbnd, uboone, icarus, miniboone by default.)")

args = parser.parse_args()

flist=[]

for subdir, dirs, files in os.walk(args.input):
    for f in files:
        if "hist_%s"%args.location in f:
            flist.append(subdir+"/"+f)

hlist={}
nFiles=0
for f in flist:
    froot=TFile(f)
    isGood=False
    for k in froot.GetListOfKeys():
        if "TH" not in k.GetClassName():
            continue
        if (k.ReadObj().GetSumOfWeights()>0):
            isGood=True
            if k.GetName() in hlist:            
                hlist[k.GetName()].Add(k.ReadObj())
            else:
                gROOT.cd()
                hlist[k.GetName()]=k.ReadObj().Clone(k.GetName())
    if (isGood):
        nFiles=nFiles+1

print("Added %i good files."%nFiles)
outfname="hist_"+args.location+".root"
if args.output:
    outfname=args.output

print("Writing histograms in %s"%outfname)
fout=TFile(outfname,"RECREATE")
keylist=hlist.keys()
sorted(keylist)
for k in keylist:
    print(k)
    hlist[k].Scale(1./float(nFiles))
    hlist[k].Write()
fout.Close()
