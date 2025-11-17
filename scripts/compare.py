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
print(sys.argv)
def getLevelOfAgreement(h1,h2):
    low1,low2,mid1,mid2,high1,high2=0,0,0,0,0,0
    for ibin in range(1,h1.GetNbinsX()+1):
        if (h1.GetBinLowEdge(ibin)<0.2):
            low1=low1+h1.GetBinContent(ibin)
            low2=low2+h2.GetBinContent(ibin)
        elif (h1.GetBinLowEdge(ibin)<2.):
            mid1=mid1+h1.GetBinContent(ibin)
            mid2=mid2+h2.GetBinContent(ibin)
        elif (h1.GetBinLowEdge(ibin)<3.):
            high1=high1+h1.GetBinContent(ibin)
            high2=high2+h2.GetBinContent(ibin)

    if low1==0:
        low1=1
    if mid1==0:
        mid1=1
    if high1==0: 
        high1=1 
    return [low2/low1,mid2/mid1,high2/high1]

def setStyle():
    gStyle.SetOptStat(0)
    style_label_font=42
    gStyle.SetLabelFont(style_label_font,"xyz")
    gStyle.SetLabelSize(0.045,"xyz")
    gStyle.SetLabelOffset(0.015,"xyz")
    gStyle.SetStatFont(style_label_font)
    gStyle.SetTitleFont(style_label_font,"xyz") # axis titles
    gStyle.SetTitleFont(style_label_font,"h")   # histogram title
    gStyle.SetTitleSize(0.065,"xyz")            # axis titles
    gStyle.SetTitleSize(0.05,"h")               # histogram title
    gStyle.SetTitleOffset(1.0,"x")
    gStyle.SetTitleOffset(0.6,"y")
    
    gStyle.SetFrameLineWidth(2)
    gStyle.SetFuncWidth(2)
    gStyle.SetHistLineWidth(2)
    
    gROOT.ForceStyle()

parser = argparse.ArgumentParser(description='Submit beam MC jobs.')

parser.add_argument("-d", "--data",
                    required=True,
                    help="Histogram file to compare to (plotted with dots). When plotting ratios this histogram is used as denominator")
parser.add_argument("-c", "--compare",
                    required=True,nargs='+',
                    help="List histogram file(s) being compared to \"data\" histogram (plotted as line).")
parser.add_argument("-o", "--output",
                    default="plots.pdf",
                    help="Output pdf file with all the plots.")


args = parser.parse_args()
hcol=[kRed, kBlue, kMagenta, kGreen+2]

for ff in [args.data]+args.compare:
    if (not os.path.isfile(ff)):
        print("%s does not exist."%ff)
        sys.exit(1)

setStyle()
fdata=TFile(args.data)
fmc={}
for ff in args.compare:
    fmc[ff]=TFile(ff)

rat={}
for inu in range(1,5):
    rat[inu]={}
    for isec in range(0,6):
        ntype="numu"
        if inu==1:
            ntype="nue"
        elif inu==2:
            ntype="nuebar"
        elif inu==4:
            ntype="numubar"

        hdata=fdata.Get("h7%i%i"%(isec,inu))
        hmc={}
        icol=0
        for ff in args.compare:
            hmc[ff]=fmc[ff].Get("h7%i%i"%(isec,inu))
            hmc[ff].Scale(500./500.)
            hdata.SetTitle(hmc[ff].GetTitle())
            hdata.SetMaximum(max(hdata.GetMaximum()*1.1,hmc[ff].GetMaximum()*1.1))
            hmc[ff].SetLineColor(hcol[ int(icol%len(hcol)) ])
            hmc[ff].SetLineStyle( int(icol/len(hcol))+1)
            icol=icol+1
            if isec==0:
                rat[inu][ff]=getLevelOfAgreement(hmc[ff],hdata)

        c=TCanvas("c%i%i"%(isec,inu),"",600,700)
        p1=TPad("p1","p1",0,0.5,1,1)
        p1.SetBottomMargin(0.);
        p1.SetLeftMargin(0.12);
        p1.Draw();
        p2=TPad("p2","p2",0,0.,1,0.5);
        p2.SetTopMargin(0);
        p2.SetBottomMargin(0.2);
        p2.SetLeftMargin(p1.GetLeftMargin());
        p2.Draw();
        p1.cd();
        hdata.SetMarkerStyle(20);
        hdata.GetXaxis().SetRangeUser(0,3)
        hdata.SetStats(0)
        hdata.GetYaxis().SetTitle("Flux")
        hdata.DrawCopy("pe")
        leg=TLegend(0.5,0.65,0.88,0.88)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.AddEntry(hdata,os.path.splitext(os.path.basename(args.data))[0])
        for ff in args.compare:
            hmc[ff].DrawCopy("samehist")
            hname=os.path.basename(ff)
            hname=os.path.splitext(hname)[0]
            print(ff,hname)
            leg.AddEntry(hmc[ff],hname,"L")
        leg.Draw()
        p2.cd()
        drawarg="hist"
        for ff in args.compare:
            hmc[ff].Divide(hdata)
            hmc[ff].SetTitle("")
            hmc[ff].GetYaxis().SetTitle("Ratio")
            hmc[ff].SetMinimum(0.4)
            hmc[ff].SetMaximum(1.6)
            hmc[ff].GetXaxis().SetRangeUser(0,3)
            hmc[ff].SetStats(0)
            hmc[ff].DrawCopy(drawarg)
            drawarg="samehist"

        l1=TLine(0,1,3,1);
        l1.Draw();
        if (inu==1 and isec==0):
            c.Print("%s("%args.output)
        elif (inu==4 and isec==5):
            c.Print("%s)"%args.output)
        else:
            c.Print(args.output)


for inu in range(1,5):
     for ff in args.compare:
         print("%40s %5.2f %5.2f %5.2f"%(ff,rat[inu][ff][0],rat[inu][ff][1],rat[inu][ff][2]))
