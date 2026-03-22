#include <iostream>
#include <sstream>
#include <string>
#include <glob.h>
#include <thread>
#include <array>

#include "TFile.h"
#include "TH1F.h"
#include "TH3F.h"
#include "TChain.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TThread.h"

#include "dk2nu/tree/dk2nu.h"
#include "dk2nu/tree/dkmeta.h"
#include "dk2nu/tree/calcLocationWeights.h"

#include <boost/program_options.hpp>

using namespace std;
using namespace boost::program_options;

struct histpackage_t
{
  std::vector<string> filelist;
  int NREDECAY;
  std::vector<double> detpos; 
  std::vector<double> prism_bins;
  double RDet;
  int rndSeed;
  bool countPOT;
  std::vector<TH3F*> hxye;
  std::vector<std::array<TH1F*, 4>> hFlux;
  std::vector<std::array<std::array<TH1F*, 4>, 4>> hparent;
  std::vector<std::array<std::array<TH1F*, 5>, 4>> hsec;
  double POT;
};

//thread function (has to be void* since root uses this to determine if it is detached or non-detached thread)  
void* FillHist(void* hp);

int main(int ac, char* av[])
{
  vector<double> detpos;
  vector<double> prism_bins;
  detpos.resize(3);
  double RDet;
  string searchpath;
  string outputfn;
  int NREDECAY;
  double userPOT=0;
  int nthread=1;
  options_description opt("Options");
  opt.add_options()
    ("help", "Print help message")
    ("input",value<string>(&searchpath),"Path pattern for input files. Put it in quotes or escape *.")
    ("output",value<string>(&outputfn)->default_value("hist.root"),"Output file name.")
    ("pot",value<double>(&userPOT),"POT used for normalization (overides counting using info in meta tree and speeds up process). \nTotal POT should be given (number of files X POT per file).")
    ("nredecays",value<int>(&NREDECAY)->default_value(1.),"Number of redecays.")
    ("detector-radius",value<double>(&RDet)->default_value(0.),"Detector radius (in cm).")
    ("detector-position",value<vector<double> >(&detpos)->multitoken(),"Detector position (in cm).")
    ("thread",value<int>(&nthread)->default_value(1),"Number of threads to run. (max set to 8)")
    ("prism-bins",value<vector<double>>(&prism_bins)->multitoken(),"PRISM bin edges (in degrees)");
    
  variables_map vm;
  
  try {
    store(parse_command_line(ac,av,opt, command_line_style::unix_style ^ command_line_style::allow_short),vm);
    notify(vm);
    if (vm.count("help")) {
      cerr<<opt<<endl;
      return 1;
    } 
    if (!vm.count("input")) {
      cerr<<"Need to provide input pattern."<<endl;
      cerr<<opt<<endl;
      return 1;
    }
    if (!vm.count("detector-position")) {
      //assume it is uboone
      detpos[0]=0;
      detpos[1]=0;
      detpos[2]=47000.;
    }
    if (!vm.count("prism-bins")) {
      // Set no PRISM
      prism_bins.push_back(0.0);
    }
  } catch (error& e) {
    cerr << e.what()<<endl<<endl;
    cerr << opt <<endl;
    return 1;
  }
  if (nthread>8) nthread=8;

  glob_t glob_result;
  cout<<"Searching "<<searchpath<<endl;
  glob(searchpath.c_str(),GLOB_TILDE,NULL,&glob_result);
  vector<string> filelist;
  for(unsigned int i=0;i<glob_result.gl_pathc;++i){
    filelist.push_back(string(glob_result.gl_pathv[i]));
  }
  globfree(&glob_result);
  cout <<"Found "<<filelist.size()<<" files"<<endl;

  //  cout<<"Total POT: "<<hp->POT<<endl;
  //cout<<"(with redecaying "<<hp->NREDECAY<<" times)"<<endl;
  //if (dkmetaTree) {
  // cout <<"POT calculated using meta data"<<endl;
  //} else {
  // cout <<"POT set using --pot option."<<endl;
  //}
  cout <<"Making histograms for detector at r=("
       <<detpos[0]<<", "
       <<detpos[1]<<", "
       <<detpos[2]<<") cm and smearing over RDet="<<RDet<<" cm"<<endl;
  cout <<"Redecaying "<<NREDECAY<<" times."<<endl;
  const int nPRISM = prism_bins.size()-1;
  std::ostringstream asts;
  asts << "Using N = " << nPRISM << " PRISM bins with edges: [";
  for( const double & pr : prism_bins ) { asts << pr << ", "; }
  asts.seekp(-2, asts.cur); asts << "] ";
  std::cout << asts.str() << std::endl;

  cout<<"Starting "<<nthread<<" threads"<<endl;
  TThread::Initialize();
  TThread* t[nthread];
  histpackage_t* hp[nthread];

  //prepare filelists
  std::vector<string> tfl[nthread];
  while (filelist.size()>0) 
    for (int i=0;i<nthread;i++) {
      if (filelist.size()>0) {
	tfl[i].push_back(filelist.back());
	filelist.pop_back();
      }
    }

  for (int i=0;i<nthread;i++) {
    hp[i]=new histpackage_t();
    hp[i]->filelist=tfl[i];
    hp[i]->NREDECAY=NREDECAY;
    hp[i]->detpos=detpos;
    hp[i]->RDet=RDet;
    hp[i]->rndSeed=i;
    hp[i]->prism_bins=prism_bins;
    if (!vm.count("pot"))
      hp[i]->countPOT=true;
    else 
      hp[i]->countPOT=false;
    t[i]=new TThread(Form("Thread_%i",i),FillHist, (void*) hp[i]);
    t[i]->Run();
  }
  TThread::Ps();
  for (int i=0;i<nthread;i++) {
    t[i]->Join();
  }

  //add histograms from all threads, for all PRISM bins
  for (int i=1;i<nthread;i++) {

    for (int p=0;p<=nPRISM;p++) {
      hp[0]->hxye[p]->Add(hp[i]->hxye[p]);
      for (int inu=0;inu<4;inu++) {
	hp[0]->hFlux[p][inu]->Add(hp[i]->hFlux[p][inu]);
	for (int ipar=0;ipar<4;ipar++) {
	  hp[0]->hparent[p][inu][ipar]->Add(hp[i]->hparent[p][inu][ipar]);
	}
	for (int isec=0;isec<5;isec++) {
	  hp[0]->hsec[p][inu][isec]->Add(hp[i]->hsec[p][inu][isec]);
	} 
      }
    } // add in all PRISM bins

    hp[0]->POT+=hp[i]->POT;
  } // sum over threads

  double totPOT=hp[0]->POT;
  if (vm.count("pot")) {
    totPOT=userPOT;
    cout <<"POT set using --pot option to "<<totPOT<<endl; 
  } else {
    cout <<"Total POT summed using meta data= "<<totPOT<<endl;
  }

  for (int p=0;p<=nPRISM;p++) {
    for (int inu=0;inu<4;inu++) {
      hp[0]->hFlux[p][inu]->Scale(1./totPOT);
      for (int ipar=0;ipar<4;ipar++) {
	hp[0]->hparent[p][inu][ipar]->Scale(1./totPOT);
      }
      for (int isec=0;isec<5;isec++) {
	hp[0]->hsec[p][inu][isec]->Scale(1./totPOT);
      } 
    }
  } // scale each PRISM bin

  //write histograms to file
  // Each collection of histograms gets its own TDirectoryFile.
  TFile fout(outputfn.c_str(),"RECREATE");
  for (int p=0;p<=nPRISM;p++) {
    std::string dirname, dirtitle;
    if( p == 0 ) {
      dirname  = std::string("inclusive");
      dirtitle = std::string("Inclusive flux");
    } else {
      dirname  = std::string(Form("prism-%02i", p));
      double plow = prism_bins[p-1]; double phigh = prism_bins[p];
      dirtitle = std::string(Form("Flux in PRISM bin [%2.1f, %2.1f]", plow, phigh));
    }
    TDirectoryFile * dirfile = new TDirectoryFile(dirname.c_str(), dirtitle.c_str());
    fout.Add(dirfile);
    dirfile->cd();
    hp[0]->hxye[p]->Write();
    for (int inu=0;inu<4;inu++) {
      hp[0]->hFlux[p][inu]->Write();
      for (int ipar=0;ipar<4;ipar++) {
	hp[0]->hparent[p][inu][ipar]->Write();
      }
    }  
    for (int inu=0;inu<4;inu++) {
      hp[0]->hFlux[p][inu]->Write(Form("h70%i",inu+1)); //same as h50x, keeping copy 
      //to be consistent with MB files
      for (int isec=0;isec<5;isec++) {
	hp[0]->hsec[p][inu][isec]->Write();
      } 
    }
    fout.cd();
  } // prism bins
  fout.cd();
  fout.Close();

  return 0;
}

void* FillHist(void* hpvoid)
{
  histpackage_t* hp=(histpackage_t*) hpvoid;
  TThread::Lock();
  TChain* dk2nuTree=new TChain("dk2nuTree");
  TChain* dkmetaTree=NULL;
  const int nPRISM = (hp->prism_bins).size()-1;
  if (hp->countPOT)
    dkmetaTree=new TChain("dkmetaTree");
  for (auto ifile : hp->filelist) {
    dk2nuTree->Add(ifile.c_str());
    if (dkmetaTree)
      dkmetaTree->Add(ifile.c_str());
  }

  bsim::Dk2Nu*  dk2nu  = new bsim::Dk2Nu;
  dk2nuTree->SetBranchAddress("dk2nu",&dk2nu);
  bsim::DkMeta* dkmeta  = new bsim::DkMeta;
  if (dkmetaTree)
    dkmetaTree->SetBranchAddress("dkmeta",&dkmeta);

  TThread::UnLock();

  Long64_t ientry=0;
  hp->POT=0;
  if (dkmetaTree) {
    while (dkmetaTree->GetEntry(ientry++)) {
      hp->POT+=dkmeta->pots;
    }
  }
  hp->POT*=hp->NREDECAY;

  string nutype[]={    "nue",        "nuebar",      "numu",         "numubar"};
  string nultx[] ={"#nu_{e}", "#bar{#nu}_{e}", "#nu_{#mu}", "#bar{#nu}_{#mu}"};
  int pdgcode[]  ={       12,             -12,          14,               -14};
  string pltx[]  ={"#mu^{#pm}","#pi^{#pm}","K^{0}_{L}","K^{#pm}"};
  string secltx[]  ={"pBe->#pi^{#pm}->...->#mu^{#pm}",
		     "pBe->#pi^{#pm}->..(not #mu^{#pm})..",
		     "pBe->K^{0}_{L}->...",
		     "pBe->K^{#pm}->...",
		     "pBe->(p or n)->..."};


  TRandom3 rndmno;
  rndmno.SetSeed(hp->rndSeed);  
  std::string suffix="";
  if (hp->rndSeed>0) 
    suffix=Form("_%i",hp->rndSeed);
  TThread::Lock();
  for (int p=0;p<=nPRISM;p++) {
    std::string psuffix = (p > 0) ? std::string(Form("_prism%02i", p)) : "";
    hp->hxye.push_back(new TH3F(Form("h_xyE%s%s",suffix.c_str(),psuffix.c_str()),
				Form("Neutrino vertices at r=(%f,%f,%f)cm;x (cm);y (cm);E (GeV)",
				     hp->detpos[0],hp->detpos[1],hp->detpos[2]),
				100,-hp->RDet,hp->RDet,
				100,-hp->RDet,hp->RDet,
				200,0,10));
    
    std::array<TH1F*, 4> hpFlux;
    for (int i=0;i<4;i++) {
      hpFlux[i] = new TH1F(Form("h50%i%s%s",i+1,suffix.c_str(),psuffix.c_str()),
			   Form("%s (all);Energy %s (GeV);#phi(%s)/50MeV/POT",
				nultx[i].c_str(),nultx[i].c_str(),nultx[i].c_str()),
			   200,0,10);
      hpFlux[i]->Sumw2();
    }
    hp->hFlux.push_back(hpFlux);
    
    std::array<std::array<TH1F*, 4>, 4> hpparent;
    std::array<std::array<TH1F*, 5>, 4> hpsec;
    for (int inu=0;inu<4;inu++) {
      for (int ipar=0;ipar<4;ipar++) {
	hpparent[inu][ipar] = new TH1F(Form("h5%i%i%s%s",ipar+1,inu+1,suffix.c_str(),psuffix.c_str()),
				       Form("...->%s->%s;Energy %s (GeV);#phi(%s)/50MeV/POT",
					    pltx[ipar].c_str(),nultx[inu].c_str(),
					    nultx[inu].c_str(),nultx[inu].c_str()),
				       200,0,10);
	hpparent[inu][ipar]->Sumw2();
      }
      for (int isec=0;isec<5;isec++) {
	hpsec[inu][isec] = new TH1F(Form("h7%i%i%s%s",isec+1,inu+1,suffix.c_str(),psuffix.c_str()),
				    Form("%s->%s;Energy %s (GeV);#phi(%s)/50MeV/POT",
					 secltx[isec].c_str(),nultx[inu].c_str(),
					 nultx[inu].c_str(),nultx[inu].c_str()),
				    200,0,10);
	hpsec[inu][isec]->Sumw2();
      }
    }
    hp->hparent.push_back(hpparent);
    hp->hsec.push_back(hpsec);
  } // for all PRISM bins, initialise

  // Note this calculation assumes all z = detpos[2].
  // For more accurate PRISM fluxes you'll need a 3D profile.
  TThread::UnLock();
  ientry=0;
  cout<<"Thread "<<hp->rndSeed<<" starting to process "<<dk2nuTree->GetNtrees()<<" files."<<endl;
  while (dk2nuTree->GetEntry(ientry++)) {
    //    if (ientry%100000==0) cout<<"Thread "<<hp->rndSeed<<" on entry "<<ientry<<endl;
    for (int ipdg=0;ipdg<4;ipdg++) {
      if (dk2nu->decay.ntype!=pdgcode[ipdg]) continue;
      
      for (int iredecay=0;iredecay<hp->NREDECAY;iredecay++) {
	double enu,wgt_xy;
	double xx=rndmno.Uniform(-hp->RDet,hp->RDet);
	double yy=rndmno.Uniform(-hp->RDet,hp->RDet);
	while (sqrt(xx*xx+yy*yy)>hp->RDet) {
	  xx=rndmno.Uniform(-hp->RDet,hp->RDet);
	  yy=rndmno.Uniform(-hp->RDet,hp->RDet);
	}
	TVector3 xyz(xx+hp->detpos[0],yy+hp->detpos[1],hp->detpos[2]);
	bsim::calcEnuWgt(dk2nu,xyz,enu,wgt_xy);
	//to compare with FluxForNuance output (MiniBooNE files)
	//normalize through whole detector area in m2
	double totwgh=wgt_xy*dk2nu->decay.nimpwt/3.14159*hp->RDet*hp->RDet*3.14159*1e-4;

	int firstInelastic=0;
	while (dk2nu->ancestor[firstInelastic].proc.find("HadronInelastic")==string::npos) firstInelastic++;

	for (int p=0;p<=nPRISM;p++){
	  hp->hxye[p]->Fill(xx,yy,enu,totwgh);
	  hp->hFlux[p][ipdg]->Fill(enu,totwgh);
	
	  if (dk2nu->decay.ptype==13 || dk2nu->decay.ptype==-13) //mu+-
	    hp->hparent[p][ipdg][0]->Fill(enu,totwgh);
	  else if (dk2nu->decay.ptype==211 || dk2nu->decay.ptype==-211) //pi+-
	    hp->hparent[p][ipdg][1]->Fill(enu,totwgh);
	  else if (dk2nu->decay.ptype==130) //K0L
	    hp->hparent[p][ipdg][2]->Fill(enu,totwgh);
	  else if (dk2nu->decay.ptype==321 || dk2nu->decay.ptype==-321) //K+-
	    hp->hparent[p][ipdg][3]->Fill(enu,totwgh);
	
	  if (fabs(dk2nu->ancestor[firstInelastic].pdg)==211 && fabs(dk2nu->decay.ptype)==13)
	    hp->hsec[p][ipdg][0]->Fill(enu,totwgh);
	  else if (fabs(dk2nu->ancestor[firstInelastic].pdg)==211)
	    hp->hsec[p][ipdg][1]->Fill(enu,totwgh);
	  else if (fabs(dk2nu->ancestor[firstInelastic].pdg)==130)
	    hp->hsec[p][ipdg][2]->Fill(enu,totwgh);
	  else if (fabs(dk2nu->ancestor[firstInelastic].pdg)==321)
	    hp->hsec[p][ipdg][3]->Fill(enu,totwgh);
	  else if (dk2nu->ancestor[firstInelastic].pdg==2212 || dk2nu->ancestor[firstInelastic].pdg==2112)
	    hp->hsec[p][ipdg][4]->Fill(enu,totwgh);
	} // loop over PRISM bins
      } // loop over redecays
    } // loop over nu pdg
  } // loop over dk2nu entries
 
  cout<<"Thread "<<hp->rndSeed<<" processed "<<ientry<<" entries. POT = "<<hp->POT<<endl;

  return NULL;
}
