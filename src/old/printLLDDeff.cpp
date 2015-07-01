#include <vector>
#include <sstream>
#include <iostream>
#include <string>
#include <time.h>

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TEntryList.h"
#include "TCanvas.h"
#include "TBranch.h"
#include "TCut.h"
#include "TGraphErrors.h"

#include "general_functions.hpp"
#include "ReadTree_comp.hpp"
#include "analyser.hpp"
#include "Lb_cuts.hpp"

using namespace std;
using namespace RooFit;


int main(int argc, char **argv)
{
	vector <string> effnames;
	effnames.push_back("reco");
	effnames.push_back("trig");
	effnames.push_back("mva");

	TString LLfilename = "/afs/cern.ch/user/p/pluca/work/Lb/Lmumu/results/Lbreleff_LL.root";	
	TString DDfilename = "/afs/cern.ch/user/p/pluca/work/Lb/Lmumu/results/Lbreleff_DD.root";

	TCanvas * cc = new TCanvas();
	
	for(unsigned s = 0; s < effnames.size(); s++)
	{
		TGraphErrors * effgrLL, * effgrDD;
		TFile * fLL = TFile::Open(LLfilename);
		fLL->GetObject(("effrel_"+effnames[s]).c_str(),effgrLL);
		TFile * fDD = TFile::Open(DDfilename);		
		fDD->GetObject(("effrel_"+effnames[s]).c_str(),effgrDD);
		
		effgrLL->SetTitle(effnames[s].c_str());
		effgrLL->Draw("AP");
		cc->Print((effnames[s]+"_LL.pdf").c_str());
		effgrDD->SetTitle(effnames[s].c_str());
		effgrDD->Draw("AP");
		cc->Print((effnames[s]+"_DD.pdf").c_str());		
	}
	
	return 0;
}





