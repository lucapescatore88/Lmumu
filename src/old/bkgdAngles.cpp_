#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraph.h"
#include <string>
#include "TString.h"
#include "Lb_cuts.hpp"

using namespace std;

int main(int argc, char **argv)
{
	TString type = "all";	
	if(argc > 1)
	{
		for(int a = 1; a < argc; a++)
		{
			string arg = argv[a];
			string str = arg.substr(2,arg.length()-2);
			if(arg.find("-t") != string::npos) type = (TString)str;
		}
	}

	
	TString baseCut = "";
	if(type == "DD") baseCut = CutsDef::DDcut;
	else if(type == "LL") baseCut = CutsDef::LLcut;


	TFile * trainFile = TFile::Open("/afs/cern.ch/work/p/pluca/weighted/Lmumu/Lb2Lmumu_CL_NBweighted.root");
	TTree * tree = (TTree *)trainFile->Get("tree");

	TCanvas * c =  new TCanvas();

	gStyle->SetOptStat(0);
	


	for(int i = 0; i < 9; i++)
	{
		TString cut = "Lb_MM>6000 && "+ (TString)Form("weight > %e",i*0.1) + baseCut;

		TString hname = Form("hL_%i",i);
		tree->Draw("cosThetaL>>"+hname+"(40,-1,1)",cut,"E");
		TH1F * hL = (TH1F*)gPad->GetPrimitive(hname);
		hL->GetXaxis()->SetTitle("cosThetaL");
		hL->Draw();
		c->Print(Form("AnglesBkg/cosThetaL_sideband_NNcut_%i.pdf",i));

		hname = Form("hB_%i",i);
		tree->Draw("cosThetaB>>"+hname+"(40,-1,1)",cut,"E");
		TH1F * hB = (TH1F*)gPad->GetPrimitive(hname);
		hB->GetXaxis()->SetTitle("cosThetaL");
		hB->Draw();
		c->Print(Form("AnglesBkg/cosThetaB_sideband_NNcut_%i.pdf",i));
	}
}



