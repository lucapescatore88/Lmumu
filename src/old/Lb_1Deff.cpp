#include <vector>
#include <sstream>
#include <iostream>
#include <string>
#include <time.h>

#include "RooRealVar.h"
#include "RooStats/SPlot.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooCBShape.h"
#include "RooArgusBG.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooDataSet.h"
#include "RooAbsData.h"
#include "RooProdPdf.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "RooArgusBG.h"
#include "RooKeysPdf.h"
#include "RooArgList.h"
#include "RooPolynomial.h"

#include "TSystem.h"
#include "TStyle.h"
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


double toteff(vector<double> effs)
{
	double tot = 1;
	for(unsigned i = 0; i < effs.size(); i++) tot*=effs[i];
	return tot;
}

double toterr(vector<double> effs, vector<double> errs)
{
	double tot = 0;
	if(errs.size()!=effs.size()) {cout << "Eff size different than err size" << endl; return -1;}
	for(unsigned i = 0; i < errs.size(); i++) tot+=TMath::Power(errs[i]/effs[i],2);
	return toteff(effs)*TMath::Sqrt(tot);
}




int main(int argc, char **argv)
{
	bool pythia6 = false;
	TString xvar = "TMath::Power(J_psi_1S_MM/1000,2)";
	string analysis = "Lb2Lmumu";
	string type = "All";
	string print = "All";
	string drawopt = "";
	bool rel = false;

	int xnbins = 8;
	double def_xbins[] = {0.1, 2.0, 4.0, 6.0, 11.0, 15.0, 16.0, 18.0, 20.0};
	double * xbins = &def_xbins[0];
	//double * xbins = new double[xnbins+1];
	//for(int i = 0; i <= xnbins; i++) xbins[i] = -1. + i*2./xnbins;
	
	if(argc > 1)
	{
		for(int a = 1; a < argc; a++)
		{
			string arg = argv[a];
			string str = arg.substr(2,arg.length()-2);
			
			if(arg.find("-t") != string::npos)
			{
				type = str;
				analysis += ("_" + type);
			}
	
			if(arg.find("-p") != string::npos) print = str;
			if(arg == "-r") rel = true;
			if(arg.find("-b") != string::npos ) xbins = decodeBinning(str,&xnbins);
			if(arg.find("-c") != string::npos ) xbins = decodeBinning(str,&xnbins,"custom");
			if(arg == "-Pythia6") pythia6 = true;
			if(arg.find("-X") != string::npos) xvar = (TString)str;
			if(arg.find("-D") != string::npos) drawopt = str;
		}
	}
	
	
	// Set trees of cancdidates previously created

	TString mctype = "MC_Pythia8_NBweighted";
	if(pythia6) mctype = "MC_Pythia6_NBweighted"; 
	TString weight = "MCnorm*(lifeTimeW > 0)*(physRate_pol0_noDecay > 0)";
	cout << "Using: " << mctype << endl;
	cout << "Plotting: " << xvar << endl;
	cout << "Binning: [";
	for(int i = 0; i < xnbins-1; i++) cout << xbins[i] << ",";
	cout << xbins[xnbins-1] << "]" << endl;
	
	TString namefileMCgeom = "/afs/cern.ch/work/p/pluca/weighted/Lmumu/Lb2Lmumu_geom"+mctype+".root";
	TFile * MCfile = TFile::Open(namefileMCgeom);
	TTree * treeMCgeom = (TTree *)MCfile->Get("MCtree");
	TString namefileMC = "/afs/cern.ch/work/p/pluca/weighted/Lmumu/Lb2Lmumu_"+mctype+".root";
	MCfile = TFile::Open(namefileMC);
	TTree * treeMC = (TTree *)MCfile->Get("tree");
	TTree * treeMCGen = (TTree *)MCfile->Get("MCtreeDecay");
	TTree * treeMCAllGen = (TTree *)MCfile->Get("MCtree");
	TString namefileMCjpsi = "/afs/cern.ch/work/p/pluca/weighted/Lmumu/Lb2JpsiL_"+mctype+".root";
	MCfile = TFile::Open(namefileMCjpsi);
	TTree * treeMCjpsi = NULL, * treeMCjpsi_Gen = NULL, * treeMCjpsi_AllGen = NULL, * treeMCjpsi_geom = NULL;
	if(rel)
	{
		treeMCjpsi = (TTree *)MCfile->Get("tree");
		treeMCjpsi_Gen = (TTree *)MCfile->Get("MCtreeDecay");
		treeMCjpsi_AllGen = (TTree *)MCfile->Get("MCtree");
	}
	TString namefileMCjpsiGeom = "/afs/cern.ch/work/p/pluca/weighted/Lmumu/Lb2JpsiL_geom"+mctype+".root";
	MCfile = TFile::Open(namefileMCjpsiGeom);
	if(rel) treeMCjpsi_geom = (TTree *)MCfile->Get("MCtree");
	
	
	TCut baseCut = CutsDef::mumuTrueID + CutsDef::avoidJpsiCut + CutsDef::massCutUnblinded;
	TCut baseJpsiCut = CutsDef::jpsiTrueID + CutsDef::jpsiCut + CutsDef::massCutUnblinded;
	if(type == "DD") { baseCut += (TCut)"pplus_TRACK_Type == 3"; baseJpsiCut += (TCut)"pplus_TRACK_Type == 3"; }
	else if(type == "LL") { baseCut += (TCut)"pplus_TRACK_Type == 5"; baseJpsiCut += (TCut)"pplus_TRACK_Type == 5"; }
	cout << "Analysisng " << type << " events" << endl;
		
	TCut geomCut = "TMath::Abs(TMath::ACos(muplus_TRUEP_Z / TMath::Sqrt(TMath::Power(muplus_TRUEP_Z,2) + TMath::Power(muplus_TRUEP_Y,2) + TMath::Power(muplus_TRUEP_X,2)))) > 0.01 && TMath::Abs(TMath::ACos(muplus_TRUEP_Z / TMath::Sqrt(TMath::Power(muplus_TRUEP_Z,2) + TMath::Power(muplus_TRUEP_Y,2) + TMath::Power(muplus_TRUEP_X,2)))) < 0.4 && TMath::Abs(TMath::ACos(muminus_TRUEP_Z / TMath::Sqrt(TMath::Power(muminus_TRUEP_Z,2) + TMath::Power(muminus_TRUEP_Y,2) + TMath::Power(muminus_TRUEP_X,2)))) > 0.01 && TMath::Abs(TMath::ACos(muminus_TRUEP_Z / TMath::Sqrt(TMath::Power(muminus_TRUEP_Z,2) + TMath::Power(muminus_TRUEP_Y,2) + TMath::Power(muminus_TRUEP_X,2)))) < 0.4";

	
	TFile * histFile = NULL;
	if(rel) histFile = new TFile(("Lbrel1Deff_"+type+".root").c_str(),"recreate");
	else histFile = new TFile(("Lbeff1D_"+type+".root").c_str(),"recreate");
		
	/**  Efficiencies  */

	TH1F * geomeff = getEff("GEO", xvar, xnbins, xbins,
			treeMCgeom, geomCut, treeMCgeom, "",
			weight+"*lifeTimeW*physRate_polp006_noDecay",
			treeMCjpsi_geom, geomCut, treeMCjpsi_geom, "",
			weight+"*lifeTimeW*physRate_pol0_noDecay",
			"-f0.5");
		
	TH1F * deteff = getEff("DET", xvar, xnbins, xbins, 
			treeMCGen, "", treeMCAllGen, "",
			weight+"*lifeTimeW*physRate_polp006_noDecay",
			treeMCjpsi_Gen, "", treeMCjpsi_AllGen, "",
			weight+"*lifeTimeW*physRate_pol0_noDecay");
		
	TH1F * recoeff = getEff("RECO", xvar, xnbins, xbins, 
			treeMC, baseCut, treeMCGen, "",
			weight+"*lifeTimeW*physRate_polp006",
			treeMCjpsi, baseJpsiCut, treeMCjpsi_Gen, "",
			weight+"*lifeTimeW*physRate_pol0");
		
	TH1F * trigeff = getEff("trig", xvar, xnbins, xbins, 
			treeMC, baseCut+CutsDef::TrigPassed, treeMC, baseCut,
			weight+"*lifeTimeW*physRate_polp006",
			treeMCjpsi, baseJpsiCut+CutsDef::TrigPassed, treeMCjpsi_Gen, baseJpsiCut,
			weight+"*lifeTimeW*physRate_pol0");
	
	TH1F * mvaeff = getEff("mva", xvar, xnbins, xbins, 
			treeMC, baseCut+CutsDef::TrigPassed+CutsDef::MVAcut, treeMC, baseCut+CutsDef::TrigPassed,
			weight+"*lifeTimeW*physRate_polp006",
			treeMCjpsi, baseJpsiCut+CutsDef::TrigPassed+CutsDef::MVAcut, treeMCjpsi_Gen, baseJpsiCut+CutsDef::TrigPassed,
			weight+"*lifeTimeW*physRate_pol0");
	
	TH1F * uppereff = getEff("upper", xvar, xnbins, xbins, 
			treeMC, baseCut+CutsDef::TrigPassed+CutsDef::MVAcut, treeMC, baseCut,
			weight+"*lifeTimeW*physRate_polp006",
			treeMCjpsi, baseJpsiCut+CutsDef::TrigPassed+CutsDef::MVAcut, treeMCjpsi_Gen, baseJpsiCut,
			weight+"*lifeTimeW*physRate_pol0");
	
	TCanvas * c = new TCanvas();
	gStyle->SetOptStat(0);
	gStyle->SetOptFit();
	
	geomeff->Fit("pol2");
	c->Print("1Deff_geom.pdf");
	deteff->Fit("pol2");
	c->Print("1Deff_det.pdf");
	if(!rel) recoeff->SetMaximum(0.2);
	recoeff->Fit("pol2");
	c->Print("1Deff_reco.pdf");
	trigeff->Fit("pol2");
	c->Print("1Deff_trig.pdf");
	mvaeff->Fit("pol2");
	c->Print("1Deff_mva.pdf");
	uppereff->Fit("pol2");
	c->Print("1Deff_upper.pdf");

	TH2F * toteff = (TH2F *)geomeff->Clone("toteff");
	toteff->Multiply(deteff);
	toteff->Multiply(recoeff);
	toteff->Multiply(uppereff);
	if(!rel) toteff->SetMaximum(0.02);
	toteff->SetTitle("Total eff");
	toteff->Fit("pol2");
	c->Print("1Deff_tot.pdf");
	
	histFile->Write();
	
	
	delete MCfile;
	delete histFile;
	return 0;
}

