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



int main(int argc, char **argv)
{
	bool pythia6 = false;

	TString xvar = "TMath::Power(J_psi_1S_MM/1000,2)";
	TString xvarname = "q2";
	TString yvar = "cosTheta";
	string analysis = "Lb2Lmumu";
	TString type = "All";
	string print = "All";
	string drawopt = "surf";
	bool rel = false;
	
	int xnbins = 9;
	double xbinstmp[] = {0.1, 2.0, 4.0, 6.0, 8.0, 11.0, 15.0, 16.0, 18.0, 20.0};
	double x2binstmp[] = {1.1, 6.0, 15.0, 20.0};
	double * xbins = &xbinstmp[0];
	int ynbins = 0;
	double * ybins  = decodeBinning("[20,-1,1]",&ynbins);
	TCut extCut = "";
	TString outName = "";

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
			if(arg == "-Pythia6") pythia6 = true;
			if(arg.find("-X") != string::npos) { xvar = str; xvarname = str; }
			if(arg.find("-x") != string::npos) 
			{
				if(str=="2bin") {xnbins = 3; xbins = &x2binstmp[0]; }
				else xbins = decodeBinning(str,&xnbins);
			}
			if(arg.find("-y") != string::npos) ybins = decodeBinning(str,&ynbins);
			if(arg.find("-Y") != string::npos) yvar = str;
			if(arg.find("-D") != string::npos) drawopt = str;
			if(arg.find("-O") != string::npos) outName = str;
			if(arg == "-highestq2") extCut = CutsDef::highestq2Cut;
			if(arg == "-highq2") extCut = CutsDef::highq2Cut;
			if(arg == "-lowq2") extCut = CutsDef::lowq2Cut;
			if(arg == "-lowestq2") extCut = CutsDef::lowestq2Cut;
		}
	}

	gStyle->SetPaintTextFormat("5.3f");
	
	// Set trees of cancdidates previously created

	TString mctype = "MC_Pythia8_NBweighted";
	if(pythia6) mctype = "MC_Pythia6_NBweighted"; 
	TString weight = "MCnorm*(lifeTimeW > 0)*(physRate_pol0_noDecay > 0)";
	cout << "Using: " << mctype << endl;
	cout << "Plotting: " << xvar << " vs " << yvar << endl;

	TString nameGeomTree = "MCtree";
	if(xvar=="cosThetaB" || yvar=="cosThetaB") nameGeomTree = "MCtreeDecay";
	TString namefileMCgeom = "/afs/cern.ch/work/p/pluca/Lmumu/weighted/Lb2Lmumu_geom"+mctype+".root";
	TFile * MCfile = TFile::Open(namefileMCgeom);
	TTree * treeMCgeom = (TTree *)MCfile->Get(nameGeomTree);
	TString namefileMC = "/afs/cern.ch/work/p/pluca/Lmumu/weighted/Lb2Lmumu_"+mctype+".root";
	MCfile = TFile::Open(namefileMC);
	TTree * treeMC = (TTree *)MCfile->Get("tree");
	TTree * treeMCGen = (TTree *)MCfile->Get("MCtreeDecay");
	TTree * treeMCAllGen = (TTree *)MCfile->Get("MCtree");
	TString namefileMCjpsi = "/afs/cern.ch/work/p/pluca/Lmumu/weighted/Lb2JpsiL_"+mctype+".root";
	MCfile = TFile::Open(namefileMCjpsi);
	TTree * treeMCjpsi = NULL, * treeMCjpsi_Gen = NULL, * treeMCjpsi_AllGen = NULL, * treeMCjpsi_geom = NULL;
	if(rel)
	{
		treeMCjpsi = (TTree *)MCfile->Get("tree");
		treeMCjpsi_Gen = (TTree *)MCfile->Get("MCtreeDecay");
		treeMCjpsi_AllGen = (TTree *)MCfile->Get("MCtree");
	}
	TString namefileMCjpsiGeom = "/afs/cern.ch/work/p/pluca/Lmumu/weighted/Lb2JpsiL_geom"+mctype+".root";
	MCfile = TFile::Open(namefileMCjpsiGeom);
	if(rel) treeMCjpsi_geom = (TTree *)MCfile->Get("MCtree");

	extCut.Print();	
	TCut baseCut = extCut + CutsDef::mumuTrueID + CutsDef::avoidJpsiCut + CutsDef::massCutUnblinded;
	TCut baseJpsiCut = extCut + CutsDef::jpsiTrueID + CutsDef::jpsiCut + CutsDef::massCutUnblinded;
	if(type == "DD") { baseCut += CutsDef::DDcut; baseJpsiCut += CutsDef::DDcut; }
	else if(type == "LL") { baseCut += CutsDef::LLcut; baseJpsiCut += CutsDef::LLcut; }
	cout << "Analysisng " << type << " events" << endl;
		
	TCut geomCut = "TMath::Abs(TMath::ACos(muplus_TRUEP_Z / TMath::Sqrt(TMath::Power(muplus_TRUEP_Z,2) + TMath::Power(muplus_TRUEP_Y,2) + TMath::Power(muplus_TRUEP_X,2)))) > 0.01 && TMath::Abs(TMath::ACos(muplus_TRUEP_Z / TMath::Sqrt(TMath::Power(muplus_TRUEP_Z,2) + TMath::Power(muplus_TRUEP_Y,2) + TMath::Power(muplus_TRUEP_X,2)))) < 0.4 && TMath::Abs(TMath::ACos(muminus_TRUEP_Z / TMath::Sqrt(TMath::Power(muminus_TRUEP_Z,2) + TMath::Power(muminus_TRUEP_Y,2) + TMath::Power(muminus_TRUEP_X,2)))) > 0.01 && TMath::Abs(TMath::ACos(muminus_TRUEP_Z / TMath::Sqrt(TMath::Power(muminus_TRUEP_Z,2) + TMath::Power(muminus_TRUEP_Y,2) + TMath::Power(muminus_TRUEP_X,2)))) < 0.4";

	
	TFile * histFile = NULL;
	if(outName!="") histFile = new TFile(outName,"recreate");
	else if(rel) histFile = new TFile("Lbrel2Deff_"+yvar+"_vs_"+xvarname+"_"+type+".root","recreate");
	else histFile = new TFile("Lbeff2D_"+yvar+"_vs_"+xvarname+"_"+type+".root","recreate");
	

	/**  Efficiencies  */

	TH2F * geomsys = getEff("GEOsys", "Lb_TRUEP_E", xvar, yvar, xnbins, xbins, ynbins, ybins,
			treeMCgeom, geomCut, treeMCgeom, "",
			weight+"*lifeTimeW*physRate_polp006_noDecay",
			treeMCjpsi_geom, geomCut, treeMCjpsi_geom, "",
			weight+"*lifeTimeW*physRate_pol0_noDecay",
			"-f0.5");

	TH2F * detsys = NULL;
	if(xvar=="cosThetaB" || yvar=="cosThetaB")
		detsys = getEff("DETsys", "Lb_TRUEP_E", xvar, yvar, xnbins, xbins, ynbins, ybins,
			treeMCGen, "", treeMCgeom, geomCut,
			weight+"*lifeTimeW*physRate_polp006_noDecay",
			treeMCjpsi_Gen, "", treeMCjpsi_geom, geomCut,
			weight+"*lifeTimeW*physRate_pol0_noDecay");
	else
		detsys = getEff("DETsys", "Lb_TRUEP_E", xvar, yvar, xnbins, xbins, ynbins, ybins,
			treeMCGen, "", treeMCAllGen, "",
			weight+"*lifeTimeW*physRate_polp006_noDecay",
			treeMCjpsi_Gen, "", treeMCjpsi_AllGen, "",
			weight+"*lifeTimeW*physRate_pol0_noDecay");
		
	TH2F * recosys = getEff("RECOsys", "Lb_TRUEP_E", xvar, yvar, xnbins, xbins, ynbins, ybins,
			treeMC, baseCut, treeMCGen, "",
			weight+"*lifeTimeW*physRate_polp006",
			treeMCjpsi, baseJpsiCut, treeMCjpsi_Gen, "",
			weight+"*lifeTimeW*physRate_pol0");
		
	TH2F * trigsys = getEff("TRIGsys", "Lb_MM", xvar, yvar, xnbins, xbins, ynbins, ybins,
			treeMC, baseCut+CutsDef::TrigPassed, treeMC, baseCut,
			weight+"*lifeTimeW*physRate_polp006",
			treeMCjpsi, baseJpsiCut+CutsDef::TrigPassed, treeMCjpsi_Gen, baseJpsiCut,
			weight+"*lifeTimeW*physRate_pol0");
	
	TH2F * mvasys = getEff("MVAsys", "Lb_MM", xvar, yvar, xnbins, xbins, ynbins, ybins,
			treeMC, baseCut+CutsDef::TrigPassed+CutsDef::MVAcut, treeMC, baseCut+CutsDef::TrigPassed,
			weight+"*lifeTimeW*physRate_polp006",
			treeMCjpsi, baseJpsiCut+CutsDef::TrigPassed+CutsDef::MVAcut, treeMCjpsi_Gen, baseJpsiCut+CutsDef::TrigPassed,
			weight+"*lifeTimeW*physRate_pol0");
	


	//////////////////////////////////////////////////////////////////////////////////////////

	weight += "*Lb_weight";

	TH2F * geomeff = getEff("GEO", "Lb_TRUEP_E", xvar, yvar, xnbins, xbins, ynbins, ybins,
			treeMCgeom, geomCut, treeMCgeom, "",
			weight+"*lifeTimeW*physRate_polp006_noDecay",
			treeMCjpsi_geom, geomCut, treeMCjpsi_geom, "",
			weight+"*lifeTimeW*physRate_pol0_noDecay",
			"-f0.5");

	TH2F * deteff = NULL;
	if(xvar=="cosThetaB" || yvar=="cosThetaB")
		deteff = getEff("DET", "Lb_TRUEP_E", xvar, yvar, xnbins, xbins, ynbins, ybins,
			treeMCGen, "", treeMCgeom, geomCut,
			weight+"*lifeTimeW*physRate_polp006_noDecay",
			treeMCjpsi_Gen, "", treeMCjpsi_geom, geomCut,
			weight+"*lifeTimeW*physRate_pol0_noDecay");
	else
		deteff = getEff("DET", "Lb_TRUEP_E", xvar, yvar, xnbins, xbins, ynbins, ybins,
			treeMCGen, "", treeMCAllGen, "",
			weight+"*lifeTimeW*physRate_polp006_noDecay",
			treeMCjpsi_Gen, "", treeMCjpsi_AllGen, "",
			weight+"*lifeTimeW*physRate_pol0_noDecay");
		
	TH2F * recoeff = getEff("RECO", "Lb_TRUEP_E", xvar, yvar, xnbins, xbins, ynbins, ybins,
			treeMC, baseCut, treeMCGen, "",
			weight+"*lifeTimeW*physRate_polp006",
			treeMCjpsi, baseJpsiCut, treeMCjpsi_Gen, "",
			weight+"*lifeTimeW*physRate_pol0");
		
	TH2F * trigeff = getEff("TRIG", "Lb_MM", xvar, yvar, xnbins, xbins, ynbins, ybins,
			treeMC, baseCut+CutsDef::TrigPassed, treeMC, baseCut,
			weight+"*lifeTimeW*physRate_polp006",
			treeMCjpsi, baseJpsiCut+CutsDef::TrigPassed, treeMCjpsi_Gen, baseJpsiCut,
			weight+"*lifeTimeW*physRate_pol0");
	
	TH2F * mvaeff = getEff("MVA", "Lb_MM", xvar, yvar, xnbins, xbins, ynbins, ybins,
			treeMC, baseCut+CutsDef::TrigPassed+CutsDef::MVAcut, treeMC, baseCut+CutsDef::TrigPassed,
			weight+"*lifeTimeW*physRate_polp006",
			treeMCjpsi, baseJpsiCut+CutsDef::TrigPassed+CutsDef::MVAcut, treeMCjpsi_Gen, baseJpsiCut+CutsDef::TrigPassed,
			weight+"*lifeTimeW*physRate_pol0");
	

	TCanvas * c = new TCanvas();
	gStyle->SetOptStat(0);
	geomeff->Draw(drawopt.c_str());
	c->Print("2Deff_geom_"+yvar+"_vs_"+xvarname+"_"+type+".pdf");
	geomeff->Write("hgeom_eff");
	deteff->Draw(drawopt.c_str());
	c->Print("2Deff_det_"+yvar+"_vs_"+xvarname+"_"+type+".pdf");
	deteff->Write("hdet_eff");
	if(!rel)recoeff->SetMaximum(0.1);
	recoeff->Draw(drawopt.c_str());
	recoeff->Write("hreco_eff");
	c->Print("2Deff_reco_"+yvar+"_vs_"+xvarname+"_"+type+".pdf");
	trigeff->Draw(drawopt.c_str());
	c->Print("2Deff_trig_"+yvar+"_vs_"+xvarname+"_"+type+".pdf");
	trigeff->Write("htrig_eff");
	mvaeff->Draw(drawopt.c_str());
	c->Print("2Deff_mva_"+yvar+"_vs_"+xvarname+"_"+type+".pdf");
	mvaeff->Write("hmva_eff");

	TH2F * uppereff = (TH2F *)recoeff->Clone("toteff");
	uppereff->Multiply(trigeff);
	uppereff->Multiply(mvaeff);
	uppereff->Draw(drawopt.c_str());
	uppereff->SetTitle(0);
	c->Print("2Deff_upper_"+yvar+"_vs_"+xvarname+"_"+type+".pdf");
	uppereff->Write("hupper_eff");

	TH2F * nodeteff = (TH2F *)uppereff->Clone("nodeteff");
	nodeteff->Multiply(geomeff);
	nodeteff->Write("hnodet_eff");

	TH2F * toteff = (TH2F *)geomeff->Clone("toteff");
	toteff->Multiply(deteff);
	toteff->Multiply(uppereff);
	if(!rel) toteff->SetMaximum(0.01);
	toteff->SetTitle("Total eff");
	toteff->Draw(drawopt.c_str());
	toteff->SetTitle(0);
	c->Print("2Deff_tot_"+yvar+"_vs_"+xvarname+"_"+type+".pdf");
	toteff->Write("htot_eff");

	TH2F * totsys = (TH2F *)geomsys->Clone("totsys");
	totsys->Multiply(detsys);
	totsys->Multiply(recosys);
	totsys->Multiply(trigsys);
	totsys->Multiply(mvasys);
	totsys->Draw(drawopt.c_str());
	totsys->SetTitle(0);
	c->Print("2Deff_totsys_"+yvar+"_vs_"+xvarname+"_"+type+".pdf");
	totsys->Write("hsys_eff");


	//histFile->Write();
	delete MCfile;
	delete histFile;
	return 0;
}

