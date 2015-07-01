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

void getAllEffSys(
	TString name, TString xvar, int xnbins, double * xbins, TString weight,
	TTree * treeUp, TCut cutUp, TTree * treeDown, TCut cutDown,
	TTree * treeJpsiUp, TCut cutJpsiUp, TTree * treeJpsiDown, TCut cutJpsiDown,
	vector <TH1F * >  * hdefault, vector <TH1F * >  * hsys, vector <TH1F * >  * hPIDsys,
	vector <TH1F * > * lfsys_plus, vector <TH1F * > * lfsys_minus,
	vector <TH1F * >  * decaysys, vector <TH1F * >  * polsys_minus, vector <TH1F * > * polsys_plus,
	vector <TH1F * >  * poljpsi1, vector <TH1F * > * poljpsi2,
	vector <TH1F * >  * poljpsi3, vector <TH1F * > * poljpsi4,
	vector <TH1F * >  * poljpsi5, vector <TH1F * > * poljpsi6,
	vector <TH1F * >  * poljpsi7, vector <TH1F * > * poljpsi8,
	vector <TH1F * >  * hwilson1, vector <TH1F * > * hwilson2, vector <TH1F * > * hwilson3,
	vector <TH1F * >  * DDsys, bool sys, string opt = "", bool jpsi = false )
{
	TString isDecay = "";
	if( (name == "GEO" && xvar!="cosThetaB") || name == "DET") isDecay = "_noDecay";

	TString polweight = "*physRate_polp006"+isDecay;
	if(jpsi) polweight = "*physRate_pol0"+isDecay;

	hsys->push_back( getEff(name+"sys", xvar, xnbins, xbins,
	treeUp, cutUp, treeDown, cutDown, weight+"*lifeTimeW"+polweight,
	treeJpsiUp, cutJpsiUp, treeJpsiDown, cutJpsiDown, weight+"*lifeTimeW*physRate_pol0"+isDecay,opt) );

	weight+="*Lb_weight";
	hdefault->push_back( getEff(name+"def", xvar, xnbins, xbins,
	treeUp, cutUp, treeDown, cutDown, weight+"*lifeTimeW"+polweight,
	treeJpsiUp, cutJpsiUp, treeJpsiDown, cutJpsiDown, weight+"*lifeTimeW*physRate_pol0"+isDecay,opt) );
	
	if(hPIDsys) hPIDsys->push_back( getEff(name+"PIDsys", xvar, xnbins, xbins,
	treeUp, cutUp, treeDown, cutDown, weight+"*lifeTimeW*pid_muplus*pid_muminus"+polweight,
	treeJpsiUp, cutJpsiUp, treeJpsiDown, cutJpsiDown, weight+"*pid_muplus*pid_muminus*lifeTimeW*physRate_pol0"+isDecay,opt) );

	if(!sys) return;

	lfsys_plus->push_back( getEff(name+"lfp", xvar, xnbins, xbins,
	treeUp, cutUp, treeDown, cutDown, weight+"*lifeTimeW_plussigma*physRate_polp006"+isDecay,
	treeJpsiUp, cutJpsiUp, treeJpsiDown, cutJpsiDown, weight+"*lifeTimeW_plussigma*physRate_pol0"+isDecay,opt) );
	lfsys_minus->push_back( getEff(name+"lfm", xvar, xnbins, xbins,
	treeUp, cutUp, treeDown, cutDown, weight+"*lifeTimeW_minussigma*physRate_polp006"+isDecay,
	treeJpsiUp, cutJpsiUp, treeJpsiDown, cutJpsiDown, weight+"*lifeTimeW_minussigma*physRate_pol0"+isDecay,opt) );	

	polsys_plus->push_back( getEff(name+"pp", xvar, xnbins, xbins,
	treeUp, cutUp, treeDown, cutDown, weight+"*lifeTimeW*physRate_polp015"+isDecay,
	treeJpsiUp, cutJpsiUp, treeJpsiDown, cutJpsiDown, weight+"*lifeTimeW*physRate_pol0"+isDecay,opt) );
	polsys_minus->push_back( getEff(name+"pm", xvar, xnbins, xbins,
	treeUp, cutUp, treeDown, cutDown, weight+"*lifeTimeW*physRate_polm003"+isDecay,
	treeJpsiUp, cutJpsiUp, treeJpsiDown, cutJpsiDown, weight+"*lifeTimeW*physRate_pol0"+isDecay,opt) );

	poljpsi1->push_back( getEff(name+"j1", xvar, xnbins, xbins,
	treeUp, cutUp, treeDown, cutDown, weight+"*lifeTimeW*physRate_polp006"+isDecay,
	treeJpsiUp, cutJpsiUp, treeJpsiDown, cutJpsiDown, weight+"*lifeTimeW*model_jpsi_1"+isDecay,opt) );
	poljpsi2->push_back( getEff(name+"j2", xvar, xnbins, xbins,
	treeUp, cutUp, treeDown, cutDown, weight+"*lifeTimeW*physRate_polp006"+isDecay,
	treeJpsiUp, cutJpsiUp, treeJpsiDown, cutJpsiDown, weight+"*lifeTimeW*model_jpsi_2"+isDecay,opt) );
	poljpsi3->push_back( getEff(name+"j3", xvar, xnbins, xbins,
	treeUp, cutUp, treeDown, cutDown, weight+"*lifeTimeW*physRate_polp006"+isDecay,
	treeJpsiUp, cutJpsiUp, treeJpsiDown, cutJpsiDown, weight+"*lifeTimeW*model_jpsi_3"+isDecay,opt) );
	poljpsi4->push_back( getEff(name+"j4", xvar, xnbins, xbins,
	treeUp, cutUp, treeDown, cutDown, weight+"*lifeTimeW*physRate_polp006"+isDecay,
	treeJpsiUp, cutJpsiUp, treeJpsiDown, cutJpsiDown, weight+"*lifeTimeW*model_jpsi_4"+isDecay,opt) );
	poljpsi5->push_back( getEff(name+"j5", xvar, xnbins, xbins,
	treeUp, cutUp, treeDown, cutDown, weight+"*lifeTimeW*physRate_polp006"+isDecay,
	treeJpsiUp, cutJpsiUp, treeJpsiDown, cutJpsiDown, weight+"*lifeTimeW*model_jpsi_5"+isDecay,opt) );
	poljpsi6->push_back( getEff(name+"j6", xvar, xnbins, xbins,
	treeUp, cutUp, treeDown, cutDown, weight+"*lifeTimeW*physRate_polp006"+isDecay,
	treeJpsiUp, cutJpsiUp, treeJpsiDown, cutJpsiDown, weight+"*lifeTimeW*model_jpsi_6"+isDecay,opt) );
	poljpsi7->push_back( getEff(name+"j7", xvar, xnbins, xbins,
	treeUp, cutUp, treeDown, cutDown, weight+"*lifeTimeW*physRate_polp006"+isDecay,
	treeJpsiUp, cutJpsiUp, treeJpsiDown, cutJpsiDown, weight+"*lifeTimeW*model_jpsi_7"+isDecay,opt) );
	poljpsi8->push_back( getEff(name+"j8", xvar, xnbins, xbins,
	treeUp, cutUp, treeDown, cutDown, weight+"*lifeTimeW*physRate_polp006"+isDecay,
	treeJpsiUp, cutJpsiUp, treeJpsiDown, cutJpsiDown, weight+"*lifeTimeW*model_jpsi_8"+isDecay,opt) );

	hwilson1->push_back( getEff(name+"wil1", xvar, xnbins, xbins,
	treeUp, cutUp, treeDown, cutDown, weight+"*lifeTimeW*physRate_pol0_wilson1"+isDecay,
	treeJpsiUp, cutJpsiUp, treeJpsiDown, cutJpsiDown, weight+"*lifeTimeW*physRate_pol0"+isDecay,opt) );
	hwilson2->push_back( getEff(name+"wil2", xvar, xnbins, xbins,
	treeUp, cutUp, treeDown, cutDown, weight+"*lifeTimeW*physRate_pol0_wilson2"+isDecay,
	treeJpsiUp, cutJpsiUp, treeJpsiDown, cutJpsiDown, weight+"*lifeTimeW*physRate_pol0"+isDecay,opt) );
	hwilson3->push_back( getEff(name+"wil3", xvar, xnbins, xbins,
	treeUp, cutUp, treeDown, cutDown, weight+"*lifeTimeW*physRate_pol0_wilson3"+isDecay,
	treeJpsiUp, cutJpsiUp, treeJpsiDown, cutJpsiDown, weight+"*lifeTimeW*physRate_pol0"+isDecay,opt) );

	decaysys->push_back( getEff(name+"decay", xvar, xnbins, xbins,
	treeUp, cutUp, treeDown, cutDown, weight+"*lifeTimeW*physRate_pol0_QCDff"+isDecay,
	treeJpsiUp, cutJpsiUp, treeJpsiDown, cutJpsiDown, weight+"*lifeTimeW*physRate_pol0"+isDecay,opt) );	
	
	if(name == "RECO")
		DDsys->push_back( getEff(name+"DD", xvar, xnbins, xbins,
		treeUp, TCut ( "DDvtx_weight*( " + (TString)cutUp + " )"), treeDown, cutDown, weight+"*lifeTimeW*physRate_polp006",
		treeJpsiUp, TCut ( "DDvtx_weight*( " + (TString)cutJpsiUp + " )"), treeJpsiDown, cutJpsiDown, weight+"*lifeTimeW*physRate_pol0") );
	else if(name == "TRIG" || name == "MVA")
		DDsys->push_back( getEff(name+"DD", xvar, xnbins, xbins,
		treeUp, cutUp, treeDown, cutDown, weight+"*lifeTimeW*physRate_polp006*DDvtx_weight",
		treeJpsiUp, cutJpsiUp, treeJpsiDown, cutJpsiDown, weight+"*lifeTimeW*physRate_pol0*DDvtx_weight") );
	else DDsys->push_back( getEff(name+"DD", xvar, xnbins, xbins,
		 treeUp, cutUp, treeDown, cutDown, weight+"*lifeTimeW*physRate_polp006_noDecay",
		 treeJpsiUp, cutJpsiUp, treeJpsiDown, cutJpsiDown, weight+"*lifeTimeW*physRate_pol0_noDecay",opt) );
}



double createSys(int j, TH1F * hdefault, TH1F * hminus, TH1F * hplus = NULL, double sum = 0)
{
	double def = hdefault->GetBinContent(j);
	double minus = def;
	double plus = def;
	if(minus) minus = hminus->GetBinContent(j);
	if(plus) plus = hplus->GetBinContent(j);
	double sys = TMath::Max(TMath::Abs(def-plus),TMath::Abs(def-minus)) / def;	
		
	if(sum) sys = TMath::Sqrt( TMath::Power(sum,2) + TMath::Power(sys,2));
	
	return sys;
}


TH1F * sqSum(TH1F * h1, TH1F * h2)
{
	TH1F * res = (TH1F *) h2->Clone("sum");
	res->Reset();
	
	for(int j = 1; j <= res->GetNbinsX(); j++)
	{
		if(h1) res->SetBinContent(j,TMath::Sqrt( TMath::Power(h1->GetBinContent(j),2) + TMath::Power(h2->GetBinContent(j),2)));
		else res->SetBinContent(j,h2->GetBinContent(j));
	}
	
	return res;
}


TH1F * getErrHist(TH1F * h, bool rel = true)
{
	TH1F * res = (TH1F *) h->Clone("errhist");
	res->Reset();
	
	for(int j = 1; j <= res->GetNbinsX(); j++)
		if(rel) res->SetBinContent(j,h->GetBinError(j)/h->GetBinContent(j));
		else res->SetBinContent(j,h->GetBinError(j));

	return res;
}



int main(int argc, char **argv)
{
	TString xvar = "TMath::Power(J_psi_1S_MM/1000,2)";
	TString xvarname = "q2";
	TString type = "All";
	TCut extCut = "";
	
	string drawopt = "";
	TString outFileName = "";
	bool rel = false, doSys = false, percent = false, jpsi = false, pythia6 = false;
	
	int xnbins = 12;
	double def_xbins[] = {0.1, 2.0, 4.0, 6.0, 8.0, 9.1, 10.1, 11.0, 12.5, 15.0, 16.0, 18.0, 20.0};
	//double def_xbins[] = {1.1, 6.0, 15.0, 20.0};
	double * xbins = &def_xbins[0];
	
	if(argc > 1)
	{
		for(int a = 1; a < argc; a++)
		{
			string arg = argv[a];
			string str = arg.substr(2,arg.length()-2);
			
			if(arg.find("-t") != string::npos) type = (TString)str;
			if(arg == "-r") rel = true;
			if(arg.find("-b") != string::npos ) xbins = decodeBinning(str,&xnbins);
			if(arg.find("-c") != string::npos ) xbins = decodeBinning(str,&xnbins,"custom");
			if(arg.find("-C") != string::npos ) extCut = (TCut)((TString)str);
			if(arg == "-highq2Cut") extCut = (TCut)"( TMath::Power(J_psi_1S_MM/1000,2) > 15. && TMath::Power(J_psi_1S_MM/1000,2) < 20. )";
			if(arg == "-Pythia6") pythia6 = true;
			if(arg.find("-X") != string::npos) { xvar = str; xvarname = xvar; }
			if(arg.find("-D") != string::npos) drawopt = str;
			if(arg.find("-sys") != string::npos) doSys = true;
			if(arg.find("-percent") != string::npos) percent = true;
			if(arg.find("-o") != string::npos) outFileName = str;
			if(arg.find("-jpsi") != string::npos) jpsi = true;
		}
	}
	
	if(jpsi) { if(xvarname=="q2") xbins = decodeBinning("[1,0,25]",&xnbins); rel = false; doSys = false; }

	// Set trees of cancdidates previously created

	TString mctype = "MC_Pythia8_NBweighted";
	if(pythia6) mctype = "MC_Pythia6_NBweighted"; 
	TString weight = "MCnorm*(lifeTimeW > 0)*(physRate_pol0_noDecay > 0)";
	cout << "Using: " << mctype << endl;
	cout << "Plotting: " << xvar << endl;
	cout << "Binning: [";
	for(int i = 0; i < xnbins; i++) cout << xbins[i] << ",";
	cout << xbins[xnbins] << "]" << endl;	

	TString namefileMCgeom = "/afs/cern.ch/work/p/pluca/Lmumu/weighted/Lb2Lmumu_geom"+mctype+".root";
	TString namefileMC = "/afs/cern.ch/work/p/pluca/Lmumu/weighted/Lb2Lmumu_"+mctype+".root";
	TString namefileMCjpsi = "/afs/cern.ch/work/p/pluca/Lmumu/weighted/Lb2JpsiL_"+mctype+".root";
	TString namefileMCjpsiGeom = "/afs/cern.ch/work/p/pluca/Lmumu/weighted/Lb2JpsiL_geom"+mctype+".root";
	TString namefileDataJpsi = "/afs/cern.ch/work/p/pluca/Lmumu/weighted/Lb2Lmumu_CL_NBweighted.root";

	TString nameGeomTree = "MCtree";
	if(xvar=="cosThetaB") nameGeomTree = "MCtreeDecay";
	TFile * MCfile = TFile::Open(namefileMCgeom);
	TTree * treeMCgeom = (TTree *)MCfile->Get(nameGeomTree);
	MCfile = TFile::Open(namefileMC);
	TTree * treeMC = (TTree *)MCfile->Get("tree");
	TTree * treeMCGen = (TTree *)MCfile->Get("MCtreeDecay");
	TTree * treeMCAllGen = (TTree *)MCfile->Get("MCtree");
	MCfile = TFile::Open("/afs/cern.ch/work/p/pluca/Lmumu/weighted/trainingSamples.root");
	TTree * treeMCmva = (TTree *)MCfile->Get("sigTestSample");

	TTree * treeMCjpsi = NULL, * treeMCjpsi_Gen = NULL, * treeMCjpsi_AllGen = NULL, * treeMCjpsi_geom = NULL;
	if(rel || jpsi)
	{
		
		MCfile = TFile::Open(namefileMCjpsi);
		if(jpsi) {
		treeMC = (TTree *)MCfile->Get("tree");
		treeMCmva = (TTree *)MCfile->Get("tree");
		treeMCGen = (TTree *)MCfile->Get("MCtreeDecay");
		treeMCAllGen = (TTree *)MCfile->Get("MCtree");	}
		else {
		treeMCjpsi = (TTree *)MCfile->Get("tree");
		treeMCjpsi_Gen = (TTree *)MCfile->Get("MCtreeDecay");
		treeMCjpsi_AllGen = (TTree *)MCfile->Get("MCtree");	}
		
		MCfile = TFile::Open(namefileMCjpsiGeom);		
		if(jpsi) treeMCgeom = (TTree *)MCfile->Get(nameGeomTree);
		else treeMCjpsi_geom = (TTree *)MCfile->Get(nameGeomTree);
	}
	
	TCut geomCut = CutsDef::geomCut;
	TCut baseCut = extCut + CutsDef::baseCutMuMu;
	if(jpsi) baseCut = extCut + CutsDef::baseCutJpsi;
	TCut baseJpsiCut = extCut + CutsDef::baseCutJpsi;
	TCut binCut = "TMath::Power(J_psi_1S_MM/1000,2) > 9.1 && TMath::Power(J_psi_1S_MM/1000,2) < 10.1";
	if(type == "DD") { baseCut += CutsDef::DDcut; baseJpsiCut += CutsDef::DDcut; }
	else if(type == "LL") { baseCut += CutsDef::LLcut; baseJpsiCut += CutsDef::LLcut; }
	cout << "Analysisng " << type << " events" << endl;

	TString myName = "Lbeff";
	if(rel) myName = "Lbreleff";
	if(doSys) myName += "AndSys";
	if(jpsi) myName += "_Jpsi"; 
	myName += ("vs"+xvarname+"_"+type+".root");
	if(outFileName!="") myName = outFileName;
	TFile * histFile = new TFile(myName,"recreate");

	/**      Calc efficiencies and systematics        */

	gROOT->ProcessLine(".x ~/work/lhcbStyle.C");

	vector <TString> effnames;
	effnames.push_back("geom");
	effnames.push_back("det");
	effnames.push_back("reco");
	effnames.push_back("mva");
	effnames.push_back("trig");

	vector <TString> complnames;
	complnames.push_back("Geometric");
	complnames.push_back("Detection");
	complnames.push_back("Reconstruction");
	complnames.push_back("MVA");
	complnames.push_back("Trigger");

	vector <TH1F * >  hdefault, hsys, hPIDsys;
	vector <TH1F * >  lfsys_plus, lfsys_minus;
	vector <TH1F * >  decaysys, DDsys;
	vector <TH1F * >  polsys_minus, polsys_plus;
	vector <TH1F * >  poljpsi1, poljpsi2, poljpsi3, poljpsi4, poljpsi5, poljpsi6, poljpsi7, poljpsi8;
	vector <TH1F * >  hwilson1, hwilson2, hwilson3;
	
	cout << "Analysing GEO sys" << endl;
	getAllEffSys("GEO", xvar, xnbins, xbins, weight,
			treeMCgeom, geomCut+extCut, treeMCgeom, extCut,
			treeMCjpsi_geom, geomCut+extCut+binCut, treeMCjpsi_geom, extCut+binCut,
			&hdefault, &hsys, NULL, &lfsys_plus, &lfsys_minus, &decaysys, &polsys_minus, &polsys_plus,
			&poljpsi1, &poljpsi2, &poljpsi3, &poljpsi4, &poljpsi5, &poljpsi6, &poljpsi7, &poljpsi8,
			&hwilson1, &hwilson2, &hwilson3,
			&DDsys, doSys, "-f0.5", jpsi);

	cout << "Analysing DET sys" << endl;
	if(xvar=="cosThetaB")
		getAllEffSys("DET", xvar, xnbins, xbins, weight,
			treeMCGen, extCut, treeMCgeom, geomCut+extCut,
			treeMCjpsi_Gen, extCut+binCut, treeMCjpsi_AllGen, geomCut+extCut+binCut,
			&hdefault, &hsys, NULL, &lfsys_plus, &lfsys_minus, &decaysys, &polsys_minus, &polsys_plus,
			&poljpsi1, &poljpsi2, &poljpsi3, &poljpsi4, &poljpsi5, &poljpsi6, &poljpsi7, &poljpsi8,
			&hwilson1, &hwilson2, &hwilson3,
			&DDsys, doSys, "", jpsi );	
	else
		getAllEffSys("DET", xvar, xnbins, xbins, weight,
			treeMCGen, extCut, treeMCAllGen, extCut,
			treeMCjpsi_Gen, extCut+binCut, treeMCjpsi_AllGen, extCut+binCut,
			&hdefault, &hsys, NULL, &lfsys_plus, &lfsys_minus, &decaysys, &polsys_minus, &polsys_plus,
			&poljpsi1, &poljpsi2, &poljpsi3, &poljpsi4, &poljpsi5, &poljpsi6, &poljpsi7, &poljpsi8,
			&hwilson1, &hwilson2, &hwilson3,
			&DDsys, doSys, "", jpsi );
	
	cout << "Analysing RECO sys" << endl;
	getAllEffSys("RECO", xvar, xnbins, xbins, weight,
			treeMC, baseCut, treeMCGen, extCut,
			treeMCjpsi, baseJpsiCut+binCut, treeMCjpsi_Gen, extCut+binCut,
			&hdefault, &hsys, NULL, &lfsys_plus, &lfsys_minus, &decaysys, &polsys_minus, &polsys_plus,
			&poljpsi1, &poljpsi2, &poljpsi3, &poljpsi4, &poljpsi5, &poljpsi6, &poljpsi7, &poljpsi8,
			&hwilson1, &hwilson2, &hwilson3,
			&DDsys, doSys, "", jpsi );
	
	cout << "Analysing MVA sys" << endl;
	getAllEffSys("MVA", xvar, xnbins, xbins, weight,
			treeMCmva, baseCut+CutsDef::MVAcut, treeMCmva, baseCut,
			treeMCjpsi, baseJpsiCut+CutsDef::MVAcut+binCut, treeMCjpsi, baseJpsiCut+binCut,
			&hdefault, &hsys, &hPIDsys, &lfsys_plus, &lfsys_minus, &decaysys, &polsys_minus, &polsys_plus,
			&poljpsi1, &poljpsi2, &poljpsi3, &poljpsi4, &poljpsi5, &poljpsi6, &poljpsi7, &poljpsi8,
			&hwilson1, &hwilson2, &hwilson3,
			&DDsys, doSys, "", jpsi );

	cout << "Analysing TRIG sys" << endl;
	getAllEffSys("TRIG", xvar, xnbins, xbins, weight,
			treeMC, baseCut+CutsDef::MVAcut+CutsDef::TrigPassed, treeMC, baseCut+CutsDef::MVAcut,
			treeMCjpsi, baseJpsiCut+CutsDef::MVAcut+CutsDef::TrigPassed+binCut, treeMCjpsi, baseJpsiCut+CutsDef::MVAcut+binCut,
			&hdefault, &hsys, &hPIDsys, &lfsys_plus, &lfsys_minus, &decaysys, &polsys_minus, &polsys_plus,
			&poljpsi1, &poljpsi2, &poljpsi3, &poljpsi4, &poljpsi5, &poljpsi6, &poljpsi7, &poljpsi8,
			&hwilson1, &hwilson2, &hwilson3,
			&DDsys, doSys, "", jpsi );
	
	TH1F * toteff_lowSel = NULL;
	if(rel || jpsi)
	{
		TString polweight = "physRate_polp006";
		if(jpsi) polweight = "physRate_pol0";
		TH1F * mva_lowSel = getEff("UPPER", xvar, xnbins, xbins,
			treeMCmva, baseCut+CutsDef::MVAcut, treeMCmva, baseCut,
			weight+"*lifeTimeW*Lb_weight*"+polweight,
			treeMCjpsi, baseJpsiCut+CutsDef::MVAcut_lowSel+binCut, treeMCjpsi, baseJpsiCut+binCut,
			weight+"*lifeTimeW*physRate_pol0*pt_weight" );
		mva_lowSel->Write("hmvaeff_lowSel");
		TH1F * uppereff_lowSel = (TH1F *)hdefault[2]->Clone("huppereff_lowSel");
		uppereff_lowSel->Multiply(hdefault[3]);
		uppereff_lowSel->Multiply(mva_lowSel);
		toteff_lowSel = (TH1F *)hdefault[0]->Clone("htoteff_lowSel");
		toteff_lowSel->Multiply(hdefault[1]);
		toteff_lowSel->Multiply(uppereff_lowSel);
		uppereff_lowSel->Write("huppereff_lowSel");
		toteff_lowSel->Write("htoteff_lowSel");

		if(jpsi)
		{
			TH1F * pid = getEff("MCPID", xvar, xnbins, xbins,
			treeMC, baseCut+"pplus_PIDp > -5", treeMC, baseCut,
			weight+"*lifeTimeW*Lb_weight*"+polweight);
			cout << pid->GetBinCenter(1) << "   " << pid->GetBinContent(1) << " +/- " << pid->GetBinError(1) << endl;
		}
	}

			
	TCanvas * c = new TCanvas();
	gStyle->SetOptStat(0);
	gStyle->SetOptFit();
	TH1F * uppereff = (TH1F *)hdefault[2]->Clone("huppereff");
	uppereff->Multiply(hdefault[3]);
	uppereff->Multiply(hdefault[4]);
	TH1F * toteff = (TH1F *)hdefault[0]->Clone("htoteff");
	if(xvar=="cosThetaB") hdefault[1]->Scale(1./hdefault[1]->Integral());
	toteff->Multiply(hdefault[1]);
	toteff->Multiply(uppereff);
	toteff->SetTitle("Total eff");
	toteff->Draw();
	//toteff->Fit("pol2");
	c->Print("effvs"+xvarname+"_"+type+"_tot.pdf");
	toteff->Write("htoteff");
	//uppereff->Fit("pol2");
	uppereff->Draw();
	c->Print("effvs"+xvarname+"_"+type+"_upper.pdf");
	uppereff->Write("huppereff");
	TH1F * tot_nodet_eff = (TH1F *)hdefault[0]->Clone("htot_nodet_eff");
	tot_nodet_eff->Multiply(uppereff);
	tot_nodet_eff->Write("htot_nodet_eff");

	TH1F * syseff = (TH1F *)hsys[0]->Clone("hsyseff");
	syseff->Multiply(hsys[1]);
	syseff->Multiply(hsys[2]);
	syseff->Multiply(hsys[3]);
	syseff->Multiply(hsys[4]);
	syseff->Write("sys_eff");

	/*
	TH1F * Kinsys = (TH1F *)syseff->Clone("Kinsys");
	Kinsys->Add(toteff,-1);
	Kinsys->GetXaxis()->SetTitle("q^{2} (GeV^{2})");
	Kinsys->GetYaxis()->SetTitle("(eff^{now} - eff^{w})");// / eff^{w}");
	Kinsys->Draw();
	c->Print("Kin_sys.pdf");
	
	TH1F * PIDsyseff = (TH1F *)hPIDsys[0]->Clone("hPIDsyseff");
	PIDsyseff->Multiply(hPIDsys[1]);

	TH1F * PIDdefeff = (TH1F *)hdefault[3]->Clone("hdefPID");
	PIDdefeff->Multiply(hdefault[4]);
	
	TH1F * PIDsys = (TH1F *)PIDsyseff->Clone("PIDsys");
	PIDsys->Add(PIDdefeff,-1);
	//PIDsys->Divide(toteff,1,"B");
	PIDsys->Write("PIDmuSys");
	PIDsys->GetXaxis()->SetTitle("q^{2} (GeV^{2})");
	PIDsys->GetYaxis()->SetTitle("(eff^{now} - eff^{w})");// / eff^{w}");
	PIDsys->Draw();
	c->Print("PIDmu_sys.pdf");
*/

	for(unsigned i = 0; i < effnames.size(); i++)
	{
		//hdefault[i]->Fit("pol2");
		if(!jpsi)
		for(int b = 0; b < hdefault[i]->GetNbinsX(); b++)
		{
			hdefault[i]->SetBinContent(hdefault[i]->GetXaxis()->FindBin(8.5),0);
			hdefault[i]->SetBinError(hdefault[i]->GetXaxis()->FindBin(8.5),0);
			hdefault[i]->SetBinContent(hdefault[i]->GetXaxis()->FindBin(10.5),0);
			hdefault[i]->SetBinError(hdefault[i]->GetXaxis()->FindBin(10.5),0);
		}
		if(rel) hdefault[i]->SetMinimum(0.5);
		if(rel) hdefault[i]->SetMaximum(1.5);
		hdefault[i]->GetXaxis()->SetTitle("q^{2} [GeV^{2}/#it{c}^{4}]");
		hdefault[i]->GetYaxis()->SetTitle("Efficiency");
		hdefault[i]->SetTitle(complnames[i]);
		hdefault[i]->Draw();
		c->Print("effvs"+xvarname+"_"+type+"_"+effnames[i]+".pdf");
		hdefault[i]->Write("h"+effnames[i]+"eff");
	}

	gStyle->SetOptFit(0);



	/** Printing out efficiencies and systematics **/

	vector <TString> sysnames;
	sysnames.push_back("Lifetime");
	sysnames.push_back("Decay Model");
	sysnames.push_back("Polarization");
	if(type=="DD") sysnames.push_back("DD vtx");
	
	/** Print efficiencies */ 

	cout << "\n\n" << xvarname << " bin " << " \t\t\t& ";
	for(unsigned s = 0; s < effnames.size(); s++) cout << effnames[s] << " \t\t\t\t& ";
	cout << "Upper \t\t\t\t& Total  \\\\" << endl;

	TGraphErrors * grtoteff = new TGraphErrors();
	TGraphErrors * grtoteff_lowSel = new TGraphErrors();
	for(int j = 1; j <= toteff->GetNbinsX(); j++)
	{
		if((xbins[j]==11 && !rel) || xbins[j]==15) continue;
		if(xbins[j]==9.1 || xbins[j-1]==10.1) continue;
		cout << fixed << setprecision(1) << "eff " << xbins[j-1] << "-" << xbins[j] << fixed << setprecision(5) << " \t & ";

		for(unsigned i = 0; i < hdefault.size(); i++)	
			cout << "$" << hdefault[i]->GetBinContent(j) << " \\pm " << hdefault[i]->GetBinError(j) << "$ \t & ";

		cout << "$" << uppereff->GetBinContent(j) << " \\pm " << uppereff->GetBinError(j) << "$ \t & ";
		cout << "$" << toteff->GetBinContent(j) << " \\pm " << toteff->GetBinError(j) << "$ \\\\ " << endl;

		grtoteff->SetPoint(j,toteff->GetBinCenter(j),toteff->GetBinContent(j));
		grtoteff->SetPointError(j,toteff->GetBinWidth(j)/2.,toteff->GetBinError(j));
		if(toteff_lowSel) grtoteff_lowSel->SetPoint(j,toteff_lowSel->GetBinCenter(j),toteff_lowSel->GetBinContent(j));
		if(toteff_lowSel) grtoteff_lowSel->SetPointError(j,toteff_lowSel->GetBinWidth(j)/2.,toteff_lowSel->GetBinError(j));
	}

	grtoteff->Write("toteff");
	if(toteff_lowSel) grtoteff_lowSel->Write("toteff_lowSel");
	if(!doSys) { delete MCfile; delete histFile; return 0; }


	/** Print sys separate in efficiency */
	
	vector < TH1F * > tmp;
	vector < vector < TH1F * > > sys_eff(4,tmp);
	vector < TH1F * > tot_sys_eff;
	vector < TH1F * > wilson1_sys_eff, wilson2_sys_eff, wilson3_sys_eff;

	cout << endl << endl << endl;
	if(!percent) cout << "\n\n" << xvarname << " bin\t\t& Value \t & Stats";
	else cout << "\n\n" << xvarname << " bin ";
	for(unsigned s = 0; s < sysnames.size(); s++) cout << "\t& " << sysnames[s];
	cout << " \\\\" << endl;

	for(int j = 1; j <= toteff->GetNbinsX(); j++)
	{
		if((xbins[j]==10.1 && !rel) || xbins[j]==15 || xbins[j]==9.1 || xbins[j]==11) continue;

		cout << "-----------------------------------------------------------------------------------------" << endl;
		cout << fixed << setprecision(1) << xbins[j-1] << "-" << xbins[j] << fixed << setprecision(3) << endl;
		cout << "-----------------------------------------------------------------------------------------" << endl;
		if(!percent) cout << fixed << setprecision(5);

		for(unsigned i = 0; i < effnames.size(); i++)
		{
			if(j==1)
			{
				sys_eff[0].push_back((TH1F*)toteff->Clone("sys_lf_"+effnames[i]));
				sys_eff[1].push_back((TH1F*)toteff->Clone("sys_decay_"+effnames[i]));
				sys_eff[2].push_back((TH1F*)toteff->Clone("sys_pol_"+effnames[i]));
				sys_eff[3].push_back((TH1F*)toteff->Clone("sys_DD_"+effnames[i]));
				tot_sys_eff.push_back((TH1F*)toteff->Clone("tot_sys_eff"));
				wilson1_sys_eff.push_back((TH1F*)toteff->Clone("wilson1_sys"));
				wilson2_sys_eff.push_back((TH1F*)toteff->Clone("wilson2_sys"));
				wilson3_sys_eff.push_back((TH1F*)toteff->Clone("wilson3_sys"));
				sys_eff[0][i]->Reset();
				sys_eff[1][i]->Reset();
				sys_eff[2][i]->Reset();
				sys_eff[3][i]->Reset();
				tot_sys_eff[i]->Reset();
				wilson1_sys_eff[i]->Reset();
				wilson2_sys_eff[i]->Reset();
				wilson3_sys_eff[i]->Reset();
			}

			double lf_sys = createSys(j,hdefault[i], lfsys_minus[i], lfsys_plus[i]);
			double decay_sys = createSys(j,hdefault[i], decaysys[i], decaysys[i]);
			double DD_sys = 0;
			double pol_sys = createSys(j, hdefault[i], polsys_minus[i], polsys_plus[i]);
			pol_sys = createSys(j, hdefault[i], poljpsi1[i], poljpsi2[i], pol_sys);
			pol_sys = createSys(j, hdefault[i], poljpsi3[i], poljpsi4[i], pol_sys);
			pol_sys = createSys(j, hdefault[i], poljpsi5[i], poljpsi6[i], pol_sys);
			pol_sys = createSys(j, hdefault[i], poljpsi7[i], poljpsi8[i], pol_sys);
			double cureff = hdefault[i]->GetBinContent(j);
			double curerr = hdefault[i]->GetBinError(j);
			double tot_eff_sys = TMath::Sqrt( TMath::Power(lf_sys,2) + TMath::Power(pol_sys,2) + TMath::Power(decay_sys,2) );

			double wilson_sys1 = createSys(j, hdefault[i], hwilson1[i], hwilson1[i]);
			double wilson_sys2 = createSys(j, hdefault[i], hwilson2[i], hwilson2[i]);
			double wilson_sys3 = createSys(j, hdefault[i], hwilson3[i], hwilson3[i]);
			
			if(type=="DD")
			{
				DD_sys = createSys(j, hdefault[i], DDsys[i], DDsys[i]);
				tot_eff_sys = TMath::Sqrt( TMath::Power(tot_eff_sys,2) + TMath::Power(DD_sys,2) );
				sys_eff[3][i]->SetBinContent(j,DD_sys);
			}

			cout << effnames[i] << " \t  & ";
			if(percent)
			{
				cout << lf_sys*100 << "\\% \t & " << decay_sys*100 << "\\% \t & " << pol_sys*100 << "\\% \t";
				if(type=="DD") cout << " & " << DD_sys*100 << " \t ";	
			}	
			else
			{
				cout << cureff << " \t & " << curerr << " \t & ";
				cout << lf_sys*cureff << " \t & " << decay_sys*cureff << " \t & " << pol_sys*cureff << " \t";
				if(type=="DD") cout << " & " << DD_sys*cureff << " \t";	
			}	
			cout << " \\\\ " << endl;
			sys_eff[0][i]->SetBinContent(j,lf_sys);
			sys_eff[1][i]->SetBinContent(j,decay_sys);
			sys_eff[2][i]->SetBinContent(j,pol_sys);
			tot_sys_eff[i]->SetBinContent(j,tot_eff_sys);
			wilson1_sys_eff[i]->SetBinContent(j,wilson_sys1);
			wilson2_sys_eff[i]->SetBinContent(j,wilson_sys2);
			wilson3_sys_eff[i]->SetBinContent(j,wilson_sys3);
		}
	}


		/** Print total sys */

	vector< TH1F * > sys;
	TH1F * tot_sys = getErrHist(toteff);
	
	for(unsigned s = 0; s < sysnames.size(); s++)
	{
		TH1F * tmp_sys = NULL;
				
		for(unsigned i = 0; i < effnames.size(); i++)
			tmp_sys = sqSum(tmp_sys,sys_eff[s][i]); 
	
		sys.push_back( tmp_sys );
		tmp_sys->Write("sys_"+sysnames[s]);
		tot_sys = sqSum(tot_sys,tmp_sys);
	}
	tot_sys->Write("sys_tot");

	TH1F * wilson_sys1 = NULL, * wilson_sys2 = NULL, * wilson_sys3 = NULL;
	for(unsigned i = 0; i < effnames.size(); i++)
	{
		wilson_sys1 = sqSum(wilson_sys1,wilson1_sys_eff[i]); 
		wilson_sys2 = sqSum(wilson_sys2,wilson2_sys_eff[i]);
		wilson_sys3 = sqSum(wilson_sys3,wilson3_sys_eff[i]);
	}

	wilson_sys1->SetBinContent(wilson_sys1->GetXaxis()->FindBin(8.5),0);
	wilson_sys1->SetBinError(wilson_sys1->GetXaxis()->FindBin(8.5),0);
	wilson_sys1->SetBinContent(wilson_sys1->GetXaxis()->FindBin(10.5),0);
	wilson_sys1->SetBinError(wilson_sys1->GetXaxis()->FindBin(10.5),0);
	wilson_sys1->SetBinContent(wilson_sys1->GetXaxis()->FindBin(13.),0);
	wilson_sys1->SetBinError(wilson_sys1->GetXaxis()->FindBin(13.),0);

	wilson_sys1->Write("rel_wilson_sys1");
	wilson_sys1->SetTitle("Relative Wilson Coeff variation sys (CNP[7]=-0.05, CNP[9]=-1.3)");
	wilson_sys1->Draw();
	c->Print("rel_wilson1_sys"+type+".pdf");

	wilson_sys2->SetBinContent(wilson_sys2->GetXaxis()->FindBin(8.5),0);
	wilson_sys2->SetBinError(wilson_sys2->GetXaxis()->FindBin(8.5),0);
	wilson_sys2->SetBinContent(wilson_sys2->GetXaxis()->FindBin(10.5),0);
	wilson_sys2->SetBinError(wilson_sys2->GetXaxis()->FindBin(10.5),0);
	wilson_sys2->SetBinContent(wilson_sys2->GetXaxis()->FindBin(13.),0);
	wilson_sys2->SetBinError(wilson_sys2->GetXaxis()->FindBin(13.),0);

	wilson_sys2->Write("rel_wilson_sys2");
	wilson_sys2->SetTitle("Relative Wilson Coeff variation sys (CNP[7]=-0.03, CNP[9]=-1.7)");
	wilson_sys2->Draw();
	c->Print("rel_wilson2_sys"+type+".pdf");

	wilson_sys3->SetBinContent(wilson_sys3->GetXaxis()->FindBin(8.5),0);
	wilson_sys3->SetBinError(wilson_sys3->GetXaxis()->FindBin(8.5),0);
	wilson_sys3->SetBinContent(wilson_sys3->GetXaxis()->FindBin(10.5),0);
	wilson_sys3->SetBinError(wilson_sys3->GetXaxis()->FindBin(10.5),0);
	wilson_sys3->SetBinContent(wilson_sys3->GetXaxis()->FindBin(13.),0);
	wilson_sys3->SetBinError(wilson_sys3->GetXaxis()->FindBin(13.),0);

	wilson_sys3->Write("rel_wilson_sys3");
	wilson_sys3->SetTitle("Relative Wilson Coeff variation sys (CNP[7]= 0.01, CNP[9]=-2.1)");
	wilson_sys3->Draw();
	c->Print("rel_wilson3_sys"+type+".pdf");

	cout << endl << endl << endl;	
	if(!percent) cout << "\n\n" << xvarname << " bin\t\t& Value \t & Stats";
	else cout << "\n\n" << xvarname << " bin ";
	for(unsigned s = 0; s < sysnames.size(); s++) cout << " \t & " << sysnames[s];
	cout << " \t & Total \\\\" << endl;


	for(unsigned i = 0; i < effnames.size(); i++)
	{
		TGraphErrors * grtot_eff = new TGraphErrors();
		cout << endl << effnames[i] << endl;

		for(int j = 1; j <= toteff->GetNbinsX(); j++)
		{
			if((xbins[j]==10.1 && !rel) || xbins[j]==15 || xbins[j] == 9.1 || xbins[j]==11) continue;
			cout << fixed << setprecision(1) << xbins[j-1] << "-" << xbins[j] << fixed << setprecision(3) << "  \t & ";
		
			double cureff = hdefault[i]->GetBinContent(j);
			if(!percent) cout << fixed << setprecision(5) << cureff << " \t & " << hdefault[i]->GetBinError(j) << "   \t & ";
		
			for(unsigned s = 0; s < sysnames.size(); s++)
				if(percent) cout << sys_eff[s][i]->GetBinContent(j)*100 << "\\% \t & ";
				else cout << sys_eff[s][i]->GetBinContent(j)*cureff << " \t & ";
			
			if(percent) cout << tot_sys_eff[i]->GetBinContent(j)*100 << "\\% ";
			else cout << tot_sys_eff[i]->GetBinContent(j)*cureff;

			//if(percent) cout << wilson_sys_eff[i]->GetBinContent(j)*100 << "\\%";
			//else cout << wilson_sys_eff[i]->GetBinContent(j)*cureff;
			cout << " \\\\ " << endl;

			grtot_eff->SetPoint(j,tot_sys_eff[i]->GetBinCenter(j),tot_sys_eff[i]->GetBinContent(j));
			grtot_eff->SetPointError(j,tot_sys_eff[i]->GetBinWidth(j)/2.,0.);
		}

		grtot_eff->Write(effnames[i]+"sys");
	}

	cout << endl << endl;	
	if(!percent) cout << "\n\n" << xvarname << " bin\t\t& Value \t\t & Stats";
	else cout << "\n\n" << xvarname << " bin ";
	for(unsigned s = 0; s < sysnames.size(); s++) cout << " \t\t & " << sysnames[s];
	cout << " \t\t & Wilson Coeff sys";
	cout << " \\\\" << endl;

	TGraphErrors * grtot = new TGraphErrors();
	for(int j = 1; j <= toteff->GetNbinsX(); j++)
	{
		if((xbins[j]==11 && !rel) || xbins[j]==15) continue;
		cout << fixed << setprecision(1) << xbins[j-1] << "-" << xbins[j] << fixed << setprecision(3);
		double cureff = toteff->GetBinContent(j);
		if(!percent) cout << fixed << setprecision(5) << cureff << " \t\t & " << toteff->GetBinError(j);
		
		for(unsigned s = 0; s < sysnames.size(); s++)
			if(percent) cout << "  \t\t & " << sys[s]->GetBinContent(j)*100 << "\\% ";
			else cout << " \t\t & " << sys[s]->GetBinContent(j)*cureff;
			
		//if(percent) cout << " \t\t & " << wilson_sys->GetBinContent(j)*100 << "\\% ";
		//else cout << " \t\t & " << wilson_sys->GetBinContent(j)*cureff;

		cout << " \\\\ " << endl;

		grtot->SetPoint(j,tot_sys->GetBinCenter(j),tot_sys->GetBinContent(j));
		grtot->SetPointError(j,tot_sys->GetBinWidth(j)/2.,0.);
	}
	
	grtot->Write("totsys");

/*
	cout << "\n\n" << xvarname << " bin\t\t& Value \\\\" << endl;
	for(int j = 1; j <= toteff->GetNbinsX(); j++)
	{
		if((xbins[j]==11 && !rel) || xbins[j]==15) continue;
		cout << fixed << setprecision(1) << xbins[j-1] << "-" << xbins[j] << fixed << setprecision(5) << " \t & ";
		cout << "$" << toteff->GetBinContent(j) << " \\pm " << tot_sys->GetBinContent(j)*toteff->GetBinContent(j) << "$ \\\\ " << endl;
	}
*/
	delete MCfile;
	delete histFile;
	return 0;
}

