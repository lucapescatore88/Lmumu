#include <vector>
#include <sstream>
#include <iostream>
#include <string>
#include <time.h>

#include "TROOT.h"
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

#include "TGraphErrors.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TEntryList.h"
#include "TCanvas.h"
#include "TBranch.h"
#include "TCut.h"

#include "general_functions.hpp"
#include "ReadTree_comp.hpp"
#include "analyser.hpp"
#include "Lb_cuts.hpp"

using namespace std;
using namespace RooFit;


int main(int argc, char **argv)
{
	string analysis = "Lb2Lmumu";
	string type = "All";
	string mother = "";
	string logPlot = "log";
	bool useSingle = false;

	if(argc > 1)
    {
          for(int a = 1; a < argc; a++)
          {
                string arg = argv[a];
                string str = arg.substr(2,arg.length()-2);

				if(arg.find("-m") != string::npos) mother = str;
                if(arg.find("-t") != string::npos) type = str;
				if(arg == "-lin") logPlot = "lin";
				if(arg == "-s") useSingle = true;
          }
    }

	gROOT->ProcessLine(".x lhcbStyle.C");

	TCut baseCut = "", singleCut = "";
	string tname = "cand";
	if(useSingle) { singleCut = "isChosenCand"; tname = "singleCand_"; cout << "Using single candidates" << endl; }

	string candfilename = "/afs/cern.ch/work/p/pluca/weighted/Lmumu/candLb_MC.root";
	string datafilename = "/afs/cern.ch/work/p/pluca/weighted/Lmumu/candLb.root";
	
	RooRealVar * vLbCons = new RooRealVar("Lb_MassConsLambda","Lb_MassConsLambda",5621.,5400.,6000.);
	RooRealVar * vLbJpsiCons = new RooRealVar("Lb_MassConsJpsiLambda","Lb_MassConsJpsiLambda",5621.5,5350.,6100.);

	Analysis * anaLbMuMu = new Analysis("Lb2Lmumu_"+type+"_MC","Lb",tname+"Lb2Lmumu",candfilename,vLbCons);
	Analysis * anaLbJpsi = new Analysis("Lb2JpsiL_"+type+"_MC","Lb",tname+"Lb2JpsiL",candfilename,vLbJpsiCons);
	Analysis * anaLbMuMu_data = new Analysis("Lb2Lmumu","Lb",tname+"Lb2Lmumu",datafilename,vLbCons);
	Analysis * anaLbJpsi_data = new Analysis("Lb2JpsiL_"+type+"_data","Lb",tname+"Lb2JpsiL",datafilename,vLbJpsiCons);
	Analysis * anaLbJpsi_data_lowSel = new Analysis("Lb2JpsiL_"+type+"_data_lowSel","Lb",tname+"Lb2JpsiL",datafilename,vLbJpsiCons,&CutsDef::cutJpsi_lowSel);


	if(mother=="Lb") baseCut += "pplus_ID>0";
	else if(mother=="Lbbar") baseCut += "pplus_ID<0";
	if(type == "DD") { baseCut += CutsDef::DDcut; }
	else if(type == "LL") { baseCut += CutsDef::LLcut; }
	

	TFile * MCFile = TFile::Open((TString)candfilename);
	TTree * BdKSmumuTree  = (TTree *)MCFile->Get("candBdKSmumu");
	TTree * BdJpsiKSTree  = (TTree *)MCFile->Get("candBdJpsiKS");
	//TTree * BuKstmumuTree  = (TTree *)MCFile->Get("candBuKstmumu");
	TTree * JpsiTailTree  = (TTree *)MCFile->Get("candJpsiTail");
	MCFile = TFile::Open("/afs/cern.ch/work/p/pluca/weighted/Lmumu/Bd2KSmumu_MC12_NBweighted.root");
	TTree * BdKSmumuTree_All  = (TTree *)MCFile->Get("EventTuple");
	MCFile = TFile::Open("/afs/cern.ch/work/p/pluca/weighted/Lmumu/Bd2JpsiKS_MC12_NBweighted.root");
	TTree * BdJpsiKSTree_All  = (TTree *)MCFile->Get("EventTuple");
	double KSmumuSelected = BdKSmumuTree->Draw("J_psi_1S_MM",baseCut);
	double JpsiKSSelected = BdJpsiKSTree->Draw("J_psi_1S_MM",baseCut);
	double effKSmumu = KSmumuSelected / (double) BdKSmumuTree_All->GetEntries();
	double effJpsiKS = JpsiKSSelected / (double) BdJpsiKSTree_All->GetEntries();
	double KS_double_eff_ratio = effKSmumu / effJpsiKS;

	TFile * histFile = new TFile(("Lbyield_"+type+".root").c_str(),"recreate");
	
	// FIT Lb data

	RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

	string model = "DCB_Sn";
	string optionsjpsi = "-quiet-range-ANDpulls-lin-stdAxis-XM(#Lambda#mu#mu) (MeV/c^{2})-noCost-nochi2-min2";

	anaLbJpsi->SetSignal((model+"-s[7]-s2[15]").c_str());
	anaLbJpsi->Initialize("");
	anaLbJpsi->Fit(5400.,5750.,200,true,"-quiet-stdAxis-XM(#Lambda#mu#mu) (MeV/c^{2})-noCost-nochi2",baseCut);
	Str2VarMap MCpars = anaLbJpsi->GetSigParams();
	PrintPars(MCpars);



	Analysis * KS = new Analysis("KS_bkg","Lb",BdJpsiKSTree,vLbJpsiCons,"DCB_OST","","-namepar");
	KS->Fit(5350.,6000.,100,true,"-quiet-noParams-stdAxis-XM(#Lambda#mu#mu) (MeV/c^{2})");
	Str2VarMap pars = KS->GetSigParams();
	RooRealVar * m_shift = new RooRealVar("m_shift","shift",2.,-10.,10.);
	if(type=="LL") m_shift = new RooRealVar("m_shift","shift",0.);
	setConstant(&pars);
	ModifyPars(&pars,"m",m_shift,"-shift");

	//Fit Jpsi
	string jpsimodel = model+"-Xn"+Form("[%f]",MCpars["n"]->getVal());
	jpsimodel += (string)"-s[7,1,12]";//+Form("[%f]",MCpars["s"]->getVal());
	jpsimodel += (string)"-s2[15,8,30]";//+Form("[%f,]",MCpars["s2"]->getVal());
	jpsimodel += (string)"-a"+Form("[%f]",MCpars["a"]->getVal());
	jpsimodel += (string)"-a2"+Form("[%f]",MCpars["a2"]->getVal());
	jpsimodel += (string)"-f"+Form("[%f]",MCpars["f"]->getVal());


	histFile->cd();
	anaLbJpsi_data->SetSignal(jpsimodel.c_str());
	anaLbJpsi_data->addBkgComponent("Comb","Exp-b[-0.003]",2.e3);
	RooRealVar * nKSjpsi = new RooRealVar("nKSjpsi","N_{K_SJ\\psi}",2.e3,0,1.e4);
	anaLbJpsi_data->addBkgComponent("JpsiKS","DCB_OST",nKSjpsi,"",pars);
	//anaLbJpsi_data->addBkgComponent("BuKst",BuKstmumuTree,1.,"-s1");
	anaLbJpsi_data->Initialize("");
	anaLbJpsi_data->Fit(5350.,6000.,150,true,optionsjpsi,baseCut+singleCut);
	anaLbJpsi_data->Print(true,(string)(((TString)optionsjpsi).ReplaceAll("-lin","-log"))+"-noParams",150);
	double valerr_jpsi;
	double val_jpsi = anaLbJpsi_data->GetSigVal(&valerr_jpsi);
	TGraphErrors * jpsihighq2gr = new TGraphErrors();
	jpsihighq2gr->SetPoint(1,0,val_jpsi);
	jpsihighq2gr->SetPointError(1,0,valerr_jpsi);
	jpsihighq2gr->Write("jpsiyield_highq2");

	Str2VarMap allpars = anaLbJpsi_data->GetParams();
	double nKSjpsi_val = nKSjpsi->getVal();
	double nKSjpsi_err = nKSjpsi->getError();

	double BRJpsiKS = (nKSjpsi_val/(effJpsiKS*0.16*2))/(val_jpsi/(0.00645));
	if(type=="LL") BRJpsiKS = (nKSjpsi_val/(effJpsiKS*0.16*2))/(val_jpsi/(0.00222));
	cout << "BR(JpsiKS)/BR(Lmumu) = " << BRJpsiKS << " +/- " << BRJpsiKS*TMath::Sqrt(TMath::Power(nKSjpsi_err/nKSjpsi_val,2) + TMath::Power(valerr_jpsi/val_jpsi,2)) << endl;
	
	anaLbJpsi_data_lowSel->SetSignal(jpsimodel.c_str());
	anaLbJpsi_data_lowSel->addBkgComponent("Comb","Exp-b[-0.003]",2.e3);
	RooRealVar * nKSjpsi2 = new RooRealVar("nKSjpsi_Lo","N_{K_SJ\\psi}",2.e3,0,1.e4);
	anaLbJpsi_data_lowSel->addBkgComponent("JpsiKS_Lo","DCB_OST",nKSjpsi2,"",pars);
	//anaLbJpsi_data_lowSel->addBkgComponent("BuKst",BuKstmumuTree,1.1,"-s1");	
	anaLbJpsi_data_lowSel->Initialize("-docuts");
	anaLbJpsi_data_lowSel->Fit(5350.,6000.,150,true,optionsjpsi,baseCut+singleCut);
	double valerr_jpsi_low;
	double val_jpsi_low = anaLbJpsi_data_lowSel->GetSigVal(&valerr_jpsi_low);
	TGraphErrors * jpsilowq2gr = new TGraphErrors();
	jpsilowq2gr->SetPoint(1,0,val_jpsi_low);
	jpsilowq2gr->SetPointError(1,0,valerr_jpsi_low);
	jpsilowq2gr->Write("jpsiyield_lowq2");



	
	TGraphErrors* q2Plot = new TGraphErrors();
	int nbins = CutsDef::nq2bins;
	double * q2min = &CutsDef::q2min_highfirst[0];
	double * q2max = &CutsDef::q2max_highfirst[0];
	//int nbins = 1;
	//double q2min[] = {6};
	//double q2max[] = {8};

	TCanvas *c = new TCanvas("c","");
	TCanvas *c2 = new TCanvas("c2","");
	TCanvas *ctemp = new TCanvas("ctemp","");
	c->Divide(2,2);
	c2->Divide(2,2);
	
	string options = "-quiet-lin-stdAxis-range-minos-XM(#Lambda#mu#mu) (MeV/c^{2})-noCost-noParams-quiet";
	vector<double> nKSmumu_v;
	m_shift->setConstant(true);
	vector< vector<double > > vals;
	vector<Str2VarMap> rarepars;
	RooRealVar * nJpsi_tail = new RooRealVar("nJpsi_tail","N_{J/#psi tail}",1.,0.,200.);

	for(int i = 0; i < nbins; i++)
	{
		TCut curq2cut = (TCut)Form("TMath::Power(J_psi_1S_MM/1000,2) >= %e  && TMath::Power(J_psi_1S_MM/1000,2) < %e",q2min[i],q2max[i]);
		TCut curcut = curq2cut+baseCut;
	
		anaLbMuMu->Reset();
		anaLbMuMu->SetName(Form("Lb2Lmumu_MC_%i",i));
		anaLbMuMu->SetSignal((model+"-s[15]-s2[30]").c_str());
		anaLbMuMu->Initialize("");
		anaLbMuMu->Fit(5400.,5750.,200,true,options,curcut);
		Str2VarMap MCmumupars = anaLbMuMu->GetSigParams();
		
		Str2VarMap rarePars = anaLbJpsi_data->GetSigParams();
		rarePars["n"] = MCmumupars["n"];
		rarePars["a"] = MCmumupars["a"];
		rarePars["a2"] = MCmumupars["a2"];
		setConstant(&rarePars);
		RooRealVar * factor = new RooRealVar(Form("factor_%i",i),"c",MCmumupars["s"]->getVal()/MCpars["s"]->getVal());
		ModifyPars(&rarePars,"s",factor);
		ModifyPars(&rarePars,"s2",factor);
			
		anaLbMuMu_data->Reset();
		if((q2min[i]==1.1 && q2max[i]==6) || (q2min[i]==15 && q2max[i]==20))
		{
			if(q2min[i]==1.1) anaLbMuMu_data->SetName(("Lb2Lmumu_"+type+"_lowQ2").c_str());
			else anaLbMuMu_data->SetName(("Lb2Lmumu_"+type+"_highQ2").c_str());
		}
		else anaLbMuMu_data->SetName(Form("q2_%i",i));
	
		TString KSname = Form("KS_bkg_"+(TString)type+"_%i",i);
		Analysis * KSmumu = new Analysis(KSname,"Lb","candBdKSmumu",candfilename,vLbCons,&curq2cut);
		double nKSmumu_MC_all = KSmumu->reducedTree->GetEntries();		
		KSmumu->SetSignal("DCB_OST",1.e4,"-namepar");
		KSmumu->Initialize("-docuts");
		double nKSmumu_MC = KSmumu->reducedTree->GetEntries();
		
		double Jpsi2mumuBr = 0.0593;
		double BRKS_mumuOverJpsi = 3.4e-7/8.73e-4;
		double nKSmumu_val = KS_double_eff_ratio*(BRKS_mumuOverJpsi/Jpsi2mumuBr)*nKSjpsi_val*nKSmumu_MC/nKSmumu_MC_all;
		RooRealVar * nKSmumu = new RooRealVar("nKSmumu","N_{K_S\\mu\\mu}",nKSmumu_val,0.,50.);
		nKSmumu->setConstant(true);
		nKSmumu_v.push_back(nKSmumu_val);

		anaLbMuMu_data->SetSignal(model.c_str(),30.,"",rarePars);
		anaLbMuMu_data->addBkgComponent("exp","Exp-b[-0.003,-0.01,0.01]",20.);
		anaLbMuMu_data->addBkgComponent("KSmumu",KSmumu->reducedTree,nKSmumu);
		//if(q2max[i] == 8.) anaLbMuMu_data->addBkgComponent("Jpsi_tail",JpsiTailTree,nJpsi_tail);
		if(q2max[i] <= 11) anaLbMuMu_data->SetBlindRegion(5500,5700);
		anaLbMuMu_data->Initialize("");
		
		if((q2min[i]==1.1 && q2max[i]==6) || (q2min[i]==15 && q2max[i]==20))
			anaLbMuMu_data->Fit(5400.,6000.,25,true,options,curcut+singleCut);
		else
		{
			RooPlot * pl = anaLbMuMu_data->Fit(5400.,6000.,25,true,options+"-noPlot-t"+(string)Form("[%3.1f,%3.1f]",q2min[i],q2max[i]),curcut+singleCut);
			if(i <= 4) c->cd(i);
			else c2->cd(i-4);
			pl->SetTitle((TString)Form("[%.1f,%.1f]",q2min[i],q2max[i]));
			pl->Draw();
			ctemp->cd();
		}
		
		rarepars.push_back(anaLbMuMu_data->GetParams());
	
		double valerr = 0;
		//double errLow = 0, errHigh = 0;
		double val = anaLbMuMu_data->GetSigVal(&valerr);
		//double val = anaLbMuMu_data->GetSigVal(&errLow,&errHigh);
		double center = (q2max[i] + q2min[i])/2.;
		double width = q2max[i] - q2min[i];
		q2Plot->SetPoint(i+1,center,val/width);
		q2Plot->SetPointError(i+1,width/2,valerr/width);
		vector<double> vv(1.,val);
		vv.push_back(valerr);
		//vv.push_back(errHigh);
		//vv.push_back(errLow);
		vals.push_back(vv);
		factor->setConstant(false);
	}
		
	c->Print(("q2_fits_"+type+"_plot1.pdf").c_str());
	c2->Print(("q2_fits_"+type+"_plot2.pdf").c_str());

	TCanvas * cq2 = new TCanvas();
	q2Plot->SetTitle(0);
	q2Plot->Draw("AP");
	q2Plot->GetXaxis()->SetTitle("q^{2}");
	cq2->Print(("q2plot_"+type+".pdf").c_str());
	q2Plot->Write("q2plot");
	
	delete cq2;
	delete c;

	cout << "KS evts prediction:" << endl;
	for(unsigned k = 0; k < nKSmumu_v.size(); k++) cout << fixed << setprecision(1) << q2min[k] << "-" << q2max[k] << " & " << nKSmumu_v[k] << endl;

	cout << "Binning min: { ";
	for(int k = 0; k < nbins; k++) cout << fixed << setprecision(1) << q2min[k] << ", ";
	cout << "}" << endl;
	cout << "Binning max: { ";
	for(int k = 0; k < nbins; k++) cout << fixed << setprecision(1) << q2max[k] << ", ";
	cout << "}" << endl;
	cout << "N_JpsiL: $" << val_jpsi << " \\pm " << valerr_jpsi << "$" << endl;
	cout << "N_JpsiL_low: $" << val_jpsi_low << " \\pm " << valerr_jpsi_low  << "$" << endl;

	cout << "N_Lmumu: {";
	for(unsigned k = 0; k < vals.size(); k++) cout << vals[k][0] << ", ";
	cout << "}" << endl;
	cout << "N_Lmumu_err: {";
	for(unsigned k = 0; k < vals.size(); k++) cout << vals[k][1] << ", ";
	cout << "}" << endl;
	if(vals[0].size()>2) for(unsigned k = 0; k < vals.size(); k++) cout << vals[k][2] << ", ";
	cout << "}" << endl;

	for(unsigned k = 0; k < vals.size(); k++) cout << q2min[k] << "-" << q2max[k] << " & $" << vals[k][0] << " \\pm " << vals[k][1] << "$ " << endl;

	cout << "Parameters: " << endl << endl;
	PrintPars(allpars,"-latext");
	cout << "----------------" << endl;
	for(unsigned k = 0; k < rarepars.size(); k++)
	{
		cout << "\\multicolumn{2}{ " << q2min[k] << "-" << q2max[k] << " } \\\\" << endl;
		PrintPars(rarepars[k],"-latex-nocost");
	}

	anaLbMuMu->PrintChi2();
	anaLbJpsi_data->PrintChi2();

	histFile->Write();
	histFile->Close();
	delete MCFile;
	delete histFile;
	
	return 0;
}





