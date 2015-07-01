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
	string analysis = "Lb2Lmumu";
	string type = "All";

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
		}
	}


	// Set trees of cancdidates previously created


	/**     PDF systematic      */

	gROOT->ProcessLine(".x lhcbStyle.C");

	string datafilename = "/afs/cern.ch/work/p/pluca/weighted/Lmumu/candLb.root";
	string candfilename = "/afs/cern.ch/work/p/pluca/weighted/Lmumu/candLb_MC.root";
		
	RooRealVar * vMuMu = new RooRealVar("Lb_MassConsLambda","Lb_MassConsLambda",5620.,5350.,6000.);
	Analysis * anaLbMuMu_MC = new Analysis("Lb2Lmumu_"+type+"_MC","Lb","candLb2Lmumu",candfilename,vMuMu,&CutsDef::cutMuMu);

	Analysis * anaLbMuMu_default = new Analysis("Lmumu_default","Lb","candLb2Lmumu",datafilename,vMuMu,&CutsDef::cutMuMu);
	Analysis * anaLbMuMu_Gauss = new Analysis("mumu_Gauss","Lb","candLb2Lmumu",datafilename,vMuMu,&CutsDef::cutMuMu);
	Analysis * anaLbMuMu_DGauss = new Analysis("mumu_DGauss","Lb","candLb2Lmumu",datafilename,vMuMu,&CutsDef::cutMuMu);
	Analysis * anaLbMuMu_CB = new Analysis("mumu_CB","Lb","candLb2Lmumu",datafilename,vMuMu,&CutsDef::cutMuMu);
	Analysis * anaLbMuMu_ChebKS = new Analysis("mumu_KS","Lb","candLb2Lmumu",datafilename,vMuMu,&CutsDef::cutMuMu);

	RooRealVar * vJpsi = new RooRealVar("Lb_MassConsJpsiLambda","Lb_MassConsJpsiLambda",5620.,5250.,6000.);
	Analysis * anaLbJpsi_default = new Analysis("Jpsi_default","Lb","candLb2JpsiL",datafilename,vJpsi,&CutsDef::cutJpsi);
	Analysis * anaLbJpsi_DGauss = new Analysis("Jpsi_DGauss","Lb","candLb2JpsiL",datafilename,vJpsi,&CutsDef::cutJpsi);
	Analysis * anaLbJpsi_Gauss = new Analysis("Jpsi_Gauss","Lb","candLb2JpsiL",datafilename,vJpsi,&CutsDef::cutJpsi);
	Analysis * anaLbJpsi_CB = new Analysis("Jpsi_CB","Lb","candLb2JpsiL",datafilename,vJpsi,&CutsDef::cutJpsi);
	Analysis * anaLbJpsi_ChebKS = new Analysis("Jpsi_KS","Lb","candLb2JpsiL",datafilename,vJpsi,&CutsDef::cutJpsi);
	Analysis * anaLbJpsi_noKst = new Analysis("Jpsi_noKst","Lb","candLb2JpsiL",datafilename,vJpsi,&CutsDef::cutJpsi);

	TFile * MCFile = TFile::Open(candfilename.c_str());
	TTree * BdJpsiKSTree = (TTree *)MCFile->Get("candBdJpsiKS");
	TTree * BuKstmumuTree  = (TTree *)MCFile->Get("candBuKstmumu");

	TFile * outFile = new TFile("Pdfsys.root","recreate");

	cout << "Computing sig PDF systematic" << endl;

	TCut highQ2cut = CutsDef::highQ2;

	double fitType = false;
	string model = "DCB_Sn";
	string optionsjpsi = "-ANDpulls-log-stdAxis-quiet-XM(#Lambda#mu#mu) (MeV/c^{2})-noCost-nochi2";
	string optionsmumu = "-stdAxis-quiet-XM(#Lambda#mu#mu) (MeV/c^{2})-noCost-nochi2";



	anaLbMuMu_MC->SetSignal((model+"-s[9]-s2[30]").c_str());
	anaLbMuMu_MC->Initialize("");
	anaLbMuMu_MC->Fit(5400.,5750.,200,true,"-stdAxis-XM(#Lambda#mu#mu) (MeV/c^{2})-noCost-nochi2",CutsDef::mumuTrueID);
	Str2VarMap MCpars = anaLbMuMu_MC->GetSigParams();

	Analysis * KS = new Analysis("KS_bkg","Lb",BdJpsiKSTree,vJpsi,"DCB_OST","","-namepar");
	KS->Fit(5300.,6000.,150,true,"-noPlot-quiet");
	Str2VarMap pars = KS->GetSigParams();
	RooRealVar * m_shift = new RooRealVar("m_shift","shift",0.,-5.,10.);
	setConstant(&pars);
	ModifyPars(&pars,"m",m_shift,"-shift");

	string jpsimodel = model+"-Xn"+Form("[%f]",MCpars["n"]->getVal());
	jpsimodel += (string)"-s[7,1,12]";//+Form("[%f]",MCpars["s"]->getVal());
	jpsimodel += (string)"-s2[15,8,30]";//+Form("[%f,]",MCpars["s2"]->getVal());
	jpsimodel += (string)"-a"+Form("[%f]",MCpars["a"]->getVal());
	jpsimodel += (string)"-a2"+Form("[%f]",MCpars["a2"]->getVal());
	jpsimodel += (string)"-f"+Form("[%f]",MCpars["f"]->getVal());

	anaLbJpsi_default->SetSignal(jpsimodel.c_str());
	RooRealVar * nKSjpsi = new RooRealVar("nKSjpsi","nKSjpsi",2.e3,0,1.e4);
	anaLbJpsi_default->addBkgComponent("JpsiKS","DCB_OST",nKSjpsi,"",pars);
	//anaLbJpsi_default->addBkgComponent("BuKst",BuKstmumuTree,350.);
	anaLbJpsi_default->Initialize();
	anaLbJpsi_default->Fit(5300.,6000.,200,true,optionsjpsi);
	double njpsiDCB = anaLbJpsi_default->GetSigVal();

	Str2VarMap jpsiSigpars = anaLbJpsi_default->GetSigParams();
	setConstant(&jpsiSigpars);
	RooRealVar * factor = new RooRealVar("factor","factor",1.,0.1,3.);
	ModifyPars(&jpsiSigpars,"s",factor);
	ModifyPars(&jpsiSigpars,"s2",factor);
	//Str2VarMap allpars = anaLbJpsi_default->GetParams();
	double nKSjpsi_val = nKSjpsi->getVal();

	Analysis * KSmumu = new Analysis("KS_mumu","Lb","candBdKSmumu",candfilename,vJpsi,&highQ2cut);
	double nKSmumu_MC_all = KSmumu->reducedTree->GetEntries();
	KSmumu->SetSignal("DCB_OST",1.e4,"-namepar");
	KSmumu->Initialize("-docuts");
	double nKSmumu_MC = KSmumu->reducedTree->GetEntries();

	double Jpsi2mumuBr = 0.0593;
	double BRKS_mumuOverJpsi = 3.4e-7/8.73e-4;
	double nKSmumu_val = (BRKS_mumuOverJpsi/Jpsi2mumuBr)*nKSjpsi_val*nKSmumu_MC/nKSmumu_MC_all;
	RooRealVar * nKSmumu = new RooRealVar("nKSmumu","nKSmumu",nKSmumu_val,0.,50.);
	nKSmumu->setConstant(true);

	anaLbMuMu_default->SetSignal(model.c_str(),30.,"-namepars",jpsiSigpars);
	anaLbMuMu_default->addBkgComponent("exp","Exp-b[0.0001]",20.);
	anaLbMuMu_default->addBkgComponent("KSmumu",KSmumu->reducedTree,nKSmumu);
	anaLbMuMu_default->Initialize("");
	anaLbMuMu_default->Fit(5400.,6000.,50,true,optionsmumu,highQ2cut);
	double nsigDCB = anaLbMuMu_default->GetSigVal();

	anaLbMuMu_Gauss->SetSignal("Gauss");
	anaLbMuMu_Gauss->addBkgComponent("exp","Exp-b[0.0001]",20.);
	anaLbMuMu_Gauss->addBkgComponent("KSmumu",KSmumu->reducedTree,nKSmumu);
	anaLbMuMu_Gauss->Initialize("");
	anaLbMuMu_Gauss->Fit(5400.,6000.,50,fitType,optionsmumu,highQ2cut);
	double nsigGauss = anaLbMuMu_Gauss->GetSigVal();

	anaLbMuMu_DGauss->SetSignal("DGauss-s1[7,5,15]-s2[15,10,30]");
	anaLbMuMu_DGauss->addBkgComponent("exp","Exp-b[0.0001]",20.);
	anaLbMuMu_DGauss->addBkgComponent("KSmumu",KSmumu->reducedTree,nKSmumu);
	anaLbMuMu_DGauss->Initialize("");
	anaLbMuMu_DGauss->Fit(5400.,6000.,50,fitType,optionsmumu,highQ2cut);
	double nsigDGauss = anaLbMuMu_DGauss->GetSigVal();

	anaLbMuMu_CB->SetSignal("CB");
	anaLbMuMu_CB->addBkgComponent("exp","Exp-b[0.0001]",20.);
	anaLbMuMu_CB->addBkgComponent("KSmumu",KSmumu->reducedTree,nKSmumu);
	anaLbMuMu_CB->Initialize("");
	anaLbMuMu_CB->Fit(5400.,6000.,50,fitType,optionsmumu,highQ2cut);
	double nsigCB = anaLbMuMu_CB->GetSigVal();

	// ----------

	anaLbJpsi_Gauss->SetSignal("Gauss",1.e4,"-namepar");
	anaLbJpsi_Gauss->addBkgComponent("JpsiKS","DCB_OST",3.e3,"",pars);
	//anaLbJpsi_Gauss->addBkgComponent("BuKst",BuKstmumuTree,350.);
	anaLbJpsi_Gauss->Initialize();
	anaLbJpsi_Gauss->Fit(5300.,6000.,200,fitType,optionsjpsi);
	double njpsiGauss = anaLbJpsi_Gauss->GetSigVal();

	anaLbJpsi_DGauss->SetSignal("DGauss-s1[7,5,15]-s2[15,10,30]",1.e4,"-namepar");
	anaLbJpsi_DGauss->addBkgComponent("JpsiKS","DCB_OST",3.e3,"",pars);
	//anaLbJpsi_DGauss->addBkgComponent("BuKst",BuKstmumuTree,350.);
	anaLbJpsi_DGauss->Initialize();
	anaLbJpsi_DGauss->Fit(5300.,6000.,200,fitType,optionsjpsi);
	double njpsiDGauss = anaLbJpsi_DGauss->GetSigVal();

	anaLbJpsi_CB->SetSignal("CB",1.e4,"-namepar");
	//anaLbJpsi_CB->addBkgComponent("JpsiKS","DCB_OST",3.e3,"",pars);
	anaLbJpsi_CB->addBkgComponent("BuKst",BuKstmumuTree,350.);
	anaLbJpsi_CB->Initialize();
	anaLbJpsi_CB->Fit(5300.,6000.,200,fitType,optionsjpsi);
	double njpsiCB = anaLbJpsi_CB->GetSigVal();

	// -----------

	anaLbJpsi_ChebKS->SetSignal((model+"-s1[7,5,15]-s2[15,10,30]").c_str());
	anaLbJpsi_ChebKS->addBkgComponent("JpsiKS",BdJpsiKSTree);
	anaLbJpsi_ChebKS->addBkgComponent("BuKst",BuKstmumuTree,350.);
	anaLbJpsi_ChebKS->Initialize();
	anaLbJpsi_ChebKS->Fit(5300.,6000.,200,fitType,optionsjpsi);
	double njpsi_ChebKS = anaLbJpsi_ChebKS->GetSigVal();

	anaLbJpsi_noKst->SetSignal((model+"-s1[7,5,15]-s2[15,10,30]").c_str());
	anaLbJpsi_noKst->addBkgComponent("JpsiKS","DCB_OST",nKSjpsi,"",pars);
	anaLbJpsi_noKst->Initialize();
	anaLbJpsi_noKst->Fit(5300.,6000.,200,fitType,optionsjpsi);
	double njpsi_noKst = anaLbJpsi_ChebKS->GetSigVal();

	anaLbMuMu_ChebKS->SetSignal((model+"-s1[7,5,15]-s2[15,10,30]").c_str(),30.,"",jpsiSigpars);
	anaLbMuMu_ChebKS->addBkgComponent("exp","Exp-b[0.0001]",20.);
	anaLbMuMu_ChebKS->Initialize("");
	anaLbMuMu_ChebKS->Fit(5400.,6000.,50,fitType,optionsmumu,highQ2cut);
	double nsigDCB_ChebKS = anaLbMuMu_ChebKS->GetSigVal();

	cout << "-----------------------------------------------------------------" << endl;	
	anaLbMuMu_default->PrintChi2();
	cout << "Nevts: " << nsigDCB << endl;
	anaLbMuMu_Gauss->PrintChi2();
	cout << "Nevts: " << nsigGauss << endl;
	anaLbMuMu_DGauss->PrintChi2();
	cout << "Nevts: " << nsigDGauss << endl;
	anaLbMuMu_CB->PrintChi2();
	cout << "Nevts: " << nsigCB << endl;
	anaLbMuMu_ChebKS->PrintChi2();
	cout << "Nevts: " << nsigDCB_ChebKS << endl;
	cout << "-----------------------------------------------------------------" << endl;		
	anaLbJpsi_default->PrintChi2();
	cout << "Nevts: " << njpsiDCB << endl;
	anaLbJpsi_DGauss->PrintChi2();
	cout << "Nevts: " << njpsiDGauss << endl;
	anaLbJpsi_Gauss->PrintChi2();
	cout << "Nevts: " << njpsiGauss << endl;
	anaLbJpsi_CB->PrintChi2();
	cout << "Nevts: " << njpsiCB << endl;
	anaLbJpsi_ChebKS->PrintChi2();
	cout << "Nevts: " << njpsi_ChebKS << endl;
	anaLbJpsi_noKst->PrintChi2();
	cout << "Nevts: " << njpsi_noKst << endl;
	cout << "-----------------------------------------------------------------" << endl;

	Analysis::SetPrintLevel("s");
	gStyle->SetOptStat("ReM");


	int nbins = 5;
	double q2min[] = {15.,11.0,15,16,18};
	double q2max[] = {20.,12.5,16,18,20};	
	
	int nexp = 50;

	vector<TH1 *> pdfsysval, bkgsysval;
	TGraph * gr = new TGraph(); 


	for(int q = 0; q < nbins; q++)
	{
		outFile->cd();
		TH1F * hsys = new TH1F(Form("hsys_%i",q),"hsys",40,-1,1);
		TH1F * hbkg = new TH1F(Form("hbkg_%i",q),"hbkg",40,-1,1);

		TCut curq2cut = (TCut)Form("TMath::Power(J_psi_1S_MM/1000,2) >= %e  && TMath::Power(J_psi_1S_MM/1000,2) < %e",q2min[q],q2max[q]);

		Analysis * myKSmumu = new Analysis("KS_mumu","Lb","candBdKSmumu",candfilename,vJpsi,&curq2cut);
		double nKSmumu_MC_all = myKSmumu->reducedTree->GetEntries();
		myKSmumu->SetSignal("DCB_OST",1.e4,"-namepar");
		myKSmumu->Initialize("-docuts");
		double nKSmumu_MC = myKSmumu->reducedTree->GetEntries();
		double nKSmumu_val = (BRKS_mumuOverJpsi/Jpsi2mumuBr)*nKSjpsi_val*nKSmumu_MC/nKSmumu_MC_all;
		RooRealVar * mynKSmumu = new RooRealVar("nKSmumu","nKSmumu",nKSmumu_val);
		mynKSmumu->setConstant(true);

		anaLbMuMu_default->Reset();
		anaLbMuMu_default->SetSignal(model.c_str(),30.,"-namepars",jpsiSigpars);
		anaLbMuMu_default->addBkgComponent("exp","Exp-b[0.0001]",20.);
		anaLbMuMu_default->addBkgComponent("KSmumu",myKSmumu->reducedTree,mynKSmumu);
		anaLbMuMu_default->Initialize("");
		anaLbMuMu_default->Fit(5400.,6000.,50,true,"-noPlot-quiet",curq2cut);
		double nevts = anaLbMuMu_default->GetSigVal();

		cout << "Bin " << q2min[q] << "-" << q2max[q] << ", nevts = " << nevts << endl;

		for(int e = 0; e < nexp; e++)
		{
			showPercentage(e,nexp);
	
			RooAbsPdf * mumupdf = (RooAbsPdf*)anaLbMuMu_default->GetModel()->Clone();
			RooAbsPdf * jpsipdf = (RooAbsPdf*)anaLbJpsi_default->GetModel()->Clone();
			Analysis * mumu = new Analysis("mumu",vMuMu,mumupdf,nevts);
			Analysis * jpsi = new Analysis("jpsi",vJpsi,jpsipdf,9965);
	
			jpsi->Reset();
			mumu->Reset();
	
			jpsi->SetSignal("DGauss-s1[7,5,15]-s2[15,10,30]");
			jpsi->addBkgComponent("JpsiKS","DCB_OST",2.e3,"",pars);
			jpsi->Initialize();
			jpsi->Fit(5350.,6000.,1000,true,"-noPlot-quiet-noCost");
			Str2VarMap jpsiPars = jpsi->GetSigParams();

			setConstant(&jpsiPars);
			RooRealVar * factor1 = new RooRealVar("factor1","factor1",1.,0.1,3.);
			ModifyPars(&jpsiPars,"s",factor1);
			ModifyPars(&jpsiPars,"s2",factor1);
			
			mumu->SetSignal("DGauss",30,"-namepars",jpsiPars);
			mumu->addBkgComponent("exp","Exp-b[-0.001,-0.1,0.1]",20.);
			mumu->addBkgComponent("KSmumu",myKSmumu->reducedTree,mynKSmumu);
			mumu->Initialize("");
			mumu->Fit(5400.,6000.,30,true,"-noPlot-quiet-noCost");
			double rel_dgauss = mumu->GetSigVal() / jpsi->GetSigVal();
	
			jpsi->Reset();
			mumu->Reset();
		
			jpsi->SetSignal(jpsimodel.c_str());
			jpsi->addBkgComponent("JpsiKS","DCB_OST",2.e3,"",pars);
			jpsi->Initialize();
			jpsi->Fit(5400.,6000.,1000,true,"-noPlot-quiet-noCost");
			Str2VarMap jpsiPars2 = jpsi->GetSigParams();

			setConstant(&jpsiPars2);
			RooRealVar * factor2 = new RooRealVar("factor2","factor2",1.,0.1,3.);
			ModifyPars(&jpsiPars2,"s",factor2);
			ModifyPars(&jpsiPars2,"s2",factor2);
			
			mumu->SetSignal(model.c_str(),30.,"-namepars",jpsiPars2);
			mumu->addBkgComponent("exp","Exp-b[-0.001,-0.1,0.1]",20.);
			mumu->addBkgComponent("KSmumu",KSmumu->reducedTree,nKSmumu_val);
			mumu->Initialize("");
			mumu->Fit(5400.,6000.,30,true,"-noPlot-noCost-quiet");
			double rel_bkg = mumu->GetSigVal() / jpsi->GetSigVal();

			jpsi->Reset();
			mumu->Reset();

			jpsi->SetSignal(jpsimodel.c_str());
			jpsi->addBkgComponent("JpsiKS","DCB_OST",2.e3,"",pars);
			jpsi->Initialize();
			jpsi->Fit(5400.,6000.,1000,true,"-noCost-noPlot-quiet");
			Str2VarMap jpsiPars3 = jpsi->GetSigParams();

			setConstant(&jpsiPars3);
			RooRealVar * factor3 = new RooRealVar("factor3","factor3",1.,0.1,3.);
			ModifyPars(&jpsiPars3,"s",factor3);
			ModifyPars(&jpsiPars3,"s2",factor3);
			
			mumu->SetSignal(model.c_str(),30.,"-namepars",jpsiPars3);
			mumu->addBkgComponent("exp","Exp-b[-0.001,-0.1,0.1]",20.);
			mumu->addBkgComponent("KSmumu",KSmumu->reducedTree,mynKSmumu);
			mumu->Initialize("");
			mumu->Fit(5400.,6000.,30,true,"-noPlot-noCost-quiet");
			double rel_def = mumu->GetSigVal() / jpsi->GetSigVal();

			if(TMath::Abs((rel_def - rel_dgauss)/rel_def) < 0.3) hsys->Fill((rel_def - rel_dgauss)/rel_def);
			if(TMath::Abs((rel_def - rel_bkg)/rel_def) < 0.3) hbkg->Fill((rel_def - rel_bkg)/rel_def);

		}

		pdfsysval.push_back(hsys);
		bkgsysval.push_back(hbkg);

		gr->SetPoint(q,(q2max[q]-q2min[q])/2.,TMath::Sqrt(TMath::Power(bkgsysval[q]->GetMean(),2) + TMath::Power(pdfsysval[q]->GetMean(),2)));
	}

	
	for(int q = 0; q < nbins; q++)
	{
		cout << "Bin " << q2min[q] << "-" << q2max[q] << endl;
		cout << "Bkg PDF sys = " << bkgsysval[q]->GetMean()*100 << "% ( error " << bkgsysval[q]->GetMeanError()*100 << ")" << endl;
		cout << "Sig PDF sys = " << pdfsysval[q]->GetMean()*100 << "% ( error " << pdfsysval[q]->GetMeanError()*100 << ")" << endl;
		cout << "Tot PDF sys = " << TMath::Sqrt(TMath::Power(bkgsysval[q]->GetMean(),2) + TMath::Power(pdfsysval[q]->GetMean(),2))*100 << "%" << endl;
	}


	MCFile->Close();
	//outFile->Write();
	gr->Write("pdfsysgr");
	outFile->Close();
	delete outFile;
	delete MCFile;

	return 0;
}

