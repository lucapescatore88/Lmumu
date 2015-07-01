#include "ReadTree_comp.hpp"
#include "multi_analyser.hpp"
#include "Lb_cuts.hpp"
#include "RooAbsReal.h"
#include "RooSimultaneous.h"
#include "RooGenericPdf.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TROOT.h"
#include "FeldmanCousins.hpp"
#include "RooKeysPdf.h"

using namespace RooFit;
using namespace std;

TString decayToDo = "Lb2Lmumu";//"Lb2JpsiL_reduced";


bool isInAllowedArea( Str2VarMap params )
{
	double a = params["afb"]->getVal();
	double f = params["fL"]->getVal();

	//cout << "isValid " << a << "  " << f << "   ->   " << (int)!((f-1)*3./4. > a || a >  -(f-1)*3./4.) << endl;
	return !((f-1)*3./4. > a || a >  -(f-1)*3./4.);
}


bool isInAllowedAreaB( Str2VarMap params )
{
	double a = params["afbB"]->getVal();
	return (TMath::Abs(a) <= 0.5);
}




void getEfficiencies(double q2min, double q2max, RooAbsPdf ** effLLpdf, RooAbsPdf ** effDDpdf, RooAbsPdf ** effLLBpdf, RooAbsPdf ** effDDBpdf, bool printeff = false)
{
	RooRealVar * cosThetaL = new RooRealVar("cosThetaL","cosThetaL",0.,-1.,1.);
	RooRealVar * cosThetaB = new RooRealVar("cosThetaB","cosThetaB",0.,-1.,1.);

	TString q2name = ((TString)Form("_q2_%4.2f_%4.2f",q2min,q2max)).ReplaceAll(".","");
	TString effbase = "/afs/cern.ch/work/p/pluca/results/";
	TH1F * effDD, * effLL, * effLLB, * effDDB;

	if(q2min == 8 && q2max == 11)
	{
		q2name = "_jpsi";
		TFile * effFile = TFile::Open(effbase+"Lbeff_JpsivscosThetaL_DD.root");
		effDD  = (TH1F *)effFile->Get("htoteff");
		effFile = TFile::Open(effbase+"Lbeff_JpsivscosThetaL_LL.root");
		effLL  = (TH1F *)effFile->Get("htoteff");
		effFile = TFile::Open(effbase+"Lbeff_JpsivscosThetaB_DD.root");
		effDDB  = (TH1F *)effFile->Get("htot_nodet_eff");
		effFile = TFile::Open(effbase+"Lbeff_JpsivscosThetaB_LL.root");
		effLLB  = (TH1F *)effFile->Get("htot_nodet_eff");
	}
	else if(q2min == 15 && q2max == 20)
	{
		TFile * effFile = TFile::Open(effbase+"LbeffvscosThetaL_DD.root");
		effDD  = (TH1F *)effFile->Get("htoteff");
		effFile = TFile::Open(effbase+"LbeffvscosThetaL_LL.root");
		effLL  = (TH1F *)effFile->Get("htoteff");
		effFile = TFile::Open(effbase+"LbeffvscosThetaB_DD.root");
		effDDB  = (TH1F *)effFile->Get("htot_nodet_eff");
		effFile = TFile::Open(effbase+"LbeffvscosThetaB_LL.root");
		effLLB  = (TH1F *)effFile->Get("htot_nodet_eff");
	}
	else
	{
		TFile * effFile = TFile::Open(effbase+"Lbeff2D_cosThetaL_vs_q2_DD.root");
		TH2F * effDD2D  = (TH2F *)effFile->Get("htot_eff");
		effDD = (TH1F*)GetSliceX(effDD2D,(q2max+q2min)/2.);
		effFile = TFile::Open(effbase+"Lbeff2D_cosThetaL_vs_q2_LL.root");
		TH2F * effLL2D  = (TH2F *)effFile->Get("htot_eff");
		effLL = (TH1F*)GetSliceX(effLL2D,(q2max+q2min)/2.);
		effFile = TFile::Open(effbase+"Lbeff2D_cosThetaB_vs_q2_DD.root");
		TH2F * effDDB2D  = (TH2F *)effFile->Get("hnodet_eff");
		effDDB = (TH1F*)GetSliceX(effDDB2D,(q2max+q2min)/2.);
		effFile = TFile::Open(effbase+"Lbeff2D_cosThetaB_vs_q2_LL.root");
		TH2F * effLLB2D  = (TH2F *)effFile->Get("hnodet_eff");
		effLLB = (TH1F*)GetSliceX(effLLB2D,(q2max+q2min)/2.);
	}

	RooDataHist * hLL = new RooDataHist("hLL","hLL",*cosThetaL,effLL);
	RooDataHist * hDD = new RooDataHist("hDD","hDD",*cosThetaL,effDD);
	(*effLLpdf) = stringToPdf("Poly2","effLL", cosThetaL);
	(*effDDpdf) = stringToPdf("Poly2","effDD", cosThetaL);
	(*effLLpdf)->fitTo(*hLL,PrintLevel(-1));
	(*effDDpdf)->fitTo(*hDD,PrintLevel(-1));
	fixParams((*effLLpdf),cosThetaL);
	fixParams((*effDDpdf),cosThetaL);	

	RooDataHist * hLLB = new RooDataHist("hLLB","hLLB",*cosThetaB,effLLB);
	RooDataHist * hDDB = new RooDataHist("hDDB","hDDB",*cosThetaB,effDDB);
	(*effLLBpdf) = stringToPdf("Poly3","effLLB", cosThetaB);
	(*effDDBpdf) = stringToPdf("Poly3","effDDB", cosThetaB);
	(*effLLBpdf)->fitTo(*hLLB,PrintLevel(-1));
	(*effDDBpdf)->fitTo(*hDDB,PrintLevel(-1));
	fixParams((*effLLBpdf),cosThetaB);
	fixParams((*effDDBpdf),cosThetaB);

	if(printeff)
	{
		TCanvas * ceff = new TCanvas();
		GetFrame(cosThetaL, hLL,(*effLLpdf),"-nochi2",0,NULL,0,"cos#theta_{l}","Tot. eff.")->Draw();
		ceff->Print("LLeffFit"+q2name+".pdf");
		GetFrame(cosThetaL, hDD,(*effDDpdf),"-nochi2",0,NULL,0,"cos#theta_{l}","Tot. eff.")->Draw();
		ceff->Print("DDeffFit"+q2name+".pdf");
		GetFrame(cosThetaB, hDDB,(*effDDBpdf),"-nochi2",0,NULL,0,"cos#theta_{#Lambda}","Tot. eff.")->Draw();
		ceff->Print("DDeffFitB"+q2name+".pdf");
		GetFrame(cosThetaB, hLLB,(*effLLBpdf),"-nochi2",0,NULL,0,"cos#theta_{#Lambda}","Tot. eff.")->Draw();
		ceff->Print("LLeffFitB"+q2name+".pdf");
		delete ceff;
	}
}


void buildBkgPdfs(double q2min, double q2max, TString name, TCut cut, RooAbsPdf ** bkg, RooAbsPdf ** bkgB, TString model = "Poly2")
{
	RooRealVar * cosThetaL = new RooRealVar("cosThetaL","cosThetaL",0.,-1.,1.);
	RooRealVar * cosThetaB = new RooRealVar("cosThetaB","cosThetaB",0.,-1.,1.);

	TString q2name = ((TString)Form("_q2_%4.2f_%4.2f",q2min,q2max)).ReplaceAll(".","");
	string datafile = "/afs/cern.ch/work/p/pluca/weighted/Lmumu/trainingSamples.root";

	string treename = "bkgTestSample";
	if(q2min == 8 && q2max == 11)
	{
		q2name = "_jpsi";
		datafile = "/afs/cern.ch/work/p/pluca/weighted/Lmumu/Lb2Lmumu_CL_NBweighted.root";
		treename = "tree";
	}
	

	Analysis * ana = new Analysis("data"+name+q2name,"Lb",treename,datafile);
	ana->AddVariable(cosThetaL);
	ana->AddVariable(cosThetaB);
	ana->AddVariable("J_psi_1S_MM");
	ana->AddVariable("Lb_MM");

	TCut sideBand = "(Lb_MM > 5690 || Lb_MM < 5550) && weight > 0.8";
	if(q2min == 8 && q2max == 11) sideBand = "Lb_MM > 6000 && weight > 0.8 && TMath::Abs(J_psi_1S_MM - 3096.916) < 92.9";
	TCut q2cut = (TCut)Form("TMath::Power(J_psi_1S_MM/1000,2) > %e && TMath::Power(J_psi_1S_MM/1000,2) < %e",q2min,q2max);
	TCut mycut = sideBand + q2cut + cut;

	ana->applyCuts(&mycut);
	RooDataSet * data = ana->GetDataSet("-recalc");

	if(model=="RooKeyPdf")
	{
		(*bkg) = new RooKeysPdf("bkg"+name,"bkg"+name,*cosThetaL,*data,RooKeysPdf::MirrorBoth,1);
		(*bkgB) = new RooKeysPdf("bkg"+name+"B","bkg"+name+"B",*cosThetaB,*data,RooKeysPdf::MirrorBoth,1);
	}
	else
	{
		(*bkg) = stringToPdf(model,"bkg"+name, cosThetaL);
		(*bkg)->fitTo(*data,PrintLevel(-1));
		fixParams((*bkg),cosThetaL);
		(*bkgB) = stringToPdf(model,"bkg"+name+"B", cosThetaB); 
		(*bkgB)->fitTo(*data,PrintLevel(-1));
		fixParams((*bkgB),cosThetaB);
	}
}




Str2VarMap getJpsiPars(TString type, TCut cut, TFile * file)
{
	string datafilename = "/afs/cern.ch/work/p/pluca/weighted/Lmumu/candLb.root";
	string candfilename = "/afs/cern.ch/work/p/pluca/weighted/Lmumu/candLb_MC.root";
	TFile * MCFile = TFile::Open(candfilename.c_str());
	TTree * BdJpsiKSTree = (TTree *)MCFile->Get("candBdJpsiKS");

	file->cd();
	RooRealVar * vMuMu = new RooRealVar("Lb_MassConsLambda","Lb_MassConsLambda",5620.,5350.,6000.);
	Analysis * anaLbMuMu_MC = new Analysis("Lb2Lmumu_"+type+"_MC","Lb","candLb2Lmumu",candfilename,vMuMu,&CutsDef::cutMuMu);
	RooRealVar * vJpsi = new RooRealVar("Lb_MassConsJpsiLambda","Lb_MassConsJpsiLambda",5620.,5250.,6000.);
	Analysis * anaLbJpsi_default = new Analysis("Jpsi_default","Lb","candLb2JpsiL",datafilename,vJpsi,&CutsDef::cutJpsi);

	string model = "DCB_Sn";
	string optionsjpsi = "-ANDpulls-log-stdAxis-quiet-XM(#Lambda#mu#mu) (MeV/c^{2})-noCost-nochi2";
	string optionsmumu = "-stdAxis-quiet-XM(#Lambda#mu#mu) (MeV/c^{2})-noCost-nochi2";

	anaLbMuMu_MC->SetSignal((model+"-s[9]-s2[30]").c_str());
	anaLbMuMu_MC->Initialize("");
	anaLbMuMu_MC->Fit(5400.,5750.,200,true,"-stdAxis-XM(#Lambda#mu#mu) (MeV/c^{2})-noCost-nochi2",CutsDef::mumuTrueID+cut);
	Str2VarMap MCpars = anaLbMuMu_MC->GetSigParams();

	Analysis * KS = new Analysis("KS_bkg","Lb",BdJpsiKSTree,vJpsi,"DCB_OST","","-namepar");
	KS->Fit(5300.,6000.,150,true,"-noPlot-quiet");
	Str2VarMap pars = KS->GetSigParams();
	RooRealVar * m_shift = new RooRealVar("m_shift","shift",0.,-5.,5.);
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
	anaLbJpsi_default->Fit(5300.,6000.,200,true,optionsjpsi,cut);
	
	Str2VarMap jpsiSigpars = anaLbJpsi_default->GetSigParams();
	setConstant(&jpsiSigpars);
	RooRealVar * factor = new RooRealVar("factor","factor",1.,0.5,3.);
	ModifyPars(&jpsiSigpars,"s",factor);
	ModifyPars(&jpsiSigpars,"s2",factor);
	return jpsiSigpars;
}



RooDataSet * getDataAndFrac(TString name, TString q2name, TreeReader * mydata, TCut cut, RooRealVar * MM, double * frac, Str2VarMap jpsiPars, double *outnsig)
{
	RooRealVar * cosThetaL = new RooRealVar("cosThetaL","cosThetaL",0.,-1.,1.);
	RooRealVar * cosThetaB = new RooRealVar("cosThetaB","cosThetaB",0.,-1.,1.);
	TCut massCut = "Lb_MassConsLambda > 5590 && Lb_MassConsLambda < 5650";

	Analysis * ana = new Analysis(name+"_mass"+q2name,"Lb",mydata,&cut,MM);
	ana->AddVariable("J_psi_1S_MM");
	ana->AddVariable(cosThetaL);
	ana->AddVariable(cosThetaB);
	RooAbsPdf * mysig = stringToPdf("Gauss","sig",MM,jpsiPars);
	RooAbsPdf * mybkg = stringToPdf("Exp","bkgM",MM);
	RooRealVar * mynsig = new RooRealVar("mynsig","mynsig",50,0,100000);
	RooRealVar * mynbkg = new RooRealVar("mynbkg","mynbkg",10,0,100000);
	RooAbsPdf * Mmodel = new RooAddPdf("MassModel","MassModel",RooArgSet(*mysig,*mybkg),RooArgSet(*mynsig,*mynbkg));
	ana->applyCuts(&cut);
	RooDataSet * data = ana->GetDataSet("-recalc");
	Mmodel->fitTo(*data,Extended(kTRUE));
	
	double sigBkg = mybkg->createIntegral(*MM,NormSet(*MM),Range("Signal"))->getVal();
	double sig = mysig->createIntegral(*MM,NormSet(*MM),Range("Signal"))->getVal();
	double nsig = mynsig->getVal();
	double nbkg = mynbkg->getVal();
	if(frac)
	{
		frac[0] = nsig*sig/(nsig*sig+nbkg*sigBkg);
		frac[1] = frac[0]*TMath::Sqrt( TMath::Power(mynsig->getError()/nsig,2) + TMath::Power(mynbkg->getError()/nbkg,2) );
	}
	TCut mycut = cut + massCut;
	ana->applyCuts(&mycut);

	TCanvas * cc = new TCanvas();
	GetFrame(MM,Mmodel,data,"-nochi2-plotAllComp",30,NULL,0,"cos#theta_{#Lambda}")->Draw();
	cc->Print("M_"+name+"_"+q2name+".pdf");
	if(*outnsig) *outnsig = nsig;
	return ana->GetDataSet("-recalc");
}





int main(int argc, char **argv)
{
	bool printeff = true;
	string fc = "none";
	
	gROOT->ProcessLine(".x lhcbStyle.C");

	if(argc > 1)
	{
		for(int a = 1; a < argc; a++)
		{
			string arg = argv[a];
			string str = arg.substr(2,arg.length()-2);

			if(arg.find("-E")!=string::npos) fc = str;
			if(arg=="-peff") printeff = true;
		}
	}
	
	int nexp = 100;
	int nbins = 6;
	double q2min[] = {8.,15.,11.0,15,16,18};
	double q2max[] = {11.,20.,12.5,16,18,20};

	TString datafilename = "/afs/cern.ch/work/p/pluca/weighted/Lmumu/candLb.root";
	TreeReader * data = new TreeReader("candLb2Lmumu");
	data->AddFile(datafilename);
	TreeReader * datajpsi = new TreeReader("candLb2JpsiL");
	datajpsi->AddFile(datafilename);

	TFile * histFile = new TFile("Afb_bkgSys.root","recreate");

	string options = "-quiet-noPlot-lin-stdAxis-XM(#Lambda#mu#mu) (MeV/c^{2})-noCost-noParams";
	Analysis::SetPrintLevel("s");

	RooRealVar * cosThetaL = new RooRealVar("cosThetaL","cosThetaL",0.,-1.,1.);
	RooRealVar * cosThetaB = new RooRealVar("cosThetaB","cosThetaB",0.,-1.,1.);
	RooRealVar * MM = new RooRealVar("Lb_MassConsLambda","Lb_MassConsLambda",5621.,5400.,6000.);
	MM->setRange("Signal",5600,5640);
	RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

	//TGraphAsymmErrors * fL_vs_q2 = new TGraphAsymmErrors();
	//TCanvas * ceff = new TCanvas();

	RooCategory * samples = new RooCategory("samples","samples");
	samples->defineType("DD");
	samples->defineType("LL");

	RooRealVar * afb = new RooRealVar("afb","afb",0.,-0.75,0.75);
	RooRealVar * fL = new RooRealVar("fL","fL",0.6,0.,1.);
	TString afbLpdf = "((3./8.)*(1.-fL)*(1 + TMath::Power(cosThetaL,2)) + afb*cosThetaL + (3./4.)*fL*(1 - TMath::Power(cosThetaL,2)))";
	RooRealVar * afbB = new RooRealVar("afbB","afbB",0.,-0.5,0.5);
	TString afbBpdf = "(1 + 2*afbB*cosThetaB)";
	RooAbsPdf * teoPdf = new RooGenericPdf("teoPdf",afbLpdf,RooArgSet(*cosThetaL,*afb,*fL));
	RooAbsPdf * teoPdfB = new RooGenericPdf("teoPdfB",afbBpdf,RooArgSet(*cosThetaB,*afbB));

	TreeReader * mydata = datajpsi;
	Str2VarMap jpsiParsLL = getJpsiPars("LL", CutsDef::LLcut, histFile);
	Str2VarMap jpsiParsDD = getJpsiPars("DD", CutsDef::DDcut, histFile);

	vector<TH1 *> fLsysh, afbsysh, afbBsysh, fLsysh_frac, afbsysh_frac, afbBsysh_frac;

	for(int i = 0; i < nbins; i++)
	{
		TString q2name = ((TString)Form("q2_%4.2f_%4.2f",q2min[i],q2max[i])).ReplaceAll(".","");
		if(i>0) { mydata = data; MM->setRange(5400,6000); }
		else { q2name = "jpsi"; MM->setRange(5500,5850); }
		TString curq2cut = Form("TMath::Power(J_psi_1S_MM/1000,2) >= %e && TMath::Power(J_psi_1S_MM/1000,2) < %e",q2min[i],q2max[i]);	
		
		cout << "------------------- q2 bin: " << q2min[i] << " - " << q2max[i] << " -----------------------" << endl;

		/**               GET AND FIT EFFICIENCIES                  **/

		RooAbsPdf * effDDpdf = NULL, * effLLpdf = NULL, * effLLBpdf = NULL, * effDDBpdf = NULL;	
		getEfficiencies(q2min[i],q2max[i],&effLLpdf,&effDDpdf,&effLLBpdf,&effDDBpdf,printeff);
		cout << "Efficiencies extracted" << endl;
		histFile->cd();


		/**                    FIT AFB                  **/


		afb->setVal(0);
		afbB->setVal(-0.37);
		fL->setVal(0.6);

		RooAbsPdf * corrPdfLL = new RooProdPdf("sigPdfLL"+q2name,"corrPdfLL",*teoPdf,*effLLpdf);
		RooAbsPdf * corrPdfDD = new RooProdPdf("sigPdfDD"+q2name,"corrPdfDD",*teoPdf,*effDDpdf);
		RooAbsPdf * corrPdfLLB = new RooProdPdf("sigPdfLLB"+q2name,"corrPdfLLB",*teoPdfB,*effLLBpdf);
		RooAbsPdf * corrPdfDDB = new RooProdPdf("sigPdfDDB"+q2name,"corrPdfDDB",*teoPdfB,*effDDBpdf);

		TCut baseCut = "";
		TCut cutLL = CutsDef::LLcut + (TCut)curq2cut + baseCut;
		TCut cutDD = CutsDef::DDcut + (TCut)curq2cut + baseCut;

		histFile->cd();
		double fracDDv[2], fracLLv[2];
		double nsigDD, nsigLL;
		RooDataSet * dataLL = getDataAndFrac("LL",q2name,mydata,cutLL,MM,&fracLLv[0],jpsiParsLL,&nsigLL);
		RooDataSet * dataDD = getDataAndFrac("DD",q2name,mydata,cutDD,MM,&fracDDv[0],jpsiParsDD,&nsigDD);
		double nevts = nsigDD+nsigLL;

		cout << fixed << setprecision(3) << fracDDv[0] << "   " << fracDDv[1] << endl;
		RooRealVar * fracLL = new RooRealVar("fracLL","fracLL",fracLLv[0]);
		RooRealVar * fracDD = new RooRealVar("fracDD","fracDD",fracDDv[0]);

		RooAbsPdf * bkgLL = NULL, * bkgLLB = NULL, * bkgDD = NULL, * bkgDDB = NULL;
		buildBkgPdfs(q2min[i],q2max[i],"LL",CutsDef::LLcut,&bkgLL,&bkgLLB);
		buildBkgPdfs(q2min[i],q2max[i],"DD",CutsDef::DDcut,&bkgDD,&bkgDDB);
	
		cout << "Backgrounds extracted" << endl;

		RooAbsPdf * modelLL = new RooAddPdf("modelLL","modelLL",RooArgSet(*corrPdfLL,*bkgLL),*fracLL);
		RooAbsPdf * modelDD = new RooAddPdf("modelDD","modelDD",RooArgSet(*corrPdfDD,*bkgDD),*fracDD);
		RooAbsPdf * modelLLB = new RooAddPdf("modelLLB","modelLLB",RooArgSet(*corrPdfLLB,*bkgLLB),*fracLL);
		RooAbsPdf * modelDDB = new RooAddPdf("modelDDB","modelDDB",RooArgSet(*corrPdfDDB,*bkgDDB),*fracDD);

		// CREATE COMBINED DATASET
		RooDataSet * combData = new RooDataSet(Form("combData_%i",i),"combined data",RooArgSet(*MM,*cosThetaL,*cosThetaB),Index(*samples),Import("DD",*dataDD),Import("LL",*dataLL));

		Str2VarMap params;
		params["fL"] = fL;
		params["afb"] = afb;	
		Str2VarMap paramsB;
		paramsB["afbB"] = afbB;

		// FIT COS LEPTON
		RooSimultaneous * combModel = new RooSimultaneous(Form("combModel_%i",i),"",*samples);
		combModel->addPdf(*modelLL,"LL");
		combModel->addPdf(*modelDD,"DD");

		RooFitResult * res = safeFit(combModel,combData,params,&isInAllowedArea);	
	
		// FIT COS HADRON
		RooSimultaneous * combModelB = new RooSimultaneous(Form("combModelB_%i",i),"",*samples);
		combModelB->addPdf(*modelLLB,"LL");
		combModelB->addPdf(*modelDDB,"DD");

		RooFitResult * resB = safeFit(combModelB,combData,paramsB,&isInAllowedAreaB);

		cout << endl << fixed << setprecision(6) << "AfbB = " << afbB->getVal() << " +/- " << afbB->getError() << endl;
		cout << "Afb = " << afb->getVal() << " +/- " << afb->getError() << endl;
		cout << "fL = " << fL->getVal() << " +/- " << fL->getError() << endl;
		cout << endl;
		cout << "lepton:  " << res->edm() << "   "  << res->covQual() << endl;
		cout << "baryon:  " << resB->edm() << "   "  << resB->covQual() << endl;
		cout << endl;

		TH1F * fLsys = new TH1F(Form("fLsys_%i",i),"fLsys",40,-1,1);
		TH1F * afbsys = new TH1F(Form("afbsys_%i",i),"afbsys",40,-1,1);
		TH1F * afbBsys = new TH1F(Form("afbBsys_%i",i),"afbBsys",40,-1,1);
		TH1F * fLsys_frac = new TH1F(Form("fLsys_frac%i",i),"fLsys",40,-1,1);
		TH1F * afbsys_frac = new TH1F(Form("afbsys_frac%i",i),"afbsys",40,-1,1);
		TH1F * afbBsys_frac = new TH1F(Form("afbBsys_frac%i",i),"afbBsys",40,-1,1);


		RooAbsPdf * mybkgDD_2 = NULL, * mybkgDDB_2 = NULL;
		buildBkgPdfs(q2min[i],q2max[i],"DD",CutsDef::DDcut,&mybkgDD_2,&mybkgDDB_2,"RooKeyPdf");

		//cout << nevts << endl;
		//TRandom3 r(0);

		for(int e = 0; e < nexp; e++)
		{
			histFile->cd();
			RooAbsPdf * toypdf = (RooAbsPdf *)modelDD->Clone();
			Analysis * toy = new Analysis("toy",cosThetaL,modelDD,nevts);
			RooAbsPdf * toypdfB = (RooAbsPdf *)modelDDB->Clone();
			Analysis * toyB = new Analysis("toyB",cosThetaB,modelDDB,nevts);
			
			afb->setVal(0);
			afbB->setVal(-0.37);
			fL->setVal(0.6);

			safeFit(toypdf,toy->GetDataSet("-recalc"),params,&isInAllowedArea);
			safeFit(toypdfB,toyB->GetDataSet("-recalc"),paramsB,&isInAllowedAreaB);
			double def_afb = afb->getVal();
			double def_fL = fL->getVal();
			double def_afbB = afbB->getVal();

			afb->setVal(0);
			afbB->setVal(-0.37);
			fL->setVal(0.6);

			RooAbsPdf * modelDD_2 = new RooAddPdf("modelDD_2","modelDD",RooArgSet(*corrPdfDD,*mybkgDD_2),*fracDD);
			RooAbsPdf * modelDDB_2 = new RooAddPdf("modelDDB_2","modelDDB",RooArgSet(*corrPdfDDB,*mybkgDDB_2),*fracDD);
			safeFit(modelDD_2,toy->GetDataSet("-recalc"),params,&isInAllowedArea);
			safeFit(modelDDB_2,toyB->GetDataSet("-recalc"),paramsB,&isInAllowedAreaB);
			double oth_afb = afb->getVal();
			double oth_fL = fL->getVal();
			double oth_afbB = afbB->getVal();

			fLsys->Fill(oth_fL-def_fL);
			afbsys->Fill(oth_afb-def_afb);
			afbBsys->Fill(oth_afbB-def_afbB);
			

			afb->setVal(0.);
			afbB->setVal(-0.37);
			fL->setVal(0.6);

			//double rdm_frac = r.Gaus(fracDDv[0],fracDDv[1]);
			double rdm_frac = fracDDv[0] + fracDDv[1];
			RooRealVar * fracDD_2 = new RooRealVar("fracDD_2","fracDD_2",rdm_frac);	
			RooAbsPdf * modelDD_3 = new RooAddPdf("modelDD_3","modelDD",RooArgSet(*corrPdfDD,*bkgDD),*fracDD_2);
			RooAbsPdf * modelDDB_3 = new RooAddPdf("modelDDB_3","modelDDB",RooArgSet(*corrPdfDDB,*bkgDDB),*fracDD_2);
			safeFit(modelDD_3,toy->GetDataSet("-recalc"),params,&isInAllowedArea);
			safeFit(modelDDB_3,toyB->GetDataSet("-recalc"),paramsB,&isInAllowedAreaB);

			double frc_afb = afb->getVal();
			double frc_fL = fL->getVal();
			double frc_afbB = afbB->getVal();

			fLsys_frac->Fill(frc_fL-def_fL);
			afbsys_frac->Fill(frc_afb-def_afb);
			afbBsys_frac->Fill(frc_afbB-def_afbB);
			
		}

		afbsysh.push_back(afbsys);
		afbBsysh.push_back(afbBsys);
		fLsysh.push_back(fLsys);
		afbsysh_frac.push_back(afbsys_frac);
		afbBsysh_frac.push_back(afbBsys_frac);
		fLsysh_frac.push_back(fLsys_frac);

	}

	
	for(int q = 0; q < nbins; q++)
	{
		cout << fixed << setprecision(2) << "-------- Bin " << q2min[q] << "-" << q2max[q] << endl;
		cout << fixed << setprecision(5) << "fL sys = " << fLsysh[q]->GetMean() << " +/- " << fLsysh[q]->GetMeanError() << endl;
		cout << "Afb sys = " << afbsysh[q]->GetMean() << " +/- " << afbsysh[q]->GetMeanError() << endl;
		cout << "AfbB sys = " << afbBsysh[q]->GetMean() << " +/- " << afbBsysh[q]->GetMeanError() << endl;
	}

	cout << "#################################################################" << endl;
	for(int q = 0; q < nbins; q++)
	{
		cout << fixed << setprecision(2) << "-------- Bin " << q2min[q] << "-" << q2max[q] << endl;
		cout << fixed << setprecision(5) << "fL sys = " << fLsysh_frac[q]->GetMean() << " +/- " << fLsysh_frac[q]->GetMeanError() << endl;
		cout << "Afb sys = " << afbsysh_frac[q]->GetMean() << " +/- " << afbsysh_frac[q]->GetMeanError() << endl;
		cout << "AfbB sys = " << afbBsysh_frac[q]->GetMean() << " +/- " << afbBsysh_frac[q]->GetMeanError() << endl;
	}

	cout << "#################################################################" << endl;
	for(int q = 0; q < nbins; q++)
	{
		cout << fixed << setprecision(2) << "-------- Bin " << q2min[q] << "-" << q2max[q] << endl;
		cout << fixed << setprecision(5) << "fL sys = " << TMath::Sqrt(TMath::Power(fLsysh_frac[q]->GetMean(),2) + TMath::Power(fLsysh[q]->GetMean(),2) )  << endl;
		cout << "Afb sys = " << TMath::Sqrt(TMath::Power(afbsysh_frac[q]->GetMean(),2) + TMath::Power(afbsysh[q]->GetMean(),2) ) << endl;
		cout << "AfbB sys = " << TMath::Sqrt(TMath::Power(afbBsysh_frac[q]->GetMean(),2) + TMath::Power(afbBsysh[q]->GetMean(),2) ) << endl;
	}

}
