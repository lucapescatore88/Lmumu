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
		effDDB  = (TH1F *)effFile->Get("htoteff");
		effFile = TFile::Open(effbase+"Lbeff_JpsivscosThetaB_LL.root");
		effLLB  = (TH1F *)effFile->Get("htoteff");
	}
	else if((q2min == 15 && q2max == 20) || (q2min == 1.1 && q2max == 6))
	{
		TFile * effFile = TFile::Open(effbase+"Lbeff2D_cosThetaL_vs_q2_DD_2bins.root");
		TH2F * effDD2D  = (TH2F *)effFile->Get("htot_eff");
		effDD = (TH1F*)GetSliceX(effDD2D,(q2max+q2min)/2.);
		effFile = TFile::Open(effbase+"Lbeff2D_cosThetaL_vs_q2_LL_2bins.root");
		TH2F * effLL2D  = (TH2F *)effFile->Get("htot_eff");
		effLL = (TH1F*)GetSliceX(effLL2D,(q2max+q2min)/2.);
		effFile = TFile::Open(effbase+"Lbeff2D_cosThetaB_vs_q2_DD_2bins.root");
		TH2F * effDDB2D  = (TH2F *)effFile->Get("htot_eff");
		effDDB = (TH1F*)GetSliceX(effDDB2D,(q2max+q2min)/2.);
		effFile = TFile::Open(effbase+"Lbeff2D_cosThetaB_vs_q2_LL_2bins.root");
		TH2F * effLLB2D  = (TH2F *)effFile->Get("htot_eff");
		effLLB = (TH1F*)GetSliceX(effLLB2D,(q2max+q2min)/2.);
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
		TH2F * effDDB2D  = (TH2F *)effFile->Get("htot_eff");
		effDDB = (TH1F*)GetSliceX(effDDB2D,(q2max+q2min)/2.);
		effFile = TFile::Open(effbase+"Lbeff2D_cosThetaB_vs_q2_LL.root");
		TH2F * effLLB2D  = (TH2F *)effFile->Get("htot_eff");
		effLLB = (TH1F*)GetSliceX(effLLB2D,(q2max+q2min)/2.);
	}

	RooDataHist * hLL = new RooDataHist("hLL","hLL",*cosThetaL,effLL);
	RooDataHist * hDD = new RooDataHist("hDD","hDD",*cosThetaL,effDD);
	(*effLLpdf) = stringToPdf("Poly1-v0-l1","effLL", cosThetaL);
	(*effDDpdf) = stringToPdf("Poly2-v0-l1","effDD", cosThetaL);
	(*effLLpdf)->fitTo(*hLL,PrintLevel(-1));
	(*effDDpdf)->fitTo(*hDD,PrintLevel(-1));
	fixParams((*effLLpdf),cosThetaL);
	fixParams((*effDDpdf),cosThetaL);	

	RooDataHist * hLLB = new RooDataHist("hLLB","hLLB",*cosThetaB,effLLB);
	RooDataHist * hDDB = new RooDataHist("hDDB","hDDB",*cosThetaB,effDDB);
	(*effLLBpdf) = stringToPdf("Poly1-v-0.2-l1","effLLB", cosThetaB);
	(*effDDBpdf) = stringToPdf("Poly2-v-0.2-l1","effDDB", cosThetaB);
	(*effLLBpdf)->fitTo(*hLLB,PrintLevel(-1));//,Range(-0.9,1));
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


void buildBkgPdfs(double q2min, double q2max, TString name, TCut cut, RooAbsPdf ** bkg, RooAbsPdf ** bkgB, TString model = "Poly2-v0-l2")
{
	RooRealVar * cosThetaL = new RooRealVar("cosThetaL","cosThetaL",0.,-1.,1.);
	RooRealVar * cosThetaB = new RooRealVar("cosThetaB","cosThetaB",0.,-1.,1.);

	TString q2name = ((TString)Form("_q2_%4.2f_%4.2f",q2min,q2max)).ReplaceAll(".","");
	string datafile = "/afs/cern.ch/work/p/pluca/weighted/Lmumu/trainingSamples.root";

	string treename = "bkgTestSample";
	int nbins = 10;
	if(q2min == 8 && q2max == 11)
	{
		q2name = "_jpsi";
		datafile = "/afs/cern.ch/work/p/pluca/weighted/Lmumu/Lb2Lmumu_CL_NBweighted.root";
		treename = "tree";
		nbins = 20;
	}
	

	Analysis * ana = new Analysis("data"+name+q2name,"Lb",treename,datafile);
	ana->AddVariable(cosThetaL);
	ana->AddVariable(cosThetaB);
	ana->AddVariable("J_psi_1S_MM");
	ana->AddVariable("Lb_MM");
	ana->AddVariable("Lb_MassConsLambda");
	ana->AddVariable("Lb_MassConsJpsiAndLambda");

	TCut sideBand = "Lb_MassConsLambda > 5670 && weight > 0.6";
	if(q2min == 8 && q2max == 11) sideBand = "Lb_MassConsLambda > 6000 && weight > 0.8 && TMath::Abs(J_psi_1S_MM - 3096.916) < 92.9";
	TCut q2cut = (TCut)Form("TMath::Power(J_psi_1S_MM/1000,2) > %e && TMath::Power(J_psi_1S_MM/1000,2) < %e",q2min,q2max);
	TCut mycut = sideBand + q2cut + cut;

	ana->applyCuts(&mycut);
	RooDataSet * data = ana->GetDataSet("-recalc");

	(*bkg) = stringToPdf(model,"bkg"+name, cosThetaL);
	(*bkg)->fitTo(*data,PrintLevel(-1));
	fixParams((*bkg),cosThetaL);
	
	//(*bkgB) = new RooKeysPdf("bkg"+name+"B","bkg"+name+"B",*cosThetaB,*data,RooKeysPdf::MirrorBoth,1);
	(*bkgB) = stringToPdf(model,"bkg"+name+"B", cosThetaB); 
	(*bkgB)->fitTo(*data,PrintLevel(-1));
	fixParams((*bkgB),cosThetaB);
	
	TCanvas * ceff = new TCanvas();
	GetFrame(cosThetaL, data,(*bkg),"-nochi2",nbins,NULL,0,"cos#theta_{l}","N bkg")->Draw();
	ceff->Print("BkgFit"+q2name+"_"+name+".pdf");
	GetFrame(cosThetaB, data,(*bkgB),"-nochi2",nbins,NULL,0,"cos#theta_{#Lambda}","N bkg")->Draw();
	ceff->Print("BkgFitB"+q2name+"_"+name+".pdf");
	delete ceff;
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
	RooRealVar * vJpsi = new RooRealVar("Lb_MassConsJpsiLambda","Lb_MassConsJpsiLambda",5620.,5400.,6000.);
	Analysis * anaLbJpsi_default = new Analysis("Jpsi_default","Lb","candLb2JpsiL",datafilename,vJpsi,&CutsDef::cutJpsi);

	string model = "DCB_Sn";
	string optionsjpsi = "-ANDpulls-log-stdAxis-quiet-XM(#Lambda#mu#mu) (MeV/c^{2})-noCost-nochi2";
	string optionsmumu = "-stdAxis-quiet-XM(#Lambda#mu#mu) (MeV/c^{2})-noCost-nochi2";

	anaLbMuMu_MC->SetSignal((model+"-s[9]-s2[30]").c_str());
	anaLbMuMu_MC->Initialize("");
	anaLbMuMu_MC->Fit(5400.,6000.,200,true,"-stdAxis-XM(#Lambda#mu#mu) (MeV/c^{2})-noCost-nochi2",CutsDef::mumuTrueID+cut);
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



RooDataSet * getDataAndFrac(TString name, TString q2name, TreeReader * mydata, TCut cut, RooRealVar * MM, double * frac, Str2VarMap jpsiPars)
{
	RooRealVar * cosThetaL = new RooRealVar("cosThetaL","cosThetaL",0.,-1.,1.);
	RooRealVar * cosThetaB = new RooRealVar("cosThetaB","cosThetaB",0.,-1.,1.);
	TCut massCut = "Lb_MassConsLambda > 5590 && Lb_MassConsLambda < 5650";// (TCut)((TString)MM->GetName()+" > 5590 && "+(TString)MM->GetName()+" < 5650");

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
	(*frac) = nsig*sig/(nsig*sig+nbkg*sigBkg);
	TCut mycut = cut + massCut;
	ana->applyCuts(&mycut);

	TCanvas * cc = new TCanvas();
	GetFrame(MM,Mmodel,data,"-nochi2-plotAllComp",30,NULL,0,"cos#theta_{#Lambda}")->Draw();
	cc->Print("M_"+name+"_"+q2name+".pdf");
	return ana->GetDataSet("-recalc");
}





int main(int argc, char **argv)
{
	bool printeff = true;
	string fc = "none";
	int dobin = -1;
	int gradbkg = 0;

	gROOT->ProcessLine(".x lhcbStyle.C");

	if(argc > 1)
	{
		for(int a = 1; a < argc; a++)
		{
			string arg = argv[a];
			string str = arg.substr(2,arg.length()-2);

			if(arg.find("-E")!=string::npos) fc = str;
			if(arg=="-peff") printeff = true;
			if(arg.find("-B")!=string::npos) dobin = ((TString)str).Atof();
			if(arg.find("-g")!=string::npos) gradbkg = ((TString)str).Atof();
		}
	}

	int start = 0;
	int nbins = 6;//CutsDef::nq2bins;
	if(dobin != -1) { start = dobin; nbins = dobin+1; }	
	double q2min[] = {8.,15.,11.0,15,16,18};//&CutsDef::q2min_highfirst[0];
	double q2max[] = {11.,20.,12.5,16,18,20};//&CutsDef::q2max_highfirst[0];

	TString datafilename = "/afs/cern.ch/work/p/pluca/weighted/Lmumu/candLb.root";
	TreeReader * data = new TreeReader("candLb2Lmumu");
	data->AddFile(datafilename);
	TreeReader * datajpsi = new TreeReader("candLb2JpsiL");
	datajpsi->AddFile(datafilename);

	TFile * histFile = new TFile("Afb_hist.root","recreate");

	string options = "-quiet-noPlot-lin-stdAxis-XM(#Lambda#mu#mu) (MeV/c^{2})-noCost-noParams";
	Analysis::SetPrintLevel("s");

	RooRealVar * cosThetaL = new RooRealVar("cosThetaL","cosThetaL",0.,-1.,1.);
	RooRealVar * cosThetaB = new RooRealVar("cosThetaB","cosThetaB",0.,-1.,1.);
	RooRealVar * MM = new RooRealVar("Lb_MassConsLambda","Lb_MassConsLambda",5621.,5400.,6000.);
	MM->setRange("Signal",5600,5640);
	RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

	TGraphAsymmErrors * Afb_vs_q2 = new TGraphAsymmErrors();
	TGraphAsymmErrors * AfbB_vs_q2 = new TGraphAsymmErrors();
	TGraphAsymmErrors * fL_vs_q2 = new TGraphAsymmErrors();
	TCanvas * ceff = new TCanvas();

	RooCategory * samples = new RooCategory("samples","samples");
	samples->defineType("DD");
	samples->defineType("LL");

	RooRealVar * afb = new RooRealVar("afb","afb",0.,-0.75,0.75);
	RooRealVar * fL = new RooRealVar("fL","fL",0.6,0.,1.);
	RooRealVar * origafb = new RooRealVar("afb","afb",0.,-0.75,0.75);
	RooRealVar * origfL = new RooRealVar("fL","fL",0.6,0.,1.);
	TString afbLpdf = "((3./8.)*(1.-fL)*(1 + TMath::Power(cosThetaL,2)) + afb*cosThetaL + (3./4.)*fL*(1 - TMath::Power(cosThetaL,2)))";
	RooRealVar * afbB = new RooRealVar("afbB","afbB",-0.3,-0.5,0.);
	RooRealVar * origafbB = new RooRealVar("afbB","afbB",-0.3,-0.5,0.);
	TString afbBpdf = "(1 + 2*afbB*cosThetaB)";
	RooAbsPdf * teoPdf = new RooGenericPdf("teoPdf",afbLpdf,RooArgSet(*cosThetaL,*afb,*fL));
	RooAbsPdf * teoPdfB = new RooGenericPdf("teoPdfB",afbBpdf,RooArgSet(*cosThetaB,*afbB));

	TreeReader * mydata = datajpsi;
	Str2VarMap jpsiParsLL = getJpsiPars("LL", CutsDef::LLcut, histFile);
	Str2VarMap jpsiParsDD = getJpsiPars("DD", CutsDef::DDcut, histFile);
	
	for(int i = start; i < nbins; i++)
	{
		TString q2name = ((TString)Form("q2_%4.2f_%4.2f",q2min[i],q2max[i])).ReplaceAll(".","");
		if(q2min[i] == 8 && q2max[i] == 11) { q2name = "jpsi"; MM->setRange(5500,5850); }
		else { mydata = data; MM->setRange(5400,6000); }
		
		TString curq2cut = Form("TMath::Power(J_psi_1S_MM/1000,2) >= %e && TMath::Power(J_psi_1S_MM/1000,2) < %e",q2min[i],q2max[i]);	
		
		cout << "------------------- q2 bin: " << q2min[i] << " - " << q2max[i] << " -----------------------" << endl;

		/**               GET AND FIT EFFICIENCIES                  **/

		RooAbsPdf * effDDpdf = NULL, * effLLpdf = NULL, * effLLBpdf = NULL, * effDDBpdf = NULL;	
		getEfficiencies(q2min[i],q2max[i],&effLLpdf,&effDDpdf,&effLLBpdf,&effDDBpdf,printeff);
		cout << "Efficiencies extracted" << endl;
		histFile->cd();


		/**                    FIT AFB                  **/


		afb->setVal(0);
		afbB->setVal(-0.3);
		fL->setVal(0.6);

		RooAbsPdf * corrPdfLL = new RooProdPdf("sigPdfLL"+q2name,"corrPdfLL",*teoPdf,*effLLpdf);
		RooAbsPdf * corrPdfDD = new RooProdPdf("sigPdfDD"+q2name,"corrPdfDD",*teoPdf,*effDDpdf);
		RooAbsPdf * corrPdfLLB = new RooProdPdf("sigPdfLLB"+q2name,"corrPdfLLB",*teoPdfB,*effLLBpdf);
		RooAbsPdf * corrPdfDDB = new RooProdPdf("sigPdfDDB"+q2name,"corrPdfDDB",*teoPdfB,*effDDBpdf);

		TCut baseCut = "";
		TCut cutLL = CutsDef::LLcut + (TCut)curq2cut + baseCut;
		TCut cutDD = CutsDef::DDcut + (TCut)curq2cut + baseCut;

		histFile->cd();
		double fracDDv = 0, fracLLv = 0;
		RooDataSet * dataLL = getDataAndFrac("LL",q2name,mydata,cutLL,MM,&fracLLv,jpsiParsLL);
		RooDataSet * dataDD = getDataAndFrac("DD",q2name,mydata,cutDD,MM,&fracDDv,jpsiParsDD);

		RooRealVar * fracLL = new RooRealVar("fracLL","fracLL",fracLLv);
		RooRealVar * fracDD = new RooRealVar("fracDD","fracDD",fracDDv);

		//RooAbsPdf * bkgLL = NULL, * bkgLLB = NULL, * bkgDD = NULL, * bkgDDB = NULL;
		//buildBkgPdfs(q2min[i],q2max[i],"LL",CutsDef::LLcut,&bkgLL,&bkgLLB);
		//buildBkgPdfs(q2min[i],q2max[i],"DD",CutsDef::DDcut,&bkgDD,&bkgDDB);

		string bkgmodel = Form("Poly%i-v0-l2",gradbkg);
		RooAbsPdf * bkgLL_teo = stringToPdf(bkgmodel.c_str(),"bkgLL_teo", cosThetaL);
		RooAbsPdf * bkgLLB_teo = stringToPdf(bkgmodel.c_str(),"bkgLLB_teo", cosThetaB);
		RooAbsPdf * bkgDD_teo = stringToPdf(bkgmodel.c_str(),"bkgDD_teo", cosThetaL); 
		RooAbsPdf * bkgDDB_teo = stringToPdf(bkgmodel.c_str(),"bkgDDB_teo", cosThetaB);
		RooAbsPdf * bkgLL = new RooProdPdf("bkgLL","",*bkgLL_teo,*effLLpdf);
		RooAbsPdf * bkgDD = new RooProdPdf("bkgDD","",*bkgDD_teo,*effDDpdf);
		RooAbsPdf * bkgLLB = new RooProdPdf("bkgLLB","",*bkgLLB_teo,*effLLBpdf);
		RooAbsPdf * bkgDDB = new RooProdPdf("bkgDDB","",*bkgDDB_teo,*effDDBpdf);


		cout << "Backgrounds extracted" << endl;

		RooAbsPdf * modelLL = new RooAddPdf("modelLL","modelLL",RooArgSet(*corrPdfLL,*bkgLL),*fracLL);
		RooAbsPdf * modelDD = new RooAddPdf("modelDD","modelDD",RooArgSet(*corrPdfDD,*bkgDD),*fracDD);
		RooAbsPdf * modelLLB = new RooAddPdf("modelLLB","modelLLB",RooArgSet(*corrPdfLLB,*bkgLLB),*fracLL);
		RooAbsPdf * modelDDB = new RooAddPdf("modelDDB","modelDDB",RooArgSet(*corrPdfDDB,*bkgDDB),*fracDD);

		// CREATE COMBINED DATASET
		RooDataSet * combData = new RooDataSet(Form("combData_%i",i),"combined data",RooArgSet(*MM,*cosThetaL,*cosThetaB),Index(*samples),Import("DD",*dataDD),Import("LL",*dataLL));


		// FIT COS LEPTON
		RooSimultaneous * combModel = new RooSimultaneous(Form("combModel_%i",i),"",*samples);
		combModel->addPdf(*modelLL,"LL");
		combModel->addPdf(*modelDD,"DD");

		Str2VarMap params;
		params["fL"] = fL;
		params["afb"] = afb;	
		Str2VarMap paramsB;
		paramsB["afbB"] = afbB;

		RooFitResult * res = safeFit(combModel,combData,params,&isInAllowedArea,"-scan",getNFreePars(combModel,RooArgSet(*cosThetaL)));
			
		int nbinsLL = 8;
		int nbinsDD = 12;
		if(q2name=="jpsi") { nbinsLL = 40; nbinsDD = 40; }

		ceff->cd();
		string opts = "-noCost-plotSigComp-fillBkg";
		//if(q2name=="jpsi") opts += "-layout 65 99 70";
		GetFrame(cosThetaL,modelLL,dataLL,opts,nbinsLL,NULL,0,"cos#theta_{l}")->Draw();
		ceff->Print("Afb_LL_"+q2name+".pdf");
		GetFrame(cosThetaL,modelDD,dataDD,opts,nbinsDD,NULL,0,"cos#theta_{l}")->Draw();
		ceff->Print("Afb_DD_"+q2name+".pdf");
	
		Afb_vs_q2->SetPoint(i,(q2max[i] + q2min[i])/2.,afb->getVal());
		fL_vs_q2->SetPoint(i,(q2max[i] + q2min[i])/2.,fL->getVal());

		// FIT COS HADRON
		RooSimultaneous * combModelB = new RooSimultaneous(Form("combModelB_%i",i),"",*samples);
		combModelB->addPdf(*modelLLB,"LL");
		combModelB->addPdf(*modelDDB,"DD");

		//RooFitResult * resB = combModelB->fitTo(*combData,Save());
		RooFitResult * resB = safeFit(combModelB,combData,paramsB,&isInAllowedAreaB,"-scan",getNFreePars(combModelB,RooArgSet(*cosThetaB)));
		
		ceff->cd();
		GetFrame(cosThetaB,modelLLB,dataLL,"-noCost-plotSigComp-fillBkg",nbinsLL,NULL,0,"cos#theta_{#Lambda}")->Draw();
		ceff->Print("AfbB_LL_"+q2name+".pdf");
		GetFrame(cosThetaB,modelDDB,dataDD,"-noCost-plotSigComp-fillBkg",nbinsDD,NULL,0,"cos#theta_{#Lambda}")->Draw();
		ceff->Print("AfbB_DD_"+q2name+".pdf");
			
		AfbB_vs_q2->SetPoint(i,(q2max[i] + q2min[i])/2.,afbB->getVal());

		cout << endl << fixed << setprecision(6) << "AfbB = " << afbB->getVal() << " +/- " << afbB->getError() << endl;
		cout << "Afb = " << afb->getVal() << " +/- " << afb->getError() << endl;
		cout << "fL = " << fL->getVal() << " +/- " << fL->getError() << endl;
		cout << endl;
		if(res) cout << "lepton:  " << res->edm() << "   "  << res->covQual() << endl;
		if(resB) cout << "baryon:  " << resB->edm() << "   "  << resB->covQual() << endl;
		cout << endl;


		if(fc!="none")
		{ 
			cout << "------------------------ FELDMAN AND COUSINS ------------------------" << endl;

			vector<double> afb_err(0,2), fL_err(0,2), afbB_err(0,2);
			vector < RooDataSet * > datas;
			vector < RooAbsPdf * > pdfs, pdfsB;
			vector < TString > cat;
			cat.push_back("LL");
			cat.push_back("DD");
			datas.push_back(dataLL);
			datas.push_back(dataDD);

			RooArgSet * origPars = new RooArgSet(*origafb,*origfL);
			pdfs.push_back(modelLL);
			pdfs.push_back(modelDD);

			//double fLval = fL->getVal(), fLerr = fL->getError();
			//double afbval = afb->getVal();
			double afbBval = afbB->getVal();
/*
			if(fc == "all" || fc == "afb")
			{
				FeldmanCousins * FC = new FeldmanCousins(q2name,cat,datas,pdfs,cosThetaL,afb,"nsig_sw");
				if(q2min[i]==11) afb_err = FC->ExtractLimits(params,origPars,-0.33,-0.09,&isInAllowedArea);
				else if( (afb->getVal()-1.4*afb->getError()) > -0.75 && (afb->getVal()+1.4*afb->getError()) < 0.75 )
					afb_err = FC->ExtractLimits(params,origPars,afb->getVal()-1.4*afb->getError(),afb->getVal()+1.4*afb->getError(),&isInAllowedArea);
				else afb_err = FC->ExtractLimits(params,origPars,-0.5,0.5,&isInAllowedArea);
				Afb_vs_q2->SetPointError(i,(q2max[i] - q2min[i])/2.,(q2max[i] - q2min[i])/2.,TMath::Abs(afb_err[0]-afbval),TMath::Abs(afb_err[1]-afbval));
			}
			if(fc == "all" || fc == "fL")
			{
				FeldmanCousins * FCfL = new FeldmanCousins(q2name,cat,datas,pdfs,cosThetaL,fL,"nsig_sw");
				if(q2min[i]==11) fL_err = FCfL->ExtractLimits(params,origPars,0.45,0.85,&isInAllowedArea);
				else fL_err = FCfL->ExtractLimits(params,origPars,fLval-1.5*fLerr,fLval+1.5*fLerr,&isInAllowedArea);
				fL_vs_q2->SetPointError(i,(q2max[i] - q2min[i])/2.,(q2max[i] - q2min[i])/2.,TMath::Abs(fL_err[0]-fLval),TMath::Abs(fL_err[1]-fLval));
			}
*/
			if(fc == "all" || fc == "afb")
			{
				FeldmanCousins * FC = new FeldmanCousins(q2name,cat,datas,pdfs,cosThetaL,afb,"",fL);
				FC->SetNExp(100);
				FC->ExtractLimits(params,origPars,-0.3,0.3,0.2,1,&isInAllowedArea);
			}

			double lims[6][2] = { {-0.5,-0.2}, {-0.47,-0.2}, {-0.49,-0.2}, {-0.4,0.1}, {-0.5,-0.2}, {-0.43,0.} };
			if(fc == "all" || fc == "afbB")
			{
				RooArgSet * origParsB = new RooArgSet(*origafbB);
				pdfsB.push_back(modelLLB);
				pdfsB.push_back(modelDDB);
				FeldmanCousins * FCB = new FeldmanCousins(q2name,cat,datas,pdfsB,cosThetaB,afbB,"nsig_sw");
				FCB->SetNExp(200);
				FCB->SetNPointsToScan(15);
				afbB_err = FCB->ExtractLimits(paramsB,origParsB,lims[i][0],lims[i][1],&isInAllowedAreaB);
				AfbB_vs_q2->SetPointError(i,(q2max[i] - q2min[i])/2.,(q2max[i] - q2min[i])/2.,TMath::Abs(afbB_err[0]-afbBval),TMath::Abs(afbB_err[1]-afbBval));
			}
		}
		else
		{
			AfbB_vs_q2->SetPointError(i,(q2max[i] - q2min[i])/2.,(q2max[i] - q2min[i])/2.,0,0);
			Afb_vs_q2->SetPointError(i,(q2max[i] - q2min[i])/2.,(q2max[i] - q2min[i])/2.,0,0);
			fL_vs_q2->SetPointError(i,(q2max[i] - q2min[i])/2.,(q2max[i] - q2min[i])/2.,0,0);
		}
	}

	ceff->cd();
	Afb_vs_q2->GetXaxis()->SetTitle("q^{2}");
	Afb_vs_q2->GetYaxis()->SetTitle("A_{FB}^{l}");
	Afb_vs_q2->SetMaximum(1);
	Afb_vs_q2->SetMinimum(-1);
	Afb_vs_q2->Draw("AP");
	ceff->Print("Afb_vs_q2.pdf");
	AfbB_vs_q2->GetXaxis()->SetTitle("q^{2}");
	AfbB_vs_q2->GetYaxis()->SetTitle("A_{FB}^{h}");
	AfbB_vs_q2->SetMaximum(1);
	AfbB_vs_q2->SetMinimum(-1);
	AfbB_vs_q2->Draw("AP");
	ceff->Print("AfbB_vs_q2.pdf");
	fL_vs_q2->GetXaxis()->SetTitle("q^{2}");
	fL_vs_q2->GetYaxis()->SetTitle("f_{L}");
	fL_vs_q2->Draw("AP");
	ceff->Print("fL_vs_q2.pdf");

	for(int bb = 0; bb < Afb_vs_q2->GetN(); bb++)
	{
		double qq, afbv, afbBv, fLv;
		Afb_vs_q2->GetPoint(bb,qq,afbv);
		AfbB_vs_q2->GetPoint(bb,qq,afbBv);
		fL_vs_q2->GetPoint(bb,qq,fLv);
		cout << fixed << setprecision(1) << q2min[bb] << " - " << q2max[bb];
		cout << fixed << setprecision(4); 

		if(fc == "all" || fc == "afb") cout << " & $" << afbv << "_{-" << Afb_vs_q2->GetErrorYlow(bb) << "}^{+" << Afb_vs_q2->GetErrorYhigh(bb)  << "} \\text{(stat)} \\pm \\text{(sys)}$ ";
		if(fc == "all" || fc == "afbB") cout << " & $" << afbBv << "_{-" << AfbB_vs_q2->GetErrorYlow(bb) << "}^{+" << AfbB_vs_q2->GetErrorYhigh(bb) << "} \\text{(stat)} \\pm \\text{(sys)}$ " ;
		if(fc == "all" || fc == "fL") cout << " & $" << fLv << "_{-" << fL_vs_q2->GetErrorYlow(bb) << "}^{+" << fL_vs_q2->GetErrorYhigh(bb)  << "} \\text{(stat)} \\pm \\text{(sys)}$ ";
		if(fc == "none") cout << " & \t" << afbv << " & \t" << afbBv << " & \t" << fLv;
		cout << "  \\\\ " << endl;
	}

}
