#ifndef RUNFUNCTIONS
#define RUNFUNCTIONS


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

namespace AngVars
{
	RooRealVar * cosThetaL = new RooRealVar("cosThetaL","cosThetaL",0.,-1.,1.);
	RooRealVar * cosThetaB = new RooRealVar("cosThetaB","cosThetaB",0.,-1.,1.);
}


void constrainParams(RooAbsPdf * pdf, RooRealVar * obs, double val, string opt = "")
{
	RooArgSet * params = pdf->getParameters(RooDataSet("dataset","",*obs));
	TIterator * it = params->createIterator();
	RooRealVar * arg;
	while( (arg = (RooRealVar*)it->Next()) )  
	{
		arg->setMin(arg->getVal()*(1.-val));
		arg->setMax(arg->getVal()*(1.+val));
		//arg->Print();
	}
	delete it;
}




vector < double > getList(string str)
{
	vector < string > v;
	vector < double > res;
	unsigned pos = 0;
	while(true)
	{
		unsigned pos2 = str.find(",",pos+1);
		v.push_back(str.substr(pos+1,pos2-pos));
		if(pos2 > 1.e9) break;
		pos = pos2;
	}
	
	for(unsigned i = 0; i < v.size(); i++) res.push_back(((TString)v[i]).Atof());
	return res;
}	



bool isInAllowedArea( RooAbsPdf * pdf )
{
	double a = GetParam(pdf,"afb")->getVal();
	double f = GetParam(pdf,"fL")->getVal();

	return !((f-1)*3./4. > a || a >  -(f-1)*3./4.);
}


bool isInAllowedAreaB( RooAbsPdf * pdf )
{
	double a = GetParam(pdf,"afbB")->getVal();
	return (TMath::Abs(a) <= 0.5);
}

bool isInAllowedArea3D( RooAbsPdf * pdf )
{
	double a = GetParam(pdf,"afb")->getVal();
	double f = GetParam(pdf,"fL")->getVal();
	double aB = GetParam(pdf,"afbB")->getVal();

	return ( !((f-1)*3./4. > a || a >  -(f-1)*3./4.) && (TMath::Abs(aB) <= 0.5)); 
}



void getEfficiencies(double q2min, double q2max, RooAbsPdf ** effLLpdf, RooAbsPdf ** effDDpdf, RooAbsPdf ** effLLBpdf, RooAbsPdf ** effDDBpdf, bool printeff = false, double def = true)
{
	TString q2name = ((TString)Form("_q2_%4.2f_%4.2f",q2min,q2max)).ReplaceAll(".","");
	TString effbase = "/afs/cern.ch/work/p/pluca/results/";
	TH1F * effDD, * effLL, * effLLB, * effDDB;
	TString hname1 = "htoteff", hname2 = "htot_eff";
	TString sname = "";
	if(!def) { hname1 = "syseff", hname2 = "hsys_eff"; sname = "_sys"; }

	if(q2min == 8 && q2max == 11)
	{
		q2name = "_jpsi";
		TFile * effFile = TFile::Open(effbase+"Lbeff_JpsivscosThetaL_DD.root");
		effDD  = (TH1F *)effFile->Get(hname1);
		effFile = TFile::Open(effbase+"Lbeff_JpsivscosThetaL_LL.root");
		effLL  = (TH1F *)effFile->Get(hname1);
		effFile = TFile::Open(effbase+"Lbeff_JpsivscosThetaB_DD.root");
		effDDB  = (TH1F *)effFile->Get(hname1);
		effFile = TFile::Open(effbase+"Lbeff_JpsivscosThetaB_LL.root");
		effLLB  = (TH1F *)effFile->Get(hname1);
	}
	else if((q2min == 15 && q2max == 20) || (q2min == 1.1 && q2max == 6))
	{
		TFile * effFile = TFile::Open(effbase+"Lbeff2D_cosThetaL_vs_q2_DD_2bins.root");
		TH2F * effDD2D  = (TH2F *)effFile->Get(hname2);
		effDD = (TH1F*)GetSliceX(effDD2D,(q2max+q2min)/2.);
		effFile = TFile::Open(effbase+"Lbeff2D_cosThetaL_vs_q2_LL_2bins.root");
		TH2F * effLL2D  = (TH2F *)effFile->Get(hname2);
		effLL = (TH1F*)GetSliceX(effLL2D,(q2max+q2min)/2.);
		effFile = TFile::Open(effbase+"Lbeff2D_cosThetaB_vs_q2_DD_2bins.root");
		TH2F * effDDB2D  = (TH2F *)effFile->Get(hname2);
		effDDB = (TH1F*)GetSliceX(effDDB2D,(q2max+q2min)/2.);
		effFile = TFile::Open(effbase+"Lbeff2D_cosThetaB_vs_q2_LL_2bins.root");
		TH2F * effLLB2D  = (TH2F *)effFile->Get(hname2);
		effLLB = (TH1F*)GetSliceX(effLLB2D,(q2max+q2min)/2.);
	}
	else if(q2min == 0 && q2max == 20)
	{
		q2name = "_All";
		TFile * effFile = TFile::Open(effbase+"LbeffvscosThetaL_DD.root");
		effDD  = (TH1F *)effFile->Get(hname1);
		effFile = TFile::Open(effbase+"LbeffvscosThetaL_LL.root");
		effLL  = (TH1F *)effFile->Get(hname1);
		effFile = TFile::Open(effbase+"LbeffvscosThetaB_DD.root");
		effDDB  = (TH1F *)effFile->Get(hname1);
		effFile = TFile::Open(effbase+"LbeffvscosThetaB_LL.root");
		effLLB  = (TH1F *)effFile->Get(hname1);
	}
	else
	{
		TFile * effFile = TFile::Open(effbase+"Lbeff2D_cosThetaL_vs_q2_DD.root");
		TH2F * effDD2D  = (TH2F *)effFile->Get(hname2);
		effDD = (TH1F*)GetSliceX(effDD2D,(q2max+q2min)/2.);
		effFile = TFile::Open(effbase+"Lbeff2D_cosThetaL_vs_q2_LL.root");
		TH2F * effLL2D  = (TH2F *)effFile->Get(hname2);
		effLL = (TH1F*)GetSliceX(effLL2D,(q2max+q2min)/2.);
		effFile = TFile::Open(effbase+"Lbeff2D_cosThetaB_vs_q2_DD.root");
		TH2F * effDDB2D  = (TH2F *)effFile->Get(hname2);
		effDDB = (TH1F*)GetSliceX(effDDB2D,(q2max+q2min)/2.);
		effFile = TFile::Open(effbase+"Lbeff2D_cosThetaB_vs_q2_LL.root");
		TH2F * effLLB2D  = (TH2F *)effFile->Get(hname2);
		effLLB = (TH1F*)GetSliceX(effLLB2D,(q2max+q2min)/2.);
	}


	effLL->Scale(10./effLL->Integral());
	effDD->Scale(10./effDD->Integral());
	effLLB->Scale(10./effLLB->Integral());
	effDDB->Scale(10./effDDB->Integral());

	RooArgList * parList_effLLB = new RooArgList("parList_effLLB");
	RooRealVar * c0_effLLB = new RooRealVar("c0_effLLB","",-0.7,-1.,0.);
	RooRealVar * c1_effLLB = new RooRealVar("c1_effLLB","",0.,-0.2,0.1);
	RooArgList * parList_effDDB = new RooArgList("parList_effDDB");
	RooRealVar * c0_effDDB = new RooRealVar("c0_effDDB","",-0.2,-0.5,0.);
	RooRealVar * c1_effDDB = new RooRealVar("c1_effDDB","",-0.2,-0.5,0.);
	parList_effLLB->add(*c0_effLLB);
	parList_effLLB->add(*c1_effLLB);
	parList_effDDB->add(*c0_effDDB);
	parList_effDDB->add(*c1_effDDB);

	(*effLLBpdf) = new RooChebychev("effLLB","", *AngVars::cosThetaB, *parList_effLLB);
	(*effDDBpdf) = new RooChebychev("effDDB","", *AngVars::cosThetaB, *parList_effDDB);

	RooDataHist * hLLB = new RooDataHist("hLLB"+sname,"hLLB",*AngVars::cosThetaB,effLLB);
	RooDataHist * hDDB = new RooDataHist("hDDB"+sname,"hDDB",*AngVars::cosThetaB,effDDB);
	(*effLLBpdf)->fitTo(*hLLB,PrintLevel(-1),Minos(true),DataError(RooAbsData::SumW2),SumW2Error(true));
	(*effDDBpdf)->fitTo(*hDDB,PrintLevel(-1),Minos(true),DataError(RooAbsData::SumW2),SumW2Error(true));
	fixParam((*effLLBpdf),AngVars::cosThetaB);
	fixParam((*effDDBpdf),AngVars::cosThetaB);
	

	RooDataHist * hLL = new RooDataHist("hLL"+sname,"hLL",*AngVars::cosThetaL,effLL);
	RooDataHist * hDD = new RooDataHist("hDD"+sname,"hDD",*AngVars::cosThetaL,effDD);
	(*effLLpdf) = stringToPdf("Cheb2-v[0.,-0.5,0.5]","effLL"+sname, AngVars::cosThetaL);
	(*effDDpdf) = stringToPdf("Cheb2-v[0.,-0.5,0.5]","effDD"+sname, AngVars::cosThetaL);
	(*effLLpdf)->fitTo(*hLL,Minos(true),DataError(RooAbsData::SumW2),PrintLevel(-1));
	(*effDDpdf)->fitTo(*hDD,Minos(true),DataError(RooAbsData::SumW2),PrintLevel(-1));
	fixParam((*effLLpdf),AngVars::cosThetaL);
	fixParam((*effDDpdf),AngVars::cosThetaL);

	if(printeff)
	{
		TCanvas * ceff = new TCanvas();

		TPaveText * tbox = new TPaveText(gStyle->GetPadLeftMargin() + 0.05,
                    0.25 - gStyle->GetPadTopMargin(),
                    gStyle->GetPadLeftMargin() + 0.30,
                    0.40 - gStyle->GetPadTopMargin(),
                    "BRNDC");
		tbox->AddText("LHCb");
		tbox->AddText("simulation");
		tbox->SetFillStyle(0);
		tbox->SetTextAlign(12);
		tbox->SetBorderSize(0);

		RooPlot * pl = GetFrame(AngVars::cosThetaL, hLL,(*effLLpdf),"-noparams-nochi2",0,"cos #theta_{l}","Tot. eff. (A.U.)");
		pl->addObject(tbox);
		pl->Draw();
		ceff->Print("LLeffFit"+q2name+".pdf");
		ceff->Print("LLeffFit"+q2name+".png");
		ceff->Print("LLeffFit"+q2name+".eps");
		ceff->Print("LLeffFit"+q2name+".C");
		pl = GetFrame(AngVars::cosThetaL, hDD,(*effDDpdf),"-noparams-nochi2",0,"cos #theta_{l}","Tot. eff. (A.U.)");
		pl->addObject(tbox);
		pl->Draw();
		ceff->Print("DDeffFit"+q2name+".pdf");
		ceff->Print("DDeffFit"+q2name+".png");
		ceff->Print("DDeffFit"+q2name+".eps");
		ceff->Print("DDeffFit"+q2name+".C");
		pl = GetFrame(AngVars::cosThetaB, hDDB,(*effDDBpdf),"-noparams-nochi2",0,"cos #theta_{h}","Tot. eff. (A.U.)");
		pl->addObject(tbox);
		pl->Draw();
		ceff->Print("DDeffFitB"+q2name+".pdf");
		ceff->Print("DDeffFitB"+q2name+".png");
		ceff->Print("DDeffFitB"+q2name+".eps");
		ceff->Print("DDeffFitB"+q2name+".C");
		pl = GetFrame(AngVars::cosThetaB, hLLB,(*effLLBpdf),"-noparams-nochi2",0,"cos #theta_{h}","Tot. eff. (A.U.)");
		pl->addObject(tbox);
		pl->Draw();
		ceff->Print("LLeffFitB"+q2name+".pdf");
		ceff->Print("LLeffFitB"+q2name+".png");
		ceff->Print("LLeffFitB"+q2name+".eps");
		ceff->Print("LLeffFitB"+q2name+".C");
		delete ceff;
	}
}



void buildBkgPdfs(double q2min, double q2max, TString name, TCut cut, RooAbsPdf ** bkg, RooAbsPdf ** bkgB, TString model = "Poly2")
{
	TString q2name = ((TString)Form("_q2_%4.2f_%4.2f",q2min,q2max)).ReplaceAll(".","");
	
	string datafile = "/afs/cern.ch/work/p/pluca/Lmumu/weighted/candLb.root";
	string treename = "candLb2JpsiL";
	if(q2min == 8 && q2max == 11)
	{
		q2name = "_jpsi";
		//datafile = "/afs/cern.ch/work/p/pluca/Lmumu/weighted/Lb2Lmumu_CL_NBweighted.root";
		//treename = "tree";
	}


	Analysis * ana = new Analysis("data"+name+q2name,"Lb",treename,datafile);
	ana->AddVariable(AngVars::cosThetaL);
	ana->AddVariable(AngVars::cosThetaB);
	ana->AddVariable("J_psi_1S_MM");
	ana->AddVariable("Lb_MM");

	TCut sideBand = "Lb_MassConsLambda > 5690";
	//if(q2min == 8 && q2max == 11) sideBand = "Lb_MassConsLambda > 5690";
	TCut q2cut = (TCut)Form("TMath::Power(J_psi_1S_MM/1000,2) > %e && TMath::Power(J_psi_1S_MM/1000,2) < %e",q2min,q2max);
	TCut mycut = sideBand + q2cut + cut;

	ana->applyCuts(mycut);
	RooDataSet * data = ana->GetDataSet("-recalc");

	if(model=="RooKeysPdf")
	{
		(*bkg) = new RooKeysPdf("Background_"+name,"bkg"+name,*AngVars::cosThetaL,*data,RooKeysPdf::MirrorBoth,2);
		(*bkgB) = new RooKeysPdf("Background_"+name+"B","bkg"+name+"B",*AngVars::cosThetaB,*data,RooKeysPdf::MirrorBoth,2);
	}
	else
	{
		(*bkg) = stringToPdf(model,"background_"+name, AngVars::cosThetaL);
		(*bkg)->fitTo(*data,PrintLevel(-1));
		fixParam((*bkg),AngVars::cosThetaL);
		(*bkgB) = stringToPdf(model,"background_"+name+"B", AngVars::cosThetaB); 
		(*bkgB)->fitTo(*data,PrintLevel(-1));
		fixParam((*bkgB),AngVars::cosThetaB);
	}
}


Str2VarMap getPars(double q2min, double q2max, RooRealVar * vv, TString type, TCut cut, TFile * file, double * range = NULL)
{
	file->cd();
	TCut q2cut = cut;
	string optionsjpsi = "-ANDpulls-log-quiet-xM(#Lambda#mu#mu) [MeV/c^{2}]-noCost-minos-noParams";
	string optionsmumu = "-quiet-xM(#Lambda#mu#mu) (MeV/c^{2})-noCost-nochi2-noParams";
	if(q2min==8 && q2max==11) q2cut += "pplus_PIDp > 10";

	string datafilename = "/afs/cern.ch/work/p/pluca/Lmumu/weighted/candLb.root";
	string candfilename = "/afs/cern.ch/work/p/pluca/Lmumu/weighted/candLb_MC.root";
	//TFile * MCFile = TFile::Open(candfilename.c_str());
	//TTree * BdJpsiKSTree = (TTree *)MCFile->Get("candBdJpsiKS");
	TreeReader * readerMC = new TreeReader("candLb2JpsiL");
	readerMC->AddFile(candfilename.c_str());
	readerMC->Initialize();
	TreeReader * readerData = new TreeReader("candLb2JpsiL");
	readerData->AddFile(datafilename.c_str());
	readerData->Initialize();
	
	file->cd();
	RooRealVar * vJpsi = new RooRealVar("Lb_MassConsJpsiLambda","Lb_MassConsJpsiLambda",5620,5400.,5900.);
	Analysis * ana_MC = new Analysis("Lb2JpsiL_"+type+"_MC","Lb",readerMC,&q2cut,vJpsi);
	Analysis * anaLbJpsi_default = new Analysis("Jpsi_default_"+type,"Lb",readerData,&q2cut,vJpsi);
	file->cd();

	file->cd();
	string model = "DCB_OST";
	ana_MC->SetSignal((model+"-Xn[2]-Xn2[2]-s[9.,3.,20.]-s2[20.,9.,100.]").c_str());
	ana_MC->Initialize("-docuts");
	ana_MC->Fit(5400.,5900.,200,false,optionsmumu);
	Str2VarMap MCpars = ana_MC->GetSigParams();
	((RooRealVar*)MCpars["a"])->setConstant();
	((RooRealVar*)MCpars["a2"])->setConstant();
	((RooRealVar*)MCpars["f"])->setConstant();

	file->cd();
	//Analysis * KS = new Analysis("KS_bkg","Lb",BdJpsiKSTree,vJpsi,"DCB_OST","","-namepar");
	//KS->Fit(5300.,6000.,150,true,"-noPlot-quiet");
	//Str2VarMap KSpars = KS->GetSigParams();
	//setConstant(&KSpars);

	file->cd();
	anaLbJpsi_default->SetSignal(model.c_str(),8000.,"",MCpars);
	//anaLbJpsi_default->addBkgComponent("JpsiKS","DCB_OST",1.,"-namepar",KSpars);
	anaLbJpsi_default->addBkgComponent("Comb","Exp-b[-0.003,-0.01,0.]",1000.);
	anaLbJpsi_default->Initialize("-docuts-namepar");
	anaLbJpsi_default->Fit(5400.,5900.,200,false,optionsjpsi);

	if(range) anaLbJpsi_default->PrintComposition(range[0],range[1]);
		
	Str2VarMap jpsiSigpars = anaLbJpsi_default->GetSigParams();
	setConstant(&jpsiSigpars);

	RooRealVar * factor = new RooRealVar("factor","c",MCpars["s"]->getVal()/jpsiSigpars["s"]->getVal(),0.5,5);
	ModifyPars(&jpsiSigpars,"s",factor);
	ModifyPars(&jpsiSigpars,"s2",factor);

	return jpsiSigpars;
}

RooDataSet * getDataAndFrac(TString name, TString q2name, TreeReader * mydata, double * range, TCut cut, RooRealVar * MM,
		double * frac, Str2VarMap jpsiPars = Str2VarMap(), double *outnsig = NULL, bool is1D = true)
{
	TCut massCut = (TCut)Form((TString)MM->GetName() + " > %e && " + (TString)MM->GetName() + " < %e",range[0],range[1]);

	//cut.Print();
	Analysis * ana = new Analysis(name+"_mass"+q2name,"Lb",mydata,&cut,MM);
	ana->AddVariable("J_psi_1S_MM");
	ana->AddVariable(AngVars::cosThetaL);
	ana->AddVariable(AngVars::cosThetaB);
	ana->SetSignal("DCB_OST",1.e4,"",jpsiPars);
	ana->addBkgComponent("bkgM","Exp",1000.,"",jpsiPars);
	ana->Initialize("-docuts");
	ana->Fit(5400,5900,50,false,"-quiet-nocost-lin");
	frac[0] = ana->GetSigFraction(range[0],range[1],&frac[1]);

	TCut mycut = cut;
	if(is1D) mycut+=massCut;
	mycut.Print();
	ana->applyCuts(mycut);
	RooDataSet * datacut = ana->GetDataSet("-recalc-docuts");

	if(outnsig) *outnsig = datacut->sumEntries();
	return datacut;
}



RooDataSet * getSideData(TString name, TString q2name, TreeReader * mydata, double * range, TCut cut, RooRealVar * MM)
{
	TCut mycut = cut + (TCut)((TString)MM->GetName() + " > 5700");
	if(range) mycut = cut + (TCut)Form((TString)MM->GetName() + " > %e",range[0]) + cut;

	Analysis * ana = new Analysis(name+"_side_mass"+q2name,"Lb",mydata,&mycut,MM);
	ana->AddVariable("J_psi_1S_MM");
	ana->AddVariable(AngVars::cosThetaL);
	ana->AddVariable(AngVars::cosThetaB);
	mycut.Print();

	return ana->GetDataSet("-recalc-docuts");
}






#endif
