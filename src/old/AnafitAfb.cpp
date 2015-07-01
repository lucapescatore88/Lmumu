#include "ReadTree_comp.hpp"
#include "multi_analyser.hpp"
#include "Lb_cuts.hpp"
#include "RooAbsReal.h"
#include "RooSimultaneous.h"
#include "RooGenericPdf.h"
#include "TGraphErrors.h"
#include "TROOT.h"
#include "FeldmanCousins.hpp"

using namespace RooFit;
using namespace std;

int main(int argc, char **argv)
{
	bool printSw = true;
	//TString massModel = "Gauss-m[5622]";
	string massModel = "DCB-m[5622]";
	TString effbase = "/afs/cern.ch/user/p/pluca/work/Lb/Lmumu/results/";
	bool printeff = false;
	TString dodata = "data";
	bool fitsingle = false;
	TString wstr = "physRate_polp006";
	TString decayToDo = "Lb2Lmumu";
	if(dodata=="genMC") wstr += "_noDecay";

	gROOT->ProcessLine(".x lhcbStyle.C");


	RooRealVar * cosThetaL = new RooRealVar("cosThetaL","cosThetaL",0.,-1.,1.);
	RooRealVar * cosThetaB = new RooRealVar("cosThetaB","cosThetaB",0.,-1.,1.);
	RooRealVar * nsig_sw = new RooRealVar("nsig_sw","nsig_sw",1,-1.e6,1.e6);
	RooRealVar * MCweight = new RooRealVar(wstr,wstr,1.,-1.e10,1.e10);
	RooRealVar * MM = new RooRealVar("Lb_MassConsLambda","Lb_MassConsLambda",5620.,5500.,5900.);
	TString datafilename = "/afs/cern.ch/user/p/pluca/work/Lb/Lmumu/candLb.root";
	if(dodata=="MC") datafilename = "/afs/cern.ch/user/p/pluca/work/Lb/Lmumu/candLb_MC.root";
	if(dodata=="genMC") datafilename = "/afs/cern.ch/work/p/pluca/weighted/Lmumu/"+(string)decayToDo+"_geomMC_Pythia8_NBweighted.root";
	TreeReader * data;
	if(dodata!="genMC") data = new TreeReader("cand"+decayToDo);
	else data = new TreeReader("MCtree");
	data->AddFile(datafilename);

	TFile * histFile = new TFile("Afb_hist.root","recreate");

	RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

	int nbins = 1;//CutsDef::nq2bins;
	double q2min[] = {15.,11.0,15,16,18};//&CutsDef::q2min_highfirst[0];
	double q2max[] = {20.,12.5,16,18,20};//&CutsDef::q2max_highfirst[0];
	
	//int nbins = CutsDef::nq2bins
	//double *q2min = &CutsDef::q2min[0];
	//double *q2max = &CutsDef::q2max[0];


	TGraphErrors * Afb_vs_q2 = new TGraphErrors();
	TGraphErrors * AfbB_vs_q2 = new TGraphErrors();
	TGraphErrors * fL_vs_q2 = new TGraphErrors();
	TCanvas * ceff = new TCanvas();

	RooCategory * samples = new RooCategory("samples","samples");
	samples->defineType("DD");
	samples->defineType("LL");

	RooRealVar * afb = new RooRealVar("afb","afb",0.,-100,100);
	RooRealVar * fL = new RooRealVar("fL","fL",0.7,-1.,10.);
	//RooRealVar * afb = new RooRealVar("afb","afb",0.,-1.,1.);
	//RooRealVar * fL = new RooRealVar("fL","fL",0.7,0.,1.);
	RooRealVar * origafb = new RooRealVar("afb","afb",0.,-1.,1.);
	RooRealVar * origfL = new RooRealVar("fL","fL",0.7,-1.,10.);
	TString afbLpdf = "((3./8.)*(1.-fL)*(1 + TMath::Power(cosThetaL,2)) + afb*cosThetaL + (3./4.)*fL*(1 - TMath::Power(cosThetaL,2)))";
	RooRealVar * afbB = new RooRealVar("afbB","afbB",0.,-100,100);
	//RooRealVar * afbB = new RooRealVar("afbB","afbB",0.,-1.,1.);
	RooRealVar * origafbB = new RooRealVar("afbB","afbB",0.,-1.,1.);
	TString afbBpdf = "(1 + 2*afbB*cosThetaB)";

	vector< vector< double > > afb_errs, afbB_errs, fL_errs;
	TList * LLlist = new TList, * DDlist = new TList;

	TCanvas * cDD = new TCanvas();
	TCanvas * cLL = new TCanvas();
	TCanvas * cDDB = new TCanvas();
	TCanvas * cLLB = new TCanvas();

	for(int i = 0; i < nbins; i++)
	{
		//if(q2min[i] < 8) continue;
		TString q2name = ((TString)Form("q2_%4.2f_%4.2f",q2min[i],q2max[i])).ReplaceAll(".","");
		TString curq2cut = Form("TMath::Power(J_psi_1S_MM/1000,2) >= %e && TMath::Power(J_psi_1S_MM/1000,2) < %e",q2min[i],q2max[i]);	
		//TString curq2cut = Form("TMath::Power(J_psi_1S_MM/1000,2) >= %e && TMath::Power(J_psi_1S_MM/1000,2) < %e && (Lb_MassConsLambda > 5680 || Lb_MassConsLambda < 5590)",q2min[i],q2max[i]); 
		cout << "------------------- q2 bin: " << q2min[i] << " - " << q2max[i] << " -----------------------" << endl;

		TFile * effFile = NULL;
		TH1F * effDD = NULL, * effLL = NULL, * effLLB = NULL, * effDDB = NULL;
		if(q2min[i] == 15 && q2max[i] == 20)
		{
			effFile = TFile::Open(effbase+"LbeffvscosThetaL_DD.root");
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

			effFile = TFile::Open(effbase+"Lbeff2D_cosThetaL_vs_q2_DD.root");
			TH2F * effDD2D  = (TH2F *)effFile->Get("htot_eff");
			effDD = (TH1F*)GetSliceX(effDD2D,(q2max[i]+q2min[i])/2.);
			effFile = TFile::Open(effbase+"Lbeff2D_cosThetaL_vs_q2_LL.root");
			TH2F * effLL2D  = (TH2F *)effFile->Get("htot_eff");
			effLL = (TH1F*)GetSliceX(effLL2D,(q2max[i]+q2min[i])/2.);
			effFile = TFile::Open(effbase+"Lbeff2D_cosThetaB_vs_q2_DD.root");
			TH2F * effDDB2D  = (TH2F *)effFile->Get("hupper_eff");
			effDDB = (TH1F*)GetSliceX(effDDB2D,(q2max[i]+q2min[i])/2.);
			effFile = TFile::Open(effbase+"Lbeff2D_cosThetaB_vs_q2_LL.root");
			TH2F * effLLB2D  = (TH2F *)effFile->Get("hupper_eff");
			effLLB = (TH1F*)GetSliceX(effLLB2D,(q2max[i]+q2min[i])/2.);
		}

		ceff->cd();

		/**                    FIT EFFICIENCY                  **/

		RooDataHist * hLL = new RooDataHist("hLL","hLL",*cosThetaL,effLL);
		RooDataHist * hDD = new RooDataHist("hDD","hDD",*cosThetaL,effDD);
		RooRealVar * c1LL = new RooRealVar("c1LL","",0.,-1.,1);
		RooRealVar * c1DD = new RooRealVar("c1DD","",0.,-1.,1);
		RooRealVar * c2LL = new RooRealVar("c2LL","",0.,-1.,1);
		RooRealVar * c2DD = new RooRealVar("c2DD","",0.,-1.,1);
		TString effLLstr = "(1 + c1LL*cosThetaL + c2LL*TMath::Power(cosThetaL,2))";
		TString effDDstr = "(1 + c1DD*cosThetaL + c2DD*TMath::Power(cosThetaL,2))";
		RooAbsPdf * effLLpdf = new RooGenericPdf("effLLpdf", "", effLLstr, RooArgSet(*cosThetaL, *c1LL, *c2LL));
		RooAbsPdf * effDDpdf = new RooGenericPdf("effDDpdf", "", effDDstr, RooArgSet(*cosThetaL, *c1DD, *c2DD));
		effLLpdf->fitTo(*hLL,PrintLevel(-1));
		effDDpdf->fitTo(*hDD,PrintLevel(-1));
		fixParams(effLLpdf,cosThetaL);
		fixParams(effDDpdf,cosThetaL);	

		RooDataHist * hLLB = new RooDataHist("hLLB","hLLB",*cosThetaB,effLLB);
		RooDataHist * hDDB = new RooDataHist("hDDB","hDDB",*cosThetaB,effDDB);
		RooRealVar * cB1LL = new RooRealVar("cB1LL","",0,-1.,1);
		RooRealVar * cB1DD = new RooRealVar("cB1DD","",0,-1.,1);
		RooRealVar * cB2LL = new RooRealVar("cB2LL","",0,-1.,1);
		RooRealVar * cB2DD = new RooRealVar("cB2DD","",0,-1.,1);
		TString effLLBstr = "(1 + cB1LL*cosThetaB + cB2LL*TMath::Power(cosThetaB,2))";
		TString effDDBstr = "(1 + cB1DD*cosThetaB + cB2DD*TMath::Power(cosThetaB,2))";
		RooAbsPdf * effLLpdfB = new RooGenericPdf("effLLpdfB", "", effLLBstr, RooArgSet(*cosThetaB, *cB1LL, *cB2LL));
		RooAbsPdf * effDDpdfB = new RooGenericPdf("effDDpdfB", "", effDDBstr, RooArgSet(*cosThetaB, *cB1DD, *cB2DD));
		effLLpdfB->fitTo(*hLLB,PrintLevel(-1));
		effDDpdfB->fitTo(*hDDB,PrintLevel(-1));
		fixParams(effLLpdfB,cosThetaB);
		fixParams(effDDpdfB,cosThetaB);

		//cout << q2min[i] << " - " << q2max[i] << " LL cosThetaL -> " << c1LL->getVal() << "  " << c2LL->getVal() << endl;
		//cout << q2min[i] << " - " << q2max[i] << " DD cosThetaL -> " << c1DD->getVal() << "  " << c2DD->getVal() << endl;
		//cout << q2min[i] << " - " << q2max[i] << " LL cosThetaB -> " << cB1LL->getVal() << "  " << cB2LL->getVal() << endl;
		//cout << q2min[i] << " - " << q2max[i] << " DD cosThetaB -> " << cB1DD->getVal() << "  " << cB2DD->getVal() << endl;


		if(printeff) {
			GetFrame(cosThetaL, hLL,effLLpdf,"-nochi2",0,NULL,0,"cos#theta_{l}","Tot. eff.")->Draw();
			ceff->Print("DDeffFit"+q2name+".pdf");
			GetFrame(cosThetaL, hDD,effDDpdf,"-nochi2",0,NULL,0,"cos#theta_{l}","Tot. eff.")->Draw();
			ceff->Print("LLeffFit"+q2name+".pdf");
			GetFrame(cosThetaB, hLLB,effLLpdfB,"-nochi2",0,NULL,0,"cos#theta_{#Lambda}","Tot. eff.")->Draw();
			ceff->Print("DDeffFitB"+q2name+".pdf");
			GetFrame(cosThetaB, hDDB,effDDpdfB,"-nochi2",0,NULL,0,"cos#theta_{#Lambda}","Tot. eff.")->Draw();
			ceff->Print("LLeffFitB"+q2name+".pdf"); }

			/**                    FIT AFB                  **/


			afb->setVal(0);
			afbB->setVal(0);
			fL->setVal(0.7);

			TString LLnorm = "1./( 1. + (2./3.)*afb*c1LL + (2./5.)*c2LL - (1./5.)*c2LL*fL )*"+effLLstr;
			TString DDnorm = "1./( 1. + (2./3.)*afb*c1DD + (2./5.)*c2DD - (1./5.)*c2DD*fL )*"+effDDstr;
			RooAbsPdf * corrPdfLL = new RooGenericPdf(Form("corrPdfLL_%i",i),LLnorm+"*"+afbLpdf,RooArgSet(*cosThetaL, *afb, *fL, *c1LL, *c2LL) );
			RooAbsPdf * corrPdfDD = new RooGenericPdf(Form("corrPdfDD_%i",i),DDnorm+"*"+afbLpdf,RooArgSet(*cosThetaL, *afb, *fL, *c1DD, *c2DD) );

			TString LLnormB = "1./( (2./3.)*( 2*afbB*cB1LL + cB2LL + 3.) )*"+effLLBstr;
			TString DDnormB = "1./( (2./3.)*( 2*afbB*cB1DD + cB2DD + 3.) )*"+effDDBstr;
			RooAbsPdf * corrPdfLLB = new RooGenericPdf(Form("corrPdfLLB_%i",i),LLnormB+"*"+afbBpdf,RooArgSet(*cosThetaB, *afbB, *cB1LL, *cB2LL) );
			RooAbsPdf * corrPdfDDB = new RooGenericPdf(Form("corrPdfDDB_%i",i),DDnormB+"*"+afbBpdf,RooArgSet(*cosThetaB, *afbB, *cB1DD, *cB2DD) );

			TCut cutLL = CutsDef::LLcut + (TCut)curq2cut;
			TCut cutDD = CutsDef::DDcut + (TCut)curq2cut;

			if(dodata=="genMC")
			{
				corrPdfLLB = new RooGenericPdf("corrPdfLL",afbBpdf,RooArgSet(*cosThetaB, *afbB, *cB1LL, *cB2LL) );
				corrPdfDDB = new RooGenericPdf("corrPdfDD",afbBpdf,RooArgSet(*cosThetaB, *afbB, *cB1DD, *cB2DD) );
				corrPdfLL = new RooGenericPdf("corrPdfLL",afbLpdf,RooArgSet(*cosThetaL, *afb, *fL, *c1LL, *c2LL) );
				corrPdfDD = new RooGenericPdf("corrPdfDD",afbLpdf,RooArgSet(*cosThetaL, *afb, *fL, *c1DD, *c2DD) );
				cutLL = (TCut)curq2cut;
				cutDD = (TCut)curq2cut;
			}

			Analysis * anaLL = new Analysis(Form("LL_mass_%i",i),"Lb",data,&cutLL,MM);
			anaLL->AddVariable(cosThetaL);
			anaLL->AddVariable(cosThetaB);
			anaLL->AddVariable("J_psi_1S_MM");
			if(dodata!="data") anaLL->SetWeight(wstr);
			RooDataSet * dataLL = anaLL->GetDataSet("-recalc-docuts");

			Analysis * anaDD = new Analysis(Form("DD_mass_%i",i),"Lb",data,&cutDD,MM);
			anaDD->AddVariable(cosThetaL);
			anaDD->AddVariable(cosThetaB);
			anaDD->AddVariable("J_psi_1S_MM");
			if(dodata!="data") anaDD->SetWeight(wstr);
			RooDataSet * dataDD = anaDD->GetDataSet("-recalc-docuts");

			RooDataSet * sdataDD, * sdataLL;

			if(dodata=="data")
			{
				sdataLL = anaLL->CalcSweight("",massModel.c_str(),"Exp");

				if(printSw) {
					GetFrame(MM,NULL,sdataLL,"-nochi2",30,NULL,0,"M(#Lambda#mu#mu) (MeV/c^{2})")->Draw();
					ceff->Print("Mass_LL_sWeighted"+q2name+".pdf");
					GetFrame(cosThetaL,NULL,sdataLL,"-nochi2",6,NULL,0,"cos#theta_{l}")->Draw();
					ceff->Print("cosThetaL_LL_sWeighted"+q2name+".pdf");
					GetFrame(cosThetaL,NULL,dataLL,"-nochi2",6,NULL,0,"cos#theta_{l}")->Draw();
					ceff->Print("cosThetaL_LL_"+q2name+".pdf");
				}

				sdataDD = anaDD->CalcSweight("",massModel.c_str(),"Exp");

				if(printSw) {
					GetFrame(MM,NULL,sdataDD,"-nochi2",30,NULL,0,"M(#Lambda#mu#mu) (MeV/c^{2})")->Draw();
					ceff->Print("Mass_DD_sWeighted"+q2name+".pdf");
					GetFrame(cosThetaL,NULL,sdataDD,"-nochi2",10,NULL,0,"cos#theta_{l}")->Draw();
					ceff->Print("cosThetaL_DD_sWeighted"+q2name+".pdf");
					GetFrame(cosThetaL,NULL,dataDD,"-nochi2",10,NULL,0,"cos#theta_{l}")->Draw();
					ceff->Print("cosThetaL_DD_"+q2name+".pdf");
				}
			}		
			else { sdataLL = dataLL; sdataDD = dataDD; }

			histFile->cd();
			TTree * LLTree = (TTree*)sdataLL->tree();
			LLTree->SetName(Form("treeLL_%i",i));
			LLlist->Add(LLTree);
			TTree * DDTree = (TTree*)sdataDD->tree();
			DDTree->SetName(Form("treeDD_%i",i));
			DDlist->Add(DDTree);


			// CREATE COMBINED DATASET
			RooDataSet * combData;
			if(dodata=="data") combData = new RooDataSet(Form("combData_%i",i),"combined data",RooArgSet(*cosThetaL,*cosThetaB,*nsig_sw),Index(*samples),Import("DD",*sdataDD),Import("LL",*sdataLL),WeightVar("nsig_sw"));
			else combData = new RooDataSet(Form("combData_%i",i),"combined data",RooArgSet(*cosThetaL,*cosThetaB,*MCweight),Index(*samples),Import("DD",*sdataDD),Import("LL",*sdataLL),WeightVar(wstr));


			// FIT COS LEPTON
			RooSimultaneous * combModel = new RooSimultaneous(Form("combModel_%i",i),"",*samples);
			combModel->addPdf(*corrPdfLL,"LL");
			combModel->addPdf(*corrPdfDD,"DD");

			combModel->fitTo(*combData,PrintLevel(-1),Verbose(kFALSE),SumW2Error(kTRUE));

			if(fitsingle) corrPdfLL->fitTo(*sdataLL,PrintLevel(-1),Verbose(kFALSE),SumW2Error(kTRUE));
			GetFrame(cosThetaL,corrPdfLL,sdataLL,"-sumW2err-nochi2-noCost",6,NULL,0,"cos#theta_{l}")->Draw();
			ceff->Print("Afb_LL_"+q2name+".pdf");
			if(fitsingle) corrPdfDD->fitTo(*sdataDD,PrintLevel(-1),Verbose(kFALSE),SumW2Error(kTRUE));		
			GetFrame(cosThetaL,corrPdfDD,sdataDD,"-sumW2err-nochi2-noCost",10,NULL,0,"cos#theta_{l}")->Draw();
			ceff->Print("Afb_DD_"+q2name+".pdf");

			Afb_vs_q2->SetPoint(i,(q2max[i] + q2min[i])/2.,afb->getVal());
			Afb_vs_q2->SetPointError(i,(q2max[i] - q2min[i])/2.,afb->getError());
			fL_vs_q2->SetPoint(i,(q2max[i] + q2min[i])/2.,fL->getVal());
			fL_vs_q2->SetPointError(i,(q2max[i] - q2min[i])/2.,fL->getError());
				
			// FIT COS HADRON
			RooSimultaneous * combModelB = new RooSimultaneous(Form("combModelB_%i",i),"",*samples);
			combModelB->addPdf(*corrPdfLLB,"LL");
			combModelB->addPdf(*corrPdfDDB,"DD");

			combModelB->fitTo(*combData,PrintLevel(-1),Verbose(kFALSE),SumW2Error(kTRUE));

			if(fitsingle) corrPdfLLB->fitTo(*sdataLL,PrintLevel(-1),Verbose(kFALSE),SumW2Error(kTRUE));
			GetFrame(cosThetaB,corrPdfLLB,sdataLL,"-sumW2err-nochi2-noCost",6,NULL,0,"cos#theta_{#Lambda}")->Draw();
			ceff->Print("AfbB_LL_"+q2name+".pdf");
			if(fitsingle) corrPdfDDB->fitTo(*sdataDD,PrintLevel(-1),Verbose(kFALSE),SumW2Error(kTRUE));		
			GetFrame(cosThetaB,corrPdfDDB,sdataDD,"-sumW2err-nochi2-noCost",10,NULL,0,"cos#theta_{#Lambda}")->Draw();
			ceff->Print("AfbB_DD_"+q2name+".pdf");

			AfbB_vs_q2->SetPoint(i,(q2max[i] + q2min[i])/2.,afbB->getVal());
			AfbB_vs_q2->SetPointError(i,(q2max[i] - q2min[i])/2.,afbB->getError());
			
			cout << endl << fixed << setprecision(6) << "AfbB = " << afbB->getVal() << " +/- " << afbB->getError() << endl;
			cout << "Afb = " << afb->getVal() << " +/- " << afb->getError() << endl;
			cout << "fL = " << fL->getVal() << " +/- " << fL->getError() << endl;
			cout << endl;
			cout << "------------------------ FELDMAN AND COUSINS ------------------------" << endl;

			vector < RooDataSet * > datas;
			vector < RooAbsPdf * > pdfs, pdfsB;
			vector < TString > cat;
			cat.push_back("LL");
			cat.push_back("DD");
			datas.push_back(sdataLL);
			datas.push_back(sdataDD);

			RooArgSet * origPars = new RooArgSet();
			origPars->add(*origafb);
			origPars->add(*origfL);
			
			pdfs.push_back(corrPdfLL);
			pdfs.push_back(corrPdfDD);

			vector< double > afb_err, afbB_err, fL_err;
/*
			double fLval = fL->getVal(), fLerr = fL->getError();
			FeldmanCousins * FC = new FeldmanCousins(q2name,cat,datas,pdfs,cosThetaL,afb,"nsig_sw");
			//FC->SetNPointsToScan(20);
			//FC->SetNExp(1000);
			if(q2min[i]==18) afb_err = FC->ExtractLimits(origPars,-0.3,0.3);
			else if( (afb->getVal()-1.4*afb->getError()) > -1 && (afb->getVal()+1.4*afb->getError()) < 1 )
		       afb_err = FC->ExtractLimits(origPars,afb->getVal()-1.4*afb->getError(),afb->getVal()+1.4*afb->getError());
		    else afb_err = FC->ExtractLimits(origPars,-0.4,0.4);

			//FeldmanCousins * FCfL = new FeldmanCousins(q2name,cat,datas,pdfs,cosThetaL,fL,"nsig_sw");
			//if(q2min[i]==11) fL_err = FCfL->ExtractLimits(origPars,0.,0.6);
			//else if (q2min[i]==18) fL_err = FCfL->ExtractLimits(origPars,0.75,0.992);
			//( (fLval-1.3*fLerr) > 0 && (fLval+1.3*fLerr) <= 1 )
			//else fL_err = FCfL->ExtractLimits(origPars,fLval-1.3*fLerr,fLval+1.3*fLerr);

			afb_errs.push_back(afb_err);
			//fL_errs.push_back(fL_err);

 			RooArgSet * origParsB = new RooArgSet();
			origParsB->add(*origafbB);
			pdfsB.push_back(corrPdfLLB);
			pdfsB.push_back(corrPdfDDB);

			FeldmanCousins * FCB = new FeldmanCousins(q2name,cat,datas,pdfsB,cosThetaB,afbB,"nsig_sw");
			if( (afbB->getVal()-1.5*afbB->getError()) > -1 && (afbB->getVal()+1.5*afbB->getError()) < 1 )
			   afbB_err = FCB->ExtractLimits(origParsB,afbB->getVal()-1.5*afbB->getError(),afbB->getVal()+1.5*afbB->getError());
			else afbB_err = FCB->ExtractLimits(origParsB,-0.4,0.4);

			afbB_errs.push_back(afbB_err);
*/
			delete effDD;
			delete effLL;
			delete effLLB;
			delete effDDB;
	}

	cDD->Print("DDeff.pdf");
	cLL->Print("LLeff.pdf");
	cDDB->Print("DDBeff.pdf");
	cLLB->Print("LLBeff.pdf");


	Afb_vs_q2->GetXaxis()->SetTitle("q^{2}");
	Afb_vs_q2->GetYaxis()->SetTitle("Afb");
	Afb_vs_q2->SetMaximum(1);
	Afb_vs_q2->SetMinimum(-1);
	Afb_vs_q2->Draw("AP");
	ceff->Print("Afb_vs_q2.pdf");
	AfbB_vs_q2->GetXaxis()->SetTitle("q^{2}");
	AfbB_vs_q2->GetYaxis()->SetTitle("AfbB");
	AfbB_vs_q2->SetMaximum(1);
	AfbB_vs_q2->SetMinimum(-1);
	AfbB_vs_q2->Draw("AP");
	ceff->Print("AfbB_vs_q2.pdf");
	fL_vs_q2->GetXaxis()->SetTitle("q^{2}");
	fL_vs_q2->GetYaxis()->SetTitle("fL");
	fL_vs_q2->Draw("AP");
	ceff->Print("fL_vs_q2.pdf");

	for(int bb = 0; bb < Afb_vs_q2->GetN(); bb++)
	{
		double qq, qqerr, afbv, afbBv, fLv;
		Afb_vs_q2->GetPoint(bb,qq,afbv);
		qqerr = Afb_vs_q2->GetErrorX(bb);
		AfbB_vs_q2->GetPoint(bb,qq,afbBv);
		fL_vs_q2->GetPoint(bb,qq,fLv);
		cout << fixed << setprecision(1) << qq-qqerr << " - " << qq+qqerr;
		cout << fixed << setprecision(4); 
		//cout << " & $" << afbv << "_{-" << TMath::Abs(afb_errs[bb][0] - afbv) << "}^{+" << TMath::Abs(afb_errs[bb][1] - afbv)  << "} \\text{(stat)} \\pm \\text{(sys)}$ ";
		//cout << " & $" << afbBv << "_{-" << TMath::Abs(afbB_errs[bb][0] - afbBv) << "}^{+" << TMath::Abs(afbB_errs[bb][1]-afbBv) << "} \\text{(stat)} \\pm \\text{(sys)}$ " ;
		//cout << " & $" << fLv << "_{-" << TMath::Abs(fL_errs[bb][0] - fLv) << "}^{+" << TMath::Abs(fL_errs[bb][1] - fLv)  << "} \\text{(stat)} \\pm \\text{(sys)}$ ";
		cout << "  \\\\ " << endl;
	}

	histFile->cd();
	TTree * finalLLtree = (TTree*)TTree::MergeTrees(LLlist);
	TTree * finalDDtree = (TTree*)TTree::MergeTrees(DDlist);
	finalLLtree->SetName("LL_data");
	finalDDtree->SetName("DD_data");
	finalLLtree->Write();
	finalDDtree->Write();

	delete ceff;
	histFile->Write();
	delete histFile;

	}
