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
#include "RooMinuit.h"

using namespace RooFit;
using namespace std;

TString decayToDo = "Lb2Lmumu";//"Lb2JpsiL_reduced";

void getEfficiencies(double q2min, double q2max, RooAbsPdf ** effDDpdf, RooAbsPdf ** effLLpdf, RooAbsPdf ** effLLBpdf, RooAbsPdf ** effDDBpdf, bool printeff = false)
{
	RooRealVar * cosThetaL = new RooRealVar("cosThetaL","cosThetaL",0.,-1.,1.);
	RooRealVar * cosThetaB = new RooRealVar("cosThetaB","cosThetaB",0.,-1.,1.);

	TString q2name = ((TString)Form("q2_%4.2f_%4.2f",q2min,q2max)).ReplaceAll(".","");
	TString effbase = "/afs/cern.ch/user/p/pluca/work/Lb/Lmumu/results/";
	if(decayToDo == "Lb2JpsiL_reduced") effbase = "/afs/cern.ch/user/p/pluca/work/Lb/Lmumu/results/Jpsi_efficiencies/";
	TH1F * effDD, * effLL, * effLLB, * effDDB;

	TFile * effFile = TFile::Open(effbase+"LbeffvscosThetaL_DD.root");
	effDD  = (TH1F *)effFile->Get("htoteff");
	effFile = TFile::Open(effbase+"LbeffvscosThetaL_LL.root");
	effLL  = (TH1F *)effFile->Get("htoteff");
	effFile = TFile::Open(effbase+"LbeffvscosThetaB_DD.root");
	effDDB  = (TH1F *)effFile->Get("htot_nodet_eff");
	effFile = TFile::Open(effbase+"LbeffvscosThetaB_LL.root");
	effLLB  = (TH1F *)effFile->Get("htot_nodet_eff");

	RooDataHist * hLL = new RooDataHist("hLL","hLL",*cosThetaL,effLL);
	RooDataHist * hDD = new RooDataHist("hDD","hDD",*cosThetaL,effDD);
	RooRealVar * c1LL = new RooRealVar("c1LL","",0.,-1.,1);
	RooRealVar * c1DD = new RooRealVar("c1DD","",0.,-1.,1);
	RooRealVar * c2LL = new RooRealVar("c2LL","",0.,-1.,1);
	RooRealVar * c2DD = new RooRealVar("c2DD","",0.,-1.,1);
	TString effLLstr = "(1 + c1LL*cosThetaL + c2LL*TMath::Power(cosThetaL,2))";
	TString effDDstr = "(1 + c1DD*cosThetaL + c2DD*TMath::Power(cosThetaL,2))";
	(*effLLpdf) = new RooGenericPdf("effLLpdf", "", effLLstr, RooArgSet(*cosThetaL, *c1LL, *c2LL));
	(*effDDpdf) = new RooGenericPdf("effDDpdf", "", effDDstr, RooArgSet(*cosThetaL, *c1DD, *c2DD));
	(*effLLpdf)->fitTo(*hLL,PrintLevel(-1));
	(*effDDpdf)->fitTo(*hDD,PrintLevel(-1));
	fixParams((*effLLpdf),cosThetaL);
	fixParams((*effDDpdf),cosThetaL);	

	RooDataHist * hLLB = new RooDataHist("hLLB","hLLB",*cosThetaB,effLLB);
	RooDataHist * hDDB = new RooDataHist("hDDB","hDDB",*cosThetaB,effDDB);
	RooRealVar * cB1LL = new RooRealVar("cB1LL","",0,-1.,1);
	RooRealVar * cB1DD = new RooRealVar("cB1DD","",0,-1.,1);
	RooRealVar * cB2LL = new RooRealVar("cB2LL","",0,-1.,1);
	RooRealVar * cB2DD = new RooRealVar("cB2DD","",0,-1.,1);
	TString effLLBstr = "(1 + cB1LL*cosThetaB + cB2LL*TMath::Power(cosThetaB,2))";
	TString effDDBstr = "(1 + cB1DD*cosThetaB + cB2DD*TMath::Power(cosThetaB,2))";
	(*effLLBpdf) = new RooGenericPdf("effLLpdfB", "", effLLBstr, RooArgSet(*cosThetaB, *cB1LL, *cB2LL));
	(*effDDBpdf) = new RooGenericPdf("effDDpdfB", "", effDDBstr, RooArgSet(*cosThetaB, *cB1DD, *cB2DD));
	(*effLLBpdf)->fitTo(*hLLB,PrintLevel(-1));
	(*effDDBpdf)->fitTo(*hDDB,PrintLevel(-1));
	fixParams((*effLLBpdf),cosThetaB);
	fixParams((*effDDBpdf),cosThetaB);
}


void findMin(RooAbsReal * nll, RooRealVar * afb, RooRealVar *fL, double r)
{
	double afbInit = afb->getVal();
	double fLInit = fL->getVal();
	double step = 2.*r/10.;
	double minLogL = 1e6, minf = 1e6, mina = 1e6;

	for(double a = afbInit - r; a <= afbInit + r; a+=step )
	   for(double f = fLInit - r; f <= fLInit + r; f+=step )
	   {
		    if( (f-1)*3./4. > a || a >  -(f-1)*3./4. ) continue;
		   
			afb->setVal(a);
			fL->setVal(f);
 
			double tmp = nll->getVal();
            if(tmp < minLogL) { minLogL = tmp; mina = a; minf = f;}
		}

		afb->setVal(mina);
		fL->setVal(minf);
		afb->setRange(mina-r,mina+r);
		fL->setRange(minf-r,minf+r);
}



RooFitResult * safeFit(RooAbsPdf * pdf, RooDataSet * data, RooRealVar * afb, RooRealVar * fL)
{
	RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
	RooAbsReal * nll = pdf->createNLL(*data);
	RooMinuit m(*nll);
	m.setPrintLevel(-1);
	RooFitResult * res = m.fit("r");

	if(res->covQual()!=3 || res->edm() > 0.1)
	{
		double  minLogL = 1e9, mina = 1e6, minf = 1e6;
		for(double a = -0.75; a <= 0.75; a+=0.01 )
		{
			for(double f = 0.; f <= 1.; f+=0.01 )
			{
				if( (f-1)*3./4. > a || a >  -(f-1)*3./4. ) continue;
				else
				{
					afb->setVal(a);
					fL->setVal(f);

					double tmp = nll->getVal();
					if(tmp < minLogL) { minLogL = tmp; mina = a; minf = f;}
				}
			}
		}

		afb->setVal(mina);
		fL->setVal(minf);

		//cout << "***************************************************************************** -> " << mina << "      " << minf << endl;
		int loop = 0;
		do
		{
			if(loop==0) res = m.fit("r");
			else if(loop==1) 
			{ 
				afb->setVal(mina);
				fL->setVal(minf);
				afb->setRange(mina-0.1,mina+0.1);
				fL->setRange(minf-0.1,minf+0.1);
				res = m.fit("r");
			}
			else
			{
				double r = 0.08;
				for(int i = 1; i < 10; i*=2) findMin(nll,afb,fL,r/i);
				afb->setRange(afb->getVal()-0.01,afb->getVal()+0.01);
				fL->setRange(fL->getVal()-0.01,fL->getVal()+0.01);
				//res = m.fit("r");
				//if(res->covQual()!=3 || res->edm() > 0.01) cout << afb->getVal() << "   " << fL->getVal() << "      CRIBBIO!!!!!!!!" << endl;
			}
			
			loop++;
		}
		while((res->covQual()!=3 || res->edm() > 0.1) && loop < 3);
	}

	return res;
}



int main(int argc, char **argv)
{
	bool printeff = false;
	double genafb = 0, genfL = 0.8;
	bool bias = false;

	gROOT->ProcessLine(".x lhcbStyle.C");

	if(argc > 1)
	{
		for(int a = 1; a < argc; a++)
		{
			string arg = argv[a];
			string str = arg.substr(2,arg.length()-2);

			if(arg=="-peff") printeff = true;
			if(arg=="-bias") bias = true;
		}
	}


	RooRealVar * cosThetaL = new RooRealVar("cosThetaL","cosThetaL",0.,-1.,1.);
	//RooRealVar * cosThetaB = new RooRealVar("cosThetaB","cosThetaB",0.,-1.,1.);

	RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

	RooRealVar * afb = new RooRealVar("afb","afb",0.,-1,1);
	RooRealVar * fL = new RooRealVar("fL","fL",0.7,0.,1.);
	TString afbLpdf = "((3./8.)*(1.-fL)*(1 + TMath::Power(cosThetaL,2)) + afb*cosThetaL + (3./4.)*fL*(1 - TMath::Power(cosThetaL,2)))";
	//RooRealVar * afbB = new RooRealVar("afbB","afbB",0.,-100,100);
	//TString afbBpdf = "(1 + 2*afbB*cosThetaB)";
	RooAbsPdf * teoPdf = new RooGenericPdf("teoPdf",afbLpdf,RooArgSet(*cosThetaL,*afb,*fL));
	TCanvas * ceff = new TCanvas();

	/**               GET AND FIT EFFICIENCIES                  **/

	RooAbsPdf * effDDpdf = NULL, * effLLpdf = NULL, * effLLBpdf = NULL, * effDDBpdf = NULL;	
	getEfficiencies(15,20,&effDDpdf,&effLLpdf,&effLLBpdf,&effDDBpdf,printeff);
	

	/**                    FIT AFB                  **/

	TH2F * fLbias = new TH2F("fLbias","",30,-0.75,0.75,20,0,1);
	TH2F * afbbias = new TH2F("afbbias","",30,-0.75,0.75,20,0,1);
	TH2F * nfailed = new TH2F("nfailed","",30,-0.75,0.75,20,0,1);

	for(double a = -0.75; a <= 0.75; a+=0.05)
	{
		for(double f = 0; f <= 1; f+=0.05)
		{
			if((f-1)*3./4. > a || a >  -(f-1)*3./4.) continue;

			int nWrongErr = 0, nWrongEDM = 0, pass = 0;
			int ntoys = 100;
			TH1F * deltaafb = new TH1F("deltaafb","",20,-0.5,0.5);
			TH1F * deltafL = new TH1F("deltafL","",20,-0.5,0.5);

			genfL = f;
			genafb = a;

			for(int i = 0; i < ntoys; i++)
			{
				afb->setVal(genafb);
				fL->setVal(genfL);

				RooAbsPdf * corrPdfLL = new RooProdPdf("corrPdfLL","corrPdfLL",*teoPdf,*effLLpdf);

				RooDataSet * data = corrPdfLL->generate(*cosThetaL,200);

				afb->setVal(0);
				fL->setVal(0.7);

				RooFitResult * res = corrPdfLL->fitTo(*data,PrintLevel(-1), Save());

				if(!bias)
				{
					if(res->covQual()!=3) nWrongErr++;
					if(res->edm()>0.01) nWrongEDM++;
					if(res->covQual()==3 && res->edm()<0.01) pass++;
					deltafL->Fill(genfL-fL->getVal());
					deltaafb->Fill(genafb-afb->getVal());
				}
				else
				{
					if(res->covQual()!=3 || res->edm()>0.01) i--;
					else
					{
						deltafL->Fill(genfL-fL->getVal());
						deltaafb->Fill(genafb-afb->getVal());
					}
				}

			}

			cout << "fL = " << genfL << "   afb = " << genafb << endl;
			cout << fixed << setprecision(4);
			cout << "EDM>0.01 -> " << (double)nWrongEDM/ntoys << endl;
			cout << "Qual!=3  -> " << (double)nWrongErr/ntoys << endl;
			cout << "No fail  -> " << (double)pass/ntoys << endl;

			nfailed->Fill(a,f,(1-pass/ntoys));
			fLbias->Fill(a,f,deltafL->GetMean());
			afbbias->Fill(a,f,deltaafb->GetMean());
		}
	}

	nfailed->Draw("colz");
	ceff->Print("nfailed.pdf");
	fLbias->Draw("colz");
	ceff->Print("fL_bias.pdf");
	afbbias->Draw("colz");
	ceff->Print("afb_bias.pdf");

	delete ceff;
	
}
