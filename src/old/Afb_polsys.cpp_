#include "ReadTree_comp.hpp"
#include "multi_analyser.hpp"
#include "Lb_cuts.hpp"
#include "RooAbsReal.h"
#include "RooSimultaneous.h"
#include "RooGenericPdf.h"
#include "TGraphErrors.h"
#include "TROOT.h"

using namespace RooFit;
using namespace std;


double maxWithSign(double null, double a, double b)
{
	double max = TMath::Max(TMath::Abs(a - null),TMath::Abs(null - b));
	if(max == TMath::Abs(a - null)) return (a - null);
	else return null - b;
}


int main(int argc, char **argv)
{
	RooRealVar * cosTheta = new RooRealVar("cosTheta","cosTheta",0,-1,1);
	RooRealVar * cosThetaL = new RooRealVar("cosThetaL","cosThetaL",0,-1,1);
	RooRealVar * cosThetaB = new RooRealVar("cosThetaB","cosThetaB",0,-1,1);


	TString effbase = "/afs/cern.ch/user/p/pluca/work/Lb/Lmumu/results/";
	TFile * effFile = TFile::Open(effbase+"Lbeff2D_cosTheta_vs_cosThetaL_All.root");
	TH1 * eff  = (TH1 *)effFile->Get("htot_eff");
	RooDataHist * heff = new RooDataHist("heff","heff",RooArgSet(*cosThetaL,*cosTheta),eff);
	RooHistPdf * effpdf = new RooHistPdf("effpdf","effpdf",RooArgSet(*cosThetaL,*cosTheta),*heff);
	TFile * effFileB = TFile::Open(effbase+"Lbeff2D_cosTheta_vs_cosThetaB_All.root");
	TH1 * effB  = (TH1 *)effFileB->Get("hupper_eff");
	RooDataHist * heffB = new RooDataHist("heffB","heffB",RooArgSet(*cosThetaB,*cosTheta),effB);
	RooHistPdf * effBpdf = new RooHistPdf("effBpdf","effBpdf",RooArgSet(*cosThetaB,*cosTheta),*heffB);

	RooRealVar * fL = new RooRealVar("fL","fL",0.82,0,1);
	RooRealVar * afb = new RooRealVar("afb","afb",0,-1,1);
	RooRealVar * afbB = new RooRealVar("afbB","afbB",0.37,-1,1);
	RooRealVar * alphaL = new RooRealVar("alphaL","alphaL",0.642);
	RooRealVar * Olp = new RooRealVar("Olp","Olp",-0.322);
	RooRealVar * Op = new RooRealVar("Op","Op",-0.461);
	RooRealVar * Ou12 = new RooRealVar("Ou12","Ou12",-0.302);
	RooRealVar * Pb = new RooRealVar("Pb","Pb",0.,-1,1);


	RooGenericPdf * distr =  new RooGenericPdf("distr","(3./8.)*(1 + cosThetaL^2)*(1-fL) + afb*cosThetaL + (3./4.)*fL*(1-cosThetaL^2) + Pb*cosTheta*(-(3./4.)*(1-cosThetaL^2)*Olp + (3./8.)*(1+cosThetaL^2)*Op - (3./8.)*cosThetaL*Ou12)",
			RooArgSet(*cosThetaL,*cosTheta,*fL,*afb,*Op,*Olp,*Ou12,*Pb) );

	RooGenericPdf * distrB =  new RooGenericPdf("distrB","1 + 2*afbB*cosThetaB + Pb*(Op-Olp)*cosTheta + alphaL*Pb*(1-2*fL)*cosTheta*cosThetaB",
			RooArgSet(*cosThetaB,*cosTheta,*afbB,*fL,*Op,*Olp,*Pb,*alphaL) );

	RooProdPdf * pdf = new RooProdPdf("pdf","pdf",*distr,*effpdf);
	RooProdPdf * pdfB = new RooProdPdf("pdfB","pdfB",*distrB,*effBpdf);



	effFile = TFile::Open(effbase+"LbeffvscosThetaL_All.root");
	TH1 * effAll  = (TH1F *)effFile->Get("htoteff");
	effFile = TFile::Open(effbase+"LbeffvscosThetaB_All.root");
	TH1 * effAllB  = (TH1F *)effFile->Get("htot_nodet_eff");

	RooDataHist * h = new RooDataHist("h","h",*cosThetaL,effAll);
	RooRealVar * c1 = new RooRealVar("c1","",0.,-100.,100);
	RooRealVar * c2 = new RooRealVar("c2","",0.,-100.,100);
	TString effstr = "(1 + c1*cosThetaL + c2*TMath::Power(cosThetaL,2))";
	RooAbsPdf * effpdfL = new RooGenericPdf("effpdfL", "", effstr, RooArgSet(*cosThetaL, *c1, *c2));
	effpdfL->fitTo(*h,PrintLevel(-1));
	fixParams(effpdfL,cosThetaL);


	RooDataHist * hB = new RooDataHist("hB","hB",*cosThetaB,effAllB);
	RooRealVar * cB1 = new RooRealVar("cB1","",0,-100.,100);
	RooRealVar * cB2 = new RooRealVar("cB2","",0,-100.,100);
	TString effBstr = "(1 + cB1*cosThetaB + cB2*TMath::Power(cosThetaB,2))";
	RooAbsPdf * effpdfB = new RooGenericPdf("effpdfB", "", effBstr, RooArgSet(*cosThetaB, *cB1, *cB2));
	effpdfB->fitTo(*hB,PrintLevel(-1));
	fixParams(effpdfB,cosThetaB);


	TString afbLpdf = "((3./8.)*(1.-fL)*(1 + TMath::Power(cosThetaL,2)) + afb*cosThetaL + (3./4.)*fL*(1 - TMath::Power(cosThetaL,2)))";
	RooAbsPdf * corrPdf = new RooGenericPdf("corrPdfLL",effstr+"*"+afbLpdf,RooArgSet(*cosThetaL, *afb, *fL, *c1, *c2) );

	TString afbBpdf = "(1 + 2*afbB*cosThetaB)";
	RooAbsPdf * corrPdfB = new RooGenericPdf("corrPdfB",effBstr+"*"+afbBpdf,RooArgSet(*cosThetaB, *afbB, *cB1, *cB2) );

	TH1F * fLsys = new TH1F("fLsys","",20,-1,1);
	TH1F * afbsys = new TH1F("afbsys","",20,-1,1);
	TH1F * afbBsys = new TH1F("afbBsys","",20,-1,1);

	RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
	TCanvas * c = new TCanvas();
	int nexp = 8000;
	int ngen = 200;

	for(int e = 0; e < nexp; e++)
	{
		showPercentage(e,nexp);

		Pb->setVal(-0.03);
		afb->setVal(0.);
		afbB->setVal(0.37);
		fL->setVal(0.82);

		RooDataSet * toy = pdf->generate(RooArgSet(*cosTheta,*cosThetaL),ngen);
		RooDataSet * toyB = pdfB->generate(RooArgSet(*cosTheta,*cosThetaB),ngen);
		/*
		   RooPlot * p = cosThetaL->frame();
		   toy->plotOn(p);
		   p->Draw();
		   c->Print("test.pdf");
		   RooPlot * pB = cosTheta->frame();
		   toyB->plotOn(pB);
		   pB->Draw();
		   c->Print("test2.pdf");
		   */
		afb->setVal(0.);
		afbB->setVal(0.);
		fL->setVal(0.7);

		corrPdf->fitTo(*toy,PrintLevel(-1));
		corrPdfB->fitTo(*toyB,PrintLevel(-1));

		double afbminus = afb->getVal();
		double afbBminus = afbB->getVal();
		double fLminus = fL->getVal();

		Pb->setVal(0.15);
		afb->setVal(0.);
		afbB->setVal(0.37);
		fL->setVal(0.82);

		RooDataSet * toy2 = pdf->generate(RooArgSet(*cosTheta,*cosThetaL),ngen);
		RooDataSet * toy2B = pdfB->generate(RooArgSet(*cosTheta,*cosThetaB),ngen);

		afb->setVal(0.);
		afbB->setVal(0.);
		fL->setVal(0.7);

		corrPdf->fitTo(*toy2,PrintLevel(-1));
		corrPdfB->fitTo(*toy2B,PrintLevel(-1));

		double afbplus = afb->getVal();
		double afbBplus = afbB->getVal();
		double fLplus = fL->getVal();
/*
		Pb->setVal(0.);
		afb->setVal(0.);
		afbB->setVal(0.37);
		fL->setVal(0.82);

		RooDataSet * toyc = pdf->generate(RooArgSet(*cosTheta,*cosThetaL),ngen);
		RooDataSet * toycB = pdfB->generate(RooArgSet(*cosTheta,*cosThetaB),ngen);

		afb->setVal(0.);
		afbB->setVal(0.);
		fL->setVal(0.7);

		corrPdf->fitTo(*toyc,PrintLevel(-1));
		corrPdfB->fitTo(*toycB,PrintLevel(-1));

		double afbnull = afb->getVal();
		double afbBnull = afbB->getVal();
		double fLnull = fL->getVal();

		fLsys->Fill(maxWithSign(fLnull,fLminus,fLplus));
		afbsys->Fill(maxWithSign(afbnull,afbminus,afbplus));
		afbBsys->Fill(maxWithSign(afbBnull,afbBminus,afbBplus));
*/
		fLsys->Fill(fLminus - fLplus);
		afbsys->Fill(afbminus - afbplus);
		afbBsys->Fill(afbBminus - afbBplus);
	}

	fLsys->GetXaxis()->SetTitle("#Delta f_{L}");
	fLsys->Draw();
	c->Print("fLsys_polarization.pdf");
	afbsys->GetXaxis()->SetTitle("#Delta A_{FB}");
	afbsys->Draw();
	c->Print("afbsys_polarization.pdf");
	afbBsys->GetXaxis()->SetTitle("#Delta A_{FB}^{h}");
	afbBsys->Draw();
	c->Print("afbBsys_polarization.pdf");

	cout << fixed << setprecision(7);
	cout << "fLsys = " << fLsys->GetMean() << " +/- " << fLsys->GetMeanError() << endl;
	cout << "afbsys = " << afbsys->GetMean() << " +/- " << afbsys->GetMeanError() << endl;
	cout << "afbBsys = " << afbBsys->GetMean() << " +/- " << afbBsys->GetMeanError() << endl;
}
