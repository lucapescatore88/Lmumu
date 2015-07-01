#include "analyser.hpp"
#include "TCanvas.h"
#include "Lb_cuts.hpp"
#include "TGraphErrors.h"
#include "RooGenericPdf.h"
#include "RooMsgService.h"
#include "RooMinuit.h"
#include "functions.hpp"

TH1F * Smear( TH1F * h )
{
	TH1F * res = (TH1F*)h->Clone();
	TRandom3 r(0);
	for(int i = 1; i <= h->GetNbinsX(); i++)
	{
		double smv = r.Gaus(h->GetBinContent(i),h->GetBinError(i));
		res->SetBinContent(i,smv);
		res->SetBinError(i,h->GetBinError(i));
	}

	return res;
}



int main(int argc, char **argv)
{
	TString var = "cosThetaL";
	TString type = "All";
	string systype = "";
	int gradbkg = 0;

	if(argc > 1)
	{
		for(int a = 1; a < argc; a++)
		{
			string arg = argv[a];
			string str = arg.substr(2,arg.length()-2);

			if(arg.find("-t") != string::npos) type = (TString)str;
			if(arg.find("-v") != string::npos) var = (TString)str;
			if(arg == "flucteff" || arg == "effsymm") systype = arg;
		}
	}


	TCut baseCut = "";
	if(type == "DD") { baseCut = CutsDef::DDcut; }
	else if(type == "LL") { baseCut = CutsDef::LLcut; }

	TString q2str = "TMath::Power(J_psi_1S_MM/1000,2)";
	int q2nbins = 5;
	double q2min[] = {15.0, 11.0, 15.0, 16.0, 18.0};
	double q2max[] = {20.0, 12.5, 16.0, 18.0, 20.0};
	double nevts[] = {200., 34., 41., 88., 72. };

	if(systype != "flucteff") q2nbins = 1; 
	vector <double> sys, fLsys;

	RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

	for(int q2 = 1; q2 < q2nbins; q2++)
	{
		TString q2name = ((TString)Form("q2_%4.2f_%4.2f",q2min[q2],q2max[q2])).ReplaceAll(".","");	

		RooRealVar * cosThetaL = new RooRealVar("cosThetaL","cosThetaL",0.,-1.,1.);
		RooRealVar * cosThetaB = new RooRealVar("cosThetaB","cosThetaB",0.,-1.,1.);

		RooAbsPdf * effDD = NULL, * effLL = NULL, * effLLB = NULL, * effDDB = NULL;	
		getEfficiencies(q2min[q2],q2max[q2],&effLL,&effDD,&effLLB,&effDDB,false);

		getEfficiencies(q2min[q2],q2max[q2],&effLL,&effDD,&effLLB,&effDDB,false);

		Analysis * normal; 
		RooRealVar * afb = new RooRealVar("afb","afb",0.04,-0.75,0.75);
		RooRealVar * fL = new RooRealVar("fL","fL",0.82,0.,1.);
		TString afbpdf = "((3./8.)*(1.-fL)*(1 + TMath::Power(cosThetaL,2)) + afb*cosThetaL + (3./4.)*fL*(1 - TMath::Power(cosThetaL,2)))";
		RooRealVar * afbB = new RooRealVar("afbB","afbB",-0.37,-1.,1.);
		TString afbBpdf = "(1 + 2*afbB*cosThetaB)";
		Analysis::SetPrintLevel("s");

		TH1F * hsys = new TH1F("hsys"+q2name,"",20,-0.5,0.5); 	
		TH1F * hfLsys = new TH1F("hfLsys"+q2name,"",20,-0.5,0.5);

		//double c1v = c1->getVal(), c2v = c2->getVal();
		//double c1_err = c1->getError(), c2_err = c2->getError();
		//double c1Bv = cB1->getVal(), c2Bv = cB2->getVal();
		//double c1B_err = cB1->getError(), c2B_err = cB2->getError();
		TRandom3 r(0);

		//c1->Print();
		//c2->Print();
		//cB1->Print();
		//cB2->Print();

		//cout <<  c1->getVal() << "  " <<  c2->getVal() << endl;
		//cout <<  c1s->getVal() << "  " <<  c2s->getVal() << endl;
		//cout <<  cB1->getVal() << "  " <<  cB2->getVal() << endl;
		//cout <<  cB1s->getVal() << "  " <<  cB2s->getVal() << endl;

		Str2VarMap params;
		params["fL"] = fL;
		params["afb"] = afb;	
		Str2VarMap paramsB;
		paramsB["afbB"] = afbB;


		for(int e = 0; e < 1000; e++)
		{
			showPercentage(e,1000);
			RooAbsPdf * corrPdf, * pdf, * corrPdfs, * bkgDD_teo, * bkgDD, * bkgDDs, * modelDD, * modelDDs;
			
			//RooRealVar * fracLL = new RooRealVar("fracLL","fracLL",fracLLv);
			RooRealVar * fracDD = new RooRealVar("fracDD","fracDD",fracDDv);

			string bkgmodel = Form("Poly%i-v0-l2",gradbkg);
			
			if(var=="cosThetaL")
			{
				pdf = new RooGenericPdf(Form("pdf_%i_",e),afbpdf,RooArgSet(*cosThetaL, *afb, *fL) );
				bkgDD_teo = stringToPdf(bkgmodel.c_str(),"bkgDD_teo", cosThetaL);

				corrPdf = new RooProdPdf(Form("corrPdf_%i",e),*pdf,*effDD);
				corrPdfs = new RooProdPdf(Form("corrPdfs_%i",e),*pdf,*effDDs);
			
				bkgDD = new RooProdPdf("bkgDD","",*bkgDD_teo,*effDD);
				bkgDDs = new RooProdPdf("bkgDDs","",*bkgDD_teo,*effDDs);

				modelDD = new RooAddPdf("modelDD","modelDD",RooArgSet(*corrPdfDD,*bkgDD),*fracDD);
				modelDDs = new RooAddPdf("modelDDs","modelDDs",RooArgSet(*corrPdfs,*bkgDDs),*fracDD);
				
				normal = new Analysis(Form("normal_%i",e),cosThetaL);

			}
			else
			{
				pdf = new RooGenericPdf(Form("pdfB_%i_",e),afbBpdf,RooArgSet(*cosThetaB, *afbB) );
				bkgDDB_teo = stringToPdf(bkgmodel.c_str(),"bkgDDB_teo", cosThetaB);

				corrPdfB = new RooProdPdf(Form("corrPdfB_%i",e),*pdf,*effDD);
				corrPdfBs = new RooProdPdf(Form("corrPdfBs_%i",e),*pdf,*effDDs);
			
				bkgDDB = new RooProdPdf("bkgDDB","",*bkgDDB_teo,*effDDB);
				bkgDDBs = new RooProdPdf("bkgDDBs","",*bkgDDB_teo,*effDDBs);

				modelDD = new RooAddPdf("modelDD","modelDD",RooArgSet(*corrPdfDD,*bkgDD),*fracDD);
				modelDDs = new RooAddPdf("modelDDs","modelDDs",RooArgSet(*corrPdfs,*bkgDDs),*fracDD);

				normal = new Analysis(Form("normal_%i",e),cosThetaB);
			}

			c1->setVal(c1v);
			c2->setVal(c2v);
			cB1->setVal(c1Bv);
			cB2->setVal(c2Bv);
			afb->setVal(0.);
			afbB->setVal(-0.37);
			fL->setVal(0.8);

			normal->SetModel(corrPdf);
			normal->Generate(nevts[q2]);

			afb->setVal(0.);
			afbB->setVal(-0.37);
			fL->setVal(0.8);

			if(var=="cosThetaL") safeFit(corrPdf,(normal->GetDataSet()),params,&isInAllowedArea);
			else safeFit(corrPdf,(normal->GetDataSet()),paramsB,&isInAllowedAreaB);

			double norm_val, noeff_val, norm_fLval = 0, noeff_fLval = 0;
			if(var=="cosThetaL") { norm_val = afb->getVal(); norm_fLval = fL->getVal(); }
			else norm_val = afbB->getVal();

			afb->setVal(0.);
			afbB->setVal(-0.37);
			fL->setVal(0.8);

			if(systype=="flucteff")
			{
				//c1->setVal(r.Gaus(c1v,c1_err));
				//c2->setVal(r.Gaus(c2v,c2_err));
				//cB1->setVal(r.Gaus(c1Bv,c1B_err));
				//cB2->setVal(r.Gaus(c2Bv,c2B_err));

				RooDataHist * hsm = NULL;
				if(var=="cosThetaL")
				{
					hsm = new RooDataHist("hsm","hsm",*cosThetaL,Smear(effLL));
					c1->setConstant(false);
					c2->setConstant(false);
					effpdf->fitTo(*hsm,PrintLevel(-1),Save())->covQual();
					fixParams(effpdf,cosThetaL);
				}
				else
				{
					hsm = new RooDataHist("hBsm","hBsm",*cosThetaB,Smear(effLLB));
					cB1->setConstant(false);
					cB2->setConstant(false);
					effpdfB->fitTo(*hsm,PrintLevel(-1),Save())->covQual();
					fixParams(effpdfB,cosThetaB);
				}

				if(var=="cosThetaL") safeFit(corrPdf,(normal->GetDataSet()),params,&isInAllowedArea);
				else safeFit(corrPdf,(normal->GetDataSet()),paramsB,&isInAllowedAreaB);

				delete hsm;
			}
			else if(systype=="effsymm")
			{
				c1->setVal(0);
				cB1->setVal(0);

				corrPdf->fitTo(*(normal->GetDataSet()),PrintLevel(-1));
			}
			else if(systype=="noeff") pdf->fitTo(*(normal->GetDataSet()),PrintLevel(-1));
			else corrPdfs->fitTo(*(normal->GetDataSet()),PrintLevel(-1));

			if(var=="cosThetaL") { noeff_val = afb->getVal(); noeff_fLval = fL->getVal(); }
			else noeff_val = afbB->getVal();

			double rel_dev = norm_val - noeff_val;
			hsys->Fill(rel_dev);
			if(var=="cosThetaL") hfLsys->Fill(norm_fLval - noeff_fLval);
		}


		sys.push_back(hsys->GetMean());
		if(var=="cosThetaL") fLsys.push_back(hfLsys->GetMean());
	}

	cout << fixed << setprecision(5);

	for(int q2 = 0; q2 < q2nbins; q2++)
	{
		if(var=="cosThetaL") cout << sys[q2] << " & " << fLsys[q2]  << endl;
		else cout << sys[q2] << endl;
	}

}
