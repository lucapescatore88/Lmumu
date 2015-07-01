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
	int nexp = 500;
	int dobin = -1;

	gROOT->ProcessLine(".x ~/work/lhcbStyle.C");

	if(argc > 1)
	{
		for(int a = 1; a < argc; a++)
		{
			string arg = argv[a];
			string str = arg.substr(2,arg.length()-2);

			if(arg.find("-B")!=string::npos) dobin = ((TString)str).Atof();
			if(arg.find("-E")!=string::npos) nexp = ((TString)str).Atof();
		}
	}

	int nbins = 2;
	int start = 1;
	if(dobin!=-1) { start = dobin; nbins = dobin+1; }
	double q2min[] = {8.,15.,11.0,15,16,18, 0.1};
	double q2max[] = {11.,20.,12.5,16,18,20, 2.};

	TString datafilename = "/afs/cern.ch/work/p/pluca/Lmumu/weighted/candLb.root";
	TreeReader * data = new TreeReader("candLb2Lmumu");
	data->AddFile(datafilename);
	TreeReader * datajpsi = new TreeReader("candLb2JpsiL");
	datajpsi->AddFile(datafilename);

	TFile * histFile = new TFile("Afb_effSys.root","recreate");

	string options = "-quiet-noPlot-lin-XM(#Lambda#mu#mu) (MeV/c^{2})-noCost-noParams";
	Analysis::SetPrintLevel("s");

	RooRealVar * cosThetaL = new RooRealVar("cosThetaL","cosThetaL",0.,-1.,1.);
	RooRealVar * cosThetaB = new RooRealVar("cosThetaB","cosThetaB",0.,-1.,1.);
	RooRealVar * MM = new RooRealVar("Lb_MassConsLambda","Lb_MassConsLambda",5621.,5450.,6000.);

	double trange[2] = {5590,5650};
	double * range = &trange[0];
	MM->setRange("Signal",5590,5650);
	MM->Print();

	RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

	RooCategory * samples = new RooCategory("samples","samples");
	samples->defineType("DD");

	TString effbase = "/afs/cern.ch/user/p/pluca/work/Lb/Lmumu/results/";
	TFile * effFile = TFile::Open(effbase+"Lbeff2D_cosThetaB_vs_cosThetaL_All.root");
	TH1 * eff  = (TH1 *)effFile->Get("htot_eff");
	RooDataHist * heff = new RooDataHist("heff","heff",RooArgSet(*cosThetaL,*cosThetaB),eff);
	RooHistPdf  * eff2Dpdf = new RooHistPdf("eff2Dpdf","eff2Dpdf",RooArgSet(*cosThetaL,*cosThetaB),*heff);

	RooRealVar * fL = new RooRealVar("fL","fL",0.7,0,1);
	RooRealVar * afb = new RooRealVar("afb","afb",0,-1,1);
	RooRealVar * afbB = new RooRealVar("afbB","afbB",-0.37,-1,1);
	RooRealVar * alphaL = new RooRealVar("alphaL","alphaL",0.642);
	RooRealVar * Op = new RooRealVar("Op","Op",-0.461);
	RooRealVar * Ou12 = new RooRealVar("Ou12","Ou12",-0.302);

	RooGenericPdf * distr =  new RooGenericPdf("distr",
			"(3./8.) + (3./8.)*(cosThetaL^2)*(1-fL) - (3./16.)*fL*(cosThetaL^2) + afb*cosThetaL"  
			"+ ((3./2.)*afbB - (3./8.)*alphaL*Op)*cosThetaB - (3./2.)*afbB*(cosThetaL^2)*cosThetaB - (3./16.)*fL"
			"+ (9./16.)*fL*(1 - cosThetaL^2) + (9./8.)*alphaL*Op*(cosThetaL^2)*cosThetaB" 
			"- (3./2.)*alphaL*Ou12*cosThetaL*cosThetaB",
			RooArgSet(*cosThetaL,*cosThetaB,*fL,*afb,*Op,*Ou12,*alphaL,*afbB) );

	RooProdPdf * pdf = new RooProdPdf("pdf","pdf",*distr,*eff2Dpdf);

	TString afbLpdf = "((3./8.)*(1.-fL)*(1 + TMath::Power(cosThetaL,2)) + afb*cosThetaL + (3./4.)*fL*(1 - TMath::Power(cosThetaL,2)))";
	TString afbBpdf = "(1 + 2*afbB*cosThetaB)";
	RooAbsPdf * teoPdf = new RooGenericPdf("teoPdf",afbLpdf,RooArgSet(*cosThetaL,*afb,*fL));
	RooAbsPdf * teoPdfB = new RooGenericPdf("teoPdfB",afbBpdf,RooArgSet(*cosThetaB,*afbB));

	for(int i = start; i < nbins; i++)
	{
		TreeReader * mydata = data;
		TString q2name = ((TString)Form("q2_%4.2f_%4.2f",q2min[i],q2max[i])).ReplaceAll(".","");
		if(q2min[i]==8 && q2max[i]==11){ mydata = datajpsi; q2name = "jpsi"; }
		TString curq2cut = Form("TMath::Power(J_psi_1S_MM/1000,2) >= %e && TMath::Power(J_psi_1S_MM/1000,2) < %e",q2min[i],q2max[i]);	

		cout << "------------------- q2 bin: " << q2min[i] << " - " << q2max[i] << " -----------------------" << endl;


		/**               GET MASS PARAMETERS                  **/


		Str2VarMap jpsiParsLL = getPars(q2min[i],q2max[i],MM,"LL", CutsDef::LLcut, histFile);
		Str2VarMap jpsiParsDD = getPars(q2min[i],q2max[i],MM,"DD", CutsDef::DDcut, histFile);


		/**               GET AND FIT EFFICIENCIES                  **/

		RooAbsPdf * effDDpdf = NULL, * effLLpdf = NULL, * effLLBpdf = NULL, * effDDBpdf = NULL;	
		getEfficiencies(q2min[i],q2max[i],&effLLpdf,&effDDpdf,&effLLBpdf,&effDDBpdf,printeff);

		cout << "Efficiencies extracted" << endl;
		histFile->cd();


		/**                    FIT AFB                  **/

		double afbgen = 0, afbBgen = -0.37, fLgen = 0.6; 

		afb->setVal(afbgen);
		afbB->setVal(afbBgen);
		fL->setVal(fLgen);

		RooAbsPdf * corrPdfLL = new RooProdPdf("sigPdfLL"+q2name,"corrPdfLL",*teoPdf,*effLLpdf);
		RooAbsPdf * corrPdfDD = new RooProdPdf("sigPdfDD"+q2name,"corrPdfDD",*teoPdf,*effDDpdf);
		RooAbsPdf * corrPdfLLB = new RooProdPdf("sigPdfLLB"+q2name,"corrPdfLLB",*teoPdfB,*effLLBpdf);
		RooAbsPdf * corrPdfDDB = new RooProdPdf("sigPdfDDB"+q2name,"corrPdfDDB",*teoPdfB,*effDDBpdf);

		TCut baseCut = "";
		TCut cutLL = CutsDef::LLcut + (TCut)curq2cut + baseCut;
		TCut cutDD = CutsDef::DDcut + (TCut)curq2cut + baseCut;

		histFile->cd();

		double fracDDv[2], fracLLv[2];
		double nsigDD, nsigLL, slopeLL, slopeDD;
		RooDataSet * dataLL = getDataAndFrac("LL",q2name,mydata,range,cutLL,MM,&fracLLv[0],&slopeLL,jpsiParsLL,&nsigLL,true);
		RooDataSet * dataDD = getDataAndFrac("DD",q2name,mydata,range,cutDD,MM,&fracDDv[0],&slopeDD,jpsiParsDD,&nsigDD,true);
		RooDataSet * sideDataLL = getSideData("LL", q2name, mydata, NULL, cutLL, MM);
		RooDataSet * sideDataDD = getSideData("DD", q2name, mydata, NULL, cutDD, MM);

		RooRealVar * fracLL = new RooRealVar("fracLL","fracLL",fracLLv[0],0.01,0.99);
		RooRealVar * fracDD = new RooRealVar("fracDD","fracDD",fracDDv[0],0.01,0.99);
		
		RooAbsPdf *bkgLL, *bkgDD, *bkgLLB, *bkgDDB;
		string bkgmodel = "Poly1-v[0.,-2,2]";
		RooAbsPdf * bkgLL_teo = stringToPdf(bkgmodel.c_str(),"bkgLL_teo", cosThetaL);
		RooAbsPdf * bkgLLB_teo = stringToPdf(bkgmodel.c_str(),"bkgLLB_teo", cosThetaB);
		RooAbsPdf * bkgDD_teo = stringToPdf(bkgmodel.c_str(),"bkgDD_teo", cosThetaL); 
		RooAbsPdf * bkgDDB_teo = stringToPdf(bkgmodel.c_str(),"bkgDDB_teo", cosThetaB);
		bkgLL = new RooProdPdf("bkgLL","bkgLL",*bkgLL_teo,*effLLpdf);
		bkgDD = new RooProdPdf("bkgDD","bkgDD",*bkgDD_teo,*effDDpdf);
		bkgLLB = new RooProdPdf("bkgLLB","bkgLLB",*bkgLLB_teo,*effLLBpdf);
		bkgDDB = new RooProdPdf("bkgDDB","bkgDDB",*bkgDDB_teo,*effDDBpdf);
	
		bkgLL->fitTo(*sideDataLL,PrintLevel(-1));
		bkgLLB->fitTo(*sideDataLL,PrintLevel(-1));
		bkgDD->fitTo(*sideDataDD,PrintLevel(-1));
		bkgDDB->fitTo(*sideDataDD,PrintLevel(-1));
		fixParam(bkgLL,cosThetaL);	
		fixParam(bkgLLB,cosThetaB);
		fixParam(bkgDDB,cosThetaB);

		RooAbsPdf * modelLL = new RooAddPdf("modelLL","modelLL",RooArgSet(*corrPdfLL,*bkgLL),*fracLL);
		RooAbsPdf * modelDD = new RooAddPdf("modelDD","modelDD",RooArgSet(*corrPdfDD,*bkgDD),*fracDD);
		RooAbsPdf * modelLLB = new RooAddPdf("modelLLB","modelLLB",RooArgSet(*corrPdfLLB,*bkgLLB),*fracLL);
		RooAbsPdf * modelDDB = new RooAddPdf("modelDDB","modelDDB",RooArgSet(*corrPdfDDB,*bkgDDB),*fracDD);

		// CREATE COMBINED DATASET
		map <string, RooDataSet *> dataMap;
		dataMap["DD"] = dataDD;
		dataMap["LL"] = dataLL;
		dataMap["sideDD"] = sideDataDD;
		dataMap["sideLL"] = sideDataLL;
		dataMap["DDB"] = dataDD;
		dataMap["LLB"] = dataLL;
		dataMap["sideDDB"] = sideDataDD;
		dataMap["sideLLB"] = sideDataLL;
		RooDataSet * combData = new RooDataSet(Form("combData_%i",i),"combined data",
			RooArgSet(*MM,*cosThetaL,*cosThetaB),Index(*samples),Import(dataMap));
		
		RooSimultaneous * combModel = new RooSimultaneous(Form("combModel_%i",i),"",*samples);
		combModel->addPdf(*modelLL,"LL");
		combModel->addPdf(*modelDD,"DD");
		combModel->addPdf(*bkgDD,"sideLL");
		combModel->addPdf(*bkgLL,"sideDD");
		RooArgSet * origPars = copyFreePars(combModel,RooArgSet(*cosThetaL));
		RooArgSet * cons = gaussianConstraints(combModel, obs);

		Str2VarMap params;
		params["fL"] = fL;
		params["afb"] = afb;	
		Str2VarMap paramsB;
		paramsB["afbB"] = afbB;
		
		RooSimultaneous * combModelB = new RooSimultaneous(Form("combModelB_%i",i),"",*samples);
		combModelB->addPdf(*modelLLB,"LL");
		combModelB->addPdf(*modelDDB,"DD");
		combModelB->addPdf(*bkgLLB,"sideLL");
		combModelB->addPdf(*bkgDDB,"sideDD");
		RooArgSet * origParsB = copyFreePars(combModelB,RooArgSet(*cosThetaB));
		RooArgSet * consB = gaussianConstraints(combModelB, obs);

		histFile->cd();
		RooFitResult * res = safeFit(combModel,combData,params,&isInAllowedArea,"-scan-fast",100,cons);
		RooFitResult * resB = safeFit(combModelB,combData,paramsB,&isInAllowedAreaB,"-scan-fast",100,consB);
		
		int ngen = dataDD->sumEntries()*fracDDv[0] + dataLL->sumEntries()*fracLLv[0]

		TH1F * fLsys = new TH1F("fLsys","",40,-1,1);
		TH1F * afbsys = new TH1F("afbsys","",40,-1,1);
		TH1F * afbBsys = new TH1F("afbBsys","",40,-1,1);

		TCanvas * c = new TCanvas();

		for(int e = 0; e < nexp; e++)
		{
			showPercentage(e,nexp);

			afb->setVal(afbgen);
			afbB->setVal(afbBgen);
			fL->setVal(fLgen);
			fracDD->setVal(fracDDv[0]);
			fracLL->setVal(fracLLv[0]);

			RooDataSet * toyDataDD = generateDataSet("DD",setL,sigDD,nsigDD*fracDDv[0],totBkgDD,nsigDD*(1.-fracDDv[0]));
			RooDataSet * toyDataDDB = generateDataSet("DDB",setB,sigDDB,nsigDD*fracDDv[0],totBkgDDB,nsigDD*(1.-fracDDv[0]));
			RooDataSet * toyDataLL = generateDataSet("LL",setL,sigLL,nsigLL*fracLLv[0],totBkgLL,nsigLL*(1.-fracLLv[0]));
			RooDataSet * toyDataLLB = generateDataSet("LLB",setB,sigLLB,nsigLL*fracLLv[0],totBkgLLB,nsigLL*(1.-fracLLv[0]));

			RooDataSet * toyData = new RooDataSet(Form("toyData_%i",i),"combined data",
					RooArgSet(*MM,*cosThetaL),Index(*samples),Import("DD",*toyDataDD),Import("LL",*toyDataLL));
			RooDataSet * toyDataB = new RooDataSet(Form("toyDataB_%i",i),"combined data",
					RooArgSet(*MM,*cosThetaB),Index(*samples),Import("DD",*toyDataDDB),Import("LL",*toyDataLLB));


		//	RooDataSet * toy = pdf->generate(RooArgSet(*cosThetaL,*cosThetaB),ngen);

			afb->setVal(afbgen);
			afbB->setVal(afbBgen);
			fL->setVal(fLgen);
			fracDD->setVal(fracDDv[0]);
			fracLL->setVal(fracLLv[0]);


			safeFit(toyData,combData,params,&isInAllowedArea,"-scan-fast",100,cons);

			double afbv = afb->getVal();
			double fLv = fL->getVal();

			afb->setVal(afbgen);
			afbB->setVal(afbBgen);
			fL->setVal(fLgen);
			fracDD->setVal(fracDDv[0]);
			fracLL->setVal(fracLLv[0]);


			safeFit(toyDataB,combData,paramsB,&isInAllowedAreaB,"-scan-fast",100,consB);

			double afbBv = afbB->getVal();

			fLsys->Fill(fLv - fLgen);
			afbsys->Fill(afbv - afbgen);
			afbBsys->Fill(afbBv - afbBgen);
		}

		fLsys->GetXaxis()->SetTitle("#Delta f_{L}");
		fLsys->Draw();
		c->Print("fLsys_efficiency"+q2name+".pdf");
		afbsys->GetXaxis()->SetTitle("#Delta A_{FB}");
		afbsys->Draw();
		c->Print("afbsys_efficiency"+q2name+".pdf");
		afbBsys->GetXaxis()->SetTitle("#Delta A_{FB}^{h}");
		afbBsys->Draw();
		c->Print("afbBsys_efficiency"+q2name+".pdf");

		cout << fixed << setprecision(7);
		cout << "fLsys = " << fLsys->GetMean() << " +/- " << fLsys->GetMeanError() << endl;
		cout << "afbsys = " << afbsys->GetMean() << " +/- " << afbsys->GetMeanError() << endl;
		cout << "afbBsys = " << afbBsys->GetMean() << " +/- " << afbBsys->GetMeanError() << endl;
	}
}
