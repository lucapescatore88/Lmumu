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
#include "functions.hpp"
#include "fitfunc.hpp"
#include <fstream>

using namespace RooFit;
using namespace std;


RooFitResult * findMinRandom( RooAbsPdf * pdf, RooDataSet * data, Str2VarMap p, ISVALIDF_PTR isValid, string opt, int nfree, RooArgSet * cons, RooAbsReal * nll = NULL)
{
	RooRealVar cosThetaL("cosThetaL","cosThetaL",0.,-1.,1.);
	RooRealVar cosThetaB("cosThetaB","cosThetaB",0.,-1.,1.);
	RooRealVar MM("Lb_MassConsLambda","Lb_MassConsLambda",5621.,5540.,5900.);
	RooArgSet obs(cosThetaL,cosThetaB,MM);

	vector < string > pnames(1,"afb");
	pnames.push_back("afbB");
	pnames.push_back("fL");
	p = GetParams(pdf,obs,pnames);

	bool afb_iscost  = ((RooRealVar*)p["afb"])->getAttribute("Constant");
	bool fL_iscost   = ((RooRealVar*)p["fL"])->getAttribute("Constant");
	bool afbB_iscost = ((RooRealVar*)p["afbB"])->getAttribute("Constant");

	((RooRealVar*)p["afb"])->setConstant();
	((RooRealVar*)p["afbB"])->setConstant();
	((RooRealVar*)p["fL"])->setConstant();
	RooArgSet * nuisances = copyFreePars(pdf,obs);
	((RooRealVar*)p["afb"])->setConstant(afb_iscost);
	((RooRealVar*)p["afbB"])->setConstant(afbB_iscost);
	((RooRealVar*)p["fL"])->setConstant(fL_iscost);


	double precision = 0.001;
	if(!nll) nll = pdf->createNLL(*data);

	const unsigned np = 200;
	double a[np], f[np], ab[np];
	double c[3], r[1] = { 1e6 };

	genPoints(np,a,f,ab);
	while (r[0] > precision)
	{
		scan(pdf,data,nll,p,nfree,opt,cons,np,a,f,ab,c,r,nuisances);
		genPointsSphere(np,c,r,a,f,ab);
	}
	return scan(pdf,data,nll,p,nfree,opt,cons,np,a,f,ab,c,r,nuisances);
}



int main(int argc, char **argv)
{
	bool printeff = true;
	string fc = "";
	int dobin = -1;
	vector < vector < double > > points;

	gROOT->ProcessLine(".x ~/work/lhcbStyle.C");

	if(argc > 3 && ((string)argv[1]).find("[")!=string::npos && ((string)argv[2]).find("[")!=string::npos)
	{
		points.push_back( getList( argv[1] ) );
		points.push_back( getList( argv[2] ) );

		for (size_t i = 0; i < points[0].size(); i++)
		{
			cout << points[0][i] << "    "  << points[1][i] << endl;
		}

		printeff = false;
	}
	else cout << "No points set" << endl;

	if(argc > 1)
	{
		for(int a = 1; a < argc; a++)
		{
			string arg = argv[a];
			string str = arg.substr(2,arg.length()-2);

			if(arg.find("-E")!=string::npos) fc = str;
			if(arg=="-peff") printeff = true;
			if(arg.find("-B")!=string::npos) dobin = ((TString)str).Atof();
		}
	}

	int start = 1;
	int nbins = 2;
	if(dobin != -1) { start = dobin; nbins = dobin+1; }	
	double q2min[] = {8 ,15.,11.0,15,16,18,0.1};
	double q2max[] = {11,20.,12.5,16,18,20,2.};

	TFile * histFile = new TFile("Afb_hist.root","recreate");

	TString datafilename = "/afs/cern.ch/work/p/pluca/Lmumu/weighted/candLb.root";
	TreeReader * data = new TreeReader("candLb2Lmumu");
	data->AddFile(datafilename);
	TreeReader * datajpsi = new TreeReader("candLb2JpsiL");
	datajpsi->AddFile(datafilename);

	histFile->cd();

	string options = "-quiet-noPlot-lin-xM(#Lambda#mu#mu) (MeV/c^{2})-noCost-noParams";
	Analysis::SetPrintLevel("s");

	RooRealVar * cosThetaL = AngVars::cosThetaL;
	RooRealVar * cosThetaB = AngVars::cosThetaB;
	RooArgSet obs(* cosThetaL,* cosThetaB);
	RooRealVar * MM = new RooRealVar("Lb_MassConsLambda","Lb_MassConsLambda",5621.,5400.,5900.);
	RooRealVar * MMjpsi = new RooRealVar("Lb_MassConsJpsiLambda","Lb_MassConsJpsiLambda",5621.,5400.,5900.);
	double range[2] = { 5580, 5660 };
	MM->setRange("Signal",range[0],range[1]);
	MMjpsi->setRange("Signal",range[0],range[1]);
	RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

	TGraphAsymmErrors * Afb_vs_q2 = new TGraphAsymmErrors();
	TGraphAsymmErrors * AfbB_vs_q2 = new TGraphAsymmErrors();
	TGraphAsymmErrors * fL_vs_q2 = new TGraphAsymmErrors();
	TCanvas * ceff = new TCanvas();

	RooCategory * samples = new RooCategory("samples","samples");
	samples->defineType("DD");
	samples->defineType("LL");
	samples->defineType("sideLL");
	samples->defineType("sideDD");
	samples->defineType("DDB");
	samples->defineType("LLB");
	samples->defineType("sideLLB");
	samples->defineType("sideDDB");

	RooRealVar * afb = new RooRealVar("afb","A_{FB}^{l}",0.,-0.75,0.75);
	RooRealVar * fL = new RooRealVar("fL","f_{L}",0.6,0.,1.);
	TString afbLpdf = "((3./8.)*(1.-fL)*(1 + TMath::Power(cosThetaL,2)) + afb*cosThetaL + (3./4.)*fL*(1 - TMath::Power(cosThetaL,2)))";
	RooRealVar * afbB = new RooRealVar("afbB","A_{FB}^{h}",-0.3,-0.5,0.);
	TString afbBpdf = "(1 + 2*afbB*cosThetaB)";
	RooAbsPdf * teoPdf = new RooGenericPdf("teoPdf",afbLpdf,RooArgSet(*cosThetaL,*afb,*fL));
	RooAbsPdf * teoPdfB = new RooGenericPdf("teoPdfB",afbBpdf,RooArgSet(*cosThetaB,*afbB));

	TreeReader * mydata = datajpsi;

	/**               GET PARAMETERS                  **/

	Str2VarMap jpsiParsLL = getPars(8, 11, MM, "LL", CutsDef::LLcut, histFile, range);
	Str2VarMap jpsiParsDD = getPars(8, 11, MM, "DD", CutsDef::DDcut, histFile, range);

	cout << "Number of rare candidates: " << data->GetEntries() << endl;

	for(int i = start; i < nbins; i++)
	{
		TString q2name = ((TString)Form("q2_%4.2f_%4.2f",q2min[i],q2max[i])).ReplaceAll(".","");
		TString curq2cut = Form("TMath::Power(J_psi_1S_MM/1000,2) >= %e && TMath::Power(J_psi_1S_MM/1000,2) < %e",q2min[i],q2max[i]);
		if(q2min[i] == 8 && q2max[i] == 11) { q2name = "jpsi"; MM = MMjpsi; curq2cut += " && pplus_PIDp > 10"; } 
		else mydata = data;

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

		RooAbsPdf * corrPdfLL = new RooProdPdf("signal_LL_"+q2name,"signal",*teoPdf,*effLLpdf);
		RooAbsPdf * corrPdfDD = new RooProdPdf("signal_DD_"+q2name,"signal",*teoPdf,*effDDpdf);
		RooAbsPdf * corrPdfLLB = new RooProdPdf("signal_LL_B_"+q2name,"signal",*teoPdfB,*effLLBpdf);
		RooAbsPdf * corrPdfDDB = new RooProdPdf("signal_DD_B_"+q2name,"signal",*teoPdfB,*effDDBpdf);

		TCut baseCut = "";
		TCut cutLL = CutsDef::LLcut + (TCut)curq2cut + baseCut;
		TCut cutDD = CutsDef::DDcut + (TCut)curq2cut + baseCut;

		histFile->cd();
		double fracDDv[2], fracLLv[2];
		RooDataSet * dataLL = getDataAndFrac("LL", q2name, mydata, range, cutLL, MM, fracLLv,jpsiParsLL,NULL,true);
		RooDataSet * dataDD = getDataAndFrac("DD", q2name, mydata, range, cutDD, MM, fracDDv,jpsiParsDD,NULL,true);
		RooDataSet * sideDataLL = getSideData("LL", q2name, mydata, NULL, cutLL, MM);
		RooDataSet * sideDataDD = getSideData("DD", q2name, mydata, NULL, cutDD, MM);

		cout << "\n ************************ \n" << fixed << setprecision(2);
		cout << "Frac LL = " << fracLLv[0] << " +/- " << fracLLv[1];
		cout << ", Frac DD = " << fracDDv[0] << " +/- " << fracDDv[1];
		cout << "\n ************************ \n" << endl;
		RooRealVar * fracLL = new RooRealVar("fracLL","f",fracLLv[0],fracLLv[0]-fracLLv[1],fracLLv[0]+fracLLv[1]);
		RooRealVar * fracDD = new RooRealVar("fracDD","f",fracDDv[0],fracDDv[0]-fracDDv[1],fracDDv[0]+fracDDv[1]);

		RooAbsPdf *bkgLL, *bkgDD, *bkgLLB, *bkgDDB;
		if(q2name == "jpsi")
		{
			string bkgmodel = "Cheb2-v[0.,-0.7,0.7]";
			bkgLL = stringToPdf(bkgmodel.c_str(),"bkg_LL_teo", cosThetaL);
			bkgDD = stringToPdf(bkgmodel.c_str(),"bkg_DD_teo", cosThetaL); 
			bkgmodel = "Cheb2-v[0.,-0.6,0.]";
			bkgLLB = stringToPdf(bkgmodel.c_str(),"bkg_LLB_teo", cosThetaB);
			bkgDDB = stringToPdf(bkgmodel.c_str(),"bkg_DDB_teo", cosThetaB);
			bkgDD->SetTitle("background");
			bkgLL->SetTitle("background");
			bkgDDB->SetTitle("background");
			bkgLLB->SetTitle("background");

			bkgLL->fitTo(*sideDataLL,PrintLevel(-1));
			bkgLLB->fitTo(*sideDataLL,PrintLevel(-1));
			bkgDD->fitTo(*sideDataDD,PrintLevel(-1));
			bkgDDB->fitTo(*sideDataDD,PrintLevel(-1));
			fixParam(bkgLL,cosThetaL);	
			fixParam(bkgLLB,cosThetaB);
			fixParam(bkgDDB,cosThetaB);
			fixParam(bkgDD,cosThetaL);
			//constrainParams(bkgLL,cosThetaL,0.5);	
			//constrainParams(bkgLLB,cosThetaB,0.5);
			//constrainParams(bkgDDB,cosThetaB,0.5);
			//constrainParams(bkgDD,cosThetaL,0.5);
			printParams(bkgLL);
			printParams(bkgDD);
			printParams(bkgLLB);
			printParams(bkgDDB);
		}
		else
		{
			string bkgmodel = "Poly1-v[0.,-2,2]";
			RooAbsPdf * bkgLL_teo = stringToPdf(bkgmodel.c_str(),"background_LL_teo", cosThetaL);
			RooAbsPdf * bkgLLB_teo = stringToPdf(bkgmodel.c_str(),"background_LLB_teo", cosThetaB);
			RooAbsPdf * bkgDD_teo = stringToPdf(bkgmodel.c_str(),"background_DD_teo", cosThetaL); 
			RooAbsPdf * bkgDDB_teo = stringToPdf(bkgmodel.c_str(),"background_DDB_teo", cosThetaB);
			bkgLL = new RooProdPdf("bkg_LL","background",*bkgLL_teo,*effLLpdf);
			bkgDD = new RooProdPdf("bkg_DD","background",*bkgDD_teo,*effDDpdf);
			bkgLLB = new RooProdPdf("bkg_LL_B","background",*bkgLLB_teo,*effLLBpdf);
			bkgDDB = new RooProdPdf("bkg_DD_B","background",*bkgDDB_teo,*effDDBpdf);

			bkgLL->fitTo(*sideDataLL,PrintLevel(-1));
			bkgLLB->fitTo(*sideDataLL,PrintLevel(-1));
			bkgDD->fitTo(*sideDataDD,PrintLevel(-1));
			bkgDDB->fitTo(*sideDataDD,PrintLevel(-1));
			fixParam(bkgLL,cosThetaL);	
			fixParam(bkgLLB,cosThetaB);
			fixParam(bkgDDB,cosThetaB);
			fixParam(bkgDD,cosThetaL);
		}

		cout << "Backgrounds extracted" << endl;

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
		Str2VarMap params_both;
		params_both["fL"] = fL;
		params_both["afb"] = afb;
		params_both["afbB"] = afbB;

		RooSimultaneous * combModelB = new RooSimultaneous(Form("combModelB_%i",i),"",*samples);
		combModelB->addPdf(*modelLLB,"LL");
		combModelB->addPdf(*modelDDB,"DD");
		combModelB->addPdf(*bkgLLB,"sideLL");
		combModelB->addPdf(*bkgDDB,"sideDD");
		RooArgSet * origParsB = copyFreePars(combModelB,RooArgSet(*cosThetaB));
		RooArgSet * consB = gaussianConstraints(combModelB, obs);

		histFile->cd();
		RooFitResult * res = NULL;
		if(fc == "") res = safeFit(combModel,combData,params,&isInAllowedArea,"-scan-fast",100,cons);

		int nbinsLL = 8;
		int nbinsDD = 8;
		if(q2name=="jpsi") { nbinsLL = 40; nbinsDD = 40; }

		ceff->cd();
		string opts = "-noCost-plotSigComp-fillBkg-noParams";
		string extraopts = "";
		//if(q2name=="jpsi") extraopts += "-layout[0.65,0.2,0.90,0.60]";

		if(points.size()==0)
		{
			TLegend * leg = new TLegend(0.6,0.29,0.88,0.50);
			leg->SetFillColor(0);
			RooPlot * p = GetFrame(cosThetaL,modelLL,dataLL,opts+extraopts,nbinsLL,"cos #theta_{l}","",leg);
			p->addObject(leg);
			p->Draw();
			ceff->Print("Afb_LL_"+q2name+".pdf");
			p = GetFrame(cosThetaL,modelDD,dataDD,opts+extraopts,nbinsDD,"cos #theta_{l}","");
			p->addObject(leg);
			p->Draw();
			ceff->Print("Afb_DD_"+q2name+".pdf");
			GetFrame(cosThetaL,bkgLL,sideDataLL,opts,nbinsLL,"cos #theta_{l}")->Draw();
			ceff->Print("Side_LL_"+q2name+".pdf");
			GetFrame(cosThetaL,bkgDD,sideDataDD,opts,nbinsDD,"cos #theta_{l}")->Draw();
			ceff->Print("Side_DD_"+q2name+".pdf");
		}

		if(i!=1) Afb_vs_q2->SetPoint(i,(q2max[i] + q2min[i])/2.,afb->getVal());
		if(i!=1) fL_vs_q2->SetPoint(i,(q2max[i] + q2min[i])/2.,fL->getVal());

		histFile->cd();
		RooFitResult * resB = NULL;
		if(fc=="") resB = safeFit(combModelB,combData,paramsB,&isInAllowedAreaB,"-scan-fast",100,consB);

		ceff->cd();

		if(points.size()==0)
		{
			TLegend * leg = new TLegend(0.6,0.65,0.9,0.9);
			leg->SetFillColor(0);
			RooPlot * p = GetFrame(cosThetaB,modelLLB,dataLL,opts,nbinsLL,"cos #theta_{h}","",leg);
			p->addObject(leg);
			p->Draw();
			ceff->Print("AfbB_LL_"+q2name+".pdf");
			p = GetFrame(cosThetaB,modelDDB,dataDD,opts,nbinsDD,"cos #theta_{h}","");
			p->addObject(leg);
			p->Draw();
			ceff->Print("AfbB_DD_"+q2name+".pdf");
			GetFrame(cosThetaB,bkgLLB,sideDataLL,opts,nbinsLL,"cos #theta_{h}")->Draw();
			ceff->Print("SideB_LL_"+q2name+".pdf");
			GetFrame(cosThetaB,bkgDDB,sideDataDD,opts,nbinsDD,"cos #theta_{h}")->Draw();
			ceff->Print("SideB_DD_"+q2name+".pdf");

			if(i<2)
			{
				//Analysis::SetPrintLevel("v");
				RooRealVar * NSigDD = new RooRealVar("NSigDD","",fracDD->getVal()*dataDD->sumEntries());
				RooRealVar * NSigLL = new RooRealVar("NSigLL","",fracLL->getVal()*dataLL->sumEntries());
				RooRealVar * NBkgDD = new RooRealVar("NBkgDD","",(1.-fracDD->getVal())*dataDD->sumEntries());
				RooRealVar * NBkgLL = new RooRealVar("NBkgLL","",(1.-fracLL->getVal())*dataLL->sumEntries());
				RooRealVar * NSigDDB = new RooRealVar("NSigDDB","",fracDD->getVal()*dataDD->sumEntries());
				RooRealVar * NSigLLB = new RooRealVar("NSigLLB","",fracLL->getVal()*dataLL->sumEntries());
				RooRealVar * NBkgDDB = new RooRealVar("NBkgDDB","",(1.-fracDD->getVal())*dataDD->sumEntries());
				RooRealVar * NBkgLLB = new RooRealVar("NBkgLLB","",(1.-fracLL->getVal())*dataLL->sumEntries());

				Analysis * anaLL = new Analysis("LL","Lb",dataLL,cosThetaL,(RooAbsPdf*)NULL);
				anaLL->SetSignal(corrPdfLL,NSigLL,"-namepar");
				anaLL->addBkgComponent("Bkg",bkgLL,NBkgLL,"-namepar");
				anaLL->Initialize("-namepar");
				Analysis * anaDD = new Analysis("DD","Lb",dataDD,cosThetaL,(RooAbsPdf*)NULL);
				anaDD->SetSignal(corrPdfDD,NSigDD,"-namepar");
				anaDD->addBkgComponent("Bkg",bkgDD,NBkgDD,"-namepar");
				anaDD->Initialize("-namepar");

				Analysis * anaLLB = new Analysis("LLB","Lb",dataLL,cosThetaB,(RooAbsPdf*)NULL);
				anaLLB->SetSignal(corrPdfLLB,NSigLLB,"-namepar");
				anaLLB->addBkgComponent("Bkg",bkgLLB,NBkgLLB,"-namepar");
				anaLLB->Initialize("-namepar");
				Analysis * anaDDB = new Analysis("DDB","Lb",dataDD,cosThetaB,(RooAbsPdf*)NULL);
				anaDDB->SetSignal(corrPdfDDB,NSigDDB,"-namepar");
				anaDDB->addBkgComponent("Bkg",bkgDDB,NBkgDDB,"-namepar");
				anaDDB->Initialize("-namepar");

				MultiAnalysis * ana = new MultiAnalysis("fit_All");
				ana->AddCategory(anaLL,"LL");
				ana->AddCategory(anaDD,"DD");

				RooPlot * pl = ana->PrintSum(Form("-Xcos #theta_{l}-nbins%i-noParams-LHCb-nototsigplot",nbinsDD),cosThetaL);
				//pl->SetMaximum(70);
				pl->Draw();
				ceff->Print("fit_"+q2name+"_All_cosThetaL.pdf");
				ceff->Print("fit_"+q2name+"_All_cosThetaL.eps");
				ceff->Print("fit_"+q2name+"_All_cosThetaL.png");
				ceff->Print("fit_"+q2name+"_All_cosThetaL.C");

				MultiAnalysis * anaB = new MultiAnalysis("fit_All");
				anaB->AddCategory(anaLLB,"LLB");
				anaB->AddCategory(anaDDB,"DDB");

				pl = anaB->PrintSum(Form("-Xcos #theta_{h}-nbins%i-noParams-LHCbDX-nototsigplot",nbinsDD),cosThetaB);
				//pl->SetMaximum(90);
				pl->Draw();
				ceff->Print("fit_"+q2name+"_All_cosThetaB.pdf");
				ceff->Print("fit_"+q2name+"_All_cosThetaB.eps");
				ceff->Print("fit_"+q2name+"_All_cosThetaB.png");
				ceff->Print("fit_"+q2name+"_All_cosThetaB.C");

				//Analysis::SetPrintLevel("s");
			}
		}

		if(i!=1) AfbB_vs_q2->SetPoint(i,(q2max[i] + q2min[i])/2.,afbB->getVal());

		cout << endl << fixed << setprecision(6) << "AfbB = " << afbB->getVal() << " +/- " << afbB->getError() << endl;
		cout << "Afb = " << afb->getVal() << " +/- " << afb->getError() << endl;
		cout << "fL = " << fL->getVal() << " +/- " << fL->getError() << endl;
		cout << endl;
		if(res) cout << "lepton:  " << res->edm() << "   "  << res->covQual() << endl;
		if(resB) cout << "baryon:  " << resB->edm() << "   "  << resB->covQual() << endl;
		cout << endl;

		RooSimultaneous * combModel_both = new RooSimultaneous(Form("combModel_%i_both",i),"",*samples);
		combModel_both->addPdf(*modelLL,"LL");
		combModel_both->addPdf(*modelDD,"DD");
		combModel_both->addPdf(*bkgDD,"sideLL");
		combModel_both->addPdf(*bkgLL,"sideDD");
		combModel_both->addPdf(*modelLLB,"LLB");
		combModel_both->addPdf(*modelDDB,"DDB");
		combModel_both->addPdf(*bkgDDB,"sideLLB");
		combModel_both->addPdf(*bkgLLB,"sideDDB");

		RooArgSet * cons_both = gaussianConstraints(combModel_both, obs);
		RooFitResult * res_both = findMinRandom(combModel_both,combData,params_both,NULL,"",100,cons_both);

		if(points.size()==0)
		{
			GetFrame(cosThetaB,modelLLB,dataLL,opts,nbinsLL,"cos #theta_{h}")->Draw();
			ceff->Print("Both_AfbB_LL_"+q2name+".pdf");
			GetFrame(cosThetaB,modelDDB,dataDD,opts,nbinsDD,"cos #theta_{h}")->Draw();
			ceff->Print("Both_AfbB_DD_"+q2name+".pdf");
			GetFrame(cosThetaB,bkgLLB,sideDataLL,opts,nbinsLL,"cos #theta_{h}")->Draw();
			ceff->Print("Both_SideB_LL_"+q2name+".pdf");
			GetFrame(cosThetaB,bkgDDB,sideDataDD,opts,nbinsDD,"cos #theta_{h}")->Draw();
			ceff->Print("Both_SideB_DD_"+q2name+".pdf");
			GetFrame(cosThetaL,modelLL,dataLL,opts,nbinsLL,"cos #theta_{l}")->Draw();
			ceff->Print("Both_Afb_LL_"+q2name+".pdf");
			GetFrame(cosThetaL,modelDD,dataDD,opts,nbinsDD,"cos #theta_{l}")->Draw();
			ceff->Print("Both_Afb_DD_"+q2name+".pdf");
			GetFrame(cosThetaL,bkgLL,sideDataLL,opts,nbinsLL,"cos #theta_{l}")->Draw();
			ceff->Print("Both_Side_LL_"+q2name+".pdf");
			GetFrame(cosThetaL,bkgDD,sideDataDD,opts,nbinsDD,"cos #theta_{l}")->Draw();
			ceff->Print("Both_Side_DD_"+q2name+".pdf");
		}

		cout << endl << fixed << setprecision(6) << "AfbB = " << afbB->getVal() << " +/- " << afbB->getError() << endl;
		cout << "Afb = " << afb->getVal() << " +/- " << afb->getError() << endl;
		cout << "fL = " << fL->getVal() << " +/- " << fL->getError() << endl;
		cout << endl;
		if(res_both) cout << "both:  " << res_both->edm() << "   "  << res_both->covQual() << endl;
		cout << endl;


		if(fc!="")
		{ 
			cout << "------------------------ FELDMAN AND COUSINS ------------------------" << endl;

			histFile->Close();
			histFile = new TFile("Afb_hist.root","recreate");

			vector<double> afb_err(0,2), fL_err(0,2), afbB_err(0,2);
			vector < RooDataSet * > datas;
			vector < RooAbsPdf * > pdfs, pdfsB;
			vector < TString > cat;
			cat.push_back("LL");
			cat.push_back("DD");
			datas.push_back(dataLL);
			datas.push_back(dataDD);
			datas.push_back(sideDataLL);
			datas.push_back(sideDataDD);

			pdfs.push_back(modelLL);
			pdfs.push_back(modelDD);
			pdfs.push_back(bkgLL);
			pdfs.push_back(bkgDD);

			pdfsB.push_back(modelLLB);
			pdfsB.push_back(modelDDB);
			pdfsB.push_back(bkgLLB);
			pdfsB.push_back(bkgDDB);

			TGraph * best = new TGraph();
			vector< vector < double > > res;

			if(fc == "all" || fc == "afbfL")
			{
				FeldmanCousins * FC = new FeldmanCousins(q2name,cat,datas,pdfs,cosThetaL,params);
				FC->AddConstraints(cons);
				FC->SetFitFunc(&safeFit);
				if(q2name=="jpsi") FC->SetNExp(500);
				else FC->SetNExp(1000);
				res = FC->ExtractLimits(points,origPars,&best,"-scan-fast");
			}
			if(fc == "afb")
			{
				FeldmanCousins * FC = new FeldmanCousins(q2name,cat,datas,pdfs,cosThetaL,afb);
				FC->AddConstraints(cons);
				FC->SetFitFunc(&safeFit);
				if(q2name=="jpsi") FC->SetNExp(500);
				else FC->SetNExp(1000);
				res = FC->ExtractLimits(points,origPars,&best,"-scan-fast");
			}
			if(fc == "fL")
			{
				FeldmanCousins * FC = new FeldmanCousins(q2name,cat,datas,pdfs,cosThetaL,fL);
				FC->AddConstraints(cons);
				FC->SetFitFunc(&safeFit);
				if(q2name=="jpsi") FC->SetNExp(500);
				else FC->SetNExp(1000);
				res = FC->ExtractLimits(points,origPars,&best,"-scan-fast");
			}

			if(fc == "all" || fc == "afbB")
			{
				FeldmanCousins * FC = new FeldmanCousins(q2name,cat,datas,pdfsB,cosThetaB,afbB);
				FC->AddConstraints(cons);
				FC->SetFitFunc(&safeFit);
				if(q2name=="jpsi") FC->SetNExp(500);
				else FC->SetNExp(1000);
				res = FC->ExtractLimits(points,origParsB,&best,"-scan-fast");
			}

			histFile->cd();
			best->Write("bestfit");

			double x,y;
			best->GetPoint(0,x,y);
			cout << "\n\n ***** Best fit: " << x << "    " << y << endl;
			for(unsigned kk = 0; kk < res.size(); kk++)
			{
				cout << "\n\n ----- Pvalues: ";
				for(unsigned ii = 0; ii < res[0].size(); ii++)
					cout << res[kk][ii] << "    ";
				cout << endl << endl;
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
	Afb_vs_q2->SetMarkerSize(1);
	Afb_vs_q2->SetMarkerStyle(20);
	Afb_vs_q2->Draw("AP");
	if(points.size()==0) ceff->Print("Afb_vs_q2.pdf");
	AfbB_vs_q2->GetXaxis()->SetTitle("q^{2}");
	AfbB_vs_q2->GetYaxis()->SetTitle("A_{FB}^{h}");
	AfbB_vs_q2->SetMaximum(1);
	AfbB_vs_q2->SetMinimum(-1);
	AfbB_vs_q2->SetMarkerSize(1);
	AfbB_vs_q2->SetMarkerStyle(20);
	AfbB_vs_q2->Draw("AP");
	if(points.size()==0) ceff->Print("AfbB_vs_q2.pdf");
	fL_vs_q2->GetXaxis()->SetTitle("q^{2}");
	fL_vs_q2->GetYaxis()->SetTitle("f_{L}");
	fL_vs_q2->Draw("AP");
	if(points.size()==0) ceff->Print("fL_vs_q2.pdf");

	for(int bb = 0; bb < Afb_vs_q2->GetN(); bb++)
	{
		double qq, afbv, afbBv, fLv;
		Afb_vs_q2->GetPoint(bb,qq,afbv);
		AfbB_vs_q2->GetPoint(bb,qq,afbBv);
		fL_vs_q2->GetPoint(bb,qq,fLv);
		cout << fixed << setprecision(1) << q2min[bb] << " - " << q2max[bb];
		cout << fixed << setprecision(4); 

		if(fc == "all" || fc == "afb" || fc == "afbfL") cout << " & $" << afbv << "_{-" << Afb_vs_q2->GetErrorYlow(bb) << "}^{+" << Afb_vs_q2->GetErrorYhigh(bb)  << "} \\text{(stat)} \\pm \\text{(sys)}$ ";
		if(fc == "all" || fc == "afbB") cout << " & $" << afbBv << "_{-" << AfbB_vs_q2->GetErrorYlow(bb) << "}^{+" << AfbB_vs_q2->GetErrorYhigh(bb) << "} \\text{(stat)} \\pm \\text{(sys)}$ " ;
		if(fc == "all" || fc == "afbfL" || fc == "fL") cout << " & $" << fLv << "_{-" << fL_vs_q2->GetErrorYlow(bb) << "}^{+" << fL_vs_q2->GetErrorYhigh(bb)  << "} \\text{(stat)} \\pm \\text{(sys)}$ ";
		if(fc == "") cout << " & \t" << afbv << " & \t" << afbBv << " & \t" << fLv;
		cout << "  \\\\ " << endl;
	}

	//histFile->Close();
}


