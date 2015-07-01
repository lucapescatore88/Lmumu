#include "TGraphErrors.h"
#include "TPaveText.h"
#include "general_functions.hpp"
#include "ReadTree_comp.hpp"
#include "analyser.hpp"
#include "Lb_cuts.hpp"
#include "multi_analyser.hpp"


using namespace std;
using namespace RooFit;

double factorLL = 0, factorDD = 0;
bool wilson = false;
string tname = "cand";
string sys = "";
TFile * histFile = NULL;
TCut singleCut = "";
ofstream KSpredict_file, effs_file;
map<double,double> pideff, pdfsys;
string model = "DCB_Sn";
map<string, vector<double> >  w1DD, w2DD, w3DD, w1LL, w2LL, w3LL;
RooRealVar * vLbJpsiCons = new RooRealVar("Lb_MassConsJpsiLambda","Lb_MassConsJpsiLambda",5621.5,5350.,6100.);
double yield_rare_LL_MC = 0, yield_rare_DD_MC = 0, yield_jpsi_LL_MC = 0, yield_jpsi_DD_MC = 0;
TCut yCut = "";
TCut pCut = "";

Analysis * fitJpsi(TString type, TCut baseCut, RooRealVar ** NJpsi, Str2VarMap * MCpars, double *KS_double_eff_ratio, double *nKSjpsi_val, bool lowQ2 = false)
{
	string candfilename = "/afs/cern.ch/work/p/pluca/Lmumu/weighted/candLb_MC.root";
	string datafilename = "/afs/cern.ch/work/p/pluca/Lmumu/weighted/candLb.root";
	histFile->cd();

	if(lowQ2) baseCut += CutsDef::cutJpsi_lowSel;
	TString jpsiname = type;
	if(lowQ2) jpsiname = "lowSel_"+type;
	Analysis * anaLbJpsi = new Analysis("Lb2JpsiL__"+jpsiname+"_MC","#Lambda_{b} #rightarrow J\\#psi#Lambda",tname+"Lb2JpsiL",candfilename,vLbJpsiCons);
	Analysis * anaLbJpsi_data = new Analysis("Lb2JpsiL__"+jpsiname+"_data","#Lambda_{b} #rightarrow J\\#psi#Lambda",tname+"Lb2JpsiL",datafilename,vLbJpsiCons);
	anaLbJpsi->SetUnits("MeV/#it{c}^{2}");
	anaLbJpsi_data->SetUnits("MeV/#it{c}^{2}");

	TFile * MCFile = TFile::Open((TString)candfilename);
	TTree * BdKSmumuTree  = (TTree *)MCFile->Get("candBdKSmumu");
	TTree * BdJpsiKSTree  = (TTree *)MCFile->Get("candBdJpsiKS");
	//TTree * BuKstmumuTree  = (TTree *)MCFile->Get("candBuKstmumu");
	MCFile = TFile::Open("/afs/cern.ch/work/p/pluca/Lmumu/weighted/Bd2KSmumu_MC12_NBweighted.root");
	TTree * BdKSmumuTree_All  = (TTree *)MCFile->Get("EventTuple");
	MCFile = TFile::Open("/afs/cern.ch/work/p/pluca/Lmumu/weighted/Bd2JpsiKS_MC12_NBweighted.root");
	TTree * BdJpsiKSTree_All  = (TTree *)MCFile->Get("EventTuple");
	double KSmumuSelected = BdKSmumuTree->Draw("J_psi_1S_MM",baseCut);
	double JpsiKSSelected = BdJpsiKSTree->Draw("J_psi_1S_MM",baseCut);
	double effKSmumu = KSmumuSelected / (double) BdKSmumuTree_All->GetEntries();
	double effJpsiKS = JpsiKSSelected / (double) BdJpsiKSTree_All->GetEntries();
	(*KS_double_eff_ratio) = effKSmumu / effJpsiKS;

	histFile->cd();

	string optionsjpsi = "-quiet-nototsigplot-range-lin-noleg-xuM(#Lambda#mu#mu)-noCost-nochi2-min2-layout[0.55,0.3,0.9,0.9]";
	anaLbJpsi->SetSignal((model+"-s[7]-s2[15]").c_str(),1.1,"-namepar");
	anaLbJpsi->Initialize("-namepar");
	anaLbJpsi->Fit(5400.,5750.,70,true,"-quiet-xuM(#Lambda#mu#mu)-noCost-noleg-nochi2-layout[0.55,0.3,0.9,0.9]",baseCut);
	(*MCpars) = anaLbJpsi->GetSigParams();
	if(type=="LL") yield_jpsi_LL_MC = anaLbJpsi->GetNSigVal();
	else if (type=="DD") yield_jpsi_DD_MC = anaLbJpsi->GetNSigVal();

	Analysis * KS = new Analysis("JpsiKS_"+type,"Lb",BdJpsiKSTree,vLbJpsiCons,"DCB_OST","","-namepar");
	KS->Fit(5350.,6000.,75,true,"-quiet-noParams-xuM(#Lambda#mu#mu)-noleg");
	Str2VarMap pars = KS->GetSigParams();
	RooRealVar * m_shift = new RooRealVar("m_shift_DD","shift",6.,-10.,10.);
	if(type=="LL") { m_shift = new RooRealVar("m_shift_LL","shift",0.,-10.,10.); m_shift->setConstant(); }
	setConstant(&pars);
	ModifyPars(&pars,"m",m_shift,"-shift");

	//Fit Jpsi
	string jpsimodel = model;
	//+"-Xn"+Form("[%f]",(*MCpars)["n"]->getVal());
	//jpsimodel += (string)"-s[7,1,15]";
	//jpsimodel += (string)"-s2[20,5,50]";
	//jpsimodel += (string)"-Xa"+Form("[%f,3,10]",(*MCpars)["a"]->getVal());
	//jpsimodel += (string)"-Xa2"+Form("[%f,0,10]",(*MCpars)["a2"]->getVal());
	//jpsimodel += (string)"-Xf"+Form("[%f,0,0.9]",(*MCpars)["f"]->getVal());
	setConstant(MCpars);
	((RooRealVar*)(*MCpars)["s"])->setConstant(0);
	((RooRealVar*)(*MCpars)["s2"])->setConstant(0);

	histFile->cd();
	anaLbJpsi_data->SetSignal(jpsimodel.c_str(),2000,"-namepar",(*MCpars));
	anaLbJpsi_data->addBkgComponent("Comb","Exp-b[-0.003,-0.01,0.]",100.,"-namepar");
	RooRealVar * nKSjpsi = new RooRealVar("nKSjpsi_"+jpsiname,"N_{K_{S}J\\psi}",1.e2,0,1.e4);
	anaLbJpsi_data->addBkgComponent("JpsiKS","DCB_OST",nKSjpsi,"-namepar",pars);
	anaLbJpsi_data->Initialize("-namepar-docuts");
	if(type=="LL") anaLbJpsi_data->Fit(5350.,6000.,75,false,optionsjpsi,baseCut + singleCut + yCut + pCut);
	else anaLbJpsi_data->Fit(5350.,6000.,75,true,optionsjpsi,baseCut + singleCut + yCut + pCut);
	anaLbJpsi_data->Print((string)(((TString)optionsjpsi).ReplaceAll("-lin","-log"))+"-noParams-noleg",150);
	
	(*NJpsi) = (RooRealVar*)anaLbJpsi_data->GetNSigPtr();
	(*nKSjpsi_val) = nKSjpsi->getVal();
	double nKSjpsi_err = nKSjpsi->getError();
	double BRJpsiKS = ((*nKSjpsi_val)/(effJpsiKS*0.16*2))/((*NJpsi)->getVal()/(0.00645));
	if(type=="LL") BRJpsiKS = ((*nKSjpsi_val)/(effJpsiKS*0.16*2))/((*NJpsi)->getVal()/(0.00222));
	cout << "BR(JpsiKS)/BR(Lmumu) = " << BRJpsiKS << " +/- " << BRJpsiKS*TMath::Sqrt(TMath::Power(nKSjpsi_err/(*nKSjpsi_val),2) + TMath::Power((*NJpsi)->getError()/(*NJpsi)->getVal(),2)) << endl;

	m_shift->setConstant(true);

	return anaLbJpsi_data;
}


map<string, vector<double> > mergeGr(TString title, TString f1name, TString f2name, bool isTH1 = false)
{
	TFile * f1 = TFile::Open(f1name);
	TGraphErrors * gr1 = new TGraphErrors();
	if(isTH1)
	{
		TH1F * hh = 0;
		f1->GetObject(title,hh);
		gr1 = new TGraphErrors(hh);
	}
	else f1->GetObject(title,gr1);
	TFile * f2 = TFile::Open(f2name);
	TGraphErrors * gr2 = new TGraphErrors();
	int init = 1;
	if(isTH1)
	{
		TH1F * hh = 0;
		f2->GetObject(title,hh);
		gr2 = new TGraphErrors(hh);
		init = 0;
	}
	else f2->GetObject(title,gr2);

	map<string, vector<double> > mymap;

	for(int i = init; i < gr1->GetN(); i++) 
	{
		double x,y;
		gr1->GetPoint(i,x,y);
		if(x<0.1) continue;
		vector<double> v(1,y);
		v.push_back(gr1->GetErrorY(i));
		mymap[Form("%4.2f",(float)x)] = v;
	}	
	for(int i = init; i < gr2->GetN(); i++) 
	{
		double x,y;
		gr2->GetPoint(i,x,y);
		if(x<0.1) continue;
		vector<double> v(1,y);
		v.push_back(gr1->GetErrorY(i));
		mymap[Form("%4.2f",(float)x)] = v;
	}
	return mymap;
}


Analysis * getAnaObject(TString type, TCut curcut, double q2min, double q2max, RooRealVar * vLbCons, RooRealVar * BR,
		map<string , vector <double> > eff, RooRealVar * Njpsi, Str2VarMap rarePars, Str2VarMap MCpars,
		double KS_double_eff_ratio, double nKSjpsi_val, string dosig = "BR")
{
	string options = "-quiet-noleg-nototsigplot-lin-range-minos-xuM(#Lambda#mu#mu)-noCost-noParams";
	string candfilename = "/afs/cern.ch/work/p/pluca/Lmumu/weighted/candLb_MC.root";
	string datafilename = "/afs/cern.ch/work/p/pluca/Lmumu/weighted/candLb.root";
	histFile->cd();

	TString q2name = ((TString)Form("_q2_%4.2f_%4.2f",q2min,q2max)).ReplaceAll(".","");
	TCut dataCut = curcut + yCut + pCut;
	Analysis * anaLbMuMu = new Analysis("Lb2Lmumu_MC_"+type+q2name,"#Lambda_{b} #rightarrow #Lambda#mu#mu",tname+"Lb2Lmumu",candfilename,vLbCons,&curcut);
	Analysis * anaLbMuMu_data = new Analysis("Lb2Lmumu_"+type+q2name,"#Lambda_{b} #rightarrow #Lambda#mu#mu",tname+"Lb2Lmumu",datafilename,vLbCons,&dataCut);
	anaLbMuMu->SetUnits("MeV/#it{c}^{2}");
	anaLbMuMu_data->SetUnits("MeV/#it{c}^{2}");

	curcut.Print();
	anaLbMuMu->SetSignal((model+"-s[15]-s2[30]").c_str(),1e3,"-namepar");
	histFile->cd();
	anaLbMuMu->Initialize("-docuts-namepar");
	anaLbMuMu->Fit(5400.,5750.,200,true,"-quiet-noleg-lin-range-xuM(#Lambda#mu#mu)-noCost-layout[0.7,0.5,0.9,0.9]");
	Str2VarMap MCmumupars = anaLbMuMu->GetSigParams();
	
	if(type=="LL") yield_rare_LL_MC = anaLbMuMu->GetNSigVal();
	else if (type=="DD") yield_rare_DD_MC = anaLbMuMu->GetNSigVal();

	rarePars["n"] = MCmumupars["n"];
	rarePars["a"] = MCmumupars["a"];
	rarePars["a2"] = MCmumupars["a2"];
	setConstant(&rarePars);
	RooRealVar * factor = new RooRealVar("factor_"+type+q2name,"c",MCmumupars["s"]->getVal()/MCpars["s"]->getVal());
	ModifyPars(&rarePars,"s",factor);
	ModifyPars(&rarePars,"s2",factor);

	if(q2min==1.1 && q2max==6) anaLbMuMu_data->SetName("Lb2Lmumu_"+type+"_lowQ2");
	else if(q2min==15 && q2max==20) anaLbMuMu_data->SetName("Lb2Lmumu_"+type+"_highQ2");

	TString KSname = "KS_bkg_"+type+"_"+q2name;
	Analysis * KSmumu = new Analysis(KSname,"Lb","candBdKSmumu",candfilename,vLbCons,&curcut);
	double nKSmumu_MC_all = KSmumu->GetReducedTree()->GetEntries();		
	KSmumu->SetSignal("DCB_OST",1.e4,"-namepar");
	KSmumu->Initialize("-docuts-namepar");
	double nKSmumu_MC = KSmumu->GetReducedTree()->GetEntries();

	double Jpsi2mumuBr = 0.0593;
	double BRKS_mumuOverJpsi = 3.4e-7/8.73e-4;
	double nKSmumu_val = KS_double_eff_ratio*(BRKS_mumuOverJpsi/Jpsi2mumuBr)*nKSjpsi_val*nKSmumu_MC/nKSmumu_MC_all;
	RooRealVar * nKSmumu = new RooRealVar("nKSmumu_"+type,"N_{K_S\\mu\\mu}",nKSmumu_val,0.,50.);
	nKSmumu->setConstant(true);

	if(type == "DD") KSpredict_file << q2min << "-" << q2max << "   &   " << nKSmumu_val << "  &   ";
	else KSpredict_file << nKSmumu_val << "\\\\" << endl << flush; 

	double curq2 = (q2max+q2min)/2.;
	string getq2 = Form("%4.2f",(float)curq2);
	double cureff = eff[getq2][0];
	if(type == "LL") cureff *= pideff[curq2];
	double njpsi = Njpsi->getVal();
	double ww = 0;
	if(wilson && type == "LL") { ww = ( w1LL[getq2][0] + w2LL[getq2][0] + w3LL[getq2][0] ) /3.; cout << "Including WC" << endl; } 
	if(wilson && type == "DD") { ww = ( w1DD[getq2][0] + w2DD[getq2][0] + w3DD[getq2][0] ) /3.; cout << "Including WC" << endl; } 

	double eff_err = TMath::Sqrt( TMath::Power(eff[getq2][1],2) + TMath::Power(cureff*ww,2) ); 
	
	double percErr = eff_err / cureff;
	if(sys == "norm_plus") njpsi += Njpsi->getError();
	else if(sys == "norm_minus") njpsi -= Njpsi->getError();
	else if(sys == "eff_plus") cureff +=  eff_err;
	else if(sys == "eff_minus") cureff -= eff_err;
	effs_file << eff[getq2][0] << "  &  " << eff_err << "  &  pce " << percErr << "  &  ";

	RooAbsReal * Nsig = NULL;
	TString formula = Form("BR * %e * %e * (%e - %e)/%e",cureff,njpsi,q2max,q2min,Jpsi2mumuBr);
	if(dosig=="BR") Nsig = new RooFormulaVar("Nsig_"+type,formula,RooArgSet(*BR));
	else if (dosig=="Ncand") Nsig = new RooRealVar("Nsig_"+type,"Nsig_"+type,100.,5.,7.e2);
	else Nsig = new RooRealVar("Nsig_"+type,"Nsig_"+type,0.,0.,0.);
	//Nsig->Print();
	if(type=="LL") { factorLL = cureff * njpsi * (q2max - q2min) / Jpsi2mumuBr; }
	else { factorDD = cureff * njpsi * (q2max - q2min) / Jpsi2mumuBr; }

	anaLbMuMu_data->SetSignal(model.c_str(),Nsig,"-namepar",rarePars);
	anaLbMuMu_data->addBkgComponent("Comb","Exp-b[-0.003,-0.01,0.01]",20.,"-namepar");
	anaLbMuMu_data->addBkgComponent("KSmumu",KSmumu->GetReducedTree(),nKSmumu,"-namepar");

	//if(q2max == 8.) anaLbMuMu_data->addBkgComponent("Jpsi_tail",JpsiTailTree,nJpsi_tail);
		
	anaLbMuMu_data->Initialize("-docuts-namepar");
	
	return anaLbMuMu_data; 
}




int main(int argc, char **argv)
{
	string analysis = "Lb2Lmumu";
	string mother = "";
	bool useSingle = false;
	string type = "all";
	string dosig = "BR";
	string part = "";
	bool doHighq2Only = false;
	string polarity = "";
	string year = "";

	gROOT->ProcessLine(".x ~/work/lhcbStyle.C");

	if(argc > 1)
	{
		for(int a = 1; a < argc; a++)
		{
			string arg = argv[a];
			string str = arg.substr(2,arg.length()-2);

			if(arg.find("-m") != string::npos) mother = str;
			if(arg.find("-p") != string::npos) polarity = str;
			if(arg.find("-y") != string::npos) year = str;
			if(arg == "-s") useSingle = true;
			if(arg == "-W") wilson = true;
			if(arg == "-nosig") dosig = "";
			if(arg == "-useNcand") dosig = "Ncand";
			if(arg == "-highQ2only") doHighq2Only = true;
			if(arg.find("-S") != string::npos) sys = str;
			if(arg.find("-t") != string::npos) type = str;
		}
	}

	cout << "Analysing " << type << endl;
	if(wilson) cout << "Using Wilson" << endl;
	if(sys!="") cout << sys << endl;

	Analysis::SetPrintLevel("v");

	pideff[17.5] = 1.00281;
	pideff[11.75] = 1.00151;
	pideff[15.5] = 1.00431;
	pideff[17.0] = 1.00215;
	pideff[19.0] = 1.00226;
	pideff[3.55] = 0.99589;
	pideff[1.05] = 0.99418;
	pideff[3.0] = 0.99523;
	pideff[5.0] = 0.99699;
	pideff[7.0] = 0.99805;

	pdfsys[17.5] = 0.01;
	pdfsys[11.75] = 0.032;
	pdfsys[15.5] = 0.028;
	pdfsys[17.0] = 0.014;
	pdfsys[19.0] = 0.025;
	pdfsys[3.55] = 0.025;
	pdfsys[1.05] = 0.025;
	pdfsys[3.0] = 0.025;
	pdfsys[5.0] = 0.025;
	pdfsys[7.0] = 0.025;

	map<string, vector<double> > effLL = mergeGr("toteff", "/afs/cern.ch/work/p/pluca/results/LbreleffAndSysvsq2_LL_7bins.root", "/afs/cern.ch/work/p/pluca/results/LbreleffAndSysvsq2_LL_2bins.root");
	map<string, vector<double> > effDD = mergeGr("toteff", "/afs/cern.ch/work/p/pluca/results/LbreleffAndSysvsq2_DD_7bins.root", "/afs/cern.ch/work/p/pluca/results/LbreleffAndSysvsq2_DD_2bins.root");
	map<string, vector<double> > effLL_low = mergeGr("toteff_lowSel", "/afs/cern.ch/work/p/pluca/results/LbreleffAndSysvsq2_LL_7bins.root", "/afs/cern.ch/work/p/pluca/results/LbreleffAndSysvsq2_LL_2bins.root");
	map<string, vector<double> > effDD_low = mergeGr("toteff_lowSel", "/afs/cern.ch/work/p/pluca/results/LbreleffAndSysvsq2_DD_7bins.root", "/afs/cern.ch/work/p/pluca/results/LbreleffAndSysvsq2_DD_2bins.root");

	w1DD = mergeGr("rel_wilson_sys1", "/afs/cern.ch/work/p/pluca/results/LbreleffAndSysvsq2_DD_7bins.root", "/afs/cern.ch/work/p/pluca/results/LbreleffAndSysvsq2_DD_2bins.root",true);
	w2DD = mergeGr("rel_wilson_sys2", "/afs/cern.ch/work/p/pluca/results/LbreleffAndSysvsq2_DD_7bins.root", "/afs/cern.ch/work/p/pluca/results/LbreleffAndSysvsq2_DD_2bins.root",true);
	w3DD = mergeGr("rel_wilson_sys3", "/afs/cern.ch/work/p/pluca/results/LbreleffAndSysvsq2_DD_7bins.root", "/afs/cern.ch/work/p/pluca/results/LbreleffAndSysvsq2_DD_2bins.root",true);
	w1LL = mergeGr("rel_wilson_sys1", "/afs/cern.ch/work/p/pluca/results/LbreleffAndSysvsq2_LL_7bins.root", "/afs/cern.ch/work/p/pluca/results/LbreleffAndSysvsq2_LL_2bins.root",true);
	w2LL = mergeGr("rel_wilson_sys2", "/afs/cern.ch/work/p/pluca/results/LbreleffAndSysvsq2_LL_7bins.root", "/afs/cern.ch/work/p/pluca/results/LbreleffAndSysvsq2_LL_2bins.root",true);
	w3LL = mergeGr("rel_wilson_sys3", "/afs/cern.ch/work/p/pluca/results/LbreleffAndSysvsq2_LL_7bins.root", "/afs/cern.ch/work/p/pluca/results/LbreleffAndSysvsq2_LL_2bins.root",true);

	TCut baseCut = "";
	if(useSingle) { singleCut = "isChosenCand"; tname = "singleCand_"; cout << "Using single candidates" << endl; }
	if(mother=="Lb") baseCut+= "pplus_ID>0";
	else if (mother=="Lbbar") baseCut+="pplus_ID<0";
	if(polarity!="")
	{
		//baseCut+=TCut(("Polarity=="+polarity).c_str());
		yCut = TCut(("Polarity=="+polarity).c_str());
	}
	if(year!="")
	{
		//baseCut+=TCut(("dataType=="+year).c_str());
		pCut = TCut(("dataType=="+year).c_str());
	}


	effs_file.open("effs.txt");
	KSpredict_file.open("KSpredictons.txt");
	ofstream yields_file, sum_yields;
	yields_file.open("yields.txt");
	sum_yields.open("sum_yields.txt");
	ofstream params_file;
	params_file.open ("params.txt");
	histFile = new TFile("Lbyield.root","recreate");

	// FIT Lb data

	RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

	RooRealVar * NJpsiLL, * NJpsiDD;
	Str2VarMap MCparsLL, MCparsDD, jpsiParsHighLL, jpsiParsHighDD;
	double KS_double_eff_ratio_LL, nKSjpsi_val_LL, KS_double_eff_ratio_DD, nKSjpsi_val_DD;
	Analysis * anaJpsiDD = fitJpsi("DD", baseCut+CutsDef::DDcut, &NJpsiDD, &MCparsDD, &KS_double_eff_ratio_DD, &nKSjpsi_val_DD, false);
	Analysis * anaJpsiLL = fitJpsi("LL", baseCut+CutsDef::LLcut, &NJpsiLL, &MCparsLL, &KS_double_eff_ratio_LL, &nKSjpsi_val_LL, false);
	Str2VarMap jpsiParsDD = anaJpsiDD->GetSigParams();
	Str2VarMap jpsiParsLL = anaJpsiLL->GetSigParams();
	jpsiParsHighLL = jpsiParsLL;
	jpsiParsHighDD = jpsiParsDD;

	double yield_rare_LL = 0, yield_rare_DD = 0;
	double yield_rare_LLerr = 0, yield_rare_DDerr = 0;
	double yield_jpsi_LL = NJpsiLL->getVal(), yield_jpsi_DD = NJpsiDD->getVal();
	double yield_jpsi_LLerr = NJpsiLL->getError(), yield_jpsi_DDerr = NJpsiDD->getError();

	vector<double> chi2DD(1,anaJpsiDD->GetProb());
	vector<double> chi2LL(1,anaJpsiLL->GetProb());

	cout << "---------- PARAMETERS ------------"<< endl;
	PrintPars(jpsiParsLL);
	cout << "------------------"<< endl;
	PrintPars(jpsiParsDD);
	cout << "------------------"<< endl;

	int nbins = CutsDef::nq2bins;
	double * q2min = &CutsDef::q2min_highfirst[0];
	double * q2max = &CutsDef::q2max_highfirst[0];

	TCanvas *cDD = new TCanvas("cDD","");
	TCanvas *c2DD = new TCanvas("c2DD","");
	TCanvas *cLL = new TCanvas("cLL","");
	TCanvas *c2LL = new TCanvas("cLL2","");
	TCanvas *ctemp = new TCanvas("ctemp","");
	cDD->Divide(2,2);
	c2DD->Divide(2,2);
	cLL->Divide(2,2);
	c2LL->Divide(2,2);
	
	TCanvas *cAll_horiz = new TCanvas("cAll_horiz","",800,400);
	cAll_horiz->Divide(4,2);
	TCanvas *cAll = new TCanvas("cAll","",550,800);
	cAll->Divide(2,4);

	TPaveText * tbox = new TPaveText(gStyle->GetPadRightMargin() + 0.63,
                   0.80 - gStyle->GetPadTopMargin(),
                   gStyle->GetPadRightMargin() + 0.83,
                   0.97 - gStyle->GetPadTopMargin(),
                   "BRNDC");
	tbox->AddText("LHCb");
	tbox->SetFillStyle(0);
	tbox->SetTextAlign(12);
	tbox->SetBorderSize(0);

	MultiAnalysis * anaJpsi = new MultiAnalysis("fitJpsi");
	anaJpsi->AddCategory(anaJpsiLL,"LL");
	anaJpsi->AddCategory(anaJpsiDD,"DD");
	ctemp->cd();
	RooPlot * jpsiplot = anaJpsi->PrintSum("-xuM(#Lambda#mu#mu)-nbins75-noParams-nototsigplot-stacked");
	jpsiplot->Draw();
	ctemp->Print("fit_Jpsi_All.pdf");
	ctemp->Print("fit_Jpsi_All.eps");
	ctemp->Print("fit_Jpsi_All.png");
	ctemp->Print("fit_Jpsi_All.C");
	ctemp->SetLogy();
	jpsiplot->addObject(tbox);
	jpsiplot->SetMinimum(10-(1e-2));
	//jpsiplot->SetMinimum(0.+(1e-2));
	jpsiplot->SetMaximum(1e4);
	jpsiplot->Draw();
	ctemp->Print("fit_Jpsi_All_log.pdf");
	ctemp->Print("fit_Jpsi_All_log.eps");
	ctemp->Print("fit_Jpsi_All_log.png");
	ctemp->Print("fit_Jpsi_All_log.C");

	ctemp->SetLogy(0);

	TGraphAsymmErrors * gr_result = new TGraphAsymmErrors();
	int p = 0;

	for(int i = 0; i < nbins; i++)
	{
		if(doHighq2Only && i > 0) break;
		if(q2min[i] == 0.1) 
		{
			Analysis *anaJpsiLL, *anaJpsiDD;
			anaJpsiDD = fitJpsi("DD", baseCut+CutsDef::DDcut, &NJpsiDD, &MCparsDD, &KS_double_eff_ratio_DD, &nKSjpsi_val_DD, true);
			anaJpsiLL = fitJpsi("LL", baseCut+CutsDef::LLcut, &NJpsiLL, &MCparsLL, &KS_double_eff_ratio_LL, &nKSjpsi_val_LL, true);
			jpsiParsDD = anaJpsiDD->GetSigParams();
			jpsiParsLL = anaJpsiLL->GetSigParams();
			effLL = effLL_low;
			effDD = effDD_low;
			MultiAnalysis * anaJpsi = new MultiAnalysis("fitJpsi_lowSel");
			anaJpsi->AddCategory(anaJpsiLL,"LL");
			anaJpsi->AddCategory(anaJpsiDD,"DD");
			ctemp->cd();
			anaJpsi->PrintSum("-xuM(#Lambda#mu#mu)-nbins75-noParams-LHCb-stacked")->Draw();
			ctemp->SetLogy();
			ctemp->Print("fit_Jpsi_All_lowSel.pdf");
			ctemp->SetLogy(0);
			chi2DD.push_back(anaJpsiDD->GetProb());
			chi2LL.push_back(anaJpsiLL->GetProb());
		}

		RooRealVar * BR = new RooRealVar("BR","BR_{\\Lambda\\mu\\mu}",1e-4,0.,1e-1);
		RooRealVar * vLbCons = new RooRealVar("Lb_MassConsLambda","Lb_MassConsLambda",5621.,5400.,6000.);
		TCut curq2cut = (TCut)Form("TMath::Power(J_psi_1S_MM/1000,2) >= %.1f  && TMath::Power(J_psi_1S_MM/1000,2) < %.1f",q2min[i],q2max[i]);
		TCut curcut = curq2cut+baseCut;
		curcut.Print();

		effs_file << fixed << setprecision(1) << q2min[i] << "-" << q2max[i] << setprecision(3) << "  &  ";

		Analysis * anaDD = getAnaObject("DD", curcut+CutsDef::DDcut, q2min[i], q2max[i], vLbCons, BR, effDD, NJpsiDD, jpsiParsDD, MCparsDD, KS_double_eff_ratio_DD, nKSjpsi_val_DD, dosig);
		Analysis * anaLL = getAnaObject("LL", curcut+CutsDef::LLcut, q2min[i], q2max[i], vLbCons, BR, effLL, NJpsiLL, jpsiParsLL, MCparsLL, KS_double_eff_ratio_LL, nKSjpsi_val_LL, dosig);

		float curq2 = (q2max[i]+q2min[i])/2.;
		effs_file << TMath::Sqrt( TMath::Power(pdfsys[curq2],2) + TMath::Power(0.0006/0.0593,2) )*100 << "\\%    \\\\" << endl << flush; 

		MultiAnalysis * ana = new MultiAnalysis("analysis");
		if(type == "LL" || type == "all") ana->AddCategory(anaLL, "LL");
		if(type == "DD" || type == "all") ana->AddCategory(anaDD, "DD");
		ana->Initialize();

		string options = "-nototsigplot-lin-range-minos-xuM(#Lambda#mu#mu)-noCost-noParams-cutLegNames";
		//-legf-leg[0.7,0.5,0.9,0.9]-cutLegNames";
		map<string, RooPlot * > plots = ana->SimultaneousFit(5400.,6000.,30,options);
		RooFitResult * fitResult = ana->GetFitResult();
		double LogL = ana->CreateLogL()->getVal();

		double errX = (q2max[i]-q2min[i])/2.;
		gr_result->SetPoint(p,(q2max[i]+q2min[i])/2.,BR->getVal());
		gr_result->SetPointError(p,errX,errX,BR->getErrorLo(),BR->getErrorHi());
		p++;

		double LLnn = anaLL->GetSigVal();
		double DDnn = anaDD->GetSigVal();
		double LLerr = anaLL->GetNSigPtr()->getPropagatedError(*fitResult);
		double DDerr = anaDD->GetNSigPtr()->getPropagatedError(*fitResult);
		//double LLerr_low = factorLL * BR->getErrorLo();
		//double LLerr = factorLL * BR->getErrorHi();
		//double DDerr_low = factorDD * BR->getErrorLo();
		//double DDerr = factorDD * BR->getErrorHi();
		if(i==0) { yield_rare_LL = LLnn; yield_rare_DD = DDnn; yield_rare_LLerr = LLerr; yield_rare_DDerr = DDerr; }
		yields_file << fixed << setprecision(1) << q2min[i] << "-" << q2max[i] << "  &  $" << LLnn << " \\pm " << LLerr << "$  &  $" << DDnn << " \\pm " << DDerr << "\\\\ " << endl << flush;
		sum_yields << fixed << setprecision(1) << q2min[i] << "-" << q2max[i] << "  &  $" << LLnn + DDnn << " \\pm ";
		//sum_yields << "^{+" << LLerr + DDerr << "}_{-";
		//sum_yields << LLerr_low + DDerr_low << "} $ & " << LogL << "\\\\ " << endl << flush;
		//sum_yields << TMath::Sqrt(TMath::Power(LLerr,2) + TMath::Power(DDerr,2)) << "$ & " << LogL << "\\\\ " << endl << flush;
		sum_yields << LLerr + DDerr << "$ & " << LogL << "\\\\ " << endl << flush;

		if((q2min[i]==1.1 && q2max[i]==6) || (q2min[i]==15 && q2max[i]==20)) ctemp->cd();
		else if(i <= 4) cDD->cd(i);
		else c2DD->cd(i-4);
		if(type == "DD" || type == "all") plots["DD"]->SetTitle((TString)Form("[%.1f,%.1f] GeV^{2} (DD events) ",q2min[i],q2max[i]));
		if(type == "DD" || type == "all") plots["DD"]->Draw();
		
		if(q2min[i]==1.1 && q2max[i]==6) ctemp->Print("fit_DD_lowQ2.pdf");
		else if(q2min[i]==15 && q2max[i]==20) ctemp->Print("fit_DD_highQ2.pdf");

		if((q2min[i]==1.1 && q2max[i]==6) || (q2min[i]==15 && q2max[i]==20)) ctemp->cd();
		else if(i <= 4) cLL->cd(i);
		else c2LL->cd(i-4);
		if(type == "LL" || type == "all") plots["LL"]->SetTitle((TString)Form("[%.1f,%.1f] GeV^{2} (LL events) ",q2min[i],q2max[i]));
		if(type == "LL" || type == "all") plots["LL"]->Draw();
	
		if(q2min[i]==1.1 && q2max[i]==6) ctemp->Print("fit_LL_lowQ2.pdf");
		else if(q2min[i]==15 && q2max[i]==20) ctemp->Print("fit_LL_highQ2.pdf");
		
		if((q2min[i]==1.1 && q2max[i]==6) || (q2min[i]==15 && q2max[i]==20)) ctemp->cd();
		else if(i <= 4) cAll_horiz->cd(i+4);
		else cAll_horiz->cd(i-4);
		
		RooPlot * pl = ana->PrintSum("-xuM(#Lambda#mu#mu)-nbins20-noParams-LHCbDX-nototsigplot");
		pl->SetMinimum(1e-3);
		//pl->SetMaximum(90.);

		if(!((q2min[i]==1.1 && q2max[i]==6) || (q2min[i]==15 && q2max[i]==20)))
		{
			//char lett = 97+i-1+4;
			//if(i > 4) lett = 97+i-1-4;
			TPaveText * text = new TPaveText(0.5,0.65,0.9,0.8,"BRNDC");
			//text->AddText("("+(TString)lett+")");
			//text->AddText(Form("%.1f < #it{q}^{2} < %.1f GeV^{2}/#it{c}^{4}",q2min[i],q2max[i]));	
			text->AddText(Form("[%.1f,%.1f] GeV^{2}/#it{c}^{4}",q2min[i],q2max[i]));
			text->SetBorderSize(0);
			text->SetFillStyle(0);
			pl->addObject(text);
		}
		pl->Draw();

		if(i <= 4) cAll->cd(i+4);
		else cAll->cd(i-4);
		if( !( (q2min[i]==1.1 && q2max[i]==6) || (q2min[i]==15 && q2max[i]==20) ) ) pl->Draw();

		if(q2min[i]==1.1 && q2max[i]==6) 
		{
			ctemp->Print("fit_All_lowQ2.pdf");
			ctemp->Print("fit_All_lowQ2.eps");
			ctemp->Print("fit_All_lowQ2.png");
			ctemp->Print("fit_All_lowQ2.C");
		}
		if(q2min[i]==15 && q2max[i]==20)
		{
			ctemp->Print("fit_All_highQ2.pdf");
			ctemp->Print("fit_All_highQ2.eps");
			ctemp->Print("fit_All_highQ2.png");
			ctemp->Print("fit_All_highQ2.C");
		}
		ctemp->cd();

		((RooRealVar *) anaDD->GetSigParams()["factor"])->setConstant(false);
		((RooRealVar *) anaLL->GetSigParams()["factor"])->setConstant(false);
		BR->setConstant(true);
		freopen ("params.txt", "a", stdout);
		cout << fixed << setprecision(2) << setprecision(2) << "\\multicolumn{5}{ " << q2min[i] << "-" << q2max[i] << "}  \\\\" << endl; 
		cout << "DD events" << endl;
		anaDD->PrintParamsTable("-nocost");
		cout << "LL events" << endl;
		anaLL->PrintParamsTable("-nocost");
		freopen ("/dev/tty", "a", stdout);
	}

	PrintPars(jpsiParsLL);
	PrintPars(jpsiParsDD);
	PrintPars(jpsiParsHighLL);
	PrintPars(jpsiParsHighDD);

	cDD->Print("q2_fits_DD_plot1.pdf");
	c2DD->Print("q2_fits_DD_plot2.pdf");
	cLL->Print("q2_fits_LL_plot1.pdf");
	c2LL->Print("q2_fits_LL_plot2.pdf");

	cAll->Print("q2_fits_All_vertical.pdf");
	cAll->Print("q2_fits_All_vertical.eps");
	cAll->Print("q2_fits_All_vertical.png");
	cAll->Print("q2_fits_All_vertical.C");
	cAll_horiz->Print("q2_fits_All_horizontal.pdf");
	cAll_horiz->Print("q2_fits_All_horizontal.eps");
	cAll_horiz->Print("q2_fits_All_horizontal.png");
	cAll_horiz->Print("q2_fits_All_horizontal.C");

	//cout << "Chi DD ->  high : " << chi2DD[0] << ",  low : " << chi2DD[1] << endl;
	//cout << "Chi LL ->  high : " << chi2LL[0] << ",  low : " << chi2LL[1] << endl;

	ofstream ratio_file;
	if (polarity=="1")
		{	ratio_file.open("log_"+mother+year+"magUp.txt"); }
	else if (polarity=="-1")
		{ ratio_file.open("log_"+mother+year+"magDown.txt"); }
	else
		{ ratio_file.open("log_"+mother+year+".txt"); }

	ratio_file << "  LL \t " << yield_rare_LL << "  " << yield_jpsi_LL << "  " << yield_rare_LL / yield_jpsi_LL << " +/- ";
	ratio_file << (yield_rare_LL / yield_jpsi_LL)*TMath::Sqrt(TMath::Power(yield_rare_LLerr/yield_rare_LL,2) + TMath::Power(yield_jpsi_LLerr/yield_jpsi_LL,2)) << endl;
	ratio_file << "  DD \t " << yield_rare_DD << "  " << yield_jpsi_DD << "  " << yield_rare_DD / yield_jpsi_DD << " +/- ";
	ratio_file << (yield_rare_DD / yield_jpsi_DD)*TMath::Sqrt(TMath::Power(yield_rare_DDerr/yield_rare_DD,2) + TMath::Power(yield_jpsi_DDerr/yield_jpsi_DD,2)) << endl;
	
	//cout << "  LL  " << yield_rare_LL << "  " << yield_jpsi_LL << "  " << yield_rare_LL / yield_jpsi_LL << " +/- ";
	//cout << (yield_rare_LL / yield_jpsi_LL)*TMath::Sqrt(TMath::Power(yield_rare_LLerr/yield_rare_LL,2) + TMath::Power(yield_jpsi_LLerr/yield_jpsi_LL,2)) << endl;
	//cout << "  DD  " << yield_rare_DD << "  " << yield_jpsi_DD << "  " << yield_rare_DD / yield_jpsi_DD << " +/- ";
	//cout << (yield_rare_DD / yield_jpsi_DD)*TMath::Sqrt(TMath::Power(yield_rare_DDerr/yield_rare_DD,2) + TMath::Power(yield_jpsi_DDerr/yield_jpsi_DD,2)) << endl;



	//cout << "rare/jpsi LL" << yield_rare_LL_MC << "  " << yield_jpsi_LL_MC << "  " << yield_rare_LL_MC / yield_jpsi_LL_MC << " +/- " << (yield_rare_LL_MC / yield_jpsi_LL_MC)*TMath::Sqrt(1./yield_rare_LL_MC + 1./yield_jpsi_LL_MC) << endl;
	//cout << "rare/jpsi DD" << yield_rare_DD_MC << "  " << yield_jpsi_DD_MC << "  " << yield_rare_DD_MC / yield_jpsi_DD_MC << " +/- " << (yield_rare_DD_MC / yield_jpsi_DD_MC)*TMath::Sqrt(1./yield_rare_DD_MC + 1./yield_jpsi_DD_MC) << endl;

	histFile->Close();
	histFile = new TFile("Lbyield.root","recreate");
	gr_result->Write("BRplot");
	histFile->Write();
	histFile->Close();
	delete histFile;

	return 0;
}





