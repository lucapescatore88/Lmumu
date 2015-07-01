#include <vector>
#include <sstream>
#include <iostream>
#include <string>
#include <time.h>

#include "RooRealVar.h"
#include "RooStats/SPlot.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooCBShape.h"
#include "RooArgusBG.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooDataSet.h"
#include "RooAbsData.h"
#include "RooProdPdf.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "RooArgusBG.h"
#include "RooKeysPdf.h"
#include "RooArgList.h"
#include "RooPolynomial.h"

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

#include "TEntryList.h"
#include "TCanvas.h"
#include "TBranch.h"
#include "TCut.h"
#include "TGraphErrors.h"

#include "general_functions.hpp"
#include "ReadTree_comp.hpp"
#include "analyser.hpp"
#include "Lb_cuts.hpp"

using namespace std;
using namespace RooFit;


vector< double > printSys(vector<string> cuts, vector<double> ndef, vector<double> errv, vector<double> nmsigma = vector<double>(), vector<double> npsigma = vector<double>(), string opt = "All", bool stats = false, string partialsys = "")
{
	vector <double >res;
	double totaleff = 1, totalerr = 0, sys = 0;
	static vector< double > totsys, partialsumsys;
		
	for(unsigned i = 1; i < cuts.size(); i++)
	{
		double def = ndef[i-1];
        double err = errv[i-1];
		totaleff *= def;
		totalerr = TMath::Sqrt( TMath::Power(totalerr,2) + TMath::Power(err/def,2) );
		double minus = def;
		double plus = def;
		if(nmsigma.size()!=0) minus = nmsigma[i-1];
		if(npsigma.size()!=0) plus = npsigma[i-1];
		sys = TMath::Max(TMath::Abs(def-plus),TMath::Abs(def-minus)) / def;
		if(totsys.size() < cuts.size()-1) totsys.push_back(err/def);
		if(totsys.size() < cuts.size()-1) partialsumsys.push_back(err/def);
        if(partialsys=="-rpartial") partialsumsys[i-1] = 0;
		if(partialsys=="-rtotal") totsys[i-1] = err/def;
		partialsumsys[i-1] = TMath::Sqrt( partialsumsys[i-1]*partialsumsys[i-1] + sys*sys );
		totsys[i-1] = TMath::Sqrt( totsys[i-1]*totsys[i-1] + sys*sys );
		
		if( opt.find("All") == string::npos && opt.find(cuts[i]) == string::npos ) continue;

		if(opt.find("None") == string::npos)
		{
			if(opt == "All")
			{
				cout << cuts[i] << "\t " << def << " +/- " << err << "\t ";
				if(partialsys.find("partial") == string::npos)
				{
					if(minus != def) cout << minus << "\t ";
					else cout << "-------\t ";
					if(plus != def) cout << plus << "\t ";
					else cout << "-------\t ";
					if(minus != def) cout << sys*100 << "% \\\\" << endl;
					else cout << "------- \\\\" << endl;
				}
				else cout << "-------\t -------\t " << partialsumsys[i-1]*100 << "% \\\\" << endl;
			}
			else if(opt.find("eff") != string::npos) cout << "\t& $" << def << " \\pm " << err << "$";
			else
			{
				if(stats) cout << "\t& " << err;
				else if(partialsys.find("-ppartial") != string::npos) cout << "\t& " << partialsumsys[i-1]*def;
				else cout << "\t& " << sys*def;
				
				if(partialsys.find("-ptotal") != string::npos) cout << "\t& " << totsys[i-1]*def;
			}
		}
		
		res.push_back(def);
		res.push_back(err/def);
		res.push_back(sys/def);
	}
	
	if(opt.find("effAll") != string::npos ) cout  << "\t & " << totaleff << " \\pm " << totalerr*totaleff << " \\\\" << endl;
	else if(opt.find("eff") != string::npos ) cout << " \\\\" << endl;

	if(partialsys.find("-ototal") != string::npos) return totsys;	
	return res;
}

vector<double > printSys(vector<string> cuts, vector<double> ndef, vector<double> errv, string opt = "All", bool stats = false, string partialsys = "")
{
	vector<double> v;
	return printSys(cuts, ndef, errv, v, v, opt, stats, partialsys);
}


double toteff(vector<double> effs)
{
	double tot = 1;
	for(unsigned i = 0; i < effs.size(); i++) tot*=effs[i];
	return tot;
}

double toterr(vector<double> effs, vector<double> errs)
{
	double tot = 0;
	if(errs.size()!=effs.size()) {cout << "Eff size different than err size" << endl; return -1;}
	for(unsigned i = 0; i < errs.size(); i++) tot+=TMath::Power(errs[i]/effs[i],2);
	return toteff(effs)*TMath::Sqrt(tot);
}



int main(int argc, char **argv)
{
	bool doPDFsys = true;
	bool doEff = true;
	bool doEffsys = true;
	bool pythia6 = false;

	string analysis = "Lb2Lmumu";
	string type = "All";
	string print = "All";
	bool rel = false;
	if(argc > 1)
	{
		for(int a = 1; a < argc; a++)
		{
			string arg = argv[a];
			string str = arg.substr(2,arg.length()-2);
			
			if(arg.find("-t") != string::npos)
			{
				type = str;
				analysis += ("_" + type);
			}
			if(arg == "-pdf") doEff = false;
			if(arg == "-eff") { doPDFsys = false; doEffsys = false;}
			if(arg == "-syseff") doPDFsys = false;
			if(arg.find("-p") != string::npos) print = str;
			if(arg == "-r") rel = true;
			if(arg == "-Pythia6") pythia6 = true;
		}
	}
	
	
	// Set trees of cancdidates previously created

	TString mctype = "MC_Pythia8";
	if(pythia6) mctype = "MC_Pythia6";

	TString weight = "MCnorm*(lifeTimeW > 0)*(physRate_pol0_noDecay > 0)";
	
	cout << "Using: " << mctype << endl;

	TString namefileMCgeom = "/afs/cern.ch/work/p/pluca/weighted/Lmumu/Lb2Lmumu_geom"+mctype+"_NBweighted.root";
	TFile * MCfile = TFile::Open(namefileMCgeom);
	TTree * treeMCgeom = (TTree *)MCfile->Get("MCtree");
	TString namefileMC = "/afs/cern.ch/work/p/pluca/weighted/Lmumu/Lb2Lmumu_"+mctype+"_NBweighted.root";
	MCfile = TFile::Open(namefileMC);
	TTree * treeMC = (TTree *)MCfile->Get("tree");
	TTree * treeMCGen = (TTree *)MCfile->Get("MCtreeDecay");
	TTree * treeMCAllGen = (TTree *)MCfile->Get("MCtree");
	TString namefileMCjpsi = "/afs/cern.ch/work/p/pluca/weighted/Lmumu/Lb2JpsiL_"+mctype+"_NBweighted.root";
	MCfile = TFile::Open(namefileMCjpsi);
	TTree * treeMCjpsi = NULL, * treeMCjpsi_Gen = NULL, * treeMCjpsi_AllGen = NULL, * treeMCjpsi_geom = NULL;
	if(rel)
	{
		treeMCjpsi = (TTree *)MCfile->Get("tree");
		treeMCjpsi_Gen = (TTree *)MCfile->Get("MCtreeDecay");
		treeMCjpsi_AllGen = (TTree *)MCfile->Get("MCtree");
	}
	TString namefileMCjpsiGeom = "/afs/cern.ch/work/p/pluca/weighted/Lmumu/Lb2JpsiL_geom"+mctype+"_NBweighted.root";
	MCfile = TFile::Open(namefileMCjpsiGeom);
	if(rel) treeMCjpsi_geom = (TTree *)MCfile->Get("MCtree");

	
	
	/**     PDF systematic      */
	
	if(doPDFsys)
	{
		TString datafilename = "/afs/cern.ch/user/p/pluca/work/Lb/Lmumu/candLb.root";
		if(TFile::Open(datafilename)==0) return 1;
	
		TreeReader * candTree_LbMuMu_data = new TreeReader("candLb2Lmumu");
		candTree_LbMuMu_data->AddFile(datafilename);
		RooRealVar * vMuMu = new RooRealVar("Lb_MassConsLambda","Lb_MassConsLambda",5620.,5250.,6000.);
		Analysis * anaLbMuMu_default = new Analysis("Lmumu_default","Lb",candTree_LbMuMu_data,&CutsDef::cutMuMu,vMuMu);
		Analysis * anaLbMuMu_Gauss = new Analysis("mumu_Gauss","Lb",candTree_LbMuMu_data,&CutsDef::cutMuMu,vMuMu);
		Analysis * anaLbMuMu_DGauss = new Analysis("mumu_DGauss","Lb",candTree_LbMuMu_data,&CutsDef::cutMuMu,vMuMu);
		Analysis * anaLbMuMu_CB = new Analysis("mumu_CB","Lb",candTree_LbMuMu_data,&CutsDef::cutMuMu,vMuMu);
		Analysis * anaLbMuMu_CBGauss = new Analysis("mumu_CBGauss","Lb",candTree_LbMuMu_data,&CutsDef::cutMuMu,vMuMu);

		TreeReader * candTree_LbJpsi_data = new TreeReader("candLb2JpsiL");
		candTree_LbJpsi_data->AddFile(datafilename);
		RooRealVar * vJpsi = new RooRealVar("Lb_MassConsJpsiLambda","Lb_MassConsLambda",5620.,5250.,6000.);
		Analysis * anaLbJpsi_default = new Analysis("Jpsi_default","Lb",candTree_LbJpsi_data,&CutsDef::cutJpsi,vJpsi);
		Analysis * anaLbJpsi_CBGauss = new Analysis("Jpsi_CBGauss","Lb",candTree_LbJpsi_data,&CutsDef::cutJpsi,vJpsi);
		Analysis * anaLbJpsi_DGauss = new Analysis("Jpsi_DGauss","Lb",candTree_LbJpsi_data,&CutsDef::cutJpsi,vJpsi);
		Analysis * anaLbJpsi_CB = new Analysis("Jpsi_CB","Lb",candTree_LbJpsi_data,&CutsDef::cutJpsi,vJpsi);
		Analysis * anaLbJpsi_ChebKS = new Analysis("Jpsi_KS","Lb",candTree_LbJpsi_data,&CutsDef::cutJpsi,vJpsi);

		TFile * MCFile = TFile::Open("candLb_MC.root");
		//TTree * LbJpsiL_MCTree = (TTree *)MCFile->Get("candLb2JpsiL");
		TTree * BdJpsiKSTree_LL = (TTree *)MCFile->Get("candBdJpsiKS_LL");
		TTree * BdJpsiKSTree_DD = (TTree *)MCFile->Get("candBdJpsiKS_DD");
		TList * list = new TList;
		list->Add(BdJpsiKSTree_LL);
		list->Add(BdJpsiKSTree_DD);
		//TTree * BdKSmumuTree = (TTree *)MCFile->Get("candBdKSmumu");
		
		TFile * outFile = new TFile("Pdfsys.root","recreate");

		cout << "Computing sig PDF systematic" << endl;
	
		TCut highQ2cut = CutsDef::highQ2;
	
		TTree * BdJpsiKSTree = TTree::MergeTrees(list);

		anaLbMuMu_default->SetSignal("DCB");
		//anaLbMuMu_default->addBkgComponent("KSmumu",BdKSmumuTree);
		anaLbMuMu_default->Initialize();
		anaLbMuMu_default->Fit(5250.,6000.,50,true,"",highQ2cut);
		double nsigDCB = anaLbMuMu_default->GetSigVal();

		anaLbMuMu_CBGauss->SetSignal("CBAndGauss-sg20");
		//anaLbMuMu_Gauss->addBkgComponent("KSmumu",BdKSmumuTree);
		anaLbMuMu_CBGauss->Initialize();
		anaLbMuMu_CBGauss->Fit(5250.,6000.,50,true,"",highQ2cut);
		double nsigCBGauss = anaLbMuMu_CBGauss->GetSigVal();
	
		anaLbMuMu_Gauss->SetSignal("Gauss");
		//anaLbMuMu_Gauss->addBkgComponent("KSmumu",BdKSmumuTree);
		anaLbMuMu_Gauss->Initialize();
		anaLbMuMu_Gauss->Fit(5250.,6000.,50,true,"",highQ2cut);
		double nsigGauss = anaLbMuMu_Gauss->GetSigVal();
		
		anaLbMuMu_DGauss->SetSignal("DGauss");
		//anaLbMuMu_Gauss->addBkgComponent("KSmumu",BdKSmumuTree);
		anaLbMuMu_DGauss->Initialize();
		anaLbMuMu_DGauss->Fit(5250.,6000.,50,true,"",highQ2cut);
		double nsigDGauss = anaLbMuMu_DGauss->GetSigVal();
		
		anaLbMuMu_CB->SetSignal("CB");
		//anaLbMuMu_CB->addBkgComponent("KSmumu",BdKSmumuTree);
		anaLbMuMu_CB->Initialize();
		anaLbMuMu_CB->Fit(5250.,6000.,50,true,"",highQ2cut);
		double nsigCB = anaLbMuMu_CB->GetSigVal();
	

		Analysis * KS = new Analysis("KS_bkg","Lb",BdJpsiKSTree,vJpsi,"DCB_OST");
		KS->Fit(5250.,6000.,100,true,"");
		Str2VarMap pars = KS->GetSigParams();
		setConstant(&pars);
		RooRealVar * m_shift = new RooRealVar("m_shift","shift",0.,-100.,100.);
		ModifyPars(&pars,"m",m_shift,"-shift");
	
		anaLbJpsi_default->SetSignal("DCB",1.e4,"-namepar");
		anaLbJpsi_default->addBkgComponent("JpsiKS","DCB_OST",3.e3,"",pars);
		anaLbJpsi_default->Initialize();
		anaLbJpsi_default->Fit(5250.,6000.,200,true,"-noCost");
		double njpsiDCB = anaLbJpsi_default->GetSigVal();	
		//Str2VarMap Jpsipars = anaLbJpsi_default->GetSigParams();
		//setConstant(&Jpsipars);
	
		anaLbJpsi_CBGauss->SetSignal("CBAndGauss",1.e4,"-namepar");
		anaLbJpsi_CBGauss->addBkgComponent("JpsiKS","DCB_OST",3.e3,"",pars);
		anaLbJpsi_CBGauss->Initialize();
		anaLbJpsi_CBGauss->Fit(5250.,6000.,200,true,"-noCost");
		double njpsiCBGauss = anaLbJpsi_CBGauss->GetSigVal();
		
		anaLbJpsi_DGauss->SetSignal("DGauss",1.e4,"-namepar");
		anaLbJpsi_DGauss->addBkgComponent("JpsiKS","DCB_OST",3.e3,"",pars);
		anaLbJpsi_DGauss->Initialize();
		anaLbJpsi_DGauss->Fit(5250.,6000.,200,true,"-noCost");
		double njpsiDGauss = anaLbJpsi_DGauss->GetSigVal();
		
		anaLbJpsi_CB->SetSignal("CB",1.e4,"-namepar");
		anaLbJpsi_CB->addBkgComponent("JpsiKS","DCB_OST",3.e3,"",pars);
		anaLbJpsi_CB->Initialize();
		anaLbJpsi_CB->Fit(5250.,6000.,200,true,"-noCost");
		double njpsiCB = anaLbJpsi_CB->GetSigVal();
		
		anaLbJpsi_ChebKS->SetSignal("DCB",1.e4,"-namepar");
		m_shift->setVal(0);
		m_shift->setConstant();
		anaLbJpsi_ChebKS->addBkgComponent("JpsiKS","DCB_OST",3.e3,"",pars);
		//anaLbJpsi_ChebKS->addBkgComponent("JpsiKS",BdJpsiKSTree);
		anaLbJpsi_ChebKS->Initialize();
		anaLbJpsi_ChebKS->Fit(5250.,6000.,200,true,"-noCost");
		double njpsi_ChebKS = anaLbJpsi_ChebKS->GetSigVal();
		
		cout << "-----------------------------------------------------------------" << endl;	
		anaLbMuMu_default->PrintChi2();
		cout << "Nevts: " << nsigDCB << endl;
		anaLbMuMu_CBGauss->PrintChi2();
		cout << "Nevts: " << nsigCBGauss << endl;
		anaLbMuMu_Gauss->PrintChi2();
		cout << "Nevts: " << nsigGauss << endl;
		anaLbMuMu_DGauss->PrintChi2();
		cout << "Nevts: " << nsigDGauss << endl;
		anaLbMuMu_CB->PrintChi2();
		cout << "Nevts: " << nsigCB << endl;
		cout << "-----------------------------------------------------------------" << endl;		
		anaLbJpsi_default->PrintChi2();
		cout << "Nevts: " << njpsiDCB << endl;
		anaLbJpsi_CBGauss->PrintChi2();
		cout << "Nevts: " << njpsiCBGauss << endl;
		anaLbJpsi_DGauss->PrintChi2();
		cout << "Nevts: " << njpsiDGauss << endl;
		anaLbJpsi_CB->PrintChi2();
		cout << "Nevts: " << njpsiCB << endl;
		anaLbJpsi_ChebKS->PrintChi2();
		cout << "Nevts: " << njpsi_ChebKS << endl;
		

		double rarePdfsys = TMath::Abs((nsigGauss-nsigDCB)/nsigDCB)*100;
		double jpsiPdfsys = TMath::Abs((njpsiCB - njpsiDCB) / njpsiDCB)*100;
		double bkgPdfsys = TMath::Abs(njpsi_ChebKS - njpsiDCB) / njpsiDCB*100;
		double relPdfsys = TMath::Sqrt(TMath::Power(rarePdfsys,2)+TMath::Power(jpsiPdfsys,2)+TMath::Power(bkgPdfsys,2));

		cout << "N_DCB = " << nsigDCB << "  N_Gauss = " << nsigGauss << "   Diff = " << TMath::Abs(nsigGauss - nsigDCB) << endl;
		cout << "Njpsi_DCB = " << njpsiDCB << "  Njpsi_CB = " << njpsiCB  << "   Diff = " << TMath::Abs(njpsiCB - njpsiDCB) << endl;
		cout << "Njpsi_Cheb = " << njpsi_ChebKS << "   " << bkgPdfsys << "%"<< endl;
		cout << fixed << setprecision(5) << "PDFsys = " << rarePdfsys << "%" << endl;
		cout << fixed << setprecision(5) << "PDFsys = " << jpsiPdfsys << "%" << endl;
		cout << fixed << setprecision(5) << "Relative tot sys = " << relPdfsys << "%" << endl;

		//cout << fixed << setprecision(5) << "Relative PDFSigsys = " << TMath::Abs(((nsigDGauss/njpsiCB)-(nsigDCB/njpsiDCB))/(nsigDCB/njpsiDCB))*100 << "%" << endl;
		//cout << fixed << setprecision(5) << "Relative PDFBkgsys = " << TMath::Abs(((nsigGauss/njpsi_ChebKS)-(nsigGauss/njpsiDCB))/(nsigGauss/njpsiDCB))*100 << "%\n\n" << endl;
		

		MCFile->Close();
		outFile->Write();
		outFile->Close();
		delete outFile;
		delete MCfile;
	}
	
	
	if(!doEff)  return 0;
	
	
	/**     Systematic on efficiency      */
	
	TCut baseCut = CutsDef::mumuTrueID + CutsDef::avoidJpsiCut + CutsDef::massCutUnblinded;
	TCut baseJpsiCut = CutsDef::jpsiTrueID + CutsDef::jpsiCut + CutsDef::massCutUnblinded;
	if(type == "DD") { baseCut += (TCut)"pplus_TRACK_Type == 3"; baseJpsiCut += (TCut)"pplus_TRACK_Type == 3"; }
	else if(type == "LL") { baseCut += (TCut)"pplus_TRACK_Type == 5"; baseJpsiCut += (TCut)"pplus_TRACK_Type == 5"; }
	cout << "Analysisng " << type << " events" << endl;
		
	TCut geomCut = "TMath::Abs(TMath::ACos(muplus_TRUEP_Z / TMath::Sqrt(TMath::Power(muplus_TRUEP_Z,2) + TMath::Power(muplus_TRUEP_Y,2) + TMath::Power(muplus_TRUEP_X,2)))) > 0.01 && TMath::Abs(TMath::ACos(muplus_TRUEP_Z / TMath::Sqrt(TMath::Power(muplus_TRUEP_Z,2) + TMath::Power(muplus_TRUEP_Y,2) + TMath::Power(muplus_TRUEP_X,2)))) < 0.4 && TMath::Abs(TMath::ACos(muminus_TRUEP_Z / TMath::Sqrt(TMath::Power(muminus_TRUEP_Z,2) + TMath::Power(muminus_TRUEP_Y,2) + TMath::Power(muminus_TRUEP_X,2)))) > 0.01 && TMath::Abs(TMath::ACos(muminus_TRUEP_Z / TMath::Sqrt(TMath::Power(muminus_TRUEP_Z,2) + TMath::Power(muminus_TRUEP_Y,2) + TMath::Power(muminus_TRUEP_X,2)))) < 0.4";

	
	vector <TCut> cuts_common;
	cuts_common.push_back(CutsDef::TrigPassed);
	cuts_common.push_back(CutsDef::MVAcut);
	vector <string> cutsnames;
	cutsnames.push_back("base");
	cutsnames.push_back("trig");
	cutsnames.push_back("mva");
	vector <string> effnames;
	effnames.push_back("base");
	effnames.push_back("geom");
	effnames.push_back("det");
	effnames.push_back("reco");
	effnames.push_back("trig");
	effnames.push_back("mva");
	
	int nbins = CutsDef::nq2bins;
	double * q2min = &CutsDef::q2min[0];
	double * q2max = &CutsDef::q2max[0];
	
	vector <TGraphErrors> gr(effnames.size()-1,TGraphErrors());
	vector <double> tmp;
	vector <vector <double> > lfsys_plus(nbins,tmp), lfsys_minus(nbins,tmp);
	vector <vector <double> > decaysys(nbins,tmp);
	vector <vector <double> > DDsys(nbins,tmp);
	vector <vector <double> > ndefault(nbins,tmp);
	vector <vector <double> > staterr(nbins,tmp);
	vector <vector <double> > polsys_minus(nbins,tmp), polsys_plus(nbins,tmp);
	vector <vector <double> > poljpsi1_minus(nbins,tmp), poljpsi1_plus(nbins,tmp), poljpsi2_minus(nbins,tmp), poljpsi2_plus(nbins,tmp), poljpsi3_minus(nbins,tmp), poljpsi3_plus(nbins,tmp), poljpsi4_minus(nbins,tmp), poljpsi4_plus(nbins,tmp);
		
		
	for(int i = 0; i < nbins; i++)
	{
		cout << "Analysing systematics for bin " << i+1 << "/" << nbins << "  ( " << q2min[i] << " - " << q2max[i] << " GeV )" << endl;
		
		vector <TCut> cuts(cuts_common);
		TString curq2cut = Form("TMath::Power(J_psi_1S_MM/1000,2) >= %e && TMath::Power(J_psi_1S_MM/1000,2) < %e",q2min[i],q2max[i]);	
		
		/**  Sys on Reco and detector interaction efficiency  */

		double staterror;
		double geomfactor = 0.5;
		if(rel) geomfactor = 1.;
		double geomeff = 0, deteff = 0, recoeff = 0;
		
		if(print.find("geom") == string::npos || print.find("All") == string::npos)
		geomeff = geomfactor*getEff(Form("nGEOdefault_%i",i), "Lb_TRUEP_E",
					treeMCgeom, (TCut)curq2cut + geomCut, treeMCgeom, (TCut)curq2cut,
					weight+"*lifeTimeW*physRate_polp006_noDecay",
					&staterror,
					treeMCjpsi_geom, geomCut, treeMCjpsi_geom, "", weight+"*lifeTimeW*physRate_pol0_noDecay");
		ndefault[i].push_back(geomeff);
		staterr[i].push_back(geomfactor*staterror);
		
		if(print.find("det") == string::npos || print.find("All") == string::npos)
		deteff = getEff(Form("nDETdefault_%i",i), "Lb_TRUEP_E",
					treeMCGen, (TCut)curq2cut, treeMCAllGen, (TCut)curq2cut,
					weight+"*lifeTimeW*physRate_polp006_noDecay",
					&staterror,
					treeMCjpsi_Gen, "", treeMCjpsi_AllGen, "",  weight+"*lifeTimeW*physRate_pol0_noDecay");
		ndefault[i].push_back(deteff);
		staterr[i].push_back(staterror);
		
		if(print.find("reco") == string::npos || print.find("All") == string::npos)
		recoeff = getEff(Form("nRECOdefault_%i",i), "Lb_TRUEP_E",
					treeMC, (TCut)curq2cut + baseCut, treeMCGen, (TCut)curq2cut,
					weight+"*lifeTimeW*physRate_polp006",
					&staterror,
					treeMCjpsi, baseJpsiCut, treeMCjpsi_Gen, "", weight+"*lifeTimeW*physRate_pol0");
		ndefault[i].push_back(recoeff);
		staterr[i].push_back(staterror);
		
		vector< vector< double > > eff = getEff(Form("effdefault_%i",i), "Lb_MM", treeMC,  baseCut + (TCut)curq2cut, cuts,
				weight+"*lifeTimeW*physRate_polp006", treeMCjpsi, baseJpsiCut, weight+"*lifeTimeW*physRate_pol0");
		ndefault[i].insert(ndefault[i].end(), eff[0].begin(), eff[0].end());
		staterr[i].insert(staterr[i].end(), eff[1].begin(), eff[1].end());
		

		for (unsigned j = 0; j < ndefault[i].size(); j++)
		{
			gr[j].SetPoint(i,(q2min[i] + q2max[i])/2.,ndefault[i][j]);
			gr[j].SetPointError(i,(q2max[i] - q2min[i])/2.,staterr[i][j]);
		}

		if(doEffsys)
		{
			/* Lifetime, Polarization and DecayModel systematics on GEOM eff (+ DD sys empty)*/

			lfsys_plus[i].push_back( geomfactor*getEff(Form("nGEOlfplus_%i",i), "Lb_TRUEP_E", treeMCgeom, (TCut)curq2cut + geomCut, treeMCgeom, (TCut)curq2cut,
					weight+"*lifeTimeW_plussigma*physRate_polp006_noDecay", NULL, treeMCjpsi_geom, geomCut, treeMCjpsi_geom, "", weight+"*lifeTimeW_plussigma*physRate_pol0_noDecay") );
			lfsys_minus[i].push_back( geomfactor*getEff(Form("nGEOlfminus_%i",i), "Lb_TRUEP_E", treeMCgeom, (TCut)curq2cut + geomCut, treeMCgeom, (TCut)curq2cut,
					weight+"*lifeTimeW_minussigma*physRate_polp006_noDecay", NULL, treeMCjpsi_geom, geomCut, treeMCjpsi_geom, "", weight+"*lifeTimeW_minussigma*physRate_pol0_noDecay") );
			
			polsys_plus[i].push_back( geomfactor*getEff(Form("nGEOpolplus_%i",i), "Lb_TRUEP_E", treeMCgeom, (TCut)curq2cut + geomCut, treeMCgeom, (TCut)curq2cut,
					weight+"*lifeTimeW*physRate_polp015_noDecay", NULL, treeMCjpsi_geom, geomCut, treeMCjpsi_geom, "", weight+"*lifeTimeW*physRate_pol0_noDecay") );
			polsys_minus[i].push_back( geomfactor*getEff(Form("nGEOpolminus_%i",i), "Lb_TRUEP_E", treeMCgeom, (TCut)curq2cut + geomCut, treeMCgeom, (TCut)curq2cut,
					weight+"*lifeTimeW*physRate_polm003_noDecay", NULL, treeMCjpsi_geom, geomCut, treeMCjpsi_geom, "", weight+"*lifeTimeW*physRate_pol0_noDecay") );
			
			poljpsi1_plus[i].push_back( geomfactor*getEff(Form("nGEOpoljpsiplus_%i",i), "Lb_TRUEP_E", treeMCgeom, (TCut)curq2cut + geomCut, treeMCgeom, (TCut)curq2cut,
					weight+"*lifeTimeW*physRate_polp006_noDecay", NULL, treeMCjpsi_geom, geomCut, treeMCjpsi_geom, "", weight+"*lifeTimeW*model_jpsi_1_noDecay") );
			poljpsi1_minus[i].push_back( geomfactor*getEff(Form("nGEOpoljpsiminus_%i",i), "Lb_TRUEP_E", treeMCgeom, (TCut)curq2cut + geomCut, treeMCgeom, (TCut)curq2cut,
					weight+"*lifeTimeW*physRate_polp006_noDecay", NULL, treeMCjpsi_geom, geomCut, treeMCjpsi_geom, "", weight+"*lifeTimeW*model_jpsi_2_noDecay") );
			
			poljpsi2_plus[i].push_back( geomfactor*getEff(Form("nGEOpoljpsiplus_%i",i), "Lb_TRUEP_E", treeMCgeom, (TCut)curq2cut + geomCut, treeMCgeom, (TCut)curq2cut,
					weight+"*lifeTimeW*physRate_polp006_noDecay", NULL, treeMCjpsi_geom, geomCut, treeMCjpsi_geom, "", weight+"*lifeTimeW*model_jpsi_3_noDecay") );
			poljpsi2_minus[i].push_back( geomfactor*getEff(Form("nGEOpoljpsiminus_%i",i), "Lb_TRUEP_E", treeMCgeom, (TCut)curq2cut + geomCut, treeMCgeom, (TCut)curq2cut,
					weight+"*lifeTimeW*physRate_polp006_noDecay", NULL, treeMCjpsi_geom, geomCut, treeMCjpsi_geom, "", weight+"*lifeTimeW*model_jpsi_4_noDecay") );
			
			poljpsi3_plus[i].push_back( geomfactor*getEff(Form("nGEOpoljpsiplus_%i",i), "Lb_TRUEP_E", treeMCgeom, (TCut)curq2cut + geomCut, treeMCgeom, (TCut)curq2cut,
					weight+"*lifeTimeW*physRate_polp006_noDecay", NULL, treeMCjpsi_geom, geomCut, treeMCjpsi_geom, "", weight+"*lifeTimeW*model_jpsi_5_noDecay") );
			poljpsi3_minus[i].push_back( geomfactor*getEff(Form("nGEOpoljpsiminus_%i",i), "Lb_TRUEP_E", treeMCgeom, (TCut)curq2cut + geomCut, treeMCgeom, (TCut)curq2cut,
					weight+"*lifeTimeW*physRate_polp006_noDecay", NULL, treeMCjpsi_geom, geomCut, treeMCjpsi_geom, "", weight+"*lifeTimeW*model_jpsi_6_noDecay") );
			
			poljpsi4_plus[i].push_back( geomfactor*getEff(Form("nGEOpoljpsiplus_%i",i), "Lb_TRUEP_E", treeMCgeom, (TCut)curq2cut + geomCut, treeMCgeom, (TCut)curq2cut,
					weight+"*lifeTimeW*physRate_polp006_noDecay", NULL, treeMCjpsi_geom, geomCut, treeMCjpsi_geom, "", weight+"*lifeTimeW*model_jpsi_7_noDecay") );
			poljpsi4_minus[i].push_back( geomfactor*getEff(Form("nGEOpoljpsiminus_%i",i), "Lb_TRUEP_E", treeMCgeom, (TCut)curq2cut + geomCut, treeMCgeom, (TCut)curq2cut,
					weight+"*lifeTimeW*physRate_polp006_noDecay", NULL, treeMCjpsi_geom, geomCut, treeMCjpsi_geom, "", weight+"*lifeTimeW*model_jpsi_8_noDecay") );
			
			decaysys[i].push_back( geomfactor*getEff(Form("nGEOdecaysys_%i",i), "Lb_TRUEP_E", treeMCgeom, (TCut)curq2cut + geomCut, treeMCgeom, (TCut)curq2cut,
					weight+"*lifeTimeW*physRate_pol0_QCDff_noDecay", NULL, treeMCjpsi_geom, geomCut, treeMCjpsi_geom, "", weight+"*lifeTimeW*physRate_pol0_noDecay") );
			
			if(type!="LL") DDsys[i].push_back( geomeff );   //For DDvtx sys
		
		
			/* Lifetime, Polarization and DecayModel systematics on DET eff (+ DD sys empty)*/
cout << "DET EFF" << endl;
			lfsys_plus[i].push_back( getEff(Form("nDETlfplus_%i",i), "Lb_TRUEP_E", treeMCGen, (TCut)curq2cut, treeMCAllGen, (TCut)curq2cut,
					weight+"*lifeTimeW_plussigma*physRate_polp006_noDecay", NULL, treeMCjpsi_Gen, "", treeMCjpsi_AllGen, "", weight+"*lifeTimeW_plussigma*physRate_pol0_noDecay") );
			lfsys_minus[i].push_back( getEff(Form("nDETlfminus_%i",i), "Lb_TRUEP_E", treeMCGen, (TCut)curq2cut, treeMCAllGen, (TCut)curq2cut,
					weight+"*lifeTimeW_minussigma*physRate_polp006_noDecay", NULL, treeMCjpsi_Gen, "", treeMCjpsi_AllGen, "", weight+"*lifeTimeW_minussigma*physRate_pol0_noDecay") );
			
			polsys_plus[i].push_back( getEff(Form("nDETpolplus_%i",i), "Lb_TRUEP_E", treeMCGen, (TCut)curq2cut, treeMCAllGen, (TCut)curq2cut,
					weight+"*lifeTimeW*physRate_polp015_noDecay", NULL, treeMCjpsi_Gen, "", treeMCjpsi_AllGen, "", weight+"*lifeTimeW*physRate_pol0_noDecay") );
			polsys_minus[i].push_back( getEff(Form("nDETpolminus_%i",i), "Lb_TRUEP_E", treeMCGen, (TCut)curq2cut, treeMCAllGen, (TCut)curq2cut,
					weight+"*lifeTimeW*physRate_polm003_noDecay", NULL, treeMCjpsi_Gen, "", treeMCjpsi_AllGen, "", weight+"*lifeTimeW*physRate_pol0_noDecay") );
			
			poljpsi1_plus[i].push_back( getEff(Form("nDETpoljpsiplus_%i",i), "Lb_TRUEP_E", treeMCGen, (TCut)curq2cut, treeMCAllGen, (TCut)curq2cut,
					weight+"*lifeTimeW*physRate_polp006_noDecay", NULL, treeMCjpsi_Gen, "", treeMCjpsi_AllGen, "", weight+"*lifeTimeW*model_jpsi_1_noDecay") );
			poljpsi1_minus[i].push_back( getEff(Form("nDETpoljpsiminus_%i",i), "Lb_TRUEP_E", treeMCGen, (TCut)curq2cut, treeMCAllGen, (TCut)curq2cut,
					weight+"*lifeTimeW*physRate_polp006_noDecay", NULL, treeMCjpsi_Gen, "", treeMCjpsi_AllGen, "", weight+"*lifeTimeW*model_jpsi_2_noDecay") );
			
			poljpsi2_plus[i].push_back( getEff(Form("nDETpoljpsiplus_%i",i), "Lb_TRUEP_E", treeMCGen, (TCut)curq2cut, treeMCAllGen, (TCut)curq2cut,
					weight+"*lifeTimeW*physRate_polp006_noDecay", NULL, treeMCjpsi_Gen, "", treeMCjpsi_AllGen, "", weight+"*lifeTimeW*model_jpsi_3_noDecay") );
			poljpsi2_minus[i].push_back( getEff(Form("nDETpoljpsiminus_%i",i), "Lb_TRUEP_E", treeMCGen, (TCut)curq2cut, treeMCAllGen, (TCut)curq2cut,
					weight+"*lifeTimeW*physRate_polp006_noDecay", NULL, treeMCjpsi_Gen, "", treeMCjpsi_AllGen, "", weight+"*lifeTimeW*model_jpsi_4_noDecay") );
			
			poljpsi3_plus[i].push_back( getEff(Form("nDETpoljpsiplus_%i",i), "Lb_TRUEP_E", treeMCGen, (TCut)curq2cut, treeMCAllGen, (TCut)curq2cut,
					weight+"*lifeTimeW*physRate_polp006_noDecay", NULL, treeMCjpsi_Gen, "", treeMCjpsi_AllGen, "", weight+"*lifeTimeW*model_jpsi_5_noDecay") );
			poljpsi3_minus[i].push_back( getEff(Form("nDETpoljpsiminus_%i",i), "Lb_TRUEP_E", treeMCGen, (TCut)curq2cut, treeMCAllGen, (TCut)curq2cut,
					weight+"*lifeTimeW*physRate_polp006_noDecay", NULL, treeMCjpsi_Gen, "", treeMCjpsi_AllGen, "", weight+"*lifeTimeW*model_jpsi_6_noDecay") );
			
			poljpsi4_plus[i].push_back( getEff(Form("nDETpoljpsiplus_%i",i), "Lb_TRUEP_E", treeMCGen, (TCut)curq2cut, treeMCAllGen, (TCut)curq2cut,
					weight+"*lifeTimeW*physRate_polp006_noDecay", NULL, treeMCjpsi_Gen, "", treeMCjpsi_AllGen, "", weight+"*lifeTimeW*model_jpsi_7_noDecay") );
			poljpsi4_minus[i].push_back( getEff(Form("nDETpoljpsiminus_%i",i), "Lb_TRUEP_E", treeMCGen, (TCut)curq2cut, treeMCAllGen, (TCut)curq2cut,
					weight+"*lifeTimeW*physRate_polp006_noDecay", NULL, treeMCjpsi_Gen, "", treeMCjpsi_AllGen, "", weight+"*lifeTimeW*model_jpsi_8_noDecay") );
		
			decaysys[i].push_back( getEff(Form("nDETdecaysys_%i",i), "Lb_TRUEP_E", treeMCGen, (TCut)curq2cut, treeMCAllGen, (TCut)curq2cut,
				weight+"*lifeTimeW*physRate_pol0_QCDff_noDecay", NULL, treeMCjpsi_Gen, "", treeMCjpsi_AllGen, "", weight+"*lifeTimeW*physRate_pol0_noDecay") );

			if(type!="LL") DDsys[i].push_back( deteff );
			
cout << "RECO EFF" << endl;
			/* Lifetime, Polarization, DD vtx eff and DecayModel systematics on RECO eff */
			
			lfsys_plus[i].push_back( getEff(Form("nRECOlfplus_%i",i), "Lb_TRUEP_E", treeMC, (TCut)curq2cut + baseCut, treeMCGen, (TCut)curq2cut,
				weight+"*lifeTimeW_plussigma*physRate_polp006", NULL, treeMCjpsi, baseJpsiCut, treeMCjpsi_Gen, "", weight+"*lifeTimeW_plussigma*physRate_pol0") );
			lfsys_minus[i].push_back( getEff(Form("nRECOlfminus_%i",i), "Lb_TRUEP_E", treeMC, (TCut)curq2cut + baseCut, treeMCGen, (TCut)curq2cut,
				weight+"*lifeTimeW_minussigma*physRate_polp006", NULL, treeMCjpsi, baseJpsiCut, treeMCjpsi_Gen, "", weight+"*lifeTimeW_minussigma*physRate_pol0") );
			
			polsys_plus[i].push_back( getEff(Form("nRECOpolplus_%i",i), "Lb_TRUEP_E", treeMC, (TCut)curq2cut + baseCut, treeMCGen, (TCut)curq2cut,
				weight+"*lifeTimeW*physRate_polp015", NULL, treeMCjpsi, baseJpsiCut, treeMCjpsi_Gen, "", weight+"*lifeTimeW*physRate_pol0") );
			polsys_minus[i].push_back( getEff(Form("nRECOpolminus_%i",i), "Lb_TRUEP_E", treeMC, (TCut)curq2cut + baseCut, treeMCGen, (TCut)curq2cut,
				weight+"*lifeTimeW*physRate_polm003", NULL, treeMCjpsi, baseJpsiCut, treeMCjpsi_Gen, "", weight+"*lifeTimeW*physRate_pol0") );
			
			poljpsi1_plus[i].push_back( getEff(Form("nRECOpoljpsiplus_%i",i), "Lb_TRUEP_E", treeMC, (TCut)curq2cut + baseCut, treeMCGen, (TCut)curq2cut,
				weight+"*lifeTimeW*physRate_polp006", NULL, treeMCjpsi, baseJpsiCut, treeMCjpsi_Gen, "", weight+"*lifeTimeW*model_jpsi_1") );
			poljpsi1_minus[i].push_back( getEff(Form("nRECOpoljpsiminus_%i",i), "Lb_TRUEP_E", treeMC, (TCut)curq2cut + baseCut, treeMCGen, (TCut)curq2cut,
				weight+"*lifeTimeW*physRate_polp006", NULL, treeMCjpsi, baseJpsiCut, treeMCjpsi_Gen, "", weight+"*lifeTimeW*model_jpsi_2") );
			
			poljpsi2_plus[i].push_back( getEff(Form("nRECOpoljpsiplus_%i",i), "Lb_TRUEP_E", treeMC, (TCut)curq2cut + baseCut, treeMCGen, (TCut)curq2cut,
				weight+"*lifeTimeW*physRate_polp006", NULL, treeMCjpsi, baseJpsiCut, treeMCjpsi_Gen, "", weight+"*lifeTimeW*model_jpsi_3") );
			poljpsi2_minus[i].push_back( getEff(Form("nRECOpoljpsiminus_%i",i), "Lb_TRUEP_E", treeMC, (TCut)curq2cut + baseCut, treeMCGen, (TCut)curq2cut,
				weight+"*lifeTimeW*physRate_polp006", NULL, treeMCjpsi, baseJpsiCut, treeMCjpsi_Gen, "", weight+"*lifeTimeW*model_jpsi_4") );
			
			poljpsi3_plus[i].push_back( getEff(Form("nRECOpoljpsiplus_%i",i), "Lb_TRUEP_E", treeMC, (TCut)curq2cut + baseCut, treeMCGen, (TCut)curq2cut,
				weight+"*lifeTimeW*physRate_polp006", NULL, treeMCjpsi, baseJpsiCut, treeMCjpsi_Gen, "", weight+"*lifeTimeW*model_jpsi_5") );
			poljpsi3_minus[i].push_back( getEff(Form("nRECOpoljpsiminus_%i",i), "Lb_TRUEP_E", treeMC, (TCut)curq2cut + baseCut, treeMCGen, (TCut)curq2cut,
				weight+"*lifeTimeW*physRate_polp006", NULL, treeMCjpsi, baseJpsiCut, treeMCjpsi_Gen, "", weight+"*lifeTimeW*model_jpsi_6") );
			
			poljpsi4_plus[i].push_back( getEff(Form("nRECOpoljpsiplus_%i",i), "Lb_TRUEP_E", treeMC, (TCut)curq2cut + baseCut, treeMCGen, (TCut)curq2cut,
				weight+"*lifeTimeW*physRate_polp006", NULL, treeMCjpsi, baseJpsiCut, treeMCjpsi_Gen, "", weight+"*lifeTimeW*model_jpsi_7") );
			poljpsi4_minus[i].push_back( getEff(Form("nRECOpoljpsiminus_%i",i), "Lb_TRUEP_E", treeMC, (TCut)curq2cut + baseCut, treeMCGen, (TCut)curq2cut,
				weight+"*lifeTimeW*physRate_polp006", NULL, treeMCjpsi, baseJpsiCut, treeMCjpsi_Gen, "", weight+"*lifeTimeW*model_jpsi_8") );
			
			
			decaysys[i].push_back( getEff(Form("nRECOdecaysys_%i",i), "Lb_TRUEP_E", treeMC, (TCut)curq2cut + baseCut, treeMCGen, (TCut)curq2cut,
				weight+"*lifeTimeW*physRate_pol0_QCDff", NULL, treeMCjpsi, baseJpsiCut, treeMCjpsi_Gen, "", weight+"*lifeTimeW*physRate_pol0") );
			
			if(type!="LL") DDsys[i].push_back( getEff(Form("nRECODDsys_%i",i), "Lb_TRUEP_E", treeMC,
				(TCut)("DDvtx_weight*("+(TString)curq2cut+" && "+(TString)baseCut+")"), treeMCGen, (TCut)curq2cut,
				weight+"*lifeTimeW*physRate_polp006", NULL, treeMCjpsi, (TCut)("DDvtx_weight*("+(TString)baseJpsiCut+")"), treeMCjpsi_Gen,"",weight+"*lifeTimeW*physRate_pol0") );
	
			/**     Lifetime systematic on Trig and MVA eff  */
			
			tmp = getEff(Form("lfminus_%i",i), "Lb_MM", treeMC,  baseCut + (TCut)curq2cut, cuts, weight+"*lifeTimeW_minussigma*physRate_polp006", treeMCjpsi, baseJpsiCut, weight+"*lifeTimeW_minussigma*physRate_pol0")[0];
			lfsys_minus[i].insert(lfsys_minus[i].end(), tmp.begin(), tmp.end());
			tmp = getEff(Form("lfplus_%i",i), "Lb_MM", treeMC,  baseCut + (TCut)curq2cut, cuts, weight+"*lifeTimeW_plussigma*physRate_polp006", treeMCjpsi, baseJpsiCut, weight+"*lifeTimeW_plussigma*physRate_pol0")[0];
			lfsys_plus[i].insert(lfsys_plus[i].end(), tmp.begin(), tmp.end());
			
			/**     Polarization systematic on Trig and MVA eff       */
			
			tmp = getEff(Form("polplus_%i",i), "Lb_MM", treeMC,  baseCut + (TCut)curq2cut, cuts, weight+"*lifeTimeW*physRate_polm003", treeMCjpsi, baseJpsiCut, weight+"*lifeTimeW*physRate_pol0")[0];
			polsys_minus[i].insert(polsys_minus[i].end(), tmp.begin(), tmp.end());
			tmp = getEff(Form("polminus_%i",i), "Lb_MM", treeMC,  baseCut + (TCut)curq2cut, cuts, weight+"*lifeTimeW*physRate_polp015", treeMCjpsi, baseJpsiCut, weight+"*lifeTimeW*physRate_pol0")[0];
			polsys_plus[i].insert(polsys_plus[i].end(), tmp.begin(), tmp.end());
			

			tmp = getEff(Form("polplus_%i",i), "Lb_MM", treeMC,  baseCut + (TCut)curq2cut, cuts, weight+"*lifeTimeW*physRate_polp006", treeMCjpsi, baseJpsiCut, weight+"*lifeTimeW*model_jpsi_1")[0];
			poljpsi1_minus[i].insert(poljpsi1_minus[i].end(), tmp.begin(), tmp.end());
			tmp = getEff(Form("polminus_%i",i), "Lb_MM", treeMC,  baseCut + (TCut)curq2cut, cuts, weight+"*lifeTimeW*physRate_polp006", treeMCjpsi, baseJpsiCut, weight+"*lifeTimeW*model_jpsi_2")[0];
			poljpsi1_plus[i].insert(poljpsi1_plus[i].end(), tmp.begin(), tmp.end());
			tmp = getEff(Form("polplus_%i",i), "Lb_MM", treeMC,  baseCut + (TCut)curq2cut, cuts, weight+"*lifeTimeW*physRate_polp006", treeMCjpsi, baseJpsiCut, weight+"*lifeTimeW*model_jpsi_3")[0];
			poljpsi2_minus[i].insert(poljpsi2_minus[i].end(), tmp.begin(), tmp.end());
			tmp = getEff(Form("polminus_%i",i), "Lb_MM", treeMC,  baseCut + (TCut)curq2cut, cuts, weight+"*lifeTimeW*physRate_polp006", treeMCjpsi, baseJpsiCut, weight+"*lifeTimeW*model_jpsi_4")[0];
			poljpsi2_plus[i].insert(poljpsi2_plus[i].end(), tmp.begin(), tmp.end());
			tmp = getEff(Form("polplus_%i",i), "Lb_MM", treeMC,  baseCut + (TCut)curq2cut, cuts, weight+"*lifeTimeW*physRate_polp006", treeMCjpsi, baseJpsiCut, weight+"*lifeTimeW*model_jpsi_5")[0];
			poljpsi3_minus[i].insert(poljpsi3_minus[i].end(), tmp.begin(), tmp.end());
			tmp = getEff(Form("polminus_%i",i), "Lb_MM", treeMC,  baseCut + (TCut)curq2cut, cuts, weight+"*lifeTimeW*physRate_polp006", treeMCjpsi, baseJpsiCut, weight+"*lifeTimeW*model_jpsi_6")[0];
			poljpsi3_plus[i].insert(poljpsi3_plus[i].end(), tmp.begin(), tmp.end());
			tmp = getEff(Form("polplus_%i",i), "Lb_MM", treeMC,  baseCut + (TCut)curq2cut, cuts, weight+"*lifeTimeW*physRate_polp006", treeMCjpsi, baseJpsiCut, weight+"*lifeTimeW*model_jpsi_7")[0];
			poljpsi4_minus[i].insert(poljpsi4_minus[i].end(), tmp.begin(), tmp.end());
			tmp = getEff(Form("polminus_%i",i), "Lb_MM", treeMC,  baseCut + (TCut)curq2cut, cuts, weight+"*lifeTimeW*physRate_polp006", treeMCjpsi, baseJpsiCut, weight+"*lifeTimeW*model_jpsi_8")[0];
			poljpsi4_plus[i].insert(poljpsi4_plus[i].end(), tmp.begin(), tmp.end());
		
			/**     Decay model systematic on Trig and MVA eff        */

			tmp = getEff(Form("decaysys_%i",i), "Lb_MM", treeMC,  baseCut + (TCut)curq2cut, cuts, weight+"*lifeTimeW*physRate_pol0_QCDff", treeMCjpsi, baseJpsiCut, weight+"*lifeTimeW*physRate_pol0")[0];
			decaysys[i].insert(decaysys[i].end(), tmp.begin(), tmp.end());
				
		
			/**     DD vtx eff systematic on Trig and MVA eff        */
		
			if(type!="LL") {
				tmp = getEff(Form("DDsys_%i",i), "Lb_MM", treeMC,  baseCut + (TCut)curq2cut, cuts, weight+"*lifeTimeW*physRate_polp006_noDecay*DDvtx_weight", treeMCjpsi, baseJpsiCut, weight+"*lifeTimeW*physRate_pol0_noDecay*DDvtx_weight")[0];
			DDsys[i].insert(DDsys[i].end(), tmp.begin(), tmp.end()); }
		}
	}
	
	
	
	
	
	
	/** Printing out efficiencies **/
	
	TFile * histFile = NULL;
	if(rel) histFile = new TFile(("Lbreleff_"+type+".root").c_str(),"recreate");
	else histFile = new TFile(("Lbeff_"+type+".root").c_str(),"recreate");
	vector <string> sysnames;
	sysnames.push_back("Lifetime");
	sysnames.push_back("Decay Model");
	if(type!="LL") sysnames.push_back("DD vtx eff");
	vector <vector < vector <double> > > plus;
	vector <vector < vector <double> > > minus;
	plus.push_back(lfsys_plus);
	minus.push_back(lfsys_minus);
	minus.push_back(decaysys);
	minus.push_back(DDsys);
	vector <TGraphErrors> sysgr(effnames.size()-1,TGraphErrors());
	
	if(doEffsys)
	{
		if(print == "All")
		{
			for(unsigned s = 0; s < sysnames.size(); s++)
			{
				cout << "\n\n" << sysnames[s] << " systematic\n" << endl;
				cout << "cut\t\t\t eff\t\t minus sig\t plus sig\t\t sys "<< endl;
		
				for(int i = 0; i < nbins; i++)
				{
					cout << fixed << setprecision(2) << q2min[i] << "-" << q2max[i] << endl;
					cout << fixed << setprecision(5);
					vector<double > v, rr;
					if(s < plus.size()) rr = printSys(effnames, ndefault[i], staterr[i], minus[s][i], plus[s][i], print,false,"-rtotal");
					else rr = printSys(effnames, ndefault[i], staterr[i], minus[s][i], v, print,false,"-rtotal");
					if(q2max[i]==20.0) cout << " ------------------------ " << endl;
				}
			}
			
			cout << "\n\nPolarization\n" << endl;
			cout << "cut\t\t\t eff\t\t minus sig\t plus sig\t\t sys "<< endl;
			for(int i = 0; i < nbins; i++)
			{
				cout << fixed << setprecision(2) << q2min[i] << "-" << q2max[i] << endl;
                                cout << fixed << setprecision(5);

				vector<double> tot;
				if(rel)
				{
					printSys(effnames, ndefault[i], staterr[i], polsys_minus[i], polsys_plus[i], "AllNone", false,  "-rpartial");
					printSys(effnames, ndefault[i], staterr[i], poljpsi1_minus[i], poljpsi1_plus[i], "AllNone", false);
					printSys(effnames, ndefault[i], staterr[i], poljpsi2_minus[i], poljpsi2_plus[i], "AllNone", false);
					printSys(effnames, ndefault[i], staterr[i], poljpsi3_minus[i], poljpsi3_plus[i], "AllNone", false);
					tot = printSys(effnames, ndefault[i], staterr[i], poljpsi4_minus[i], poljpsi4_plus[i], "All", false, "-ppartial-ototal");
				}
				else tot = printSys(effnames, ndefault[i], staterr[i], polsys_minus[i], polsys_plus[i], "All", false, "-ototal");

				for(unsigned t = 0; t < tot.size(); t++) sysgr[t].SetPoint(i,(q2min[i] + q2max[i])/2,tot[t]);
			}
		}
		else
		{
			unsigned nloop = 1;
			if(print=="all") nloop = effnames.size();
			for( unsigned e = 1; e < nloop; e++ )
			{
				if(nloop > 1) print = effnames[e];
				cout << "\n\nEfficiency: " << print << endl;
				
				cout << "\n\nq2 bin\t\t& Statistics";
				for(unsigned s = 0; s < sysnames.size(); s++) cout << "\t& " << sysnames[s];
				cout << "\t& Polarization\t& Total Sys \\\\ \n \\hline" << endl;
				
				for(int i = 0; i < nbins; i++)
				{
					vector<double > myeff, v, tot;
					cout << fixed << setprecision(2) << q2min[i] << "-" << q2max[i] << fixed << setprecision(5);
					double totalerror = 0;
				
					for(unsigned s = 0; s < sysnames.size(); s++)
					{
						if(s == 0)
						{
							myeff = printSys(effnames, ndefault[i], staterr[i], print, true,"-rtotal");
							totalerror *= TMath::Power(myeff[1],2);
						}
						if(s < plus.size()) myeff = printSys(effnames, ndefault[i], staterr[i], minus[s][i], plus[s][i], print);
						else myeff = printSys(effnames, ndefault[i], staterr[i], minus[s][i], v, print);
						totalerror += TMath::Power(myeff[2],2);
					}
			
					if(rel)
					{
						printSys(effnames, ndefault[i], staterr[i], polsys_minus[i], polsys_plus[i], print+"None", false,  "-rpartial");
						printSys(effnames, ndefault[i], staterr[i], poljpsi1_minus[i], poljpsi1_plus[i], print+"None", false);
						printSys(effnames, ndefault[i], staterr[i], poljpsi2_minus[i], poljpsi2_plus[i], print+"None", false);
						printSys(effnames, ndefault[i], staterr[i], poljpsi3_minus[i], poljpsi3_plus[i], print+"None", false);
						tot = printSys(effnames, ndefault[i], staterr[i], poljpsi4_minus[i], poljpsi4_plus[i], print, false, "-ppartial-ptotal-ototal");
					}
					else tot = printSys(effnames, ndefault[i], staterr[i], polsys_minus[i], polsys_plus[i], print, false,"-ptotal-ototal");
					cout << "  \\\\" << endl;
					if(q2max[i]==20.0) cout << " \\hline" << endl;
					for(unsigned t = 0; t < tot.size(); t++) sysgr[t].SetPoint(i,(q2min[i] + q2max[i])/2,tot[t]);
				}
			}
		}
	}
	else
	{
		if(print == "All")
		{
			cout << "\n\nq2 bin";
			for(unsigned e = 1; e < effnames.size(); e++) cout << "\t\t\t& " << effnames[e];
			cout << "\t\t & Total  \\\\ \\hline" << endl;
		}
		
		for(int i = 0; i < nbins; i++)
		{
			cout << fixed << setprecision(2) << q2min[i] << "-" << q2max[i] << fixed << setprecision(5);
			vector < double > rr = printSys(effnames, ndefault[i], staterr[i], "eff"+print);
			if(q2max[i]==20.0) cout << " \\hline" << endl;
		}
	}
	

	if(sysgr.size()>0)
	{
		for(unsigned e = 0; e < sysgr.size(); e++)
		{
			string grname = "sys_";
			if(rel) grname = "sysrel_";
			grname += (string)(((TString)effnames[e+1]).ReplaceAll(" ",""));
			sysgr[e].Write(grname.c_str());
		}
	}
	
	if(gr.size()>0)
	{
		TCanvas * cc = new TCanvas();
		for(unsigned e = 0; e < gr.size(); e++)
		{
			string grname = "eff_";
			if(rel) grname = "effrel_";
			grname += (string)(((TString)effnames[e+1]).ReplaceAll(" ",""));
			gr[e].Write(grname.c_str());
			gr[e].SetMarkerStyle(22);
			gr[e].SetMarkerSize(0.9);
			gr[e].SetMarkerColor(1);
			gr[e].Draw("AP");
			gr[e].GetXaxis()->SetTitle("q2");
			gr[e].GetYaxis()->SetTitle((effnames[e+1]).c_str());
			cc->Print((grname+"_"+type+".pdf").c_str());
		}

		delete cc;
	}
	
	histFile->Write();
	histFile->Close();
	
	delete MCfile;
	delete histFile;
	return 0;
}

