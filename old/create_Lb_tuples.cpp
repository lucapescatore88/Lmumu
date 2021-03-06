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
#include "RooPolynomial.h"

#include "TLorentzVector.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TEntryList.h"
#include "TCanvas.h"
#include "TBranch.h"
#include "TCut.h"

#include "ReadTree_comp.hpp"
#include "analyser.hpp"
#include "Lb_cuts.hpp"
#include "DPHelpers.hpp"

using namespace std;
using namespace RooFit;

void addVariables(TreeReader * reader, TTree * newTree, bool reset)
{
	static bool MC = false;
	static float const0, const1, const2;
	static float cosTheta, cosThetaL, cosThetaB, phiL, phiB, dphi;
	static float cosTheta_TRUE, cosThetaL_TRUE, cosThetaB_TRUE, phiL_TRUE, phiB_TRUE, dphi_TRUE;

	if(reset)
	{
		newTree->Branch("Lb_MassCons",&const0,"Lb_MassCons/F");
		newTree->Branch("Lb_MassConsLambda",&const1,"Lb_MassConsLambda/F");
		newTree->Branch("Lb_MassConsJpsiLambda",&const2,"Lb_MassConsJpsiLambda/F");

		newTree->Branch("cosTheta",&cosTheta,"cosTheta/F");
		newTree->Branch("cosThetaL",&cosThetaL,"cosThetaL/F");
		newTree->Branch("cosThetaB",&cosThetaB,"cosThetaB/F");
		newTree->Branch("phiL",&phiL,"phiL/F");
		newTree->Branch("phiB",&phiB,"phiB/F");
		newTree->Branch("dphi",&dphi,"dphi/F");

		if(reader->HasVar("Lb_TRUEP_X"))
		{
			MC = true;
			newTree->Branch("cosTheta_TRUE",&cosTheta_TRUE,"cosTheta_TRUE/F");
			newTree->Branch("cosThetaL_TRUE",&cosThetaL_TRUE,"cosThetaL_TRUE/F");
			newTree->Branch("cosThetaB_TRUE",&cosThetaB_TRUE,"cosThetaB_TRUE/F");
			newTree->Branch("phiL_TRUE",&phiL_TRUE,"phiL_TRUE/F");
			newTree->Branch("phiB_TRUE",&phiB_TRUE,"phiB_TRUE/F");
			newTree->Branch("dphi_TRUE",&dphi_TRUE,"dphi_TRUE/F");
		}
	}
	else
	{
		int year = reader->GetValue("dataType");
		float perg = 4.e6;
		if(year == 2011) perg = 3.5e6;
		float pmom = TMath::Sqrt(perg*perg-938.3*938.3);
		TLorentzVector initialProton(0.,0.,pmom,perg);

		const0 = reader->GetValue("Lb_MassCons_M",0);
		const1 = reader->GetValue("Lb_MassConsLambda_M",0);
		const2 = reader->GetValue("Lb_MassConsJpsiAndLambda_M",0);

		TLorentzVector Lambda, pion, proton, Lb, Jpsi, mup, mum;

		Lb.SetPxPyPzE(reader->GetValue("Lb_PX"),
				reader->GetValue("Lb_PY"),
				reader->GetValue("Lb_PZ"),
				reader->GetValue("Lb_PE"));
		Jpsi.SetPxPyPzE(reader->GetValue("J_psi_1S_PX"),
				reader->GetValue("J_psi_1S_PY"),
				reader->GetValue("J_psi_1S_PZ"),
				reader->GetValue("J_psi_1S_PE"));
		mup.SetPxPyPzE(reader->GetValue("muplus_PX"),
				reader->GetValue("muplus_PY"),
				reader->GetValue("muplus_PZ"),
				reader->GetValue("muplus_PE"));
		mum.SetPxPyPzE(reader->GetValue("muminus_PX"),
				reader->GetValue("muminus_PY"),
				reader->GetValue("muminus_PZ"),
				reader->GetValue("muminus_PE"));
		Lambda.SetPxPyPzE(reader->GetValue("Lambda0_PX"),
				reader->GetValue("Lambda0_PY"),
				reader->GetValue("Lambda0_PZ"),
				reader->GetValue("Lambda0_PE"));
		proton.SetPxPyPzE(reader->GetValue("pplus_PX"),
				reader->GetValue("pplus_PY"),
				reader->GetValue("pplus_PZ"),
				reader->GetValue("pplus_PE"));
		pion.SetPxPyPzE(reader->GetValue("piminus_PX"),
				reader->GetValue("piminus_PY"),
				reader->GetValue("piminus_PZ"),
				reader->GetValue("piminus_PE"));

		int pcharge = 1;
		if(reader->GetValue("pplus_ID") < 0) pcharge = -1;

		double tcosTheta, tcosThetaL, tcosThetaB, tphiL, tphiB, tdphi;
		DPHelpers::LbPsiRAngles(initialProton,Lb,Jpsi,Lambda,mup,mum,proton,pion,pcharge,
				tcosTheta,tcosThetaL,tcosThetaB,tphiL,tphiB,tdphi);
		cosTheta = (float)tcosTheta, cosThetaL = (float)tcosThetaL, cosThetaB = (float)tcosThetaB, phiL = (float)tphiL, phiB = (float)tphiB, dphi = (float)tdphi;

		//TRUE angular variables only for MC
		if(MC)
		{
			TLorentzVector Lambda_TRUE, pion_TRUE, proton_TRUE, Lb_TRUE, Jpsi_TRUE, mup_TRUE, mum_TRUE;

			Lb_TRUE.SetPxPyPzE(reader->GetValue("Lb_TRUEP_X"),
					reader->GetValue("Lb_TRUEP_Y"),
					reader->GetValue("Lb_TRUEP_Z"),
					reader->GetValue("Lb_TRUEP_E"));
			mup_TRUE.SetPxPyPzE(reader->GetValue("muplus_TRUEP_X"),
					reader->GetValue("muplus_TRUEP_Y"),
					reader->GetValue("muplus_TRUEP_Z"),
					reader->GetValue("muplus_TRUEP_E"));
			mum_TRUE.SetPxPyPzE(reader->GetValue("muminus_TRUEP_X"),
					reader->GetValue("muminus_TRUEP_Y"),
					reader->GetValue("muminus_TRUEP_Z"),
					reader->GetValue("muminus_TRUEP_E"));
			Jpsi_TRUE = mup_TRUE + mum_TRUE;
			Lambda_TRUE.SetPxPyPzE(reader->GetValue("Lambda0_TRUEP_X"),
					reader->GetValue("Lambda0_TRUEP_Y"),
					reader->GetValue("Lambda0_TRUEP_Z"),
					reader->GetValue("Lambda0_TRUEP_E"));
			proton_TRUE.SetPxPyPzE(reader->GetValue("pplus_TRUEP_X"),
					reader->GetValue("pplus_TRUEP_Y"),
					reader->GetValue("pplus_TRUEP_Z"),
					reader->GetValue("pplus_TRUEP_E"));
			pion_TRUE.SetPxPyPzE(reader->GetValue("piminus_TRUEP_X"),
					reader->GetValue("piminus_TRUEP_Y"),
					reader->GetValue("piminus_TRUEP_Z"),
					reader->GetValue("piminus_TRUEP_E"));

			int pcharge_TRUE = 1;
			if(reader->GetValue("pplus_TRUEID") < 0) pcharge_TRUE = -1;

			double tcosTheta_TRUE, tcosThetaL_TRUE, tcosThetaB_TRUE, tphiL_TRUE, tphiB_TRUE, tdphi_TRUE;
			DPHelpers::LbPsiRAngles(initialProton,Lb_TRUE,Jpsi_TRUE,Lambda_TRUE,mup_TRUE,mum_TRUE,proton_TRUE,pion_TRUE,pcharge_TRUE,
					tcosTheta_TRUE,tcosThetaL_TRUE,tcosThetaB_TRUE,tphiL_TRUE,tphiB_TRUE,tdphi_TRUE);
			cosTheta_TRUE = (float)tcosTheta_TRUE, cosThetaL_TRUE = (float)tcosThetaL_TRUE, cosThetaB_TRUE = (float)tcosThetaB_TRUE, phiL_TRUE = (float)tphiL_TRUE, phiB_TRUE = (float)tphiB_TRUE, dphi_TRUE = (float)tdphi_TRUE;
		}
	}
}


void RenameMass(TreeReader * reader, TTree * newTree, bool reset)
{
	static float mass1, mass2, mass3, mass4;
	if(reset)
	{
		newTree->Branch("Lb_MassCons",&mass1,"Lb_MassCons/F");
		newTree->Branch("Lb_MassConsLambda",&mass2,"Lb_MassConsLambda/F");
		newTree->Branch("Lb_MassConsJpsiLambda",&mass3,"Lb_MassConsJpsiLambda/F");
		newTree->Branch("J_psi_1S_MM",&mass4,"J_psi_1S_MM/F");
	}
	else
	{
		mass1 = reader->GetValue("Lb_MASS");
		mass2 = reader->GetValue("Lb_MASS");
		mass3 = reader->GetValue("Lb_MASS");
		mass4 = reader->GetValue("J_psi_1S_MASS");
	}
}



int main(int argc, char **argv)
{
	TString analysis = "Lb2Lmumu";
	bool MC = false;
	TString base = "/afs/cern.ch/work/p/pluca/weighted/Lmumu/";

	if(argc > 1)
	{
		string arg = argv[1];
		if(arg == "MC") MC = true;
	}

	vector< string > novar;
	novar.push_back("Lb_MassCons");
	novar.push_back("Lb_MassConsLambda");
	novar.push_back("Lb_MassConsJpsiLambda");
	novar.push_back("cosTheta");
	novar.push_back("cosThetaL");
	novar.push_back("cosThetaB");
	novar.push_back("phiL");
	novar.push_back("phiB");
	novar.push_back("dphi");
	novar.push_back("cosTheta_TRUE");
	novar.push_back("cosThetaL_TRUE");
	novar.push_back("cosThetaB_TRUE");
	novar.push_back("phiL_TRUE");
	novar.push_back("phiB_TRUE");
	novar.push_back("dphi_TRUE");

	TCut cutJpsi = CutsDef::cutJpsi;
	TCut cutMuMu = CutsDef::cutMuMu_veto;

	TreeReader* treeReader = new TreeReader("tree");

	TString namefile = base + "candLb";
	if(MC) namefile += "_MC";
	namefile += ".root";

	TFile * candFile = new TFile(namefile,"recreate");

	if(!MC) treeReader->AddFile(base+analysis+"_CL_NBweighted.root");
	else treeReader->AddFile(base+analysis+"_MC_Pythia8_NBweighted.root");
	treeReader->Initialize(novar,"except");

	Analysis * anaLbMuMu = new Analysis("Lb2Lmumu","Lb",treeReader,&cutMuMu);

	candFile->cd();
	TTree * candLbMuMu = anaLbMuMu->applyCuts(&addVariables);
	candLbMuMu->Write();
	string tnameMuMu = candLbMuMu->GetName();

	candFile->Close();
	candFile = TFile::Open(namefile,"update");
	TTree * singleCand_LbMuMu = anaLbMuMu->checkMultiple("weight",namefile,tnameMuMu,&randomKill);
	singleCand_LbMuMu->Write();



	if(MC)
	{
		treeReader = new TreeReader("tree");
		treeReader->AddFile(base+"Lb2JpsiL_MC_Pythia8_NBweighted.root");
		treeReader->Initialize(novar,"except");
	}
	Analysis * anaLbJpsi = new Analysis("Lb2JpsiL","Lb",treeReader,&cutJpsi);

	candFile->cd();
	TTree * candLbJpsi = anaLbJpsi->applyCuts(&addVariables);
	candLbJpsi->Write();
	string tnameJpsi = candLbJpsi->GetName();

	candFile->Close();
	candFile = TFile::Open(namefile,"update");	
	TTree * singleCand_LbJpsi = anaLbJpsi->checkMultiple("weight",namefile,tnameJpsi,&randomKill);
	singleCand_LbJpsi->Write();

	candFile->cd();
	TTree * candLbJpsi_reduced = anaLbJpsi->applyCuts(&addVariables,300);
	candLbJpsi_reduced->SetName("candLb2JpsiL_reduced");
	candLbJpsi_reduced->Write();



	if(MC)
	{
		candFile->cd();
		TCut jpsiSwap = cutJpsi + CutsDef::jpsiSwapID;
		TCut mumuSwap = cutMuMu + CutsDef::mumuSwapID;
		TTree * mumuSwapTree = anaLbMuMu->applyCuts(&mumuSwap, false,&addVariables);
		mumuSwapTree->SetName("candLmumuSwap");
		mumuSwapTree->Write();
		TTree * jpsiSwapTree = anaLbJpsi->applyCuts(&jpsiSwap, false, &addVariables);
		jpsiSwapTree->SetName("candJpsiLSwap");
		jpsiSwapTree->Write();

		TreeReader * KSReader = new TreeReader("tree");
		KSReader->AddFile(base+"Bd2JpsiKS_MC12_NBweighted.root");
		KSReader->Initialize(novar,"except");
		TCut cutBdLL = cutJpsi + CutsDef::LLcut;
		TCut cutBdDD = cutJpsi + CutsDef::DDcut;
		Analysis * KSAnalysis_LL = new Analysis("BdJpsiKS_LL","B0",KSReader,&cutBdLL);
		TTree *KSTree_LL = KSAnalysis_LL->applyCuts(&addVariables);		
		KSTree_LL->Write();
		Analysis * KSAnalysis_DD = new Analysis("BdJpsiKS_DD","B0",KSReader,&cutBdDD);
		TTree *KSTree_DD = KSAnalysis_DD->applyCuts(&addVariables);		
		KSTree_DD->Write();
		Analysis * KSAnalysis_all = new Analysis("BdJpsiKS","B0",KSReader,&cutJpsi);
		TTree *KSTree = KSAnalysis_all->applyCuts(&addVariables);		
		KSTree->Write();

		candFile->cd();
		TreeReader * KstmumuReader = new TreeReader("tree");
		KstmumuReader->AddFile(base+"Bu2Kstmumu_MC12_NBweighted.root");
		KstmumuReader->Initialize(novar,"except");
		Analysis * KstmumuAnalysis = new Analysis("BuKstmumu","B0",KstmumuReader,&cutMuMu);
		TTree *KstmumuTree = KstmumuAnalysis->applyCuts(&addVariables);
		KstmumuTree->Write();

		candFile->cd();
		TreeReader * KSmumuReader = new TreeReader("tree");
		KSmumuReader->AddFile(base+"Bd2KSmumu_MC12_NBweighted.root");
		KSmumuReader->Initialize(novar,"except");
		Analysis * KSmumuAnalysis = new Analysis("BdKSmumu","B0",KSmumuReader,&cutMuMu);
		TTree *KSmumuTree = KSmumuAnalysis->applyCuts(&addVariables);
		KSmumuTree->Write();

		candFile->cd();
		TreeReader * JpsiGenReader = new TreeReader("tree");
		JpsiGenReader->AddFile("/afs/cern.ch/work/k/kreps/public/LbLMuMuAna/generatorLevel/LbJpsiLGenOnlyDaughInAccForRadiativeTail.root");
		JpsiGenReader->Initialize();
		TCut JpsiTailCut = "TMath::Power(J_psi_1S_MASS/1000,2) < 8 && Lb_MASS > 5300 && Lambda0_MASS > 1105 && Lambda0_MASS < 1125";
		Analysis * JpsiTailAnalysis = new Analysis("JpsiTail","Lb",JpsiGenReader,&JpsiTailCut);
		TTree *JpsiTailTree = JpsiTailAnalysis->applyCuts(&RenameMass, 0.1);
		JpsiTailTree->Write();
	}

	candFile->Close();
	delete candFile;

	return 0;
}





