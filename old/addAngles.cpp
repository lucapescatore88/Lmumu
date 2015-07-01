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

		if(reader->HasVar("Lb_TRUEP_X"));
		{
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
	}
}



int main(int argc, char **argv)
{
	TString analysis = "Lb2Lmumu";

	vector< string > novar;
	novar.push_back("cosTheta");
	novar.push_back("cosThetaL");
	novar.push_back("cosThetaB");
	novar.push_back("phiL");
	novar.push_back("phiB");
	novar.push_back("dphi");

	TreeReader* treeReader = new TreeReader("tree");
	if(argc>3) treeReader = new TreeReader(argv[3]);

	TString namefile = argv[1];

	TFile * candFile = new TFile(namefile,"recreate");

	treeReader->AddFile(argv[2]);
	treeReader->Initialize(novar,"except");

	TCut noCut = "";
	Analysis * anaLbMuMu = new Analysis("Lb2Lmumu","Lb",treeReader,&noCut);

	TTree * candLbMuMu = anaLbMuMu->applyCuts(&addVariables);
	if(argc>3) candLbMuMu->SetName(argv[3]);
	else candLbMuMu->SetName("tree");
	candLbMuMu->Write();

	candFile->Close();
	delete candFile;

	return 0;
}





