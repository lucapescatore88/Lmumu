#ifndef ADD_LB_ANGLES_HPP
#define ADD_LB_ANGLES_HPP


#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooStats/SPlot.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooCBShape.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooChi2Var.h"
#include "RooDataHist.h"
#include "RooAbsData.h"
#include "RooProdPdf.h"
#include "RooPolynomial.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "RooAbsPdf.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooConstVar.h"
#include "RooNumConvPdf.h"
#include "RooFFTConvPdf.h"
#include "RooKeysPdf.h"
#include "RooHistPdf.h"

#include "TCanvas.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TEntryList.h"
#include "TCanvas.h"
#include "TBranch.h"
#include "TCut.h"
#include "TMath.h"
#include "TIterator.h"

#include <vector>
#include <sstream>
#include <iostream>
#include <string>
#include <time.h>
#include <iomanip>

#include "ReadTree_comp.hpp"
#include "DPHelpers.hpp"
#include "general_functions.hpp"



void addAngles(TreeReader * reader, TTree * newTree, bool reset)
{
	static float const0, const1, const2;
	static float cosTheta, cosThetaL, cosThetaB, phiL, phiB, dphi;

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
	}
	else
	{
		TLorentzVector initialProton(0.,0.,3500000.,TMath::Sqrt(3500000.*3500000.+938.*938.));
		
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
	
		DPHelpers::LbPsiRAngles(initialProton,Lb,Jpsi,Lambda,mup,mum,proton,pion,pcharge,
					cosTheta,cosThetaL,cosThetaB,phiL,phiB,dphi);
	}
}


#endif
