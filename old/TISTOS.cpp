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


int main(int argc, char **argv)
{
	TString base = "/afs/cern.ch/work/p/pluca/weighted/Lmumu/";

	TCut cutJpsi = "TMath::Abs(J_psi_1S_MM - 3096.9) < 30 && TMath::Abs(Lb_MassConsLambda_M[0] - 5621) < 30";

	TFile * f = TFile::Open(base+"Lb2Lmumu_CL_NBweighted.root");	
	TTree * candLbJpsi = (TTree*)f->Get("tree");

	TCut L0Passed = "(Lb_L0MuonDecision_TIS || Lb_L0DiMuonDecision_TIS)";
	TCut Hlt1Passed = "(Lb_Hlt1TrackAllL0Decision_TIS || Lb_Hlt1TrackMuonDecision_TIS || Lb_Hlt1MuTrackDecision_TIS || Lb_Hlt1DiMuonHighMassDecision_TIS)";
	TCut Hlt2Passed = "(Lb_Hlt2Topo2BodyBBDTDecision_TIS || Lb_Hlt2Topo3BodyBBDTDecision_TIS || Lb_Hlt2Topo4BodyBBDTDecision_TIS || Lb_Hlt2TopoMu2BodyBBDTDecision_TIS || Lb_Hlt2TopoMu3BodyBBDTDecision_TIS || Lb_Hlt2TopoMu4BodyBBDTDecision_TIS || Lb_Hlt2SingleMuonDecision_TIS || Lb_Hlt2DiMuonDetachedDecision_TIS)";
	TCut TIStriggers = (TCut)( (TString)L0Passed + " && " + (TString)Hlt1Passed + " && " + (TString)Hlt2Passed );



	TCut baseCut = CutsDef::MVAcut;
	int nbins = 3;
	double bins[] = {3000,1.5e4,3.5e4,2e5};
	TString var = "muminus_P";

	for(int i = 0; i <= nbins; i++)
	{

		if(i>0) 
		{
			baseCut = CutsDef::MVAcut;
			baseCut += (TCut)(var + Form(">%f && ",bins[i-1]) + var + Form("<%f ",bins[i]));
			baseCut.Print();
		}
		candLbJpsi->Draw("Lb_MM>>hTOS",cutJpsi + CutsDef::TrigPassed + TIStriggers + baseCut,"E");
		TH1F * hTOS = (TH1F*)gPad->GetPrimitive("hTOS");

		candLbJpsi->Draw("Lb_MM>>hTIS",cutJpsi + TIStriggers + baseCut,"E");
		TH1F * hTIS = (TH1F*)gPad->GetPrimitive("hTIS");

		double ratio = hTOS->GetEntries() / hTIS->GetEntries();
		cout << ratio << " +/-" << TMath::Sqrt(ratio*(1-ratio)/hTIS->GetEntries()) << endl;

		TFile * fMC = TFile::Open(base+"Lb2JpsiL_MC_Pythia8_NBweighted.root");
		TTree * candLbJpsiMC = (TTree*)fMC->Get("tree");


		TString weight = "MCnorm*(lifeTimeW > 0)*(physRate_pol0_noDecay > 0)*lifeTimeW*physRate_pol0"; 
		candLbJpsiMC->Draw("Lb_MM>>hMCTOS",weight + "*(" + (TString)( baseCut + CutsDef::baseCutJpsi + CutsDef::massCutUnblinded + cutJpsi + CutsDef::TrigPassed) + ")","E");
		TH1F * hMCTOS = (TH1F*)gPad->GetPrimitive("hMCTOS");
		candLbJpsiMC->Draw("Lb_MM>>hMCAll",weight + "*(" + (TString)( baseCut + CutsDef::baseCutJpsi + CutsDef::massCutUnblinded + cutJpsi) + ")","E");
		TH1F * hMCAll = (TH1F*)gPad->GetPrimitive("hMCAll");

		double ratioMC = hMCTOS->Integral() / hMCAll->Integral();
		cout << ratioMC << " +/-" << TMath::Sqrt(ratioMC*(1-ratioMC)/hMCAll->GetEntries()) << endl;	
	}

	return 0;
}





