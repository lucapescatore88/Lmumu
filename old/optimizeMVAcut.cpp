#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TMath.h"
#include "TEntryList.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TF1.h"
#include "TLatex.h"

#include "general_functions.hpp"
#include "Lb_cuts.hpp"

using namespace std;

int main(int argc, char** argv)
{
	/*
	Options:
	-f [sig file] [bkg file]
	-c [cut]
	-s [steps]
	-w [name weight]
	-t [analysis type] //DD or LL
	*/
	
	TString weight = "MCnorm*lifeTimeW*(lifeTimeW > 0)*physRate_pol0*(physRate_pol0 > 0)";
	TString analysis = "Lmumu";
	TString type = "all";
	TString namefilesig = "/afs/cern.ch/work/p/pluca/weighted/Lmumu/Lb2Lmumu_MC_Pythia8_NBweighted.root";
	TString namefilebkg = "/afs/cern.ch/work/p/pluca/weighted/Lmumu/Lb2Lmumu_CL_NBweighted.root";
	TString nameweight = "weight";
	int nsteps = 100;
	double defCut = 0;
	
	
	if(argc > 1)
	{
		for(int a = 1; a < argc; a++)
		{
			string arg = argv[a];
			string str = arg.substr(2,arg.length()-2);
			
			if(arg.find("-f") != string::npos && argc >= (a+3))
			{
				namefilesig = argv[a+1];
				namefilebkg = argv[a+2];
			}
			if(arg.find("-c") != string::npos)
			{
				defCut = ((TString)str).Atof();
				nsteps = 0;
			}
			if(arg.find("-s") != string::npos) nsteps = ((TString)str).Atoi();
			if(arg.find("-w") != string::npos) nameweight = str;
			if(arg.find("-t") != string::npos)
			{
				string type = str;
				analysis += ("_" + type);
			}
		}
	}

	TFile * fsig = TFile::Open(namefilesig);
	TTree * treeSigTmp = (TTree *)fsig->Get("tree");
	TTree * treeSig = (TTree *)treeSigTmp->Clone("treeSig");
	TFile * fbkg = TFile::Open(namefilebkg);
	TTree * treeBkgTmp = (TTree *)fbkg->Get("tree");
	TTree * treeBkg = (TTree *)treeBkgTmp->Clone("treeBkg");
	
	TCut sideBandCut = "Lb_MassConsLambda_M[0] > 6000";
	TCut trigAndJpsiCut = CutsDef::TrigPassed + CutsDef::avoidJpsiCut;
	if(type == "DD") trigAndJpsiCut += (TCut)"pplus_TRACK_Type == 3";
	else if(type == "LL") trigAndJpsiCut += (TCut)"pplus_TRACK_Type == 5";
	
	TCanvas *c0 = new TCanvas();
	treeBkg->Draw("Lb_MassConsLambda_M[0]>>hsideband(100,5900,7100)",trigAndJpsiCut + sideBandCut,"E");
	TH1 * hh = (TH1*)gPad->GetPrimitive("hsideband");
	TF1 * f = new TF1("f","[0]*TMath::Exp([1]*x)",5400,8000);
	//f->SetParameter(0,100);
	f->SetParameter(0,1000);
	f->SetParameter(1,-0.005);
	hh->Fit(f,"","",6100,6900);
	double bkgNorm = f->Integral(5619.4 - 30.,5619.4 + 30.) / f->Integral(6000,7000);
	c0->Print("Lbsideband.pdf");
	
	double JpsimumuBr = 0.0593;
	double Lb2LmumuBr_Over_Lb2LjpsiBr = 1.54e-3; //From our Lmumu article
	
	/*
	double sigma_pp2bbbarX = 284.e9; //(femtobarn => microbarn)
    double lumi = luminosity(namefilebkg);
	double Lb2LjpsiBr_timesfL = 5.8e-5; //PDG gives this quantity
	double Lb2LmumuBr_timesfL = Lb2LmumuBr_Over_Lb2LjpsiBr * Lb2LjpsiBr_timesfL;
	double LppiBr = 0.639;
	double acc = 0.1875;
	
	TreeReader * MCreader = new TreeReader("tree",namefilesig);
	MCreader->GetEntriy(1);
	double NGen = MCreader->GetValue("NGen");
	
	double norm_lumi = 1.1 * sigma_pp2bbbarX * lumi * Lb2LmumuBr_timesfL * acc / (double) NGen;
	*/
	
	TCut lowQ2 = "TMath::Power(J_psi_1S_MM/1000,2) < 8";
	TCut highQ2 = "TMath::Power(J_psi_1S_MM/1000,2) > 15";
	
	double nJpsiMuMuAfterStrip = 23419;
	treeSig->Draw("Lb_MM>>hhTot",weight+"*("+ (TString)(trigAndJpsiCut + CutsDef::mumuTrueID + highQ2) + ")");
	TH1 *hhTot = (TH1 *)gPad->GetPrimitive("hhTot");
	double nTotMCLmumu = hhTot->Integral();
	double norm = ( nJpsiMuMuAfterStrip * Lb2LmumuBr_Over_Lb2LjpsiBr / JpsimumuBr ) / nTotMCLmumu;
	treeSig->Draw("Lb_MM>>hhTotLow",weight+"*("+ (TString)(trigAndJpsiCut + CutsDef::mumuTrueID + lowQ2) + ")");
	TH1 *hhTotLow = (TH1 *)gPad->GetPrimitive("hhTotLow");
	double nTotLowMCLmumu = hhTotLow->Integral();
	double ratioLowOverHigh = nTotMCLmumu / nTotLowMCLmumu;
	double normLow = ( nJpsiMuMuAfterStrip * Lb2LmumuBr_Over_Lb2LjpsiBr / ( JpsimumuBr * ratioLowOverHigh) ) / nTotLowMCLmumu;
	
	cout << "Signal norm = " << norm << ",   Bkg norm = " << bkgNorm << endl;
		
	double optimalWlowQ2 = defCut;
	double optimalWhighQ2 = defCut;

	if(defCut == 0)
	{
		optimalWlowQ2 = optimizeCut((const char*)(analysis+"_lowQ2"), "Lb", treeSig, treeBkg, trigAndJpsiCut + lowQ2, CutsDef::mumuTrueID, sideBandCut, normLow, bkgNorm, weight, nsteps, "significance",(const char *)nameweight);
		optimalWhighQ2 = optimizeCut((const char*)(analysis+"_highQ2"), "Lb", treeSig, treeBkg, trigAndJpsiCut + highQ2, CutsDef::mumuTrueID, sideBandCut, norm, bkgNorm, weight, nsteps, "significance", (const char *)nameweight);
	}
	else cout << "Cut high/low Q2 = " << defCut << endl;


	treeSig->Draw("Lb_MM>>hhLowQ2",weight+"*("+(TString)(trigAndJpsiCut + CutsDef::mumuTrueID + lowQ2) + ")");
	TH1 *hhLowQ2 = (TH1 *)gPad->GetPrimitive("hhLowQ2");
	double SallLowQ2 = norm * hhLowQ2->Integral();
	treeSig->Draw("Lb_MM>>hhHighQ2",weight+"*("+(TString)(trigAndJpsiCut + CutsDef::mumuTrueID + highQ2) + ")");
	TH1 *hhHighQ2 = (TH1 *)gPad->GetPrimitive("hhHighQ2");
	double SallHighQ2 = norm * hhHighQ2->Integral();
	
	treeSig->Draw("Lb_MM>>hh0",weight+"*("+(TString)(CutsDef::mumuTrueID) + ")");
	TH1 *hh0 = (TH1 *)gPad->GetPrimitive("hh0");
	double Snotrig = norm * hh0->Integral();

	treeSig->Draw("Lb_MM>>hh5",weight+"*("+(TString)(CutsDef::TrigPassed + CutsDef::mumuTrueID) + ")");
	TH1 *hh5 = (TH1 *)gPad->GetPrimitive("hh5");
	double StrigCut = norm * hh5->Integral();
	
	TString selectLowQ2((TString)nameweight + Form(" > %e",optimalWlowQ2) +" && " + (TString)(trigAndJpsiCut + lowQ2));
	TString selectHighQ2((TString)nameweight + Form(" > %e",optimalWhighQ2) + " && " + (TString)(trigAndJpsiCut + highQ2));
	TString selectLowQ2MC(weight+"*("+(TString)nameweight + Form(" > %e",optimalWlowQ2) +" && " + (TString)(trigAndJpsiCut  + CutsDef::mumuTrueID + lowQ2) + ")");
        TString selectHighQ2MC(weight+"*("+(TString)nameweight + Form(" > %e",optimalWhighQ2) + " && " + (TString)(trigAndJpsiCut  + CutsDef::mumuTrueID + highQ2) + ")");

	treeSig->Draw("Lb_MM>>hhLowQ2pass", selectLowQ2MC);
	TH1 *hhLowQ2pass = (TH1 *)gPad->GetPrimitive("hhLowQ2pass");
	double SpassLowQ2 = norm * hhLowQ2pass->Integral();
	treeSig->Draw("Lb_MM>>hhHighQ2pass", selectHighQ2MC);
	TH1 *hhHighQ2pass = (TH1 *)gPad->GetPrimitive("hhHighQ2pass");
	double SpassHighQ2 = norm * hhHighQ2pass->Integral();
	
	treeBkg->Draw("Lb_MM>>hhLowQ2bkg", (TCut)selectLowQ2 + sideBandCut);
	TH1 *hhLowQ2bkg = (TH1 *)gPad->GetPrimitive("hhLowQ2bkg");
	double BpassLowQ2 = bkgNorm * hhLowQ2bkg->Integral();
	treeBkg->Draw("Lb_MM>>hhHighQ2bkg", (TCut)selectHighQ2 + sideBandCut);
	TH1 *hhHighQ2bkg = (TH1 *)gPad->GetPrimitive("hhHighQ2bkg");
	double BpassHighQ2 = bkgNorm * hhHighQ2bkg->Integral();
	
	double SigLowQ2 = SpassLowQ2 / TMath::Sqrt(BpassLowQ2 + SpassLowQ2);
	double SigHighQ2 = SpassHighQ2 / TMath::Sqrt(BpassHighQ2 + SpassHighQ2);
	
	cout << "Trigger efficiency = " << StrigCut / Snotrig * 100 << "%" << endl;
	cout << "Jpsi cut (on Lmumu) efficiency = " << (SallLowQ2+SallHighQ2) / StrigCut * 100 << "%" << endl;
	cout << "Predicted signal events lowQ2 = " << SpassLowQ2 << ",  bkg events = " << BpassLowQ2 << " (in |mLb - m| < 30)\n   =>   Significance = " << SigLowQ2 << endl;
	cout << "Predicted signal events highQ2 = " << SpassHighQ2 << ",  bkg events = " << BpassHighQ2 << " (in |mLb - m| < 30)\n   =>   Significance = " << SigHighQ2 << endl;
	
	return 0;
}
