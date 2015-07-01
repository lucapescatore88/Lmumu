#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraph.h"

int main()
{
	TFile * trainFile = TFile::Open("/afs/cern.ch/work/p/pluca/weighted/Lmumu/trainingSamples.root");
	TTree * sigTestS = (TTree *)trainFile->Get("sigTestSample");
	TTree * sigTrainS = (TTree *)trainFile->Get("sigTrainSample");
	TTree * bkgTrainS = (TTree *)trainFile->Get("bkgTrainSample");
	TTree * bkgTestS = (TTree *)trainFile->Get("bkgTestSample");

	trainFile = TFile::Open("/afs/cern.ch/work/p/pluca/weighted/Lmumu/Lb2Lmumu_CL_NBweighted.root");
	TTree * bkgTestS_prof = (TTree *)trainFile->Get("tree");

	TCanvas * c =  new TCanvas();

	gStyle->SetOptStat(0);

	bkgTestS_prof->Draw("weight:Lb_MM>>bkgProf(100,5700,7000)","TMath::Abs(J_psi_1S_MM-3096)>50","prof");
	TH1F * bkgProf = (TH1F*)gPad->GetPrimitive("bkgProf");
	bkgProf->GetYaxis()->SetTitle("<NNout>");
	bkgProf->GetXaxis()->SetTitle("m(#Lambda#mu#mu) (MeV/c^{2})");
	bkgProf->SetTitle(0);
	bkgProf->SetMinimum(0);
	bkgProf->SetMaximum(1);
	bkgProf->Draw();
	c->Print("NNout_profile_vs_LbMM_bkgData.pdf");

	sigTestS->Draw("weight:Lb_MM>>sigProf(100,5600,5630)","","prof");
	TH1F * sigProf = (TH1F*)gPad->GetPrimitive("sigProf");
	sigProf->GetYaxis()->SetTitle("<NNout>");
	sigProf->GetXaxis()->SetTitle("m(#Lambda#mu#mu) (MeV/c^{2})");
	sigProf->SetTitle(0);
	sigProf->SetMinimum(0);
	sigProf->SetMaximum(1);
	sigProf->Draw();
	c->Print("NNout_profile_vs_LbMM_MCsignal.pdf");

	bkgTestS->Draw("weight>>hBTest(50,0,1)","","E");
	TH1F * bkgTestW = (TH1F*)gPad->GetPrimitive("hBTest");
	bkgTrainS->Draw("weight>>hBTrain(50,0,1)","","");
	TH1F * bkgTrainW = (TH1F*)gPad->GetPrimitive("hBTrain");
	sigTrainS->Draw("weight>>hSTrain(50,0,1)","","");
	TH1F * sigTrainW = (TH1F*)gPad->GetPrimitive("hSTrain");
	sigTestS->Draw("weight>>hSTest(50,0,1)","","E");
	TH1F * sigTestW = (TH1F*)gPad->GetPrimitive("hSTest");

	c->SetLogy();

	bkgTrainW->SetMarkerColor(2);
	bkgTestW->SetMarkerColor(2);
	sigTrainW->SetMarkerColor(4);
	sigTestW->SetMarkerColor(4);
	bkgTrainW->SetMarkerSize(0.9);
	sigTrainW->SetMarkerSize(0.9);
	sigTestW->SetMarkerSize(0.9);
	bkgTestW->SetMarkerSize(0.9);
	bkgTrainW->SetMarkerStyle(20);
	sigTrainW->SetMarkerStyle(22);
	sigTestW->SetMarkerStyle(22);
	bkgTestW->SetMarkerStyle(20);

	bkgTrainW->SetFillColor(2);
	bkgTestW->SetFillColor(2);
	sigTrainW->SetFillColor(4);
	sigTestW->SetFillColor(4);
	bkgTrainW->SetFillStyle(3353);
	sigTrainW->SetFillStyle(3335);
	sigTestW->SetFillStyle(1001);
	bkgTestW->SetFillStyle(1001);
	//bkgTestW->SetFillColorAlpha(2, 0.35);
	//sigTestW->SetFillColorAlpha(2, 0.35);

	bkgTrainW->Scale(1./bkgTrainW->Integral());
	sigTrainW->Scale(1./sigTrainW->Integral());
	sigTestW->Scale(1./sigTestW->Integral());
	bkgTestW->Scale(1./(bkgTestW->Integral()));

	//sigTrainW->GetYaxis()->SetTitle("");
	sigTestW->GetXaxis()->SetTitle("NNout");
	sigTestW->SetTitle(0);

	sigTestW->Draw("");
	bkgTestW->Draw("same");
	bkgTrainW->Draw("same");
	sigTrainW->Draw("same");
	c->Print("TrainAndTest.pdf");

	TFile * hoptFile = TFile::Open("optimise_Lmumu_highQ2.root");
	TGraph * highq2ROC = (TGraph *)hoptFile->Get("ROC");
	TFile * loptFile = TFile::Open("optimise_Lmumu_lowQ2.root");
	TGraph * lowq2ROC = (TGraph *)loptFile->Get("ROC");
	
	highq2ROC->SetMarkerSize(0.4);
	lowq2ROC->SetMarkerSize(0.4);
	highq2ROC->SetMarkerStyle(21);
	lowq2ROC->SetMarkerStyle(20);
	highq2ROC->SetMarkerColor(1);
	lowq2ROC->SetMarkerColor(2);

	c->SetLogy(0);

	lowq2ROC->GetYaxis()->SetTitle("sig. eff.");
	lowq2ROC->GetXaxis()->SetTitle("bkg. rej.");
	lowq2ROC->Draw("AP");
	highq2ROC->Draw("Psame");
	c->Print("ROC.pdf");


	
}



