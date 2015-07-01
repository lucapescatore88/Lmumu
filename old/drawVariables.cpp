#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TF1.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TF2.h>
#include <TChain.h>
#include <TString.h>
#include <TRandom.h>
#include <TMath.h>
#include <TLegend.h>
#include <TCut.h>
#include <TString.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TPaveStats.h>
#include <TEfficiency.h>

#include <limits>
#include <vector>
#include <sstream>

#include "ReadTree_comp.hpp"
#include "analyser.hpp"
#include "Lb_cuts.hpp"
#include "general_functions.hpp"

using namespace std;


int main(int argc, char **argv)
{
	TString type = "all";	
	if(argc > 1)
	{
		for(int a = 1; a < argc; a++)
		{
			string arg = argv[a];
			string str = arg.substr(2,arg.length()-2);
			if(arg.find("-t") != string::npos) type = str;
		}
	}

	cout << "Analysing " << type << endl;

	TCut baseCut = "";//CutsDef::TrigPassed;//"TMath::Abs(J_psi_1S_MM - 3096) < 30"+CutsDef::cutJpsi;
	TCut dataCut = "TMath::Abs(J_psi_1S_MM - 3096) < 30";
	TCut MCcut = "";
	if(type == "DD") baseCut += CutsDef::DDcut;
	else if(type == "LL") baseCut += CutsDef::LLcut;



	//Attach data files	

	TCanvas *c = new TCanvas("c","");



	TString MCFileName = "/afs/cern.ch/work/p/pluca/weighted/Lmumu/Lb2Lmumu_MC_Pythia8_NBweighted.root";
	TString MCJpsiFileName = "/afs/cern.ch/work/p/pluca/weighted/Lmumu/Lb2JpsiL_MC_Pythia8_NBweighted.root";
	TString DataFileName = "/afs/cern.ch/work/p/pluca/weighted/Lmumu/Lb2Lmumu_CL_NBweighted.root";
	TString MCFileGeomName = "/afs/cern.ch/work/p/pluca/weighted/Lmumu/Lb2Lmumu_geomMC_Pythia8_NBweighted.root";
	TString MCFileJpsiGeomName = "/afs/cern.ch/work/p/pluca/weighted/Lmumu/Lb2JpsiL_geomMC_Pythia8_NBweighted.root";
	TFile * MCFile = TFile::Open(MCFileName);
	TTree * MCtree = (TTree *)MCFile->Get("tree");
	TFile * MCJpsiFile = TFile::Open(MCJpsiFileName);
	TTree * MCJpsitree = (TTree *)MCJpsiFile->Get("tree");
	TFile * DataFile = TFile::Open(DataFileName);
	TTree * dataTree = (TTree *)DataFile->Get("tree");
	TFile * MCFileGeom = TFile::Open(MCFileGeomName);
	TTree * MCtreeGeom = (TTree *)MCFileGeom->Get("MCtree");
	TFile * MCFileJpsiGeom = TFile::Open(MCFileJpsiGeomName);
	TTree * MCtreeJpsiGeom = (TTree *)MCFileJpsiGeom->Get("MCtree");
	
	/*
	MCtree->Draw("TMath::Power(J_psi_1S_MM/1000,2)>>hq2nw","","E");
	TH1F * hq2nw = (TH1F*)gPad->GetPrimitive("hq2nw");
	hq2nw->Scale(1./hq2nw->Integral());
	hq2nw->SetMarkerColor(1);
	hq2nw->SetMarkerSize(0.8);
	hq2nw->SetMarkerStyle(20);
	MCtree->Draw("TMath::Power(J_psi_1S_MM/1000,2)>>hq2w","physRate_pol0","E");
	TH1F * hq2w = (TH1F*)gPad->GetPrimitive("hq2w");
	hq2w->Scale(1./hq2w->Integral());
	hq2w->SetMarkerColor(2);
	hq2w->SetMarkerSize(0.8);
	hq2w->SetMarkerStyle(20);
	MCJpsitree->Draw("TMath::Power(J_psi_1S_MM/1000,2)>>hq2nw_JPsi","","E");
	TH1F * hq2nw_JPsi = (TH1F*)gPad->GetPrimitive("hq2nw_JPsi");
	hq2nw_JPsi->Scale(1./hq2nw_JPsi->Integral());
	hq2nw_JPsi->SetMarkerColor(1);
	hq2nw_JPsi->SetMarkerSize(0.8);
	hq2nw_JPsi->SetMarkerStyle(21);
	MCJpsitree->Draw("TMath::Power(J_psi_1S_MM/1000,2)>>hq2w_JPsi","physRate_pol0","E");
	TH1F * hq2w_JPsi = (TH1F*)gPad->GetPrimitive("hq2w_JPsi");
	hq2w_JPsi->Scale(1./hq2w_JPsi->Integral());
	hq2w_JPsi->SetMarkerColor(2);
	hq2w_JPsi->SetMarkerSize(0.8);
	hq2w_JPsi->SetMarkerStyle(21);

	hq2w->Draw();
	hq2nw->Draw("same");
	hq2nw_JPsi->Draw("same");
	hq2w_JPsi->Draw("same");
	c->Print("q2_compare.pdf");
*/

	dataTree->Draw("Lb_MassConsJpsiAndLambda_M>>hsideband(100,5500,7000)",dataCut+baseCut,"E");
	TH1 * hh = (TH1*)gPad->GetPrimitive("hsideband");
	TF1 * f = new TF1("f","[0]*TMath::Exp([1]*x)+[2]*TMath::Gaus(x,[3],[4])",5000,7100);
	TF1 * exp = new TF1("exp","[0]*TMath::Exp([1]*x)",5000,7100);
	f->SetParameter(0,5000);
	f->SetParameter(1,-0.005);
	if(type=="DD") f->SetParameter(2,10000);
	else f->SetParameter(2,1000);
	f->SetParameter(3,5621);
	f->SetParameter(4,7);
	hh->Fit(f,"","",5570,7000);
	double pars[5];
	f->GetParameters(pars);
	exp->SetParameter(0,pars[0]);
	exp->SetParameter(1,pars[1]);
	double factor = exp->Integral(5605,5635) / exp->Integral(5800,7000);
	cout << factor << endl;
	c->Print("side_pre"+type+".pdf");

	

	gStyle->SetOptStat(0);

	int nvars = 26;
	string vars[] = {"Lambda0_ENDVERTEX_Z","muplus_PIDmu","muplus_ProbNNmu","muplus_IPCHI2_OWNPV","J_psi_1S_IPCHI2_OWNPV","Lb_TAU",
		"Lb_P","Lb_PT","Lambda0_PT","Lambda0_IPCHI2_OWNPV","Lambda0_FD_OWNPV","Lb_DIRA_OWNPV","muplus_TRACK_CHI2NDOF",
		"TMath::Min(muplus_IPCHI2_OWNPV,muminus_IPCHI2_OWNPV)","TMath::Max(muplus_IPCHI2_OWNPV,muminus_IPCHI2_OWNPV)",
		"Lb_MassConsLambda_chi2[0]/Lb_MassConsLambda_nDOF[0]","pplus_IPCHI2_OWNPV","piminus_IPCHI2_OWNPV","weight",
		"piminus_P","piminus_PT","cosThetaB","cosThetaL","Lambda0_P","pplus_P","pplus_PT"};
	double min[] = {-200, -4, 0, 0,    0,    0,    0,     0,   0,   0,   0,    0.9999, 0, 0,   0,   0,   0,    0,    0, 1500, 0,  0, 0, 7e3, 0, 0};
	double max[] = {4e3,  14, 1, 50e3, 20e3, 0.05, 400e3, 5e4, 10e3, 1e3, 4000, 1.,     3, 5e4, 5e4, 120, 3e3, 10e3, 1, 15e3, 5e3, 1, 1, 12e4, 5e5, 1e4};
	bool log[] = {false,false,false,true,true,true,false,true,false,true,false,true,false,true,true,true,true,true,true,false,true,false,false,true,true,true};

	//Begin computing

	TLegend * leg = new TLegend(0.7,0.7,0.99,0.99);

	for(int v = 0; v < nvars; v++)
	{
		c->SetLogy(log[v]);

		TString namevar = (TString)vars[v];
		namevar.ReplaceAll(":","");
		namevar.ReplaceAll(",","");
		namevar.ReplaceAll("(","");
		namevar.ReplaceAll(")","");
		namevar.ReplaceAll("[","");
		namevar.ReplaceAll("]","");
		namevar.ReplaceAll("/","");
		namevar.ReplaceAll("TMath","");

		TString binning = Form("(40,%e,%e)",min[v],max[v]);
		TString nameMCh = "hMC"+namevar+Form("_%i",v);
		TString nameMCJpsih = "hMCJpsi"+namevar+Form("_%i",v);
		TString nameMC_Prw_h = "hMCPrw"+namevar+Form("_%i",v);
		TString nameMC_q2rw_h = "hMCq2rw"+namevar+Form("_%i",v);
		TString nameMC_fullrw_h = "hMCfullrw"+namevar+Form("_%i",v);
		TString nameDatah = "hData"+namevar+Form("_%i",v);
		TString nameSideh = "hSide"+namevar+Form("_%i",v);
		TH1F * hMC = NULL, * hMCJpsi = NULL, * hMC_q2rw = NULL, * hMC_Prw = NULL, * hMC_fullrw = NULL;
		if((MCcut+baseCut)!="")
		{
			MCtree->Draw(vars[v]+">>"+nameMCh+binning,"MCnorm*lifeTimeW*("+(TString)(MCcut+baseCut)+")","E");
			hMC = (TH1F*)gPad->GetPrimitive(nameMCh);
			MCtree->Draw(vars[v]+">>"+nameMC_q2rw_h+binning,"physRate_polp006*MCnorm*lifeTimeW*("+(TString)(MCcut+baseCut)+")","E");
			hMC_q2rw = (TH1F*)gPad->GetPrimitive(nameMC_q2rw_h);
			//MCtree->Draw(vars[v]+">>"+nameMC_Prw_h+binning,"Lb_weight*MCnorm*lifeTimeW*("+(TString)(MCcut+baseCut)+")","E");
			//hMC_Prw = (TH1F*)gPad->GetPrimitive(nameMC_Prw_h);
			//MCtree->Draw(vars[v]+">>"+nameMC_fullrw_h+binning,"Lb_weight*physRate_polp006*MCnorm*lifeTimeW*("+(TString)(MCcut+baseCut)+")","E");
			//hMC_fullrw = (TH1F*)gPad->GetPrimitive(nameMC_fullrw_h);
			MCtree->Draw(vars[v]+">>"+nameMC_Prw_h+binning,"pt_weight*MCnorm*lifeTimeW*("+(TString)(MCcut+baseCut)+")","E");
			hMC_Prw = (TH1F*)gPad->GetPrimitive(nameMC_Prw_h);
			MCtree->Draw(vars[v]+">>"+nameMC_fullrw_h+binning,"pt_weight*physRate_polp006*MCnorm*lifeTimeW*("+(TString)(MCcut+baseCut)+")","E");
			hMC_fullrw = (TH1F*)gPad->GetPrimitive(nameMC_fullrw_h);
			MCJpsitree->Draw(vars[v]+">>"+nameMCJpsih+binning,"MCnorm*lifeTimeW*("+(TString)(MCcut+baseCut)+")","E");
			hMCJpsi = (TH1F*)gPad->GetPrimitive(nameMCJpsih);
		}
		else
		{
			MCtree->Draw(vars[v]+">>"+nameMCh+binning,"MCnorm*lifeTimeW","E");
			hMC = (TH1F*)gPad->GetPrimitive(nameMCh);
			MCtree->Draw(vars[v]+">>"+nameMC_q2rw_h+binning,"physRate_polp006*MCnorm*lifeTimeW","E");
			hMC_q2rw = (TH1F*)gPad->GetPrimitive(nameMC_q2rw_h);
			//MCtree->Draw(vars[v]+">>"+nameMC_Prw_h+binning,"Lb_weight*MCnorm*lifeTimeW","E");
			//hMC_Prw = (TH1F*)gPad->GetPrimitive(nameMC_Prw_h);
			//MCtree->Draw(vars[v]+">>"+nameMC_fullrw_h+binning,"Lb_weight*physRate_polp006*MCnorm*lifeTimeW","E");
			//hMC_fullrw = (TH1F*)gPad->GetPrimitive(nameMC_fullrw_h);
			MCtree->Draw(vars[v]+">>"+nameMC_Prw_h+binning,"pt_weight*MCnorm*lifeTimeW","E");
			hMC_Prw = (TH1F*)gPad->GetPrimitive(nameMC_Prw_h);
			MCtree->Draw(vars[v]+">>"+nameMC_fullrw_h+binning,"pt_weight*physRate_polp006*MCnorm*lifeTimeW","E");
			hMC_fullrw = (TH1F*)gPad->GetPrimitive(nameMC_fullrw_h);
			MCJpsitree->Draw(vars[v]+">>"+nameMCJpsih+binning,"MCnorm*lifeTimeW","E");
			hMCJpsi = (TH1F*)gPad->GetPrimitive(nameMCJpsih);
		}	

		dataTree->Draw(vars[v]+">>"+nameDatah+binning,dataCut+baseCut+"Lb_MassConsJpsiAndLambda_M > 5605 && Lb_MassConsJpsiAndLambda_M < 5635","E");
		TH1F * hData = (TH1F*)gPad->GetPrimitive(nameDatah);
		dataTree->Draw(vars[v]+">>"+nameSideh+binning,dataCut+baseCut+"Lb_MassConsJpsiAndLambda_M > 5800 && Lb_MassConsJpsiAndLambda_M < 7000","E");
		TH1F * hSide = (TH1F*)gPad->GetPrimitive(nameSideh);

		hData->Add(hSide,-factor);

		hMC->SetMarkerStyle(20);
		hMC->SetMarkerColor(1);
		hMC->SetMarkerSize(0.8);
		hMC_fullrw->SetMarkerStyle(20);
		hMC_fullrw->SetMarkerColor(8);
		hMC_fullrw->SetMarkerSize(0.7);
		hMC_Prw->SetMarkerStyle(20);
		hMC_Prw->SetMarkerColor(6);
		hMC_Prw->SetMarkerSize(0.7);
		hMC_q2rw->SetMarkerStyle(20);
		hMC_q2rw->SetMarkerColor(3);
		hMC_q2rw->SetMarkerSize(0.7);
		hMCJpsi->SetMarkerStyle(20);
		hMCJpsi->SetMarkerColor(5);
		hMCJpsi->SetMarkerSize(0.7);


		hData->SetMarkerStyle(21);
		hData->SetMarkerColor(2);
		hData->SetMarkerSize(0.8);
		hSide->SetMarkerStyle(22);
		hSide->SetMarkerColor(4);
		hSide->SetMarkerSize(0.8);

		hMC->Scale(1./hMC->Integral());
		hMC_q2rw->Scale(1./hMC_q2rw->Integral());
		hMC_Prw->Scale(1./hMC_Prw->Integral());
		hMC_fullrw->Scale(1./hMC_fullrw->Integral());
		hMCJpsi->Scale(1./hMCJpsi->Integral());
		hData->Scale(1./hData->Integral());
		hSide->Scale(1./hSide->Integral());



		if(v==0)
		{
			leg->AddEntry(hData,"Data","P");
			leg->AddEntry(hSide,"Side band","P");
			leg->AddEntry(hMC,"MC","P");
			//leg->AddEntry(hMC_q2rw,"MC dec mod W","P");
			//leg->AddEntry(hMC_Prw,"MC Lb mom W","P");
			leg->AddEntry(hMC_fullrw,"MC fully W","P");
			//leg->AddEntry(hMCJpsi,"MC JpsiL","P");
		}

		hData->SetTitle(0);
		hData->GetXaxis()->SetTitle(vars[v].c_str());
		hData->Draw();
		hMC->Draw("same");
		hSide->Draw("same");
		hMC_fullrw->Draw("same");
		//hMC_q2rw->Draw("same");
		//hMC_Prw->Draw("same");
		//hMCJpsi->Draw("same");
		leg->Draw();
		c->Print("MC_data_comp/"+namevar+"_plot"+type+".pdf");
	}


	//--------------------------------------------------------------------------------



	TString binning = "(10,-1,1)";
	MCJpsitree->Draw("cosThetaB>>hMCLL"+binning,MCcut+baseCut+"pplus_TRACK_Type==3","E");
	TH1F * hMCLL = (TH1F*)gPad->GetPrimitive("hMCLL");
	dataTree->Draw("cosThetaB>>hDataLL"+binning,dataCut+baseCut+"pplus_TRACK_Type==3 && Lb_MassConsJpsiAndLambda_M > 5605 && Lb_MassConsJpsiAndLambda_M < 5635","E");
	TH1F * hDataLL = (TH1F*)gPad->GetPrimitive("hDataLL");
	dataTree->Draw("cosThetaB>>hSideLL"+binning,dataCut+baseCut+"pplus_TRACK_Type==3 && Lb_MassConsJpsiAndLambda_M > 5800 && Lb_MassConsJpsiAndLambda_M < 7000","E");
	TH1F * hSideLL = (TH1F*)gPad->GetPrimitive("hSideLL");

	hDataLL->Add(hSideLL,-factor);
	hMCLL->Scale(1./hMCLL->Integral());
	hDataLL->Scale(1./hDataLL->Integral());


	MCJpsitree->Draw("cosThetaB>>hMCDD"+binning,MCcut+baseCut+"pplus_TRACK_Type==5","E");
	TH1F * hMCDD = (TH1F*)gPad->GetPrimitive("hMCDD");
	dataTree->Draw("cosThetaB>>hDataDD"+binning,dataCut+baseCut+"pplus_TRACK_Type==5 && Lb_MassConsJpsiAndLambda_M > 5605 && Lb_MassConsJpsiAndLambda_M < 5635","E");
	TH1F * hDataDD = (TH1F*)gPad->GetPrimitive("hDataDD");
	dataTree->Draw("cosThetaB>>hSideDD"+binning,dataCut+baseCut+"pplus_TRACK_Type==5 && Lb_MassConsJpsiAndLambda_M > 5800 && Lb_MassConsJpsiAndLambda_M < 7000","E");
	TH1F * hSideDD = (TH1F*)gPad->GetPrimitive("hSideDD");

	hDataDD->Add(hSideDD,-factor);
	hMCDD->Scale(1./hMCDD->Integral());
	hDataDD->Scale(1./hDataDD->Integral());

	hDataDD->SetMarkerSize(0.8);
	hDataDD->SetMarkerColor(1);
	hDataDD->SetMarkerStyle(20);
	hMCDD->SetMarkerSize(0.8);
	hMCDD->SetMarkerColor(2);
	hMCDD->SetMarkerStyle(20);
	hDataLL->SetMarkerSize(0.8);
	hDataLL->SetMarkerColor(1);
	hDataLL->SetMarkerStyle(21);
	hMCLL->SetMarkerSize(0.8);
	hMCLL->SetMarkerColor(2);
	hMCLL->SetMarkerStyle(21);

	TLegend * leg2 = new TLegend(0.7,0.7,0.99,0.99);
	leg2->AddEntry(hDataLL,"Data LL","P");
	leg2->AddEntry(hDataDD,"Data DD","P");
	leg2->AddEntry(hMCLL,"MC LL","P");
	leg2->AddEntry(hMCDD,"MC DD","P");

	c->SetLogy(0);
	hDataDD->GetXaxis()->SetTitle("cos#theta_{#Lambda}");
	hDataLL->Draw();
	hDataDD->Draw("same");
	hMCDD->Draw("same");
	hMCLL->Draw("same");
	leg2->Draw("same");
	c->Print("costhetaB_distribs.pdf");
	c->SetLogy();

	TH1F * dataRatio = (TH1F *)hDataDD->Clone("dataRatio");
	dataRatio->GetXaxis()->SetTitle("cos#theta_{#Lambda}");
	dataRatio->Divide(hDataLL);
	TH1F * mcRatio = (TH1F *)hMCDD->Clone("mcRatio");
	mcRatio->GetXaxis()->SetTitle("cos#theta_{#Lambda}");
	mcRatio->Divide(hMCLL);

	TH1F * ratio = (TH1F *)dataRatio->Clone("ratio");
	ratio->Divide(mcRatio);
	ratio->GetXaxis()->SetTitle("cos#theta_{#Lambda}");
	ratio->Draw();
	c->SetLogy(0);
	c->Print("DD_Over_LL_MC_Over_Data_cosThetaB.pdf");

	TLegend * leg3 = new TLegend(0.7,0.7,0.99,0.99);
	leg3->AddEntry(dataRatio,"Data","P");
	leg3->AddEntry(mcRatio,"MC","P");

	dataRatio->SetMarkerSize(0.8);
	dataRatio->SetMarkerColor(1);
	dataRatio->SetMarkerStyle(20);
	mcRatio->SetMarkerSize(0.8);
	mcRatio->SetMarkerColor(2);
	mcRatio->SetMarkerStyle(21);
	mcRatio->Draw();
	dataRatio->Draw("same");
	leg3->Draw("same");
	c->Print("DD_Over_LL_MC_And_Data_cosThetaB.pdf");

	//-----------------------------------------------------------------------------

	TFile * outFile = new TFile("MC_reweight.root","recreate");

	double p[] = {1e4,6e4,7e4,8e4,9e4,10e4,12e4,15e4,100e4};
	double pt[] = {0.,2e3,3e3,4e3,5e3,6e3,7e3,8e3,10e3,100e3};

	TH2F * hMC2D = new TH2F("hMC2D","hMC2D",8,pt,7,p);
	TH2F * hData2D = new TH2F("hData2D","hData2D",8,pt,7,p);
	TH2F * hSide2D = new TH2F("hSide2D","hSide2D",8,pt,7,p);
	MCtree->Draw("Lb_P:Lb_PT>>hMC2D","MCnorm*lifeTimeW","E");
	dataTree->Draw("Lb_P:Lb_PT>>hData2D",dataCut+baseCut+"Lb_MassConsJpsiAndLambda_M > 5605 && Lb_MassConsJpsiAndLambda_M < 5635","E");
	dataTree->Draw("Lb_P:Lb_PT>>hSide2D",dataCut+baseCut+"Lb_MassConsJpsiAndLambda_M > 5800 && Lb_MassConsJpsiAndLambda_M < 7000","E");

	hData2D->Add(hSide2D,-factor);
	hData2D->Scale(1./hData2D->Integral());
	hMC2D->Scale(1./hMC2D->Integral());
	TH2F * ratio2D = (TH2F *)hData2D->Clone("ratio2D");
	ratio2D->Divide(hMC2D);
	ratio2D->Draw("colz text");
	c->Print("ratio_Lb_p_pt.pdf");

	ratio2D->Write();

	double Lbpt[] = {0.,3e3,4e3,5e3,6e3,7e3,8e3,10e3,100e3};
	double Lpt[] = {0.,1.8e3,2.1e3,2.4e3,2.8e3,3.2e3,100e3};
	
	TH2F * hMC2D_2 = new TH2F("hMC2D_2","hMC2D",7,Lbpt,5,Lpt);
	TH2F * hData2D_2 = new TH2F("hData2D_2","hData2D",7,Lbpt,5,Lpt);
	TH2F * hSide2D_2 = new TH2F("hSide2D_2","hSide2D",7,Lbpt,5,Lpt);
	//MCtree->Draw("Lambda0_TRUEPT:Lb_TRUEPT>>hMC2D_2","MCnorm*lifeTimeW*physRate_polp006*("+(TString)baseCut+")","E");
	MCtree->Draw("Lambda0_TRUEPT:Lb_TRUEPT>>hMC2D_2","MCnorm*lifeTimeW*physRate_polp006","E");
	dataTree->Draw("Lambda0_PT:Lb_PT>>hData2D_2",dataCut+baseCut+"Lb_MassConsJpsiAndLambda_M > 5605 && Lb_MassConsJpsiAndLambda_M < 5635","E");
	dataTree->Draw("Lambda0_PT:Lb_PT>>hSide2D_2",dataCut+baseCut+"Lb_MassConsJpsiAndLambda_M > 5800 && Lb_MassConsJpsiAndLambda_M < 7000","E");

	hData2D_2->Add(hSide2D_2,-factor);
	hData2D_2->Scale(1./hData2D_2->Integral());
	hMC2D_2->Scale(1./hMC2D_2->Integral());
	TH2F * ratio2D_2 = (TH2F *)hData2D_2->Clone("ratio2D_pt");
	ratio2D_2->Divide(hMC2D_2);
	ratio2D_2->Draw("colz text");
	c->Print("ratio_Lb_Lambda0_pt.pdf");

	ratio2D_2->Write();

	outFile->Close();

	//--------------------------------------------------------------
	
	dataTree->Draw("Lb_MM>>hbefore(100,5000,7000)",CutsDef::avoidJpsiCut_large);
	TH1F * hb = (TH1F*)gPad->GetPrimitive("hbefore");
	hb->SetTitle(0);
	hb->GetXaxis()->SetTitle("M(#Lambda#mu#mu) (MeV)");
	c->Print("Lb_MM_beforeMVAcut.pdf");
	dataTree->Draw("Lb_MM>>hafter(100,5000,7000)",CutsDef::MVAcut+CutsDef::avoidJpsiCut_large);
	TH1F * ha = (TH1F*)gPad->GetPrimitive("hafter");
	ha->SetTitle(0);
	ha->GetXaxis()->SetTitle("M(#Lambda#mu#mu) (MeV)");
	c->Print("Lb_MM_afterMVAcut.pdf");



	TFile * SelMCFile = TFile::Open("/afs/cern.ch/work/p/pluca/weighted/Lmumu/candLb_MC.root");
	TTree * BuKstmumuTree  = (TTree *)SelMCFile->Get("candBuKstmumu");
	TTree * JpsiTailTree  = (TTree *)SelMCFile->Get("candJpsiTail");
	
	BuKstmumuTree->Draw("Lb_MassConsLambda>>hKst(40,5400,6000");
	TH1F * hKst = (TH1F*)gPad->GetPrimitive("hKst");
	hKst->SetTitle(0);
	hKst->GetXaxis()->SetTitle("M(#Lambda#mu#mu) (MeV)");
	c->Print("Kst_plus_distrib.pdf");

	JpsiTailTree->Draw("Lb_MassConsLambda>>hTail(100,5400,6000");
	TH1F * hTail = (TH1F*)gPad->GetPrimitive("hTail");
	hTail->SetTitle(0);
	hTail->GetXaxis()->SetTitle("M(#Lambda#mu#mu) (MeV)");
	c->Print("Jpsi_tail_distrib.pdf");


	TFile * SelDataFile = TFile::Open("/afs/cern.ch/work/p/pluca/weighted/Lmumu/candLb.root");
	TTree * Lb2LmumuTree  = (TTree *)SelDataFile->Get("candLb2Lmumu");

	Lb2LmumuTree->Draw("Lambda0_MM>>hL0");
	TH1F * hL0 = (TH1F*)gPad->GetPrimitive("hL0");
	hL0->SetTitle(0);
	hL0->GetXaxis()->SetTitle("M(p#pi) (MeV)");
	c->Print("Lambda0_mass.pdf");

	Lb2LmumuTree->Draw("J_psi_1S_MM>>hq2");
	TH1F * hq2 = (TH1F*)gPad->GetPrimitive("hq2");
	hq2->SetTitle(0);
	hq2->GetXaxis()->SetTitle("M(#mu#mu) (MeV)");
	c->Print("mumu_mass.pdf");

	
//-----------------------------------------------------------------------------------------------


	TCut binCut = "TMath::Power(J_psi_1S_MM/1000,2) > 9.1 && TMath::Power(J_psi_1S_MM/1000,2) < 10.1";
	//MCTree->Draw("muplus_P>>hmumu_Pmu(40,0,1e4)",binCut,"E");
	//TH1F * hmumu_Pmu = (TH1F*)gPad->GetPrimitive("hmumu_Pmu",binCut,"E");
	//MCTree->Draw("pplus_P>>hmumu_Pp(40,0,1.5e4)",binCut,"E");
	//TH1F * hmumu_Pp = (TH1F*)gPad->GetPrimitive("hmumu_Pp",binCut,"E");
	//MCTree->Draw("piminus_P>>hmumu_Ppi(40,0,1.5e4)");
	//TH1F * hmumu_Ppi = (TH1F*)gPad->GetPrimitive("hmumu_Ppi",binCut,"E");
	MCtreeGeom->Draw("muplus_TRUEPT>>hmumu_PTmu(40,0,1e4)",binCut,"E");
	TH1F * hmumu_PTmu = (TH1F*)gPad->GetPrimitive("hmumu_PTmu");
	MCtreeGeom->Draw("pplus_TRUEPT>>hmumu_PTp(40,0,1.5e4)",binCut,"E");
	TH1F * hmumu_PTp = (TH1F*)gPad->GetPrimitive("hmumu_PTp");
	MCtreeGeom->Draw("piminus_TRUEPT>>hmumu_PTpi(40,0,1.5e4)",binCut,"E");
	TH1F * hmumu_PTpi = (TH1F*)gPad->GetPrimitive("hmumu_PTpi");


	//MCTreeJpsi->Draw("muplus_P>>hjpsi_Pmu(40,0,1e4)");
	//TH1F * hjpsi_Pmu = (TH1F*)gPad->GetPrimitive("hmumu_Pmu",binCut,"E");
	//MCTreeJpsi->Draw("pplus_P>>hjpsi_Pp(40,0,1.5e4)");
	//TH1F * hjpsi_Pp = (TH1F*)gPad->GetPrimitive("hmumu_Pp",binCut,"E");
	//MCTreeJpsi->Draw("piminus_P>>hjpsi_Ppi(40,0,1.5e4)");
	//TH1F * hjpsi_Ppi = (TH1F*)gPad->GetPrimitive("hmumu_Ppi",binCut,"E");
	MCtreeJpsiGeom->Draw("muplus_TRUEPT>>hjpsi_PTmu(40,0,1e4)",binCut,"E");
	TH1F * hjpsi_PTmu = (TH1F*)gPad->GetPrimitive("hjpsi_PTmu");
	MCtreeJpsiGeom->Draw("pplus_TRUEPT>>hjpsi_PTp(40,0,1.5e4)",binCut,"E");
	TH1F * hjpsi_PTp = (TH1F*)gPad->GetPrimitive("hjpsi_PTp");
	MCtreeJpsiGeom->Draw("piminus_TRUEPT>>hjpsi_PTpi(40,0,1.5e4)",binCut,"E");
	TH1F * hjpsi_PTpi = (TH1F*)gPad->GetPrimitive("hjpsi_PTpi");

	hmumu_PTmu->Scale(1./hmumu_PTmu->Integral());
	hjpsi_PTmu->Scale(1./hjpsi_PTmu->Integral());
	hmumu_PTpi->Scale(1./hmumu_PTpi->Integral());
	hjpsi_PTpi->Scale(1./hjpsi_PTpi->Integral());
	hmumu_PTp->Scale(1./hmumu_PTp->Integral());
	hjpsi_PTp->Scale(1./hjpsi_PTp->Integral());

	TH1F * ratio_mu = (TH1F*)hmumu_PTmu->Clone("ratio_mu");
	ratio_mu->Divide(hjpsi_PTmu);
	ratio_mu->GetXaxis()->SetTitle("PT muon");
	ratio_mu->GetYaxis()->SetTitle("mumu / jpsi");
	ratio_mu->Draw();
	c->Print("ratio_mu_PT.pdf");
	TH1F * ratio_pi = (TH1F*)hmumu_PTpi->Clone("ratio_pi");
	ratio_pi->Divide(hjpsi_PTpi);
	ratio_pi->GetXaxis()->SetTitle("PT pion");
	ratio_pi->GetYaxis()->SetTitle("mumu / jpsi");
	ratio_pi->Draw();
	c->Print("ratio_pi_PT.pdf");
	TH1F * ratio_p = (TH1F*)hmumu_PTp->Clone("ratio_p");
	ratio_p->Divide(hjpsi_PTp);
	ratio_p->GetXaxis()->SetTitle("PT proton");
	ratio_p->GetYaxis()->SetTitle("mumu / jpsi");
	ratio_p->Draw();
	c->Print("ratio_p_PT.pdf");

	return 0;
}


/*
	   string candfilename = "/afs/cern.ch/work/p/pluca/weighted/Lmumu/candLb_MC.root";
	   string datafilename = "/afs/cern.ch/work/p/pluca/weighted/Lmumu/candLb.root";

	   RooRealVar * vLbJpsiCons = new RooRealVar("Lb_MassConsJpsiLambda","Lb_MassConsJpsiLambda",5621.5,5350.,7000.);
	   vLbJpsiCons->setRange("Signal",5600,5625);
	   vLbJpsiCons->setRange("Sideband",5700,7000);
	   Analysis * anaLbJpsi = new Analysis("Lb2JpsiL_"+type+"_MC","Lb","candLb2JpsiL",candfilename,vLbJpsiCons);
	   Analysis * anaLbJpsi_data = new Analysis("Lb2JpsiL_"+type+"_data","Lb","candLb2JpsiL",datafilename,vLbJpsiCons);

	   MCFile = TFile::Open((TString)candfilename);
	   TTree * BdJpsiKSTree  = (TTree *)MCFile->Get("candBdJpsiKS");

	   TFile * outFile = new TFile("out.root","recreate");

	   string model = "DCB_Sn";
	   string optionsjpsi = "-quiet-range-ANDpulls-lin-stdAxis-XM(#Lambda#mu#mu) (MeV/c^{2})-noCost-nochi2-min2";

	   anaLbJpsi->SetSignal((model+"-s[7]-s2[15]").c_str());
	   anaLbJpsi->Initialize("");
	   anaLbJpsi->Fit(5400.,5750.,200,true,"-quiet-stdAxis-XM(#Lambda#mu#mu) (MeV/c^{2})-noCost-nochi2",baseCut);
	   Str2VarMap MCpars = anaLbJpsi->GetSigParams();

	   Analysis * KS = new Analysis("KS_bkg","Lb",BdJpsiKSTree,vLbJpsiCons,"DCB_OST","","-namepar");
	   KS->Fit(5350.,6000.,100,true,"-quiet-noParams-stdAxis-XM(#Lambda#mu#mu) (MeV/c^{2})");
	   Str2VarMap pars = KS->GetSigParams();
	   RooRealVar * m_shift = new RooRealVar("m_shift","shift",2.,-10.,10.);
	   if(type=="LL") m_shift = new RooRealVar("m_shift","shift",0.);
	   setConstant(&pars);
	   ModifyPars(&pars,"m",m_shift,"-shift");

	   string jpsimodel = model+"-Xn"+Form("[%f]",MCpars["n"]->getVal());
	   jpsimodel += (string)"-s[7,1,12]";//+Form("[%f]",MCpars["s"]->getVal());
	   jpsimodel += (string)"-s2[15,8,30]";//+Form("[%f,]",MCpars["s2"]->getVal());
	   jpsimodel += (string)"-a"+Form("[%f]",MCpars["a"]->getVal());
	   jpsimodel += (string)"-a2"+Form("[%f]",MCpars["a2"]->getVal());
	   jpsimodel += (string)"-f"+Form("[%f]",MCpars["f"]->getVal());

	   anaLbJpsi_data->SetSignal(jpsimodel.c_str());
	   anaLbJpsi_data->SetBkgMode( true );
	   anaLbJpsi_data->addBkgComponent("JpsiKS","DCB_OST",2.e3,"",pars);
	   anaLbJpsi_data->addBkgComponent("Comb","Exp-b[-0.003]",2.e3);
	   anaLbJpsi_data->Initialize("");
	   anaLbJpsi_data->Fit(5350.,6000.,150,true,optionsjpsi,baseCut);

	//double sigBkg = anaLbJpsi_data->GetTotBkg()->createIntegral(*vLbJpsiCons,Range("Signal"))->getVal();
	//double sideBkg = anaLbJpsi_data->GetTotBkg()->createIntegral(*vLbJpsiCons,Range("Sideband"))->getVal();
	//double nbkg = anaLbJpsi_data->GetTotNBkg()->getVal();
	//double sig = anaLbJpsi_data->GetSig()->createIntegral(*vLbJpsiCons,NormSet(*vLbJpsiCons),Range("Signal"))->getVal();
	//double nsig = anaLbJpsi_data->GetNSigPtr()->getVal();
	//double factor = 1. - nsig*sig/(nsig*sig+nbkg*sigBkg);
	cout << sigBkg << "   " << sideBkg << endl; 
	double factor = sigBkg/sideBkg;
	cout << fixed << setprecision(3) << "FACTOR = " << factor << endl;
	*/


		
