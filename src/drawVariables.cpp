#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TF1.h>
#include <TH2D.h>
#include <TH2Poly.h>
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
#include "functions.hpp"

using namespace std;



int main(int argc, char **argv)
{
	TString type = "all";	
	bool dovars = false;
	if(argc > 1)
	{
		for(int a = 1; a < argc; a++)
		{
			string arg = argv[a];
			string str = arg.substr(2,arg.length()-2);
			if(arg.find("-t") != string::npos) type = str;
			if(arg.find("-v") != string::npos) dovars = true;
		}
	}

	cout << "Analysing " << type << endl;

	TCut mvaCut = CutsDef::MVAcut;
	TCut baseCut = mvaCut + CutsDef::TrigPassed + CutsDef::Lmasscut;
	TCut dataCut = "TMath::Abs(J_psi_1S_MM - 3096) < 30";
	TCut MCcut = "TMath::Abs(J_psi_1S_MM - 3096) < 30";
	if(type == "DD") baseCut += CutsDef::DDcut;
	else if(type == "LL") baseCut += CutsDef::LLcut;



	//Attach data files	

	TCanvas *c = new TCanvas("c","");

	TString kinw = "Lb_weight";

	TString MCFileName = "/afs/cern.ch/work/p/pluca/Lmumu/weighted/Lb2Lmumu_MC_Pythia8_NBweighted.root";
	TString MCJpsiFileName = "/afs/cern.ch/work/p/pluca/Lmumu/weighted/Lb2JpsiL_MC_Pythia8_NBweighted.root";
	TString DataFileName = "/afs/cern.ch/work/p/pluca/Lmumu/weighted/Lb2Lmumu_CL_NBweighted.root";
	TString MCFileGeomName = "/afs/cern.ch/work/p/pluca/Lmumu/weighted/Lb2Lmumu_geomMC_Pythia8_NBweighted.root";
	TString MCFileJpsiGeomName = "/afs/cern.ch/work/p/pluca/Lmumu/weighted/Lb2JpsiL_geomMC_Pythia8_NBweighted.root";
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

	TFile * outFile = new TFile("MC_reweight"+type+".root","recreate");


	gStyle->SetOptStat(0);

	if(dovars)
	{

		int nvars = 28;
		string vars[] = {"Lambda0_ENDVERTEX_Z","muplus_PIDmu","muplus_ProbNNmu","muplus_IPCHI2_OWNPV","J_psi_1S_IPCHI2_OWNPV","Lb_TAU",
			"Lb_P","Lb_PT","Lambda0_PT","Lambda0_IPCHI2_OWNPV","Lambda0_FD_OWNPV","Lb_DIRA_OWNPV","muplus_TRACK_CHI2NDOF",
			"TMath::Min(muplus_IPCHI2_OWNPV,muminus_IPCHI2_OWNPV)","TMath::Max(muplus_IPCHI2_OWNPV,muminus_IPCHI2_OWNPV)",
			"Lb_MassConsLambda_chi2[0]/Lb_MassConsLambda_nDOF[0]","pplus_IPCHI2_OWNPV","piminus_IPCHI2_OWNPV","weight",
			"piminus_P","piminus_PT","cosThetaB","cosThetaL","Lambda0_P","pplus_P","pplus_PT","muplus_P","muplus_PT"};
		double min[] = {-200, -4, 0, 0,    0,    0,    0,     0,   0,   0,   0,    0.9999, 0, 0,   0,   0,   0,    0,    0, 1500, 0,  0, 0, 7e3, 0, 0,0,0};
		double max[] = {4e3,  14, 1, 50e3, 20e3, 0.05, 400e3, 5e4, 10e3, 1e3, 4000, 1.,     3, 5e4, 5e4, 120, 3e3, 10e3, 1, 15e3, 5e3, 1, 1, 12e4, 5e5, 1e4,6e5,12e4};
		bool log[] = {false,false,false,true,true,true,false,false,false,true,false,true,false,true,true,true,true,true,true,false,false,false,false,true,false,false,true,true};

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
			TString nameMC_MuMufullrw_h = "hMCMuMufullrw"+namevar+Form("_%i",v);
			TString nameDatah = "hData"+namevar+Form("_%i",v);
			TString nameSideh = "hSide"+namevar+Form("_%i",v);
			TH1F * hMC = NULL, * hMCJpsi = NULL, * hMC_q2rw = NULL, * hMC_Prw = NULL, * hMC_fullrw = NULL, * hMC_MuMufullrw = NULL;

			if((MCcut+baseCut)!="")
			{
				MCtree->Draw(vars[v]+">>"+nameMCh+binning,"MCnorm*lifeTimeW*("+(TString)(MCcut+baseCut)+")","E");
				hMC = (TH1F*)gPad->GetPrimitive(nameMCh);
				MCJpsitree->Draw(vars[v]+">>"+nameMC_q2rw_h+binning,"physRate_pol0*MCnorm*lifeTimeW*("+(TString)(MCcut+baseCut)+")","E");
				hMC_q2rw = (TH1F*)gPad->GetPrimitive(nameMC_q2rw_h);
				MCJpsitree->Draw(vars[v]+">>"+nameMC_Prw_h+binning,kinw+"*MCnorm*lifeTimeW*("+(TString)(MCcut+baseCut)+")","E");
				hMC_Prw = (TH1F*)gPad->GetPrimitive(nameMC_Prw_h);
				MCJpsitree->Draw(vars[v]+">>"+nameMC_fullrw_h+binning,kinw+"*physRate_pol0*MCnorm*lifeTimeW*("+(TString)(MCcut+baseCut)+")","E");
				hMC_fullrw = (TH1F*)gPad->GetPrimitive(nameMC_fullrw_h);
				MCtree->Draw(vars[v]+">>"+nameMC_MuMufullrw_h+binning,kinw+"*physRate_polp006*MCnorm*lifeTimeW*("+(TString)(MCcut+baseCut)+")","E");
				hMC_MuMufullrw = (TH1F*)gPad->GetPrimitive(nameMC_MuMufullrw_h);
				MCJpsitree->Draw(vars[v]+">>"+nameMCJpsih+binning,"MCnorm*lifeTimeW*("+(TString)(MCcut+baseCut)+")","E");
				hMCJpsi = (TH1F*)gPad->GetPrimitive(nameMCJpsih);
			}
			else
			{
				MCtree->Draw(vars[v]+">>"+nameMCh+binning,"MCnorm*lifeTimeW","E");
				hMC = (TH1F*)gPad->GetPrimitive(nameMCh);
				MCJpsitree->Draw(vars[v]+">>"+nameMC_q2rw_h+binning,"physRate_pol0*MCnorm*lifeTimeW","E");
				hMC_q2rw = (TH1F*)gPad->GetPrimitive(nameMC_q2rw_h);
				MCJpsitree->Draw(vars[v]+">>"+nameMC_Prw_h+binning,kinw+"*MCnorm*lifeTimeW","E");
				hMC_Prw = (TH1F*)gPad->GetPrimitive(nameMC_Prw_h);
				MCJpsitree->Draw(vars[v]+">>"+nameMC_fullrw_h+binning,kinw+"*physRate_pol0*MCnorm*lifeTimeW","E");
				hMC_fullrw = (TH1F*)gPad->GetPrimitive(nameMC_fullrw_h);
				MCtree->Draw(vars[v]+">>"+nameMC_MuMufullrw_h+binning,kinw+"*physRate_polp006*MCnorm*lifeTimeW","E");
				hMC_MuMufullrw = (TH1F*)gPad->GetPrimitive(nameMC_MuMufullrw_h);
				MCJpsitree->Draw(vars[v]+">>"+nameMCJpsih+binning,"MCnorm*lifeTimeW","E");
				hMCJpsi = (TH1F*)gPad->GetPrimitive(nameMCJpsih);
			}	

			dataTree->Draw(vars[v]+">>"+nameDatah+binning,dataCut+baseCut+"Lb_MassConsJpsiAndLambda_M > 5605 && Lb_MassConsJpsiAndLambda_M < 5635","E");
			TH1F * hData = (TH1F*)gPad->GetPrimitive(nameDatah);
			dataTree->Draw(vars[v]+">>"+nameSideh+binning,dataCut+baseCut+"Lb_MassConsJpsiAndLambda_M > 5800 && Lb_MassConsJpsiAndLambda_M < 7000","E");
			TH1F * hSide = (TH1F*)gPad->GetPrimitive(nameSideh);
			hData->Add(hSide,-factor);

			hMC->SetMarkerStyle(24);
			hMC->SetMarkerColor(1);
			hMC->SetMarkerSize(0.8);
			hMC_MuMufullrw->SetMarkerStyle(24);
			hMC_MuMufullrw->SetMarkerColor(8);
			hMC_MuMufullrw->SetMarkerSize(0.7);
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
			hMCJpsi->SetMarkerColor(1);
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
			hMC_MuMufullrw->Scale(1./hMC_MuMufullrw->Integral());
			hData->Scale(1./hData->Integral());
			hSide->Scale(1./hSide->Integral());



			if(v==0)
			{
				leg->AddEntry(hData,"Data (JpsiL)","P");
				leg->AddEntry(hSide,"Side band","P");
				//leg->AddEntry(hMC,"MC Lmumu","P");
				//leg->AddEntry(hMC_q2rw,"MC JpsiL dec mod W","P");
				//leg->AddEntry(hMC_Prw,"MC JpsiL Lb mom W","P");
				leg->AddEntry(hMC_fullrw,"MC Jpsi fully W","P");
				//leg->AddEntry(hMC_MuMufullrw,"MC Lmumu fully W","P");
				leg->AddEntry(hMCJpsi,"MC JpsiL","P");
			}

			hData->SetTitle(0);
			hData->GetXaxis()->SetTitle(vars[v].c_str());
			hData->Draw();
			//hMC->Draw("same");
			hSide->Draw("same");
			hMC_fullrw->Draw("same");
			//hMC_MuMufullrw->Draw("same");
			//hMC_q2rw->Draw("same");
			//hMC_Prw->Draw("same");
			hMCJpsi->Draw("same");
			leg->Draw();
			c->Print("MC_data_comp/"+namevar+"_plot"+type+".pdf");
		}
	}

	//--------------------------------------------------------------------------------


	if(type=="DD" || type == "LL")
	{
		TString binning = "(10,-1,1)";
		MCJpsitree->Draw("cosThetaB>>hMCLL"+binning,"("+(TString)(MCcut+baseCut)+")*"+kinw+"*physRate_pol0*MCnorm*lifeTimeW","E");
		TH1F * hMCLL = (TH1F*)gPad->GetPrimitive("hMCLL");
		dataTree->Draw("cosThetaB>>hDataLL"+binning,dataCut+baseCut+"Lb_MassConsJpsiAndLambda_M > 5605 && Lb_MassConsJpsiAndLambda_M < 5635","E");
		TH1F * hDataLL = (TH1F*)gPad->GetPrimitive("hDataLL");
		dataTree->Draw("cosThetaB>>hSideLL"+binning,dataCut+baseCut+"Lb_MassConsJpsiAndLambda_M > 5800 && Lb_MassConsJpsiAndLambda_M < 7000","E");
		TH1F * hSideLL = (TH1F*)gPad->GetPrimitive("hSideLL");

		hDataLL->SetMarkerSize(0.8);
		hDataLL->SetMarkerColor(1);
		if(type=="DD") hDataLL->SetMarkerStyle(21);
		else hDataLL->SetMarkerStyle(20);
		hMCLL->SetMarkerSize(0.8);
		hMCLL->SetMarkerColor(2);
		if(type=="DD") hMCLL->SetMarkerStyle(21);
		else hMCLL->SetMarkerStyle(20);

		hDataLL->Add(hSideLL,-factor);
		hMCLL->Scale(1./hMCLL->Integral());
		hDataLL->Scale(1./hDataLL->Integral());
		hDataLL->Write("cosThetaB_"+type+"_data");
		hMCLL->Write("cosThetaB_"+type+"_MC");
	}

	//-----------------------------------------------------------------------------

	c->SetLogy(0);

	double p[] = {0,5e4,6e4,7e4,8e4,10e4,100e4};
	double pt[] = {0.,3e3,4e3,5e3,6e3,7e3,8e3,10e3,100e3};
	/*	
		TH2Poly * hMC2D = new TH2Poly("hMC2D","",pt[0],pt[8],p[0],p[8]);
		TH2Poly * hData2D = new TH2Poly("hData2D","",pt[0],pt[8],p[0],p[8]);
		TH2Poly * hSide2D = new TH2Poly("hSide2D","",pt[0],pt[8],p[0],p[8]);
		for(int i = 0; i < 7; i++)
		for(int j = 0; j < 8; j++)
		{
		hMC2D->AddBin(pt[i],p[j],pt[i+1],p[j+1]);
		hData2D->AddBin(pt[i],p[j],pt[i+1],p[j+1]);
		hSide2D->AddBin(pt[i],p[j],pt[i+1],p[j+1]);
		}

		hMC2D->AddBin(pt[7],p[0],pt[8],p[1]);
		hData2D->AddBin(pt[7],p[0],pt[8],p[1]);
		hSide2D->AddBin(pt[7],p[0],pt[8],p[1]);
		hMC2D->AddBin(pt[7],p[1],pt[8],p[2]);
		hData2D->AddBin(pt[7],p[1],pt[8],p[2]);
		hSide2D->AddBin(pt[7],p[1],pt[8],p[2]);
		hMC2D->AddBin(pt[7],p[2],pt[8],p[3]);
		hData2D->AddBin(pt[7],p[2],pt[8],p[3]);
		hSide2D->AddBin(pt[7],p[2],pt[8],p[3]);
		hMC2D->AddBin(pt[7],p[3],pt[8],p[4]);
		hData2D->AddBin(pt[7],p[3],pt[8],p[4]);
		hSide2D->AddBin(pt[7],p[3],pt[8],p[4]);
		hMC2D->AddBin(pt[7],p[4],pt[8],p[5]);
		hData2D->AddBin(pt[7],p[4],pt[8],p[5]);
		hSide2D->AddBin(pt[7],p[4],pt[8],p[5]);
		hMC2D->AddBin(pt[7],p[5],pt[8],p[8]);
		hData2D->AddBin(pt[7],p[5],pt[8],p[8]);
		hSide2D->AddBin(pt[7],p[5],pt[8],p[8]);
		*/


	TH2F * hMC2D = new TH2F("hMC2D","",7,pt,5,p);
	TH2F * hData2D = new TH2F("hData2D","",7,pt,5,p);
	TH2F * hSide2D = new TH2F("hSide2D","",7,pt,5,p);
	TString var = "TMath::Sqrt(TMath::Power(Lb_TRUEP_X,2)+TMath::Power(Lb_TRUEP_Y,2)+TMath::Power(Lb_TRUEP_Z,2)):Lb_TRUEPT";
	MCJpsitree->Draw(var+">>hMC2D","MCnorm*lifeTimeW*physRate_pol0","E");
	dataTree->Draw("Lb_P:Lb_PT>>hData2D",dataCut+baseCut+"Lb_MassConsJpsiAndLambda_M > 5605 && Lb_MassConsJpsiAndLambda_M < 5635","E");
	dataTree->Draw("Lb_P:Lb_PT>>hSide2D",dataCut+baseCut+"Lb_MassConsJpsiAndLambda_M > 5800 && Lb_MassConsJpsiAndLambda_M < 7000","E");

	hData2D->Add(hSide2D,-factor);
	hData2D->Scale(1./hData2D->Integral());
	hMC2D->Scale(1./hMC2D->Integral());
	TH2F * ratio2D = (TH2F *)hData2D->Clone("ratio2D");
	ratio2D->Divide(hMC2D);
	ratio2D->GetXaxis()->SetTitle("p_{T}(#Lambda_{b}) (MeV/c)");
	ratio2D->GetYaxis()->SetTitle("p(#Lambda_{b}) (MeV/c)");
	ratio2D->Draw("colz text");
	c->Print("ratio_Lb_p_pt.pdf");

	ratio2D->Write();

	//double Lbpt[] = {0.,3e3,4e3,5e3,6e3,7e3,8e3,10e3,100e3};
	double pt2[] = {0.,5e3,6e3,7e3,8e3,11e3,100e3};
	double Lpt[] = {0.,2.1e3,2.4e3,2.8e3,3.2e3,100e3};

	TH2F * hMC2D_2 = new TH2F("hMC2D_2","hMC2D",6,pt2,5,Lpt);
	TH2F * hData2D_2 = new TH2F("hData2D_2","hData2D",6,pt2,5,Lpt);
	TH2F * hSide2D_2 = new TH2F("hSide2D_2","hSide2D",6,pt2,5,Lpt);
	MCJpsitree->Draw("Lambda0_TRUEPT:Lb_TRUEPT>>hMC2D_2","MCnorm*lifeTimeW*physRate_pol0","E");
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

	double pidB[] = {-4,4,6,9,11.5,16};

	TH1F * hMC2D_3 = new TH1F("hMC2D_3","hMC2D",5,pidB);
	TH1F * hData2D_3 = new TH1F("hData2D_3","hData2D",5,pidB);
	TH1F * hSide2D_3 = new TH1F("hSide2D_3","hSide2D",5,pidB);
	MCJpsitree->Draw("muplus_PIDmu>>hMC2D_3","MCnorm*lifeTimeW*physRate_pol0","E");
	dataTree->Draw("muplus_PIDmu>>hData2D_3",dataCut+baseCut+"Lb_MassConsJpsiAndLambda_M > 5605 && Lb_MassConsJpsiAndLambda_M < 5635","E");
	dataTree->Draw("muplus_PIDmu>>hSide2D_3",dataCut+baseCut+"Lb_MassConsJpsiAndLambda_M > 5800 && Lb_MassConsJpsiAndLambda_M < 7000","E");

	hData2D_3->Add(hSide2D_3,-factor);
	hData2D_3->Scale(1./hData2D_3->Integral());
	hMC2D_3->Scale(1./hMC2D_3->Integral());
	TH2F * ratio2D_3 = (TH2F *)hData2D_3->Clone("ratio_PIDmu");
	ratio2D_3->Divide(hMC2D_3);
	ratio2D_3->Draw("colz text");
	c->Print("ratio_PIDmu.pdf");

	ratio2D_3->Write();

	outFile->Close();

	//--------------------------------------------------------------
	/*
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



	   TFile * SelMCFile = TFile::Open("/afs/cern.ch/work/p/pluca/Lmumu/weighted/candLb_MC.root");
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


	   TFile * SelDataFile = TFile::Open("/afs/cern.ch/work/p/pluca/Lmumu/weighted/candLb.root");
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
	   */

	//-----------------------------------------------------------------------------------------------


	TCut binCut = "TMath::Power(J_psi_1S_MM/1000,2) > 9.1 && TMath::Power(J_psi_1S_MM/1000,2) < 10.1";
	MCtreeGeom->Draw("muplus_TRUEPT>>hmumu_PTmu(40,0,1e4)",binCut,"E");
	TH1F * hmumu_PTmu = (TH1F*)gPad->GetPrimitive("hmumu_PTmu");
	MCtreeGeom->Draw("pplus_TRUEPT>>hmumu_PTp(40,0,1.5e4)",binCut,"E");
	TH1F * hmumu_PTp = (TH1F*)gPad->GetPrimitive("hmumu_PTp");
	MCtreeGeom->Draw("piminus_TRUEPT>>hmumu_PTpi(40,0,1.5e4)",binCut,"E");
	TH1F * hmumu_PTpi = (TH1F*)gPad->GetPrimitive("hmumu_PTpi");


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

	//-------------------------------------------------------------------------------------

	MCtree->Draw("TMath::Power(J_psi_1S_MM/1000,2)>>hw1(40,0,20)","physRate_pol0_wilson1_noDecay","E");
	TH1F * hw1 = (TH1F*)gPad->GetPrimitive("hw1");
	MCtree->Draw("TMath::Power(J_psi_1S_MM/1000,2)>>hw2(40,0,20)","physRate_pol0_wilson2_noDecay","E");
	TH1F * hw2 = (TH1F*)gPad->GetPrimitive("hw2");
	MCtree->Draw("TMath::Power(J_psi_1S_MM/1000,2)>>hw3(40,0,20)","physRate_pol0_wilson3_noDecay","E");
	TH1F * hw3 = (TH1F*)gPad->GetPrimitive("hw3");
	MCtree->Draw("TMath::Power(J_psi_1S_MM/1000,2)>>hdef(40,0,20)","physRate_polp006_noDecay","E");
	TH1F * hdef = (TH1F*)gPad->GetPrimitive("hdef");

	//hdef->Scale(1./hdef->Integral());
//	hw1->Scale(hdef->Integral(15,20)/hw1->Integral(15,20));
//	hw2->Scale(hdef->Integral(15,20)/hw2->Integral(15,20));
//	hw3->Scale(hdef->Integral(15,20)/hw3->Integral(15,20));

	hdef->SetMarkerSize(0.8);
	hdef->SetMarkerColor(1);
	hdef->SetMarkerStyle(21);
	hw1->SetMarkerSize(0.8);
	hw1->SetMarkerColor(2);
	hw1->SetMarkerStyle(22);
	hw2->SetMarkerSize(0.8);
	hw2->SetMarkerColor(3);
	hw2->SetMarkerStyle(22);
	hw3->SetMarkerSize(0.8);
	hw3->SetMarkerColor(4);
	hw3->SetMarkerStyle(22);

	hdef->SetTitle(0);
	hdef->Draw();
	hw1->Draw("same");
	hw2->Draw("same");
	hw3->Draw("same");

	TLegend * legw = new TLegend(0.1,0.65,0.5,0.9);
	legw->AddEntry(hdef, "SM");
	legw->AddEntry(hw1, "C_{7}^{NP} = -0.05, C_{9}^{NP} = -1.3");
	legw->AddEntry(hw2, "C_{7}^{NP} = -0.03, C_{9}^{NP} = -1.7");
	legw->AddEntry(hw3, "C_{7}^{NP} = 0.01, C_{9}^{NP} = -2.1");
	legw->Draw();
	c->Print("wilson_q2.pdf");


	//-------------------------------------------------------------------------------------



	if(type=="DD" || type == "LL")
	{
		double min = 1.1, max = 6;
		TString selDataFileName = "/afs/cern.ch/work/p/pluca/Lmumu/weighted/candLb.root";
		TFile * selDataFile = TFile::Open(selDataFileName);
		TTree * SelDataTree = (TTree *)selDataFile->Get("candLb2Lmumu");

		RooRealVar * cosThetaL = new RooRealVar("cosThetaL","cosThetaL",0.,-1.,1.);
		RooRealVar * cosThetaB = new RooRealVar("cosThetaB","cosThetaB",0.,-1.,1.);

		RooAbsPdf * eff = NULL, * effB = NULL, * dummy = NULL, * dummyB = NULL;	
		if(type=="DD") getEfficiencies(15,20,&dummy,&eff,&dummyB,&effB);
		else getEfficiencies(min,max,&eff,&dummy,&effB,&dummyB);

		SelDataTree->Draw("cosThetaB>>hcosThetaB(50,-1,1)","Lb_MassConsLambda > 5700 && "+(TString)Form("TMath::Power(J_psi_1S_MM/1000,2) >= %e && TMath::Power(J_psi_1S_MM/1000,2) < %e",min,max),"E");
		TH1F * hcosThetaB = (TH1F*)gPad->GetPrimitive("hcosThetaB");
		SelDataTree->Draw("cosThetaL>>hcosThetaL(50,-1,1)","Lb_MassConsLambda > 5700 && "+(TString)Form("TMath::Power(J_psi_1S_MM/1000,2) >= %e && TMath::Power(J_psi_1S_MM/1000,2) < %e",min,max),"E");
		TH1F * hcosThetaL = (TH1F*)gPad->GetPrimitive("hcosThetaL");

		RooDataHist * hhcosThetaL = new RooDataHist("hhcosThetaL","",*cosThetaL,hcosThetaL);
		RooDataHist * hhcosThetaB = new RooDataHist("hhcosThetaB","",*cosThetaB,hcosThetaB);

		GetFrame(cosThetaB,effB,hhcosThetaB,(string)"-noCost-t"+(string)type+" events",50,NULL,vector<string>(),"cos#theta_{#Lambda}")->Draw();
		c->Print("cosThetaB_bkg_"+type+".pdf");
		GetFrame(cosThetaL,eff,hhcosThetaL,(string)"-noCost-t"+(string)type+" events",50,NULL,vector<string>(),"cos#theta_{l}")->Draw();
		c->Print("cosThetaL_bkg_"+type+".pdf");

		//RooAbsPdf * bkgKey = NULL, * bkgBKey = NULL;
		//buildBkgPdfs(min,max,(TString)type,baseCut,&bkgKey,&bkgBKey);

		//GetFrame(cosThetaB,bkgBKey,hhcosThetaB,(string)"-noCost-t"+(string)type+" events",50,NULL,0,"cos#theta_{#Lambda}")->Draw();
		//c->Print("cosThetaB_bkg_"+type+"_Key.pdf");
		//GetFrame(cosThetaL,bkgKey,hhcosThetaL,(string)"-noCost-t"+(string)type+" events",50,NULL,0,"cos#theta_{l}")->Draw();
		//c->Print("cosThetaL_bkg_"+type+"_Key.pdf");
	}

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



