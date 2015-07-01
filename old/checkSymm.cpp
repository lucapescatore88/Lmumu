#include "ReadTree_comp.hpp"
#include "multi_analyser.hpp"
#include "Lb_cuts.hpp"
#include "RooAbsReal.h"
#include "RooSimultaneous.h"
#include "RooGenericPdf.h"
#include "TGraphErrors.h"
#include "TROOT.h"
#include "FeldmanCousins.hpp"
#include "TLegend.h"

using namespace RooFit;
using namespace std;

int main(int argc, char **argv)
{
	TString effbase = "/afs/cern.ch/user/p/pluca/work/Lb/Lmumu/results/";
	TCanvas * cDD = new TCanvas();
	TCanvas * cLL = new TCanvas();
	TCanvas * cDDB = new TCanvas();
	TCanvas * cLLB = new TCanvas();

	TString cut = "";
	if(argc > 1)
	{
		for(int a = 1; a < argc; a++)
		{
			string arg = argv[a];
			string str = arg.substr(2,arg.length()-2);
			
			if(arg == "-highestq2") cut = "_"+arg.substr(1,arg.length()-1);
			if(arg == "-highq2") cut = "_"+arg.substr(1,arg.length()-1);
			if(arg == "-lowq2") cut = "_"+arg.substr(1,arg.length()-1);
			if(arg == "-lowestq2") cut = "_"+arg.substr(1,arg.length()-1);
		}
	}
	

	double step = 0.2;
	gStyle->SetOptStat(0);
	int p = 0;

	TLegend * l1 = new TLegend(0.55,0.65,0.99,0.99);
	TLegend * l2 = new TLegend(0.55,0.65,0.99,0.99);

	for(int i = 0; i < 10; i+=3)
	{
		cout << (-1+(i+0.5)*step) << endl;
		TFile * effFile = NULL;
		TH1F * effDD = NULL, * effLL = NULL, * effLLB = NULL, * effDDB = NULL;
		effFile = TFile::Open(effbase+"Lbeff_cosThetaL_vs_cosThetaB_DD"+cut+".root");
		TH2F * effDD2D  = (TH2F *)effFile->Get("hupper_eff");
		effDD = (TH1F*)GetSliceX(effDD2D,(double)(-1+(i+0.5)*step));
		effFile = TFile::Open(effbase+"Lbeff_cosThetaL_vs_cosThetaB_LL"+cut+".root");
		TH2F * effLL2D  = (TH2F *)effFile->Get("hupper_eff");
		effLL = (TH1F*)GetSliceX(effLL2D,(double)(-1+(i+0.5)*step));
		effFile = TFile::Open(effbase+"Lbeff_cosThetaB_vs_cosThetaL_DD"+cut+".root");
		TH2F * effDDB2D  = (TH2F *)effFile->Get("hupper_eff");
		effDDB = (TH1F*)GetSliceX(effDDB2D,(double)(-1+(i+0.5)*step));
		effFile = TFile::Open(effbase+"Lbeff_cosThetaB_vs_cosThetaL_LL"+cut+".root");
		TH2F * effLLB2D  = (TH2F *)effFile->Get("hupper_eff");
		effLLB = (TH1F*)GetSliceX(effLLB2D,(double)(-1+(i+0.5)*step));
		effLL->Scale(1./effLL->Integral());
		effDD->Scale(1./effDD->Integral());
		effLLB->Scale(1./effLLB->Integral());
		effDDB->Scale(1./effDDB->Integral());
		effDDB->SetMinimum(0.);
		effLLB->SetMinimum(0.);
		effLL->SetMinimum(0.);
		effDD->SetMinimum(0.);
		//effDDB->SetMaximum(0.2);
		//effLLB->SetMaximum(0.2);
		//effLL->SetMaximum(0.2);
		//effDD->SetMaximum(0.2);
		effDD->SetTitle("DD");
		effDDB->SetTitle("DD");
		effLL->SetTitle("LL");
		effLLB->SetTitle("LL");
		
		l1->AddEntry(effDD,Form("%4.2f < cosThetaB < %4.2f",-1+i*step,-1+(i+1)*step),"P");
		l2->AddEntry(effDDB,Form("%4.2f < cosThetaL < %4.2f",-1+i*step,-1+(i+1)*step),"P");

		cDD->cd();
		effDD->SetMarkerStyle(20+p);
		effDD->SetMarkerSize(1.3);
		effDD->SetMarkerColor(1+p);
		if(i==0) effDD->Draw(); 
		else effDD->Draw("same");
		cLL->cd();
		effLL->SetMarkerStyle(20+p);
		effLL->SetMarkerSize(1.3);
		effLL->SetMarkerColor(1+p);
		if(i==0) effLL->Draw(); 
		else effLL->Draw("same");
		cDDB->cd();
		effDDB->SetMarkerStyle(20+p);
		effDDB->SetMarkerSize(1.3);
		effDDB->SetMarkerColor(1+p);
		if(i==0) effDDB->Draw(); 
		else effDDB->Draw("same");
		cLLB->cd();
		effLLB->SetMarkerStyle(20.+p);
		effLLB->SetMarkerSize(1.3);
		effLLB->SetMarkerColor(1+p);
		if(i==0) effLLB->Draw(); 
		else effLLB->Draw("same");

		effDD->SetMaximum(0.3);
		effLL->SetMaximum(0.3);
		effLLB->SetMaximum(0.3);
		effDDB->SetMaximum(0.3);


		p++;
	}

	cDD->cd();
	l1->Draw("P");
	cLL->cd();
	l1->Draw("P");
	cDDB->cd();
	l2->Draw("P");
	cLLB->cd();
	l2->Draw("P");

    cDD->Print("DDeff"+cut+".pdf");
	cLL->Print("LLeff"+cut+".pdf");
	cDDB->Print("DDBeff"+cut+".pdf");
	cLLB->Print("LLBeff"+cut+".pdf");

/*
	TFile * effFile = TFile::Open(effbase+"Lbeff2D_cosThetaB_vs_cosThetaL_DD.root");
	TH2F * eff2D  = (TH2F *)effFile->Get("hupper_eff");
	effFile = TFile::Open(effbase+"LbeffvscosThetaB_DD.root");
	TH2F * effB  = (TH2F *)effFile->Get("huppereff");
	effFile = TFile::Open(effbase+"LbeffvscosThetaL_DD.root");
	TH2F * effL  = (TH2F *)effFile->Get("huppereff");
	TH2F * res = (TH2F*)eff2D->Clone();
	res->SetName("ratio");

	for(int x = 0; x < eff2D->GetNbinsX(); x++)
		for(int y = 0; y < eff2D->GetNbinsY(); y++)
		{
			double cont = eff2D->GetBinContent(x,y);
			double vx = eff2D->GetXaxis()->GetBinCenter(x);
			double vy = eff2D->GetYaxis()->GetBinCenter(y);
			
			double binL = effL->FindBin(vx);
			double binB = effB->FindBin(vy);
			double contL = effL->GetBinContent(binL);
			double contB = effB->GetBinContent(binB);
			
			res->SetBinContent(x,y,contL*contB/cont);
		}

	TCanvas * c = new TCanvas();
	res->Draw("colz");
	c->Print("ratio.pdf");
*/
}

