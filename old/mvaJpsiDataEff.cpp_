#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TCut.h"
#include "Lb_cuts.hpp"
#include <string>
#include "analyser.hpp"

using namespace std;

int main(int argc, char** argv)
{
	TString namefilebkg = "/afs/cern.ch/work/p/pluca/weighted/Lmumu/Lb2Lmumu_CL_NBweighted.root";
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

	TFile * fbkg = TFile::Open(namefilebkg);
	TTree * treeBkgTmp = (TTree *)fbkg->Get("tree");
	TTree * treeBkg = (TTree *)treeBkgTmp->Clone("treeBkg");

	double maxSig = 5635, minSig = 5605;//, minSide = 5700;

	TCut jpsiCut = "TMath::Abs(J_psi_1S_MM - 3096) < 20";
	TCut baseCut = CutsDef::Lmasscut + jpsiCut;
	if(type == "DD") baseCut += CutsDef::DDcut;
	else if(type == "LL") baseCut += CutsDef::LLcut;

	TCanvas * d = new TCanvas();
	treeBkg->Draw("Lb_MassConsJpsiLambda>>hsideband(72,5530,7000)",baseCut,"E");
	TH1 * hh = (TH1*)gPad->GetPrimitive("hsideband");
	TF1 * f = new TF1("f","[0]*TMath::Exp([1]*x)+[2]*TMath::Gaus(x,[3],[4])",5000,7100);
	TF1 * exp = new TF1("exp","[0]*TMath::Exp([1]*x)",5000,7100);
	f->SetParameter(0,1000);
	f->SetParameter(1,-0.005);
	f->SetParameter(2,10000);
	f->SetParameter(3,5620);
	f->SetParameter(4,7);
	hh->Fit(f,"","",5560,7000);
	double pars[5];
	f->GetParameters(pars);
	exp->SetParameter(0,pars[0]);
	exp->SetParameter(1,pars[1]);
	double bkgNum = exp->Integral(minSig,maxSig) / f->Integral(minSig,maxSig);
	cout << bkgNum << endl;

	TLine l1(minSig,0,minSig,10000);
	TLine l2(maxSig,0,maxSig,10000);
	l1.Draw("same");
	l2.Draw("same");
	d->Print("side_pre"+type+".pdf");


	treeBkg->Draw("Lb_MassConsJpsiLambda>>hsidepost(72,5560,7000)",baseCut + CutsDef::MVAcut,"E");
	TH1 * hh2 = (TH1*)gPad->GetPrimitive("hsidepost");
	TF1 * f2 = new TF1("f","[0]*TMath::Exp([1]*x)+[2]*TMath::Gaus(x,[3],[4])",5000,7100);
	TF1 * exp2 = new TF1("f","[0]*TMath::Exp([1]*x)",5000,7100);
	f2->SetParameter(0,100);
	f2->SetParameter(1,-0.001);
	f2->SetParameter(2,1000);
	f2->SetParameter(3,5620);
	f2->SetParameter(4,15);
	hh2->Fit(f2,"","",5560,7000);
	double pars2[5];
	f2->GetParameters(pars2);
	exp2->SetParameter(0,pars2[0]);
	exp2->SetParameter(1,pars2[1]);
	double bkgNum_post = exp2->Integral(minSig,maxSig) / f2->Integral(minSig,maxSig);
	cout << bkgNum_post << endl;
	l1.Draw("same");
	l2.Draw("same");
	d->Print("side_post"+type+".pdf");



	TCut LbMassCut = Form("Lb_MassConsJpsiLambda > %e && Lb_MassConsJpsiLambda < %e",minSig,maxSig);
	treeBkg->Draw("Lb_MassConsJpsiLambda>>hpre",baseCut + LbMassCut,"E");
	TH1 * hh3 = (TH1*)gPad->GetPrimitive("hpre");
	double preMVA = hh3->GetEntries();
	treeBkg->Draw("Lb_MassConsJpsiLambda>>hpost",baseCut + CutsDef::MVAcut + LbMassCut,"E");
	TH1 * hh4 = (TH1*)gPad->GetPrimitive("hpost");
 	double postMVA = hh4->GetEntries();

	cout << "nBkg = " << bkgNum << "   nBkg post = " << bkgNum_post << endl;
	cout << "ntot = " << preMVA << "   ntot post = " << postMVA << endl;
	double eff = (postMVA*(1. - bkgNum_post))/(preMVA*(1. - bkgNum));
	double ntot = preMVA*(1. - bkgNum);
	cout << "MVA eff from Jpsi data -> " << eff << " \\pm " << TMath::Sqrt( TMath::Power(0.05,2) + eff*(1-eff)/ntot ) << endl;
}
