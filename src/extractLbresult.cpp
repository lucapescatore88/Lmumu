#include <vector>
#include <sstream>
#include <iostream>
#include <string>
#include <time.h>

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TEntryList.h"
#include "TCanvas.h"
#include "TBranch.h"
#include "TCut.h"
#include "TGraphErrors.h"
#include "TROOT.h"
#include "TLegend.h"

#include "general_functions.hpp"
#include "ReadTree_comp.hpp"
#include "analyser.hpp"
#include "Lb_cuts.hpp"


using namespace std;
using namespace RooFit;


TGraphErrors * mergeGr(TString title, TString f1name, TString f2name)
{
	TFile * f1 = TFile::Open(f1name);
	TGraphErrors * gr1 = new TGraphErrors();
	f1->GetObject(title,gr1);
	TFile * f2 = TFile::Open(f2name);
	TGraphErrors * gr2 = new TGraphErrors();
	f2->GetObject(title,gr2);

	TGraphErrors * gr = new TGraphErrors();

	unsigned p = 1;
	double x,y;
	for(int i = 1; i < gr1->GetN(); i++) 
	{ 
		gr1->GetPoint(i,x,y);
		gr->SetPoint(p,x,y);
		gr->SetPointError(p,gr1->GetErrorX(i),gr1->GetErrorY(i));
		p++;
	}
	for(int i = 1; i < gr2->GetN(); i++)
	{ 
		gr2->GetPoint(i,x,y);
		gr->SetPoint(p,x,y);
		gr->SetPointError(p,gr2->GetErrorX(i),gr2->GetErrorY(i));
		p++;
	}

	return gr;
}



int main(int argc, char **argv)
{
	string print = "yield";
	if(argc > 1)
		for(int a = 1; a < argc; a++)
		{
			string arg = argv[a];
			string str = arg.substr(2,arg.length()-2);

			if(arg.find("-p") != string::npos) print = str;
		}


	gROOT->ProcessLine(".x lhcbStyle.C");

	TCanvas * cc = new TCanvas();
	TFile * histFile = new TFile("Lbresults.root","recreate");

	vector <TString> effnames;
	effnames.push_back("geom");
	effnames.push_back("det");
	effnames.push_back("reco");
	effnames.push_back("trig");
	effnames.push_back("mva");

	double Jpsi2mumuBr = 0.0593;
	double Jpsi2mumuBr_err = 0.0006;
	double jpsiyield, jpsiyield_lowSel;
	map<double,double> pdfsys, pideff;

	pdfsys[17.5] = 0.01;
	pdfsys[11.75] = 0.032;
	pdfsys[15.5] = 0.028;
	pdfsys[17.0] = 0.014;
	pdfsys[19.0] = 0.025;
	pdfsys[3.55] = 0.025;
	pdfsys[1.05] = 0.025;
	pdfsys[3.0] = 0.025;
	pdfsys[5.0] = 0.025;
	pdfsys[7.0] = 0.025;

	pideff[17.5] = 1.00281;
	pideff[11.75] = 1.00151;
	pideff[15.5] = 1.00431;
	pideff[17.0] = 1.00215;
	pideff[19.0] = 1.00226;
	pideff[3.55] = 0.99589;
	pideff[1.05] = 0.99418;
	pideff[3.0] = 0.99523;
	pideff[5.0] = 0.99699;
	pideff[7.0] = 0.99805;

	TGraphErrors * combresgr = new TGraphErrors();
	TGraphErrors * corr_sys_gr = new TGraphErrors();
	TGraphErrors * LL_uncorrplot = new TGraphErrors();
	TGraphErrors * DD_uncorrplot = new TGraphErrors();
	TGraphErrors * LL_eff = new TGraphErrors();
	TGraphErrors * DD_eff = new TGraphErrors();
	TGraphErrors * LL_stat = new TGraphErrors();
	TGraphErrors * DD_stat = new TGraphErrors();
	TGraphErrors * LL_ncsys = new TGraphErrors();
	TGraphErrors * DD_ncsys = new TGraphErrors();
	TGraphErrors * LL_norm = new TGraphErrors();
	TGraphErrors * DD_norm = new TGraphErrors();

	TString type = "LL";
	TGraphErrors * DDres = NULL, * LLres = NULL;

	for(int t = 0; t < 2; t++)
	{
		if(t==0) type = "DD";
		else  type = "LL"; 	

		cout << "Event type: " << type << endl;

		TString datafilename = "/afs/cern.ch/work/p/pluca/results/Lbyield_"+type+".root";
		TFile * fyield = TFile::Open(datafilename);
		TGraphErrors * yield, * jpsiyieldgr, * jpsiyieldgr_lowSel;
		fyield->GetObject("q2plot",yield);
		fyield->GetObject("jpsiyield_highq2",jpsiyieldgr);
		fyield->GetObject("jpsiyield_lowq2",jpsiyieldgr_lowSel);
		double dummy = 0;
		jpsiyieldgr->GetPoint(1,dummy,jpsiyield);
		double jpsierr = jpsiyieldgr->GetErrorY(1);
		jpsiyieldgr_lowSel->GetPoint(1,dummy,jpsiyield_lowSel);
		double jpsierr_lowSel = jpsiyieldgr_lowSel->GetErrorY(1);

		cout << "Njpsi = " << jpsiyield << " \\pm " << jpsierr << endl; 
		cout << "Njpsi lowSel = " << jpsiyield_lowSel << " \\pm " << jpsierr_lowSel << endl;

		TString f1name = "/afs/cern.ch/work/p/pluca/results/LbreleffAndSysvsq2_"+type+"_7bins.root";
		TString f2name = "/afs/cern.ch/work/p/pluca/results/LbreleffAndSysvsq2_"+type+"_2bins.root";

		TGraphErrors * resultgr = new TGraphErrors();
		TGraphErrors * toteffgr = mergeGr("toteff", f1name, f2name);
		TGraphErrors * toteffgr_lowSel = mergeGr("toteff_lowSel", f1name, f2name);
		vector < TGraphErrors * > sysgr;
		for(unsigned s = 0; s < effnames.size(); s++) sysgr.push_back(mergeGr(effnames[s]+"sys", f1name, f2name));

		unsigned p = 1;
		for(int i = 1; i < toteffgr->GetN(); i++)	
		{
			double corr_err = 0, nocorr_err = 0;
			double errq2 = 0, q2, tmpq2, nn, errnn = 0, toteff;

			toteffgr->GetPoint(i,tmpq2,toteff);
			double erreff = toteffgr->GetErrorY(i-1);

			if( tmpq2 < 1.e-10 || tmpq2 > 1.e9 ) continue;
			for(int h = 0; h < yield->GetN(); h++) 
			{
				yield->GetPoint(h,q2,nn);
				errq2 = yield->GetErrorX(h);
				errnn = yield->GetErrorY(h);
				if(TMath::Abs(q2 - tmpq2) < 1e-5) break;
			}
			if( TMath::Abs(q2 - tmpq2) > 1e-5 ) continue;
			if(type=="LL") toteff*=pideff[q2];

			for(unsigned s = 0; s < effnames.size(); s++)
			{
				double sys;
				sysgr[s]->GetPoint(i,tmpq2,sys);

				if(s < 2) corr_err = TMath::Sqrt(TMath::Power(corr_err,2) + TMath::Power(sys,2));
				else nocorr_err = TMath::Sqrt(TMath::Power(nocorr_err,2) + TMath::Power(sys,2));
			}

			double cur_jpsiyield = jpsiyield;
			if(q2 < 9) 
			{
				cur_jpsiyield = jpsiyield_lowSel;
				toteffgr_lowSel->GetPoint(i,tmpq2,toteff); 
				erreff = toteffgr_lowSel->GetErrorY(i-1);
			}
			nocorr_err = TMath::Sqrt(TMath::Power(nocorr_err,2) + TMath::Power(erreff/toteff,2));
			double rel_corr_yield = Jpsi2mumuBr * nn / cur_jpsiyield / toteff;
			double staterr = errnn/nn;
			double normerr = jpsierr/jpsiyield;
			double totstaterr = TMath::Sqrt(TMath::Power(staterr,2) + TMath::Power(normerr,2));
			double tot_nocorr_err = TMath::Sqrt(TMath::Power(totstaterr,2) + TMath::Power(nocorr_err,2));
			double tot_corr_err = TMath::Sqrt(TMath::Power(Jpsi2mumuBr_err/Jpsi2mumuBr,2) + TMath::Power(pdfsys[q2],2) + TMath::Power(corr_err,2));
			double resulterr = rel_corr_yield * TMath::Sqrt(TMath::Power(tot_nocorr_err,2) + TMath::Power(tot_corr_err,2));

			if(q2 > 9)
			{
				if(print=="yield") cout << fixed << setprecision(1) << q2-errq2 << "-" << q2+errq2 << fixed << setprecision(3) << " \t & $" << rel_corr_yield*1.e3 << " \\pm " << rel_corr_yield*1.e3*staterr << "$ (stat) $\\pm" << rel_corr_yield*1.e3*nocorr_err << "$ (uncorr) $\\pm" << rel_corr_yield*1.e3*tot_corr_err << "$ (corr) \\\\" << endl;
				else cout << q2-errq2 << "-" << q2+errq2 << ": nevts = " << nn << "  +/-  " << errnn << endl;

				resultgr->SetPoint(p,q2,rel_corr_yield);
				resultgr->SetPointError(p,errq2,resulterr);
			}

			if(type=="LL") 
			{
				LL_uncorrplot->SetPoint(p,q2,rel_corr_yield);
				LL_uncorrplot->SetPointError(p,q2,rel_corr_yield*tot_nocorr_err);
				LL_norm->SetPoint(p,q2,normerr*rel_corr_yield);
				LL_stat->SetPoint(p,q2,staterr*rel_corr_yield);
				LL_ncsys->SetPoint(p,q2,nocorr_err*rel_corr_yield);
				LL_eff->SetPoint(p,q2,toteff);
				LL_eff->SetPointError(p,q2,toteff*nocorr_err);
				corr_sys_gr->SetPoint(p,q2,corr_err);
				corr_sys_gr->SetPointError(p,errq2,0.);
			}
			else
			{
				DD_uncorrplot->SetPoint(p,q2,rel_corr_yield);
				DD_uncorrplot->SetPointError(p,q2,rel_corr_yield*tot_nocorr_err);
				DD_norm->SetPoint(p,q2,normerr*rel_corr_yield);
				DD_stat->SetPoint(p,q2,staterr*rel_corr_yield);
				DD_ncsys->SetPoint(p,q2,nocorr_err*rel_corr_yield);
				DD_eff->SetPoint(p,q2,toteff);
				DD_eff->SetPointError(p,q2,toteff*nocorr_err);	
			}
			p++;
		}

		histFile->cd();
		cc->cd();
		resultgr->Write("corrq2_"+type);
		resultgr->SetMarkerStyle(22);
		resultgr->SetMarkerSize(0.9);
		resultgr->SetMarkerColor(1);
		resultgr->Draw("AP");
		resultgr->GetYaxis()->SetTitle("(1 / B(#Lambda_{b} #rightarrow J/#psi #Lambda)) dB(#Lambda_{b} #rightarrow #Lambda #mu #mu) / dq^{2}");
		resultgr->GetYaxis()->SetTitleSize(0.04);
		resultgr->GetYaxis()->SetTitleOffset(1.4);
		resultgr->GetXaxis()->SetTitle("q^{2} (GeV^{2})");	
		resultgr->GetYaxis()->SetLabelSize(0.04);
		resultgr->SetMinimum(0.);
		cc->Print("q2result_"+type+".pdf");
		if(type=="LL") LLres = resultgr;
		else  DDres = resultgr;
	}


	LLres->SetMarkerColor(4);
	DDres->SetMarkerColor(2);

	TLegend * leg = new TLegend(0.18,0.9,0.3,0.75);
	leg->AddEntry(LLres,"LL","P");
	leg->AddEntry(DDres,"DD","P");

	DDres->Draw("AP");
	LLres->Draw("Psame");
	//leg->Draw("same");
	cc->Print("q2result_both.pdf");

	if(print=="yield") cout << "\n\nq2 bin \t\t & Combined yield \\\\" << endl; 
	else cout << "\n\nq2 bin \t\t & Eff DD \t & Uncorr err \t & Eff LL \t & Uncorr err \t & Corr err \\\\" << endl;

	double chi2 = 0;
	int point = 1;

	for(int i = 1; i < corr_sys_gr->GetN(); i++)
	{
		double q2, errq2, corr_sys, DDstat, LLstat;
		double DDncsys, LLncsys, DDnorm, LLnorm;
		double DDyield, LLyield, DDerr, LLerr;
		double DDeff, LLeff, DDeff_err, LLeff_err;

		DD_ncsys->GetPoint(i,q2,DDncsys);
		LL_ncsys->GetPoint(i,q2,LLncsys);
		DD_stat->GetPoint(i,q2,DDstat);
		LL_stat->GetPoint(i,q2,LLstat);
		LL_norm->GetPoint(i,q2,LLnorm);
		DD_norm->GetPoint(i,q2,DDnorm);
		corr_sys_gr->GetPoint(i,q2,corr_sys);
		errq2 = corr_sys_gr->GetErrorX(i);
		DD_uncorrplot->GetPoint(i,q2,DDyield);
		DDerr = DD_uncorrplot->GetErrorY(i);
		LL_uncorrplot->GetPoint(i,q2,LLyield);
		LLerr = LL_uncorrplot->GetErrorY(i);
		DD_eff->GetPoint(i,q2,DDeff);
		DDeff_err = DD_eff->GetErrorY(i);
		LL_eff->GetPoint(i,q2,LLeff);
		LLeff_err = LL_eff->GetErrorY(i);
		double sumW = 1./(DDerr*DDerr) + 1./(LLerr*LLerr);
		double combYield = (DDyield/(DDerr*DDerr) + LLyield/(LLerr*LLerr)) / sumW;
		double combErr = 1./TMath::Sqrt(sumW);
		double combStat = 1./TMath::Sqrt(1./(DDstat*DDstat) + 1./(LLstat*LLstat));
		double combNorm = 1./TMath::Sqrt(1./(DDnorm*DDnorm) + 1./(LLnorm*LLnorm));
		double noCorrSys = 1./TMath::Sqrt(1./TMath::Power(DDyield*DDncsys,2) + 1./TMath::Power(LLyield*LLncsys,2));
		double totCorrErr = TMath::Sqrt(TMath::Power(corr_sys,2) + TMath::Power(pdfsys[q2],2) + TMath::Power(Jpsi2mumuBr_err/Jpsi2mumuBr,2));
		double combErrPlusUncorr = combYield * TMath::Sqrt(TMath::Power(totCorrErr,2) + TMath::Power(combErr/combYield,2));
		double totSys = TMath::Sqrt(TMath::Power(combYield*totCorrErr,2) + TMath::Power(noCorrSys,2));

		if(q2 > 1.e10 || q2 < 1.e-6 || combYield < 1.e-6 || combYield > 1.e10) continue;
		else
		{
			if(print!="yield") cout << fixed << setprecision(1) << q2-errq2 << "-" << q2+errq2 << fixed << setprecision(3) << " \t & " << DDeff << " \t & " << DDeff_err << " \t & " << LLeff << " \t & " << LLeff_err << " \t & " << totCorrErr*100 << "\\% \\\\" << endl;
			
			if(q2>10)
			{
				if(print=="yield") cout << fixed << setprecision(1) << q2-errq2 << "-" << q2+errq2 << fixed << setprecision(5) << " \t & $" << combYield*1.e3 << " \\pm " << combStat*1.e3 << "$ (stat) $\\pm " << combNorm*1.e3 << "$ (norm) $\\pm " << 1.e3*totSys << "$ (sys) \\\\" << endl; 
				
				if(q2!=17.5) chi2 += TMath::Power(DDyield - LLyield,2)/(DDerr*DDerr + LLerr*LLerr);
				
				combresgr->SetPoint(point,q2,combYield);
				combresgr->SetPointError(point,errq2,combErrPlusUncorr);
				point++;
			}
		}
	}

	cout << endl << endl << "Consistency Chi2 = " << chi2/4. << "  pvalue = " << TMath::Prob(chi2,4.) << endl;

	histFile->cd();
	cc->cd();
	combresgr->Write("combresult");
	combresgr->SetMarkerStyle(22);
	combresgr->SetMarkerSize(0.9);
	combresgr->SetMarkerColor(1);
	combresgr->Draw("AP");
	combresgr->GetYaxis()->SetTitle("(1 / B(#Lambda_{b} #rightarrow J/#psi #Lambda))  dB(#Lambda_{b} #rightarrow #Lambda #mu #mu) / dq^{2}");
	combresgr->GetYaxis()->SetTitleSize(0.04);
	combresgr->GetYaxis()->SetTitleOffset(1.6);
	combresgr->GetXaxis()->SetTitle("q^{2} (GeV^{2})");
	cc->Print("combined_result.pdf");

	histFile->Write();
	histFile->Close();	
	delete histFile;

	return 0;
}





