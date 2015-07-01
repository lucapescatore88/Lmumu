#include "analyser.hpp"
#include "TCanvas.h"
#include "Lb_cuts.hpp"
#include "TGraphErrors.h"
#include "RooGenericPdf.h"
#include "RooMsgService.h"

int main(int argc, char **argv)
{
	TCanvas * c = new TCanvas();
	TString var = "cosThetaL";
	TString type = "All";
	bool dosys = false;

	if(argc > 1)
	{
		for(int a = 1; a < argc; a++)
		{
			string arg = argv[a];
			string str = arg.substr(2,arg.length()-2);

			if(arg.find("-t") != string::npos) type = (TString)str;
			if(arg.find("-v") != string::npos) var = (TString)str;
			if(arg=="-sys") dosys = true;
		}
	}

	int min = -1;
	int max = 1;
	double nstep = 10;

	TCut baseCut = "";
	if(type == "DD") { baseCut = CutsDef::DDcut; }
	else if(type == "LL") { baseCut = CutsDef::LLcut; }

	TFile *_file0 = TFile::Open("/afs/cern.ch/user/p/pluca/work/Lb/Lmumu/candLb_MC.root");
	TTree * MCtree = (TTree*)_file0->Get("candLb2Lmumu");

	TString q2str = "TMath::Power(J_psi_1S_MM/1000,2)";
	//int q2nbins = CutsDef::nq2bins;
	//double * q2min = &CutsDef::q2min_highfirst[0];
	//double * q2max = &CutsDef::q2max_highfirst[0];
	TGraphErrors * q2Res = new TGraphErrors();

	int q2nbins = 5;
	double q2min[] = {11.0, 15.0, 16.0, 18.0, 15.0};
	double q2max[] = {12.5, 16.0, 18.0, 20.0, 20.0};


	if(!dosys)
	{
		gStyle->SetOptStat(0);
		MCtree->Draw(var+":"+var+"_TRUE",baseCut);
		c->Print(var+"_vs_"+var+"_TRUE_"+type+".pdf");
		MCtree->Draw(var+"-"+var+"_TRUE:"+var,baseCut);
		c->Print("RmT_vs_"+var+"_"+type+".pdf");
		MCtree->Draw(var+"-"+var+"_TRUE:TMath::Power(J_psi_1S_MM/1000,2)>>hRmTvsq2",baseCut);
		TH1 * hRmTvsq2 = (TH1*)gPad->GetPrimitive("hRmTvsq2");
		hRmTvsq2->GetXaxis()->SetTitle("q^{2}");
		c->Print("RmT"+var+"_vs_q2_"+type+".pdf");

		gStyle->SetOptStat("mer");
		TGraphErrors * cosThetaRes = new TGraphErrors();

		for(int k = 0; k < nstep; k++)
		{
			TString hname = Form("h%i",k);
			double mymin = min + k*(max-min)/nstep;
			double mymax = min + (k+1.)*(max-min)/nstep;
			MCtree->Draw(var+"-"+var+"_TRUE>>"+hname+"(50,-0.25,0.25)",(TCut)(var+Form(" < %f && ",mymax) + var + Form(" > %f",mymin))+baseCut);
			TH1 * h = (TH1*)gPad->GetPrimitive(hname);
			TF1 * mygaus = new TF1("mygaus","gaus");
			mygaus->FixParameter(1,0.);
			h->Fit("mygaus","QE");
			TF1 *fit = h->GetFunction("mygaus");
			TString myminstr = ((TString)Form("%4.2f",mymin)).ReplaceAll(".","");
			TString mymaxstr = ((TString)Form("%4.2f",mymax)).ReplaceAll(".","");
			myminstr.ReplaceAll("-","m");
			mymaxstr.ReplaceAll("-","m");
			h->GetXaxis()->SetTitle(var);
			c->Print("RmT"+var+"_"+myminstr+"-"+mymaxstr+"_"+type+".pdf");
			cosThetaRes->SetPoint(k,(mymax+mymin)/2.,fit->GetParameter(2));
			cosThetaRes->SetPointError(k,(mymax-mymin)/2.,fit->GetParError(2));
			delete mygaus;
			delete fit;
			delete h;
		}

		cosThetaRes->Draw("AP");
		cosThetaRes->GetXaxis()->SetTitle(var);
		cosThetaRes->GetYaxis()->SetTitle(var + " resolution");
		cosThetaRes->SetMarkerStyle(22);
		cosThetaRes->SetMarkerColor(2);
		cosThetaRes->SetMarkerSize(1);
		cosThetaRes->SetMinimum(0.);
		c->Print("resolution_"+var+"_"+type+".pdf");


		double q2binning[] = {0.1, 2.0, 4.0, 6.0, 11.0, 12.5, 15.0, 16.0, 18.0, 20.0};
		double cosThetabinning[] = {-1,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0};

		TH2F * h2D = new TH2F("2dresolution","2Dresolution",9,q2binning,10,cosThetabinning);
		for(int q2 = 0; q2 < 9; q2++)
		{
			if(q2binning[q2] == 6.0 || q2binning[q2] == 12.5) continue;
			for(int k = 0; k < nstep; k++)
			{
				TString hname = Form("h2d%i",q2);
				double mymin = q2binning[q2];
				double mymax = q2binning[q2+1];
				double myminvar = min + k*(max-min)/nstep;
				double mymaxvar = min + (k+1.)*(max-min)/nstep;
				MCtree->Draw(var+"-"+var+"_TRUE>>"+hname+"(50,-0.25,0.25)",(TCut)(q2str+Form(" < %f && ",mymax) + q2str + Form(" > %f",mymin))+(TCut)(var+Form(" < %f && ",mymaxvar) + q2str + Form(" > %f",myminvar))+baseCut);
				TH1 * h = (TH1*)gPad->GetPrimitive(hname);
				TF1 * mygaus = new TF1("mygaus","gaus");
				mygaus->FixParameter(1,0.);
				h->Fit("mygaus","QE");
				TF1 *fit = h->GetFunction("mygaus");
				h2D->Fill((mymin+mymax)/2.,(myminvar+mymaxvar)/2.,fit->GetParameter(2));
				delete mygaus;
				delete fit;
				delete h;
			}
		}

		gStyle->SetOptStat(0);
		//h2D->SetMinimum(0.);
		h2D->Draw("colz");
		h2D->GetXaxis()->SetTitle("q^{2}");
		h2D->GetYaxis()->SetTitle(var);
		h2D->SetTitle(var + " resolution");
		c->Print("resolution2D_"+var+"_"+type+".pdf");

		TF1 * mygaus = new TF1("mygaus","gaus");
		MCtree->Draw(var+"-"+var+"_TRUE>>havgResol",baseCut);
		TH1 * h = (TH1*)gPad->GetPrimitive("havgResol");
		h->Fit("mygaus","QE");
		TF1 *fit = h->GetFunction("mygaus");
		double resol = fit->GetParameter(2);

		cout << "************************** Average "+var+ " resolution ("+type+") = " << resol << endl;
	}

	//Resoltution systematic

	if(dosys)
	{
		cout << fixed << setprecision(5);
		TFile * outfile = new TFile("Afbsys_"+var+".root","recreate");
		double nevts[] = { 34., 41., 88., 72., 200 };

		TGraphErrors * hq2sys = new TGraphErrors();
		TGraphErrors * hq2sysfL = new TGraphErrors();

		RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
		TString effbase = "/afs/cern.ch/user/p/pluca/work/Lb/Lmumu/results/";
		TFile * effFile = TFile::Open(effbase+"LbeffvscosThetaL_LL.root");
		TH1F * effLL  = (TH1F *)effFile->Get("htoteff");
		effFile = TFile::Open(effbase+"LbeffvscosThetaB_LL.root");
		TH1F * effLLB  = (TH1F *)effFile->Get("htot_nodet_eff");

		RooRealVar * cosThetaL = new RooRealVar("cosThetaL","cosThetaL",0.,-1.,1.);
		RooRealVar * cosThetaB = new RooRealVar("cosThetaB","cosThetaB",0.,-1.,1.);

		RooDataHist * hLL = new RooDataHist("hLL","hLL",*cosThetaL,effLL);
		RooRealVar * c1LL = new RooRealVar("c1LL","",0.,-100.,100);
		RooRealVar * c2LL = new RooRealVar("c2LL","",0.,-100.,100);
		TString effLLstr = "(1 + c1LL*cosThetaL + c2LL*TMath::Power(cosThetaL,2))";
		RooAbsPdf * effLLpdf = new RooGenericPdf("effLLpdf", "", effLLstr, RooArgSet(*cosThetaL, *c1LL, *c2LL));
		effLLpdf->fitTo(*hLL,PrintLevel(-1));
		fixParams(effLLpdf,cosThetaL);

		RooDataHist * hLLB = new RooDataHist("hLLB","hLLB",*cosThetaB,effLLB);
		RooRealVar * cB1LL = new RooRealVar("cB1LL","",0,-100.,100);
		RooRealVar * cB2LL = new RooRealVar("cB2LL","",0,-100.,100);
		TString effLLBstr = "(1 + cB1LL*cosThetaB + cB2LL*TMath::Power(cosThetaB,2))";
		RooAbsPdf * effLLpdfB = new RooGenericPdf("effLLpdfB", "", effLLBstr, RooArgSet(*cosThetaB, *cB1LL, *cB2LL));
		effLLpdfB->fitTo(*hLLB,PrintLevel(-1));
		fixParams(effLLpdfB,cosThetaB);

		Analysis * normal, * smeared; 
		RooRealVar * afb = new RooRealVar("afb","afb",0.04,-100.,100.);
		RooRealVar * fL = new RooRealVar("fL","fL",0.82,-100.,100.);
		TString afbpdf = "((3./8.)*(1.-fL)*(1 + TMath::Power(cosThetaL,2)) + afb*cosThetaL + (3./4.)*fL*(1 - TMath::Power(cosThetaL,2)))";
		RooRealVar * afbB = new RooRealVar("afbB","afbB",-0.37,-100.,100.);
		TString afbBpdf = "(1 + 2*afbB*cosThetaB)";
		Analysis::SetPrintLevel("s");

		int p = 0;
		for(int q2bin = 0; q2bin < q2nbins; q2bin++)
		{
			TString hname = Form("hh%i",q2bin);
			double mymin = q2min[q2bin];
			double mymax = q2max[q2bin];
			MCtree->Draw(var+"-"+var+"_TRUE>>"+hname+"(50,-0.25,0.25)",(TCut)(q2str+Form(" < %f && ",mymax) + q2str + Form(" > %f",mymin))+"pplus_TRACK_Type == 5");
			TH1 * h = (TH1*)gPad->GetPrimitive(hname);
			TF1 * mygaus = new TF1("mygaus","gaus");
			mygaus->FixParameter(1,0.);
			h->Fit("mygaus","QE");
			TF1 *fit = h->GetFunction("mygaus");
			TString myminstr = ((TString)Form("%4.2f",mymin)).ReplaceAll(".","");
			TString mymaxstr = ((TString)Form("%4.2f",mymax)).ReplaceAll(".","");
			myminstr.ReplaceAll("-","m");
			mymaxstr.ReplaceAll("-","m");
			h->GetXaxis()->SetTitle("q_{2}");
			c->Print("RmT"+var+"_q2_"+myminstr+"-"+mymaxstr+"_"+type+".pdf");
			q2Res->SetPoint(q2bin,(mymax+mymin)/2.,fit->GetParameter(2));
			q2Res->SetPointError(q2bin,(mymax-mymin)/2.,fit->GetParError(2));

			TString smearval = Form("%4.2f",fit->GetParameter(2));
			cout << "************************************ Examining: " << mymin << " - " << mymax << "  Nevts = " << nevts[q2bin] << endl;
			TH1F * hsys = new TH1F(Form("hsys_%i",q2bin),"",20,-0.5,0.5); 	
			TH1F * hfLsys = new TH1F(Form("hfLsys_%i",q2bin),"",20,-0.5,0.5);

			for(int e = 0; e < 1000; e++)
			{
				showPercentage(e,1000);
				RooAbsPdf * corrPdf;
				if(var=="cosThetaL")
				{
					corrPdf = new RooGenericPdf(Form("corrPdf%i_%i",e,q2bin),afbpdf+"*"+effLLstr,RooArgSet(*cosThetaL, *afb, *fL, *c1LL, *c2LL) );
					normal = new Analysis(Form("normal%i_%i",e,q2bin),cosThetaL);
					smeared = new Analysis(Form("smeared%i_%i",e,q2bin),cosThetaL);
				}
				else
				{
					corrPdf = new RooGenericPdf(Form("corrPdf%i_%i",e,q2bin),afbBpdf+"*"+effLLBstr,RooArgSet(*cosThetaB, *afbB, *cB1LL, *cB2LL) );
					normal = new Analysis(Form("normal%i_%i",e,q2bin),cosThetaB);
					smeared = new Analysis(Form("smeared%i_%i",e,q2bin),cosThetaB);
				}

				afb->setVal(0.0429);
				afbB->setVal(-0.37);
				fL->setVal(0.82);

				normal->SetModel(corrPdf);
				normal->Generate(nevts[q2bin],Form("-seed%i-subtree",e+1));

				smeared->SetModel(corrPdf);
				smeared->Generate(nevts[q2bin],Form("-seed%i-smear"+smearval+"-subtree",e+1));

				afb->setVal(0.);
				afbB->setVal(0.);
				fL->setVal(0.7);

				corrPdf->fitTo(*(normal->GetDataSet()),PrintLevel(-1));	
				double norm_val, smear_val, smear_fLval = 0, norm_fLval = 0;
				if(var=="cosThetaL") { norm_val = afb->getVal(); norm_fLval = fL->getVal(); }
				else norm_val = afbB->getVal();

				afb->setVal(0.);
				afbB->setVal(0.);
				fL->setVal(0.7);

				corrPdf->fitTo(*(smeared->GetDataSet()),PrintLevel(-1));
				if(var=="cosThetaL") { smear_val = afb->getVal(); smear_fLval = fL->getVal(); }
				else smear_val = afbB->getVal();

				double rel_dev = norm_val - smear_val;
				hsys->Fill(rel_dev);
				if(var=="cosThetaL") hfLsys->Fill(norm_fLval - smear_fLval);
				//cout << fixed << setprecision(5) << norm_val << "   " << smear_val << endl;
			}

			hq2sys->SetPoint(p,(q2min[q2bin]+q2max[q2bin])/2.,hsys->GetMean());
			if(var=="cosThetaL") hq2sysfL->SetPoint(p,(q2max[q2bin]+q2min[q2bin])/2.,hfLsys->GetMean());
			hq2sys->SetPointError(p,(q2max[q2bin]-q2min[q2bin])/2.,hsys->GetMeanError());
			if(var=="cosThetaL") hq2sysfL->SetPointError(p,(q2max[q2bin]-q2min[q2bin])/2.,hfLsys->GetMeanError());
			p++;
		}

		q2Res->Draw("AP");
		q2Res->GetXaxis()->SetTitle("q^{2}");
		q2Res->GetYaxis()->SetTitle(var + " resolution");
		q2Res->SetMarkerStyle(22);
		q2Res->SetMarkerColor(2);
		q2Res->SetMarkerSize(1);
		q2Res->SetMinimum(0.);
		c->Print("resolution_"+var+"vsq2_"+type+".pdf");

		for(int qq = 0; qq < q2Res->GetN(); qq++)
		{ 
			double x,y,ex,afbsys,fLsys = 0; 
			ex = q2Res->GetErrorX(qq);
			q2Res->GetPoint(qq,x,y);
			hq2sys->GetPoint(qq,x,afbsys);
			if(var=="cosThetaL") hq2sysfL->GetPoint(qq,x,fLsys);
			cout << fixed << setprecision(1) << x-ex << "-"<< x+ex << " & "  << fixed << setprecision(5) << y << " & " << afbsys;
			if(var=="cosThetaL") cout << " & " << fLsys;
			cout << " \\\\" << endl;
		}

		outfile->cd();
		hq2sys->Draw("AP");
		if(var=="cosThetaL") hq2sys->GetYaxis()->SetTitle("A_{FB}^{l}");
		hq2sys->GetYaxis()->SetTitle("A_{FB}^{B}");
		hq2sys->GetXaxis()->SetTitle("q^{2}");
		c->Print("Afbsys_resolution_"+var+".pdf");
		hq2sys->Write();
		if(var=="cosThetaL") {
			hq2sysfL->Draw("AP");
			hq2sysfL->GetYaxis()->SetTitle("A_{FB}^{B}");
			hq2sysfL->GetXaxis()->SetTitle("q^{2}");
		    c->Print("fLsys_resolution.pdf");
		    hq2sysfL->Write();
		}

		outfile->Close();
		delete outfile;
	}

	delete c;
	delete _file0;
}
