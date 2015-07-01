#include "TCanvas.h"
#include "TH1F.h"
#include "TCut.h"
#include "TFile.h"
#include "ReadTree_comp.hpp"
#include "TH1F.h"
#include "general_functions.hpp"
#include "TStyle.h"
#include "TMath.h"

using namespace std;

int main()
{
	gStyle->SetOptStat(0);	

	vector<TString> trig;
	trig.push_back( "Lb_Hlt2Topo2BodyBBDTDecision_TOS" );
	trig.push_back(	"Lb_Hlt2Topo3BodyBBDTDecision_TOS" );
	trig.push_back(	"Lb_Hlt2Topo4BodyBBDTDecision_TOS" );
	trig.push_back(	"Lb_Hlt2TopoMu2BodyBBDTDecision_TOS" );
	trig.push_back(	"Lb_Hlt2TopoMu3BodyBBDTDecision_TOS" );
	trig.push_back(	"Lb_Hlt2TopoMu4BodyBBDTDecision_TOS" );
	trig.push_back(	"Lb_Hlt2SingleMuonDecision_TOS" );
	trig.push_back(	"Lb_Hlt2DiMuonDetachedDecision_TOS" );

	TreeReader * reader = new TreeReader("tree");
	reader->AddFile("/afs/cern.ch/work/p/pluca/weighted/Lmumu/Lb2Lmumu_MC_Pythia8_NBweighted.root");
	reader->Initialize();

	TH1F * hout = new TH1F("trig","trig",trig.size(),0,trig.size());
	int ntot = reader->GetEntries();


	for(size_t t = 0; t < trig.size(); t++)	
	{
		TString label = trig[t];
		label.ReplaceAll("Lb_","");
		label.ReplaceAll("_TOS","");
		label.ReplaceAll("Decision","");
		hout->GetXaxis()->SetBinLabel(t+1,label);
	}

	int nevts = 0;
	for(int i = 0; i < ntot; i++)
	{
		showPercentage(i,ntot);
		reader->GetEntry(i);
		if(TMath::Power(reader->GetValue("J_psi_1S_MM")/1000,2) < 15 ) continue;

		nevts++;
		for(size_t t = 0; t < trig.size(); t++)
		{
			if(reader->GetValue(trig[t])) hout->Fill(t);
		}
	}

	TCanvas * c = new TCanvas();
	hout->GetYaxis()->SetTitle("single trig. eff.");
	hout->GetYaxis()->SetTitleOffset(1.25);
	hout->SetTitle(0);
	hout->Scale(1./(double)nevts);
	hout->Draw();
	c->Print("trig_highq2.pdf");
}


