#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TList.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom3.h"

#include "general_functions.hpp"
#include "ReadTree_comp.hpp"
#include "NBfunctions.hpp"


#include "DPHelpers.hpp"
#include "LbLmumuAng.hh"
#include "LbJpsiL2.hh"
#include "LbJpsiLAng.hh"

#include "NeuroBayesExpert.hh"

#include <iostream>
#include <string>
#include <cstring>


TLorentzVector get4momentum(TreeReader * reader, string pname);

TTree * appWieght(TreeReader * reader, vector<float> MCnorm = vector<float>(), string tname = "",  bool physW = false, bool mva = true, bool Lnodecay = false);

TH2F * weightHist = NULL;
TH2F * PTweightHist = NULL;
TH1F * pidWhist = NULL;

int main(int argc, char** argv)
{
	// Arguments:
	//  1 - treename
	//  2 - output file name
	//  3 - list of input files
	
	TreeReader * reader = new TreeReader(argv[1]);

	for(int a = 3; a < argc; a++)
	{
		reader->AddFile(argv[a]);
		cout << "Reading: " << argv[a] << endl;
	}
	cout << "Done" << endl;

	TFile * weightFile = TFile::Open("/afs/cern.ch/user/p/pluca/work/Lb/Lmumu/MC_reweight.root");
	weightHist = (TH2F *)weightFile->Get("ratio2D"); 
	PTweightHist = (TH2F *)weightFile->Get("ratio2D_pt");
	pidWhist = (TH1F *)weightFile->Get("ratio_PIDmu");

	vector<string> novar(1,"weight");
	novar.push_back("Lb_weight");
	novar.push_back("pid_muplus");
	novar.push_back("pid_muminus");
	novar.push_back("pt_weight");

	reader->Initialize(novar,"except");

	TFile out(argv[2],"recreate");

	TTree* newTree = new TTree(argv[1],"");
	reader->BranchNewTree(newTree);

	string expertfile = "expert_12.nb";
	Expert* nbExpert = new Expert(expertfile.c_str());

	double isMC = reader->HasVar("Lb_TRUEPT");

	float NNOut, Lb_weight, pt_weight, pidWeight_plus, pidWeight_minus;
	newTree->Branch("weight",&NNOut,"weight/F");
	if(isMC) newTree->Branch("Lb_weight",&Lb_weight,"Lb_weight/F");
	if(isMC) newTree->Branch("pt_weight",&pt_weight,"pt_weight/F");
	if(isMC) newTree->Branch("pid_muplus",&pidWeight_plus,"pid_muplus/F");
	if(isMC) newTree->Branch("pid_muminus",&pidWeight_minus,"pid_muminus/F");

	int ntot = reader->GetEntries();
	for( int i = 0; i < ntot; ++i )
	{
		showPercentage(i,ntot);
		reader->GetEntry(i);

		string filename = (string)reader->GetChain()->GetFile()->GetName();

		if(filename.find("CL")!=string::npos && reader->GetValue("Lb_MM") < 5000) continue;

		float inputs[100];
		fillInputArray(reader,inputs);
		if(isnan(inputs[0])) continue;

		NNOut = nbExpert->nb_expert(inputs);
		NNOut = 0.5*(NNOut+1);

		float pidplus = reader->GetValue("muplus_PIDmu");
		float pidminus = reader->GetValue("muminus_PIDmu");
		pidWeight_plus = pidWhist->GetBinContent(pidWhist->GetXaxis()->FindBin(pidplus));
		pidWeight_minus = pidWhist->GetBinContent(pidWhist->GetXaxis()->FindBin(pidminus));

		if(isMC)
		{
			//Lb momentum reweight
				
			float Lbpx = reader->GetValue("Lb_TRUEP_X");
			float Lbpy = reader->GetValue("Lb_TRUEP_Y");
			float Lbpz = reader->GetValue("Lb_TRUEP_Z");
			float Lbp = TMath::Sqrt(TMath::Power(Lbpx,2) + TMath::Power(Lbpy,2) + TMath::Power(Lbpz,2));
			float Lbpt = reader->GetValue("Lb_TRUEPT");
			Lb_weight = weightHist->GetBinContent(weightHist->GetXaxis()->FindBin(Lbpt),weightHist->GetYaxis()->FindBin(Lbp));
			
			float Lambda0pt = reader->GetValue("Lambda0_TRUEPT");
			pt_weight = PTweightHist->GetBinContent(PTweightHist->GetXaxis()->FindBin(Lbpt),PTweightHist->GetYaxis()->FindBin(Lambda0pt));
		}

		newTree->Fill();
	}

	newTree->Write();
	out.Close();

	delete newTree;
}



TLorentzVector get4momentum(TreeReader * reader, string pname)
{
	TLorentzVector p;
	string base = pname + "_TRUEP";
	p.SetPxPyPzE(reader->GetValue((base+"_X").c_str()),
			reader->GetValue((base+"_Y").c_str()),
			reader->GetValue((base+"_Z").c_str()),
			reader->GetValue((base+"_E").c_str()));
	return p;
}





