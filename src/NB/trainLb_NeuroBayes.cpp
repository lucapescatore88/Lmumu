#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TRandom3.h"

#include "NeuroBayesTeacher.hh"
#include "NeuroBayesExpert.hh"

#include <iostream>
#include <string>

#include "NBfunctions.hpp"
#include "ReadTree_comp.hpp"
#include "general_functions.hpp"


int main(int argc, char** argv)
{
	string dataType = "12";
	//if(argc > 1) if((string)argv[1] == "11") dataType = "11";
	
	NeuroBayesTeacher* nb = NeuroBayesTeacher::Instance();
	
	nb->NB_DEF_TASK("CLASSIFICATION");
	
	//setup network topology
	int nvar = 20;  // Set this to number of inputs to your NN
	

	char ** varnames = new char*[nvar];
	varnames[0]  = "chi2_DTF";
	varnames[1]  = "Lb_TAU";
	varnames[2]  = "Lb_DIRA_OWNPV";
	varnames[3]  = "Lb_IPCHI2_OWNPV";
	varnames[4]  = "max_mu_IPCHI2_OWNPV";
	varnames[5]  = "min_mu_TRACKCHI2";
	varnames[6]  = "min_mu_PID";
	varnames[7]  = "min_mu_PID";
	varnames[8]  = "LL_Lambda0_IPCHI2_OWNPV";
	varnames[9] = "LL_Lambda0_FDCHI2_OWNPV";
	varnames[10] = "LL_Lambda0_PT";
	varnames[11] = "DD_Lambda0_IPCHI2_OWNPV";
	varnames[12] = "DD_Lambda0_FDCHI2_OWNPV";
	varnames[13] = "DD_Lambda0_PT";
	varnames[14]  = "DD_pplus_IPCHI2_OWNPV";
	varnames[15]  = "DD_piminus_IPCHI2_OWNPV";
	varnames[16]  = "DD_piminus_PT";
	varnames[17]  = "LL_pplus_IPCHI2_OWNPV";
	varnames[18]  = "LL_piminus_IPCHI2_OWNPV";
	varnames[19]  = "LL_piminus_PT";

	nb->NB_DEF_NODE1(nvar+1);
	nb->NB_DEF_NODE2(nvar);      // nodes in hidden layer 
	nb->NB_DEF_NODE3(1);       // nodes in output layer

	nb->NB_DEF_TASK("CLA");    // binominal classification

	nb->NB_DEF_PRE(822);
	//  nb->NB_DEF_PRE(812);
	nb->NB_DEF_REG("REG");           // 'OFF','REG' (def) ,'ARD','ASR','ALL'
	nb->NB_DEF_LOSS("ENTROPY");      // 'ENTROPY'(def),'QUADRATIC'
	nb->NB_DEF_METHOD("BFGS");
	nb->NB_DEF_SHAPE("DIAG");
	nb->NB_DEF_LEARNDIAG(1);

	nb->NB_DEF_RTRAIN(1.0);          // use 70% of events for training
	//  nb->NB_DEF_EPOCH(200);           // weight update after n events

	nb->NB_DEF_SPEED(2.0);           // multiplicative factor to enhance global learning speed
	nb->NB_DEF_MAXLEARN(1.0);        // multiplicative factor to limit the global learning speed in any direction, this number should be smaller than NB_DEF_SPEED

	nb->NB_DEF_ITER(100);             // number of training iteration
	//nb->NB_DEF_ITER(0);             // number of training iteration

    //int i = 4701;
	//int j = 29; 
	//nb->NB_RANVIN(i,j,2);            // random number seed initialisation, i has to be an odd number, the third argument is a debugging flag

	nb->SetOutputFile(("expert_"+dataType+".nb").c_str());  // expert file
	SetupNNPrepro(nb);
	
	// MC
	TreeReader* reader = new TreeReader("tree");
	reader->AddFile("/afs/cern.ch/work/p/pluca/Lmumu/weighted/Lb2Lmumu_MC_Pythia8_NBweighted_new.root");
	reader->Initialize();

	
	// We take all signal and 20% of background
	
	nb->SetTarget(1);
	int ntot = reader->GetEntries();
	int npassedMC = 0;
	
	cout << "Read in " << ntot << " events" << endl;
	
	int nstepMC = 5;
	//if(dataType=="11") nstepMC = 5;

	TFile ofile("/afs/cern.ch/work/p/pluca/Lmumu/weighted/samplesMVA_"+(TString)dataType+".root","recreate");
	TTree * sigTrainSample = new TTree("sigTrainSample","");
	reader->BranchNewTree(sigTrainSample);
	TTree * sigTestSample = new TTree("sigTestSample","");
	reader->BranchNewTree(sigTestSample);

	for(int event = 0; event < ntot; event++)
	{
		reader->GetEntry(event);
		
		if( TrueID(reader) && TriggerPassed(reader))
		{
			if( event%nstepMC==0 && npassedMC <= 4e4 )
			{
				npassedMC++;
				float InputArray[nvar+1];
				fillInputArray(reader,InputArray);
				if(isnan(InputArray[0])) continue;
				nb->SetWeight(reader->GetValue("Lb_weight"));
				nb->SetNextInput(nvar,InputArray);
				sigTrainSample->Fill();
			}
			else sigTestSample->Fill();
		}
	}


	// Data
	TreeReader* reader2 = new TreeReader("tree");
	reader2->AddFile("/afs/cern.ch/work/p/pluca/Lmumu/weighted/Lb2Lmumu_CL_NBweighted.root");
	reader2->Initialize();

	TTree * bkgTrainSample = new TTree("bkgTrainSample","");
	reader2->BranchNewTree(bkgTrainSample);
	TTree * bkgTestSample = new TTree("bkgTestSample","");
	reader2->BranchNewTree(bkgTestSample);

	nb->SetTarget(0);
	int ntot2 = reader2->GetEntries();
	int npassed = 0;
	
	cout << "Read in " << ntot2 << " events" << endl;
	int nstep = 2;
	//if(dataType=="11") nstep = 2;
	
	for(int event = 0; event < ntot2; event++)
	{
		reader2->GetEntry(event);
		
		double massLb = reader2->GetValue("Lb_MassConsLambda_M",0);
		double massJpsi = reader2->GetValue("J_psi_1S_MM");
		if(massLb > 6000 && TMath::Abs(massJpsi - 3096) > 100 && TMath::Abs(massJpsi - 3686) > 90 && TriggerPassed(reader2))
		{
			if(event%nstep==0 && npassed <=4e4)
			{
				float InputArray[100];
				npassed++;
				fillInputArray(reader2,InputArray);
				if(isnan(InputArray[0])) continue;
				nb->SetWeight(1.);
				nb->SetNextInput(nvar,InputArray);
				bkgTrainSample->Fill();
			}
			else bkgTestSample->Fill();
		}
	}
	
	bkgTrainSample->Write();
	bkgTestSample->Write();
	sigTestSample->Write();
	sigTrainSample->Write();
	ofile.Close();

	cout << "\nData used = " << npassed << ", MC used = " << npassedMC << endl;
	
	cout << "Train the Network\n" << endl;
	nb->TrainNet();
	nb->nb_correl_signi(varnames,"correl_signi.txt","correl_signi.html");

	cout << "\nData used = " << npassed << ", MC used = " << npassedMC << endl;
	
	return 0;
}

