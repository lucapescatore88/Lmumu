#include "ReadTree_comp.hpp"
#include "NeuroBayesTeacher.hh"
#include <cmath>

void SetupNNPrepro(NeuroBayesTeacher* nb)
{
	nb->SetIndividualPreproFlag(0,14);
	nb->SetIndividualPreproFlag(1,14);
	nb->SetIndividualPreproFlag(2,14);
	nb->SetIndividualPreproFlag(3,14);
	nb->SetIndividualPreproFlag(4,14);
	nb->SetIndividualPreproFlag(5,14);
	nb->SetIndividualPreproFlag(6,14);
	nb->SetIndividualPreproFlag(7,14);

	nb->SetIndividualPreproFlag(8,34);
	nb->SetIndividualPreproFlag(9,34);
	nb->SetIndividualPreproFlag(10,34);
	nb->SetIndividualPreproFlag(11,34);
	nb->SetIndividualPreproFlag(12,34);
	nb->SetIndividualPreproFlag(13,34);
	nb->SetIndividualPreproFlag(14,34);
	nb->SetIndividualPreproFlag(15,34);
	nb->SetIndividualPreproFlag(16,34);
	nb->SetIndividualPreproFlag(17,34);
	nb->SetIndividualPreproFlag(18,34);
	nb->SetIndividualPreproFlag(19,34);
}



void fillInputArray(TreeReader* reader, float* array)
{
	float chi2DTF = reader->GetValue( "Lb_MassConsLambda_chi2",0) / reader->GetValue( "Lb_MassConsLambda_nDOF",0);
	if(!std::isfinite(chi2DTF)) chi2DTF = 1.e6;
	array[0] = chi2DTF;
	array[1] = reader->GetValue( "Lb_TAU");
	array[2] = reader->GetValue( "Lb_DIRA_OWNPV");	
	array[3] = reader->GetValue( "Lb_IPCHI2_OWNPV");
	
	float mu1IPchi2 = reader->GetValue( "muminus_IPCHI2_OWNPV");
	float mu2IPchi2 = reader->GetValue( "muplus_IPCHI2_OWNPV");
	
	array[4] = TMath::Max(mu1IPchi2,mu2IPchi2);
	array[5] = TMath::Min(mu1IPchi2,mu2IPchi2);

	//float pidMu1 = reader->GetValue( "muminus_ProbNNmu");
	//float pidMu2 = reader->GetValue( "muplus_ProbNNmu");
	
	float pidMu1 = reader->GetValue( "muminus_PIDmu");
	float pidMu2 = reader->GetValue( "muplus_PIDmu");

	array[6] = TMath::Min(pidMu1,pidMu2);
	array[7] = TMath::Max(pidMu1,pidMu2);

	int trackType = reader->GetValue( "pplus_TRACK_Type");
	if(trackType == 3)
	{
		array[8] = reader->GetValue( "Lambda0_IPCHI2_OWNPV");
		array[9] = reader->GetValue( "Lambda0_FD_OWNPV");
		array[10] = reader->GetValue( "Lambda0_PT");
		array[11] = -999;
		array[12] = -999;
		array[13] = -999;
		array[14] = reader->GetValue( "piminus_IPCHI2_OWNPV");
		array[15] = reader->GetValue( "pplus_IPCHI2_OWNPV");
		array[16] = reader->GetValue( "piminus_PT");
		array[17] = -999;
		array[18] = -999;
		array[19] = -999;
	}
	else
	{
		array[8] = -999;
		array[9] = -999;
		array[10] = -999;
		array[11] = reader->GetValue( "Lambda0_IPCHI2_OWNPV");
		array[12] = reader->GetValue( "Lambda0_FD_OWNPV");
		array[13] = reader->GetValue( "Lambda0_PT");
		array[14] = -999;
		array[15] = -999;
		array[16] = -999;
		array[17] = reader->GetValue( "piminus_IPCHI2_OWNPV");
		array[18] = reader->GetValue( "pplus_IPCHI2_OWNPV");
		array[19] = reader->GetValue( "piminus_PT");
	}
}


bool TrueID(TreeReader *reader, string type = "mumu", bool fullid = true);

bool TrueID(TreeReader *reader, string type,  bool fullid )
{
	Int_t Lambda0_mother = TMath::Abs(reader->GetValue("Lambda0_MC_MOTHER_ID"));
	Int_t piminus_mother = TMath::Abs(reader->GetValue("piminus_MC_MOTHER_ID"));
	Int_t pplus_mother = TMath::Abs(reader->GetValue("pplus_MC_MOTHER_ID"));
	Int_t muminus_mother = TMath::Abs(reader->GetValue("muminus_MC_MOTHER_ID"));
	Int_t muplus_mother = TMath::Abs(reader->GetValue("muplus_MC_MOTHER_ID"));
	Int_t dimuon_mother = TMath::Abs(reader->GetValue("J_psi_1S_MC_MOTHER_ID"));
	
	Int_t Lb_ID = TMath::Abs(reader->GetValue("Lb_TRUEID"));
	Int_t Lambda0_ID = TMath::Abs(reader->GetValue("Lambda0_TRUEID"));
	Int_t piminus_ID = TMath::Abs(reader->GetValue("piminus_TRUEID"));
	Int_t pplus_ID = TMath::Abs(reader->GetValue("pplus_TRUEID"));
	Int_t muminus_ID = TMath::Abs(reader->GetValue("muminus_TRUEID"));
	Int_t muplus_ID = TMath::Abs(reader->GetValue("muplus_TRUEID"));
	Int_t dimuon_ID = TMath::Abs(reader->GetValue("J_psi_1S_TRUEID"));

	bool res;
	if(type=="jpsi")
	{
		res = (Lambda0_mother == 5122 && piminus_mother == 3122 && pplus_mother == 3122 && muminus_mother == 443 && muplus_mother == 443 && dimuon_mother == 5122);
		if (fullid) res = (res && Lb_ID == 5122 && Lambda0_ID == 3122 && piminus_ID == 211 && pplus_ID == 2212 && muminus_ID == 13 && muplus_ID == 13 && dimuon_ID == 443);
	}
	else
	{
		res = (Lambda0_mother == 5122 && piminus_mother == 3122 && pplus_mother == 3122 && muminus_mother == 5122 && muplus_mother == 5122);
		if (fullid) res = (res && Lb_ID == 5122 && Lambda0_ID == 3122 && piminus_ID == 211 && pplus_ID == 2212 && muminus_ID == 13 && muplus_ID == 13);
	}
	
	return res;
}


bool isCharm(TreeReader *reader)
{
	double mass = reader->GetValue("J_psi_1S_MM");
	if(TMath::Abs(mass - 3096.916) < 100 || TMath::Abs(mass - 3686.1) < 90)
		return true;
	else return false;
}


bool TriggerPassed(TreeReader *reader)
{
	Bool_t L0Muon_TOS = reader->GetValue("Lb_L0MuonDecision_TOS");
	Bool_t L0DiMuon_TOS = reader->GetValue("Lb_L0DiMuonDecision_TOS");

	Bool_t HLT1All_TOS = reader->GetValue("Lb_Hlt1TrackAllL0Decision_TOS");
	Bool_t HLT1Mu_TOS = reader->GetValue("Lb_Hlt1TrackMuonDecision_TOS");
	Bool_t HLT1DiMu_TOS = reader->GetValue("Lb_Hlt1DiMuonHighMassDecision_TOS");

	Bool_t HLT2Mu2_TOS = reader->GetValue("Lb_Hlt2TopoMu2BodyBBDTDecision_TOS");
	Bool_t HLT2Mu3_TOS = reader->GetValue("Lb_Hlt2TopoMu3BodyBBDTDecision_TOS");
	Bool_t HLT2Mu4_TOS = reader->GetValue("Lb_Hlt2TopoMu4BodyBBDTDecision_TOS");
	Bool_t HLT2Topo2_TOS = reader->GetValue("Lb_Hlt2Topo2BodyBBDTDecision_TOS");
	Bool_t HLT2Topo3_TOS = reader->GetValue("Lb_Hlt2Topo3BodyBBDTDecision_TOS");
	Bool_t HLT2Topo4_TOS = reader->GetValue("Lb_Hlt2Topo4BodyBBDTDecision_TOS");
	Bool_t HLT2SM_TOS = reader->GetValue("Lb_Hlt2SingleMuonDecision_TOS");
	Bool_t HLT2DM_TOS = reader->GetValue("Lb_Hlt2DiMuonDetachedDecision_TOS");
	
	
	if (L0Muon_TOS == true || L0DiMuon_TOS == true)
	{
		if(HLT1DiMu_TOS == true || HLT1All_TOS == true || HLT1Mu_TOS == true)	
		{
			if(HLT2Topo2_TOS == true || HLT2Topo3_TOS == true || HLT2Topo4_TOS == true || HLT2Mu2_TOS == true
						|| HLT2Mu3_TOS == true || HLT2Mu4_TOS == true || HLT2SM_TOS == true || HLT2DM_TOS == true)
			{
				return true;
			}
			else return false;
		}
		else return false;
	}
	else return false;
}
