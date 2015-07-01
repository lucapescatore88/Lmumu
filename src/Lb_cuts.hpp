#ifndef LB_CUTS_H
#define LB_CUTS_H

#include "TCut.h"
#include <iostream>

namespace CutsDef
{
	int nq2bins = 10;
	double q2min[] = {0.1, 2.0, 4.0, 6.0, 11.0, 15.0, 16.0, 18.0, 1.1, 15.0};
	double q2max[] = {2.0, 4.0, 6.0, 8.0, 12.5, 16.0, 18.0, 20.0, 6.0, 20.0};
	double q2min_highfirst[] = {15.0, 11.0, 15.0, 16.0, 18.0, 0.1, 2.0, 4.0, 6.0, 1.1};
	double q2max_highfirst[] = {20.0, 12.5, 16.0, 18.0, 20.0, 2.0, 4.0, 6.0, 8.0, 6.0};

	TCut highQ2 = "TMath::Power(J_psi_1S_MM/1000,2) > 9.";
	TCut lowQ2 = !highQ2;

	TCut highq2Cut = "TMath::Power(J_psi_1S_MM/1000,2) > 15 && TMath::Power(J_psi_1S_MM/1000,2) < 20";
	TCut lowq2Cut = "TMath::Power(J_psi_1S_MM/1000,2) > 1.1 && TMath::Power(J_psi_1S_MM/1000,2) < 6";
	TCut lowestq2Cut = "TMath::Power(J_psi_1S_MM/1000,2) > 1.1 && TMath::Power(J_psi_1S_MM/1000,2) < 3";
	TCut highestq2Cut = "TMath::Power(J_psi_1S_MM/1000,2) > 18 && TMath::Power(J_psi_1S_MM/1000,2) < 20";

	TCut geomCut = "TMath::Abs(TMath::ACos(muplus_TRUEP_Z / TMath::Sqrt(TMath::Power(muplus_TRUEP_Z,2) + TMath::Power(muplus_TRUEP_Y,2) + TMath::Power(muplus_TRUEP_X,2)))) > 0.01 && TMath::Abs(TMath::ACos(muplus_TRUEP_Z / TMath::Sqrt(TMath::Power(muplus_TRUEP_Z,2) + TMath::Power(muplus_TRUEP_Y,2) + TMath::Power(muplus_TRUEP_X,2)))) < 0.4 && TMath::Abs(TMath::ACos(muminus_TRUEP_Z / TMath::Sqrt(TMath::Power(muminus_TRUEP_Z,2) + TMath::Power(muminus_TRUEP_Y,2) + TMath::Power(muminus_TRUEP_X,2)))) > 0.01 && TMath::Abs(TMath::ACos(muminus_TRUEP_Z / TMath::Sqrt(TMath::Power(muminus_TRUEP_Z,2) + TMath::Power(muminus_TRUEP_Y,2) + TMath::Power(muminus_TRUEP_X,2)))) < 0.4";

	TCut massCutUnblinded = "Lb_MassConsLambda_M[0] > 5200 && Lb_MassConsLambda_M[0] < 6300";
	TCut Lmasscut = "Lambda0_MM > 1105 && Lambda0_MM < 1125";

	TCut L0Passed = "(Lb_L0MuonDecision_TOS || Lb_L0DiMuonDecision_TOS)";
	TCut Hlt1Passed = "(Lb_Hlt1TrackAllL0Decision_TOS || Lb_Hlt1TrackMuonDecision_TOS || Lb_Hlt1DiMuonHighMassDecision_TOS)";
	TCut Hlt2Passed = "(Lb_Hlt2Topo2BodyBBDTDecision_TOS || Lb_Hlt2Topo3BodyBBDTDecision_TOS || Lb_Hlt2Topo4BodyBBDTDecision_TOS || Lb_Hlt2TopoMu2BodyBBDTDecision_TOS || Lb_Hlt2TopoMu3BodyBBDTDecision_TOS || Lb_Hlt2TopoMu4BodyBBDTDecision_TOS || Lb_Hlt2SingleMuonDecision_TOS || Lb_Hlt2DiMuonDetachedDecision_TOS)";
	TCut TrigPassed = L0Passed + Hlt1Passed + Hlt2Passed;

	TCut avoidJpsiCut = "!( TMath::Abs(J_psi_1S_MM - 3096) < 100 || TMath::Abs(J_psi_1S_MM - 3686.1) < 90 )";
	TCut avoidJpsiCut_large = "!( ( J_psi_1S_MM < 3196 && J_psi_1S_MM > 2900 ) || TMath::Abs(J_psi_1S_MM - 3686.1) < 90 )";
	TCut jpsiCut = "TMath::Abs(J_psi_1S_MM - 3096.916) < 92.9";
	
	//TCut MVAcut = "( (weight > 0.732 && TMath::Power(J_psi_1S_MM/1000,2) > 9.590888711) || (weight > 0.969 && TMath::Power(J_psi_1S_MM/1000,2) <= 9.590888711) )";
	TCut MVAcut = "( (weight > 0.66 && TMath::Power(J_psi_1S_MM/1000,2) > 9) || (weight > 0.96 && TMath::Power(J_psi_1S_MM/1000,2) <= 9) )";
	// TCut MVAcut = "( (weight > 0.60 && TMath::Power(J_psi_1S_MM/1000,2) > 9.590888711) || (weight > 0.90 && TMath::Power(J_psi_1S_MM/1000,2) <= 9.590888711) )";
	TCut MVAcut_lowSel = "( (weight > 0.66 && TMath::Power(J_psi_1S_MM/1000,2) > 11) || (weight > 0.96 && TMath::Power(J_psi_1S_MM/1000,2) <= 11) )";

	TCut mumuTrueID = "( (TMath::Abs(Lambda0_MC_MOTHER_ID) == 5122 && TMath::Abs(piminus_MC_MOTHER_ID) == 3122 && TMath::Abs(pplus_MC_MOTHER_ID) == 3122 && TMath::Abs(muminus_MC_MOTHER_ID) == 5122 && TMath::Abs(muplus_MC_MOTHER_ID) == 5122) && (TMath::Abs(Lb_TRUEID) == 5122 && TMath::Abs(Lambda0_TRUEID) == 3122 && TMath::Abs(piminus_TRUEID) == 211 && TMath::Abs(pplus_TRUEID) == 2212 && TMath::Abs(muminus_TRUEID) == 13 && TMath::Abs(muplus_TRUEID) == 13) )";
	TCut jpsiTrueID = "((TMath::Abs(Lambda0_MC_MOTHER_ID) == 5122 && TMath::Abs(piminus_MC_MOTHER_ID) == 3122 && TMath::Abs(pplus_MC_MOTHER_ID) == 3122 && TMath::Abs(muminus_MC_MOTHER_ID) == 443 && TMath::Abs(muplus_MC_MOTHER_ID) == 443 && TMath::Abs(J_psi_1S_MC_MOTHER_ID) == 5122) && (TMath::Abs(Lb_TRUEID) == 5122 && TMath::Abs(Lambda0_TRUEID) == 3122 && TMath::Abs(piminus_TRUEID) == 211 && TMath::Abs(pplus_TRUEID) == 2212 && TMath::Abs(muminus_TRUEID) == 13 && TMath::Abs(muplus_TRUEID) == 13 && TMath::Abs(J_psi_1S_TRUEID) == 443) )";

	TCut DDcut = "pplus_TRACK_Type == 5";
	TCut LLcut = "pplus_TRACK_Type == 3";
	
	TCut mumuSwapID = !mumuTrueID;
	TCut jpsiSwapID = !jpsiTrueID;
	
	TCut upperCut = TrigPassed + MVAcut;
	TCut upperCut_lowSel = TrigPassed + MVAcut_lowSel;

	TCut cutJpsi = upperCut + massCutUnblinded + Lmasscut + jpsiCut;
	TCut cutJpsi_lowSel = upperCut_lowSel + massCutUnblinded + Lmasscut + jpsiCut;
	TCut cutMuMu = upperCut + massCutUnblinded + Lmasscut;
	TCut cutMuMu_notrig = MVAcut + massCutUnblinded + Lmasscut;
	TCut cutMuMu_veto = upperCut + massCutUnblinded + Lmasscut + avoidJpsiCut;
	TCut cutMuMu_noMVA = TrigPassed + massCutUnblinded + Lmasscut;

	TCut baseCutJpsi = massCutUnblinded + Lmasscut + jpsiTrueID + jpsiCut;
	TCut baseCutMuMu = massCutUnblinded + Lmasscut + mumuTrueID;
}

#endif
