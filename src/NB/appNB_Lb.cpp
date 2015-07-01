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

TTree * appWeight(TreeReader * reader, vector<float> MCnorm = vector<float>(), string tname = "",  bool physW = false, bool mva = true, bool Lnodecay = false);

TTree * appWeight_LbAndMVA(TreeReader * reader, vector<float> MCnorm = vector<float>(), string tname = "",  bool physW = false, bool mva = true, bool Lnodecay = false);

TH2F * weightHist = NULL;
TH2F * PTweightHist = NULL;
TH1 * vtxeffhist = NULL;
TH1F * pidWhist = NULL;
bool isMCGeom = false;

int main(int argc, char** argv)
{
	// Arguments:
	//  1 - MC
	//  2 - output file name
	//  3 - list of input files

	string dotree = "all";
	bool origTree = false;
	bool isMC = false;
	if( argc > 3 )
	{
		if(((string)argv[1]) == "MC") isMC = true;
		else if (((string)argv[1]) == "MCGeom") isMCGeom = true;

		if(((string)argv[2]) == "orig") origTree = true;
		dotree = argv[3];
	}
	if( argc < 6 ) cout << "Missing arguments, please check source for details\n";


	TString treeName = "Lb2Lmumu_Tuple/DecayTree";
	TString treeNameMCgeom = "Lb1_Tuple/MCDecayTree";
	TString treeNameMCgeomDecay = "Lb1_Tuple/MCDecayTree";
	TString treeNameMC = "Lb2Lmumu_MCTuple/MCDecayTree";
	TString treeNameMCdecay = "Lb2Lppimumu_MCTuple/MCDecayTree";

	if(!origTree)
	{
		treeName = "tree";
		treeNameMCgeom = "MCtree";
		treeNameMCgeomDecay = "MCtreeDecay";
		treeNameMC = "MCtree";
		treeNameMCdecay = "MCtreeDecay";
	}


	TreeReader * reader = new TreeReader(treeName);
	TreeReader * MCreaderDecay = NULL, * MCreader = NULL;
	if(isMCGeom)
	{
		MCreader = new TreeReader(treeNameMCgeom);
		MCreaderDecay = new TreeReader(treeNameMCgeomDecay);
	}
	else
	{
		MCreaderDecay = new TreeReader(treeNameMCdecay);
		MCreader = new TreeReader(treeNameMC);
	}

	TList * list = new TList;
	double NGen11 = 0, NGen12 = 0, NGen_PreJune12 = 0;

	for(int a = 5; a < argc; a++)
	{
		if(!isMCGeom) reader->AddFile(argv[a]);
		cout << "Reading: " << argv[a] << endl;

		TFile * f = TFile::Open(argv[a]);
		TTree * tevt = 0;

		if(isMC || isMCGeom || ((string)argv[a]).find("/B")!=string::npos)
		{
			if(((string)argv[a]).find("/B")==string::npos)
			{
				MCreader->AddFile(argv[a]);
				MCreaderDecay->AddFile(argv[a]);
			}

			if(origTree)
			{
				if(isMC || ((string)argv[a]).find("/B")!=string::npos) f->GetObject("EventTuple/EventTuple",tevt);
				else f->GetObject("Lb1_Tuple/MCDecayTree",tevt);
				double nentries = tevt->GetEntries();
				if(((string)argv[a]).find("11")!=string::npos) NGen11 += nentries;
				else if(((string)argv[a]).find("early2012")!=string::npos) NGen_PreJune12+=nentries;
				else NGen12 += nentries;
			}
			else f->GetObject("EventTuple",tevt);
		}
		else
		{
			if(origTree) f->GetObject("GetIntegratedLuminosity/LumiTuple",tevt);
			else f->GetObject("LumiTuple",tevt);
		}

		if(tevt != 0 && !isMCGeom)
		{
			cout << "Adding " << tevt->GetName() << endl;
			list->Add(tevt);
		}
		cout << "Done" << endl;
	}

	vector<float> MCnorm;
	if(!origTree && (isMC || isMCGeom) ) MCnorm.push_back(0);
	else 
	{
		if(isMC)
		{
			MCnorm.push_back(1.e6 * 0.944 / NGen11); // 2011
			MCnorm.push_back(1.e6 * 0.517 / NGen_PreJune12); // 2012 PreJune
			MCnorm.push_back(1.e6 * 1.442 / (2 * NGen12) ); // 2012 PostJone
			cout << MCnorm[0] << "   " << MCnorm[1] << "    " << MCnorm[2] << endl;
			cout << NGen11 << "   " << NGen_PreJune12 << "    " << NGen12 << endl;
		}
		else if(isMCGeom)
		{
			MCnorm.push_back(1.e6 * 0.944 / NGen11); // 2011
			MCnorm.push_back(0.); // 2012 PreJune
			MCnorm.push_back(1.e6 * 1.959 / NGen12); // 2012 PostJone	
		}
	}

	if(!isMCGeom) reader->Initialize();
	vector<string > ex;
	ex.push_back("DLL");
	if(isMC || isMCGeom) MCreader->Initialize(ex,"exceptcontains");
	if(isMC || isMCGeom) MCreaderDecay->Initialize(ex,"exceptcontains");


	TFile * VtxEffFile = TFile::Open("/afs/cern.ch/user/p/pluca/work/Lb/Lmumu/vertex_eff_ratio_DD.root");
	TH1 * tmpvtxeff = (TH1 *)VtxEffFile->Get("vertex_eff_ratio");
	vtxeffhist = (TH1*)tmpvtxeff->Clone();

	TFile * weightFile = TFile::Open("/afs/cern.ch/user/p/pluca/work/Lb/Lmumu/MC_reweight.root");
	weightHist = (TH2F *)weightFile->Get("ratio2D");
	PTweightHist = (TH2F *)weightFile->Get("ratio2D_pt");
	pidWhist = (TH1F *)weightFile->Get("ratio_PIDmu");

	TFile out(argv[4],"recreate");

	cout << "Weighting " << dotree << endl;

	if(!isMCGeom)
	{
		TTree * lumitree = TTree::MergeTrees(list);
		cout << "NGen = " << lumitree->GetEntries() << endl;
		lumitree->Write();

		TTree * newTree = NULL;

		if(dotree=="all" || dotree=="tree")
		{
			if(origTree)
			{
				if(isMC) newTree = appWeight(reader,MCnorm,"tree",true);
				else newTree = appWeight(reader);
			}
			else
			{
				if(isMC) newTree = appWeight_LbAndMVA(reader,MCnorm,"tree",true);
				else newTree = appWeight_LbAndMVA(reader);
			}
			newTree->Write();
		}
	}
	
	if(isMC || isMCGeom)
	{
		if(origTree)
		{
			if(dotree=="all" || dotree=="MCtree") 
			{ TTree * newMCTree = appWeight(MCreader,MCnorm,"MCtree",true,false,true); newMCTree->Write(); }
			if(dotree=="all" || dotree=="MCDecaytree") 
			{ TTree * newMCTreeDecay = appWeight(MCreaderDecay,MCnorm,"MCtreeDecay",true,false); newMCTreeDecay->Write(); }
		}
		else
		{
			if(dotree=="all" || dotree=="MCtree")
			{ TTree * newMCTree = appWeight_LbAndMVA(MCreader,MCnorm,"MCtree",true,false,true); newMCTree->Write(); }
			if(dotree=="all" || dotree=="MCDecaytree") 
			{ TTree * newMCTreeDecay = appWeight_LbAndMVA(MCreaderDecay,MCnorm,"MCtreeDecay",true,false); newMCTreeDecay->Write(); }
		}
	}

	out.Close();
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


TTree * appWeight_LbAndMVA(TreeReader * reader, vector<float> MCnorm, string tname, bool physW, bool mva, bool Lnodecay)
{
	TString nametree = "tree";
	if(tname != "") nametree = tname;
	TTree* newTree = new TTree(nametree,"");
	reader->BranchNewTree(newTree);

	float NNOut1;
	float Lb_weight, pt_weight, pidWeight_plus, pidWeight_minus, pplus_eta;
	float physRate_pol0_wilson1_noDecay, physRate_pol0_wilson2_noDecay, physRate_pol0_wilson3_noDecay, physRate_pol0_wilson1, physRate_pol0_wilson2, physRate_pol0_wilson3;

	string type = "mumu";

	if(mva)
	{
		newTree->Branch("weight",&NNOut1,"weight/F");
		if(MCnorm.size() > 0) newTree->Branch("pid_muplus",&pidWeight_plus,"pid_muplus/F");
		if(MCnorm.size() > 0) newTree->Branch("pid_muminus",&pidWeight_minus,"pid_muminus/F");
	}
	if(MCnorm.size() > 0 || isMCGeom)
	{
		if(mva) newTree->Branch("pplus_TRACK_Eta",&pplus_eta,"pplus_TRACK_Eta/F");
		newTree->Branch("Lb_weight",&Lb_weight,"Lb_weight/F");
		newTree->Branch("pt_weight",&pt_weight,"pt_weight/F");

		if(((string)reader->GetChain()->GetFile()->GetName()).find("Jpsi")!=string::npos) type = "jpsi";
		if(physW && type == "mumu")
		{		
			newTree->Branch("physRate_pol0_wilson1_noDecay",&physRate_pol0_wilson1_noDecay,"physRate_pol0_wilson1_noDecay/F");
			newTree->Branch("physRate_pol0_wilson2_noDecay",&physRate_pol0_wilson2_noDecay,"physRate_pol0_wilson2_noDecay/F");
			newTree->Branch("physRate_pol0_wilson3_noDecay",&physRate_pol0_wilson3_noDecay,"physRate_pol0_wilson3_noDecay/F");

			if(!Lnodecay)
			{
				newTree->Branch("physRate_pol0_wilson1",&physRate_pol0_wilson1,"physRate_pol0_wilson1/F");
				newTree->Branch("physRate_pol0_wilson2",&physRate_pol0_wilson2,"physRate_pol0_wilson2/F");
				newTree->Branch("physRate_pol0_wilson3",&physRate_pol0_wilson3,"physRate_pol0_wilson3/F");		
			}
		}
	}

	LbLmumuAng model_wilson1(0.06,1), model_wilson2(0.06,2), model_wilson3(0.06,3);

	string expertfile = "expert_12.nb";
	Expert* nbExpert = new Expert(expertfile.c_str());

	int ntot = reader->GetEntries();
	cout << tname << "  " << ntot << endl;

	for( int i = 0; i < ntot; ++i )
	{
		showPercentage(i,ntot);
		reader->GetEntry(i);

		if(MCnorm.size() > 0 || isMCGeom)
		{
			//Lb momentum reweight

			float Lbpz = reader->GetValue("Lb_TRUEP_Z");
			float Lbpt = reader->GetValue("Lb_TRUEPT");
			float Lbp = TMath::Sqrt( TMath::Power(Lbpz,2) + TMath::Power(Lbpt,2) );
			Lb_weight = weightHist->GetBinContent(weightHist->GetXaxis()->FindBin(Lbpt),weightHist->GetYaxis()->FindBin(Lbp));

			float Lambda0pt = reader->GetValue("Lambda0_TRUEPT");
			pt_weight = PTweightHist->GetBinContent(PTweightHist->GetXaxis()->FindBin(Lbpt),PTweightHist->GetYaxis()->FindBin(Lambda0pt));

			float p_pz = reader->GetValue("pplus_TRUEP_Z");
			float p_pt = reader->GetValue("pplus_TRUEPT");
			float p_p = TMath::Sqrt( TMath::Power(p_pz,2) + TMath::Power(p_pt,2) );
			pplus_eta = 0.5*TMath::Log((p_p + p_pz)/(p_p - p_pz));

			if(type=="mumu")
			{
				bool is11 = (reader->GetValue("dataType") == 2011);

				float perg = 4.e6;
				if(is11) perg = 3.5e6;
				float pmom = TMath::Sqrt(perg*perg-938.3*938.3);
				TLorentzVector initialProton(0.,0.,pmom,perg);
				TLorentzVector pion, proton;
				TLorentzVector Lb = get4momentum(reader,"Lb");
				if(!Lnodecay)
				{
					pion = get4momentum(reader,"piminus");
					proton = get4momentum(reader,"pplus");
				}
				else
				{
					pion.SetPxPyPzE(0,0,0,0);
					proton.SetPxPyPzE(0,0,0,0);
				}
				TLorentzVector Lambda = get4momentum(reader,"Lambda0");
				TLorentzVector mup = get4momentum(reader,"muplus");
				TLorentzVector mum = get4momentum(reader,"muminus");
				TLorentzVector Jpsi = mum+mup;

				int pcharge = 1;
				if(reader->GetValue("Lambda0_MC_MOTHER_ID") < 0) pcharge = -1;


				float q2 = TMath::Power(Jpsi.M()/1000.,2);

				if(physW)
				{
					double cosTheta,cosThetaL,cosThetaB,phiL,phiB,dphi;
					DPHelpers::LbPsiRAngles(initialProton,Lb,Jpsi,Lambda,mup,mum,proton,pion,pcharge,
							cosTheta,cosThetaL,cosThetaB,phiL,phiB,dphi);

					if(!Lnodecay)
					{
						physRate_pol0_wilson1 = model_wilson1.physicsRate(q2, cosTheta, cosThetaL, cosThetaB, phiL, phiB);
						if(isnan(physRate_pol0_wilson1)) { cout << "NAN   " << endl; continue; } 
						physRate_pol0_wilson2 = model_wilson2.physicsRate(q2, cosTheta, cosThetaL, cosThetaB, phiL, phiB);
						physRate_pol0_wilson3 = model_wilson3.physicsRate(q2, cosTheta, cosThetaL, cosThetaB, phiL, phiB);
					}

					physRate_pol0_wilson1_noDecay = model_wilson1.physicsRate(q2, cosTheta, cosThetaL, phiL);
					if(isnan(physRate_pol0_wilson1_noDecay)) { cout << "noDecay NAN" << endl; continue; } 
					physRate_pol0_wilson2_noDecay = model_wilson2.physicsRate(q2, cosTheta, cosThetaL, phiL);
					physRate_pol0_wilson3_noDecay = model_wilson3.physicsRate(q2, cosTheta, cosThetaL, phiL);
				}
			}

		}

		if(mva)
		{
			if(MCnorm.size() > 0)
			{
				float pidplus = reader->GetValue("muplus_PIDmu");
				float pidminus = reader->GetValue("muminus_PIDmu");
				pidWeight_plus = pidWhist->GetBinContent(pidWhist->GetXaxis()->FindBin(pidplus));
				pidWeight_minus = pidWhist->GetBinContent(pidWhist->GetXaxis()->FindBin(pidminus));
			}

			float inputs[100];
			fillInputArray(reader,inputs);
			if(isnan(inputs[0])) continue;

			NNOut1 = nbExpert->nb_expert(inputs);
			NNOut1 = 0.5*(NNOut1+1);
		}

		newTree->Fill();
	}

	return newTree;
}




TTree * appWeight(TreeReader * reader, vector<float> MCnorm, string tname, bool physW, bool mva, bool Lnodecay)
{
	TString nametree = "tree";
	if(tname != "") nametree = tname;
	TTree* newTree = new TTree(nametree,"");
	reader->BranchNewTree(newTree);

	reader->GetEntry(0);
	string type = "mumu";
	if(((string)reader->GetChain()->GetFile()->GetName()).find("Jpsi")!=string::npos) type = "jpsi";

	float NNOut1;
	int year = 2012;
	float myMCnorm = 1., pplus_eta;
	float myq2, DDvtx_weight,pidWeight_plus, pidWeight_minus;;
	float lfTimeW, lfTimeW_psigma, lfTimeW_msigma;
	float physRate_pol0, physRate_polm003, physRate_polp015, physRate_polp006, physRate_pol0_QCDff;
	float modelRate_jpsi_1, modelRate_jpsi_2, modelRate_jpsi_3, modelRate_jpsi_4, modelRate_jpsi_5, modelRate_jpsi_6, modelRate_jpsi_7, modelRate_jpsi_8;
	float physRate_pol0_noDecay, physRate_polm003_noDecay, physRate_polp015_noDecay, physRate_polp006_noDecay, physRate_pol0_QCDff_noDecay;
	float modelRate_jpsi_1_noDecay, modelRate_jpsi_2_noDecay, modelRate_jpsi_3_noDecay, modelRate_jpsi_4_noDecay, modelRate_jpsi_5_noDecay, modelRate_jpsi_6_noDecay, modelRate_jpsi_7_noDecay, modelRate_jpsi_8_noDecay;
	float physRate_pol0_wilson1_noDecay, physRate_pol0_wilson2_noDecay, physRate_pol0_wilson3_noDecay, physRate_pol0_wilson1, physRate_pol0_wilson2, physRate_pol0_wilson3;
	double cosTheta, cosThetaL, cosThetaB, phiL, phiB, dphi;
	float Lb_weight, pt_weight;

	newTree->Branch("dataType",&year,"dataType/I");
	if(mva)
	{
		newTree->Branch("weight",&NNOut1,"weight/F");
		if(MCnorm.size() > 0) newTree->Branch("DDvtx_weight",&DDvtx_weight,"DDvtx_weight/F");
		if(MCnorm.size() > 0) newTree->Branch("pid_muplus",&pidWeight_plus,"pid_muplus/F");
		if(MCnorm.size() > 0) newTree->Branch("pid_muminus",&pidWeight_minus,"pid_muminus/F");
	}
	if(MCnorm.size() > 0)
	{
		if(!reader->HasVar("J_psi_1S_MM")) newTree->Branch("J_psi_1S_MM",&myq2,"J_psi_1S_MM/F");
		if(mva) newTree->Branch("pplus_TRACK_Eta",&pplus_eta,"pplus_TRACK_Eta/F");
		newTree->Branch("MCnorm",&myMCnorm,"MCnorm/F");
		newTree->Branch("Lb_weight",&Lb_weight,"Lb_weight/F");
		newTree->Branch("pt_weight",&pt_weight,"pt_weight/F");
		newTree->Branch("lifeTimeW",&lfTimeW,"lifeTimeW/F");
		newTree->Branch("lifeTimeW_plussigma",&lfTimeW_psigma,"lifeTimeW_plussigma/F");
		newTree->Branch("lifeTimeW_minussigma",&lfTimeW_msigma,"lifeTimeW_minussigma/F");

		if(physW)
		{
			newTree->Branch("cosTheta",&cosTheta,"cosTheta/D");
			newTree->Branch("cosThetaL",&cosThetaL,"cosThetaL/D");
			newTree->Branch("cosThetaB",&cosThetaB,"cosThetaB/D");
			newTree->Branch("phiL",&phiL,"phiL/D");
			newTree->Branch("phiB",&phiB,"phiB/D");
			newTree->Branch("dphi",&dphi,"dphi/D");

			newTree->Branch("physRate_pol0_noDecay",&physRate_pol0_noDecay,"physRate_pol0_noDecay/F");

			if(type == "mumu")
			{
				newTree->Branch("physRate_polp015_noDecay",&physRate_polp015_noDecay,"physRate_polp015_noDecay/F");
				newTree->Branch("physRate_polm003_noDecay",&physRate_polm003_noDecay,"physRate_polm003_noDecay/F");
				newTree->Branch("physRate_polp006_noDecay",&physRate_polp006_noDecay,"physRate_polp006_noDecay/F");
				newTree->Branch("physRate_pol0_QCDff_noDecay",&physRate_pol0_QCDff_noDecay,"physRate_pol0_QCDff_noDecay/F");
				newTree->Branch("physRate_pol0_wilson1_noDecay",&physRate_pol0_wilson1_noDecay,"physRate_pol0_wilson1_noDecay/F");
				newTree->Branch("physRate_pol0_wilson2_noDecay",&physRate_pol0_wilson2_noDecay,"physRate_pol0_wilson2_noDecay/F");
				newTree->Branch("physRate_pol0_wilson3_noDecay",&physRate_pol0_wilson3_noDecay,"physRate_pol0_wilson3_noDecay/F");
			}
			else
			{
				newTree->Branch("model_jpsi_1_noDecay",&modelRate_jpsi_1_noDecay,"model_jpsi_1_noDecay/F");
				newTree->Branch("model_jpsi_2_noDecay",&modelRate_jpsi_2_noDecay,"model_jpsi_2_noDecay/F");
				newTree->Branch("model_jpsi_3_noDecay",&modelRate_jpsi_3_noDecay,"model_jpsi_3_noDecay/F");
				newTree->Branch("model_jpsi_4_noDecay",&modelRate_jpsi_4_noDecay,"model_jpsi_4_noDecay/F");
				newTree->Branch("model_jpsi_5_noDecay",&modelRate_jpsi_5_noDecay,"model_jpsi_5_noDecay/F");
				newTree->Branch("model_jpsi_6_noDecay",&modelRate_jpsi_6_noDecay,"model_jpsi_6_noDecay/F");
				newTree->Branch("model_jpsi_7_noDecay",&modelRate_jpsi_7_noDecay,"model_jpsi_7_noDecay/F");
				newTree->Branch("model_jpsi_8_noDecay",&modelRate_jpsi_8_noDecay,"model_jpsi_8_noDecay/F");
			}


			if(!Lnodecay)
			{
				newTree->Branch("physRate_pol0",&physRate_pol0,"physRate_pol0/F");

				if(type == "mumu")
				{
					newTree->Branch("physRate_polp015",&physRate_polp015,"physRate_polp015/F");
					newTree->Branch("physRate_polm003",&physRate_polm003,"physRate_polm003/F");
					newTree->Branch("physRate_polp006",&physRate_polp006,"physRate_polp006/F");
					newTree->Branch("physRate_pol0_QCDff",&physRate_pol0_QCDff,"physRate_pol0_QCDff/F");
					newTree->Branch("physRate_pol0_wilson1",&physRate_pol0_wilson1,"physRate_pol0_wilson1/F");
					newTree->Branch("physRate_pol0_wilson2",&physRate_pol0_wilson2,"physRate_pol0_wilson2/F");
					newTree->Branch("physRate_pol0_wilson3",&physRate_pol0_wilson3,"physRate_pol0_wilson3/F");
				}
				else
				{
					newTree->Branch("model_jpsi_1",&modelRate_jpsi_1,"model_jpsi_1/F");
					newTree->Branch("model_jpsi_2",&modelRate_jpsi_2,"model_jpsi_2/F");
					newTree->Branch("model_jpsi_3",&modelRate_jpsi_3,"model_jpsi_3/F");
					newTree->Branch("model_jpsi_4",&modelRate_jpsi_4,"model_jpsi_4/F");
					newTree->Branch("model_jpsi_5",&modelRate_jpsi_5,"model_jpsi_5/F");
					newTree->Branch("model_jpsi_6",&modelRate_jpsi_6,"model_jpsi_6/F");
					newTree->Branch("model_jpsi_7",&modelRate_jpsi_7,"model_jpsi_7/F");
					newTree->Branch("model_jpsi_8",&modelRate_jpsi_8,"model_jpsi_8/F");
				}
			}
		}
	}

	double genLife = 0.001424702;
	double pdgLife = 0.001482;
	double lfsigma = 0.000021;

	LbLmumuAng model_pol0(0.), model_polp015(0.15), model_polm003(-0.03), model_polp006(0.06), model_pol0_QCDff(0.);
	LbLmumuAng model_wilson1(0.06,1), model_wilson2(0.06,2), model_wilson3(0.06,3);
	LbJpsiL2 model_jpsi_0, model_jpsi_1, model_jpsi_2, model_jpsi_3, model_jpsi_4, model_jpsi_5, model_jpsi_6, model_jpsi_7, model_jpsi_8;
	model_pol0_QCDff.setFFtype(LbLmumuAng::LQCD);
	model_jpsi_1.setVariation(1);
	model_jpsi_2.setVariation(2);
	model_jpsi_3.setVariation(3);
	model_jpsi_4.setVariation(4);
	model_jpsi_5.setVariation(5);
	model_jpsi_6.setVariation(6);
	model_jpsi_7.setVariation(7);
	model_jpsi_8.setVariation(8);

	string expertfile = "expert_12.nb";
	Expert* nbExpert = new Expert(expertfile.c_str());

	int ntot = reader->GetEntries();
	cout << tname << "  " << ntot << endl;
	//ntot = 1000;
	for( int i = 0; i < ntot; ++i )
	{
		showPercentage(i,ntot);
		reader->GetEntry(i);
		string filename = (string)reader->GetChain()->GetFile()->GetName();
		bool is11 = (filename.find("11") != string::npos);

		if(MCnorm.size() > 0 || isMCGeom)
		{
			//Calc Lb P/PT weight

			float Lbpx = reader->GetValue("Lb_TRUEP_X");
			float Lbpy = reader->GetValue("Lb_TRUEP_Y");
			float Lbpz = reader->GetValue("Lb_TRUEP_Z");
			float Lbp = TMath::Sqrt(TMath::Power(Lbpx,2) + TMath::Power(Lbpy,2) + TMath::Power(Lbpz,2));
			float Lbpt = reader->GetValue("Lb_TRUEPT");
			Lb_weight = weightHist->GetBinContent(weightHist->GetXaxis()->FindBin(Lbpt),weightHist->GetYaxis()->FindBin(Lbp));

			float Lambda0pt = reader->GetValue("Lambda0_TRUEPT");
			pt_weight = PTweightHist->GetBinContent(PTweightHist->GetXaxis()->FindBin(Lbpt),PTweightHist->GetYaxis()->FindBin(Lambda0pt));

			//Calc lifetime weight

			double decTau = reader->GetValue("Lb_TRUETAU");
			if( decTau >=0 )
			{
				lfTimeW = TMath::Exp(-decTau/pdgLife)/TMath::Exp(-decTau/genLife); 
				lfTimeW_psigma = TMath::Exp(-decTau/(pdgLife+lfsigma))/TMath::Exp(-decTau/genLife);
				lfTimeW_msigma = TMath::Exp(-decTau/(pdgLife-lfsigma))/TMath::Exp(-decTau/genLife);
			}
			else { lfTimeW = lfTimeW_psigma = lfTimeW_msigma = -1; }


			//Calc physics weight

			float perg = 4.e6;
			if(is11) perg = 3.5e6;
			float pmom = TMath::Sqrt(perg*perg-938.3*938.3);
			TLorentzVector initialProton(0.,0.,pmom,perg);
			TLorentzVector Jpsi, pion, proton;
			TLorentzVector Lb = get4momentum(reader,"Lb");
			if(!Lnodecay)
			{
				pion = get4momentum(reader,"piminus");
				proton = get4momentum(reader,"pplus");
			}
			else
			{
				pion.SetPxPyPzE(0,0,0,0);
				proton.SetPxPyPzE(0,0,0,0);
			}
			TLorentzVector Lambda = get4momentum(reader,"Lambda0");
			TLorentzVector mup = get4momentum(reader,"muplus");
			TLorentzVector mum = get4momentum(reader,"muminus");

			int pcharge = 1;
			if(reader->GetValue("Lambda0_MC_MOTHER_ID") < 0) pcharge = -1;

			if(type == "jpsi")
			{
				if(isMCGeom) Jpsi = get4momentum(reader,"Jpsi");
				else Jpsi = get4momentum(reader,"J_psi_1S");

				if(physW)
				{
					LbJpsiL2::LbPsiRAngles(initialProton, Lb, Jpsi, Lambda, mup, mum, proton, pion, pcharge,
							cosTheta, cosThetaL, cosThetaB, phiL, phiB, dphi);

					if(!Lnodecay)
					{
						physRate_pol0 = model_jpsi_0.physicsRate(cosTheta, cosThetaL, cosThetaB, phiL, phiB);
						if(isnan(physRate_pol0)) { cout << "\nNAN:   " <<  cosTheta << "   " <<  cosThetaL << "   "  <<  phiL << endl; continue;}
						modelRate_jpsi_1 = model_jpsi_1.physicsRate(cosTheta, cosThetaL, cosThetaB, phiL, phiB);
						modelRate_jpsi_2 = model_jpsi_2.physicsRate(cosTheta, cosThetaL, cosThetaB, phiL, phiB);
						modelRate_jpsi_3 = model_jpsi_3.physicsRate(cosTheta, cosThetaL, cosThetaB, phiL, phiB);
						modelRate_jpsi_4 = model_jpsi_4.physicsRate(cosTheta, cosThetaL, cosThetaB, phiL, phiB);
						modelRate_jpsi_5 = model_jpsi_5.physicsRate(cosTheta, cosThetaL, cosThetaB, phiL, phiB);
						modelRate_jpsi_6 = model_jpsi_6.physicsRate(cosTheta, cosThetaL, cosThetaB, phiL, phiB);
						modelRate_jpsi_7 = model_jpsi_7.physicsRate(cosTheta, cosThetaL, cosThetaB, phiL, phiB);
						modelRate_jpsi_8 = model_jpsi_8.physicsRate(cosTheta, cosThetaL, cosThetaB, phiL, phiB);
					}

					physRate_pol0_noDecay = model_jpsi_0.physicsRate(cosTheta, cosThetaL, phiL);
					if(isnan(physRate_pol0_noDecay)) { cout << "\nnoDecay NAN:   " <<  cosTheta << "   " <<  cosThetaL << "   "  <<  phiL << endl; continue;}
					modelRate_jpsi_1_noDecay = model_jpsi_1.physicsRate(cosTheta, cosThetaL, phiL);
					modelRate_jpsi_2_noDecay = model_jpsi_2.physicsRate(cosTheta, cosThetaL, phiL);
					modelRate_jpsi_3_noDecay = model_jpsi_3.physicsRate(cosTheta, cosThetaL, phiL);
					modelRate_jpsi_4_noDecay = model_jpsi_4.physicsRate(cosTheta, cosThetaL, phiL);
					modelRate_jpsi_5_noDecay = model_jpsi_5.physicsRate(cosTheta, cosThetaL, phiL);
					modelRate_jpsi_6_noDecay = model_jpsi_6.physicsRate(cosTheta, cosThetaL, phiL);
					modelRate_jpsi_7_noDecay = model_jpsi_7.physicsRate(cosTheta, cosThetaL, phiL);
					modelRate_jpsi_8_noDecay = model_jpsi_8.physicsRate(cosTheta, cosThetaL, phiL);
				}
			}
			else
			{
				Jpsi = mup+mum;
				float q2 = TMath::Power(Jpsi.M()/1000.,2);

				if(physW)
				{
					DPHelpers::LbPsiRAngles(initialProton,Lb,Jpsi,Lambda,mup,mum,proton,pion,pcharge,
							cosTheta,cosThetaL,cosThetaB,phiL,phiB,dphi);

					if(!Lnodecay)
					{
						physRate_pol0 = model_pol0.physicsRate(q2, cosTheta, cosThetaL, cosThetaB, phiL, phiB);
						if(isnan(physRate_pol0)) { cout << "NAN:   " << q2 << "   " <<  cosTheta << "   " <<  cosThetaL << "   "  <<  phiL << endl; continue;}
						physRate_polm003 = model_polm003.physicsRate(q2, cosTheta, cosThetaL, cosThetaB, phiL, phiB);
						physRate_polp015 = model_polp015.physicsRate(q2, cosTheta, cosThetaL, cosThetaB, phiL, phiB);
						physRate_polp006 = model_polp006.physicsRate(q2, cosTheta, cosThetaL, cosThetaB, phiL, phiB);
						physRate_pol0_QCDff = model_pol0_QCDff.physicsRate(q2, cosTheta, cosThetaL, cosThetaB, phiL, phiB);
						physRate_pol0_wilson1 = model_wilson1.physicsRate(q2, cosTheta, cosThetaL, cosThetaB, phiL, phiB);
						physRate_pol0_wilson2 = model_wilson2.physicsRate(q2, cosTheta, cosThetaL, cosThetaB, phiL, phiB);
						physRate_pol0_wilson3 = model_wilson3.physicsRate(q2, cosTheta, cosThetaL, cosThetaB, phiL, phiB);
					}

					physRate_pol0_noDecay = model_pol0.physicsRate(q2, cosTheta, cosThetaL, phiL);
					if(isnan(physRate_pol0_noDecay)) { cout << "noDecay NAN:   " << q2 << "   " <<  cosTheta << "   " <<  cosThetaL << "   "  <<  phiL << endl; continue;}
					physRate_polm003_noDecay = model_polm003.physicsRate(q2, cosTheta, cosThetaL, phiL);
					physRate_polp015_noDecay = model_polp015.physicsRate(q2, cosTheta, cosThetaL, phiL);
					physRate_polp006_noDecay = model_polp006.physicsRate(q2, cosTheta, cosThetaL, phiL);
					physRate_pol0_QCDff_noDecay = model_pol0_QCDff.physicsRate(q2, cosTheta, cosThetaL, phiL);
					physRate_pol0_wilson1_noDecay = model_wilson1.physicsRate(q2, cosTheta, cosThetaL, phiL);
					physRate_pol0_wilson2_noDecay = model_wilson2.physicsRate(q2, cosTheta, cosThetaL, phiL);
					physRate_pol0_wilson3_noDecay = model_wilson3.physicsRate(q2, cosTheta, cosThetaL, phiL);
				}
			}

			myq2 = Jpsi.M();
		}

		if( is11 )
		{
			year = 2011;
			if(MCnorm.size() > 0) myMCnorm = MCnorm[0];
		}
		else if(filename.find("early2012")!=string::npos && MCnorm.size() > 0) myMCnorm = MCnorm[1]; 
		else if(MCnorm.size() > 0) myMCnorm = MCnorm[2];

		if(mva)
		{
			if(MCnorm.size() > 0)
			{
				float p_pz = reader->GetValue("pplus_TRUEP_Z");
				float p_pt = reader->GetValue("pplus_TRUEPT");
				float p_p = TMath::Sqrt( TMath::Power(p_pz,2) + TMath::Power(p_pt,2) );
				pplus_eta = 0.5*TMath::Log((p_p + p_pz)/(p_p - p_pz));

				DDvtx_weight = 1.;
				if(reader->GetValue("pplus_TRACK_Type") == 5)
				{
					float Lp = reader->GetValue("Lambda0_P");
					if(Lp > 100e3) Lp = 99e3;
					if(Lp < 10e3) Lp = 10e3;
					DDvtx_weight = vtxeffhist->GetBinContent(vtxeffhist->GetXaxis()->FindBin(Lp));
				}
				float pidplus = reader->GetValue("muplus_PIDmu");
				float pidminus = reader->GetValue("muminus_PIDmu");
				pidWeight_plus = pidWhist->GetBinContent(pidWhist->GetXaxis()->FindBin(pidplus));
				pidWeight_minus = pidWhist->GetBinContent(pidWhist->GetXaxis()->FindBin(pidminus));
			}

			float inputs[100];
			fillInputArray(reader,inputs);
			if(isnan(inputs[0])) continue;

			NNOut1 = nbExpert->nb_expert(inputs);
			NNOut1 = 0.5*(NNOut1+1);
		}

		newTree->Fill();
	}

	return newTree;
}





