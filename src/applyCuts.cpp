#include "ReadTree_comp.hpp"
#include "general_functions.hpp"

using namespace std;


TTree * applyCuts(TString name, TreeReader * dataReader, TCut * _cuts, void (*addFunc)(TreeReader *,  TTree *, bool) = NULL, double frac = -1)
{
	TString skimname("skim");
	if(_cuts) dataReader->GetChain()->Draw(">>" + skimname,*_cuts,"entrylist");
	else dataReader->GetChain()->Draw(">>" + skimname,"","entrylist");
	TEntryList *skim = (TEntryList*)gDirectory->Get(skimname);
	if(skim == NULL) return NULL;
	int nEntries = skim->GetN();
	dataReader->SetEntryList(skim);

	TTree * newTree = new TTree("cand"+name,"");
	dataReader->BranchNewTree(newTree);
	if(addFunc) addFunc(dataReader,newTree,true);
	
	cout << "N candidates = " << nEntries << endl;
	if(frac > 0 && frac < 1) { nEntries *= frac; cout << "Using only " << 100*frac << "% of the entries" << endl; }
	else if (frac > 1) { nEntries = frac; cout << "Using only " << frac << " entries" << endl; }
	
	for(Long64_t i = 0 ; i < nEntries ; i++)
	{
		showPercentage(i,nEntries);

		dataReader->GetEntry(i,skim);
		if(addFunc) addFunc(dataReader,newTree,false);
		
		newTree->Fill();
	}
	
	dataReader->SetEntryList(0);

	delete skim;
	return newTree;
}



int main(int argc, char **argv)
{
	TString nametree, namefile, ofilename;
	TCut cuts;

	if(argc >= 4)
	{
		nametree = argv[1];
		namefile = argv[2];
		cuts = argv[3];

		cout << nametree << "  " << namefile << endl;
		cuts.Print();

		ofilename = namefile;
		ofilename = ofilename.ReplaceAll(".root","")+"_cuts.root";

		for(int a = 1; a < argc; a++)
		{
			string arg = argv[a];
			string str = arg.substr(2,arg.length()-2);

			if(arg.find("-O")!=string::npos) ofilename = str;
		}
	}
	else { cout << "Too few arguments!" << endl; return 1; }

	TreeReader * reader = new TreeReader(nametree);
	reader->AddFile(namefile);
	reader->Initialize();

	TFile * ofile = new TFile(ofilename,"recreate");
	
	TTree * newTree = applyCuts("",reader,&cuts);

	newTree->Write();
	ofile->Close();
}


