#include <general_functions.hpp>

int main(int argc, char ** argv)
{
	double lumi = 0;
	for(int i = 1; i < argc; i++) 
	{
		cout << argv[i] << endl;
		lumi += luminosity((TString)argv[i],"");
	}
	cout << "Tot lumi = " << lumi << endl;
}
