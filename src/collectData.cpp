#include "TFile.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TList.h"
#include "TObjArray.h"
#include "TGraph.h"
#include "TROOT.h"
#include <iostream>

using namespace std;

int main()
{
	TFile * file = new TFile("out.root","update");
	TH2F * grid =  (TH2F *)file->Get("pvalues_grid");

	TCanvas * c = new TCanvas();
	double cont_vals[1] = { 0.32 };
	grid->SetContour(1, cont_vals);
	grid->SetMinimum(0);
	grid->SetMaximum(1);
	grid->SetLineWidth(1.5);

	grid->Draw("CONT LIST");

	c->Update();

	TObjArray *conts = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
	TList *lcontour1 = (TList*)conts->At(0);
	TGraph *gc1 = (TGraph*)lcontour1->First();

	double f0 = 0, a0 = 0;
	for(int i = 0; i < gc1->GetN(); i++)
	{
		double a, f;
		gc1->GetPoint(i,a,f);
		if(i==0) { a0=a; f0=f; } 
		if( (f-1)*3./4. > a )
		{
			double m = 4./3.;
			double q = 1.;
			gc1->SetPoint(i,(m*a+f-q)/(2.*m),(m*a+f+q)/2.);
		}
		else if( a > -(f-1)*3./4. )
		{
			double m = -4./3.;
			double q = 1.;
			gc1->SetPoint(i,(m*a+f-q)/(2.*m),(m*a+f+q)/2.);
		}
	}
	gc1->SetPoint(gc1->GetN(),a0,f0);

	gc1->Write("corrected_grid");
	file->Close();
}

