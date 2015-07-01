#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"

int main()
{
	TF1 * pdf = new TF1("pdf","(3./8.)*(1.-[1])*(1 + TMath::Power(x,2)) + [0]*x + (3./4.)*[1]*(1 - TMath::Power(x,2))",-1.,1.);

	int p = 1;
	double afb = -2;
	double maxafb = 2;
	double maxfL = 2;
	double afbstep = (maxafb - afb)/500.;
	double fL = 0;
	double fLstep = (maxfL - fL)/500.;

	TGraph * gr = new TGraph();

	for(; afb <= maxafb; afb += afbstep)
	{
		for(fL = 0; fL <= maxfL; fL += fLstep)
		{
			pdf->SetParameter(0,afb);
			pdf->SetParameter(1,fL);

			bool isOK = true;
			for(double x = -1; x <= 1; x+=0.001)
				if(pdf->Eval(x) < 0 ) { isOK = false; break; }
			
			if(isOK) { gr->SetPoint(p,afb,fL); p++; }
		}
	}

	TCanvas * c = new TCanvas();
	gr->Draw("AP");
	gr->GetXaxis()->SetTitle("afb");
	gr->GetYaxis()->SetTitle("fL");
	gr->Draw("AP");
	c->Print("scan.pdf");
}
