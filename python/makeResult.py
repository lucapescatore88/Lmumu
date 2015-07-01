from ROOT import *
from string import *
import sys

gROOT.ProcessLine(".x ~/work/lhcbStyle.C");

printeff = True

from ROOT import TPaveText, gStyle

lhcbName = TPaveText(gStyle.GetPadLeftMargin() + 0.05,
                        0.87 - gStyle.GetPadTopMargin(),
                        gStyle.GetPadLeftMargin() + 0.20,
                        0.95 - gStyle.GetPadTopMargin(),
                        "BRNDC");

lhcbName.AddText("LHCb");
lhcbName.SetFillColor(0);
lhcbName.SetTextAlign(12);
lhcbName.SetBorderSize(0);


def makeResult(base, type = 'rel') :	

	file = TFile(base+"Lb_yield/Lbyield.root")
	central = file.Get("BRplot")
	file = TFile(base+"Lb_yield_norm_plus/Lbyield.root")
	norm_plus = file.Get("BRplot")
	file = TFile(base+"Lb_yield_norm_minus/Lbyield.root")
	norm_minus = file.Get("BRplot")
	file = TFile(base+"Lb_yield_eff_plus/Lbyield.root")
	eff_plus = file.Get("BRplot")
	file = TFile(base+"Lb_yield_eff_minus/Lbyield.root")
	eff_minus = file.Get("BRplot")

	#central.Print()
	result = TGraphAsymmErrors()

	pdfsys = { 17.5 : 0.01, 11.75 : 0.032, 15.5 : 0.028, 17.0 : 0.014, 19.0 : 0.025, 3.55 : 0.042, 1.05 : 0.034, 3.0 : 0.038, 5.0 : 0.066, 7.0 : 0.020}

	Jpsi2mumuBr = 0.0593
	Jpsi2mumuBr_err = 0.0006
	Jpsi2mumuBr_relerr = Jpsi2mumuBr_err/Jpsi2mumuBr


	print "\n\nq2 bin \t\t & Combined yield \\\\"

	for i in range(0,central.GetN()) :
		x, c = Double(0), Double(0)
		np, nm = Double(0), Double(0)
		ep, em = Double(0), Double(0)
		central.GetPoint(i,x,c)
		#if(x < 9) :
		#	continue
	
		errlow = central.GetErrorYlow(i)
		errhigh = central.GetErrorYhigh(i)
		errX = central.GetErrorX(i)
		norm_plus.GetPoint(i,x,np)
		norm_minus.GetPoint(i,x,nm)
		eff_plus.GetPoint(i,x,ep)
		eff_minus.GetPoint(i,x,em)
	
		sys_np = TMath.Abs(float(np) - float(c))
		sys_nm = TMath.Abs(float(nm) - float(c))
		sys_ep = TMath.Sqrt( TMath.Power( TMath.Abs(float(ep) - float(c)),2) + TMath.Power(pdfsys[x]*float(c),2) + TMath.Power(float(c)*Jpsi2mumuBr_relerr,2) )
		sys_em = TMath.Sqrt( TMath.Power( TMath.Abs(float(em) - float(c)),2) + TMath.Power(pdfsys[x]*float(c),2) + TMath.Power(float(c)*Jpsi2mumuBr_relerr,2) )
	
		#print TMath.Abs(float(em) - float(c))*1000, pdfsys[x]*float(c)*1000, 
		if type == 'abs' :
			JpsiLBR = 6.2e-4;
			JpsiLBrrelErr = 1.4/6.2
			
			absErrLow  = TMath.Sqrt( TMath.Power(errlow,2) + TMath.Power(sys_nm,2) + TMath.Power(sys_em,2)  )
			absErrHigh = TMath.Sqrt( TMath.Power(errhigh,2) + TMath.Power(sys_np,2) + TMath.Power(sys_ep,2) )
			normerr    = c*JpsiLBR*JpsiLBrrelErr
			absErrLow  = TMath::Sqrt( TMath.Power(absErrLow * JpsiLBR / c,2) + TMath.Power(normerr,2) )
			absErrHigh = TMath::Sqrt( TMath.Power(absErrHigh * JpsiLBR / c,2) + TMath.Power(normerr,2) )

			result.SetPoint(i,x,c*JpsiLBR)	
			result.SetPointError(i,errX,errX,absErrLow,absErrHigh)
	
			out = '{:5.1f}-{:5.1f}  & ${:10.4f}'.format(float(x)-errX,float(x)+errX, float(c)*1.e3)
			out += '^[+{0:10.4f}]'.format(errhigh*1.e3)
			out += '_[{0:5.4f}] \\text[(stat)] '.format(errlow*1.e3)
			out += '^[+{0:10.4f}]'.format(TMath::Sqrt( TMath::Power(sys_ep,2) + TMath::Power(sys_np,2) )*1.e3)
			out += '_[-{0:5.4f}] \\text[(sys)] '.format(TMath::Sqrt( TMath::Power(sys_em,2) + TMath::Power(sys_nm,2) )*1.e3)
			out += '{0:5.4f} \\text[(norm)] $ \\\\'.format(normerr) 
			print out.replace("[","{").replace("]","}").replace("-    ","-").replace("+    ","+")
		else
			result.SetPoint(i,x,c)
			result.SetPointError(i,errX,errX,
				TMath.Sqrt(TMath.Power(errlow,2) + TMath.Power(sys_nm,2) + TMath.Power(sys_em,2) ),
				TMath.Sqrt(TMath.Power(errhigh,2) + TMath.Power(sys_np,2) + TMath.Power(sys_ep,2) ) )
	
			out = '{:5.1f}-{:5.1f}  & ${:10.4f}'.format(float(x)-errX,float(x)+errX, float(c)*1.e3)
			out += '^[+{0:10.4f}]'.format(errhigh*1.e3)
			out += '_[{0:5.4f}] \\text[(stat)] '.format(errlow*1.e3)
			out += '^[+{0:10.4f}]'.format(sys_ep*1.e3)
			out += '_[-{0:5.4f}] \\text[(sys)] '.format(sys_em*1.e3)
			out += '^[+{0:10.4f}]'.format(sys_np*1.e3)
			out += '_[-{0:5.4f}] \\text[(norm)] $ \\\\'.format(sys_nm*1.e3) 
			print out.replace("[","{").replace("]","}").replace("-    ","-").replace("+    ","+")


	print "\n\n"
	for i in range(0,result.GetN()) :
		x, c = Double(0), Double(0)
		result.GetPoint(i,x,c)
		error = result.GetErrorYlow(i)
		errX = result.GetErrorXlow(i)
		#print c, error
		print '{:5.1f}-{:5.1f}  & ${:10.4f}'.format(float(x)-errX,float(x)+errX,c/error)

	return result;

if __name__ == "__main__" :

	h = TH1F("hh","",1000,0,21)
	gStyle.SetOptStat(0)
	h.SetMinimum(0)
	h.SetMaximum(0.00028)
	h.GetYaxis().SetTitle("(1 / B(#Lambda_{b} #rightarrow J/#psi #Lambda))  dB(#Lambda_{b} #rightarrow #Lambda #mu #mu) / dq^{2}");
	h.GetYaxis().SetTitleSize(0.04)
	h.GetYaxis().SetLabelSize(0.045)
	h.GetYaxis().SetTitleOffset(1.4)
	h.GetXaxis().SetTitle("q^{2} [GeV^{2}/#it{c}^{4}])")
	c = TCanvas()

	if(len(sys.argv)>1) :
		print "############################## DD events ######################################"
		hDD = makeResult("/afs/cern.ch/work/p/pluca/jobs/yield/DD/")
		print "############################## LL events ######################################"
		hLL = makeResult("/afs/cern.ch/work/p/pluca/jobs/yield/LL/")
		h.Draw()
		hDD.SetMarkerColor(2)
		hDD.Draw("P SAME")
		hLL.SetMarkerColor(4)
		hLL.SetMarkerStyle(22)
		hLL.Draw("P SAME")
		leg = TLegend(0.2,0.7,0.5,0.9);
		leg.AddEntry(hLL,"LL","P")
		leg.AddEntry(hDD,"DD","P")
		leg.SetFillColor(kWhite);
		leg.Draw()
		lhcbName.Draw()
		c.Print("q2result_both.pdf")
		
	else :
		hh = makeResult("/afs/cern.ch/work/p/pluca/jobs/yield/comb/")
		h.Draw()
		hh.Draw("P SAME")
		lhcbName.Draw()
		c.Print("combined_result.pdf")
		c.Print("combined_result.eps")

		hh = makeResult("/afs/cern.ch/work/p/pluca/jobs/yield/comb/",'abs')
		h.GetYaxis().SetTitle("dB(#Lambda_{b} #rightarrow #Lambda #mu #mu) / dq^{2}");
		h.Draw()
		hh.Draw("P SAME")
		lhcbName.Draw()
		c.Print("combined_result_absolute.pdf")
		c.Print("combined_result_absolute.eps")

