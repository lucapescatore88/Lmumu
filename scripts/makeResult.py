from ROOT import *
from string import *
import sys

gROOT.ProcessLine(".x ~/work/lhcbStyle.C");

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


def makeResult(base, type = 'rel', err = 'tot') :	

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


	if(type!="sys") :
		print "\n\nq2 bin \t\t & Combined yield \\\\"

	p = 0

	for i in range(0,central.GetN()) :
		x, c = Double(0), Double(0)
		np, nm = Double(0), Double(0)
		ep, em = Double(0), Double(0)
		central.GetPoint(i,x,c)

		if x < 1e-9 and c < 1e-9 :
			continue

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
		sys_toteff_p = TMath.Abs(float(ep) - float(c))/c
		sys_toteff_m = TMath.Abs(float(em) - float(c))/c

		#print TMath.Abs(float(em) - float(c))*1000, pdfsys[x]*float(c)*1000, 
		if type == 'abs' :
			JpsiLBR = 6.2e-4;
			JpsiLBrrelErr = 1.4/6.2
			
			absErrLow  = TMath.Sqrt( TMath.Power(errlow,2) + TMath.Power(sys_nm,2) + TMath.Power(sys_em,2)  )
			absErrHigh = TMath.Sqrt( TMath.Power(errhigh,2) + TMath.Power(sys_np,2) + TMath.Power(sys_ep,2) )
			normerr    = c*JpsiLBR*JpsiLBrrelErr
			if err=='tot':
				absErrLow  = TMath.Sqrt( TMath.Power(absErrLow * JpsiLBR,2) + TMath.Power(normerr,2) )
				absErrHigh = TMath.Sqrt( TMath.Power(absErrHigh * JpsiLBR,2) + TMath.Power(normerr,2) )
			if err=='norm' :
				absErrLow  = normerr
				absErrHigh = normerr
			if err=='rare':
				absErrLow  = absErrLow * JpsiLBR
				absErrHigh = absErrHigh * JpsiLBR

			if x != 17.5 and x != 3.55 :
				
				result.SetPoint(p,x,c*JpsiLBR*1e6)	
				result.SetPointError(p,errX,errX,absErrLow*1e6,absErrHigh*1e6)
				p+=1

			out = '{:.1f}--{:.1f}  & ${:^10.2f}'.format(float(x)-errX,float(x)+errX, float(c*JpsiLBR)*1.e7)
			out += '\\;^[+{0:.2f}]'.format(errhigh*JpsiLBR*1.e7)
			out += '_[{0:.2f}]\\; \\text[(stat)]\\; '.format(errlow*JpsiLBR*1.e7)
			
			out += '^[+{0:.2f}]'.format(TMath.Sqrt( TMath.Power(sys_ep*JpsiLBR,2) + TMath.Power(sys_np*JpsiLBR,2) )*1.e7)
			out += '_[-{0:.2f}]\\; \\text[(sys)]\\; '.format(TMath.Sqrt( TMath.Power(sys_em*JpsiLBR,2) + TMath.Power(sys_nm*JpsiLBR,2) )*1.e7)
	
			out += '^[+{0:.2f}]'.format(normerr*1.e7)
			out += '_[-{0:.2f}]\\; \\text[(norm)] $ \\\\'.format(normerr*1.e7)

			print out.replace("[","{").replace("]","}").replace("-    ","-").replace("+    ","+")
		elif type == 'sys' :
			if(i==0) :
				print "q^{2} bin \t& efficiency \t\t & Yields "
			out = '{:.1f}--{:.1f}  \t & '.format(float(x)-errX,float(x)+errX)
			out+= '[-{:.1f},{:.1f}]\\% \t & '.format(sys_toteff_m*100,sys_toteff_p*100)
			out+= '{:.1f}\\% \\\\'.format(pdfsys[x]*100)
			print out

		else :
			
			if x != 17.5 and x != 3.55 :
	
				result.SetPoint(p,x,c*1e3)
				if err=='tot' :
					result.SetPointError(p,errX,errX,
						TMath.Sqrt(TMath.Power(errlow,2) + TMath.Power(sys_nm,2) + TMath.Power(sys_em,2) )*1e3,
						TMath.Sqrt(TMath.Power(errhigh,2) + TMath.Power(sys_np,2) + TMath.Power(sys_ep,2) )*1e3 )
				elif err=='sys' :
					result.SetPointError(p,errX,errX,
                         TMath.Sqrt(TMath.Power(sys_nm,2) + TMath.Power(sys_em,2) )*1e3,
                         TMath.Sqrt(TMath.Power(sys_np,2) + TMath.Power(sys_ep,2) )*1e3 )
				elif err=='stat' :
					result.SetPointError(p,errX,errX,errlow*1e3,errhigh*1e3)
				p+=1

		
			out = '{:.1f}--{:.1f}  & ${:^10.2f}\\;'.format(float(x)-errX,float(x)+errX, float(c)*1.e4)
			out += '^[+{0:.2f}]'.format(errhigh*1.e4)
			out += '_[{0:.2f}]\\; \\text[(stat)]\\; '.format(errlow*1.e4)
			
			out += '^[+{0:.2f}]'.format( TMath.Sqrt( TMath.Power(sys_ep,2) + TMath.Power(sys_np,2) )*1.e4 )
			out += '_[-{0:.2f}]\\; \\text[(sys)]  $ \\\\'.format( TMath.Sqrt( TMath.Power(sys_em,2) + TMath.Power(sys_nm,2) )*1.e4 )

			#out += '^[+{0:10.4f}]'.format(sys_ep*1.e3)
			#out += '_[-{0:5.4f}] \\text[(sys)] '.format(sys_em*1.e3)
			#out += '^[+{0:10.4f}]'.format(sys_np*1.e3)
			#out += '_[-{0:5.4f}] \\text[(norm)] $ \\\\'.format(sys_nm*1.e3)
			print out.replace("[","{").replace("]","}").replace("-  ","-").replace("+  ","+")

	print "\n\n"
	#for i in range(0,result.GetN()) :
	#	x, c = Double(0), Double(0)
	#	result.GetPoint(i,x,c)
	#	error = result.GetErrorYlow(i)
	#	errX = result.GetErrorXlow(i)
	#	print c, error
	#	print '{:5.1f}-{:5.1f}  & ${:10.4f}'.format(float(x)-errX,float(x)+errX,c/error)

	return result;

if __name__ == "__main__" :

	h = TH1F("hh","",1000,0,21)
	gStyle.SetOptStat(0)
	h.SetMinimum(0)
	h.SetMaximum(0.28)
	h.GetYaxis().SetTitle("(1 / B(#Lambda_{b} #rightarrow J/#psi #Lambda))  dB(#Lambda_{b} #rightarrow #Lambda #mu #mu) / d#it{q}^{2} #upoint 10^{-3}");
	h.GetYaxis().SetTitleSize(0.04)
	h.GetYaxis().SetLabelSize(0.045)
	h.GetYaxis().SetTitleOffset(1.4)
	h.GetXaxis().SetTitle("#it{q}^{2} [GeV^{2}/#it{c}^{4}]")
	c = TCanvas()

	gStyle.SetEndErrorSize(4)

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
		hh = makeResult("/afs/cern.ch/work/p/pluca/jobs/yield_local/comb/")
		h.Draw()
		hh.Draw("P SAME")
		hh_stat = makeResult("/afs/cern.ch/work/p/pluca/jobs/yield_local/comb/",'rel','sys')
		hh_stat.Draw("P SAME")
		lhcbName.Draw()
		c.Print("combined_result_2err.pdf")
		c.Print("combined_result_2err.eps")

		hh_tot = makeResult("/afs/cern.ch/work/p/pluca/jobs/yield_local/comb/",'abs')
		h.SetMaximum(0.18)
		h.GetYaxis().SetTitle("dB(#Lambda_{b} #rightarrow #Lambda #mu #mu) / d#it{q}^{2} #upoint 10^{-6}");
		h.Draw()

		hh_tot.Draw("P SAME")
		hh_norm = makeResult("/afs/cern.ch/work/p/pluca/jobs/yield_local/comb/",'abs','norm')
		hh_norm.Draw("P SAME")
		lhcbName.Draw()
		c.Print("combined_result_absolute_2err.pdf")
		c.Print("combined_result_absolute_2err.eps")

		makeResult("/afs/cern.ch/work/p/pluca/jobs/yield_local/comb/",'sys')

