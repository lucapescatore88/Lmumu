from ROOT import *
from string import *
import sys

gROOT.ProcessLine(".x ~/work/lhcbStyle.C");

from ROOT import TPaveText, gStyle

lhcbName = TPaveText(gStyle.GetPadLeftMargin() + 0.05,
                        0.82 - gStyle.GetPadTopMargin(),
                        gStyle.GetPadLeftMargin() + 0.25,
                        0.95 - gStyle.GetPadTopMargin(),
                        "BRNDC");

lhcbName.AddText("LHCb");
lhcbName.SetFillColor(0);
lhcbName.SetTextAlign(12);
lhcbName.SetBorderSize(0);

lhcbName_low = TPaveText(gStyle.GetPadLeftMargin() + 0.58,
                        0.25 - gStyle.GetPadTopMargin(),
                        gStyle.GetPadLeftMargin() + 0.8,
                        0.37 - gStyle.GetPadTopMargin(),
                        "BRNDC");

lhcbName_low.AddText("LHCb");
lhcbName_low.SetFillColor(0);
lhcbName_low.SetTextAlign(12);
lhcbName_low.SetBorderSize(0);

## Stefan values
BrTheo = { 17.5 : 0.76, 11.75 : 1.004, 15.5 : 1.035, 17.0 : 0.956, 19.0 : 0.656, 3.55 : 0.52, 1.05 : 0.653, 3.0 : 0.556, 5.0 : 0.655, 7.0 : 0.817}
BrTheo_err = { 17.5 : 0.18, 11.75 : 0.36, 15.5 : 0.29, 17.0 : 0.24, 19.0 : 0.14, 3.55 : 0.27, 1.05 : 0.38, 3.0 : 0.30, 5.0 : 0.32, 7.0 : 0.37}

## Michal values 
#BrTheo = { 17.5 : 0.76, 11.75 : 0.89, 15.5 : 0.92, 17.0 : 0.85, 19.0 : 0.58, 3.55 : 0.52, 1.05 : 0.58, 3.0 : 0.49, 5.0 : 0.58, 7.0 : 0.72}
#BrTheo_err = { 17.5 : 0.18, 11.75 : 0.32, 15.5 : 0.25, 17.0 : 0.21, 19.0 : 0.12, 3.55 : 0.27, 1.05 : 0.33, 3.0 : 0.26, 5.0 : 0.29, 7.0 : 0.33}

q2_err = { 17.5 : 2.5, 11.75 : 0.75, 15.5 : 0.5, 17.0 : 1., 19.0 : 1., 3.55 : 2.95, 1.05 : 0.95, 3.0 : 1., 5.0 : 1., 7.0 : 1.}
BRTheo_gr = TGraphErrors()

i = 0
for k,v in BrTheo.iteritems() :
	if k != 17.5 and k != 3.55 :
		BRTheo_gr.SetPoint(i,k,v)
		i+=1
i = 0
for k,v in BrTheo_err.iteritems() :
	if k != 17.5 and k != 3.55 :
		BRTheo_gr.SetPointError(i,q2_err[k],v)
		i+=1

def makeResult(base, type = 'rel', err = 'tot') :	

	file = TFile(base+"Lb_yield/Lbyield.root")
	central = file.Get("BRplot")
	file = TFile(base+"Lb_yield_norm_plus/Lbyield.root")
	norm_plus = file.Get("BRplot")
	file = TFile(base+"Lb_yield_norm_minus/Lbyield.root")
	norm_minus = file.Get("BRplot")
	base=base.replace("yield_local","yield")
	file = TFile(base+"Lb_yield_eff_plus_noWC/Lbyield.root")
	eff_plus = file.Get("BRplot")
	file = TFile(base+"Lb_yield_eff_minus_noWC/Lbyield.root")
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
			#JpsiLBR = 6.2e-4
			#JpsiLBrrelErr = 1.4/6.2
			JpsiLBR = 6.3e-4;
			JpsiLBrrelErr = 1.3/6.3
			
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
				
				result.SetPoint(p,x,c*JpsiLBR*1e7)	
				result.SetPointError(p,errX,errX,absErrLow*1e7,absErrHigh*1e7)
				p+=1

			out = '{:.1f}--{:.1f}  & ${:^10.2f} \\;\\;\\;'.format(float(x)-errX,float(x)+errX, float(c*JpsiLBR)*1.e7)
			out += '^[+\\,{0:.2f}]'.format(errhigh*JpsiLBR*1.e7)
			out += '_[-\\,{0:.2f}]\\;\\;\\;'.format(TMath.Abs(errlow)*JpsiLBR*1.e7)
			
			out += '^[+\\,{0:.2f}]'.format(TMath.Sqrt( TMath.Power(sys_ep*JpsiLBR,2) + TMath.Power(sys_np*JpsiLBR,2) )*1.e7)
			out += '_[-\\,{0:.2f}]\\;\\;\\;'.format(TMath.Sqrt( TMath.Power(sys_em*JpsiLBR,2) + TMath.Power(sys_nm*JpsiLBR,2) )*1.e7)
	
			out += '\\pm\\,{0:.2f}\\; \\\\'.format(normerr*1.e7)

			print out.replace("[","{").replace("]","}").replace("-    ","-").replace("+    ","+")
		elif type == 'sys' :
			if(i==0) :
				print "q^{2} bin \t& Yields \t\t & Efficiency "
			out  = '{:.1f}--{:.1f}  \t & '.format(float(x)-errX,float(x)+errX)
			out += '{:.1f} '.format(pdfsys[x]*100)
			#out+= '[-{:.1f},{:.1f}]\\% \t & '.format(sys_toteff_m*100,sys_toteff_p*100)
			out += ' &\t $_[-{:.1f}],^[{:.1f}]$ \t '.format(sys_toteff_m*100,sys_toteff_p*100)
			out += ' \\\\'
			print out.replace("[","{").replace("]","}")

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

			out = '{:.1f}--{:.1f}  & ${:^10.2f} \\;\\;\\;'.format(float(x)-errX,float(x)+errX, float(c)*1.e4)
			out += '^[+\\,{0:.2f}]'.format(errhigh*1.e4)
			out += '_[-\\,{0:.2f}]\\;\\;\\; '.format(TMath.Abs(errlow*1.e4))
			
			out += '^[+\\,{0:.2f}]'.format( TMath.Sqrt( TMath.Power(sys_ep,2) + TMath.Power(sys_np,2) )*1.e4 )
			out += '_[-\\,{0:.2f}]\\; $ \\\\'.format( TMath.Sqrt( TMath.Power(sys_em,2) + TMath.Power(sys_nm,2) )*1.e4 )

			print out.replace("[","{").replace("]","}").replace("-  ","-").replace("+  ","+")

	print "\n\n"

	return result;

if __name__ == "__main__" :

	h = TH1F("hh","",1000,0,21)
	gStyle.SetOptStat(0)
	h.SetMinimum(1e-3)
	h.SetMaximum(0.28)
	h.GetYaxis().SetTitle("(1 / #it{B}(#Lambda_{b} #rightarrow J/#psi #Lambda))  d#it{B}(#Lambda_{b} #rightarrow #Lambda #mu #mu) / d#it{q}^{2} #upoint 10^{-3}");
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
		h.GetYaxis().SetLabelSize(0.055)
		h.GetYaxis().SetTitleOffset(1.5)
		h.Draw()
		hh.Draw("P SAME")
		hh_stat = makeResult("/afs/cern.ch/work/p/pluca/jobs/yield_local/comb/",'rel','sys')
		hh_stat.Draw("P SAME")
		lhcbName.Draw()
		c.Print("combined_result_2err.pdf")
		c.Print("combined_result_2err.eps")
		c.Print("combined_result_2err.png")
		c.Print("combined_result_2err.C")

		#BRTheo_gr.Draw("AP")
		#c.Print("prediction.pdf")

		hh_tot = makeResult("/afs/cern.ch/work/p/pluca/jobs/yield_local/comb/",'abs')
		h.SetMaximum(1.8)
		h.GetYaxis().SetLabelSize(0.055)
		h.GetYaxis().SetTitleSize(0.05)
		h.GetYaxis().SetTitleOffset(1.2)
		h.GetYaxis().SetTitle("d#it{B}(#Lambda_{b} #rightarrow #Lambda #mu #mu) / d#it{q}^{2} [10^{-7}(GeV^{2}/#it{c}^{4})^{-1}]");
		h.Draw()
		
		## Theory prediction
		BRTheo_gr.SetMarkerSize(0)
		BRTheo_gr.SetMarkerColor(4)
		BRTheo_gr.SetLineWidth(0)
		BRTheo_gr.SetLineColor(7)
		BRTheo_gr.SetFillColor(7)
		BRTheo_gr.Draw("P SAME 2")
		BRTheo_gr_tmp = BRTheo_gr.Clone()
		for p in xrange(BRTheo_gr.GetN()) :
			BRTheo_gr_tmp.SetPointError(p,BRTheo_gr_tmp.GetErrorX(p),0)
		BRTheo_gr_tmp.SetMarkerSize(0)
		BRTheo_gr_tmp.SetLineColor(4)
		BRTheo_gr_tmp.SetLineStyle(kDashed)

		BRTheo_gr_tmp.Draw("P SAME")
		hh_tot.Draw("P SAME")
		hh_norm = makeResult("/afs/cern.ch/work/p/pluca/jobs/yield_local/comb/",'abs','norm')
		hh_norm.Draw("P SAME")
		lhcbName_low.Draw()
		
		leg_abs = TLegend(0.2,0.65,0.5,0.9)
		leg_abs.AddEntry(BRTheo_gr,"SM prediction","F")
		leg_abs.AddEntry(hh_tot,"Data","PL")
		leg_abs.SetFillColor(0)
		leg_abs.Draw()
		
		gPad.RedrawAxis()
		
		c.Print("combined_result_absolute_2err.pdf")
		c.Print("combined_result_absolute_2err.eps")
		c.Print("combined_result_absolute_2err.png")
		c.Print("combined_result_absolute_2err.C")

		makeResult("/afs/cern.ch/work/p/pluca/jobs/yield_local/comb/",'sys')

	## Absolute with swapped errors

		hh_tot = makeResult("/afs/cern.ch/work/p/pluca/jobs/yield_local/comb/",'abs')
		h.SetMaximum(1.8)
		h.GetYaxis().SetLabelSize(0.055)
		h.GetYaxis().SetTitleSize(0.05)
		h.GetYaxis().SetTitleOffset(1.2)
		h.GetYaxis().SetTitle("d#it{B}(#Lambda_{b} #rightarrow #Lambda #mu #mu) / d#it{q}^{2} [10^{-7}(GeV^{2}/#it{c}^{4})^{-1}]");
		h.Draw()
		
		## Theory prediction
		BRTheo_gr.SetMarkerSize(0)
		BRTheo_gr.SetMarkerColor(4)
		BRTheo_gr.SetLineWidth(0)
		BRTheo_gr.SetLineColor(7)
		BRTheo_gr.SetFillColor(7)
		BRTheo_gr.Draw("P SAME 2")
		BRTheo_gr_tmp = BRTheo_gr.Clone()
		for p in xrange(BRTheo_gr.GetN()) :
			BRTheo_gr_tmp.SetPointError(p,BRTheo_gr_tmp.GetErrorX(p),0)
		BRTheo_gr_tmp.SetMarkerSize(0)
		BRTheo_gr_tmp.SetLineColor(4)
		BRTheo_gr_tmp.SetLineStyle(kDashed)

		BRTheo_gr_tmp.Draw("P SAME")
		hh_tot.Draw("P SAME")
		hh_norm = makeResult("/afs/cern.ch/work/p/pluca/jobs/yield_local/comb/",'abs','rare')
		hh_norm.Draw("P SAME")
		lhcbName_low.Draw()
	
		chi2 = 0
		chi2_low = 0
		chi2_high = 0
		ndf = 0

		for i in xrange(hh_tot.GetN()) :
			x, v = Double(0), Double(0)
			ep, em = Double(0), Double(0)
			hh_tot.GetPoint(i,x,v)
			errlow = hh_tot.GetErrorYlow(i)
			errhigh = hh_tot.GetErrorYhigh(i)

			for j in xrange(BRTheo_gr.GetN()) :
				xt, vt = Double(0), Double(0)
				ept, emt = Double(0), Double(0)
				BRTheo_gr.GetPoint(j,xt,vt)
				errlow_t = BRTheo_gr.GetErrorYlow(j)
				errhigh_t = BRTheo_gr.GetErrorYhigh(j)

				if(abs(xt - x) < 0.1) :
					ndf+=1
					chi2+=(vt-v)*(vt-v)/(TMath.Power((errlow+errhigh)/2.,2) + TMath.Power((errlow_t+errhigh_t)/2.,2))
					
					if errlow > errhigh :
						chi2_high+=(vt-v)*(vt-v)/(TMath.Power(errlow,2) + TMath.Power((errlow_t+errhigh_t)/2.,2))
						chi2_low+=(vt-v)*(vt-v)/(TMath.Power(errhigh,2) + TMath.Power((errlow_t+errhigh_t)/2.,2))
					else :
						chi2_high+=(vt-v)*(vt-v)/(TMath.Power(errhigh,2) + TMath.Power((errlow_t+errhigh_t)/2.,2))
						chi2_low+=(vt-v)*(vt-v)/(TMath.Power(errlow,2) + TMath.Power((errlow_t+errhigh_t)/2.,2))

		print "Chi2 theory and data -> ", chi2, '(ndf = {:3})'.format(ndf), '    ----->  chi2/ndf = ', chi2/ndf, '   Prob = ', TMath.Prob(chi2,ndf)	
		print "Chi2 theory and data (using smallest error) -> ", chi2_low, '(ndf = {:3})'.format(ndf), '    ----->  chi2/ndf = ', chi2_low/ndf, '   Prob = ', TMath.Prob(chi2_low,ndf)
		print "Chi2 theory and data (using biggest error) -> ", chi2_high, '(ndf = {:3})'.format(ndf), '    ----->  chi2/ndf = ', chi2_high/ndf, '   Prob = ', TMath.Prob(chi2_high,ndf)

		leg_abs = TLegend(0.2,0.7,0.4,0.9)
		leg_abs.AddEntry(BRTheo_gr,"SM prediction","F")
		leg_abs.AddEntry(hh_tot,"Data","PL")
		leg_abs.SetFillColor(0)
		leg_abs.Draw()
		
		gPad.RedrawAxis()
		
		c.Print("combined_result_absolute_2err_rareErrIn.pdf")
		c.Print("combined_result_absolute_2err_rareErrIn.eps")
		c.Print("combined_result_absolute_2err_rareErrIn.png")
		c.Print("combined_result_absolute_2err_rareErrIn.C")

		makeResult("/afs/cern.ch/work/p/pluca/jobs/yield_local/comb/",'sys')



