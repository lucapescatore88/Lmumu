from ROOT import *
from string import *
import sys

gROOT.ProcessLine(".x ~/work/lhcbStyle.C");

from ROOT import TPaveText, gStyle

lhcbName = TPaveText(gStyle.GetPadLeftMargin() + 0.05,
                        0.80 - gStyle.GetPadTopMargin(),
                        gStyle.GetPadLeftMargin() + 0.30,
                        0.95 - gStyle.GetPadTopMargin(),
                        "BRNDC");

lhcbName.AddText("LHCb");
lhcbName.AddText("simulation");
lhcbName.SetFillColor(0);
lhcbName.SetTextAlign(12);
lhcbName.SetBorderSize(0);

base = "$WORK/results/"
file = TFile(base+"LbreleffAndSysvsq2_DD_2bins.root")
effDD_2bin = file.Get("toteff")
file = TFile(base+"LbreleffAndSysvsq2_LL_2bins.root")
effLL_2bin = file.Get("toteff")
file = TFile(base+"LbreleffAndSysvsq2_DD_7bins.root")
effDD_7bin = file.Get("toteff")
file = TFile(base+"LbreleffAndSysvsq2_LL_7bins.root")
effLL_7bin = file.Get("toteff")
file = TFile(base+"LbreleffAndSysvsq2_DD_2bins.root")
effDD_2bin_low = file.Get("toteff_lowSel")
file = TFile(base+"LbreleffAndSysvsq2_LL_2bins.root")
effLL_2bin_low = file.Get("toteff_lowSel")
file = TFile(base+"LbreleffAndSysvsq2_DD_7bins.root")
effDD_7bin_low = file.Get("toteff_lowSel")
file = TFile(base+"LbreleffAndSysvsq2_LL_7bins.root")
effLL_7bin_low = file.Get("toteff_lowSel")


grDD = TGraphAsymmErrors()
grLL = TGraphAsymmErrors()
p = 0

#for i in range(0,effDD_2bin.GetN()) :
#	x, c = Double(0), Double(0)
#	effDD_2bin.GetPoint(i,x,c)
#	errlow = effDD_2bin.GetErrorYlow(i)
#	errhigh = effDD_2bin.GetErrorYhigh(i)
#	errX = effDD_2bin.GetErrorX(i)

#	if(x > 8 and c > 1e-6) :
#		grDD.SetPoint(p,x,c)
#		grDD.SetPointError(p,errX,errX,errlow,errhigh)
#		p+=1

for i in range(0,effDD_7bin.GetN()) :
	x, c = Double(0), Double(0)
	effDD_7bin.GetPoint(i,x,c)
	errlow = effDD_7bin.GetErrorYlow(i)
	errhigh = effDD_7bin.GetErrorYhigh(i)
	errX = effDD_7bin.GetErrorX(i)

	if(x > 10 and c > 1e-6) :
		grDD.SetPoint(p,x,c)
		grDD.SetPointError(p,errX,errX,errlow,errhigh)
		p+=1

#for i in range(0,effDD_2bin_low.GetN()) :
#	x, c = Double(0), Double(0)
#	effDD_2bin_low.GetPoint(i,x,c)
#	errlow = effDD_2bin_low.GetErrorYlow(i)
#	errhigh = effDD_2bin_low.GetErrorYhigh(i)
#	errX = effDD_2bin_low.GetErrorX(i)

#	if(x < 8 and c > 1e-6) :
#		grDD.SetPoint(p,x,c)
#		grDD.SetPointError(p,errX,errX,errlow,errhigh)
#		p+=1

for i in range(0,effDD_7bin_low.GetN()) :
	x, c = Double(0), Double(0)
	effDD_7bin_low.GetPoint(i,x,c)
	errlow = effDD_7bin_low.GetErrorYlow(i)
	errhigh = effDD_7bin_low.GetErrorYhigh(i)
	errX = effDD_7bin_low.GetErrorX(i)

	if(x < 8 and c > 1e-6) :
		grDD.SetPoint(p,x,c)
		grDD.SetPointError(p,errX,errX,errlow,errhigh)
		p+=1


p = 0

#for i in range(0,effLL_2bin.GetN()) :
#	x, c = Double(0), Double(0)
#	effLL_2bin.GetPoint(i,x,c)
#	errlow = effLL_2bin.GetErrorYlow(i)
#	errhigh = effLL_2bin.GetErrorYhigh(i)
#	errX = effLL_2bin.GetErrorX(i)

#	if(x > 10 and c > 1e-6) :
#		grLL.SetPoint(p,x,c)
#		grLL.SetPointError(p,errX,errX,errlow,errhigh)
#		p+=1

for i in range(0,effLL_7bin.GetN()) :
	x, c = Double(0), Double(0)
	effLL_7bin.GetPoint(i,x,c)
	errlow = effLL_7bin.GetErrorYlow(i)
	errhigh = effLL_7bin.GetErrorYhigh(i)
	errX = effLL_7bin.GetErrorX(i)

	if(x > 10 and c > 1e-6) :
		grLL.SetPoint(p,x,c)
		grLL.SetPointError(p,errX,errX,errlow,errhigh)
		p+=1
		
#for i in range(0,effLL_2bin_low.GetN()) :
#	x, c = Double(0), Double(0)
#	effLL_2bin_low.GetPoint(i,x,c)
#	errlow = effLL_2bin_low.GetErrorYlow(i)
#	errhigh = effLL_2bin_low.GetErrorYhigh(i)
#	errX = effLL_2bin_low.GetErrorX(i)

#	if(x < 8 and c > 1e-6) :
#		grLL.SetPoint(p,x,c)
#		grLL.SetPointError(p,errX,errX,errlow,errhigh)
#		p+=1

for i in range(0,effLL_7bin_low.GetN()) :
	x, c = Double(0), Double(0)
	effLL_7bin_low.GetPoint(i,x,c)
	errlow = effLL_7bin_low.GetErrorYlow(i)
	errhigh = effLL_7bin_low.GetErrorYhigh(i)
	errX = effLL_7bin_low.GetErrorX(i)

	if(x < 8 and c > 1e-6) :
		grLL.SetPoint(p,x,c)
		grLL.SetPointError(p,errX,errX,errlow,errhigh)
		p+=1



c = TCanvas()

grDD.SetMaximum(1.4)
grDD.SetMinimum(0.4)
grDD.SetMarkerStyle(20)
grDD.SetMarkerSize(1)
grDD.SetMarkerColor(1)
grLL.SetMarkerStyle(25)
grLL.SetMarkerSize(1)
grLL.SetMarkerColor(1)

leg = TLegend(0.65,0.2,0.9,0.4)
leg.AddEntry(grLL,"long","P")
leg.AddEntry(grDD,"downstream","P")
leg.SetFillColor(0)

grDD.GetXaxis().SetTitle("#it{q}^{2} [GeV^{2}/#it{c}^{4}]")
grDD.GetYaxis().SetTitle("Relative efficiency")

grDD.Draw("AP")
grLL.Draw("P same")
leg.Draw()
lhcbName.Draw()

c.Print("Efficiency.pdf")
c.Print("Efficiency.eps")
c.Print("Efficiency.png")
c.Print("Efficiency.C")





