from ROOT import *
import time
import sys
import os
import subprocess
from array import array

gROOT.ProcessLine(".x ~/work/lhcbStyle.C")

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



c = TCanvas()

file = TFile("/afs/cern.ch/work/p/pluca/results/LbreleffAndSysvsq2_DD_7bins.root")
h_7binsDD = file.Get("toteff");
file = TFile("/afs/cern.ch/work/p/pluca/results/LbreleffAndSysvsq2_DD_2bins.root")
h_2binsDD = file.Get("toteff");
file = TFile("/afs/cern.ch/work/p/pluca/results/LbreleffAndSysvsq2_LL_7bins.root")
h_7binsLL = file.Get("toteff");
file = TFile("/afs/cern.ch/work/p/pluca/results/LbreleffAndSysvsq2_LL_2bins.root")
h_2binsLL = file.Get("toteff");

h_7binsDD.SetMarkerColor(1)
h_7binsDD.SetMarkerStyle(20)
h_7binsDD.SetMarkerSize(1)
h_2binsDD.SetMarkerColor(1)
h_2binsDD.SetMarkerStyle(20)
h_2binsDD.SetMarkerSize(1)

h_7binsLL.SetMarkerColor(2)
h_7binsLL.SetMarkerStyle(21)
h_7binsLL.SetMarkerSize(1)
h_2binsLL.SetMarkerColor(2)
h_2binsLL.SetMarkerStyle(21)
h_2binsLL.SetMarkerSize(1)

h_7binsDD.GetXaxis().SetTitle("q^{2} [GeV^{2}/#it{c}^{4}]")
h_7binsDD.GetYaxis().SetTitleOffset(0.8)
h_7binsDD.GetYaxis().SetTitle("Total relative efficiency")

h_7binsDD.Draw("AP")
h_2binsDD.Draw("P same")
h_7binsLL.Draw("P same")
h_2binsLL.Draw("P same")

leg = TLegend(0.7,0.2,0.9,0.4)
leg.AddEntry(h_7binsDD,"DD","P")
leg.AddEntry(h_7binsLL,"LL","P")
leg.SetFillStyle(0)
leg.Draw()
lhcbName.Draw()

c.Print("Efficiency.pdf")
c.Print("Efficiency.eps")


def getAngGr(type) :
	
	c = TCanvas()

	out = subprocess.check_output("grep res res_"+type+" | tr \"res\" \"   \" | tr \"_{\" \"  \" | tr \"}^{\" \"   \" | tr \"}\" \" \" ", shell=True)

	x = [17.5, 11.75, 15.5, 17, 19, 1.05]
	errx = [2.5, 0.75, 0.5, 1, 1, 0.95]
	vals = []
	errp = []
	errm = []

	for l in out.split("\n") :
		vv = l.split()
		if len(vv) > 0 :
			vals.append(float(vv[0]))
			errm.append(TMath.Abs(float(vv[1])))
			errp.append(TMath.Abs(float(vv[2])))

	gr = TGraphAsymmErrors(6,array('f',x),array('f',vals),array('f',errx),array('f',errx),array('f',errm),array('f',errp))

	gr.GetXaxis().SetTitle("q^{2} [GeV^{2}/#it{c}^{4}]")

	if type == "afb" :
		gr.GetYaxis().SetTitle("A_{FB}^{l}")
		gr.GetYaxis().SetRangeUser(-0.75,0.75)
		gr.Draw("AP")
		lhcbName.Draw()
		c.Print("Afb.pdf")
		c.Print("Afb.eps")

	if type == "afbB" :
		gr.GetYaxis().SetRangeUser(-0.5,0.5)
		gr.GetYaxis().SetTitle("A_{FB}^{h}")
		gr.Draw("AP")
		lhcbName.Draw()
		c.Print("Afbh.pdf")
		c.Print("Afbh.eps")

	if type == "fL" :
		gr.GetYaxis().SetRangeUser(0.,1.)
		gr.GetYaxis().SetTitle("f_{L}")
		gr.Draw("AP")
		lhcbName.Draw()
		c.Print("fL.pdf")
		c.Print("fL.eps")



getAngGr("afb")
getAngGr("afbB")
getAngGr("fL")
