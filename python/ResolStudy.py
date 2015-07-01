#! /usr/bin/env python

from ROOT import *
from sys import argv
from array import array
from string import *
import os

weight = "Lb_weight*MCnorm*physRate_polp006*lifeTimeW"
#file = TFile(argv[1])
#tree = file.Get(argv[2])
file = TFile("/afs/cern.ch/work/p/pluca/Lmumu/weighted/Lb2Lmumu_MC_Pythia8_NBweighted.root")
tree = file.Get("tree")

bins = [11,12.5,15,16,18,20]

c = TCanvas()

truev = "( TMath::Power(muplus_TRUEP_E + muminus_TRUEP_E,2) - TMath::Power(muplus_TRUEP_X + muminus_TRUEP_X,2) - TMath::Power(muplus_TRUEP_Y + muminus_TRUEP_Y,2)  - TMath::Power(muplus_TRUEP_Z + muminus_TRUEP_Z,2) )/1000000 "
recov = "TMath::Power(J_psi_1S_MM/1000,2)"

vars = [recov,"cosThetaB","cosThetaL"]

tree.Draw(truev + " : " + recov + ">>resmatrix","","colz")
hresm = gPad.GetPrimitive("resmatrix")

hresm.GetXaxis().SetTitle( "reco" )
hresm.GetYaxis().SetTitle( "true" )
hresm.SetTitle( recov )
	
c.Print("resolution_matrix_vs.pdf")





for v in vars :

	c = TCanvas()
	gStyle.SetOptStat(0)
	vname = v.replace("(","").replace(")","").replace(":","").replace("/","").replace(",2","_2")

	tree.Draw("(" + truev + " - " + recov + " ) / " + recov + " : " + v + ">>hdiff"+vname,"","colz")
	hdiff = gPad.GetPrimitive("hdiff"+vname)

	hdiff.GetXaxis().SetTitle( v )
	hdiff.GetYaxis().SetTitle("(reco - true)/reco")
	hdiff.SetTitle("")

	c.Print("spread_vs_"+vname+".pdf")


hres = TGraph()
habsres = TGraph()
hbinmig = TGraph()

for b in range(0,len(bins)-1) :

	c = TCanvas()
	vname = recov.replace("(","").replace(")","").replace(":","").replace("/","").replace(",2","_2")

	if weight == "" :
		tree.Draw( "(" + truev + " - " + recov + " ) / " + recov + ">>hrms", recov + " > " + str(bins[b]) + " && " + recov + " < " + str(bins[b+1]) , "colz" )
	else :
		tree.Draw( "(" + truev + " - " + recov + " ) / " + recov + ">>hrms", weight + "*(" + recov + " > " + str(bins[b]) + " && " + recov + " < " + str(bins[b+1]) + ")" , "colz" )

	hdiff = gPad.GetPrimitive("hrms")
	hdiff.GetXaxis().SetTitle("(reco - true)/reco" )
	hdiff.Fit("gaus")
	c.Print("resolution_"+str(bins[b])+"_"+str(bins[b+1])+".pdf")

	func = hdiff.GetFunction("gaus")
	res = func.GetParameter(2)

	hres.SetPoint(b,(bins[b+1]+bins[b])/2,res*100)
	habsres.SetPoint(b,(bins[b+1]+bins[b])/2,res*(bins[b+1]+bins[b])/2)

	c = TCanvas()
	print recov + " > " + str(bins[b]) + " && " + recov + " < " + str(bins[b+1])
	
	in_inside = 0
	in_outside = 0
	out_inside = 0
	
	cut_in_inside = recov + " > " + str(bins[b]) + " && " + recov + " < " + str(bins[b+1]) + " && " + truev + " > " + str(bins[b]) + " && " + truev + " < " + str(bins[b+1])
	print "In inside: ", cut_in_inside
	if weight == "" :
		tree.Draw(truev + ">>hin_inside(1000,0,23)", cut_in_inside, "colz" )
	else :
		tree.Draw(truev + ">>hin_inside(1000,0,23)",weight + "*(" + cut_in_inside + ")", "colz" )

	hin_inside = gPad.GetPrimitive("hin_inside")
	in_inside = hin_inside.Integral()


	cut_in_all = truev + " > " + str(bins[b]) + " && " + truev + " < " + str(bins[b+1])
	print "In all: ", cut_in_all
	if weight == "" :
		tree.Draw(truev + ">>hin_all(1000,0,23)", cut_in_all , "colz" )
	else :
		tree.Draw(truev + ">>hin_all(1000,0,23)",weight + "*(" + cut_in_all + ")", "colz" )

	hin_all = gPad.GetPrimitive("hin_all")
	in_all = hin_all.Integral()
	in_outside = in_all - in_inside;


	cut_out_inside = "(" + truev + " < " + str(bins[b]) + " || " + truev + " > " + str(bins[b+1]) + ") && " + recov + " > " + str(bins[b]) + " && " + recov + " < " + str(bins[b+1])
	print "Out inside : ", cut_out_inside

	if weight == "" :
		tree.Draw(truev + ">>hout_inside(1000,0,23)", cut_out_inside , "colz" )
	else :
		tree.Draw(truev + ">>hout_inside(1000,0,23)",weight + "*(" + cut_out_inside + ")", "colz" )

	hout_inside = gPad.GetPrimitive("hout_inside")
	out_inside = hout_inside.Integral()


	binmigration = (out_inside - in_outside) / in_all
	hbinmig.SetPoint(b,(bins[b+1]+bins[b])/2,binmigration*100)

habsres.SetMarkerColor(1)
habsres.SetMarkerStyle(21)
habsres.SetMarkerSize(1)

hres.SetMarkerColor(1)
hres.SetMarkerStyle(21)
hres.SetMarkerSize(1)

hbinmig.SetMarkerColor(1)
hbinmig.SetMarkerStyle(21)
hbinmig.SetMarkerSize(1)

habsres.GetXaxis().SetTitle(recov)
habsres.GetYaxis().SetTitle("resolution")
habsres.Draw("AP")
c.Print("absolute_resolution_vs_"+vname+".pdf")
hres.GetXaxis().SetTitle(recov)
hres.GetYaxis().SetTitle("resolution (%)")
hres.Draw("AP")
c.Print("resolution_vs_"+vname+".pdf")
hbinmig.GetYaxis().SetTitle("bin migration (%)")
hbinmig.GetXaxis().SetTitle(recov)
hbinmig.Draw("AP")
c.Print("binmigration_vs_"+vname+".pdf")








