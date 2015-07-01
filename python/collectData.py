#! /usr/bin/env python

import sys
from ROOT import *
import os
from array import array
from string import find
import subprocess

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



dict = { 'B1' : '1500_2000', 'B2' : '1100_1250', 'B3' : '1500_1600', 'B4' : '1600_1800', 'B5' : '1800_2000', 'B6' : '010_200' }
tdict = { 'B1' : '15.0 < q^{2} < 20.0', 'B2' : '11.0 < q^{2} < 12.5', 'B3' : '15.0 < q^{2} < 16.0', 'B4' : '16.0 < q^{2} < 18.0', 'B5' : '18.0 < q^{2} < 20.0', 'B6' : '0.1 < q^{2} < 2.0' }

def extractLimits(gr, CL) :

	entr = int(gr.GetNbinsX())
	
	max, pvmax = -1e6, -1e6
	xmin, xmax, pvxmax, pvxmin = 0, 0, 0, 0
	pv, x = -1, -1
	for i in range(1,entr+1) :
		pv = gr.GetBinContent(i)
		x = gr.GetBinCenter(i)
		if(i==1) :
			xmin = x
			pvxmin = pv
		if(i==(entr)) :
			xmax = x
			pvxmax = pv
		if(pv > pvmax)  :
			pvmax = pv
			max = x
		
	lp, rp = -1, -1;
	for i in range(1,entr+1) :
		pv = gr.GetBinContent(i);
		x = gr.GetBinCenter(i)
		if(x < max and pv >= (1 - CL) and lp==-1.) :
			lp = i;
		if(x > max and pv < (1 - CL) and rp==-1.) :
			rp = i;
		
	res = []
	x1, p1, x2, p2 = 0, 0, 0, 0
	if(max != xmin and pvxmin < (1.-CL)) :
		p1 = gr.GetBinContent(lp-1);
		p2 = gr.GetBinContent(lp);
		x1 = gr.GetBinCenter(lp-1);
		x2 = gr.GetBinCenter(lp);
		res.append(x1 + (x2-x1)*((1.-CL)-p1)/(p2-p1))
	else :
		res.append(xmin)
	if(max != xmax and pvxmax < (1.-CL) ) :
		p1 = gr.GetBinContent(rp-1);
		p2 = gr.GetBinContent(rp);
		x1 = gr.GetBinCenter(rp-1);
		x2 = gr.GetBinCenter(rp);
		res.append(x1 + (x2-x1)*((1.-CL)-p1)/(p2-p1));
	else :
		res.append(xmax);
	
	return res;




def run_command(command):
	return subprocess.Popen(command,stdout=subprocess.PIPE).communicate()

def getFloat(s) :
	l = []
	for t in s.split():
		try:
			l.append(float(t))
		except ValueError:
			pass
	return l

def drange(start, stop, nbins):
	res = []
	r = start
	step = (stop - start)/nbins
	if start==0 :
		step = 0.1
	
	while r <= stop:
		res.append(r)
		r += step
	return res



##################################################################

maxAfb = 0.75
minAfb = -0.75
maxfL = 1
minfL = 0

gRoot->ProcessLine(".x ~/work/lhcbStyle.C")
base = "/afs/cern.ch/work/p/pluca/jobs/"
#base = "/afs/cern.ch/user/p/pluca/eos/lhcb/user/p/pluca/jobs/"

folder = sys.argv[1]
obs = sys.argv[2]
namefile = "FC_"+sys.argv[2]+"_"+folder[-2:]+".dat"
if(sys.argv[1]=="-d" and len(sys.argv) > 2) :
	print "Using already existing data"
	folder = sys.argv[2]
	obs = sys.argv[3]
	namefile = "FC_"+obs+"_"+sys.argv[2][-2:]+".dat"
else :
	os.system( "grep Pvalues: " + base+sys.argv[1] + "/*/out | awk '{ print $4, $5, $6 }' > " + namefile )

stepfile = open(base+folder+"/steps.txt","r")
line1 = stepfile.readline()
line2 = stepfile.readline()
stepAfb = getFloat(line1)[0]
stepfL = getFloat(line2)[0]

nfLbins = (maxfL - minfL)/stepfL + 1
nAfbins = (maxAfb - minAfb)/stepAfb + 1

is1D = False
if find(folder,"1D") > -1 :
	is1D = True

if obs == "afbB" :
	minAfb = -0.5
	maxAfb = 0.5
	nAfbins = (maxAfb - minAfb)/stepAfb + 1
	print "Afbh step = ", stepAfb, "  (", nAfbins, ")"
elif obs == "afb" :
	minAfb = -0.75
	maxAfb = 0.75
	nAfbins = (maxAfb - minAfb)/stepAfb + 1
	print "Afb step = ", stepAfb, "  (", nAfbins, ")"
elif obs == "fL" :
	minAfb = 0.
	maxAfb = 1.
	nAfbins = (maxAfb - minAfb)/stepAfb + 1
	print "fL step = ", stepAfb, "  (", nAfbins, ")"
else :
	print "Afb step = ", stepAfb, ",  fL step = ", stepfL

best = run_command(["/bin/grep","Best",base+folder+"/1/out"])[0].replace("\n","").replace(" ***** Best fit:","") 
best_coord = getFloat(best)
grpoint = TGraph()
grpoint.SetPoint(0,best_coord[0],best_coord[1])

file = TFile("out.root","recreate")

tree = TTree("dataTree","")
tree.ReadFile(namefile,"afb:fL:pvalue")
tree.Write("pvalues")

c = TCanvas()

if obs == "afbfL" :
	var = "fL:afb>>hh("+str(nAfbins)+","+str(minAfb-stepAfb/2)+","+str(maxAfb+stepAfb/2)+","+str(nfLbins)+","+str(minfL-stepfL/2)+","+str(maxfL+stepfL/2)+")"
else :
	var = "afb>>hh("+str(nAfbins)+","+str(minAfb-stepAfb/2)+","+str(maxAfb+stepAfb/2)+")"

#print var
tree.Draw(var,"pvalue","colz")
gr = gPad.GetPrimitive("hh")

if obs != "afbfL" :
	#gr.Draw("HIST C P");
	#c.Print("pvalue_"+obs+"_"+dict[folder[-2:]]+".pdf");
	error = extractLimits(gr,0.68)
	string = 'res  {:5.4}_[{:5.4}]^[{:5.4}]'.format(best_coord[0], error[0]-best_coord[0], error[1]-best_coord[0])
	print string.replace("[","{").replace("]","}")

	gr.SetMarkerStyle(20);
	gr.SetMarkerSize(0.8);
	gr.SetMarkerColor(1);
	gr.SetLineColor(1);
	gr.SetTitle("");
	if obs == "afb" :
		gr.GetXaxis().SetTitle("A_{FB}^{l}");
	elif obs == "afbB" :
		gr.GetXaxis().SetTitle("A_{FB}^{h}");
	elif obs == "fL" :
		gr.GetXaxis().SetTitle("f_{L}");
	gr.GetYaxis().SetTitle("pvalue");
	gr.Draw("HIST C P");
	line = TLine(minAfb,0.32,maxAfb*0.99,0.32);
	line.SetLineColor(kRed);
	line.Draw();

	line2 = TLine(error[0],0,error[0],1);
	line2.SetLineColor(kRed);
	line2.Draw();
	line3 = TLine(error[1],0.,error[1],1);
	line3.SetLineColor(kRed);
	line3.Draw();
	c.Print("pvalue_"+obs+"_"+dict[folder[-2:]]+".pdf");
	sys.exit()

gr.Write("pvalues_grid")
file.Close();

os.system("./collectData.out")
file = TFile("out.root","update")
myconts = file.Get("corrected_grid")
conts = TGraph(myconts)

errAm, errAp, errFm, errFp = 1e6, -1e6, 1e6, -1e6
for i in range(0,conts.GetN()) :
	x, y = Double(0), Double(0)
	conts.GetPoint(i,x,y);
	if ( x < errAm ) :
		errAm = x;
	if ( x > errAp ) :
		errAp = x;
	if ( y < errFm ) :
		errFm = y;
	if ( y > errFp ) :
		errFp = y;

string = 'A_[FB]^[\ell] = {:5.4}_[{:5.4}]^[{:5.4}]'.format(best_coord[0], errAm-best_coord[0], errAp-best_coord[0])
print string.replace("[","{").replace("]","}")
string = 'f_[L] = {:5.4}_[{:5.4}]^[{:5.4}]'.format(best_coord[1], errFm-best_coord[1], errFp-best_coord[1])
print string.replace("[","{").replace("]","}")



nbins = 10
xx = array( 'f' , drange(minAfb*1.05,maxAfb*1.05,nbins) )
yy = array( 'f' , drange(minfL-0.05,maxfL+0.05,nbins) )
htmp = TH2F("tmp","",nbins,xx,nbins,yy)
htmp.GetXaxis().SetTitle("A_{FB}^{l}")
htmp.GetYaxis().SetTitle("f_{L}")
htmp.SetTitle(tdict[folder[-2:]]+" GeV^{2}/#it{c}^{4}")
htmp.Draw()

title = TPaveText(0.45,0.9,0.55,1.,"BRNDC")
title.SetFillStyle(0)
title.SetBorderSize(0)
title.AddText(tdict[folder[-2:]]+" GeV^{2}/#it{c}^{4}")
#title.Draw()

#conts.SetLineColor(1)
#conts.SetLineStyle(1)
#conts.SetLineWidth(1)
gStyle.SetOptStat(0)

xt = array( 'f' , [ -0.75,0.,0.75 ] )
yt = array( 'f' , [ 0.,1.,0. ] )
pline = TPolyLine(3,xt,yt)
pline.SetLineColor(1)
pline.SetLineWidth(1)
pline.Draw("f")

gStyle.SetOptStat(0)
conts.Draw("L SAME")
#gr.Draw("COLZ TEXT")

if best_coord[0] > -(best_coord[1]-1)*3./4. :
	m = -4./3.
	q = 1.
	grpoint.SetPoint(0,(m*best_coord[0]+best_coord[1]-q)/(2*m),(m*best_coord[0]+best_coord[1]+q)/2.)

grpoint.SetMarkerColor(4)
grpoint.SetMarkerSize(1.1)
grpoint.SetMarkerStyle(29)
grpoint.Draw("P same")

lhcbName.Draw()

if(is1D) :
	c.Print("contours_"+dict[folder[-2:]]+"_1D.pdf")
	c.Print("contours_"+dict[folder[-2:]]+"_1D.eps")
else :
	c.Print("contours_"+dict[folder[-2:]]+"_3D.pdf")

file.Close()






