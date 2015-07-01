from ROOT import *

gROOT.ProcessLine(".x ~/work/lhcbStyle.C");

from ROOT import TPaveText, gStyle

lhcbName = TPaveText(gStyle.GetPadLeftMargin() + 0.05,
                        0.83 - gStyle.GetPadTopMargin(),
                        gStyle.GetPadLeftMargin() + 0.20,
                        0.95 - gStyle.GetPadTopMargin(),
                        "BRNDC");

lhcbName.AddText("LHCb");
lhcbName.SetFillColor(0);
lhcbName.SetTextAlign(12);
lhcbName.SetBorderSize(0);

line = TLine(0,0,21,0)
line.SetLineStyle(8)
line.SetLineColor(16)


def getPredictions(name = "Lambdab_Lambda_mu_mu_AFB_binned.dat"):
	
	lines = open("../Lambdab_Lambda_mu_mu_SM_predictions/"+name).readlines()

	i=0
	gr = TGraphErrors()
	for l in lines :
		el = map(float,l.split())
		q2 = (el[1]+el[0])/2.
		q2_err = (el[1]-el[0])/2.
		gr.SetPoint(i,q2,el[2])
		gr.SetPointError(i,q2_err,el[3])
		#print q2, q2_err, '{:.4} +- {:.4}'.format(el[2],el[3])
		i+=1

	return gr;


afb_sys = { 1.05 : 0.034, 11.75 : 0.056, 15.5 : 0.034, 17.0 : 0.043, 19.0 : 0.035, 17.5 : 0.034 }
afbB_sys = { 1.05 : 0.151, 11.75 : 0.043, 15.5 : 0.032, 17.0 : 0.031, 19.0 : 0.032, 17.5 : 0.026 }
fL_sys = { 1.05 : 0.083, 11.75 : 0.058, 15.5 : 0.048, 17.0 : 0.046, 19.0 : 0.042, 17.5 : 0.032 }

f_afb = open("res_afb").readlines()
f_afbB = open("res_afbB").readlines()
f_fL = open("res_fL").readlines()

q2 = []
q2err = []
afb = []

for l in f_afb :
	l = l.split()

	## get q2 bin
	q2lim = l[0].replace("-"," ").split()
	q2.append( ( float(q2lim[1])+float(q2lim[0]) ) / 2. )
	q2err.append( ( float(q2lim[1])-float(q2lim[0]) ) / 2. )

	## get value and errors
	vals = l[2].replace("_{"," ").replace("}^{"," ").replace("}","").split()
	afb.append( { "val" : float(vals[0]), "errm" : float(vals[1]), "errp" : float(vals[2]) } )

afbB = []

for l in f_afbB :
	l = l.split()

	## get value and errors
	vals = l[2].replace("_{"," ").replace("}^{"," ").replace("}","").split()
	afbB.append( { "val" : float(vals[0]), "errm" : float(vals[1]), "errp" : float(vals[2]) } )

fL = []

for l in f_fL :
	l = l.split()

	## get value and errors
	vals = l[2].replace("_{"," ").replace("}^{"," ").replace("}","").split()
	fL.append( { "val" : float(vals[0]), "errm" : float(vals[1]), "errp" : float(vals[2]) } )


gr_afb = TGraphAsymmErrors()
gr_afbB = TGraphAsymmErrors()
gr_afb_int = TGraphAsymmErrors()
gr_afbB_int = TGraphAsymmErrors()

p = 0
for q,qe,a,aB in zip(q2,q2err,afb,afbB) :

	if q == 17.5 : continue

	gr_afb.SetPoint(p,q,a["val"])
	gr_afbB.SetPoint(p,q,aB["val"])
	
	gr_afb.SetPointError(p,qe,qe,
			TMath.Sqrt( TMath.Power(a["errm"],2) + TMath.Power(afb_sys[q],2)), TMath.Sqrt( TMath.Power(a["errp"],2) + TMath.Power(afb_sys[q],2)) )
	gr_afbB.SetPointError(p,qe,qe,
			TMath.Sqrt( TMath.Power(aB["errm"],2) + TMath.Power(afbB_sys[q],2)), TMath.Sqrt( TMath.Power(aB["errp"],2) + TMath.Power(afbB_sys[q],2)) )
	p+=1

gr_afb_int.SetPoint(0,17.5,afb[5]["val"])
gr_afbB_int.SetPoint(0,17.5,afbB[5]["val"])
	
gr_afb_int.SetPointError(0,qe,qe,
		TMath.Sqrt( TMath.Power(afb[5]["errm"],2) + TMath.Power(afb_sys[17.5],2)), TMath.Sqrt( TMath.Power(afbB[5]["errp"],2) + TMath.Power(afbB_sys[17.5],2)) )
gr_afbB_int.SetPointError(0,qe,qe,
		TMath.Sqrt( TMath.Power(afbB[5]["errm"],2) + TMath.Power(afbB_sys[17.5],2)), TMath.Sqrt( TMath.Power(afbB[5]["errp"],2) + TMath.Power(afbB_sys[17.5],2)) )
	

c = TCanvas()
gr_afb.SetMinimum(-1)
gr_afb.SetMaximum(1)
gr_afb.GetXaxis().SetTitle("#it{q}^{2} [GeV^{2}/#it{c}^{4}]")
gr_afb.GetYaxis().SetTitle("A_{FB}^{#it{l}}")
gr_afbB.SetMinimum(-0.5)
gr_afbB.SetMaximum(0.5)
gr_afbB.GetXaxis().SetTitle("#it{q}^{2} [GeV^{2}/#it{c}^{4}]")
gr_afbB.GetYaxis().SetTitle("A_{FB}^{#it{h}}")

gr_afb.SetMarkerStyle(21)
gr_afb.SetMarkerSize(1)
gr_afbB.SetMarkerStyle(21)
gr_afbB.SetMarkerSize(1)

gr_afb_int.SetMarkerStyle(22)
gr_afb_int.SetMarkerSize(1.3)
gr_afb_int.SetMarkerColor(2)
gr_afbB_int.SetMarkerStyle(22)
gr_afbB_int.SetMarkerSize(1.3)
gr_afbB_int.SetMarkerColor(2)


Theo_gr = getPredictions()
Theo_gr.SetMarkerSize(0)
Theo_gr.SetMarkerColor(4)
Theo_gr.SetLineWidth(0)
Theo_gr.SetLineColor(7)
Theo_gr.SetFillColor(7)
Theo_gr_tmp = Theo_gr.Clone()
for p in xrange(Theo_gr.GetN()) :
	Theo_gr_tmp.SetPointError(p,Theo_gr_tmp.GetErrorX(p),0)
Theo_gr_tmp.SetMarkerSize(0)
Theo_gr_tmp.SetLineColor(4)
Theo_gr_tmp.SetLineStyle(kDashed)

TheoB_gr = getPredictions("Lambdab_Lambda_mu_mu_AFBh_binned.dat")
TheoB_gr.SetMarkerSize(0)
TheoB_gr.SetMarkerColor(4)
TheoB_gr.SetLineWidth(0)
TheoB_gr.SetLineColor(7)
TheoB_gr.SetFillColor(7)
TheoB_gr_tmp = TheoB_gr.Clone()
for p in xrange(Theo_gr.GetN()) :
	TheoB_gr_tmp.SetPointError(p,TheoB_gr_tmp.GetErrorX(p),0)
TheoB_gr_tmp.SetMarkerSize(0)
TheoB_gr_tmp.SetLineColor(4)
TheoB_gr_tmp.SetLineStyle(kDashed)

leg_abs = TLegend(0.65,0.65,0.9,0.9)
leg_abs.AddEntry(Theo_gr,"SM prediction","F")
leg_abs.AddEntry(gr_afb,"Data","PL")
leg_abs.SetFillColor(0)

predict = True
predstr = ""
if predict : 
	predstr = "_pred"

gr_afb.Draw("AP")
line.Draw()
if predict : Theo_gr.Draw("P SAME 2")
if predict : Theo_gr_tmp.Draw("P SAME")
lhcbName.Draw()
if predict : leg_abs.Draw()
gr_afb.Draw("P SAME")
c.Print("Afb_vs_q2"+predstr+".pdf")
c.Print("Afb_vs_q2"+predstr+".eps")
c.Print("Afb_vs_q2"+predstr+".png")
c.Print("Afb_vs_q2"+predstr+".C")

gr_afbB.Draw("AP")
line.Draw()
if predict : TheoB_gr.Draw("P SAME 2")
if predict : TheoB_gr_tmp.Draw("P SAME")
lhcbName.Draw()
if predict : leg_abs.Draw()
gr_afbB.Draw("P SAME")
c.Print("AfbB_vs_q2"+predstr+".pdf")
c.Print("AfbB_vs_q2"+predstr+".eps")
c.Print("AfbB_vs_q2"+predstr+".png")
c.Print("AfbB_vs_q2"+predstr+".C")

gr_afb.Draw("AP")
line.Draw()
if predict : Theo_gr.Draw("P SAME 2")
if predict : Theo_gr_tmp.Draw("P SAME")
gr_afb_int.Draw("P SAME")
lhcbName.Draw()
if predict : leg_abs.Draw()
gr_afb.Draw("P SAME")
gr_afb_int.Draw("P SAME")
c.Print("Afb_vs_q2_int"+predstr+".pdf")
c.Print("Afb_vs_q2_int"+predstr+".eps")
c.Print("Afb_vs_q2_int"+predstr+".png")
c.Print("Afb_vs_q2_int"+predstr+".C")


gr_afbB.Draw("AP")
line.Draw()
if predict : TheoB_gr.Draw("P SAME 2")
if predict : TheoB_gr_tmp.Draw("P SAME")
gr_afbB_int.Draw("P SAME")
lhcbName.Draw()
if predict : leg_abs.Draw()
gr_afbB.Draw("P SAME")
gr_afbB_int.Draw("P SAME")
c.Print("AfbB_vs_q2_int"+predstr+".pdf")
c.Print("AfbB_vs_q2_int"+predstr+".eps")
c.Print("AfbB_vs_q2_int"+predstr+".png")
c.Print("AfbB_vs_q2_int"+predstr+".C")


for q,qe,a,f in zip(q2,q2err,afb,fL) :

	q2str = '{:.1f}--{:.1f} & '.format(q-qe,q+qe)
	afbstr = '${:.2f} \\; ^[+{:.2f}]_[{:.2f}] \\; \\pm \\; {:.2f}$'.format(a["val"],a["errp"],a["errm"],afb_sys[q])
	flstr = '${:.2f} \\; ^[+{:.2f}]_[{:.2f}] \\; \\pm \\; {:.2f}$'.format(f["val"],f["errp"],f["errm"],fL_sys[q])
	print q2str+afbstr.replace("[","{").replace("]","}")+"  &  ",flstr.replace("[","{").replace("]","}")+" \\\\"

print
print

for q,qe,aB in zip(q2,q2err,afbB) :

	q2str = '{:.1f}--{:.1f} & '.format(q-qe,q+qe)
	afbstr = '${:.2f} \\; ^[+{:.2f}]_[{:.2f}] \\; \\pm \\; {:.2f}$'.format(aB["val"],aB["errp"],aB["errm"],afbB_sys[q])
	print q2str+afbstr.replace("[","{").replace("]","}")+" \\\\"




chi2 = 0
chi2_low = 0
chi2_high = 0
ndf = 0

for i in xrange(gr_afb.GetN()) :
	x, v = Double(0), Double(0)
	ep, em = Double(0), Double(0)
	gr_afb.GetPoint(i,x,v)
	errlow = gr_afb.GetErrorYlow(i)
	errhigh = gr_afb.GetErrorYhigh(i)
			
	for j in xrange(Theo_gr.GetN()) :
		xt, vt = Double(0), Double(0)
		ept, emt = Double(0), Double(0)
		Theo_gr.GetPoint(j,xt,vt)
		errlow_t = Theo_gr.GetErrorYlow(j)
		errhigh_t = Theo_gr.GetErrorYhigh(j)

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


