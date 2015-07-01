from ROOT import *
from string import *
import sys

gROOT.ProcessLine(".x ~/work/lhcbStyle.C");

base = "/afs/cern.ch/work/p/pluca/Lmumu/weighted/"
fBuKst = TFile(base+"Bu2Kstmumu_MC12_NBweighted.root")
tBuKst = fBuKst.Get("tree")
fBdKS = TFile(base+"Bd2JpsiKS_MC12_NBweighted.root")
tBdKS = fBdKS.Get("tree")
fJpsiL = TFile(base+"Lb2JpsiL_MC_Pythia8_NBweighted.root")
tJpsiL = fJpsiL.Get("tree")
fLmumu = TFile(base+"Lb2Lmumu_MC_Pythia8_NBweighted.root")
tLmumu = fLmumu.Get("tree")

c = TCanvas()

tBuKst.Draw("Lb_MassConsLambda_M>>h(50,5400,6000)")
h = gPad.GetPrimitive("h")
h.GetXaxis().SetTitle("m(#Lambda#mu#mu) [MeV/#it{c}^{2}]")
h.GetYaxis().SetTitle("A.U.")
h.SetTitle("")

c.Print("Bu2Kstplus_mass.pdf")

#tBdKS.Draw("Lb_MM>>h(50,5000,7000)")

#c.Print("Bd2JpsiKS_mass.pdf")

tLmumu.Draw("TMath::Power(J_psi_1S_MM/1000,2)>>h(50,0,20)")
h = gPad.GetPrimitive("h")
h.GetXaxis().SetTitle("q^{2} [GeV^{2}/#it{c}^{4}]")
h.GetYaxis().SetTitle("A.U.")
h.SetTitle("")
h.Draw()

c.Print("Q2_beforemodel.pdf")

tLmumu.Draw("TMath::Power(J_psi_1S_MM/1000,2)>>h2(50,0,20)","physRate_pol0")
h2 = gPad.GetPrimitive("h2")
h2.GetXaxis().SetTitle("q^{2} [GeV^{2}/#it{c}^{4}]")
h2.GetYaxis().SetTitle("A.U.")
h2.SetTitle("")
h2.Draw()

c.Print("Q2_aftermodel.pdf")


tJpsiL.Draw("Lb_MassConsLambda_M>>h(50,5400,6000)","TMath::Power(J_psi_1S_MM/1000,2) < 8")
h = gPad.GetPrimitive("h")
h.GetXaxis().SetTitle("m(#Lambda#mu#mu) [MeV/#it{c}^{2}]")
h.GetYaxis().SetTitle("A.U.")
h.SetTitle("")
h.Draw()

c.Print("JpsiL_leakage_mass.pdf")



tJpsiL.Draw("Lb_MassConsJpsiAndLambda_M>>h(100,5500,5700)")
h = gPad.GetPrimitive("h")
h.GetXaxis().SetTitle("m(#Lambda#mu#mu) [MeV/#it{c}^{2}]")
h.GetYaxis().SetTitle("A.U.")
h.SetLineColor(1)
h.SetTitle("")

tJpsiL.Draw("Lb_Fit_M>>h2(100,5500,5700)")
h2 = gPad.GetPrimitive("h2")
h2.SetLineColor(2)
h2.SetTitle("")

tJpsiL.Draw("Lb_M>>h3(100,5500,5700)")
h3 = gPad.GetPrimitive("h3")
h3.SetLineColor(4)
h3.SetTitle("")

h.Draw()
#h2.Draw("same")
h3.Draw("same")

leg = TLegend(0.2,0.7,0.5,0.9)
leg.SetFillColor(0)
leg.AddEntry(h3,"Not contrained","l")
#leg.AddEntry(h2,"Basic contraints","l")
leg.AddEntry(h,"Contrained","l")
leg.Draw()

c.Print("DTF_performance.pdf")




