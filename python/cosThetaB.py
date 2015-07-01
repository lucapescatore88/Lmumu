from ROOT import *

c = TCanvas();
DD = TFile("MC_reweightDD.root");
hDataDD = DD.Get("cosThetaB_DD_data");
hMCDD = DD.Get("cosThetaB_DD_MC");
LL = TFile("MC_reweightLL.root");
hDataLL = LL.Get("cosThetaB_LL_data");
hMCLL = LL.Get("cosThetaB_LL_MC");

leg2 = TLegend(0.7,0.7,0.99,0.99);
leg2.AddEntry(hDataLL,"Data LL","P");
leg2.AddEntry(hDataDD,"Data DD","P");
leg2.AddEntry(hMCLL,"MC LL","P");
leg2.AddEntry(hMCDD,"MC DD","P");

c.SetLogy(0);
hDataDD.GetXaxis().SetTitle("cos#theta_{#Lambda}");
hDataLL.Draw();
hDataDD.Draw("same");
hMCDD.Draw("same");
hMCLL.Draw("same");
leg2.Draw("same");
c.Print("costhetaB_distribs.pdf");
c.SetLogy();

dataRatio = hDataDD.Clone("dataRatio");
dataRatio.GetXaxis().SetTitle("cos#theta_{#Lambda}");
dataRatio.Divide(hDataLL);
mcRatio = hMCDD.Clone("mcRatio");
mcRatio.GetXaxis().SetTitle("cos#theta_{#Lambda}");
mcRatio.Divide(hMCLL);

ratio = dataRatio.Clone("ratio");
ratio.Divide(mcRatio);
ratio.GetXaxis().SetTitle("cos#theta_{#Lambda}");
ratio.GetYaxis().SetTitle("(DD / LL)_{data} / (DD / LL)_{MC}");
ratio.Draw();
c.SetLogy(0);
c.Print("DD_Over_LL_MC_Over_Data_cosThetaB.pdf");

leg3 = TLegend(0.7,0.7,0.99,0.99);
leg3.AddEntry(dataRatio,"Data","P");
leg3.AddEntry(mcRatio,"MC","P");

dataRatio.SetMarkerSize(0.8);
dataRatio.SetMarkerColor(1);
dataRatio.SetMarkerStyle(20);
mcRatio.SetMarkerSize(0.8);
mcRatio.SetMarkerColor(2);
mcRatio.SetMarkerStyle(21);
mcRatio.GetYaxis().SetTitle("DD / LL");
mcRatio.Draw();
dataRatio.Draw("same");
leg3.Draw("same");
c.Print("DD_Over_LL_MC_And_Data_cosThetaB.pdf");


DDRatio = hDataDD.Clone("DDRatio");
DDRatio.GetXaxis().SetTitle("cos#theta_{#Lambda}");
DDRatio.Divide(hMCDD);
LLRatio = hDataLL.Clone("LLRatio");
LLRatio.GetXaxis().SetTitle("cos#theta_{#Lambda}");
LLRatio.Divide(hMCLL);

leg4 = TLegend(0.7,0.7,0.99,0.99);
leg4.AddEntry(LLRatio,"LL","P");
leg4.AddEntry(DDRatio,"DD","P");

LLRatio.SetMarkerSize(0.8);
LLRatio.SetMarkerColor(1);
LLRatio.SetMarkerStyle(20);
DDRatio.SetMarkerSize(0.8);
DDRatio.SetMarkerColor(2);
DDRatio.SetMarkerStyle(21);
LLRatio.GetYaxis().SetTitle("data / MC");
LLRatio.Draw();
DDRatio.Draw("same");
leg4.Draw("same");
c.Print("Data_Over_MC_DD_And_LL_cosThetaB.pdf");
