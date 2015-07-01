/*
 * Plot type definies what output format:
 *   0 = eps
 *   1 = pdf
 *   2 = png
 *   3 = gif
 */

void purityPlots(std::string filename, bool zeroIt = false, 
   int plotType = 0 )
{

//  gSystem->CompileMacro("~/work/NeuroBayes/cpptools//analysis.C");

  gStyle->SetOptStat(0);

  std::string extension;
  switch (plotType)
  {
    case (0): extension="eps";
              break;
    case (1): extension="pdf";
              break;
    case (2): extension="png";
              break;
    case (3): extension="gif";
              break;
    default: extension="eps";
  }

  TH1F** PurEffGraph = new TH1F*[5*1];
  
  TFile file(filename.c_str());
  TH1F* h101,*h201;
  if ( zeroIt )
  {
    h101=(TH1F*)file.Get("h9102");
    h201=(TH1F*)file.Get("h9202");
  }
  else
  {
    h101=(TH1F*)file.Get("h101");
    h201=(TH1F*)file.Get("h201");
  }
  
  // Calculate purity-efficiency plots
  PurityEfficiency(h101, h201, &(PurEffGraph[0]),1,0,false);

  TCanvas c1;
  c1.SetGrid();
  c1.Draw();

  // purity vs signal eff

  PurEffGraph[2]->SetMinimum(0);
  PurEffGraph[2]->SetMaximum(1);
  PurEffGraph[2]->GetYaxis()->SetTitle("Signal purity");
  PurEffGraph[2]->GetYaxis()->SetTitleOffset(1.0);
  PurEffGraph[2]->GetXaxis()->SetTitle("Signal efficiency");
  PurEffGraph[2]->GetXaxis()->SetTitleOffset(1.0);
  PurEffGraph[2]->GetXaxis()->SetNdivisions(510);
  PurEffGraph[2]->Draw();
  PurEffGraph[3]->Draw("same");
  c1.Update();
  std::string name="signal_eff_pur.";
  c1.Print((name+extension).c_str());

  c1.Clear();

  // Signal eff vs eff

  PurEffGraph[4]->SetMinimum(0);
  PurEffGraph[4]->SetMaximum(1);
  PurEffGraph[4]->GetYaxis()->SetTitle("Signal efficiency");
  PurEffGraph[4]->GetYaxis()->SetTitleOffset(1.0);
  PurEffGraph[4]->GetXaxis()->SetTitle("Efficiency");
  PurEffGraph[4]->GetXaxis()->SetTitleOffset(1.0);
  PurEffGraph[4]->GetXaxis()->SetNdivisions(510);
  PurEffGraph[4]->Draw("l");
  TLine l1(0,0,1,1);
  l1.SetLineColor(4);
  l1.Draw();
  int nBins = PurEffGraph[2]->GetNbinsX(); 
  TLine l2(0,0,PurEffGraph[2]->GetBinContent(nBins),1);
  l2.SetLineColor(2);
  l2.Draw();

  c1.Update();
  name="signal_eff_eff.";
  c1.Print((name+extension).c_str());

  // Purity vs NNout

  h101->Sumw2();
  h201->Sumw2();

  TH1F* hsum = (TH1F*) h201->Clone("hsum");
  hsum->Add(h101);
  TH1F* hpur =(TH1F*) h201->Clone("hpur");
  plotPurity(hpur,hsum,1,1);

  hpur->GetXaxis()->SetTitleFont(42);
  hpur->GetXaxis()->SetLabelFont(42);
  hpur->GetYaxis()->SetTitleFont(42);
  hpur->GetYaxis()->SetLabelFont(42);
  hpur->GetXaxis()->SetTitleSize(0.055);
  hpur->GetXaxis()->SetLabelSize(0.05);
  hpur->GetYaxis()->SetTitleSize(0.055);
  hpur->GetYaxis()->SetLabelSize(0.05);
  hpur->GetYaxis()->SetTitleOffset(1.4);
  hpur->GetXaxis()->SetTitleOffset(1.0);
  hpur->GetXaxis()->SetNdivisions(505);

  c1.Update();
  name="purity_NN.";
  c1.Print((name+extension).c_str());

  // Lets add also NN output distribution to have all in one go
  int bkgMax=h101->GetMaximum();
  int sigMax=h201->GetMaximum();
 
  if (sigMax>bkgMax)
    h101->SetMaximum(sigMax);

  h101->SetTitle("");
  h101->GetXaxis()->SetTitle("Network output");
  h101->GetYaxis()->SetTitle("");
  h101->GetXaxis()->SetTitleFont(42);
  h101->GetXaxis()->SetLabelFont(42);
  h101->GetYaxis()->SetTitleFont(42);
  h101->GetYaxis()->SetLabelFont(42);
  h101->GetXaxis()->SetTitleSize(0.055);
  h101->GetXaxis()->SetLabelSize(0.05);
  h101->GetYaxis()->SetTitleSize(0.055);
  h101->GetYaxis()->SetLabelSize(0.05);
  h101->GetYaxis()->SetTitleOffset(1.4);
  h101->GetXaxis()->SetTitleOffset(1.0);
  h101->GetXaxis()->SetNdivisions(505);
  h101->SetLineWidth(2);
  h101->Draw("hist");
  h201->SetLineColor(2);
  h201->SetLineWidth(2);
  h201->Draw("histsame");
  c1.Update();
  name="NN_dist.";
  c1.Print((name+extension).c_str());

  // Lets try to get correlation matrix to eps file
  c1.SetLeftMargin(0.08);
  c1.SetRightMargin(0.15);
  c1.SetGrid(0,0);
  correlations(&file,1);
  TH1F* h200=(TH1F*)gROOT->FindObject("h2000");
  h200->GetXaxis()->SetLabelFont(42);
  h200->GetXaxis()->SetLabelSize(0.05);
  h200->GetYaxis()->SetLabelFont(42);
  h200->GetYaxis()->SetLabelSize(0.05);
  h200->GetZaxis()->SetLabelFont(42);
  h200->GetZaxis()->SetLabelSize(0.04);
  c1.Clear();
  h200->Draw("colz");
  
  c1.Update();
  name="correlation.";
  c1.Print((name+extension).c_str());

  file.Close();
  
}

