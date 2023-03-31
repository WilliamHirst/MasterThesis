{

  //Set ATLAS style
  gROOT->LoadMacro("utils.h");
  gROOT->SetStyle("ATLAS");
  atlasStyle->SetTitleSize(0.06,"Y");
  gROOT->ForceStyle();

  Int_t Nmasses;

  ifstream input("limits.txt");

  input >> Nmasses;
  cout << "Number of mass points: " << Nmasses << endl;

  Double_t* mass = new Double_t[Nmasses];
  Double_t* xsec = new Double_t[Nmasses];
  Double_t* exclusion = new Double_t[Nmasses];
  Double_t* exclusionM2S = new Double_t[Nmasses];
  Double_t* exclusionM1S = new Double_t[Nmasses];
  Double_t* exclusionExp = new Double_t[Nmasses];
  Double_t* exclusionP1S = new Double_t[Nmasses];
  Double_t* exclusionP2S = new Double_t[Nmasses];

  for (Int_t i = 0; i < Nmasses; i++) {
    
    input >> mass[i];
    input >> xsec[i];
    input >> exclusion[i];
    input >> exclusionM2S[i];
    input >> exclusionM1S[i];
    input >> exclusionExp[i];
    input >> exclusionP1S[i];
    input >> exclusionP2S[i];
  }

  TGraph* xsecG = new TGraph(Nmasses,mass,xsec);
  TGraph* exclusionG = new TGraph(Nmasses,mass,exclusion);
  TGraph* exclusionM2SG = new TGraph(Nmasses,mass,exclusionM2S);
  TGraph* exclusionM1SG = new TGraph(Nmasses,mass,exclusionM1S);
  TGraph* exclusionExpG = new TGraph(Nmasses,mass,exclusionExp);
  TGraph* exclusionP1SG = new TGraph(Nmasses,mass,exclusionP1S);
  TGraph* exclusionP2SG = new TGraph(Nmasses,mass,exclusionP2S);

  exclusionM2SG->SetLineColor(5);
  exclusionM1SG->SetLineColor(3);
  exclusionExpG->SetLineStyle(2);
  exclusionP2SG->SetLineColor(5);
  exclusionP1SG->SetLineColor(3);

  exclusionG->SetMarkerStyle(8);

  TGraph* fillM2S = new TGraph(2*Nmasses);
  for (Int_t i = 0; i < Nmasses; i++) {
    fillM2S->SetPoint(i,mass[i],exclusionM1S[i]);
    fillM2S->SetPoint(Nmasses+i,mass[Nmasses-i-1],exclusionM2S[Nmasses-i-1]);
  }
  fillM2S->SetFillColor(5);

  TGraph* fill1S = new TGraph(2*Nmasses);
  for (Int_t i = 0; i < Nmasses; i++) {
    fill1S->SetPoint(i,mass[i],exclusionP1S[i]);
    fill1S->SetPoint(Nmasses+i,mass[Nmasses-i-1],exclusionM1S[Nmasses-i-1]);
  }
  fill1S->SetFillColor(3);

  TGraph* fillP2S = new TGraph(2*Nmasses);
  for (Int_t i = 0; i < Nmasses; i++) {
    fillP2S->SetPoint(i,mass[i],exclusionP2S[i]);
    fillP2S->SetPoint(Nmasses+i,mass[Nmasses-i-1],exclusionP1S[Nmasses-i-1]);
  }
  fillP2S->SetFillColor(5);

  xsecG->GetHistogram()->SetYTitle("#sigma [fb]");
  xsecG->GetHistogram()->SetXTitle("W' mass [GeV]");
  
  xsecG->GetHistogram()->GetYaxis()->SetRangeUser(5.0e-2,5.0e3);
  
  TCanvas* c = new TCanvas;
  c->SetLogy();

  xsecG->Draw("AC*");

  fillM2S->Draw("f");
  fill1S->Draw("f");
  fillP2S->Draw("f");

  exclusionM2SG->Draw("L same");
  exclusionM1SG->Draw("L same");
  exclusionP1SG->Draw("L same");
  exclusionP2SG->Draw("L same");
  exclusionExpG->Draw("L same");
  exclusionG->Draw("LP same");
  xsecG->Draw("C*");
  
  TLegend* leg = new TLegend(0.533046,0.677966,0.899425,0.900424);
  leg->AddEntry(xsecG,"#sigma(W' #rightarrow l#nu)","PL");
  leg->AddEntry(exclusionG,"Excluded at 95% CL","PL");
  leg->AddEntry(exclusionExpG,"Expected limit","L");
  leg->AddEntry(fill1S,"#pm1#sigma","F");
  leg->AddEntry(fillP2S,"#pm2#sigma","F");
  leg->Draw();

  c->SaveAs("./limitplot.png")
}
