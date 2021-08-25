void ana3Dprint3x3x1(){


   TCanvas *Tlva = new TCanvas("Tlva","Tlva",500,500);
   Tlva->SetGrid();
   Tlva->DrawFrame(0,0,1,1);
   const char *longstring = "K_{S}... K^{*0}... #frac{2s}{#pi#alpha^{2}} #frac{d#sigma}{dcos#theta} (e^{+}e^{-} #rightarrow f#bar{f} ) = #left| #frac{1}{1 - #Delta#alpha} #right|^{2} (1+cos^{2}#theta)";
 
  
  TGraph *gLYvsDist = new TGraph(); // Distance from MPPC vs Landau MPV (ADC)

  // // Measurement #1 (w/ Teflon white) 
  // gLYvsDist->SetPoint(gLYvsDist->GetN(),1.,887.);   
  // gLYvsDist->SetPoint(gLYvsDist->GetN(),1.5,876); 
  // gLYvsDist->SetPoint(gLYvsDist->GetN(),2.,879);  
  // gLYvsDist->SetPoint(gLYvsDist->GetN(),2.5,857.);
  // gLYvsDist->SetPoint(gLYvsDist->GetN(),3.,838.);  
  // gLYvsDist->SetPoint(gLYvsDist->GetN(),3.5,833.); 
  // gLYvsDist->SetPoint(gLYvsDist->GetN(),4.,830);  
  // gLYvsDist->SetPoint(gLYvsDist->GetN(),4.5,824.);   

  // Measurement #2 (w/ Teflon white - MPV) 
  //gLYvsDist->SetPoint(gLYvsDist->GetN(),1.,);   
  //gLYvsDist->SetPoint(gLYvsDist->GetN(),1.5,814); 
  //gLYvsDist->SetPoint(gLYvsDist->GetN(),2.,792);  
  //gLYvsDist->SetPoint(gLYvsDist->GetN(),2.5,766);  
  //gLYvsDist->SetPoint(gLYvsDist->GetN(),3.,742); 
  //gLYvsDist->SetPoint(gLYvsDist->GetN(),3.5,727);  
  //gLYvsDist->SetPoint(gLYvsDist->GetN(),4.,715);  
  //gLYvsDist->SetPoint(gLYvsDist->GetN(),4.5,697);    

  // Measurement #2 (w/ Teflon white - Mean) 
  //gLYvsDist->SetPoint(gLYvsDist->GetN(),1.,);  
  //gLYvsDist->SetPoint(gLYvsDist->GetN(),1.5,972); 
  //gLYvsDist->SetPoint(gLYvsDist->GetN(),2.,914); 
  //gLYvsDist->SetPoint(gLYvsDist->GetN(),2.5,866); 
  //gLYvsDist->SetPoint(gLYvsDist->GetN(),3.,829); 
  //gLYvsDist->SetPoint(gLYvsDist->GetN(),3.5,797); 
  //gLYvsDist->SetPoint(gLYvsDist->GetN(),4.,767); 
  //gLYvsDist->SetPoint(gLYvsDist->GetN(),4.5,750);   

  // Measurement 2021-7-22
  gLYvsDist->SetPoint(0,4.5,415.47);
  gLYvsDist->SetPoint(1,4.0,421.13);
  gLYvsDist->SetPoint(2,3.5,431.40);
  gLYvsDist->SetPoint(3,3.0,443.81);
  gLYvsDist->SetPoint(4,2.5,453.50);
  //gLYvsDist->SetPoint(5,2.0,453.50); // Not used same as 2.5 measurement
  gLYvsDist->SetPoint(5,1.5,477.08);

  gLYvsDist->Fit("expo");

  double slope = gLYvsDist->GetFunction("expo")->GetParameter(1);
  cout << "slope = " << slope << endl;
  double attleng = TMath::Abs(1/slope);
  cout << "attenuation length = " << attleng << " cm" << endl;
  

  TCanvas *cLYvsDist = new TCanvas("cLYvsDist","cLYvsDist");
  gLYvsDist->SetMarkerSize(1);
  gLYvsDist->SetMarkerStyle(21);
  gLYvsDist->GetXaxis()->SetTitle("Distance from SiPM (cm)");
  gLYvsDist->GetYaxis()->SetTitle("Light Yield (ADC)");
  gLYvsDist->GetXaxis()->SetLabelSize(.045);
  gLYvsDist->GetXaxis()->SetLabelSize(.045);
  gLYvsDist->GetXaxis()->SetTitleSize(.05);
  gLYvsDist->GetYaxis()->SetTitleSize(.05);
  gLYvsDist->GetXaxis()->SetTitleOffset(.95);
  gLYvsDist->GetYaxis()->SetTitleOffset(.95);
  gLYvsDist->Draw("ap");


  //TLatex latex;
  //latex.SetTextSize(0.025);
  //latex.SetTextAlign(13);  //align at top
  //latex.DrawLatex(.2,.9,"K_{S}");
  //latex.DrawLatex(.3,.9,"K^{*0}");
  //latex.DrawLatex(.2,.8,longstring);

  //string text = "Attenuation Length = " + std::to_string(attleng);;
  char ctext[200];
  sprintf(ctext,"Attenuation Length = %.1f cm",attleng);
  
  TPaveText *pt = new TPaveText(0.5,0.7,0.88,0.85,"NDC");
  pt->SetBorderSize(1);
  pt->SetTextSize(0.04);
  //pt->SetFillColor(0);
  pt->SetFillColor(19);
  pt->SetTextAlign(12);
  pt->SetFillStyle(12);
  pt->AddText(ctext);
  pt->Draw();

  
  return;
}
