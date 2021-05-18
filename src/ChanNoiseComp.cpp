void ChanNoiseComp(){

  TString fin_name;
  TString fout_name;

  const int fChanNum = 18;
  TH1D *noise_esti1[fChanNum];
  TH1D *noise_esti2[fChanNum];
  TH1D *noise_esti3[fChanNum];
  TH1D *noise_esti4[fChanNum];
  TString name;

  // Read plots from each file
  fin_name = "../../results/ChanNoise_SFGDCubes_GluedFiber_May04.root";
  TFile *fin1 = new TFile(fin_name.Data(),"read");
  fin1->cd();

  for(int i = 0; i < fChanNum; i++){
    name.Form("channel%i_noise",i+1);
    noise_esti1[i] = (TH1D*)fin1->Get(name);
  }

  fin_name = "../../results/ChanNoise_SFGDCubes_GluedFiber_May05.root";
  TFile *fin2 = new TFile(fin_name.Data(),"read");
  fin2->cd();

  for(int i = 0; i < fChanNum; i++){
    name.Form("channel%i_noise",i+1);
    noise_esti2[i] = (TH1D*)fin2->Get(name);
  }

  fin_name = "../../results/ChanNoise_SFGDCubes_GluedFiber_May06.root";
  TFile *fin3 = new TFile(fin_name.Data(),"read");
  fin3->cd();

  for(int i = 0; i < fChanNum; i++){
    name.Form("channel%i_noise",i+1);
    noise_esti3[i] = (TH1D*)fin3->Get(name);
  }

  fin_name = "../../results/ChanNoise_SFGDCubes_GluedFiber_May07.root";
  TFile *fin4 = new TFile(fin_name.Data(),"read");
  fin4->cd();

  for(int i = 0; i < fChanNum; i++){
    name.Form("channel%i_noise",i+1);
    noise_esti4[i] = (TH1D*)fin4->Get(name);
  }

  // Normalize each plot and set some draw options
  int color[4] = {2,4,28,30};
  double scale_factor;
  double upp_limit = 0.15;

  for(int i = 0; i < fChanNum; i++){
    scale_factor = noise_esti1[i]->Integral();
    noise_esti1[i]->Scale(1./scale_factor);
    noise_esti1[i]->GetYaxis()->SetRangeUser(0,upp_limit);
    noise_esti1[i]->SetLineColor(2);

    scale_factor = noise_esti2[i]->Integral();
    noise_esti2[i]->Scale(1./scale_factor);
    noise_esti2[i]->GetYaxis()->SetRangeUser(0,upp_limit);
    noise_esti2[i]->SetLineColor(4);

    scale_factor = noise_esti3[i]->Integral();
    noise_esti3[i]->Scale(1./scale_factor);
    noise_esti3[i]->GetYaxis()->SetRangeUser(0,upp_limit);
    noise_esti3[i]->SetLineColor(28);

    scale_factor = noise_esti4[i]->Integral();
    noise_esti4[i]->Scale(1./scale_factor);
    noise_esti4[i]->GetYaxis()->SetRangeUser(0,upp_limit);
    noise_esti4[i]->SetLineColor(30);
  } 

  // Draw the plots
  TLegend *lg = new TLegend(0.67,0.5,0.86,0.87) ;
  lg->SetLineWidth(0);
  lg->SetFillStyle(0);
  lg->AddEntry(noise_esti1[0],"May 4th","l");
  lg->AddEntry(noise_esti2[0],"May 5th","l");
  lg->AddEntry(noise_esti3[0],"May 6th","l");
  lg->AddEntry(noise_esti4[0],"May 7th","l");

  gStyle->SetOptStat(0);

  // MPPC channel face 1 + 3
  TCanvas *c1 = new TCanvas("noisecomp_xz","noisecomp_xz",1200,1200);
  c1->Divide(3,3);

  c1->cd(1);
  noise_esti1[12]->Draw("hist");
  noise_esti2[12]->Draw("hist same");
  noise_esti3[12]->Draw("hist same");
  noise_esti4[12]->Draw("hist same");
  lg->Draw("same");

  c1->cd(2);
  noise_esti1[3]->Draw("hist");
  noise_esti2[3]->Draw("hist same");
  noise_esti3[3]->Draw("hist same");
  noise_esti4[3]->Draw("hist same");
  lg->Draw("same");

  c1->cd(3);
  noise_esti1[13]->Draw("hist");
  noise_esti2[13]->Draw("hist same");
  noise_esti3[13]->Draw("hist same");
  noise_esti4[13]->Draw("hist same");
  lg->Draw("same");

  c1->cd(4);
  noise_esti1[2]->Draw("hist");
  noise_esti2[2]->Draw("hist same");
  noise_esti3[2]->Draw("hist same");
  noise_esti4[2]->Draw("hist same");
  lg->Draw("same");

  c1->cd(5);
  noise_esti1[11]->Draw("hist");
  noise_esti2[11]->Draw("hist same");
  noise_esti3[11]->Draw("hist same");
  noise_esti4[11]->Draw("hist same");
  lg->Draw("same");

  c1->cd(6);
  noise_esti1[1]->Draw("hist");
  noise_esti2[1]->Draw("hist same");
  noise_esti3[1]->Draw("hist same");
  noise_esti4[1]->Draw("hist same");
  lg->Draw("same");

  c1->cd(7);
  noise_esti1[9]->Draw("hist");
  noise_esti2[9]->Draw("hist same");
  noise_esti3[9]->Draw("hist same");
  noise_esti4[9]->Draw("hist same");
  lg->Draw("same");

  c1->cd(8);
  noise_esti1[0]->Draw("hist");
  noise_esti2[0]->Draw("hist same");
  noise_esti3[0]->Draw("hist same");
  noise_esti4[0]->Draw("hist same");
  lg->Draw("same");

  c1->cd(9);
  noise_esti1[10]->Draw("hist");
  noise_esti2[10]->Draw("hist same");
  noise_esti3[10]->Draw("hist same");
  noise_esti4[10]->Draw("hist same");
  lg->Draw("same");

  c1->Update();

  // MPPC channel face 2 + 4
  TCanvas *c2 = new TCanvas("noisecomp_yz","noisecomp_yz",1200,1200);
  c2->Divide(3,3);

  c2->cd(1);
  noise_esti1[8]->Draw("hist");
  noise_esti2[8]->Draw("hist same");
  noise_esti3[8]->Draw("hist same");
  noise_esti4[8]->Draw("hist same");
  lg->Draw("same");

  c2->cd(2);
  noise_esti1[17]->Draw("hist");
  noise_esti2[17]->Draw("hist same");
  noise_esti3[17]->Draw("hist same");
  noise_esti4[17]->Draw("hist same");
  lg->Draw("same");

  c2->cd(3);
  noise_esti1[7]->Draw("hist");
  noise_esti2[7]->Draw("hist same");
  noise_esti3[7]->Draw("hist same");
  noise_esti4[7]->Draw("hist same");
  lg->Draw("same");

  c2->cd(4);
  noise_esti1[15]->Draw("hist");
  noise_esti2[15]->Draw("hist same");
  noise_esti3[15]->Draw("hist same");
  noise_esti4[15]->Draw("hist same");
  lg->Draw("same");

  c2->cd(5);
  noise_esti1[6]->Draw("hist");
  noise_esti2[6]->Draw("hist same");
  noise_esti3[6]->Draw("hist same");
  noise_esti4[6]->Draw("hist same");
  lg->Draw("same");

  c2->cd(6);
  noise_esti1[16]->Draw("hist");
  noise_esti2[16]->Draw("hist same");
  noise_esti3[16]->Draw("hist same");
  noise_esti4[16]->Draw("hist same");
  lg->Draw("same");

  c2->cd(7);
  noise_esti1[5]->Draw("hist");
  noise_esti2[5]->Draw("hist same");
  noise_esti3[5]->Draw("hist same");
  noise_esti4[5]->Draw("hist same");
  lg->Draw("same");

  c2->cd(8);
  noise_esti1[14]->Draw("hist");
  noise_esti2[14]->Draw("hist same");
  noise_esti3[14]->Draw("hist same");
  noise_esti4[14]->Draw("hist same");
  lg->Draw("same");

  c2->cd(9);
  noise_esti1[4]->Draw("hist");
  noise_esti2[4]->Draw("hist same");
  noise_esti3[4]->Draw("hist same");
  noise_esti4[4]->Draw("hist same");
  lg->Draw("same");

  c2->Update();

  // Save the plots into output file
  TString prefix = "../../../plots/scintillator_cube/";
  TString type;
  TString suffix;

  type = "sfgdcubes_gluedfiber/";

  suffix = prefix + type + "noisecomp_xz.png";
  c1->SaveAs(suffix);

  suffix = prefix + type + "noisecomp_yz.png";
  c2->SaveAs(suffix);

  fout_name = "../../results/NoiseComp_SFGDCubes_GluedFiber.root"; 

  TFile *fout = new TFile(fout_name.Data(),"recreate");
  fout->cd();

  c1->Write();
  c2->Write();

  fout->Close();
  fin1->Close();
  fin2->Close();
  fin3->Close();
  fin4->Close();

}
