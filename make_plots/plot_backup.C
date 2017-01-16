Double_t myVoight(Double_t *x, Double_t *par)
{
  Float_t xx =x[0];
  Double_t f=0.;
  Double_t a = 1.;
  
  if(xx>0.){
    a=TMath::Sqrt(xx/par[1]);
  }
  f = par[0]*TMath::Voigt(xx-par[1],par[2],par[3]*a,4);
  
  return f;
}

Double_t myFuncBWG(Double_t *x, Double_t *par)
{
  Float_t xx =x[0];
  Double_t f=0.;
  Double_t a = 1.;
  
  if(xx>0.){
    a=TMath::Sqrt(xx/par[1]);
  }
  f = par[0]*TMath::Voigt(xx-par[1],par[2],par[3]*a,4) + par[4] + par[5]*xx + par[6]*xx*xx;// + par[7]*xx*xx*xx ;
  
  return f;
}


void fitMKVoight(TH1 *h33 , Double_t low, Double_t high, Double_t p0, Double_t p1, Double_t p2, Double_t p3, Double_t initialPar, Double_t width, Double_t Gamma, Double_t factor, Int_t draw_opt){
  
  //double nEnt = h33->GetEntries();
  double nEnt = h33->GetMaximum();
  double in_par[8] = {100, initialPar, width , Gamma, p0, p1, p2}; //{A,mean,sigma,Gamma}
  TF1 *fitter = new TF1("fitter",myFuncBWG,low,high,8);
  fitter->SetParameters(&in_par[0]);
  fitter->SetParLimits(0,10,nEnt); fitter->SetParLimits(1,initialPar-width,initialPar+width);
  fitter->SetLineColor(kBlue);
  h33->Fit("fitter","REM+","same");
  //h33->Draw("E");
  TF1 *backFcn = new TF1("backFcn", "pol2",low,high);
  TF1 *signalFcn = new TF1("signalFcn", myVoight,low,high,4);
  signalFcn->SetLineColor(kBlue);
  signalFcn->SetLineWidth(2);
  Double_t par[8];
  fitter->GetParameters(par);
  backFcn->SetParameters(&par[4]);
  backFcn->SetLineStyle(2);
  backFcn->SetLineColor(6);
  backFcn->SetLineWidth(1);
  
  
  signalFcn->SetParameters(par);
  signalFcn->SetLineStyle(2);
  signalFcn->SetLineColor(4);
  signalFcn->SetLineWidth(1);
  // backFcn->Draw("same");
  //signalFcn->Draw("same");
  //Double_t Intg = abs(signalFcn->Integral(par[1]-factor*par[2],par[1]+factor*par[2]));
  ROOT::Math::GSLIntegrator ig(1.E-6,1.E-6,1000);
  ROOT::Math::WrappedTF1 wf(*fitter);
  ig.SetFunction(wf);
  double Intg = ig.Integral(par[1]-factor*par[2],par[1]+factor*par[2]);
  Double_t Intb = abs(backFcn->Integral(par[1]-factor*par[2],par[1]+factor*par[2]));
  
  
  Double_t binw = h33->GetBinWidth(1);
  Int_t yield = Intg/binw;
  Int_t bckgd = Intb/binw;
  //Double_t ratio = double(yield)/TMath::Sqrt(double(bckgd));
  Double_t ratio = double(yield)/double(yield+bckgd);
  
  cout << yield << "\t" << ratio << endl;
  TAxis *x=h33->GetXaxis();
  TAxis *y=h33->GetYaxis();
  
  Double_t startx=x->GetXmin()+0.45*(x->GetXmax()-x->GetXmin());
  /*
   Double_t starty0=0.5*h33->GetMaximum();
   Double_t starty1=0.56*h33->GetMaximum();
   Double_t starty2=0.62*h33->GetMaximum();
   Double_t starty3=0.68*h33->GetMaximum();
   Double_t starty4=0.74*h33->GetMaximum();
   Double_t starty5=0.8*h33->GetMaximum();
   */
  Double_t starty0=0.35*h33->GetMaximum();
  Double_t starty1=0.43*h33->GetMaximum();
  Double_t starty2=0.49*h33->GetMaximum();
  Double_t starty3=0.56*h33->GetMaximum();
  Double_t starty4=0.63*h33->GetMaximum();
  Double_t starty5=0.7*h33->GetMaximum();
  
  
  
  double meanError = fitter->GetParError(1);
  double sigmaError = fitter->GetParError(2);
  if (draw_opt ==1) {
    
    TLatex *sum = new TLatex(startx*0.93, starty5,Form("Yield: %i",yield));
    TLatex *sum12 = new TLatex(startx*0.93, starty4,Form("Background: %i",bckgd));
    TLatex *sum0=new TLatex(startx*0.93, starty3,Form("Range: #pm %2.1f #sigma",factor));
    TLatex *sum2=new TLatex(startx*0.93, starty2,Form("Mean: %4.4f #pm %.4f GeV",par[1], meanError));
    TLatex *sum3=new TLatex(startx*0.93, starty1,Form("#sigma: %5.4f #pm %.4f GeV",par[2], sigmaError));
    //TLatex *ra = new TLatex(startx*0.93, starty0,Form("#frac{S}{#sqrt{B}}= %.1f", ratio));
    TLatex *ra = new TLatex(startx*0.93, starty0,Form("#frac{S}{S+B} = %1.4f", ratio));
    
    
    sum->SetTextSize(0.04);
    sum->SetTextColor(kBlue);//2
    sum->Draw("same");
    sum12->SetTextSize(0.04);
    sum12->SetTextColor(kBlue);//6
    sum12->Draw("same");
    
    ra->SetTextSize(0.04);
    ra->SetTextColor(kBlue);//was black before;
    ra->Draw("same");
    
    
    sum0->SetTextSize(0.04);
    sum0->SetTextColor(kBlue);//2
    sum0->Draw("same");
    sum2->SetTextSize(0.04);
    sum2->SetTextColor(kBlue);//4
    sum2->Draw("same");
    sum3->SetTextSize(0.04);
    sum3->SetTextColor(kBlue);//4
    sum3->Draw("same");
    
    
  }
}



void plot_backup(){
  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  gROOT->ForceStyle();
  TFile *in = TFile::Open("IncLusive_Plots_vmdtest.root");
  
  TH1D *hIVEpEmGam = (TH1D*)in->Get("hIVEpEmGam");
  TH1D *hIVEpEmGam_cut = (TH1D*)in->Get("hIVEpEmGam_cut");
  
  hIVEpEmGam->SetLineColor(kBlack);
  hIVEpEmGam->SetStats(false);
  hIVEpEmGam_cut->SetStats(false);
  hIVEpEmGam->GetXaxis()->SetTitle("M(e^{+}e^{-}#gamma) [GeV]");
  hIVEpEmGam->GetXaxis()->SetTitleSize(0.05);
  hIVEpEmGam->GetXaxis()->SetTitleOffset(0.8);
  hIVEpEmGam->GetYaxis()->SetTitle("Counts / 35.5 MeV");
  hIVEpEmGam->GetYaxis()->SetTitleSize(0.05);
  hIVEpEmGam->GetYaxis()->SetTitleOffset(1.0);
  
  hIVEpEmGam->SetLineWidth(2);
  
  hIVEpEmGam_cut->SetLineColor(kRed);
  hIVEpEmGam_cut->SetFillColor(kRed);
  hIVEpEmGam_cut->SetFillStyle(3001);
  
  TH1D *hMMPEmX = (TH1D*)in->Get("hMMPEmX");
  TH1D *hMMPEmX_cut = (TH1D*)in->Get("hMMPEmX_cut");
  
  hMMPEmX->SetLineColor(kBlack);
  hMMPEmX->SetLineWidth(2);
  
  hMMPEmX->GetXaxis()->SetTitle("M_{x}(pe^{-}) [GeV]");
  hMMPEmX->SetStats(false);
  hMMPEmX_cut->SetStats(false);
  hMMPEmX->GetXaxis()->SetTitleSize(0.05);
  hMMPEmX->GetXaxis()->SetTitleOffset(0.8);
  hMMPEmX->GetYaxis()->SetTitle("Counts / 35.5 MeV");
  hMMPEmX->GetYaxis()->SetTitleSize(0.05);
  hMMPEmX->GetYaxis()->SetTitleOffset(1.0);
  
  hMMPEmX_cut->SetLineColor(kRed);
  hMMPEmX_cut->SetFillColor(kRed);
  hMMPEmX_cut->SetFillStyle(3001);
  
  
  TH1D *IVrest = (TH1D*)hIVEpEmGam->Clone();
  int critIVbin = IVrest->GetXaxis()->FindFixBin(1.5);
  
  TH1D *MMrest = (TH1D*)hMMPEmX->Clone();
  int critMMbin = MMrest->GetXaxis()->FindFixBin(1.2);
  
  for(int i=1;i<= IVrest->GetNbinsX();i++){
    if(i<critIVbin){
      IVrest->SetBinContent(i,0);
    }
  }
  IVrest->SetFillStyle(3001);
  IVrest->SetFillColor(8);
  IVrest->SetLineColor(8);
  
  for(int i=1;i<= MMrest->GetNbinsX();i++){
    if(i>critMMbin){
      MMrest->SetBinContent(i,0);
    }
  }
  MMrest->SetFillStyle(3001);
  MMrest->SetFillColor(8);
  MMrest->SetLineColor(8);
  
  TLegend *leg1 = new TLegend(0.25,0.7,0.9,0.9);
  leg1->SetTextSize(0.055);

  leg1->SetFillColor(0);
  leg1->AddEntry(IVrest,"Select e^{-}' from e^{-}p#rightarrow e^{-}'pX");
  leg1->AddEntry(hIVEpEmGam_cut,"Select e^{-} from #eta'#rightarrow e^{+}e^{-}#gamma");
  
  TLegend *leg2 = new TLegend(0.25,0.7,0.9,0.9);
  leg2->SetTextSize(0.055);

  leg2->SetFillColor(0);
  leg2->AddEntry(MMrest,"Select e^{-}' from e^{-}p#rightarrow e^{-}'pX");
  leg2->AddEntry(hMMPEmX_cut,"Select e^{-} from #eta'#rightarrow e^{+}e^{-}#gamma");
  
  TCanvas *c = new TCanvas("c","",1200,500);
  c->Divide(2);
  c->cd(1);
  IVrest->GetYaxis()->SetRangeUser(0,570.e03);
  IVrest->Draw();
  
      fitMKVoight(hIVEpEmGam, 0.6, 2.1, 0, 0, 0, 0, 0.957, 0.01, 0.001, 2.5, 1);
  
  
      hIVEpEmGam_cut->Draw("same");
      hIVEpEmGam->Draw("same");
      leg1->Draw("same");
      c->cd(2);
  
      MMrest->Draw();
      fitMKVoight(hMMPEmX, 0.6, 2.1, 0, 0, 0, 0, 0.957, 0.01, 0.001, 2.5, 1);
  
      hMMPEmX_cut->Draw("same");
      hMMPEmX->Draw("same");
      leg2->Draw("same");
      c->cd();
  
      TH1D *hIVEpEm = (TH1D*)in->Get("hIVEpEm");
      TH1D *hIVEpEm_cut = (TH1D*)in->Get("hIVEpEm_cut");
      TH1D *hEpEm_contam = (TH1D*)in->Get("hEpEm_contam");
  
  
      hIVEpEm->GetXaxis()->SetTitle("M(e^{+}e^{-}) [GeV]");
      hIVEpEm->GetXaxis()->SetTitleSize(0.05);
      hIVEpEm->GetXaxis()->SetTitleOffset(0.8);
  
      hIVEpEm->GetYaxis()->SetTitle("Counts / 40 MeV");
      hIVEpEm->GetYaxis()->SetTitleSize(0.05);
      hIVEpEm->GetYaxis()->SetTitleOffset(0.8);
  
      hIVEpEm->SetLineColor(kBlack);
      hIVEpEm->SetLineWidth(2);
  
  
      hIVEpEm_cut->SetLineColor(kRed);
      hIVEpEm_cut->SetFillColor(kRed);
      hIVEpEm_cut->SetFillStyle(3001);
  
  
      hEpEm_contam->SetLineColor(8);
      hEpEm_contam->SetFillColor(8);
      hEpEm_contam->SetFillStyle(3001);
  
      TLegend *leg3 = new TLegend(0.35,0.77,0.9,0.9);
  leg3->SetTextSize(0.055);

      leg3->SetFillColor(0);
      //leg3->AddEntry((TObject*)0,"","");

      leg3->AddEntry(hIVEpEm_cut,"Select e^{-} from #eta'#rightarrowe^{+}e^{-}#gamma");
      leg3->AddEntry(hEpEm_contam,"Select e^{-}' from e^{-}p#rightarrow e^{-}'pX");
  
      TCanvas *c2 = new TCanvas("c2","",1);
      hIVEpEm->Draw();
      hIVEpEm_cut->Rebin(4);
      hIVEpEm_cut->Draw("same");
      hEpEm_contam->Rebin(4);
      hEpEm_contam->Draw("same");
      hIVEpEm->Draw("same");
      leg3->Draw("same");
      c2->SetLogy();
  
  
      TH1D *hEmP = (TH1D*)in->Get("hEmP");
      TH1D *hEmP_cut = (TH1D*)in->Get("hEmP_cut");
  
      hEmP->GetXaxis()->SetTitle("Momentum of e^{-} [GeV/c]");
      hEmP->GetXaxis()->SetTitleSize(0.05);
      hEmP->GetXaxis()->SetTitleOffset(0.8);
  
      hEmP->GetYaxis()->SetTitle("Counts / 110 MeV");
      hEmP->GetYaxis()->SetTitleSize(0.05);
      hEmP->GetYaxis()->SetTitleOffset(0.8);
  
      hEmP->SetLineColor(kBlack);
      hEmP->SetLineWidth(2);
  
      hEmP_cut->SetLineColor(kRed);
      hEmP_cut->SetFillColor(kRed);
      hEmP_cut->SetFillStyle(3001);
  
      TLegend *leg4 = new TLegend(0.6,0.6,0.9,0.9);
      leg4->SetFillStyle(0);
      leg4->AddEntry(hEmP_cut,"Select e^{-} from #eta'#rightarrow e^{+}e^{-}#gamma");
  
     // TCanvas *c3 = new TCanvas("c3","",1);
//      hEmP->Draw();
//      hEmP_cut->Draw("same");
//      leg4->Draw("same");
      //c3->cd();
  
  
  
  
  
  
}
