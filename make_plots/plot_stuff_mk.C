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
Double_t Pole_FFII(Double_t *x, Double_t *par) {
  
  Double_t Lambdasqr= par[1];
  Double_t gammasqr = par[2];
  Double_t qsqr = (x[0]*x[0]);
  Double_t f = 0.;
  Double_t Norm = par[0];
  
  Double_t nom = Lambdasqr*(Lambdasqr + gammasqr);
  Double_t demon = pow((Lambdasqr-qsqr),2) + Lambdasqr*gammasqr;
  f = nom/demon;
  return Norm*f;
  //return Norm*abs(f*f);
  
  
}

Double_t Eval_Kroll_wada(Double_t *x, Double_t *par) {
  Double_t structure_constant = 1./(137.035999074);
  Double_t pi = TMath::Pi();
  Double_t Melectron = 0.000510999; //Electron Mass
  Double_t M = par[0];
  Double_t fpi = 0.0924; //MeV
  Double_t f8pi = 1.3*fpi;
  Double_t f0pi = 1.04*fpi;
  Double_t theta_mix = 10.6*pi/180.;
  
  Double_t Mp_const = structure_constant/(pi*fpi)*(1./sqrt(3.0))*(fpi/f8pi*sin(theta_mix) + 2.0*sqrt(2.0)*fpi/f0pi*cos(theta_mix));
  
  Double_t Norm = 1./(64.*pi*4.)*pow(M,3)*pow(Mp_const,2); //factor 4 unknown
  Double_t qsqr = (x[0]*x[0]);
  Double_t f=0.;
  
  f = (1./Norm)*(4.0*structure_constant/(3.*pi)*(1./sqrt(qsqr))*sqrt(1 - 4.*(pow(Melectron,2)/qsqr)) * (1. + 2.0 *(pow(Melectron,2)/qsqr)) * (pow((1. - qsqr/(pow(M,2))),3)) ) ;
  
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



void plot_stuff_mk(){

  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  gROOT->ForceStyle();
  TFile *in = TFile::Open("IncLusive_Plots_vmdtest.root");
  TFile *inflat = TFile::Open("IncLusive_Plots_flat.root");

  
  //MK stuff

  TCanvas *cmkII =  new TCanvas("cmkII","cmkII",1200,500);
  //cmkII->Divide(1,2);
  cmkII->cd();
  cmkII->SetLogy();
  TH1D *hIVEpEm_cut_clone = (TH1D*)in->Get("hIVEpEm_cut_clone");
  hIVEpEm_cut_clone->SetLineColor(kRed);
  hIVEpEm_cut_clone->SetYTitle("Expected Counts / 10 MeV");
  hIVEpEm_cut_clone->GetYaxis()->SetTitleSize(0.05);
  hIVEpEm_cut_clone->GetYaxis()->SetTitleOffset(0.8);
  hIVEpEm_cut_clone->GetXaxis()->SetTitle("M(e^{+}e^{-}) [GeV]");
  hIVEpEm_cut_clone->GetXaxis()->SetTitleSize(0.05);
  hIVEpEm_cut_clone->GetXaxis()->SetTitleOffset(0.8);
  TH1D *hEpEm_corrected = (TH1D*)in->Get("hEpEm_corrected");
  hEpEm_corrected->SetYTitle("Counts / 10 MeV");
  hEpEm_corrected->GetYaxis()->SetTitleSize(0.05);
  hEpEm_corrected->GetYaxis()->SetTitleOffset(0.8);
  hEpEm_corrected->GetXaxis()->SetTitle("M(e^{+}e^{-}) [GeV]");
  hEpEm_corrected->GetXaxis()->SetTitleSize(0.05);
  hEpEm_corrected->Draw("EP");
  hIVEpEm_cut_clone->Draw("EP same");

  TLegend *legmkII = new TLegend(0.45,0.7,0.9,0.9);
  legmkII->SetTextSize(0.05);
  legmkII->SetFillColor(0);
  legmkII->AddEntry(hIVEpEm_cut_clone,"Expected counts in 80 days (N_{tot}=52785)","l");

  legmkII->AddEntry(hEpEm_corrected,"Acceptance corrected counts in 80 days","l");
  legmkII->Draw("same");
  
  
  TCanvas *cmkIII =  new TCanvas("cmkIII","cmkIII",1200,500);
  cmkIII->cd();
  TH1D *hEpEm_acceptance = (TH1D*)in->Get("hEpEm_acceptance");
  hEpEm_acceptance->SetLineColor(kRed);
  hEpEm_acceptance->SetYTitle("Acceptance / 10 MeV");
  hEpEm_acceptance->GetYaxis()->SetTitleSize(0.05);
  hEpEm_acceptance->GetYaxis()->SetTitleOffset(0.8);
  hEpEm_acceptance->GetXaxis()->SetTitle("M(e^{+}e^{-}) [GeV]");
  hEpEm_acceptance->GetXaxis()->SetTitleSize(0.05);
  hEpEm_acceptance->GetXaxis()->SetTitleOffset(0.8);
  TH1D *hEpEm_acceptance_flat = (TH1D*)inflat->Get("hEpEm_acceptance");
  hEpEm_acceptance_flat->SetYTitle("Acceptance / 10 MeV");
  hEpEm_acceptance_flat->GetYaxis()->SetTitleSize(0.05);
  hEpEm_acceptance_flat->GetYaxis()->SetTitleOffset(0.8);
  hEpEm_acceptance_flat->GetXaxis()->SetTitle("M(e^{+}e^{-}) [GeV]");
  hEpEm_acceptance_flat->GetXaxis()->SetTitleSize(0.05);
  hEpEm_acceptance_flat->GetXaxis()->SetTitleOffset(0.8);
  hEpEm_acceptance->GetYaxis()->SetRangeUser(0,0.16);
  hEpEm_acceptance->Draw("EP");
  hEpEm_acceptance_flat->Draw("EP same");
  
  TLegend *legmkIII = new TLegend(0.11,0.7,0.45,0.9);
  legmkIII->SetTextSize(0.05);
  legmkIII->SetFillColor(0);
  legmkIII->AddEntry(hEpEm_acceptance,"QED+VMD M(e^{+}e^{-}) Acceptance","l");
  
  legmkIII->AddEntry(hEpEm_acceptance_flat,"Flat M(e^{+}e^{-}) Acceptance","l");
  legmkIII->Draw("same");
  
  
  //
  double QED_par[8] = {0.957}; //{Mass}
  TF1 *QED_norm = new TF1("QED_norm",Eval_Kroll_wada,0.,1.,1);
  QED_norm->SetParameters(&QED_par[0]);
  
  TH1D *hEpEm_corrected_flat = new TH1D("hEpEm_corrected_flat","hEpEm_corrected_flat",100,0.0,1);
  TH1D *hEpEm_QEDnorm_flat = new TH1D("hEpEm_QEDnorm_flat","hEpEm_QEDnorm_flat",100,0.0,1);
  TH1D *hEpEm_corrected_sys = new TH1D("hEpEm_corrected_sys","hEpEm_corrected_sys",100,0.0,1);

  TF1 *fsys = new TF1("fsys","0.05*x",0,1);

  for (int i = 1; i<hIVEpEm_cut_clone->GetNbinsX(); i++) {
    
    Double_t bin_factor = 1.65; //this needs to be solved at somepoint
    Double_t intotal_acceptance = hEpEm_acceptance_flat->GetBinContent(i)*bin_factor;
    Double_t intotal_events_upper;
    
    if (intotal_acceptance == 0) {
      intotal_events_upper = 0;
    }else{
      intotal_events_upper = hIVEpEm_cut_clone->GetBinContent(i)/intotal_acceptance;
    }
    
    hEpEm_corrected_flat->SetBinContent(i,intotal_events_upper);
    hEpEm_corrected_flat->SetBinError(i,sqrt(intotal_events_upper));
    
    Double_t QED_factor = QED_norm->Eval(hIVEpEm_cut_clone->GetBinCenter(i));//2.0e-06;
    hEpEm_QEDnorm_flat->SetBinContent(i,intotal_events_upper/QED_factor);
    //hEpEm_QEDnorm->SetBinError(i,sqrt(intotal_events_upper/QED_factor));
    

    //For systematic
    Double_t num_sys = fsys->Eval(hIVEpEm_cut_clone->GetBinCenter(i));
    Double_t intotal_acceptance_sys = hEpEm_acceptance->GetBinContent(i)*bin_factor*(1. + num_sys);
    Double_t intotal_events_upper_sys;
    
    if (intotal_acceptance_sys == 0) {
      intotal_events_upper_sys = 0;
    }else{
      intotal_events_upper_sys = hIVEpEm_cut_clone->GetBinContent(i)/intotal_acceptance_sys;
    }
    hEpEm_corrected_sys->SetBinContent(i,intotal_events_upper_sys/QED_factor);
    //hEpEm_QEDnorm->SetBinError(i,sqrt(intotal_events_upper/QED_factor));
    //cout<<bin_factor*num_sys*hEpEm_acceptance->GetBinContent(i)<<"  "<<intotal_events_upper_sys<<"  "<<i<<endl;

  }
  cout<<"Total & "<<TMath::Nint(hIVEpEm_cut_clone->Integral())<<" & "<<TMath::Nint(sqrt(hIVEpEm_cut_clone->Integral("")))<<"\\\\"<<endl;

  
  
  //
  
  double pole_par[3] = {10, 0.59,0.0144}; //{A,Lambda,Gamma}
  
  TF1 *FF_fittepole = new TF1("FF_fittepole", Pole_FFII,0.02,0.92,3);
  FF_fittepole->SetParameters(&pole_par[0]);
  FF_fittepole->SetParLimits(0,0.,100.);
  FF_fittepole->SetParLimits(1,pole_par[1] - pole_par[1]*0.1,pole_par[1] + pole_par[1]*0.1);
  FF_fittepole->SetParLimits(2,pole_par[2] - pole_par[2]*0.5,pole_par[2] + pole_par[2]*0.5);
  
  //for flat
  double pole_parII[3] = {10, 0.59,0.0144}; //{A,Lambda,Gamma}
  
  TF1 *FF_fittepoleII = new TF1("FF_fittepoleII", Pole_FFII,0.02,0.85,3);
  FF_fittepoleII->SetParameters(&pole_parII[0]);
  FF_fittepoleII->SetParLimits(0,0.,100.);
  FF_fittepoleII->SetParLimits(1,pole_parII[1] - pole_parII[1]*0.1,pole_parII[1] + pole_parII[1]*0.1);
  FF_fittepoleII->SetParLimits(2,pole_parII[2] - pole_parII[2]*0.5,pole_parII[2] + pole_parII[2]*0.5);
  FF_fittepoleII->SetLineColor(kBlue);
  //TF1 *FF_fitterGaus = new TF1("FF_fitterGaus","gaus",0.55,0.92);
  
  TCanvas *cmkI =  new TCanvas("cmkI","cmkI",1200,500);
  
  TH1D *hEpEm_QEDnorm = (TH1D*)in->Get("hEpEm_QEDnorm");
  cmkI->cd();
  hEpEm_QEDnorm->SetTitle("Expected distribution of |F(q^{2})|^{2}");
  hEpEm_QEDnorm->SetLineColor(kBlack);
  hEpEm_QEDnorm->SetYTitle("|F(q^{2})|^{2}");
  hEpEm_QEDnorm->GetYaxis()->SetTitleSize(0.05);
  hEpEm_QEDnorm->GetYaxis()->SetTitleOffset(0.8);
  hEpEm_QEDnorm->GetXaxis()->SetTitle("M(e^{+}e^{-}) [GeV]");
  hEpEm_QEDnorm->GetXaxis()->SetTitleSize(0.05);
  hEpEm_QEDnorm->GetXaxis()->SetTitleOffset(0.8);
  hEpEm_QEDnorm->Draw("EP");
  hEpEm_QEDnorm_flat->SetLineColor(8);
  hEpEm_QEDnorm_flat->Draw("EP same");
  hEpEm_QEDnorm_flat->Fit("FF_fittepoleII","REM");
  hEpEm_QEDnorm->Fit("FF_fittepole","REM+");

  Double_t Lambda = FF_fittepole->GetParameter(1);
  Double_t Lambdaerr = FF_fittepole->GetParError(1);
  Double_t bn = 1./Lambda;
  Double_t bn_err = Lambdaerr/(bn*bn);
  
  
  Double_t Lambda_flat = FF_fittepoleII->GetParameter(1);
  Double_t Lambdaerr_flat = FF_fittepoleII->GetParError(1);
  Double_t bn_flat = 1./Lambda_flat;
  Double_t bn_err_flat = Lambdaerr_flat/(bn_flat*bn_flat);
  
  
  TString sLambda = Form("#Lambda^{2}_{fit}  = %2.4f #pm %2.4f  ", Lambda,Lambdaerr);
  TString sbn = Form("b_{n} = %2.4f #pm %2.4f  ", bn,bn_err);
  
  Double_t Lambda_gen = 0.5776;
  Double_t bn_gen = 1./0.5776;
  TString sLambda_gen = Form("#Lambda^{2}_{gen}  = %2.4f  ", 0.5776);
  TString sbn_gen = Form("b_{n gen} = %2.4f  ", 1./0.5776);
  
  TString sLambda_flat = Form("#Lambda^{2}_{fit}  = %2.4f #pm %2.4f  ", Lambda_flat,Lambdaerr_flat);
  TString sbn_flat = Form("b_{n} = %2.4f #pm %2.4f  ", bn_flat,bn_err_flat);


  
  TLegend *legmkI = new TLegend(0.12,0.55,0.37,0.9);
  legmkI->SetTextSize(0.04);

  legmkI->SetHeader("\t QED+VMD M(e^{+}e^{-}) Acceptance");
  legmkI->SetFillColor(0);
  legmkI->AddEntry(FF_fittepole,sLambda,"l");
  legmkI->AddEntry((TObject*)0,sbn,"");
  legmkI->AddEntry((TObject*)0,sLambda_gen,"");
  legmkI->AddEntry((TObject*)0,sbn_gen,"");
  legmkI->Draw("same");
  
  TLegend *legmkI_I = new TLegend(0.37,0.55,0.57,0.9);
  legmkI_I->SetTextSize(0.04);

  legmkI_I->SetHeader("\t Flat M(e^{+}e^{-}) Acceptance");
  legmkI_I->SetFillColor(0);
  legmkI_I->AddEntry(FF_fittepoleII,sLambda_flat,"l");
  legmkI_I->AddEntry((TObject*)0,sbn_flat,"");
  legmkI_I->AddEntry((TObject*)0,"","");
  legmkI_I->AddEntry((TObject*)0,"","");

  legmkI_I->Draw("same");

  
  TCanvas *cmksys =  new TCanvas("cmksys","cmksys",1200,500);
  cmksys->cd();
  double pole_parsys[3] = {10, 0.59,0.0144}; //{A,Lambda,Gamma}
  
  TF1 *FF_fittepolesys = new TF1("FF_fittepolesys", Pole_FFII,0.02,0.92,3);
  FF_fittepolesys->SetParameters(&pole_parsys[0]);
  FF_fittepolesys->SetParLimits(0,0.,100.);
  FF_fittepolesys->SetParLimits(1,pole_parsys[1] - pole_parsys[1]*0.1,pole_parsys[1] + pole_parsys[1]*0.1);
  FF_fittepolesys->SetParLimits(2,pole_parsys[2] - pole_parsys[2]*0.5,pole_parsys[2] + pole_parsys[2]*0.5);
  FF_fittepolesys->SetLineColor(kBlue);
  //hEpEm_QEDnorm->Draw("EP");
  hEpEm_corrected_sys->SetTitle("Expected distribution of |F(q^{2})|^{2} with 5% systematic on acceptance");
  hEpEm_corrected_sys->SetLineColor(kBlack);
  hEpEm_corrected_sys->SetYTitle("|F(q^{2})|^{2}");
  hEpEm_corrected_sys->GetYaxis()->SetTitleSize(0.05);
  hEpEm_corrected_sys->GetYaxis()->SetTitleOffset(0.8);
  hEpEm_corrected_sys->GetXaxis()->SetTitle("M(e^{+}e^{-}) [GeV]");
  hEpEm_corrected_sys->GetXaxis()->SetTitleSize(0.05);
  hEpEm_corrected_sys->GetXaxis()->SetTitleOffset(0.8);

  hEpEm_corrected_sys->Draw("EP");
  hEpEm_corrected_sys->Fit("FF_fittepolesys","REM");

  Double_t Lambda_sys = FF_fittepolesys->GetParameter(1);
  Double_t Lambdaerr_sys = FF_fittepolesys->GetParError(1);
  Double_t bn_sys = 1./Lambda_sys;
  Double_t bn_err_sys = Lambdaerr_sys/(bn_sys*bn_sys);
 
  TString sLambda_sys = Form("#Lambda^{2}_{fit}  = %2.4f #pm %2.4f  ", Lambda_sys,Lambdaerr_sys);
  TString sbn_sys = Form("b_{n} = %2.4f #pm %2.4f  ", bn_sys,bn_err_sys);

  legmkI->Draw("same");
  
  TLegend *legmksys = new TLegend(0.37,0.55,0.63,0.9);
  legmksys->SetTextSize(0.04);
  
  legmksys->SetHeader("\t QED+VMD M(e^{+}e^{-}) Acceptance ");
  legmksys->SetFillColor(0);
  legmksys->AddEntry((TObject*)0,"with 5% systematic","");

  legmksys->AddEntry(FF_fittepolesys,sLambda_sys,"l");
  legmksys->AddEntry((TObject*)0,sbn_sys,"");
  legmksys->AddEntry((TObject*)0,sLambda_gen,"");
  legmksys->AddEntry((TObject*)0,sbn_gen,"");
  
  legmksys->Draw("same");

  //lets calculate percent error and percent difference
  Double_t Lambda_percent_diff = abs(Lambda_sys - Lambda)/(0.5*(Lambda_sys + Lambda))*100.;
  Double_t bn_percent_diff = abs(bn_sys - bn)/(0.5*(bn_sys + bn))*100.;
  
  Double_t Lambda_percent_err = abs((Lambda_gen - Lambda)/Lambda_gen)*100.;
  Double_t bn_percent_err = abs((bn_gen - bn)/bn_gen)*100.;

  Double_t Lambda_percent_errsys = abs((Lambda_gen - Lambda_sys)/Lambda_gen)*100.;
  Double_t bn_percent_errsys = abs((bn_gen - bn_sys)/bn_gen)*100.;
  
  TString sLambda_perdif = Form("#Lambda^{2}_{sys}  = %2.4f #pm %2.4f  ", Lambda_percent_diff,bn_percent_diff);

  
  cout<<Lambda<<"  "<<bn<<endl;
  cout<<Lambda_sys<<"  "<<bn_sys<<endl;
  cout<<Lambda_gen<<"  "<<bn_gen<<endl;

  
  cout<<"Percent Diff "<<Lambda_percent_diff<<"  "<<bn_percent_diff<<endl;
  cout<<"Percent Error "<<Lambda_percent_err<<"  "<<bn_percent_err<<endl;
  cout<<"Percent Error sys "<<Lambda_percent_errsys<<"  "<<bn_percent_errsys<<endl;

  //cmkI->Print("../figures/counts/result.pdf");
  //cmkIII->Print("../figures/counts/acceptance.pdf");
  //cmkII->Print("../figures/counts/counts.pdf");
  //cmksys->Print("../figures/counts/sys.pdf");
  
  cmkI->Print("../figures/counts/result.C");
  cmkIII->Print("../figures/counts/acceptance.C");
  cmkII->Print("../figures/counts/counts.C");
  cmksys->Print("../figures/counts/sys.C");

}





