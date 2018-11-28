#include "TROOT.h"
#include "TRint.h"
#include "TH1.h"
#include "TF1.h"
#include "TLatex.h"

#include "TH1F.h"
#include "TH2.h"

#include "TCanvas.h"
#include "TChain.h"
#include "TVector3.h"
#include "TMath.h"
#include "TFile.h"
#include "TSystem.h"
#include "TGStatusBar.h"
#include "TSystem.h"
#include "TXMLEngine.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TNtuple.h"
#include "TString.h"
#include "TClonesArray.h"

#include "Math/GSLIntegrator.h"
#include "Math/WrappedTF1.h"

#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <unistd.h>

void fitMKGaus(TH1 * , Double_t , Double_t , Double_t , Double_t , Double_t , Double_t , Double_t , Double_t , Double_t , Int_t);
void fitMKVoight(TH1 *, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Int_t);

Double_t Eval_Kroll_wada(Double_t *, Double_t *);
Double_t VMD_FF(Double_t *, Double_t *);
Double_t Pole_FF(Double_t *, Double_t *);
Double_t Pole_FFII(Double_t *, Double_t *);

Double_t myVoight(Double_t *, Double_t *);

TString XSection_Hist_num(Double_t);

Double_t OneBoostComparison(TLorentzVector daughter1, TLorentzVector parent){
  
  daughter1.Boost(-parent.BoostVector());
  
  Double_t CosTheta = daughter1.CosTheta();
  
  if(TMath::IsNaN(CosTheta)){return -1000.;}
  else{return CosTheta;}
  
}


void Inclusive_Proposal_Plots(){
  
  //gStyle->SetPaintTextFormat("4.2f");
  TFile *fhist = new TFile("/Users/michaelkunkel/WORK/CLAS/CLAS12/CODES/SIMUALTION/ETAPRIME/ElectroProduction/ROOT_FILES/Hists_etaP.root","READ");
  
  TChain *Reaction = new TChain("etaP_MC");
  Reaction->Add("VMD/etaP_fastMC_vmd_EpEm.root");
  //Reaction->Add("FLAT/etaP_fastMC_flat_EpEm.root");
  
  //  Reaction->SetBranchStatus("*",0);
  //  Reaction->SetBranchStatus("EmPEmEpGam_acc",1);
  //  Reaction->SetBranchStatus("PEmEpGam_acc",1);
  //  Reaction->SetBranchStatus("EmEmEpGam_acc",1);
  //  Reaction->SetBranchStatus("EmPEmEpGam_acc",1);
  //  Reaction->SetBranchStatus("EmPEmEpGam_acc",1);

  TLorentzVector vEm_gen, vEp_gen, vP_gen, vGamma_gen, vScatter_gen;
  TLorentzVector vEm, vEp, vP, vGamma, vScatter;
  Double_t M_el = 0.000511;
  Double_t M_P = 0.938272;
  TLorentzVector vT(0.0, 0.0, 0.0, M_P);
  
  Float_t EmPEmEpGam_acc; // All detected
  Float_t PEmEpGam_acc;   // missing scattered
  Float_t EmEmEpGam_acc;  // missing proton
  Float_t EmPEpGam_acc;  // missing dalitz electron
  Float_t AnyEmPEpGam_acc;       // any electron P Ep Gam
  
  Reaction->SetBranchAddress("AnyEmPEpGam_acc",&AnyEmPEpGam_acc);
  Reaction->SetBranchAddress("PEmEpGam_acc",&PEmEpGam_acc);
  Reaction->SetBranchAddress("EmPEpGam_acc",&EmPEpGam_acc);
  
  //Generated Quantities
  Double_t P_P_gen, P_Theta_gen, P_Phi_gen;
  Double_t Ep_P_gen, Ep_Theta_gen, Ep_Phi_gen, Ep_E_gen;
  Double_t Em_P_gen, Em_Theta_gen, Em_Phi_gen, Em_E_gen;
  Double_t Gam_P_gen, Gam_Theta_gen, Gam_Phi_gen;
  Double_t Em_s_P_gen, Em_s_Theta_gen, Em_s_Phi_gen, Em_s_E_gen;
  
  Reaction->SetBranchAddress("P_P_gen",&P_P_gen);
  Reaction->SetBranchAddress("P_Theta_gen",&P_Theta_gen);
  Reaction->SetBranchAddress("P_Phi_gen",&P_Phi_gen);
  Reaction->SetBranchAddress("Ep_P_gen",&Ep_P_gen);
  Reaction->SetBranchAddress("Ep_Theta_gen",&Ep_Theta_gen);
  Reaction->SetBranchAddress("Ep_Phi_gen",&Ep_Phi_gen);
  Reaction->SetBranchAddress("Ep_E_gen",&Ep_E_gen);
  Reaction->SetBranchAddress("Em_P_gen",&Em_P_gen);
  Reaction->SetBranchAddress("Em_Theta_gen",&Em_Theta_gen);
  Reaction->SetBranchAddress("Em_Phi_gen",&Em_Phi_gen);
  Reaction->SetBranchAddress("Em_E_gen",&Em_E_gen);
  Reaction->SetBranchAddress("Gam_P_gen",&Gam_P_gen);
  Reaction->SetBranchAddress("Gam_Theta_gen",&Gam_Theta_gen);
  Reaction->SetBranchAddress("Gam_Phi_gen",&Gam_Phi_gen);
  Reaction->SetBranchAddress("Em_s_P_gen",&Em_s_P_gen);
  Reaction->SetBranchAddress("Em_s_Theta_gen",&Em_s_Theta_gen);
  Reaction->SetBranchAddress("Em_s_Phi_gen",&Em_s_Phi_gen);
  Reaction->SetBranchAddress("Em_s_E_gen",&Em_s_E_gen);
  
  
  //Reconstructed quantites
  Double_t P_P, P_Theta, P_Phi;
  Double_t Ep_P, Ep_Theta, Ep_Phi, Ep_E;
  Double_t Em_P, Em_Theta, Em_Phi, Em_E;
  Double_t Gam_P, Gam_Theta, Gam_Phi;
  Double_t Em_s_P, Em_s_Theta, Em_s_Phi, Em_s_E;
  Reaction->SetBranchAddress("P_P",&P_P);
  Reaction->SetBranchAddress("P_Theta",&P_Theta);
  Reaction->SetBranchAddress("P_Phi",&P_Phi);
  Reaction->SetBranchAddress("Ep_P",&Ep_P);
  Reaction->SetBranchAddress("Ep_Theta",&Ep_Theta);
  Reaction->SetBranchAddress("Ep_Phi",&Ep_Phi);
  Reaction->SetBranchAddress("Ep_E",&Ep_E);
  Reaction->SetBranchAddress("Em_P",&Em_P);
  Reaction->SetBranchAddress("Em_Theta",&Em_Theta);
  Reaction->SetBranchAddress("Em_Phi",&Em_Phi);
  Reaction->SetBranchAddress("Em_E",&Em_E);
  Reaction->SetBranchAddress("Gam_P",&Gam_P);
  Reaction->SetBranchAddress("Gam_Theta",&Gam_Theta);
  Reaction->SetBranchAddress("Gam_Phi",&Gam_Phi);
  Reaction->SetBranchAddress("Em_s_P",&Em_s_P);
  Reaction->SetBranchAddress("Em_s_Theta",&Em_s_Theta);
  Reaction->SetBranchAddress("Em_s_Phi",&Em_s_Phi);
  Reaction->SetBranchAddress("Em_s_E",&Em_s_E);
  
  Double_t IV_EpEm, IV_EpEm_gen, Qsq_gen, Qsq;
  Reaction->SetBranchAddress("IV_EpEm",&IV_EpEm);
  Reaction->SetBranchAddress("IV_EpEm_gen",&IV_EpEm_gen);
  Reaction->SetBranchAddress("Qsq_gen",&Qsq_gen);
  Reaction->SetBranchAddress("Qsq",&Qsq);
  
  //Lets get total number of entries
  Int_t nevents = Reaction->GetEntries();
  cout<<nevents<<endl;
  TH1D *hIVEpEmGam = new TH1D("hIVEpEmGam","hIVEpEmGam",100,0.45,4.);
  TH1D *hIVEpEm = new TH1D("hIVEpEm","hIVEpEm",100,0.0,4.);
  
  TH1D *hEmP = new TH1D("hEmP","hEmP",100,0.0,11.);
  TH1D *hEm1P = new TH1D("hEm1P","hEm1P",100,0.0,11.);
  TH1D *hEm2P = new TH1D("hEm2P","hEm2P",100,0.0,11.);
  
  TH1D *hIVEpEm1Gam = new TH1D("hIVEpEm1Gam","hIVEpEm1Gam",100,0.45,4.);
  TH1D *hIVEpEm2Gam = new TH1D("hIVEpEm2Gam","hIVEpEm2Gam",100,0.45,4.);
  
  TH1D *hMMPEmX= new TH1D("hMMPEmX","hMMPEmX",100,0.45,4.);
  TH1D *hMMPEm1X= new TH1D("hMMPEm1X","hMMPEm1X",100,0.45,4.);
  TH1D *hMMPEm2X= new TH1D("hMMPEm2X","hMMPEm2X",100,0.45,4.);
  
  TH1D *hMM2PEmEpGam = new TH1D("hMM2PEmEpGam","hMM2PEmEpGam",100,-1.4,1.4);
  TH1D *hQsqr_inclusive = new TH1D("hQsqr_inclusive","hQsqr_inclusive",100,0,22.4);

  TH1D *hQsqr_recon = new TH1D("hQsqr_recon","hQsqr_recon",100,0,22.4);

  //cut
  TH1D *hIVEpEmGam_cut = new TH1D("hIVEpEmGam_cut","hIVEpEmGam_cut",100,0.45,4.);
  TH1D *hIVEpEm_cut = new TH1D("hIVEpEm_cut","hIVEpEm_cut",100,0.0,1);
  
  TH1D *hEmP_cut = new TH1D("hEmP_cut","hEmP_cut",100,0.0,11.);
  TH1D *hEm1P_cut = new TH1D("hEm1P_cut","hEm1P_cut",100,0.0,11.);
  TH1D *hEm2P_cut = new TH1D("hEm2P_cut","hEm2P_cut",100,0.0,11.);
  
  TH1D *hIVEpEm1Gam_cut = new TH1D("hIVEpEm1Gam_cut","hIVEpEm1Gam_cut",100,0.45,4.);
  TH1D *hIVEpEm2Gam_cut = new TH1D("hIVEpEm2Gam_cut","hIVEpEm2Gam_cut",100,0.45,4.);
  
  TH1D *hMMPEmX_cut = new TH1D("hMMPEmX_cut","hMMPEmX_cut",100,0.45,4.);
  TH1D *hMMPEm1X_cut = new TH1D("hMMPEm1X_cut","hMMPEm1X_cut",100,0.45,4.);
  TH1D *hMMPEm2X_cut = new TH1D("hMMPEm2X_cut","hMMPEm2X_cut",100,0.45,4.);
  
  TH1D *hMM2PEmEpGam_cut = new TH1D("hMM2PEmEpGam_cut","hMM2PEmEpGam_cut",100,-1.4,1.4);
  TH1D *hQsqr_inclusive_cut = new TH1D("hQsqr_inclusive_cut","hQsqr_inclusive_cut",100,0,22.4);

  TH1D *hQsqr_recon_cut = new TH1D("hQsqr_recon_cut","hQsqr_recon_cut",100,0,22.4);

  //
  TH1D *hIVEpEm_excl = new TH1D("hIVEpEm_excl","hIVEpEm_excl",100,0.0,1);
  TH1D *hIVEpEm_gen = new TH1D("hIVEpEm_gen","hIVEpEm_gen",100,0.0,1);
  
  TH2D *hIVEpEm_Qsq_gen = new TH2D("hIVEpEm_Qsq_gen","hIVEpEm_Qsq_gen",100,0.0,1,90,0,18);
  TH2D *hIVEpEm_Qsq = new TH2D("hIVEpEm_Qsq","hIVEpEm_Qsq",100,0.0,1,90,0,18);
  TH2D *hEm_sTheta_Qsq = new TH2D("hEm_sTheta_Qsq","hEm_sTheta_Qsq",100,0.0,50.,100,0,18);
  TH2D *hEm_sE_Qsq = new TH2D("hEm_sE_Qsq","hEm_sE_Qsq",100,0.0,11.,100,0,18);
  TH2D *hEm_sE_Theta = new TH2D("hEm_sE_Theta","hEm_sE_Theta",100,0.0,11.,100,0,50);

  TH1D *hQsqr_gen = new TH1D("hQsqr_gen","hQsqr_gen",100,0,22.4);
  TH1D *hQsqr_inclusive_gen = new TH1D("hQsqr_inclusive_gen","hQsqr_inclusive_gen",100,0,22.4);
  TH1D *hQsqr_inclusive_gen_acc = new TH1D("hQsqr_inclusive_gen_acc","hQsqr_inclusive_gen_acc",100,0,22.4);
  TH1D *hQsqr_inclusive_gen_acccut = new TH1D("hQsqr_inclusive_gen_acccut","hQsqr_inclusive_gen_acccut",100,0,22.4);
  TH1D *hMMPEmEpGam_gen = new TH1D("hMMPEmEpGam_gen","hMMPEmEpGam_gen",100,-0.004,0.004);

  
  
  TH1D *hIVEpEm_AvXsection = new TH1D("hIVEpEm_AvXsection","hIVEpEm_AvXsection",100,0.0,1);
  hIVEpEm_AvXsection->Sumw2();
  
  TH2D *W_v_CosTh = new TH2D("W_v_CosTh","W_v_CosTh",100,-1,1,100,1.4,5.4);
  //Test histo
  TH2D *hscat_reconE_theta = new TH2D("hscat_reconE_theta","hscat_reconE_theta",100,0,50,100,0,50);
  TH2D *hscat_reconE_Qsqr = new TH2D("hscat_reconE_Qsqr","hscat_reconE_Qsqr",100,0,18,100,0,18);

  for (int i=0; i<nevents/100; i++)//nevents
  {
    //f->cd();
    if(!(i%500000)) std::cout << "\r done " << i << " out of " << nevents << " ==> " << double(i)*100.0/double(nevents) << "%" << flush;
    //if (processed_event>=10000) {break;}
    Reaction->GetEntry(i);
    //gen
    Double_t P_px_gen = (P_P_gen*sin(P_Theta_gen)*cos(P_Phi_gen));
    Double_t P_py_gen = (P_P_gen*sin(P_Theta_gen)*sin(P_Phi_gen));
    Double_t P_pz_gen = (P_P_gen*cos(P_Theta_gen));
    
    Double_t Em_px_gen = Em_P_gen*sin(Em_Theta_gen)*cos(Em_Phi_gen);
    Double_t Em_py_gen = Em_P_gen*sin(Em_Theta_gen)*sin(Em_Phi_gen);
    Double_t Em_pz_gen = Em_P_gen*cos(Em_Theta_gen);
    
    Double_t Ep_px_gen = Ep_P_gen*sin(Ep_Theta_gen)*cos(Ep_Phi_gen);
    Double_t Ep_py_gen = Ep_P_gen*sin(Ep_Theta_gen)*sin(Ep_Phi_gen);
    Double_t Ep_pz_gen = Ep_P_gen*cos(Ep_Theta_gen);
    
    Double_t Em_s_px_gen = Em_s_P_gen*sin(Em_s_Theta_gen)*cos(Em_s_Phi_gen);
    Double_t Em_s_py_gen = Em_s_P_gen*sin(Em_s_Theta_gen)*sin(Em_s_Phi_gen);
    Double_t Em_s_pz_gen = Em_s_P_gen*cos(Em_s_Theta_gen);
    
    Double_t Gam_px_gen = Gam_P_gen*sin(Gam_Theta_gen)*cos(Gam_Phi_gen);
    Double_t Gam_py_gen = Gam_P_gen*sin(Gam_Theta_gen)*sin(Gam_Phi_gen);
    Double_t Gam_pz_gen = Gam_P_gen*cos(Gam_Theta_gen);
    

    
    vP_gen.SetPxPyPzE(P_px_gen, P_py_gen, P_pz_gen, sqrt(P_P_gen*P_P_gen + M_P*M_P));
    vEp_gen.SetPxPyPzE(Ep_px_gen, Ep_py_gen, Ep_pz_gen, sqrt(Ep_P_gen*Ep_P_gen + M_el*M_el));
    vGamma_gen.SetPxPyPzE(Gam_px_gen, Gam_py_gen, Gam_pz_gen, Gam_P_gen);
    vEm_gen.SetPxPyPzE(Em_px_gen, Em_py_gen, Em_pz_gen, sqrt(Em_P_gen*Em_P_gen + M_el*M_el));
    vScatter_gen.SetPxPyPzE(Em_s_px_gen, Em_s_py_gen, Em_s_pz_gen, sqrt(Em_s_P_gen*Em_s_P_gen + M_el*M_el));
    
    Double_t beam_mass_pseudo = ((vScatter_gen + vP_gen + vEm_gen + vEp_gen +vGamma_gen) - vT).M();
    TLorentzVector vEg(0.0, 0.0, 11.0, sqrt(pow(11.0,2) + pow(beam_mass_pseudo,2))); //why beam_mass_pseudo??? Fuck if I know!!

    
    TLorentzVector vMMPEmEpGam_gen, vQsqr_inc_gen;
    vMMPEmEpGam_gen = vEg + vT - (vP_gen + vEm_gen + vEp_gen +vGamma_gen);
    vQsqr_inc_gen = vEg - vMMPEmEpGam_gen;
    
    hIVEpEm_gen->Fill(IV_EpEm_gen);
    hIVEpEm_Qsq_gen->Fill(IV_EpEm_gen,-Qsq_gen);
    hQsqr_gen->Fill(-Qsq_gen);
    hQsqr_inclusive_gen->Fill(-vQsqr_inc_gen.M2());
    hMMPEmEpGam_gen->Fill(vMMPEmEpGam_gen.M());
    //First lets do an analysis with either lepton detected
    if (AnyEmPEpGam_acc ==1.) {//AnyEmPEpGam_acc
      Double_t P_px = (P_P*sin(P_Theta)*cos(P_Phi));
      Double_t P_py = (P_P*sin(P_Theta)*sin(P_Phi));
      Double_t P_pz = (P_P*cos(P_Theta));
      
      Double_t Ep_px = Ep_P*sin(Ep_Theta)*cos(Ep_Phi);
      Double_t Ep_py = Ep_P*sin(Ep_Theta)*sin(Ep_Phi);
      Double_t Ep_pz = Ep_P*cos(Ep_Theta);
      
      Double_t Gam_px = Gam_P*sin(Gam_Theta)*cos(Gam_Phi);
      Double_t Gam_py = Gam_P*sin(Gam_Theta)*sin(Gam_Phi);
      Double_t Gam_pz = Gam_P*cos(Gam_Theta);
      //Lets set the lorentz vectors
      vP.SetPxPyPzE(P_px, P_py, P_pz, sqrt(P_P*P_P + M_P*M_P));
      vEp.SetPxPyPzE(Ep_px, Ep_py, Ep_pz, sqrt(Ep_P*Ep_P + M_el*M_el));
      vGamma.SetPxPyPzE(Gam_px, Gam_py, Gam_pz, Gam_P);
      
      //check which electron
      if (PEmEpGam_acc ==1. && EmPEpGam_acc ==0.) { //Dalitz electron accepted only
        //cout<<PEmEpGam_acc<<"  "<<EmPEpGam_acc<<endl;

        Double_t Em_px = Em_P*sin(Em_Theta)*cos(Em_Phi);
        Double_t Em_py = Em_P*sin(Em_Theta)*sin(Em_Phi);
        Double_t Em_pz = Em_P*cos(Em_Theta);
        
        TLorentzVector vEmX, vIV_EpEmGam, vIV_EpEm; //When one electron is detects, both
        TLorentzVector vMM_PEmX;
        TLorentzVector vMMPEmEpGam, vQsqr_inc;

        vEmX.SetPxPyPzE(Em_px, Em_py, Em_pz, sqrt(Em_P*Em_P + M_el*M_el));
        
        vIV_EpEmGam = (vEp + vEmX + vGamma);
        vIV_EpEm = (vEp + vEmX);
        vMM_PEmX = vEg + vT - (vP + vEmX);
        vMMPEmEpGam = vEg + vT - (vP + vEmX + vEp +vGamma);
        vQsqr_inc = vEg - vMMPEmEpGam;
        
        hIVEpEmGam->Fill(vIV_EpEmGam.M());
        hMMPEmX->Fill(vMM_PEmX.M());
        hIVEpEm->Fill(vIV_EpEm.M());
        
        hMMPEm1X->Fill(vMM_PEmX.M());
        
        hEmP->Fill(vEmX.Vect().Mag());
        hEm1P->Fill(vEmX.Vect().Mag());
        
        hMM2PEmEpGam->Fill(vMMPEmEpGam.M2());
        
        hQsqr_inclusive_gen_acc->Fill(-vQsqr_inc_gen.M2());

        if (abs(vIV_EpEmGam.M() - 0.9577) < 2.5*0.0365 && vEmX.E() > 0.5) { //&& -vQsqr_inc.M2() < 11.5
          hIVEpEmGam_cut->Fill(vIV_EpEmGam.M());
          hMMPEmX_cut->Fill(vMM_PEmX.M());
          hIVEpEm_cut->Fill(vIV_EpEm.M());
          hMM2PEmEpGam_cut->Fill(vMMPEmEpGam.M2());
          hQsqr_inclusive_cut->Fill(-vQsqr_inc.M2());
          hQsqr_inclusive_gen_acccut->Fill(-vQsqr_inc_gen.M2());
          hEm_sTheta_Qsq->Fill(vQsqr_inc_gen.Theta()*180./TMath::Pi(),-vQsqr_inc_gen.M2());
          hEm_sE_Qsq->Fill(vMMPEmEpGam.E(),-vQsqr_inc_gen.M2());
          hEm_sE_Theta->Fill(vMMPEmEpGam.E(),vQsqr_inc_gen.Theta()*180./TMath::Pi());

          hEmP_cut->Fill(vEmX.Vect().Mag());
          hIVEpEm_excl->Fill(IV_EpEm);
          hIVEpEm_Qsq->Fill(vIV_EpEm.M(),-vQsqr_inc.M2());
          hQsqr_inclusive->Fill(-vQsqr_inc.M2());

          
        }
        
      }
      
      if (!PEmEpGam_acc && EmPEpGam_acc) { //scattered electron accepted only
        Double_t Em_s_px = Em_s_P*sin(Em_s_Theta)*cos(Em_s_Phi);
        Double_t Em_s_py = Em_s_P*sin(Em_s_Theta)*sin(Em_s_Phi);
        Double_t Em_s_pz = Em_s_P*cos(Em_s_Theta);
        
        TLorentzVector vEmX, vIV_EpEmGam, vIV_EpEm; //When one electron is detects, both
        TLorentzVector vMM_PEmX;
        TLorentzVector vMMPEmEpGam, vQsqr_inc;
        
        hQsqr_recon->Fill(-Qsq);
        vEmX.SetPxPyPzE(Em_s_px, Em_s_py, Em_s_pz, sqrt(Em_s_P*Em_s_P + M_el*M_el));
        vIV_EpEmGam = (vEp + vEmX + vGamma);
        vIV_EpEm = (vEp + vEmX);
        vMM_PEmX = vEg + vT - (vP + vEmX);
        vMMPEmEpGam = vEg + vT - (vP + vEmX + vEp +vGamma);
        vQsqr_inc = vEg - vMMPEmEpGam;

        hIVEpEmGam->Fill(vIV_EpEmGam.M());
        hIVEpEm->Fill(vIV_EpEm.M());
        
        hMMPEmX->Fill(vMM_PEmX.M());
        hMMPEm2X->Fill(vMM_PEmX.M());
        
        hEmP->Fill(vEmX.Vect().Mag());
        hEm2P->Fill(vEmX.Vect().Mag());
        hMM2PEmEpGam->Fill(vMMPEmEpGam.M2());


        if (abs(vIV_EpEmGam.M() - 0.9577) < 2.5*0.0365 && vEmX.E() > 0.5 ) {//&& -vQsqr_inc.M2() < 11.5
          hQsqr_recon_cut->Fill(-Qsq);

          hIVEpEmGam_cut->Fill(vIV_EpEmGam.M());
          hMMPEmX_cut->Fill(vMM_PEmX.M());
          hIVEpEm_cut->Fill(vIV_EpEm.M());
          hMM2PEmEpGam_cut->Fill(vMMPEmEpGam.M2());
          hQsqr_inclusive_cut->Fill(-vQsqr_inc.M2());

          hEmP_cut->Fill(vEmX.Vect().Mag());
          hIVEpEm_Qsq->Fill(vIV_EpEm.M(),-vQsqr_inc.M2());
          hQsqr_inclusive->Fill(-vQsqr_inc.M2());

        }
        
      }
      /*
      if (PEmEpGam_acc && EmPEpGam_acc) { //Both electrons accepted
        
        TLorentzVector vEm1X, vEm2X, vIV_EpEm1Gam, vIV_EpEm2Gam; //When boths electrons is detects, both

        
        Double_t Em_px = Em_P*sin(Em_Theta)*cos(Em_Phi);
        Double_t Em_py = Em_P*sin(Em_Theta)*sin(Em_Phi);
        Double_t Em_pz = Em_P*cos(Em_Theta);
        
        Double_t Em_s_px = Em_s_P*sin(Em_s_Theta)*cos(Em_s_Phi);
        Double_t Em_s_py = Em_s_P*sin(Em_s_Theta)*sin(Em_s_Phi);
        Double_t Em_s_pz = Em_s_P*cos(Em_s_Theta);
        
        vEm1X.SetPxPyPzE(Em_px, Em_py, Em_pz, sqrt(Em_P*Em_P + M_el*M_el));
        vEm2X.SetPxPyPzE(Em_s_px, Em_s_py, Em_s_pz, sqrt(Em_s_P*Em_s_P + M_el*M_el));
        
        TLorentzVector vMMPEmEpGam_exc, vQsqr_exl, vQsqr_exl_scatter;
        vMMPEmEpGam_exc = vEg + vT - (vP + vEm1X + vEp +vGamma);
        vQsqr_exl = vEg - vMMPEmEpGam_exc;
        vQsqr_exl_scatter = vEg - vEm2X;
        
        hscat_reconE_theta->Fill(vEm2X.Theta()*180./TMath::Pi(),vMMPEmEpGam_exc.Theta()*180./TMath::Pi());
        hscat_reconE_Qsqr->Fill(-vQsqr_exl_scatter.M2(),-vQsqr_exl.M2());

        
        vIV_EpEm1Gam = (vEp + vEm1X + vGamma);
        vIV_EpEm2Gam = (vEp + vEm2X + vGamma);
        
        hIVEpEm1Gam->Fill(vIV_EpEm1Gam.M());
        hIVEpEm2Gam->Fill(vIV_EpEm2Gam.M());
        
        //        hEm1P->Fill(vEm1X.Vect().Mag());
        //        hEm2P->Fill(vEm2X.Vect().Mag());
        
      }
      */
      
    }
    
    
  }
  
  TCanvas *c1 = new TCanvas("c1","c1",800,1000);
  c1->Divide(2,2);
  c1->cd(1);
  hIVEpEmGam->Draw();
  fitMKVoight(hIVEpEmGam, 0.6, 2.1, 0, 0, 0, 0, 0.957, 0.01, 0.001, 2.5, 1);
  c1->cd(2);
  hMMPEmX->Draw();
  fitMKVoight(hMMPEmX, 0.6, 2.1, 0, 0, 0, 0, 0.957, 0.01, 0.001, 2.5, 1);
  c1->cd(3)->SetLogy();
  hIVEpEm->Draw();
  c1->cd(4);
  hEmP->Draw();
  
  TCanvas *c2 = new TCanvas("c2","c2",800,1000);
  c2->Divide(2,2);
  c2->cd(1);
  hIVEpEmGam_cut->Draw();
  //fitMKVoight(hIVEpEmGam, 0.6, 2.1, 0, 0, 0, 0, 0.957, 0.01, 0.001, 2.5, 1);
  c2->cd(2);
  hMMPEmX_cut->Draw();
  //fitMKVoight(hMMPEmX, 0.6, 2.1, 0, 0, 0, 0, 0.957, 0.01, 0.001, 2.5, 1);
  c2->cd(3)->SetLogy();
  hIVEpEm_cut->Draw();
  c2->cd(4);
  hEmP_cut->Draw();
  
  TCanvas *c3 = new TCanvas("c3","c3",800,1000);
  c3->Divide(2,2);
  c3->cd(1);
  TH1D *hEpEm_contam = new TH1D("hEpEm_contam","hEpEm_contam",100,0.0,1);
  //hEpEm_contam->Sumw2();
  hEpEm_contam->Add(hIVEpEm_excl,hIVEpEm_cut,-1,1);
  c3->cd(1);
  hEpEm_contam->SetDirectory(0);
  hEpEm_contam->Draw("EP");
  c3->cd(2);
  TH1D *hEpEm_acceptance = new TH1D("hEpEm_acceptance","hEpEm_acceptance",100,0.0,1);
  hEpEm_acceptance->Sumw2();
  hEpEm_acceptance->Divide(hIVEpEm_cut,hIVEpEm_gen);
  hEpEm_acceptance->Draw("EP");
  c3->cd(3);
  TH2D *hIVEpEm_Qsq_gen_acceptance = new TH2D("hIVEpEm_Qsq_gen_acceptance","hIVEpEm_Qsq_gen_acceptance",100,0.0,1,90,0,18);
  hIVEpEm_Qsq_gen_acceptance->Sumw2();
  hIVEpEm_Qsq_gen_acceptance->Divide(hIVEpEm_Qsq,hIVEpEm_Qsq_gen);
  hIVEpEm_Qsq_gen_acceptance->Draw("colz");
  
  
  TCanvas *c4 = new TCanvas("c4","c4",800,1000);
  c4->Divide(2,2);
  c4->cd(1);
  hMM2PEmEpGam->Draw();
  c4->cd(2);
  hMM2PEmEpGam_cut->Draw();
  c4->cd(3);
  hQsqr_inclusive->Draw();
  c4->cd(4);
  hQsqr_inclusive_cut->Draw();
  
  
  
  TCanvas *c5 = new TCanvas("c5","c5",800,1000);
  c5->cd(1);

  int dannybin_1 = hIVEpEm_Qsq_gen_acceptance->GetYaxis()->FindFixBin(0.5);
  int dannybin_2 = hIVEpEm_Qsq_gen_acceptance->GetYaxis()->FindFixBin(3.5);
  int dannybin_3 = hIVEpEm_Qsq_gen_acceptance->GetYaxis()->FindFixBin(7.5);
  int dannybin_4 = hIVEpEm_Qsq_gen_acceptance->GetYaxis()->FindFixBin(11.5);
  int dannybin_5 = hIVEpEm_Qsq_gen_acceptance->GetYaxis()->FindFixBin(15.5);
  
  TH1D *hsameacc = (TH1D*)hIVEpEm_Qsq_gen_acceptance->ProjectionX("hsameacc",dannybin_1,dannybin_2,"");
  TH1D *hsameacc_2 = (TH1D*)hIVEpEm_Qsq_gen_acceptance->ProjectionX("hsameacc_2",dannybin_2,dannybin_3,"");
  TH1D *hsameacc_3 = (TH1D*)hIVEpEm_Qsq_gen_acceptance->ProjectionX("hsameacc_3",dannybin_3,dannybin_4,"");
  TH1D *hsameacc_4 = (TH1D*)hIVEpEm_Qsq_gen_acceptance->ProjectionX("hsameacc_4",dannybin_4,dannybin_5,"");

  hsameacc->Draw();
  hsameacc_2->Draw("same");
  hsameacc_3->Draw("same");
  hsameacc_4->Draw("same");
  
  
  TCanvas *c6 = new TCanvas("c6","c6",800,1000);

  c6->Divide(2,2);
  c6->cd(1);
  hMMPEmEpGam_gen->Draw("");
  c6->cd(2);
  hscat_reconE_theta->Draw("colz");
  c6->cd(3);
  hEm_sE_Qsq->Draw("colz");
  c6->cd(4);
  hEm_sE_Theta->Draw("colz");
//  hQsqr_inclusive_gen_acc->Draw();
//  c5->cd(3);
//  hQsqr_inclusive_gen_acccut->Draw();
//  c5->cd(4);
//  hEm_sTheta_Qsq->Draw("colz");
  //hQsqr_inclusive_gen->Draw();
  
  //
  
  //Lets now calculate the number of expected results
  Double_t etaP_eegam_BR = 5.13e-04; //  BR from PRELIM CLAS result
  
  Double_t time = 6.912e06; //this is in seconds 80 days
  //Derek
  Double_t total_etaP_expected = 2.05897e+09;//5.0e07;
  Double_t total_etaP_Dalitzexpected = total_etaP_expected * etaP_eegam_BR;
  cout<<total_etaP_Dalitzexpected<<" total_etaP_Dalitzexpected / total_etaP_Dalitzexpected with aceptance--> "<<total_etaP_Dalitzexpected*hEpEm_acceptance->Integral("width")<<endl;
  //
  
  Double_t InclIntegrated_Events_upper = 0.0;
  Double_t spread_events = total_etaP_Dalitzexpected/hIVEpEm_cut->GetNbinsX();//number seen per bin
  
  for (int i = 1; i<hEpEm_acceptance->GetNbinsX(); i++) {
    Double_t intotal_acceptance = hEpEm_acceptance->GetBinContent(i);
    Double_t intotal_events_upper;
    if (intotal_acceptance == 0) {
      intotal_events_upper = 0;
    }else{
      intotal_events_upper = (spread_events*intotal_acceptance);  //THIS IS HOW Ntot was calculated
    }
    InclIntegrated_Events_upper = InclIntegrated_Events_upper + intotal_events_upper;
  }
  cout<<InclIntegrated_Events_upper<<" Events Expected  52134.4 from previous production analysis"<<endl;
  //Done counting
  //Now lets make plots of what the spectrum will look like with and without acceptance and QED normalized
  //Before acceptance
  Double_t Nincl_accepted = hIVEpEm_cut->GetEntries();
  TH1D *hIVEpEm_cut_clone = (TH1D*)hIVEpEm_cut->Clone();
  hIVEpEm_cut_clone->SetName("hIVEpEm_cut_clone");
  hIVEpEm_cut_clone->Scale(InclIntegrated_Events_upper/Nincl_accepted);
  Double_t clone_sum = 0.0;
  cout<<" IN clone "<<hIVEpEm_cut_clone->Integral()<<endl;
  //Now lets normalozed by QED
  //For QED normalization
  //
  
  double QED_par[8] = {0.957}; //{Mass}
  TF1 *QED_norm = new TF1("QED_norm",Eval_Kroll_wada,0.,1.,1);
  QED_norm->SetParameters(&QED_par[0]);
  
  TH1D *hEpEm_corrected= new TH1D("hEpEm_corrected","hEpEm_corrected",100,0.0,1);
  TH1D *hEpEm_QEDnorm= new TH1D("hEpEm_QEDnorm","hEpEm_QEDnorm",100,0.0,1);
  //hEpEm_QEDnorm->Sumw2();
  
  for (int i = 1; i<hIVEpEm_cut_clone->GetNbinsX(); i++) {
    
    //first I want to make sure the spectrum of the data has the proper errors
    
    hIVEpEm_cut_clone->SetBinError(i,sqrt(hIVEpEm_cut_clone->GetBinContent(i)));
    
    
    Double_t bin_factor = 1.65; //this needs to be solved at somepoint
    Double_t intotal_acceptance = hEpEm_acceptance->GetBinContent(i)*bin_factor;
    clone_sum = clone_sum + hIVEpEm_cut_clone->GetBinContent(i);
    Double_t intotal_events_upper;
    if (intotal_acceptance == 0) {
      intotal_events_upper = 0;
    }else{
      intotal_events_upper = hIVEpEm_cut_clone->GetBinContent(i)/intotal_acceptance;
    }
    hEpEm_corrected->SetBinContent(i,intotal_events_upper);
    hEpEm_corrected->SetBinError(i,sqrt(intotal_events_upper));
    
    Double_t QED_factor = QED_norm->Eval(hIVEpEm_cut_clone->GetBinCenter(i));//2.0e-06;
    
    hEpEm_QEDnorm->SetBinContent(i,intotal_events_upper/QED_factor);
    //hEpEm_QEDnorm->SetBinError(i,sqrt(intotal_events_upper/QED_factor));
    
  }
  cout<<"hEpEm_corrected "<<hEpEm_corrected->Integral()<<endl;
  TCanvas *cQED = new TCanvas("cQED","cQED",800,1000);
  cQED->Divide(1,3);
  cQED->cd(1);
  hIVEpEm_cut_clone->Draw("EP");
  cQED->cd(2);
  hEpEm_corrected->Draw("EP");
  cQED->cd(3);
  hEpEm_QEDnorm->Draw("EP");
  double FF_par[2] = {1.,0.77}; //{Normalization,Lambda}
  //double FF_par[1] = {0.77}; //{Normalization,Lambda}
  TF1 *FF_fitter = new TF1("FF_fitter",VMD_FF,0.2,0.76,2);
  FF_fitter->SetParameters(&FF_par[0]);
  FF_fitter->SetParLimits(0,0.,10.);
  FF_fitter->SetParLimits(1,FF_par[0] - FF_par[0]*0.5,FF_par[0] + FF_par[0]*0.5);
  //hEpEm_QEDnorm->Fit("FF_fitter","REM");

  //other side
  FF_fitter->GetParameters(FF_par);
  //double FF_parII[2] = {1.,0.77}; //{Normalization,Lambda}
  TF1 *FF_fitterII = new TF1("FF_fitterII",VMD_FF,0.775,0.9,2);
  FF_fitterII->SetParameter(0,FF_par[0]);
  FF_fitterII->SetParameter(1,FF_par[1]);
  FF_fitterII->SetParLimits(1,FF_par[0] - FF_par[0]*0.5,FF_par[0] + FF_par[0]*0.5);

  //FF_fitterII->SetParameters(&FF_parII[0]);
  //FF_fitterII->SetParLimits(0,0.,10.);
  //FF_fitterII->SetParLimits(1,FF_parII[0] - FF_parII[0]*0.5,FF_parII[0] + FF_parII[0]*0.5);
  //hEpEm_QEDnorm->Fit("FF_fitterII","REM+");

  //Lets try gauss fit
  double Voight_par[4] = {10, 0.77, 0.01 , 0.001}; //{A,mean,sigma,Gamma}

  TF1 *FF_fitterGaus = new TF1("FF_fitterGaus", myVoight,0.55,0.92,4);
  FF_fitterGaus->SetParameters(&Voight_par[0]);

  //TF1 *FF_fitterGaus = new TF1("FF_fitterGaus","gaus",0.55,0.92);
  //hEpEm_QEDnorm->Fit("FF_fitterGaus","REM");
  
  
  //now lets try BW pole
  
  double pole_par[3] = {10, 0.81,0.12}; //{A,Lambda,Gamma}
  //double pole_par[3] = {10, 0.49,0.0144}; //{A,Lambda,Gamma}

  TF1 *FF_fittepole = new TF1("FF_fittepole", Pole_FF,0.02,0.92,3);
  FF_fittepole->SetParameters(&pole_par[0]);
  FF_fittepole->SetParLimits(0,0.,100.);
  FF_fittepole->SetParLimits(1,pole_par[1] - pole_par[1]*0.01,pole_par[1] + pole_par[1]*0.01);
  FF_fittepole->SetParLimits(2,pole_par[2] - pole_par[2]*0.5,pole_par[2] + pole_par[2]*0.5);

  //TF1 *FF_fitterGaus = new TF1("FF_fitterGaus","gaus",0.55,0.92);
  hEpEm_QEDnorm->Fit("FF_fittepole","REM");
  cout<<hEpEm_acceptance->Integral("width")<<" Acceptance Integral"<<endl;
  
  c1->Print("Incl_plots.pdf(");
  c2->Print("Incl_plots.pdf");
  c3->Print("Incl_plots.pdf)");
  
  //TFile fplot("FLAT/IncLusive_Plots_flat.root","recreate");
  TFile fplot("VMD/IncLusive_Plots_vmdtest.root","recreate");

  
  hIVEpEmGam->Write();
  hMMPEmX->Write();
  hIVEpEm->Write();
  hEmP->Write();
  hMM2PEmEpGam->Write();
  hQsqr_inclusive->Write();
  
  hIVEpEmGam_cut->Write();
  hMMPEmX_cut->Write();
  hIVEpEm_cut->Write();
  hEmP_cut->Write();
  hMM2PEmEpGam_cut->Write();
  hQsqr_inclusive_cut->Write();
  
  hEpEm_contam->Write();
  hEpEm_acceptance->Write();
  hIVEpEm_Qsq_gen->Write();
  hIVEpEm_Qsq->Write();
  hIVEpEm_Qsq_gen_acceptance->Write();
  hIVEpEm_cut_clone->Write();
  hEpEm_corrected->Write();
  hEpEm_QEDnorm->Write();
  
  
  
  fplot.Write();
  fplot.Close();
  
}

void fitMKGaus(TH1 *h33 , Double_t low, Double_t high, Double_t p0, Double_t p1, Double_t p2, Double_t p3, Double_t initialPar, Double_t width, Double_t factor, Int_t draw_opt){
  
  //double nEnt = h33->GetEntries();
  double nEnt = h33->GetMaximum();
  
  TF1 *fitter = new TF1("fitter","gaus + [3] + [4]*x + [5]*x*x + [6]*x*x*x",low,high);
  fitter->SetParameters(100, initialPar, width, p0, p1, p2, p3); fitter->SetParLimits(0,10,nEnt); fitter->SetParLimits(1,initialPar-width,initialPar+width);
  h33->Fit("fitter","REM","same");
  //h33->Draw("E");
  TF1 *backFcn = new TF1("backFcn", "pol3",low,high);
  TF1 *signalFcn = new TF1("signalFcn", "gaus",low,high);
  signalFcn->SetLineColor(2);
  signalFcn->SetLineWidth(2);
  Double_t par[7];
  fitter->GetParameters(par);
  backFcn->SetParameters(&par[3]);
  backFcn->SetLineStyle(2);
  backFcn->SetLineColor(6);
  backFcn->SetLineWidth(1);
  
  
  signalFcn->SetParameters(par);
  signalFcn->SetLineStyle(2);
  signalFcn->SetLineColor(4);
  signalFcn->SetLineWidth(1);
  backFcn->Draw("same");
  signalFcn->Draw("same");
  Double_t Intg = abs(signalFcn->Integral(par[1]-factor*par[2],par[1]+factor*par[2]));
  Double_t Intb = abs(backFcn->Integral(par[1]-factor*par[2],par[1]+factor*par[2]));
  
  
  Double_t binw = h33->GetBinWidth(1);
  Int_t yield = Intg/binw;
  Int_t bckgd = Intb/binw;
  //Double_t ratio = double(yield)/TMath::Sqrt(double(bckgd));
  Double_t ratio = double(yield)/double(yield+bckgd);
  
  cout << yield << "\t" << ratio << endl;
  TAxis *x=h33->GetXaxis();
  TAxis *y=h33->GetYaxis();
  
  Double_t startx=x->GetXmin()+0.75*(x->GetXmax()-x->GetXmin());
  
  Double_t starty0=0.35*h33->GetMaximum();
  Double_t starty1=0.45*h33->GetMaximum();
  Double_t starty2=0.55*h33->GetMaximum();
  Double_t starty3=0.65*h33->GetMaximum();
  Double_t starty4=0.75*h33->GetMaximum();
  Double_t starty5=0.85*h33->GetMaximum();
  
  
  double meanError = fitter->GetParError(1);
  double sigmaError = fitter->GetParError(2);
  if (draw_opt ==1) {
    
    TLatex *sum = new TLatex(startx*0.93, starty5,Form("Yield: %i",yield));
    TLatex *sum12 = new TLatex(startx*0.93, starty4,Form("Background: %i",bckgd));
    TLatex *sum0=new TLatex(startx*0.93, starty3,Form("Range: #pm %2.1f #sigma",factor));
    TLatex *sum2=new TLatex(startx*0.93, starty2,Form("Mean:%4.4f #pm %.4f GeV",par[1], meanError));
    TLatex *sum3=new TLatex(startx*0.93, starty1,Form("#sigma:%5.4f #pm %.4f GeV",par[2], sigmaError));
    //TLatex *ra = new TLatex(startx*0.93, starty0,Form("#frac{S}{#sqrt{B}}= %.1f", ratio));
    TLatex *ra = new TLatex(startx*0.93, starty0,Form("#frac{S}{S+B}= %1.1f", ratio));
    
    
    sum->SetTextSize(0.04);
    sum->SetTextColor(2);
    sum->Draw("same");
    sum12->SetTextSize(0.04);
    sum12->SetTextColor(6);
    sum12->Draw("same");
    
    ra->SetTextSize(0.04);
    ra->Draw("same");
    
    
    sum0->SetTextSize(0.04);
    sum0->SetTextColor(2);
    sum0->Draw("same");
    sum2->SetTextSize(0.04);
    sum2->SetTextColor(4);
    sum2->Draw("same");
    sum3->SetTextSize(0.04);
    sum3->SetTextColor(4);
    sum3->Draw("same");
    
    
  }
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
  
  f = (Norm)*(2.0*structure_constant/(3.*pi)*(1./sqrt(qsqr))*sqrt(1 - 4.*(pow(Melectron,2)/qsqr)) * (1. + 2.0 *(pow(Melectron,2)/qsqr)) * (pow((1. - qsqr/(pow(M,2))),3)) ) ;
  
  return f;
}

Double_t VMD_FF(Double_t *x, Double_t *par) {
  
//  if (x[0] > 0.742 && x[0] < 0.779) {
//    TF1::RejectPoint();
//    return 0;
//  }else{
//    
//    Double_t Lambda= par[1]*par[1];
//    Double_t qsqr = (x[0]*x[0]);
//    Double_t f = 0.;
//    Double_t Norm = par[0];
//    Double_t demon = 1. - qsqr/Lambda;
//    f = 1./demon;
//    return Norm*abs(f*f);
//  }
  Double_t Lambda= par[1]*par[1];
  Double_t qsqr = (x[0]*x[0]);
  Double_t f = 0.;
  Double_t Norm = par[0];
  Double_t demon = 1. - qsqr/Lambda;
  f = 1./demon;
  return Norm*abs(f*f);
  //return Norm*abs(f*f);
  
}
Double_t Pole_FF(Double_t *x, Double_t *par) {

  Double_t Lambda= par[1]*par[1];
  Double_t gamma = par[2];
  Double_t qsqr = (x[0]*x[0]);
  Double_t f = 0.;
  Double_t Norm = par[0];
  
  Double_t nom = Lambda*Lambda*(Lambda*Lambda + gamma*gamma);
  Double_t demon = pow((Lambda*Lambda-qsqr),2) + Lambda*Lambda*gamma*gamma;
  f = nom/demon;
  return Norm*f;
  //return Norm*abs(f*f);
  
  
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

void fitMKVoight(TH1 *h33 , Double_t low, Double_t high, Double_t p0, Double_t p1, Double_t p2, Double_t p3, Double_t initialPar, Double_t width, Double_t Gamma, Double_t factor, Int_t draw_opt){
  
  //double nEnt = h33->GetEntries();
  double nEnt = h33->GetMaximum();
  double in_par[8] = {100, initialPar, width , Gamma, p0, p1, p2}; //{A,mean,sigma,Gamma}
  TF1 *fitter = new TF1("fitter",myFuncBWG,low,high,8);
  fitter->SetParameters(&in_par[0]);
  fitter->SetParLimits(0,10,nEnt); fitter->SetParLimits(1,initialPar-width,initialPar+width);
  h33->Fit("fitter","REM","same");
  //h33->Draw("E");
  TF1 *backFcn = new TF1("backFcn", "pol2",low,high);
  TF1 *signalFcn = new TF1("signalFcn", myVoight,low,high,4);
  signalFcn->SetLineColor(2);
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
  backFcn->Draw("same");
  signalFcn->Draw("same");
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
  
  Double_t starty0=0.5*h33->GetMaximum();
  Double_t starty1=0.56*h33->GetMaximum();
  Double_t starty2=0.62*h33->GetMaximum();
  Double_t starty3=0.68*h33->GetMaximum();
  Double_t starty4=0.74*h33->GetMaximum();
  Double_t starty5=0.8*h33->GetMaximum();
  
  
  double meanError = fitter->GetParError(1);
  double sigmaError = fitter->GetParError(2);
  if (draw_opt ==1) {
    
    TLatex *sum = new TLatex(startx*0.93, starty5,Form("Yield: %i",yield));
    TLatex *sum12 = new TLatex(startx*0.93, starty4,Form("Background: %i",bckgd));
    TLatex *sum0=new TLatex(startx*0.93, starty3,Form("Range: #pm %2.1f #sigma",factor));
    TLatex *sum2=new TLatex(startx*0.93, starty2,Form("Mean:%4.4f #pm %.4f GeV",par[1], meanError));
    TLatex *sum3=new TLatex(startx*0.93, starty1,Form("#sigma:%5.4f #pm %.4f GeV",par[2], sigmaError));
    //TLatex *ra = new TLatex(startx*0.93, starty0,Form("#frac{S}{#sqrt{B}}= %.1f", ratio));
    TLatex *ra = new TLatex(startx*0.93, starty0,Form("#frac{S}{S+B}= %1.4f", ratio));
    
    
    sum->SetTextSize(0.04);
    sum->SetTextColor(2);
    sum->Draw("same");
    sum12->SetTextSize(0.04);
    sum12->SetTextColor(6);
    sum12->Draw("same");
    
    ra->SetTextSize(0.04);
    ra->Draw("same");
    
    
    sum0->SetTextSize(0.04);
    sum0->SetTextColor(2);
    sum0->Draw("same");
    sum2->SetTextSize(0.04);
    sum2->SetTextColor(4);
    sum2->Draw("same");
    sum3->SetTextSize(0.04);
    sum3->SetTextColor(4);
    sum3->Draw("same");
    
    
  }
}


TString XSection_Hist_num(Double_t sqrt_s){
  
  //if(sqrt_s > 1.6 && sqrt_s <=1.91 ){ return "hist1";}
  if(sqrt_s <=1.91 ){ return "hist1";}
  else if(sqrt_s > 1.91 && sqrt_s <=1.92 ){ return "hist2";}
  else if(sqrt_s > 1.92 && sqrt_s <=1.93 ){ return "hist3";}
  else if(sqrt_s > 1.93 && sqrt_s <=1.94 ){ return "hist4";}
  else if(sqrt_s > 1.94 && sqrt_s <=1.95 ){ return "hist5";}
  else if(sqrt_s > 1.95 && sqrt_s <=1.96 ){ return "hist6";}
  else if(sqrt_s > 1.96 && sqrt_s <=1.97 ){ return "hist7";}
  else if(sqrt_s > 1.97 && sqrt_s <=1.98 ){ return "hist8";}
  else if(sqrt_s > 1.98 && sqrt_s <=1.99 ){ return "hist9";}
  else if(sqrt_s > 1.99 && sqrt_s <=2 ){ return "hist10";}
  else if(sqrt_s > 2 && sqrt_s <=2.01 ){ return "hist11";}
  else if(sqrt_s > 2.01 && sqrt_s <=2.02 ){ return "hist12";}
  else if(sqrt_s > 2.02 && sqrt_s <=2.03 ){ return "hist13";}
  else if(sqrt_s > 2.03 && sqrt_s <=2.04 ){ return "hist14";}
  else if(sqrt_s > 2.04 && sqrt_s <=2.05 ){ return "hist15";}
  else if(sqrt_s > 2.05 && sqrt_s <=2.06 ){ return "hist16";}
  else if(sqrt_s > 2.06 && sqrt_s <=2.07 ){ return "hist17";}
  else if(sqrt_s > 2.07 && sqrt_s <=2.08 ){ return "hist18";}
  
  
  else if(sqrt_s > 2.08 && sqrt_s <=2.09 ){ return "hist19";}
  else if(sqrt_s > 2.09 && sqrt_s <=2.1 ){ return "hist20";}
  
  else if(sqrt_s > 2.1 && sqrt_s <=2.12 ){ return "hist21";}
  else if(sqrt_s > 2.12 && sqrt_s <=2.14 ){ return "hist22";}
  else if(sqrt_s > 2.14 && sqrt_s <=2.16 ){ return "hist23";}
  else if(sqrt_s > 2.16 && sqrt_s <=2.18 ){ return "hist24";}
  else if(sqrt_s > 2.18 && sqrt_s <=2.2 ){ return "hist25";}
  else if(sqrt_s > 2.2 && sqrt_s <=2.22 ){ return "hist26";}
  else if(sqrt_s > 2.22 && sqrt_s <=2.24 ){ return "hist27";}
  else if(sqrt_s > 2.24 && sqrt_s <=2.26 ){ return "hist28";}
  else if(sqrt_s > 2.26 && sqrt_s <=2.28 ){ return "hist29";}
  else if(sqrt_s > 2.28 && sqrt_s <=2.3 ){ return "hist30";}
  else if(sqrt_s > 2.3 && sqrt_s <=2.32 ){ return "hist31";}
  else if(sqrt_s > 2.32 && sqrt_s <=2.34 ){ return "hist32";}
  else if(sqrt_s > 2.34 && sqrt_s <=2.36 ){ return "hist33";}
  else if(sqrt_s > 2.36 && sqrt_s <=2.4 ){ return "hist34";}
  else if(sqrt_s > 2.4 && sqrt_s <=2.44 ){ return "hist35";}
  else if(sqrt_s > 2.44 && sqrt_s <=2.48 ){ return "hist36";}
  else if(sqrt_s > 2.48 && sqrt_s <=2.52 ){ return "hist37";}
  else if(sqrt_s > 2.52 && sqrt_s <=2.56 ){ return "hist38";}
  else if(sqrt_s > 2.56 && sqrt_s <=2.6 ){ return "hist39";}
  else if(sqrt_s > 2.6 && sqrt_s <=2.64 ){ return "hist40";}
  else if(sqrt_s > 2.64 && sqrt_s <=2.68 ){ return "hist41";}
  else if(sqrt_s > 2.68 && sqrt_s <=2.73 ){ return "hist42";}
  else if(sqrt_s > 2.73){ return "hist42";}
  else {
    //cout<<"NOT IN BOUNDS WITH z VERTEX"<<endl;
    return "";
  }
  //cout<<"I am at the end "<<z<<"  "<<p<<endl;
  
  
}






