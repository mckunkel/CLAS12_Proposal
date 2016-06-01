

TH1D *hIVEPEm = (TH1D*)in->Get("hIVEpEm");
TH1D *hIVEPEm_cut = (TH1D*)in->Get("hIVEpEm_cut");
TH1D *hEpEm_contam = (TH1D*)in->Get("hEpEm_contam");

hIVEpEm->GetXaxis()->SetTitle("M(e^{+}e^{-}) [GeV]");
hIVEpEm->GetXaxis()->SetTitleSize(0.05);
hIVEpEm->GetXaxis()->SetTitleOffset(0.8);

hIVEpEm->GetYaxis()->SetTitle("Counts");
hIVEpEm->GetYaxis()->SetTitleSize(0.05);
hIVEpEm->GetYaxis()->SetTitleOffset(1.0);

TCanvas *c2 = new TCanvas("c2","",1);
c2->cd();

hIVEpEm->Draw();

