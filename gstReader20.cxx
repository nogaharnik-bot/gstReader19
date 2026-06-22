#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TMath.h>
#include <TCanvas.h>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

void gstReader20()
{
    string File_Name       = "11_1000060120_EM_5_98636GeV_G18_02a_00_000_Q2_0_4.gst.root";
    string File_TTree_Name = "gst";

    cout << "Opening file..." << endl;
    TFile* Root_File = TFile::Open(File_Name.c_str());
    if (!Root_File || Root_File->IsZombie()) { cout << "File not found!" << endl; return; }

    TTree* Root_Tree = (TTree*)Root_File->Get(File_TTree_Name.c_str());
    if (!Root_Tree) { cout << "Tree not found!" << endl; return; }
    cout << "Entries: " << Root_Tree->GetEntries() << endl;

    Root_Tree->SetBranchStatus("*",    0);
    Root_Tree->SetBranchStatus("cthl", 1);
    Root_Tree->SetBranchStatus("nf",   1);
    Root_Tree->SetBranchStatus("pdgf", 1);
    Root_Tree->SetBranchStatus("pf",   1);
    Root_Tree->SetBranchStatus("pxf",  1);
    Root_Tree->SetBranchStatus("pyf",  1);
    Root_Tree->SetBranchStatus("pzf",  1);

    double cthl = 0.0;
    int    nf   = 0;
    const  int kMax = 200;
    int    pdgf[kMax];
    double pf[kMax], pxf[kMax], pyf[kMax], pzf[kMax];

    Root_Tree->SetBranchAddress("cthl", &cthl);
    Root_Tree->SetBranchAddress("nf",   &nf);
    Root_Tree->SetBranchAddress("pdgf",  pdgf);
    Root_Tree->SetBranchAddress("pf",    pf);
    Root_Tree->SetBranchAddress("pxf",   pxf);
    Root_Tree->SetBranchAddress("pyf",   pyf);
    Root_Tree->SetBranchAddress("pzf",   pzf);

    double tMax = TMath::Pi();
    double pMax = TMath::Pi();

    TH2D* hKp  = new TH2D("hKp",  "K^{+};#theta (rad);#phi (rad)",   100, 0, tMax, 100, -pMax, pMax);
    TH2D* hKm  = new TH2D("hKm",  "K^{-};#theta (rad);#phi (rad)",   100, 0, tMax, 100, -pMax, pMax);
    TH2D* hPip = new TH2D("hPip", "#pi^{+};#theta (rad);#phi (rad)", 100, 0, tMax, 100, -pMax, pMax);
    TH2D* hPim = new TH2D("hPim", "#pi^{-};#theta (rad);#phi (rad)", 100, 0, tMax, 100, -pMax, pMax);

    Long64_t nEntries = Root_Tree->GetEntries();

    for (Long64_t i = 0; i < nEntries; i++)
    {
        if (i % 1000000 == 0)
            cout << "Processing event " << i << " / " << nEntries << endl;

        Root_Tree->GetEntry(i);
        if (cthl < -1 || cthl > 1) continue;

        for (int j = 0; j < nf; j++)
        {
            double p = pf[j];
            if (p <= 0) continue;

            double theta = TMath::ACos(pzf[j] / p);
            double phi   = TMath::ATan2(pyf[j], pxf[j]);

            if      (pdgf[j] ==  321) hKp->Fill(theta, phi);
            else if (pdgf[j] == -321) hKm->Fill(theta, phi);
            else if (pdgf[j] ==  211) hPip->Fill(theta, phi);
            else if (pdgf[j] == -211) hPim->Fill(theta, phi);
        }
    }

    cout << "Drawing..." << endl;

    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    TCanvas* c = new TCanvas("c", "#phi vs #theta", 1600, 800);
    c->Divide(2, 2);

    vector<TH2D*> hists = {hKp, hKm, hPip, hPim};

    for (int i = 0; i < 4; i++)
    {
        c->cd(i+1);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);
        gPad->SetRightMargin(0.15);
        hists[i]->Draw("colz");
    }

    c->SaveAs("phi_vs_theta.png");
    cout << "Done!" << endl;

    Root_File->Close();
}