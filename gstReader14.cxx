#include <TFile.h>
#include <TTree.h>
#include <TLegend.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <map>
#include <iostream>
#include <vector>
#include <string>

using namespace std;

void gstReader14()
{
    string File_Path          = "";
    string File_Name          = "11_1000060120_EM_5_98636GeV_G18_02a_00_000_Q2_0_4.gst.root";
    string File_TTree_Name    = "gst";
    string File_Path_and_Name = File_Path + File_Name;

    cout << "Opening file..." << endl;
    TFile* Root_File = TFile::Open(File_Path_and_Name.c_str());
    if (!Root_File || Root_File->IsZombie()) { cout << "File not found!" << endl; return; }
    cout << "File opened" << endl;

    TTree* Root_Tree = (TTree*)Root_File->Get(File_TTree_Name.c_str());
    if (!Root_Tree) { cout << "Tree not found!" << endl; return; }
    cout << "Tree found, entries: " << Root_Tree->GetEntries() << endl;

    Root_Tree->SetBranchStatus("*",    0);
    Root_Tree->SetBranchStatus("Ev",   1);
    Root_Tree->SetBranchStatus("El",   1);
    Root_Tree->SetBranchStatus("cthl", 1);
    Root_Tree->SetBranchStatus("nf",   1);
    Root_Tree->SetBranchStatus("pdgf", 1);
    Root_Tree->SetBranchStatus("pf",   1);
    Root_Tree->SetBranchStatus("pxf",  1);
    Root_Tree->SetBranchStatus("pyf",  1);
    Root_Tree->SetBranchStatus("pzf",  1);

    int    nf = 0;
    const  int kMax = 200;
    int    pdgf[kMax];
    double pf[kMax], pxf[kMax], pyf[kMax], pzf[kMax];
    double gst_Ev = 0.0, gst_El = 0.0, cthl = 0.0;

    Root_Tree->SetBranchAddress("Ev",   &gst_Ev);
    Root_Tree->SetBranchAddress("El",   &gst_El);
    Root_Tree->SetBranchAddress("cthl", &cthl);
    Root_Tree->SetBranchAddress("nf",   &nf);
    Root_Tree->SetBranchAddress("pdgf",  pdgf);
    Root_Tree->SetBranchAddress("pf",    pf);
    Root_Tree->SetBranchAddress("pxf",   pxf);
    Root_Tree->SetBranchAddress("pyf",   pyf);
    Root_Tree->SetBranchAddress("pzf",   pzf);

    // momentum histograms
    TH1D* hKpMom      = new TH1D("hKpMom",      "K^{+} Momentum;p (GeV);Events",       100, 0, 5);
    TH1D* hKmMom      = new TH1D("hKmMom",      "K^{-} Momentum;p (GeV);Events",       100, 0, 5);
    TH1D* hPipMom     = new TH1D("hPipMom",     "#pi^{+} Momentum;p (GeV);Events",     100, 0, 5);
    TH1D* hPimMom     = new TH1D("hPimMom",     "#pi^{-} Momentum;p (GeV);Events",     100, 0, 5);
    TH1D* hProtonMom  = new TH1D("hProtonMom",  "Proton Momentum;p (GeV);Events",      100, 0, 5);
    TH1D* hNeutronMom = new TH1D("hNeutronMom", "Neutron Momentum;p (GeV);Events",     100, 0, 5);

    // angle histograms
    TH1D* hKpTheta      = new TH1D("hKpTheta",      "K^{+} Angle;#theta (rad);Events",       100, 0, TMath::Pi());
    TH1D* hKmTheta      = new TH1D("hKmTheta",      "K^{-} Angle;#theta (rad);Events",       100, 0, TMath::Pi());
    TH1D* hPipTheta     = new TH1D("hPipTheta",     "#pi^{+} Angle;#theta (rad);Events",     100, 0, TMath::Pi());
    TH1D* hPimTheta     = new TH1D("hPimTheta",     "#pi^{-} Angle;#theta (rad);Events",     100, 0, TMath::Pi());
    TH1D* hProtonTheta  = new TH1D("hProtonTheta",  "Proton Angle;#theta (rad);Events",      100, 0, TMath::Pi());
    TH1D* hNeutronTheta = new TH1D("hNeutronTheta", "Neutron Angle;#theta (rad);Events",     100, 0, TMath::Pi());

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

            if      (pdgf[j] == 321)  { hKpMom->Fill(p);      hKpTheta->Fill(theta);      }
            else if (pdgf[j] == -321) { hKmMom->Fill(p);      hKmTheta->Fill(theta);      }
            else if (pdgf[j] == 211)  { hPipMom->Fill(p);     hPipTheta->Fill(theta);     }
            else if (pdgf[j] == -211) { hPimMom->Fill(p);     hPimTheta->Fill(theta);     }
            else if (pdgf[j] == 2212) { hProtonMom->Fill(p);  hProtonTheta->Fill(theta);  }
            else if (pdgf[j] == 2112) { hNeutronMom->Fill(p); hNeutronTheta->Fill(theta); }
        }
    }

    cout << "Event loop done, drawing..." << endl;

    gStyle->SetOptStat(0);

    vector<int> colors = {kBlack, kRed, kBlue, kGreen+2, kViolet+1, kCyan+2};

    // momentum canvas
    TCanvas* cMom = new TCanvas("cMom", "Particle Momenta", 1600, 1200);
    cMom->Divide(3, 2);

    vector<TH1D*> momHists = {hKpMom, hKmMom, hPipMom, hPimMom, hProtonMom, hNeutronMom};
    vector<TString> momTitles = {"K^{+}", "K^{-}", "#pi^{+}", "#pi^{-}", "Proton", "Neutron"};

    for (int i = 0; i < (int)momHists.size(); i++)
    {
        cMom->cd(i+1);
        gPad->SetLogy();
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);
        momHists[i]->SetLineColor(colors[i]);
        momHists[i]->SetLineWidth(2);
        momHists[i]->Draw("hist");

        TLegend* leg = new TLegend(0.6, 0.7, 0.9, 0.9);
        leg->AddEntry(momHists[i], momTitles[i], "l");
        leg->SetBorderSize(0);
        leg->SetTextSize(0.04);
        leg->Draw();
    }
    cMom->SaveAs("particle_momenta.png");

    // angle canvas
    TCanvas* cTheta = new TCanvas("cTheta", "Particle Theta", 1600, 1200);
    cTheta->Divide(3, 2);

    vector<TH1D*> thetaHists = {hKpTheta, hKmTheta, hPipTheta, hPimTheta, hProtonTheta, hNeutronTheta};

    for (int i = 0; i < (int)thetaHists.size(); i++)
    {
        cTheta->cd(i+1);
        gPad->SetLogy();
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);
        thetaHists[i]->SetLineColor(colors[i]);
        thetaHists[i]->SetLineWidth(2);
        thetaHists[i]->Draw("hist");

        TLegend* leg = new TLegend(0.6, 0.7, 0.9, 0.9);
        leg->AddEntry(thetaHists[i], momTitles[i], "l");
        leg->SetBorderSize(0);
        leg->SetTextSize(0.04);
        leg->Draw();
    }
    cTheta->SaveAs("particle_theta.png");

    // overlay momentum canvas
    TCanvas* cMomOverlay = new TCanvas("cMomOverlay", "Momenta Overlay", 800, 600);
    gPad->SetLogy();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);

    TLegend* legMom = new TLegend(0.7, 0.5, 0.95, 0.9);
    legMom->SetBorderSize(1);
    legMom->SetTextSize(0.03);

double maxMom = 0;                                          
for (int i = 0; i < (int)momHists.size(); i++)             
    maxMom = max(maxMom, momHists[i]->GetMaximum()); 

for (int i = 0; i < (int)momHists.size(); i++)    
    {
        TH1D* hc = (TH1D*)momHists[i]->Clone(TString::Format("momclone_%d", i));
        hc->SetLineColor(colors[i]);
        hc->SetLineWidth(2);
        hc->SetMaximum(maxMom * 10);   
        hc->SetMinimum(1);
	hc->SetTitle("Momentum by Particle");
        if (i == 0) hc->Draw("hist");
        else        hc->Draw("hist same");
        legMom->AddEntry(hc, momTitles[i], "l");
    }
    legMom->Draw();
    cMomOverlay->SaveAs("momenta_overlay.png");

    // overlay angle canvas
    TCanvas* cThetaOverlay = new TCanvas("cThetaOverlay", "Theta Overlay", 800, 600);
    gPad->SetLogy();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);

    TLegend* legTheta = new TLegend(0.7, 0.5, 0.95, 0.9);
    legTheta->SetBorderSize(1);
    legTheta->SetTextSize(0.03);

double maxThet = 0;                                          
for (int i = 0; i < (int)momHists.size(); i++)             
    maxThet = max(maxThet, thetaHists[i]->GetMaximum());
for (int i = 0; i < (int)thetaHists.size(); i++)
    {
        TH1D* hc = (TH1D*)thetaHists[i]->Clone(TString::Format("thetaclone_%d", i));
        hc->SetLineColor(colors[i]);
        hc->SetLineWidth(2);
	hc->SetMaximum(maxThet * 10);
	hc->SetMinimum(1);
        hc->SetTitle("Theta by Particle");
	if (i == 0) hc->Draw("hist");
        else        hc->Draw("hist same");
        legTheta->AddEntry(hc, momTitles[i], "l");
    }
    legTheta->Draw();
    cThetaOverlay->SaveAs("theta_overlay.png");

    cout << "Done!" << endl;

    Root_File->Close();
}
