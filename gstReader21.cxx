#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TMath.h>
#include <TCanvas.h>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

void gstReader21()
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
    Root_Tree->SetBranchStatus("nip",  1);
    Root_Tree->SetBranchStatus("nin",  1);
    Root_Tree->SetBranchStatus("nfp",  1);
    Root_Tree->SetBranchStatus("nfn",  1);

    double cthl = 0.0;
    int    nf   = 0;
    int    nip, nin, nfp, nfn;
    const  int kMax = 200;
    int    pdgf[kMax];
    double pf[kMax];

    Root_Tree->SetBranchAddress("cthl", &cthl);
    Root_Tree->SetBranchAddress("nf",   &nf);
    Root_Tree->SetBranchAddress("pdgf",  pdgf);
    Root_Tree->SetBranchAddress("pf",    pf);
    Root_Tree->SetBranchAddress("nip",  &nip);
    Root_Tree->SetBranchAddress("nin",  &nin);
    Root_Tree->SetBranchAddress("nfp",  &nfp);
    Root_Tree->SetBranchAddress("nfn",  &nfn);

    // FSI conserved
    TH1D* hKpCons  = new TH1D("hKpCons",  "K^{+} FSI Conserved;p (GeV);Events",  100, 0, 5);
    TH1D* hKmCons  = new TH1D("hKmCons",  "K^{-} FSI Conserved;p (GeV);Events",  100, 0, 5);

    // FSI not conserved
    TH1D* hKpNotCons  = new TH1D("hKpNotCons",  "K^{+} FSI Not Conserved;p (GeV);Events",  100, 0, 5);
    TH1D* hKmNotCons  = new TH1D("hKmNotCons",  "K^{-} FSI Not Conserved;p (GeV);Events",  100, 0, 5);

    // FSI subtypes — proton absorbed, neutron absorbed, charge exchange, knockout
    TH1D* hKpAbsorb   = new TH1D("hKpAbsorb",   "K^{+} Nucleon Absorbed;p (GeV);Events",   100, 0, 5);
    TH1D* hKpExchange = new TH1D("hKpExchange",  "K^{+} Charge Exchange;p (GeV);Events",   100, 0, 5);
    TH1D* hKpKnockout = new TH1D("hKpKnockout",  "K^{+} Nucleon Knockout;p (GeV);Events",  100, 0, 5);

    TH1D* hKmAbsorb   = new TH1D("hKmAbsorb",   "K^{-} Nucleon Absorbed;p (GeV);Events",   100, 0, 5);
    TH1D* hKmExchange = new TH1D("hKmExchange",  "K^{-} Charge Exchange;p (GeV);Events",   100, 0, 5);
    TH1D* hKmKnockout = new TH1D("hKmKnockout",  "K^{-} Nucleon Knockout;p (GeV);Events",  100, 0, 5);

    Long64_t nEntries = Root_Tree->GetEntries();

    for (Long64_t i = 0; i < nEntries; i++)
    {
        if (i % 1000000 == 0)
            cout << "Processing event " << i << " / " << nEntries << endl;

        Root_Tree->GetEntry(i);
        if (cthl < -1 || cthl > 1) continue;

        bool conserved  = (nip == nfp && nin == nfn);
        bool exchange   = (!conserved && (nip+nin) == (nfp+nfn));
        bool knockout   = (nfp+nfn > nip+nin);
        bool absorbed   = (nfp+nfn < nip+nin);

        for (int j = 0; j < nf; j++)
        {
            double p = pf[j];
            if (p <= 0) continue;

            if (pdgf[j] == 321)
            {
                if (conserved)       hKpCons->Fill(p);
                else                 hKpNotCons->Fill(p);
                if (absorbed)        hKpAbsorb->Fill(p);
                else if (exchange)   hKpExchange->Fill(p);
                else if (knockout)   hKpKnockout->Fill(p);
            }
            else if (pdgf[j] == -321)
            {
                if (conserved)       hKmCons->Fill(p);
                else                 hKmNotCons->Fill(p);
                if (absorbed)        hKmAbsorb->Fill(p);
                else if (exchange)   hKmExchange->Fill(p);
                else if (knockout)   hKmKnockout->Fill(p);
            }
        }
    }

    cout << "Drawing..." << endl;
    gStyle->SetOptStat(0);

    vector<int> colors = {kBlack, kRed, kBlue, kGreen+2, kOrange+7};

    auto DrawOverlay = [&](vector<TH1D*> hists, vector<TString> titles, TString canvasName, TString saveName)
    {
        TCanvas* c = new TCanvas(canvasName, canvasName, 800, 600);
        gPad->SetLogy();
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);

        double maxVal = 0;
        for (auto h : hists)
            maxVal = max(maxVal, h->GetMaximum());

        TLegend* leg = new TLegend(0.55, 0.6, 0.95, 0.9);
        leg->SetBorderSize(1);
        leg->SetTextSize(0.03);

        for (int i = 0; i < (int)hists.size(); i++)
        {
            TH1D* hc = (TH1D*)hists[i]->Clone(TString::Format("clone_%s_%d", canvasName.Data(), i));
            hc->SetLineColor(colors[i]);
            hc->SetLineWidth(2);
            hc->SetMaximum(maxVal * 10);
            hc->SetMinimum(1);
            if (i == 0) hc->Draw("hist");
            else        hc->Draw("hist same");
            leg->AddEntry(hc, titles[i], "l");
        }
        leg->Draw();
        c->SaveAs(saveName);
    };

    // K+ conserved vs not conserved
    DrawOverlay(
        {hKpCons, hKpNotCons},
        {"Conserved", "Not Conserved"},
        "cKpConsVsNot", "Kp_cons_vs_not.png"
    );

    // K- conserved vs not conserved
    DrawOverlay(
        {hKmCons, hKmNotCons},
        {"Conserved", "Not Conserved"},
        "cKmConsVsNot", "Km_cons_vs_not.png"
    );

    // K+ FSI subtypes
    DrawOverlay(
        {hKpCons, hKpAbsorb, hKpExchange, hKpKnockout},
        {"Conserved", "Absorbed", "Charge Exchange", "Knockout"},
        "cKpSubtypes", "Kp_FSI_subtypes.png"
    );

    // K- FSI subtypes
    DrawOverlay(
        {hKmCons, hKmAbsorb, hKmExchange, hKmKnockout},
        {"Conserved", "Absorbed", "Charge Exchange", "Knockout"},
        "cKmSubtypes", "Km_FSI_subtypes.png"
    );

    cout << "Done!" << endl;
    Root_File->Close();
}