// root -l -b -q gstReader19.cxx

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TLegend.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TString.h>
#include <map>

#include <iostream>
#include <vector>
#include <string>

using namespace std;

struct Topology {
    TString name;
    TString title;
    int nProton;
    int nNeutron;
    int nPip;
    int nPim;
    int nKp;
    int nKm;
};

vector<TH1D*> MakeHists(const TString& name, const TString& title)
{
    return {
        new TH1D("et_" + name,    title +";Energy Transfer [GeV];Events", 100, 0, 6),
        new TH1D("theta_" + name, title +";#theta_{e} [rad];Events",      100, 0, 0.5),
        new TH1D("Q2_" + name,    title +";Q^{2} [GeV^{2}];Events",       100, 0, 1),
        new TH1D("W_" + name,     title +";W [GeV/c^{2}];Events",         100, 0, 3.5)
    };
}

bool MatchTopology(const Topology& tp,
                   int nProton,
                   int nNeutron,
                   int nPip,
                   int nPim,
                   int nKp,
                   int nKm,)
{
    if (tp.nProton >= 0 && nProton != tp.nProton) return false;
    if (tp.nNeutron >= 0 && nNeutron != tp.nNeutron) return false;
    if (tp.nPip >= 0 && nPip != tp.nPip) return false;
    if (tp.nPim >= 0 && nPim != tp.nPim) return false;
    if (tp.nKp >= 0 && nKp != tp.nKp) return false;
    if (tp.nKm >= 0 && nKm != tp.nKm) return false;
    return true;
}

void FillHists(vector<TH1D*>& histos, double et, double theta, double Q2, double W)
{
    histos[0]->Fill(et);
    histos[1]->Fill(theta);
    histos[2]->Fill(Q2);
    histos[3]->Fill(W);
}

vector<vector<TH1D*>> MakeHistoSet(const vector<Topology>& list)
{
    vector<vector<TH1D*>> histos;
    histos.reserve(list.size());
    for (const auto& tp : list)
        histos.push_back(MakeHists(tp.name, tp.title));
    return histos;
}

void DrawTopologies(const vector<vector<TH1D*>>& histos,
                    const vector<Topology>& topos,
                    const TString& prefix,
                    const TString& canvasTitle,
                    int minEntries = 10,
                    bool logY = false,
                    const TString& outputDir = ".")
{
    TString canvasName = prefix + (logY ?"_LOG" :"");
    TCanvas* c = new TCanvas(canvasName, canvasTitle, 1600, 1200);
    c->Divide(2, 2);

    for (int var = 0; var < 4; ++var) {
        c->cd(var + 1);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.15);
        gPad->SetBottomMargin(0.15);
        gPad->SetLogy(logY ? 1 : 0);
        gStyle->SetOptStat(0);

        double maxVal = 0;
        for (size_t t = 0; t < topos.size(); ++t) {
            if (histos[t][var]->GetEntries() < minEntries) continue;
            TH1D* clone = (TH1D*)histos[t][var]->Clone();
            maxVal = max(maxVal, clone->GetMaximum());
            delete clone;
        }

        TLegend* leg = new TLegend(0.87, 0.1, 0.99, 0.9);
        leg->SetTextSize(0.022);
        leg->SetBorderSize(1);

        bool first = true;
        int colorIdx = 0;
        const vector<int> colors = {
            kBlack, kBlue, kRed, kGreen + 2, kOrange + 7,
            kViolet + 1, kCyan + 2, kMagenta + 2, kTeal + 3, kGray + 2
        };
        const vector<int> lineStyles = {1,2,3,4,5,6,7,8,9};

        for (size_t t = 0; t < topos.size(); ++t) {
            if (histos[t][var]->GetEntries() < minEntries) continue;
            TH1D* hc = (TH1D*)histos[t][var]->Clone(TString::Format("clone_%s_%d", topos[t].name.Data(), var));
            hc->SetLineColor(colors[colorIdx % colors.size()]);
            hc->SetLineStyle(lineStyles[colorIdx % lineStyles.size()]);
            hc->SetLineWidth(2);
            hc->SetFillStyle(0);
            hc->GetYaxis()->SetTitle("Events");
            hc->SetTitle(canvasTitle);
            hc->SetMaximum(maxVal * 1.1);
            hc->SetMinimum(logY ? 1 : 0);

            if (first) {
                hc->Draw("hist");
                first = false;
            } else {
                hc->Draw("hist same");
            }
            leg->AddEntry(hc, topos[t].title,"l");
            colorIdx++;
        }

        leg->Draw();
    }

    TString fileSuffix = outputDir + "/" + prefix + (logY ?"_log" :"");
    c->SaveAs(fileSuffix +".png");
}

void gstReader19()
{
    // Create 6GeV output directory
    gSystem->Exec("mkdir -p 6GeV");
    
    const string File_Path = "/pnfs/genie/persistent/users/asportes/2N_Analysis_Samples/C12/GEM21_11a_00_000/5986MeV_Q2_0_40/master-routine_validation_01-eScattering/";
    const string File_TTree_Name = "gst";
    const vector<string> File_Names = {
        "e_on_1000060120_5986MeV_1.gst.root",

        // "e_on_1000060120_5986MeV_*.gst.root"
        // add more ROOT file names here
    };

    TChain* Root_Tree3 = new TChain(File_TTree_Name.c_str());
    int filesAdded = 0;
    for (const auto& fileName : File_Names) {
        const string fullPath = File_Path + fileName;
        filesAdded += Root_Tree3->Add(fullPath.c_str());
    }

    if (filesAdded == 0) {
        cerr << "Error adding files to TChain for tree: " << File_TTree_Name << endl;
        return;
    }

    Root_Tree3->SetBranchStatus("*", 0);
    Root_Tree3->SetBranchStatus("Ev", 1);
    Root_Tree3->SetBranchStatus("El", 1);
    Root_Tree3->SetBranchStatus("cthl", 1);
    Root_Tree3->SetBranchStatus("nf", 1);
    Root_Tree3->SetBranchStatus("pdgf", 1);
    Root_Tree3->SetBranchStatus("pf", 1);
    Root_Tree3->SetBranchStatus("qel", 1);

    int nf = 0;
    const int kMax = 200;
    int pdgf[kMax];
    double pf[kMax];
    double gst_Ev = 0.0;
    double gst_El = 0.0;
    double cthl = 0.0;
    bool qel = false;

    const double kBaryonMomentumThreshold = 0.25; // 250 MeV/c
    const double kMesonMomentumThreshold = 0.15;  // 150 MeV/c

    Root_Tree3->SetBranchAddress("Ev", &gst_Ev);
    Root_Tree3->SetBranchAddress("El", &gst_El);
    Root_Tree3->SetBranchAddress("cthl", &cthl);
    Root_Tree3->SetBranchAddress("nf", &nf);
    Root_Tree3->SetBranchAddress("pdgf", pdgf);
    Root_Tree3->SetBranchAddress("pf", pf);
    Root_Tree3->SetBranchAddress("qel", &qel);

    // 1kp //
    vector<Topology> topos1 = {
        
    };

    // 1p + 1k //
    vector<Topology> topos2 = {
        {"topo_1_1_0_0_0_1_0_0","1p 1n 1K^{+}", 1, 1, 0, 0, 1, 0},
        {"topo_1_0_0_0_0_1_0_0","1p 1K^{+}", 1, 0, 0, 0, 1, 0},
        {"topo_1_1_0_1_0_1_0_0","1p 1n 1#pi^{-} 1K^{+}", 1, 1, 0, 1, 1, 0},
        {"topo_1_1_1_1_0_1_0_0","1p 1n 1#pi^{+} 1#pi^{-} 1K^{+}", 1, 1, 1, 1, 1, 0},
        {"topo_1_0_1_1_0_1_0_0","1p 1#pi^{+} 1#pi^{-} 1K^{+}", 1, 0, 1, 1, 1, 0},
        {"topo_1_2_0_0_0_1_0_0","1p 2n 1K^{+}", 1, 2, 0, 0, 1, 0},
        {"topo_1_0_1_2_0_1_0_0","1p 1#pi^{+} 2#pi^{-}1K^{+}", 1, 0, 1, 2, 1, 0},
        {"topo_1_2_1_1_0_1_0_0","1p 2n 1#pi^{+} 1#pi^{-}1K^{+}", 1, 2, 1, 1, 1, 0},
        {"topo_1_1_1_2_0_1_0_0","1p 1n 1#pi^{+} 2#pi^{-} 1K^{+}", 1, 1, 1, 2, 1, 0},
        {"topo_1_2_0_1_0_1_0_0","1p 2n 1#pi^{-} 1K^{+}", 1, 2, 0, 1, 1, 0},
        {"topo_1_0_0_2_0_1_0_0","1p 2#pi^{-} 1K^{+}", 1, 0, 0, 2, 1, 0},
        {"topo_1_3_0_0_0_1_0_0","1p 3n 1K^{+}", 1, 3, 0, 0, 1, 0},
        {"topo_1_2_1_0_0_1_0_0","1p 2n 1#pi^{+} 1K^{+}", 1, 2, 1, 0, 1, 0},
        {"topo_1_3_1_1_0_1_0_0","1p 3n 1#pi^{+} 1#pi^{-} 1K^{+}", 1, 3, 1, 1, 1, 0},
        {"topo_1_4_0_0_0_1_0_0","1p 4n 1K^{+}", 1, 4, 0, 0, 1, 0},
        {"topo_1_0_1_0_0_0_1_0","1p 1#pi^{+} 0K^{+} 1K^{-}", 1, 0, 1, 0, 0, 1},
    };

    // 1p + 2K//
     vector<Topology> topos3 = {
        {"topo_1_0_0_0_0_1_1_0","1p 1K^{+} 1K^{-}", 1, 0, 0, 0, 1, 1},
        {"topo_1_1_0_0_0_1_1_0","1p 1n 1K^{+} 1K^{-}", 1, 1, 0, 0, 1, 1},
        {"topo_1_0_0_1_0_1_1_0","1p 1#pi^{-} 1K^{+} 1K^{-}", 1, 0, 0, 1, 1, 1},
        {"topo_1_0_1_1_0_1_1_0","1p 1#pi^{+} 1#pi^{-} 1K^{+} 1K^{-}", 1, 0, 1, 1, 1, 1},
        {"topo_1_2_0_0_0_1_1_0","1p 2n 1K^{+} 1K^{-}", 1, 2, 0, 0, 1, 1},
        {"topo_1_1_1_0_0_1_1_0","1p 1n 1#pi^{+} 1K^{+} 1K^{-}", 1, 1, 1, 0, 1, 1},
        {"topo_1_1_1_1_0_1_1_0","1p 1n 1#pi^{+} 1#pi^{-} 1K^{+} 1K^{-}", 1, 1, 1, 1, 1, 1},
        {"topo_1_2_1_0_0_1_1_0","1p 2n 1#pi^{+} 1K^{+} 1K^{-}", 1, 2, 1, 0, 1, 1},
        {"topo_1_1_0_1_0_1_1_0","1p 1n 1#pi^{-} 1K^{+} 1K^{-}", 1, 1, 0, 1, 1, 1},
        {"topo_1_0_1_0_0_0_1_1","1p 1#pi^{+} 1K^{-} ", 1, 0, 1, 0, 0, 1},
        {"topo_1_0_0_0_0_0_1_1","1p 1K^{-}", 1, 0, 0, 0, 0, 1}
    };

    // 1n + 1k //
     vector<Topology> topos4 = {
        {"topo_0_1_0_0_0_1_0_0", "1n 1K^{+}", 0, 1, 0, 0, 1, 0},
        {"topo_0_1_1_1_0_1_0_0", "1n 1#pi^{+} 1#pi^{-} 1K^{+}", 0, 1, 1, 1, 1, 0},
        {"topo_0_1_0_1_0_1_0_0", "1n 1#pi^{-} 1K^{+}", 0, 1, 0, 1, 1, 0},
        {"topo_1_1_0_0_0_1_0_0", "1p 1n 1K^{+}", 1, 1, 0, 0, 1, 0},
        {"topo_1_1_0_1_0_1_0_0", "1p 1n 1#pi^{-} 1K^{+}", 1, 1, 0, 1, 1, 0},
        {"topo_2_1_0_0_0_1_0_0", "2p 1n 1K^{+}", 2, 1, 0, 0, 1, 0},
        {"topo_2_1_0_1_0_1_0_0", "2p 1n 1#pi^{-} 1K^{+}", 2, 1, 0, 1, 1, 0},
        {"topo_0_1_1_2_0_1_0_0", "1n 1#pi^{+} 2#pi^{-} 1K^{+}", 0, 1, 1, 2, 1, 0},
        {"topo_1_1_1_1_0_1_0_0", "1p 1n 1#pi^{+} 1#pi^{-} 1K^{+}", 1, 1, 1, 1, 1, 0},
        {"topo_0_1_1_0_0_1_0_0", "1n 1#pi^{+} 1K^{+}", 0, 1, 1, 0, 1, 0},
        {"topo_1_1_1_2_0_1_0_0", "1p 1n 1#pi^{+} 2#pi^{-} 1K^{+}", 1, 1, 1, 2, 1, 0},
        {"topo_2_1_1_1_0_1_0_0", "2p 1n 1#pi^{+} 1#pi^{-} 1K^{+}", 2, 1, 1, 1, 1, 0},
        {"topo_3_1_0_1_0_1_0_0", "3p 1n 1#pi^{-} 1K^{+}", 3, 1, 0, 1, 1, 0},
        {"topo_3_1_0_0_0_1_0_0", "3p 1n 1K^{+}", 3, 1, 0, 0, 1, 0},
        {"topo_2_1_1_2_0_1_0_0", "2p 1n 1#pi^{+} 2#pi^{-} 1K^{+}", 2, 1, 1, 2, 1, 0},
        {"topo_0_1_2_1_0_1_0_0", "1n 2#pi^{+} 1#pi^{-} 1K^{+}", 0, 1, 2, 1, 1, 0},
    };

    // 1n + 2K //
     vector<Topology> topos5 = {
        {"topo_0_1_0_0_0_1_1_0", "1n 1K^{+} 1K^{-}", 0, 1, 0, 0, 1, 1},
        {"topo_0_1_1_0_0_1_1_0", "1n 1#pi^{+} 1K^{+} 1K^{-}", 0, 1, 1, 0, 1, 1},
        {"topo_1_1_0_0_0_1_1_0", "1p 1n 1K^{+} 1K^{-}", 1, 1, 0, 0, 1, 1},
        {"topo_2_1_0_0_0_1_1_0", "2p 1n 1K^{+} 1K^{-}", 2, 1, 0, 0, 1, 1},
        {"topo_0_1_1_1_0_1_1_0", "1n 1#pi^{+} 1#pi^{-} 1K^{+} 1K^{-}", 0, 1, 1, 1, 1, 1},
        {"topo_1_1_1_0_0_1_1_0", "1p 1n 1#pi^{+} 1K^{+} 1K^{-}", 1, 1, 1, 0, 1, 1},
        {"topo_1_1_1_1_0_1_1_0", "1p 1n 1#pi^{+} 1#pi^{-} 1K^{+} 1K^{-}", 1, 1, 1, 1, 1, 1},
        {"topo_3_1_0_0_0_1_1_0", "3p 1n 1K^{+} 1K^{-}", 3, 1, 0, 0, 1, 1},
        {"topo_1_1_0_1_0_1_1_0", "1p 1n 1#pi^{-} 1K^{+} 1K^{-}", 1, 1, 0, 1, 1, 1},
        {"topo_2_1_0_1_0_1_1_0", "2p 1n 1#pi^{-} 1K^{+} 1K^{-}", 2, 1, 0, 1, 1, 1},
        {"topo_0_1_2_1_0_1_1_0", "1n 2#pi^{+} 1#pi^{-} 1K^{+} 1K^{-}", 0, 1, 2, 1, 1, 1},
        {"topo_2_1_1_1_0_1_1_0", "2p 1n 1#pi^{+} 1#pi^{-} 1K^{+} 1K^{-}", 2, 1, 1, 1, 1, 1},
    };

    // 1p + 1pi //
    vector<Topology> topos6 = {
        {"1p0n0pip1pim", "1p 0n 0#pi^{+} 1#pi^{-}", 1, 0, 0, 1, 0, 0},
        {"1p1n0pip1pim", "1p 1n 0#pi^{+} 1#pi^{-}", 1, 1, 0, 1, 0, 0},
        {"1p2n0pip1pim", "1p 2n 0#pi^{+} 1#pi^{-}", 1, 2, 0, 1, 0, 0},
        {"1p0n0pip1pim1kp", "1p 0n 0#pi^{+} 1#pi^{-} 1K^{+}", 1, 0, 0, 1, 1, 0},
        {"1p0n0pip1pim", "1p 0n 0#pi^{+} 1#pi^{-} ", 1, 0, 0, 1, 0, 0},
        {"1p3n0pip1pim", "1p 3n 0#pi^{+} 1#pi^{-}", 1, 3, 0, 1, 0, 0},
        {"1p1n0pip1pim1kp", "1p 1n 0#pi^{+} 1#pi^{-} 1K^{+}", 1, 1, 0, 1, 1, 0},
        {"1p4n0pip1pim", "1p 4n 0#pi^{+} 1#pi^{-}", 1, 4, 0, 1, 0, 0},
        {"1p0n0pip1pim1kp1km", "1p 0n 0#pi^{+} 1#pi^{-} 1K^{+} 1K^{-}", 1, 0, 0, 1, 1, 1},
        {"1p5n0pip1pim", "1p 5n 0#pi^{+} 1#pi^{-}", 1, 5, 0, 1, 0, 0},
        {"1p1n0pip1pim", "1p 1n 0#pi^{+} 1#pi^{-} ", 1, 1, 0, 1, 0, 0},
        {"1p2n0pip1pim1kp", "1p 2n 0#pi^{+} 1#pi^{-} 1K^{+}", 1, 2, 0, 1, 1, 0},
        {"1p1n0pip1pim1kp1km", "1p 1n 0#pi^{+} 1#pi^{-} 1K^{+} 1K^{-}", 1, 1, 0, 1, 1, 1},
        {"1p1n1pip", "1p 1n 1#pi^{+}", 1, 1, 1, 0, 0, 0},
        {"1p2n1pip", "1p 2n 1#pi^{+}", 1, 2, 1, 0, 0, 0},
        {"1p3n1pip", "1p 3n 1#pi^{+}", 1, 3, 1, 0, 0, 0},
        {"1p4n1pip", "1p 4n 1#pi^{+}", 1, 4, 1, 0, 0, 0},
        {"1p5n1pip", "1p 5n 1#pi^{+}", 1, 5, 1, 0, 0, 0},
        {"1p6n1pip", "1p 6n 1#pi^{+}", 1, 6, 1, 0, 0, 0},
        {"1p1n1pip", "1p 1n 1#pi^{+} ", 1, 1, 1, 0, 0, 0},
        {"1p0n1pip", "1p 0n 1#pi^{+} ", 1, 0, 1, 0, 0, 0},
        {"1p7n1pip", "1p 7n 1#pi^{+}", 1, 7, 1, 0, 0, 0},
        {"1p2n1pip", "1p 2n 1#pi^{+} ", 1, 2, 1, 0, 0, 0},
        {"1p0n1pip1km", "1p 0n 1#pi^{+} 1K^{-} ", 1, 0, 1, 0, 0, 1},
        {"1p1n1pip1kp1km", "1p 1n 1#pi^{+} 1K^{+} 1K^{-}", 1, 1, 1, 0, 1, 1},
        {"1p2n1pip1kp1km", "1p 2n 1#pi^{+} 1K^{+} 1K^{-}", 1, 2, 1, 0, 1, 1},
        {"1p3n1pip", "1p 3n 1#pi^{+} ", 1, 3, 1, 0, 0, 0},
        {"1p2n1pip1kp", "1p 2n 1#pi^{+} 1K^{+}", 1, 2, 1, 0, 1, 0}
    };

    // All Topo //
    vector<Topology> topos7 = {
        {"1p0n0pi",         "1p 0#pi",                    1,  0,   0,   0,   0,  0},
        {"1p0n1pim",        "1p 1#pi^{-}",                1,  0,   0,   1,   0,  0},
        {"1p0n1pip1pim",    "1p 1#pi^{+} 1#pi^{-}",       1,  0,   1,   1,   0,  0},
        {"1p0n1pip2pim",    "1p 1#pi^{+} 2#pi^{-}",       1,  0,   1,   2,   0,  0},
        {"1p0n2pip2pim",    "1p 2#pi^{+} 2#pi^{-}",       1,  0,   2,   2,   0,  0},
        {"1p1n0pi",         "1p 1n 0#pi",                 1,  1,   0,   0,   0,  0},
        {"1p1n1pim",        "1p 1n 1#pi^{-}",             1,  1,   0,   1,   0,  0},
        {"1p1n1pip1pim",    "1p 1n 1#pi^{+} 1#pi^{-}",    1,  1,   1,   1,   0,  0},
        {"1p2n0pi",         "1p 2n 0#pi",                 1,  2,   0,   0,   0,  0},
        {"1p3n0pi",         "1p 3n 0#pi",                 1,  3,   0,   0,   0,  0},
        {"2p0n0pi",         "2p 0#pi",                    2,  0,   0,   0,   0,  0},
        {"2p0n1pim",        "2p 1#pi^{-}",                2,  0,   0,   1,   0,  0},
        {"2p1n0pi",         "2p 1n 0#pi",                 2,  1,   0,   0,   0,  0},
        {"2p1n1pim",        "2p 1n 1#pi^{-}",             2,  1,   0,   1,   0,  0},
        {"2p2n0pi",         "2p 2n 0#pi",                 2,  2,   0,   0,   0,  0},
        {"3p0n0pi",         "3p 0#pi",                    3,  0,   0,   0,   0,  0},
        {"3p1n0pi",         "3p 1n 0#pi",                 3,  1,   0,   0,   0,  0},
        {"3p2n0pi",         "3p 2n 0#pi",                 3,  2,   0,   0,   0,  0},
        {"0p1n0pi",         "0p 1n 0#pi",                 0,  1,   0,   0,   0,  0},
        {"0p1n1pip",        "0p 1n 1#pi^{+}",             0,  1,   1,   0,   0,  0},
        {"0p1n1pip1pim",    "0p 1n 1#pi^{+} 1#pi^{-}",    0,  1,   1,   1,   0,  0},
        {"0p2n0pi",         "0p 2n 0#pi",                 0,  2,   0,   0,   0,  0},
        {"0p2n1pip",        "0p 2n 1#pi^{+}",             0,  2,   1,   0,   0,  0},
    };

    // 1p + 2 pi //
    vector<Topology> topos8 = {
        {"1p0n1pip1pim",           "1p 0n 1#pi^{+} 1#pi^{-}",           1,  0,   1,   1,   0,  0},
        {"1p1n1pip1pim",           "1p 1n 1#pi^{+} 1#pi^{-}",           1,  1,   1,   1,   0,  0},
        {"1p2n1pip1pim",           "1p 2n 1#pi^{+} 1#pi^{-}",           1,  2,   1,   1,   0,  0},
        {"1p3n1pip1pim",           "1p 3n 1#pi^{+} 1#pi^{-}",           1,  3,   1,   1,   0,  0},
        {"1p4n1pip1pim",           "1p 4n 1#pi^{+} 1#pi^{-}",           1,  4,   1,   1,   0,  0},
        {"1p0n1pip1pim",           "1p 0n 1#pi^{+} 1#pi^{-} ",          1,  0,   1,   1,   0,  0},
        {"1p5n1pip1pim",           "1p 5n 1#pi^{+} 1#pi^{-}",           1,  5,   1,   1,   0,  0},
        {"1p1n1pip1pim",           "1p 1n 1#pi^{+} 1#pi^{-} ",          1,  1,   1,   1,   0,  0},
        {"1p1n1pip1pim1kp",        "1p 1n 1#pi^{+} 1#pi^{-} 1K^{+}",    1,  1,   1,   1,   1,  0},
        {"1p6n1pip1pim",           "1p 6n 1#pi^{+} 1#pi^{-}",           1,  6,   1,   1,   0,  0},
        {"1p0n1pip1pim1kp",        "1p 0n 1#pi^{+} 1#pi^{-} 1K^{+}",    1,  0,   1,   1,   1,  0},
        {"1p0n1pip1pim1kp1km",     "1p 0n 1#pi^{+} 1#pi^{-} 1K^{+} 1K^{-}",  1,  0,   1,   1,   1,  1},
        {"1p2n1pip1pim1kp",        "1p 2n 1#pi^{+} 1#pi^{-} 1K^{+}",    1,  2,   1,   1,   1,  0},
        {"1p1n1pip1pim1kp1km",     "1p 1n 1#pi^{+} 1#pi^{-} 1K^{+} 1K^{-}",  1,  1,   1,   1,   1,  1},
        {"1p2n1pip1pim",           "1p 2n 1#pi^{+} 1#pi^{-} ",          1,  2,   1,   1,   0,  0},
        {"1p3n1pip1pim1kp",        "1p 3n 1#pi^{+} 1#pi^{-} 1K^{+}",    1,  3,   1,   1,   1,  0},
        {"1p0n0pip2pim1kp",        "1p 0n 0#pi^{+} 2#pi^{-} 1K^{+}",    1,  0,   0,   2,   1,  0},
        {"1p2n2pip",               "1p 2n 2#pi^{+}",                    1,  2,   2,   0,   0,  0},
        {"1p3n2pip",               "1p 3n 2#pi^{+}",                    1,  3,   2,   0,   0,  0},
        {"1p4n2pip",               "1p 4n 2#pi^{+}",                    1,  4,   2,   0,   0,  0},
        {"1p5n2pip",               "1p 5n 2#pi^{+}",                    1,  5,   2,   0,   0,  0},
        {"1p6n2pip",               "1p 6n 2#pi^{+}",                    1,  6,   2,   0,   0,  0},
        {"1p7n2pip",               "1p 7n 2#pi^{+}",                    1,  7,   2,   0,   0,  0},
    };

    // 1km //
    vector<Topology> topos9 = {
        
    };

    auto histos1 = MakeHistoSet(topos1);
    auto histos2 = MakeHistoSet(topos2);
    auto histos3 = MakeHistoSet(topos3);
    auto histos4 = MakeHistoSet(topos4);
    auto histos5 = MakeHistoSet(topos5);
    auto histos6 = MakeHistoSet(topos6);
    auto histos7 = MakeHistoSet(topos7);
    auto histos7_qel = MakeHistoSet(topos7);
    auto histos8 = MakeHistoSet(topos8);
    auto histos9 = MakeHistoSet(topos9);

    Long64_t nEntries = Root_Tree3->GetEntries();
    map<tuple<int,int,int,int,int,int,int,int>, int> topoCount;

    cout<<"Total number of events: " << nEntries << endl;

    for (Long64_t i = 0; i < nEntries; ++i) {
        if (i+1 % 1000000 == 0)
            cout <<"Processing event " << i <<" /" << nEntries << endl;

        Root_Tree3->GetEntry(i);
        if (cthl < -1 || cthl > 1) continue;

        int nProton = 0;
        int nNeutron = 0;
        int nPip = 0;
        int nPim = 0;
        int nPi0 = 0;
        int nKp = 0;
        int nKm = 0;
        int nK0 = 0;

        for (int j = 0; j < nf; ++j) {
            const bool aboveBaryon = pf[j] > kBaryonMomentumThreshold;
            const bool aboveMeson = pf[j] > kMesonMomentumThreshold;

            if      (pdgf[j] == 2212 && aboveBaryon) nProton++;
            else if (pdgf[j] == 2112 && aboveBaryon) nNeutron++;
            else if (pdgf[j] == 211  && aboveMeson)  nPip++;
            else if (pdgf[j] == -211 && aboveMeson)  nPim++;
            else if (pdgf[j] == 321  && aboveMeson)  nKp++;
            else if (pdgf[j] == -321 && aboveMeson)  nKm++;
        }

        double et = gst_Ev - gst_El;
        double theta = TMath::ACos(cthl);
        double Q2 = 2 * gst_Ev * gst_El * (1 - cthl);
        double W = TMath::Sqrt(TMath::Max(0.0, 0.938272 * 0.938272 + 2 * 0.938272 * et - Q2));

        topoCount[{nProton, nNeutron, nPip, nPim, nPi0, nKp, nKm, nK0}]++;

        for (size_t t = 0; t < topos1.size(); ++t) {
            if (MatchTopology(topos1[t], nProton, nNeutron, nPip, nPim, nPi0, nKp, nKm, nK0))
                FillHists(histos1[t], et, theta, Q2, W);
        }

        for (size_t t = 0; t < topos2.size(); ++t) {
            if (MatchTopology(topos2[t], nProton, nNeutron, nPip, nPim, nPi0, nKp, nKm, nK0))
                FillHists(histos2[t], et, theta, Q2, W);
        }

        for (size_t t = 0; t < topos3.size(); ++t) {
            if (MatchTopology(topos3[t], nProton, nNeutron, nPip, nPim, nPi0, nKp, nKm, nK0))
                FillHists(histos3[t], et, theta, Q2, W);
        }
        for (size_t t = 0; t < topos4.size(); ++t) {
            if (MatchTopology(topos4[t], nProton, nNeutron, nPip, nPim, nPi0, nKp, nKm, nK0))
                FillHists(histos4[t], et, theta, Q2, W);
        }

        for (size_t t = 0; t < topos5.size(); ++t) {
            if (MatchTopology(topos5[t], nProton, nNeutron, nPip, nPim, nPi0, nKp, nKm, nK0))
                FillHists(histos5[t], et, theta, Q2, W);
        }

        for (size_t t = 0; t < topos6.size(); ++t) {
            if (MatchTopology(topos6[t], nProton, nNeutron, nPip, nPim, nPi0, nKp, nKm, nK0))
                FillHists(histos6[t], et, theta, Q2, W);
        }

        for (size_t t = 0; t < topos7.size(); ++t) {
            if (MatchTopology(topos7[t], nProton, nNeutron, nPip, nPim, nPi0, nKp, nKm, nK0))
                FillHists(histos7[t], et, theta, Q2, W);
        }

        // Process topos7 only for quasielastic events
        if (qel) {
            for (size_t t = 0; t < topos7.size(); ++t) {
                if (MatchTopology(topos7[t], nProton, nNeutron, nPip, nPim, nPi0, nKp, nKm, nK0))
                    FillHists(histos7_qel[t], et, theta, Q2, W);
            }
        }

        for (size_t t = 0; t < topos8.size(); ++t) {
            if (MatchTopology(topos8[t], nProton, nNeutron, nPip, nPim, nPi0, nKp, nKm, nK0))
                FillHists(histos8[t], et, theta, Q2, W);
        }

        for (size_t t = 0; t < topos9.size(); ++t) {
            if (MatchTopology(topos9[t], nProton, nNeutron, nPip, nPim, nPi0, nKp, nKm, nK0))
                FillHists(histos9[t], et, theta, Q2, W);
        }
    }
    

    DrawTopologies(histos1, topos1,"K+","K^{+}", 10, false, "6GeV");
    DrawTopologies(histos1, topos1,"K+","K^{+}", 10, true, "6GeV");
    DrawTopologies(histos2, topos2,"1nP1nK","1nP 1nK", 10, false, "6GeV");
    DrawTopologies(histos2, topos2,"1nP1nK","1nP 1nK", 10, true, "6GeV");
    DrawTopologies(histos3, topos3,"1nP2nK","1nP 2nK", 10, false, "6GeV");
    DrawTopologies(histos3, topos3,"1nP2nK","1nP 2nK", 10, true, "6GeV");
    DrawTopologies(histos4, topos4,"1nN1nK","1nN 1nK", 10, false, "6GeV");
    DrawTopologies(histos4, topos4,"1nN1nK","1nN 1nK", 10, true, "6GeV");
    DrawTopologies(histos5, topos5,"1nN2nK","1nN 2nK", 10, false, "6GeV");
    DrawTopologies(histos5, topos5,"1nN2nK","1nN 2nK", 10, true, "6GeV");
    DrawTopologies(histos6, topos6,"1nP1npi","1nP 1n#pi", 10, false, "6GeV");
    DrawTopologies(histos6, topos6,"1nP1npi","1nP 1n#pi", 10, true, "6GeV");
    DrawTopologies(histos7, topos7,"AllTopologies","All Topologies", 10, false, "6GeV");
    DrawTopologies(histos7, topos7,"AllTopologies","All Topologies", 10, true, "6GeV");
    DrawTopologies(histos7_qel, topos7,"AllTopologies_QEL","All Topologies (QEL only)", 10, false, "6GeV");
    DrawTopologies(histos7_qel, topos7,"AllTopologies_QEL","All Topologies (QEL only)", 10, true, "6GeV");
    DrawTopologies(histos8, topos8,"1nP2npi","1nP 2n#pi", 10, false, "6GeV");
    DrawTopologies(histos8, topos8,"1nP2npi","1nP 2n#pi", 10, true, "6GeV");
    DrawTopologies(histos9, topos9,"1Km","1K^{-}", 10, false, "6GeV");
    DrawTopologies(histos9, topos9,"1Km","1K^{-}", 10, true, "6GeV");
}
