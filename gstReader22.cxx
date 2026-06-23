#include <TFile.h>
#include <TTree.h>
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
    int nPi0;
    int nKp;
    int nKm;
    int nK0;
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
                   int nPi0,
                   int nKp,
                   int nKm,
                   int nK0)
{
    if (tp.nProton >= 0 && nProton != tp.nProton) return false;
    if (tp.nNeutron >= 0 && nNeutron != tp.nNeutron) return false;
    if (tp.nPip >= 0 && nPip != tp.nPip) return false;
    if (tp.nPim >= 0 && nPim != tp.nPim) return false;
    if (tp.nPi0 >= 0 && nPi0 != tp.nPi0) return false;
    if (tp.nKp >= 0 && nKp != tp.nKp) return false;
    if (tp.nKm >= 0 && nKm != tp.nKm) return false;
    if (tp.nK0 >= 0 && nK0 != tp.nK0) return false;
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
                    bool logY = false)
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

    TString fileSuffix = prefix + (logY ?"_log" :"");
    c->SaveAs(fileSuffix +".png");
}

void gstReader19()
{
    const string File_Path ="";
    const string File_Name ="11_1000060120_EM_5_98636GeV_G18_02a_00_000_Q2_0_4.gst.root";
    const string File_TTree_Name ="gst";
    const string File_Path_and_Name = File_Path + File_Name;

    TFile* Root_File3 = TFile::Open(File_Path_and_Name.c_str());
    if (!Root_File3 || Root_File3->IsZombie()) {
        cerr <<"Error opening file:" << File_Path_and_Name << endl;
        return;
    }

    TTree* Root_Tree3 = static_cast<TTree*>(Root_File3->Get(File_TTree_Name.c_str()));
    if (!Root_Tree3) {
        cerr <<"Error getting tree:" << File_TTree_Name << endl;
        Root_File3->Close();
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

    vector<Topology> topos1 = {
       
        {"1p0n0pip1pim", "1p 0n 0#pi^{+} 1#pi^{-}", 1, 0, 0, 1, 0, 0, 0, 0},
        {"1p1n0pip1pim", "1p 1n 0#pi^{+} 1#pi^{-}", 1, 1, 0, 1, 0, 0, 0, 0},
        {"1p2n0pip1pim", "1p 2n 0#pi^{+} 1#pi^{-}", 1, 2, 0, 1, 0, 0, 0, 0},
        {"1p0n0pip1pim1kp", "1p 0n 0#pi^{+} 1#pi^{-} 1K^{+}", 1, 0, 0, 1, 0, 1, 0, 0},
        {"1p0n0pip1pim1k0", "1p 0n 0#pi^{+} 1#pi^{-} 1K^{0}", 1, 0, 0, 1, 0, 0, 0, 1},
        {"1p3n0pip1pim", "1p 3n 0#pi^{+} 1#pi^{-}", 1, 3, 0, 1, 0, 0, 0, 0},
        {"1p1n0pip1pim1kp", "1p 1n 0#pi^{+} 1#pi^{-} 1K^{+}", 1, 1, 0, 1, 0, 1, 0, 0},
        {"1p4n0pip1pim", "1p 4n 0#pi^{+} 1#pi^{-}", 1, 4, 0, 1, 0, 0, 0, 0},
        {"1p0n0pip1pim1kp1km", "1p 0n 0#pi^{+} 1#pi^{-} 1K^{+} 1K^{-}", 1, 0, 0, 1, 0, 1, 1, 0},
        {"1p5n0pip1pim", "1p 5n 0#pi^{+} 1#pi^{-}", 1, 5, 0, 1, 0, 0, 0, 0},
        {"1p1n0pip1pim1k0", "1p 1n 0#pi^{+} 1#pi^{-} 1K^{0}", 1, 1, 0, 1, 0, 0, 0, 1},
        {"1p2n0pip1pim1kp", "1p 2n 0#pi^{+} 1#pi^{-} 1K^{+}", 1, 2, 0, 1, 0, 1, 0, 0},
        {"1p1n0pip1pim1kp1km", "1p 1n 0#pi^{+} 1#pi^{-} 1K^{+} 1K^{-}", 1, 1, 0, 1, 0, 1, 1, 0},
        {"1p1n1pip", "1p 1n 1#pi^{+}", 1, 1, 1, 0, 0, 0, 0, 0},
        {"1p2n1pip", "1p 2n 1#pi^{+}", 1, 2, 1, 0, 0, 0, 0, 0},
        {"1p3n1pip", "1p 3n 1#pi^{+}", 1, 3, 1, 0, 0, 0, 0, 0},
        {"1p4n1pip", "1p 4n 1#pi^{+}", 1, 4, 1, 0, 0, 0, 0, 0},
        {"1p5n1pip", "1p 5n 1#pi^{+}", 1, 5, 1, 0, 0, 0, 0, 0},
        {"1p6n1pip", "1p 6n 1#pi^{+}", 1, 6, 1, 0, 0, 0, 0, 0},
        {"1p1n1pip1k0", "1p 1n 1#pi^{+} 1K^{0}", 1, 1, 1, 0, 0, 0, 0, 1},
        {"1p0n1pip1k0", "1p 0n 1#pi^{+} 1K^{0}", 1, 0, 1, 0, 0, 0, 0, 1},
        {"1p7n1pip", "1p 7n 1#pi^{+}", 1, 7, 1, 0, 0, 0, 0, 0},
        {"1p2n1pip1k0", "1p 2n 1#pi^{+} 1K^{0}", 1, 2, 1, 0, 0, 0, 0, 1},
        {"1p0n1pip1km1k0", "1p 0n 1#pi^{+} 1K^{-} 1K^{0}", 1, 0, 1, 0, 0, 0, 1, 1},
        {"1p1n1pip1kp1km", "1p 1n 1#pi^{+} 1K^{+} 1K^{-}", 1, 1, 1, 0, 0, 1, 1, 0},
        {"1p2n1pip1kp1km", "1p 2n 1#pi^{+} 1K^{+} 1K^{-}", 1, 2, 1, 0, 0, 1, 1, 0},
        {"1p3n1pip1k0", "1p 3n 1#pi^{+} 1K^{0}", 1, 3, 1, 0, 0, 0, 0, 1},
        {"1p2n1pip1kp", "1p 2n 1#pi^{+} 1K^{+}", 1, 2, 1, 0, 0, 1, 0, 0}
    };


    


    auto histos1 = MakeHistoSet(topos1);



    Long64_t nEntries = Root_Tree3->GetEntries();
    map<tuple<int,int,int,int,int,int,int,int>, int> topoCount;

    for (Long64_t i = 0; i < nEntries; ++i) {
        if (i % 1000000 == 0)
            cout <<"Processing event" << i <<" /" << nEntries << endl;

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
            else if (pdgf[j] == 111  && aboveMeson)  nPi0++;
            else if (pdgf[j] == 321  && aboveMeson)  nKp++;
            else if (pdgf[j] == -321 && aboveMeson)  nKm++;
            else if (pdgf[j] == 311  && aboveMeson)  nK0++;
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
    }

    DrawTopologies(histos1, topos1,"p+pi","p+pi (w/threshholds)", 10, false);
    DrawTopologies(histos1, topos1,"p+pi","p+pi (w/threshholds)", 10, true);


    Root_File3->Close();
}
