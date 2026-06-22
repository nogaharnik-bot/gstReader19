#include <TFile.h>
#include <TTree.h>
#include <TLegend.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TDirectory.h>
#include <TChain.h>
#include <THStack.h>
#include <map>

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <fstream>
#include <stdlib.h>

using namespace std;

void gstReader12()
{
    string File_Path          = "";
    string File_Name          = "11_1000060120_EM_5_98636GeV_G18_02a_00_000_Q2_0_4.gst.root";
    string File_TTree_Name    = "gst";
    string File_Path_and_Name = File_Path + File_Name;
    
    TFile* Root_File3 = TFile::Open(File_Path_and_Name.c_str());
    TTree* Root_Tree3 = (TTree*)Root_File3->Get(File_TTree_Name.c_str());

    Root_Tree3->SetBranchStatus("*",    0);
    Root_Tree3->SetBranchStatus("Ev",   1);
    Root_Tree3->SetBranchStatus("El",   1);
    Root_Tree3->SetBranchStatus("cthl", 1);
    Root_Tree3->SetBranchStatus("nf",   1);
    Root_Tree3->SetBranchStatus("pdgf", 1);
    Root_Tree3->SetBranchStatus("pf",   1);
    Root_Tree3->SetBranchStatus("qel", 1);

    int    nf = 0;
    const  int kMax = 200;
    int    pdgf[kMax];
    double pf[kMax];
    double gst_Ev=0.0, gst_El=0.0, cthl=0.0;
    bool   qel;

    Root_Tree3->SetBranchAddress("Ev",   &gst_Ev);
    Root_Tree3->SetBranchAddress("El",   &gst_El);
    Root_Tree3->SetBranchAddress("cthl", &cthl);
    Root_Tree3->SetBranchAddress("nf",   &nf);
    Root_Tree3->SetBranchAddress("pdgf",  pdgf);
    Root_Tree3->SetBranchAddress("pf",    pf);
    Root_Tree3->SetBranchAddress("qel", &qel);

    auto MakeHists = [](TString name, TString title) -> vector<TH1D*>
    {
        return {
            new TH1D("et_"+name,    title+";Energy Transfer (GeV);Events", 100, 0,   6),
            new TH1D("theta_"+name, title+";#theta_{e} (rad);Events",      100, 0,   .5),
            new TH1D("Q2_"+name,    title+";Q^{2} (GeV^{2});Events",       100, 0,   1),
            new TH1D("W_"+name,     title+";W (GeV);Events",               100, 0, 3.5)
        };
    };

    struct Topology {
        TString name, title;
        int nProton, nNeutron, nPip, nPim, nPi0, nKp, nKm, nK0, nPi, nK;
    };

    vector<Topology> topos = {
        // name           title            nP  nN   nPip nPim nPi0 nKp nKm nK0 nPi nK
    // nP>0 nKp=1, nKm=1
{"1p0n1pim_Kp1Km0",       "1p 0n 1#pi^{-} K^{+}",              1,  0,  0,  1,  0,  1,  0,  0,  0,  0},
{"1p1n0pi_Kp1Km0",        "1p 1n 0#pi K^{+}",                  1,  1,  0,  0,  0,  1,  0,  0,  0,  0},
{"1p0n0pi_Kp1Km0",        "1p 0n 0#pi K^{+}",                  1,  0,  0,  0,  0,  1,  0,  0,  0,  0},
{"1p1n1pim_Kp1Km0",       "1p 1n 1#pi^{-} K^{+}",              1,  1,  0,  1,  0,  1,  0,  0,  0,  0},
{"2p1n0pi_Kp1Km0",        "2p 1n 0#pi K^{+}",                  2,  1,  0,  0,  0,  1,  0,  0,  0,  0},
{"2p1n1pim_Kp1Km0",       "2p 1n 1#pi^{-} K^{+}",              2,  1,  0,  1,  0,  1,  0,  0,  0,  0},
{"1p1n1pip1pim_Kp1Km0",   "1p 1n 1#pi^{+}1#pi^{-} K^{+}",      1,  1,  1,  1,  0,  1,  0,  0,  0,  0},
{"1p0n1pip1pim_Kp1Km0",   "1p 0n 1#pi^{+}1#pi^{-} K^{+}",      1,  0,  1,  1,  0,  1,  0,  0,  0,  0},
{"2p0n1pim_Kp1Km0",       "2p 0n 1#pi^{-} K^{+}",              2,  0,  0,  1,  0,  1,  0,  0,  0,  0},
{"1p2n0pi_Kp1Km0",        "1p 2n 0#pi K^{+}",                  1,  2,  0,  0,  0,  1,  0,  0,  0,  0},
{"1p0n1pip2pim_Kp1Km0",   "1p 0n 1#pi^{+}2#pi^{-} K^{+}",      1,  0,  1,  2,  0,  1,  0,  0,  0,  0},
{"2p2n0pi_Kp1Km0",        "2p 2n 0#pi K^{+}",                  2,  2,  0,  0,  0,  1,  0,  0,  0,  0},
{"1p2n1pip1pim_Kp1Km0",   "1p 2n 1#pi^{+}1#pi^{-} K^{+}",      1,  2,  1,  1,  0,  1,  0,  0,  0,  0},
{"1p1n1pip2pim_Kp1Km0",   "1p 1n 1#pi^{+}2#pi^{-} K^{+}",      1,  1,  1,  2,  0,  1,  0,  0,  0,  0},
{"2p1n1pip1pim_Kp1Km0",   "2p 1n 1#pi^{+}1#pi^{-} K^{+}",      2,  1,  1,  1,  0,  1,  0,  0,  0,  0},
{"1p2n1pim_Kp1Km0",       "1p 2n 1#pi^{-} K^{+}",              1,  2,  0,  1,  0,  1,  0,  0,  0,  0},
{"1p0n2pim_Kp1Km0",       "1p 0n 2#pi^{-} K^{+}",              1,  0,  0,  2,  0,  1,  0,  0,  0,  0},
{"3p1n1pim_Kp1Km0",       "3p 1n 1#pi^{-} K^{+}",              3,  1,  0,  1,  0,  1,  0,  0,  0,  0},
{"2p2n1pim_Kp1Km0",       "2p 2n 1#pi^{-} K^{+}",              2,  2,  0,  1,  0,  1,  0,  0,  0,  0},
{"1p3n0pi_Kp1Km0",        "1p 3n 0#pi K^{+}",                  1,  3,  0,  0,  0,  1,  0,  0,  0,  0},
{"2p2n1pip1pim_Kp1Km0",   "2p 2n 1#pi^{+}1#pi^{-} K^{+}",      2,  2,  1,  1,  0,  1,  0,  0,  0,  0},
{"2p3n0pi_Kp1Km0",        "2p 3n 0#pi K^{+}",                  2,  3,  0,  0,  0,  1,  0,  0,  0,  0},
{"1p2n1pip_Kp1Km0",       "1p 2n 1#pi^{+} K^{+}",              1,  2,  1,  0,  0,  1,  0,  0,  0,  0},
{"3p2n0pi_Kp1Km0",        "3p 2n 0#pi K^{+}",                  3,  2,  0,  0,  0,  1,  0,  0,  0,  0},
{"2p0n2pim_Kp1Km0",       "2p 0n 2#pi^{-} K^{+}",              2,  0,  0,  2,  0,  1,  0,  0,  0,  0},
{"3p1n0pi_Kp1Km0",        "3p 1n 0#pi K^{+}",                  3,  1,  0,  0,  0,  1,  0,  0,  0,  0},
{"2p1n2pim_Kp1Km0",       "2p 1n 1#pi^{+}2#pi^{-} K^{+}",      2,  1,  1,  2,  0,  1,  0,  0,  0,  0},
{"3p2n1pim_Kp1Km0",       "3p 2n 1#pi^{-} K^{+}",              3,  2,  0,  1,  0,  1,  0,  0,  0,  0},
{"1p3n1pip1pim_Kp1Km0",   "1p 3n 1#pi^{+}1#pi^{-} K^{+}",      1,  3,  1,  1,  0,  1,  0,  0,  0,  0},
{"1p4n0pi_Kp1Km0",        "1p 4n 0#pi K^{+}",                  1,  4,  0,  0,  0,  1,  0,  0,  0,  0},
};    

    vector<vector<TH1D*>> histos;
    for (auto& t : topos)
        histos.push_back(MakeHists(t.name, t.title));

//    vector<TH2D*> h2ds;
//    for (auto& t : topos)
//    {
//        histos.push_back(MakeHists(t.name, t.title));
//        h2ds.push_back(new TH2D("h2d_"+t.name, t.title+";W (GeV);Q^{2} (GeV^{2})",
//                                 100, 0.8, 2, 100, 0, 1));
//    }

    auto FillHists = [](vector<TH1D*>& h, double et, double theta, double Q2, double W)
    {
        h[0]->Fill(et);
        h[1]->Fill(theta);
        h[2]->Fill(Q2);
        h[3]->Fill(W);
    };

    Long64_t nEntries = Root_Tree3->GetEntries();
    map<tuple<int,int,int,int,int,int,int,int>, int> topoCount;

    for (Long64_t i = 0; i < nEntries; i++)
    {
        if (i % 1000000 == 0)
            cout << "Processing event " << i << " / " << nEntries << endl;

        Root_Tree3->GetEntry(i);
        if (cthl < -1 || cthl > 1) continue;

        int nProton=0, nNeutron=0, nPip=0, nPim=0, nPi0=0, nKp=0, nKm=0, nK0=0;

        for (int j = 0; j < nf; j++)
        {
            if      (pdgf[j] == 2212) nProton++;
            else if (pdgf[j] == 2112) nNeutron++;
            else if (pdgf[j] == 211)  nPip++;
            else if (pdgf[j] == -211) nPim++;
            else if (pdgf[j] == 111)  nPi0++;
            else if (pdgf[j] == 321)  nKp++;
            else if (pdgf[j] == -321) nKm++;
            else if (pdgf[j] == 311)  nK0++;
        }

        int nPi = nPip + nPim + nPi0;
        int nK  = nKp  + nKm  + nK0;

        double et    = gst_Ev - gst_El;
        double theta = TMath::ACos(cthl);
        double Q2    = 2 * gst_Ev * gst_El * (1 - cthl);
        double W     = TMath::Sqrt(TMath::Max(0.0, 0.938272*0.938272 + 2*0.938272*et - Q2));

        topoCount[{nProton, nNeutron, nPip, nPim, nPi0, nKp, nKm, nK0}]++;

        for (int t = 0; t < (int)topos.size(); t++)
        {
            auto& tp = topos[t];
	bool protonMatch  = (tp.nProton  == -1) ? true : (nProton  == tp.nProton);
	bool neutronMatch = (tp.nNeutron == -1) ? true : (nNeutron == tp.nNeutron);
	bool pionMatch    = (tp.nPi > 0)  ? (nPi == tp.nPi)  :
                        (tp.nPi == -1) ? true :
                        (nPip == tp.nPip && nPim == tp.nPim && nPi0 == tp.nPi0);
	bool kaonMatch    = (tp.nK  > 0)  ? (nK  == tp.nK)   :
                        (tp.nK == -1)  ? true :
                        (nKp == tp.nKp && nKm == tp.nKm && nK0 == tp.nK0);
	if (nProton == tp.nProton &&
        	neutronMatch &&
        	pionMatch   &&
		kaonMatch)
            {
                FillHists(histos[t], et, theta, Q2, W);
//                h2ds[t]->Fill(W, Q2);
            }
        }
    }

    //  topology counts 
//    cout << "\nnP nN nPip nPim nPi0 nKp nKm nK0 : count : qel" << endl;
//    for (auto& entry : topoCount)
//    {
//        auto& k = entry.first;
//        cout << get<0>(k) << "  " << get<1>(k) << "  " << get<2>(k) << "  "
//             << get<3>(k) << "  " << get<4>(k) << "  " << get<5>(k) << "  "
//             << get<6>(k) << "  " << get<7>(k) << "  :  " << entry.second << endl;
//    }

//    cout << "\nTopology entries:" << endl;
//    for (int t = 0; t < (int)topos.size(); t++)
//        cout << topos[t].name << ": " << histos[t][0]->GetEntries() << endl;

    vector<int> colors = {
        kBlack, kRed, kBlue, kGreen+2, kOrange+7,
        kViolet+1, kCyan+2, kMagenta+1, kYellow+3, kRed+3,
        kBlue+3, kGreen+4, kOrange-3, kPink+9, kAzure+7,
        kTeal+3, kSpring-6, kViolet-4, kOrange+3, kRed-7,
        kBlue-7, kGreen-6, kMagenta-7, kCyan-6
    };

    vector<int> lineStyles = {1,1,1,1,1,1,1,1, 2,2,2,2,2, 3,3,3,3,3, 7,7,7,7,7, 9,9,9,9,9};

    int minEntries = 10;

//auto DrawTopo = [](vector<TH1D*>& h, TString name)
//    {
//        TCanvas* c = new TCanvas("c_"+name, name, 1600, 1200);
//        c->Divide(2, 2);
//        for (int pad = 1; pad <= 4; pad++)
//        {
//            c->cd(pad);
//            gPad->SetLeftMargin(0.15);
//            gPad->SetRightMargin(0.15);
//	    gPad->SetBottomMargin(0.15);
//            h[pad-1]->Draw("hist");
//        }
//        c->SaveAs(name + ".png");
//    };

//    for (int t = 0; t < (int)topos.size(); t++)
//        DrawTopo(histos[t], topos[t].name);

    TCanvas* cALL = new TCanvas("cALL", "nP>01nKp0nKm Topologies", 1600, 1200);
    cALL->Divide(2, 2);

    for (int var = 0; var < 4; var++)
    {
        cALL->cd(var+1);
        gPad->SetLeftMargin(0.15);
	gPad->SetRightMargin(0.15);
        gPad->SetBottomMargin(0.15);
	gStyle->SetOptStat(0);
        double maxVal = 0;
        for (int t = 0; t < (int)topos.size(); t++)
        {
            if (histos[t][var]->GetEntries() < minEntries) continue;
            TH1D* hc = (TH1D*)histos[t][var]->Clone();
            maxVal = max(maxVal, hc->GetMaximum());
            delete hc;
        }

        TLegend* leg = new TLegend(0.87, 0.1, 0.99, 0.9);
        leg->SetTextSize(0.022);
        leg->SetBorderSize(1);

        bool first = true;
        int colorIdx = 0;
        for (int t = 0; t < (int)topos.size(); t++)
        {
            if (histos[t][var]->GetEntries() < minEntries) continue;

            TH1D* hc = (TH1D*)histos[t][var]->Clone(
                TString::Format("clone_%d_%d", t, var));

            hc->SetLineColor(colors[colorIdx % colors.size()]);
            hc->SetLineStyle(lineStyles[colorIdx % lineStyles.size()]);
            hc->SetLineWidth(2);
            hc->SetFillStyle(0);
            hc->GetYaxis()->SetTitle("Events");
            hc->SetTitle("nP>0 1nKp 0nKm Topologies");
            hc->SetMaximum(maxVal * 1.1);
            hc->SetMinimum(0);

            if (first) { hc->Draw("hist"); first = false; }
            else         hc->Draw("hist same");
            leg->AddEntry(hc, topos[t].title, "l");
            colorIdx++;
        }
        leg->Draw();
    }
    cALL->SaveAs("nP>01nkp0nkm_topologies.png");

    TCanvas* cLOG = new TCanvas("cLOG", "Log nP>01nKp0nKm Topologies", 1600, 1200);
    cLOG->Divide(2, 2);

    for (int var = 0; var < 4; var++)
    {
        cLOG->cd(var+1);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.15);
        gPad->SetBottomMargin(0.15);
	gPad->SetLogy(1);
        gStyle->SetOptStat(0);
        double maxVal = 0;
        for (int t = 0; t < (int)topos.size(); t++)
        {
            if (histos[t][var]->GetEntries() < minEntries) continue;
            TH1D* hc = (TH1D*)histos[t][var]->Clone();
            maxVal = max(maxVal, hc->GetMaximum());
            delete hc;
        }

        TLegend* leg  = new TLegend(0.87, 0.1, 0.99, 0.9);
        leg->SetTextSize(0.022);
        leg->SetBorderSize(1);

        bool first = true;
        int colorIdx = 0;
        for (int t = 0; t < (int)topos.size(); t++)
        {
            if (histos[t][var]->GetEntries() < minEntries) continue;

            TH1D* hc = (TH1D*)histos[t][var]->Clone(
                TString::Format("clone_%d_%d", t, var));

            hc->SetLineColor(colors[colorIdx % colors.size()]);
            hc->SetLineStyle(lineStyles[colorIdx % lineStyles.size()]);
            hc->SetLineWidth(2);
            hc->SetFillStyle(0);
            hc->GetYaxis()->SetTitle("Events");
            hc->SetTitle("nP>0 1nKp 0nKm Topologies");
            hc->SetMaximum(maxVal * 1.1);
            hc->SetMinimum(1);

            if (first) { hc->Draw("hist"); first = false; }
            else         hc->Draw("hist same");
            leg->AddEntry(hc, topos[t].title, "l");
            colorIdx++;
        }
        leg->Draw();
	    }
	    cLOG->SaveAs("nP>01nkp0nkm_topologies->log.png");
    // 2D 
//    gStyle->SetPalette(kBird);
//    for (int t = 0; t < (int)topos.size(); t++)
//    {
//        if (h2ds[t]->GetEntries() < minEntries) continue;
//        TCanvas* c2d = new TCanvas("c2d_"+topos[t].name, topos[t].name, 800, 600);
//        c2d->SetLeftMargin(0.15);
//        c2d->SetBottomMargin(0.15);
//        h2ds[t]->Draw("colz");
///        c2d->SaveAs("h2d_"+topos[t].name+".png");
//    }

    Root_File3->Close();
}
