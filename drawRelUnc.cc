#include <fstream>
#include <sstream>
#include <iostream>
#include <map>
#include <vector>
#include <cmath>
#include "include/helperFunctions.h"


#include "TFile.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1.h"
#include "TStyle.h"
#include "TLine.h"





int main(int argc, char* argv[]){


    SetAtlasStyle();
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    gStyle->SetTitleSize(0.04, "XY");
    gStyle->SetLabelSize(0.04, "XY");

    gStyle->SetEndErrorSize(0); 

    std::string errFileName(argv[1]);
    std::ifstream infile(errFileName);
    std::string line;

    std::vector<SystEntry> entries;
    while (std::getline(infile, line)) {

        if (line.empty() || line[0] == '#') continue;

        std::stringstream ss(line);
        SystEntry e;

        ss >> e.group >> e.up >> e.nominal >> e.down >> e.groupID >> e.type >> e.colorHex;
        entries.push_back(e);
    }

    TH1D* h_pythia = loadHist("/data/ldelagra/LundPlane/Unfolding/syst_PrimaryLJP_2026/", "nominal", "UData");


    std::map<std::string, TH1D*> groupAccum;
    std::map<std::string, int> groupColor;
    for (auto& e : entries) {

        // std::cout << "au moins je suis ici .. "<< std::endl;
        TH1D* h_nom  = loadHist("/data/ldelagra/LundPlane/Unfolding/syst_PrimaryLJP_2026/", e.nominal, "UData");
        TH1D* h_up   = loadHist("/data/ldelagra/LundPlane/Unfolding/syst_PrimaryLJP_2026/", e.up, "UData");
        TH1D* h_down = loadHist("/data/ldelagra/LundPlane/Unfolding/syst_PrimaryLJP_2026/", e.down, "UData");
        // get the relative uncertainty 
        TH1D* h_rel = computeRel(h_nom, h_up, h_down, e.type);

        
        if (groupAccum.find(e.group) == groupAccum.end()) {
            groupAccum[e.group] = (TH1D*)h_rel->Clone(("h_" + e.group).c_str());
            groupAccum[e.group]->Reset();
        }

        for (int i = 1; i <= h_rel->GetNbinsX(); ++i) {
            double val = groupAccum[e.group]->GetBinContent(i);
            double add = h_rel->GetBinContent(i);

            groupAccum[e.group]->SetBinContent(i, val + add * add);
        }

        // add some colors to histos per group
        int col = TColor::GetColor(e.colorHex.c_str());
        groupColor[e.group] = col;

    }

    // Get group syst
    for (auto& [name, hist] : groupAccum) {
        for (int i = 1; i <= hist->GetNbinsX(); ++i) {
            hist->SetBinContent(i, std::sqrt(hist->GetBinContent(i)));
            hist->SetLineColor(groupColor[name]);
            hist->SetMarkerColor(groupColor[name]);
            hist->SetFillColorAlpha(groupColor[name], 0.35);

        }
    }


    // Get total syst
    TH1D* h_total = (TH1D*)h_pythia->Clone("h_total");
    h_total->Reset();
    for (int i = 1; i <= h_total->GetNbinsX(); ++i) {

        double sum2 = 0;
        for (auto& [name, hist] : groupAccum) {
            double v = hist->GetBinContent(i);
            sum2 += v * v;
        }
        h_total->SetBinContent(i, std::sqrt(sum2));
    }


    DrawTH1DPlots(h_total, groupAccum["TRACKING"], groupAccum["PURW"], groupAccum["JES"], groupAccum["MODELING"], -5, 5);

    
}
