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
    TH1D* h_powheg = loadHist("/data/ldelagra/LundPlane/Unfolding/syst_PrimaryLJP_2026/", "POWPY", "UData");
    TH1D* h_sherpa = loadHist("/data/ldelagra/LundPlane/Unfolding/syst_PrimaryLJP_2026/", "SHERPA", "UData");
    TH1D* h_herwig = loadHist("/data/ldelagra/LundPlane/Unfolding/syst_PrimaryLJP_2026/", "HERWIG", "UData");
    TH1D* h_data   = loadHist("/data/ldelagra/LundPlane/Unfolding/nominal_PrimaryLJP/", "Final_Files_16-03-2026", "UData_TotStat");
    TH1D* h_unfold = loadHist("/data/ldelagra/LundPlane/Unfolding/nominal_PrimaryLJP/", "", "biasIDS");
    // applyStatUncertainty(h_data, "/data/ldelagra/LundPlane/Unfolding/nominal_PrimaryLJP/", "Final_Files_16-03-2026"); // fonction obselete

    auto kt_pythia = extract_kT_distributions(h_pythia, "nominal");
    auto kt_powheg = extract_kT_distributions(h_powheg, "POWPY");
    auto kt_sherpa = extract_kT_distributions(h_sherpa, "SHERPA");
    auto kt_herwig = extract_kT_distributions(h_herwig, "HERWIG");
    auto kt_data   = extract_kT_distributions(h_data, "DATA");

    // draw 1D histos (2 pads, top: true, middle: mc/data)
    DrawTH1DPlots(kt_data, kt_pythia, kt_powheg, kt_sherpa, kt_herwig);

    // draw 1D histos (3 pads, top: true, middle: mc/data, bottom: rel. unc.)
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
        sum2 += h_unfold->GetBinContent(i) * h_unfold->GetBinContent(i); // add the unfolding uncertainty to the total
        for (auto& [name, hist] : groupAccum) {
            double v = hist->GetBinContent(i);
            sum2 += v * v;
        }
        h_total->SetBinContent(i, std::sqrt(sum2));
    }

    // these are already rel uncertainties, not 1D distributions from the ROOT files
    auto kt_purw     = extract_kT_distributions(groupAccum["PURW"], "purw");
    auto kt_tracking = extract_kT_distributions(groupAccum["TRACKING"], "tracking");
    auto kt_JES      = extract_kT_distributions(groupAccum["JES"], "JES");
    auto kt_Modeling = extract_kT_distributions(groupAccum["MODELING"], "modeling");
    auto kt_total    = extract_kT_distributions(h_total, "total");
    auto kt_unfold   = extract_kT_distributions(h_unfold, "biasIDS");

    DrawTH1DPlotsWithError( kt_data, kt_pythia, kt_powheg, kt_sherpa, kt_herwig, kt_total, kt_tracking, kt_purw, kt_JES, kt_Modeling, kt_unfold);
    

}