
#include "include/Utils.h"
#include "include/PlottingTools.h"
#include "TMath.h"
#include "TStyle.h"
// #include "include/helperFunctionsFinal.h"

std::map<TString, TH2D*>  plot_LJP_2D(TString MCSample, TH1D* h1 ) {
    
    

    // TString inputPath = "/data/ldelagra/LundPlane/Unfolding//nominal_PrimaryLJP/Final_Files_16-03-2026//4iter/StatUnc.root";
    // TString inputPath = "/data/ldelagra/LundPlane/Unfolding/syst_PrimaryLJP_2026/nominal//4iter/Udata.root"; 
    TString outputDirStat = "/home/jmsardain/LJP/myplotting/2DPlots/"+MCSample+"/Stat/";
    TString outputDirLJP = "/home/jmsardain/LJP/myplotting/2DPlots/"+MCSample+"/";
    std::cout << " " << outputDirLJP << std::endl;
    std::cout << " " << outputDirStat << std::endl;
    gSystem->mkdir(outputDirStat, true);
    gSystem->mkdir(outputDirLJP, true);

    // Define binning
    std::vector<double> bins_kT = {-10, -3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 10};
    std::vector<double> bins_dR;
    for (double r = 0; r <= 5.0; r += 0.5) bins_dR.push_back(r);
    std::vector<double> bins_pT = {75, 118, 183, 277, 415, 657, 1034, 2011, 4256};

    const int nBins_kT = bins_kT.size() - 1;
    const int nBins_dR = bins_dR.size() - 1;
    const int nBins_pT = bins_pT.size() - 1;
    
    std::map<TString, TH2D*> hLJP_map;
    TCanvas*cLJP = new TCanvas("", "", 1200, 1200); 
    TCanvas*cStat = new TCanvas("", "", 1200, 1200); 
    TString outputFileLJP = Form("%sLundSlice%s_pTbins", outputDirLJP.Data(), MCSample.Data());
    TString outputFileStat = Form("%sLundSlice%s_pTbins", outputDirStat.Data(), MCSample.Data());

    cLJP->SaveAs(outputFileLJP+".pdf[");
    cStat->SaveAs(outputFileStat+".pdf[");
    // Loop over pT bins
    for (int ipt = 0; ipt < nBins_pT; ++ipt) {
        TString nameStat = Form("Stat_kTvsDeltaR_pTbin%d", ipt);
        TString nameLJP = Form("Lund_kTvsDeltaR_pTbin%d", ipt);

        TString title = Form("%.0f < p_{T} < %.0f GeV", bins_pT[ipt], bins_pT[ipt + 1]);
        TH2D* hStat = new TH2D(nameStat, title,
                            nBins_dR, bins_dR.data(),
                            nBins_kT, bins_kT.data());

        TH2D* hLJP = new TH2D(nameLJP, title,
                            nBins_dR, bins_dR.data(),
                            nBins_kT, bins_kT.data());

        for (int idr = 0; idr < nBins_dR; ++idr) {
            for (int ikt = 0; ikt < nBins_kT; ++ikt) {
                int globalBin = ikt + nBins_kT * (idr + nBins_dR * ipt) + 1; // ROOT bins start at 1
                double content = h1->GetBinContent(globalBin);
                double error = h1->GetBinError(globalBin);

                double kt_center = 0.5 * (bins_kT[ikt] + bins_kT[ikt+1]);
                double dr_center = 0.5 * (bins_dR[idr] + bins_dR[idr+1]);
                if (content > 0.) { 
                    hStat->Fill(dr_center, kt_center, error/content);
                    hLJP->Fill(dr_center, kt_center, content);
                }
                else { 
                    if (ipt==1) {
                        if ( bins_kT[ikt] >0 && bins_kT[ikt] < 4) {
                            if ( bins_dR[idr] ==0.5 ) {
                                std::cout << " global bin : " << globalBin << std::endl; 
                            }
                        }
                    }

                    hStat->Fill(dr_center, kt_center, 0.); 
                    hLJP->Fill(dr_center, kt_center, 0.);
                }
            }
        }



        Draw2DLund(cStat, hStat, outputFileStat, ipt, title, "Statistical Uncertainties");
        Draw2DLund(cLJP, hLJP, outputFileLJP, ipt, title, "Primary Lund Jet Plane");

        TString key = Form("%s_pTbin%d", MCSample.Data(), ipt);
        hLJP_map[key] = hLJP;
        delete hStat;
        // delete hLJP;
    }
    cLJP->SaveAs(outputFileLJP+".pdf]");
    cStat->SaveAs(outputFileStat+".pdf]");
    return hLJP_map;

}

std::map<TString, TH2D*> make_ratio_map( const std::map<TString, TH2D*>& numerator,
                                        const std::map<TString, TH2D*>& denominator,
                                        TString denominatorSample,
                                        TString ratioName)
{
    std::map<TString, TH2D*> ratioMap;

    
    for (const auto& [key, hNum] : numerator) {
        TString pTbin = key;
        pTbin.ReplaceAll("Pythia_", "");
        // the key for the alt. sample is altSample_pTbin not Pythia_pTbin .... if i don't do this, break seg violation
        TString denKey = denominatorSample + "_" + pTbin;

        auto it = denominator.find(denKey);
        TH2D* hDen = it->second;
        TH2D* hRatio = (TH2D*)hNum->Clone(Form("%s_%s", ratioName.Data(), key.Data()));

        hRatio->SetDirectory(0);
        hRatio->SetTitle( Form("%s, %s", ratioName.Data(), key.Data()));
        hRatio->Divide(hDen);
        ratioMap[key] = hRatio;
    }

    return ratioMap;
}

void save_ratio_maps( const std::map<TString, TH2D*>& ratioMap,
                        TString outputDir,
                        TString ratioTitle, 
                        TString sampleName) 
{
    gSystem->mkdir(outputDir, true);

    TCanvas*cRatio = new TCanvas("", "", 1200, 1200); 
    // TCanvas*cStat = new TCanvas("", "", 1200, 1200); 
    TString outputFile = outputDir + "Ratio_Pythia_"+sampleName+"_pTbin";

    cRatio->SaveAs(outputFile+".pdf[");

    for (const auto& [key, hist] : ratioMap) {
        TString pTbin = key;
        pTbin.ReplaceAll("Pythia_pTbin", "");
        int ipt = pTbin.Atoi();
        // TString outputFile = outputDir + Form("Ratio_pTbin%d", ipt);

        Draw2DLund(cRatio,  hist, outputFile, ipt, hist->GetTitle(), ratioTitle);
    }
    cRatio->SaveAs(outputFile+".pdf]");
}

int main(int argc, char* argv[]){
    
    SetAtlasStyle();
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    gStyle->SetTitleSize(0.04, "XY");
    gStyle->SetLabelSize(0.04, "XY");

    gStyle->SetEndErrorSize(0); 

    // std::string mcSample(argv[1]);
    // TString MCSample = TString(mcSample); 


    TH1D* h_pythia = loadHist("/data/ldelagra/LundPlane/Unfolding/syst_PrimaryLJP_2026/", "nominal", "UData");
    TH1D* h_powheg = loadHist("/data/ldelagra/LundPlane/Unfolding/syst_PrimaryLJP_2026/", "POWPY", "UData");
    TH1D* h_sherpa = loadHist("/data/ldelagra/LundPlane/Unfolding/syst_PrimaryLJP_2026/", "SHERPA", "UData");
    TH1D* h_herwig = loadHist("/data/ldelagra/LundPlane/Unfolding/syst_PrimaryLJP_2026/", "HERWIG", "UData");
    TH1D* h_data   = loadHist("/data/ldelagra/LundPlane/Unfolding/nominal_PrimaryLJP/", "Final_Files_16-03-2026", "UData_TotStat");
    TH1D* h_unfold = loadHist("/data/ldelagra/LundPlane/Unfolding/nominal_PrimaryLJP/Final_Files_16-03-2026/BiasIDS/", "", "biasIDS");
    
    
    std::map<TString, TH2D*> hLJP_pythia = plot_LJP_2D("Pythia", h_pythia);
    std::map<TString, TH2D*> hLJP_powheg = plot_LJP_2D("POWPY", h_powheg); 
    std::map<TString, TH2D*> hLJP_sherpa = plot_LJP_2D("SHERPA", h_sherpa);
    std::map<TString, TH2D*> hLJP_herwig = plot_LJP_2D("HERWIG", h_herwig);
    std::map<TString, TH2D*> hLJP_data   = plot_LJP_2D("Data", h_data);
    std::map<TString, TH2D*> hLJP_unfold = plot_LJP_2D("UNFOLD", h_unfold);

    std::map<TString, TH2D*> ratio_pythia_sherpa = make_ratio_map(hLJP_pythia, hLJP_sherpa, "SHERPA",  "Pythia_Sherpa");
    std::map<TString, TH2D*> ratio_pythia_powheg = make_ratio_map(hLJP_pythia, hLJP_powheg, "POWPY", "Pythia_Powheg");
    std::map<TString, TH2D*> ratio_pythia_herwig = make_ratio_map(hLJP_pythia, hLJP_herwig, "HERWIG", "Pythia_Herwig");
    
    save_ratio_maps( ratio_pythia_sherpa, "/home/jmsardain/LJP/myplotting/2DPlots/Ratios/Pythia_Sherpa/", "Pythia / Sherpa", "Sherpa");
    save_ratio_maps( ratio_pythia_powheg, "/home/jmsardain/LJP/myplotting/2DPlots/Ratios/Pythia_Powheg/", "Pythia / POWPY", "POWPY");
    save_ratio_maps( ratio_pythia_herwig, "/home/jmsardain/LJP/myplotting/2DPlots/Ratios/Pythia_Herwig/", "Pythia / HERWIG", "HERWIG");
    return 0;
}
