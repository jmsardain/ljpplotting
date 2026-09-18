// g++ `root-config --cflags --ldflags --libs --glibs` plot_LJP_2D.cpp -o plot_LJP_2D.exe

#include "utils/Utils.h"
#include "utils/PlottingTools.h"
#include "TMath.h"

void plot_LJP_2D() {
    SetAtlasStyle();

    TString inputPath = "/data/ldelagra/LundPlane/Unfolding//nominal_PrimaryLJP/Final_Files_16-03-2026//4iter/StatUnc.root";
    TString histName = "UData_TotStat";
    TString outputDirStat = "/data/ldelagra/LundPlane/Unfolding//nominal_PrimaryLJP/Final_Files_16-03-2026//4iter/plots/Stat";
    TString outputDirLJP = "/data/ldelagra/LundPlane/Unfolding//nominal_PrimaryLJP/Final_Files_16-03-2026//4iter/plots//";

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

    TFile* file = TFile::Open(inputPath, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Could not open input file: " << inputPath << std::endl;
        return;
    }

    TH1D* h1 = (TH1D*)file->Get(histName);
    if (!h1) {
        std::cerr << "Histogram " << histName << " not found!" << std::endl;
        file->Close();
        return;
    }

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

                //std::cout << "globalBin " << globalBin << std::endl;
                //std::cout << "Filling bin (" << dr_center << ", "<< kt_center << ") with " << content << std::endl;
                
                //if ( error/content < 0.) {std::cout << "negative entry for bin (" << dr_center << ", "<< kt_center << ") with " << error/content << std::endl;}
                if (content > 0.) { 
                    //std::cout << "=====================================" << std::endl;
                    //std::cout << "global bin " << globalBin << std::endl;
                    //std::cout << "Filling bin (" << dr_center << ", "<< kt_center << ")" << std::endl;
                    //std::cout << "error " << error << std::endl;
                    //std::cout << "content " << content << std::endl;
                    //std::cout << "error/content " << error/content << std::endl;
                    hStat->Fill(dr_center, kt_center, error/content);
                    hLJP->Fill(dr_center, kt_center, content);
                }
                else { 
                    hStat->Fill(dr_center, kt_center, 0.); 
                    hLJP->Fill(dr_center, kt_center, 0.);
                }
            }
        }

        TString outputFileLJP = Form("%sLundSlice_pTbin%d", outputDirLJP.Data(), ipt);
        TString outputFileStat = Form("%sLundSlice_pTbin%d", outputDirStat.Data(), ipt);

        Draw2DLund(hStat, outputFileStat, ipt, title, "Statistical Uncertainties");
        Draw2DLund(hLJP, outputFileLJP, ipt, title, "Primary Lund Jet Plane");

        delete hStat;
        delete hLJP;
    }

    file->Close();
}

int main() {
    plot_LJP_2D();
    return 0;
}
