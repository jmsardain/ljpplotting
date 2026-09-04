#ifndef HELPERFUNCTIONS_H
#define HELPERFUNCTIONS_H

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <vector> 
#include <TLegend.h>
#include <iostream>
#include <TColor.h>
#include <TMath.h>
#include <utility>
#include <cmath>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>
#include <AtlasLabels.C>
#include <AtlasStyle.C>
#include <AtlasUtils.C>
#include <TFile.h>
#include <TGraphAsymmErrors.h>


struct SystEntry {
    std::string group;
    std::string up;
    std::string nominal;
    std::string down;
    int groupID;
    std::string type;
    std::string colorHex;
};


inline TLegend* makeLegends(const std::vector<TH1D*>& hists,
                     const std::vector<std::string>& errNames,
                     double padHeight,
                     double xPos,
                     double yPos,
                     double dX = 0.05,
                     int nColumns = 3)
{
    double Ystart = yPos;

    double spacing = 0.011 / padHeight;

    int nRows = std::ceil((double)hists.size() / nColumns);
    double dY = nRows * spacing / padHeight;
    dX = 0.97 - xPos;
    // double spacing = 0.05; 

    // int nRows = std::ceil((double)hists.size() / nColumns);
    // double dY = nRows * spacing;

    double yLow = Ystart - dY;
    if (yLow < 0.2)  yLow = 0.16;

    double textSize = 0.031 / padHeight;

    TLegend* leg = new TLegend(xPos, yLow, xPos + dX, Ystart, "", "brNDC");
    // TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9, "", "brNDC");

    leg->SetNColumns(nColumns);
    leg->SetBorderSize(0);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    // leg->SetTextSize(textSize);

    for (size_t i = 0; i < hists.size(); ++i) {

        std::string title = hists[i]->GetTitle();
        const std::string& errName = errNames[i];
        leg->AddEntry(hists[i], errName.c_str(), "l");
    }

    return leg;

}


inline void configureHist(TH1D* h, int color, int markerStyle = 20, double fillAlpha = 0.0)
{
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    h->SetMarkerStyle(markerStyle);
    h->SetMarkerSize(1.2);
    h->SetLineWidth(2);

    if (fillAlpha > 0)
        h->SetFillColorAlpha(color, fillAlpha);
}

inline void configureHistFromOther(TH1D* h, TH1D* ref)
{
    h->SetLineColor(ref->GetLineColor());
    h->SetMarkerColor(ref->GetMarkerColor());
    h->SetMarkerStyle(ref->GetMarkerStyle());
    h->SetMarkerSize(ref->GetMarkerSize());
    h->SetLineWidth(ref->GetLineWidth());
}

inline void configurePlot(TH1D* h, double padHeight)
{
    h->GetYaxis()->SetLabelSize(0.035 / padHeight);
    h->GetYaxis()->SetTitleSize(0.035 / padHeight);
    h->GetYaxis()->SetTitleOffset(2.0 * padHeight);
    h->GetYaxis()->SetTickLength(0.02);

    h->GetXaxis()->SetLabelSize(0.035 / padHeight);
    h->GetXaxis()->SetTitleSize(0.035 / padHeight);
    h->GetXaxis()->SetTitleOffset(1.3);
    h->GetXaxis()->SetTickLength(0.02 / padHeight);
}

inline void configurePlot(TGraphAsymmErrors* h, double padHeight)
{
    h->GetYaxis()->SetLabelSize(0.035 / padHeight);
    h->GetYaxis()->SetTitleSize(0.035 / padHeight);
    h->GetYaxis()->SetTitleOffset(2.0 * padHeight);
    h->GetYaxis()->SetTickLength(0.02);

    h->GetXaxis()->SetLabelSize(0.035 / padHeight);
    h->GetXaxis()->SetTitleSize(0.035 / padHeight);
    h->GetXaxis()->SetTitleOffset(1.3);
    h->GetXaxis()->SetTickLength(0.02 / padHeight);
}

inline void makeLabels(const std::vector<std::string>& labels,
                const std::string& atlasLabel,
                double padHeight,
                double xPos)
{
    double textSize = 0.03 / padHeight;

    gStyle->SetTextSize(textSize);

    double Y = 0.82;
    double delta = textSize * 1.1;

    // ATLAS label
    // ATLASLabel(xPos, Y, kBlack, 0.12, textSize, atlasLabel.c_str());
    ATLASLabel(xPos, Y,"Preliminary"); 

    Y -= delta;

    // additional labels
    for (const auto& label : labels) {
        // myText(xPos, Y, kBlack, textSize, label.c_str());
        myText(xPos, Y, kBlack, label.c_str());
        Y -= delta;
    }
}


inline TH1D* loadHist(const std::string& path, const std::string& tag, TString histName) {

    std::string fullPath = path + tag + "/4iter/Udata.root";
    if (histName == "UData_TotStat") { 
        fullPath = path + tag + "/4iter/StatUnc.root";
    }
    std::cout << "Loading histogram " << histName << " from: " << path << "/" << tag << std::endl;
    TFile* f = TFile::Open(fullPath.c_str());
    TH1D* h = (TH1D*)f->Get(histName);
    h->SetDirectory(0);
    f->Close();
    return (TH1D*)h->Clone();
}

inline void applyStatUncertainty(TH1D* h_data, const std::string& path, const std::string& tag) {
    /* cette fonction est obselete ! */
    std::string fullPath = path + tag + "/4iter/StatUnc.root";
    std::cout << "Loading stat uncertainty from: " << fullPath << std::endl;

    TFile* f = TFile::Open(fullPath.c_str());
    TH1D* h_err = (TH1D*)f->Get("UData_TotStat");
    h_err->SetDirectory(0);
    f->Close();
    for (int b = 1; b <= h_data->GetNbinsX(); ++b) {
        h_data->SetBinError(b, h_err->GetBinContent(b));
    }

    delete h_err;
}

inline TH1D* computeRel(TH1D* h_nom, TH1D* h_up, TH1D* h_down, const std::string& type) {

    TH1D* h_rel = (TH1D*)h_nom->Clone();
    h_rel->Reset();

    for (int i = 1; i <= h_nom->GetNbinsX(); ++i) {

        double nom = h_nom->GetBinContent(i);
        if (nom == 0) continue;

        double unc = 0;

        if (type == "ENVELOPE") {
            double up = h_up->GetBinContent(i);
            unc = std::fabs(up - nom) / nom;
        }
        else if (type == "SYMM") {
            double up = h_up->GetBinContent(i);
            double down = h_down->GetBinContent(i);
            unc = std::fabs(up - down) / (2.0 * nom);
        }

        h_rel->SetBinContent(i, unc);
    }

    return h_rel;
}

inline std::vector<std::vector<TH1D*>> extract_kT_distributions(TH1D* h_global, TString type = "nominal") {

    std::vector<double> kT_bins = {-10,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,4.0,10.};

    std::vector<double> pt_bins = {75,118,183,277,415,657,1034,2011,4256};
    std::vector<double> dr_bins = {0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5};

    
    std::cout << "Extracting kT distributions for type: " << type << std::endl;
    
    int NkT = kT_bins.size() - 1;

    int NdeltaR = 10;
    int NpT = 8;

    std::vector<std::vector<TH1D*>> result;

    for (int ipt = 0; ipt < NpT; ++ipt) {

        std::vector<TH1D*> vec_dR;

        for (int idr = 0; idr < NdeltaR; ++idr) {

            std::string name = "h_kT_pt" + std::to_string(ipt) + "_dR" + std::to_string(idr);
            TString finalName = type + "_" + name.c_str();
            std::string title = Form("%.0f<p_{T}<%.0f, %.1f<#DeltaR<%.1f",
                            pt_bins[ipt], pt_bins[ipt+1],
                            dr_bins[idr], dr_bins[idr+1]);

            TH1D* h_kT = new TH1D(finalName, "", NkT, &kT_bins[0]);

            for (int ikt = 0; ikt < NkT; ++ikt) {

                int global_bin = ipt * (NdeltaR * NkT)
                               + idr * NkT
                               + ikt + 1; 

                double val = h_global->GetBinContent(global_bin);
                double err = h_global->GetBinError(global_bin);
                h_kT->SetBinContent(ikt + 1, val);
                h_kT->SetBinError(ikt + 1, err);
                // if (type=="DATA") {std::cout << " " << err << std::endl;}

                // h_kT->SetTitle(title.c_str());
                h_kT->GetXaxis()->SetTitle("ln(k_{T} / GeV)");
                // h_kT->GetYaxis()->SetTitle("Uncertainty");
                

            }

            vec_dR.push_back(h_kT);
        }

        result.push_back(vec_dR);
    }

    return result;
}

// convert a TH1D into a TGraphAsymmErrors (symmetric errors from bin errors)
inline TGraphAsymmErrors* TH1DtoTGraphAsymmErrors(TH1D* h) {
    int n = h->GetNbinsX();
    TGraphAsymmErrors* g = new TGraphAsymmErrors(n);
    for (int b = 1; b <= n; ++b) {
        double x   = h->GetBinCenter(b);
        double y   = h->GetBinContent(b);
        double exl = h->GetBinWidth(b) / 2.;
        double exh = exl;
        double eyl = h->GetBinError(b);
        double eyh = eyl;
        g->SetPoint(b - 1, x, y);
        g->SetPointError(b - 1, exl, exh, eyl, eyh);
    }
    return g;
}

// ratio graph g_num / g_den, errors propagated from g_num only 
inline TGraphAsymmErrors* divideGraphs(TGraphAsymmErrors* g_num, TGraphAsymmErrors* g_den) {
    int n = g_num->GetN();
    TGraphAsymmErrors* g_ratio = new TGraphAsymmErrors(n);
    for (int i = 0; i < n; ++i) {
        double x, y_num, dummy, y_den;
        g_num->GetPoint(i, x, y_num);
        g_den->GetPoint(i, dummy, y_den);

        double exl = g_num->GetErrorXlow(i);
        double exh = g_num->GetErrorXhigh(i);
        double eyl = g_num->GetErrorYlow(i);
        double eyh = g_num->GetErrorYhigh(i);

        double ratio = (y_den != 0) ? y_num / y_den : 0;
        double reyl  = (y_den != 0) ? eyl / y_den   : 0;
        double reyh  = (y_den != 0) ? eyh / y_den   : 0;

        g_ratio->SetPoint(i, x, ratio);
        g_ratio->SetPointError(i, exl, exh, reyl, reyh);
    }
    return g_ratio;
}

// stat-error band for data, centered at 1, relative error from data itself 
inline TGraphAsymmErrors* makeStatBand(TGraphAsymmErrors* g_data) {
    int n = g_data->GetN();
    TGraphAsymmErrors* g_band = new TGraphAsymmErrors(n);
    for (int i = 0; i < n; ++i) {
        double x, y;
        g_data->GetPoint(i, x, y);
        double exl = g_data->GetErrorXlow(i);
        double exh = g_data->GetErrorXhigh(i);
        double eyl = g_data->GetErrorYlow(i);
        double eyh = g_data->GetErrorYhigh(i);

        double rel_l = (y != 0) ? eyl / y : 0;
        double rel_h = (y != 0) ? eyh / y : 0;

        g_band->SetPoint(i, x, 1.0);
        g_band->SetPointError(i, exl, exh, rel_l, rel_h);
    }
    return g_band;
}

inline std::pair<TGraphAsymmErrors*, TGraphAsymmErrors*> makeSymmetric(TGraphAsymmErrors* graph, TString color) {
    int n = graph->GetN();

    TGraphAsymmErrors* g_up   = new TGraphAsymmErrors(n);
    TGraphAsymmErrors* g_down = new TGraphAsymmErrors(n);

    for (int i = 0; i < n; ++i) {
        double x, y;
        graph->GetPoint(i, x, y);
        double exl = graph->GetErrorXlow(i);
        double exh = graph->GetErrorXhigh(i);
        double eyl = graph->GetErrorYlow(i);
        double eyh = graph->GetErrorYhigh(i);

        double posVal =  y; 
        double negVal = -y; 
        g_up->SetPoint(i, x, posVal);
        g_up->SetPointError(i, exl, exh, 0, 0);
        g_down->SetPoint(i, x, negVal);
        g_down->SetPointError(i, exl, exh, 0, 0);
    }
    g_up->SetLineColor(TColor::GetColor(color));
    g_down->SetLineColor(TColor::GetColor(color));

    g_up->SetLineColor(TColor::GetColor(color));
    g_down->SetLineColor(TColor::GetColor(color));

    g_up->SetMarkerColor(TColor::GetColor(color));
    g_down->SetMarkerColor(TColor::GetColor(color));
    

    return std::make_pair(g_up, g_down);
}

inline std::pair<TH1D*, TH1D*> makeSymmetric(TH1D*h, TString type){ 
    
    TH1D* h_up = (TH1D*)h->Clone();
    TH1D* h_down = (TH1D*)h->Clone();

    for (int b = 1; b < h->GetNbinsX()+1; ++b) {
        double val = h->GetBinContent(b);
        h_up->SetBinContent(b, +val);
        h_down->SetBinContent(b, -val);
    }

    TString color = ""; 
    int linestyle = -1;
    
    if (type == "total")   {  color = "#cccccc"; linestyle = 1;  }
    if (type == "tracking"){  color = "#d62728"; linestyle = 2;  }
    if (type == "purw")    {  color = "#8B4000"; linestyle = 3;  }
    if (type == "modeling"){  color = "#FFA62B"; linestyle = 4;  }
    if (type == "jes")     {  color = "#1f77b4"; linestyle = 5;  }
    
    h_up->SetLineColor(TColor::GetColor(color));
    h_down->SetLineColor(TColor::GetColor(color));

    h_up->SetLineColor(TColor::GetColor(color));
    h_down->SetLineColor(TColor::GetColor(color));

    h_up->SetMarkerColor(TColor::GetColor(color));
    h_down->SetMarkerColor(TColor::GetColor(color));

    h_up->SetLineStyle(linestyle); 
    h_down->SetLineStyle(linestyle); 
    if (type == "total") { 
        h_up->SetFillColorAlpha(kGray+1, 0.4);
        h_down->SetFillColorAlpha(kGray+1, 0.4);
    }
    return std::make_pair(h_up, h_down);

}




inline void DrawTH1DPlots(
    const std::vector<std::vector<TH1D*>>& kt_data,
    const std::vector<std::vector<TH1D*>>& kt_pythia,
    const std::vector<std::vector<TH1D*>>& kt_powheg,
    const std::vector<std::vector<TH1D*>>& kt_sherpa,
    const std::vector<std::vector<TH1D*>>& kt_herwig
)
{

    TCanvas*c = new TCanvas("", "", 600, 600); 
    c->SetLogy(); 
    c->SaveAs("kt_1D_dataMC.pdf[");


    std::vector<double> kT_bins = {-10,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,4.0,10.};
    std::vector<double> pt_bins = {75,118,183,277,415,657,1034,2011,4256};
    std::vector<double> dr_bins = {0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5};

    int NkT = kT_bins.size() - 1;
    int NdeltaR = 10;
    int NpT = 8;

    TLegend* legUpper = new TLegend(0.8, 0.6, 0.88, 0.9);
    legUpper->SetBorderSize(0);
    legUpper->SetFillStyle(0);

    double ratioHeight = 0.3;

    for (int ipt = 1; ipt < NpT; ++ipt) {
        for (int idr = 0; idr < NdeltaR; ++idr) {
           
            TCanvas* c1 = new TCanvas(
                Form("c_pt%d_dR%d", ipt, idr),
                "", 800, 800
            );

            TPad* upperPad = new TPad("pad1","",0.0,1.2*ratioHeight,1.0,1.0);
            TPad* ratioPad = new TPad("pad2","",0.0,0.1,1.0,1.2*ratioHeight);


            upperPad->SetBottomMargin(0);
            upperPad->SetLeftMargin(0.15);
            upperPad->SetRightMargin(0.05);

            ratioPad->SetBottomMargin(0);
            ratioPad->SetTopMargin(0);
            ratioPad->SetLeftMargin(0.15);
            ratioPad->SetRightMargin(0.05);

            upperPad->Draw();
            ratioPad->Draw();


            upperPad->cd();

            
            TGraphAsymmErrors* g_data   = TH1DtoTGraphAsymmErrors(kt_data[ipt][idr]);
            TGraphAsymmErrors* g_pythia = TH1DtoTGraphAsymmErrors(kt_pythia[ipt][idr]);
            TGraphAsymmErrors* g_powheg = TH1DtoTGraphAsymmErrors(kt_powheg[ipt][idr]);
            TGraphAsymmErrors* g_sherpa = TH1DtoTGraphAsymmErrors(kt_sherpa[ipt][idr]);
            TGraphAsymmErrors* g_herwig = TH1DtoTGraphAsymmErrors(kt_herwig[ipt][idr]);

            g_pythia->SetLineColor(TColor::GetColor("#5520c7"));
            g_powheg->SetLineColor(TColor::GetColor("#473868"));
            g_sherpa->SetLineColor(TColor::GetColor("#FF6F61"));
            g_herwig->SetLineColor(TColor::GetColor("#015457"));
            g_data->SetLineColor(TColor::GetColor("#000000"));
            
            
            g_pythia->SetMarkerStyle(21);
            g_powheg->SetMarkerStyle(46);
            g_sherpa->SetMarkerStyle(22);
            g_herwig->SetMarkerStyle(34);
            g_data->SetMarkerStyle(20);

            g_pythia->SetMarkerColor(TColor::GetColor("#5520c7"));
            g_powheg->SetMarkerColor(TColor::GetColor("#473868"));
            g_sherpa->SetMarkerColor(TColor::GetColor("#FF6F61"));
            g_herwig->SetMarkerColor(TColor::GetColor("#015457"));
            g_data->SetMarkerColor(TColor::GetColor("#000000"));

            if (ipt == 0 && idr == 0) {
                legUpper->AddEntry(g_data,   "Data",   "lp");
                legUpper->AddEntry(g_pythia, "Pythia", "p");
                legUpper->AddEntry(g_powheg, "Powheg", "p");
                legUpper->AddEntry(g_sherpa, "Sherpa", "p");
                legUpper->AddEntry(g_herwig, "Herwig", "p");
            }

            g_data->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#sigma}{dln(k_{T})}");
            g_data->GetXaxis()->SetTitle("ln(k_{T} / GeV)");
            g_data->GetXaxis()->SetLimits(-10, 10);
            // g_data->SetMinimum(1e-4);   
            // g_data->SetMaximum(0.3);     
            double errorPadHeight = ratioHeight + 0.1;
            configurePlot(g_data, errorPadHeight);

            g_data->Draw("AP");
            g_pythia->Draw("P SAME");
            g_powheg->Draw("P SAME");
            g_sherpa->Draw("P SAME");
            g_herwig->Draw("P SAME");
            legUpper->Draw("SAME");

            std::vector<std::string> labels;
            labels.push_back("#sqrt{s} = 13 TeV, 140 fb^{-1}");
            labels.push_back(Form("%.1f < ln(R/#DeltaR) < %.1f", dr_bins[idr], dr_bins[idr+1]));
            labels.push_back(Form("%.0f < p_{T} / GeV < %.0f", pt_bins[ipt], pt_bins[ipt+1]));
            makeLabels(labels, "Internal", 1.0 - 2*ratioHeight, 0.2);

    

            ratioPad->cd();
            
            TGraphAsymmErrors* g_ratio_pythia = divideGraphs(g_pythia, g_data);
            TGraphAsymmErrors* g_ratio_powheg = divideGraphs(g_powheg, g_data);
            TGraphAsymmErrors* g_ratio_sherpa = divideGraphs(g_sherpa, g_data);
            TGraphAsymmErrors* g_ratio_herwig = divideGraphs(g_herwig, g_data);

            g_ratio_pythia->SetMarkerStyle(21); g_ratio_pythia->SetMarkerColor(TColor::GetColor("#5520c7")); g_ratio_pythia->SetLineColor(TColor::GetColor("#5520c7"));
            g_ratio_powheg->SetMarkerStyle(46); g_ratio_powheg->SetMarkerColor(TColor::GetColor("#473868")); g_ratio_powheg->SetLineColor(TColor::GetColor("#473868"));
            g_ratio_sherpa->SetMarkerStyle(22); g_ratio_sherpa->SetMarkerColor(TColor::GetColor("#FF6F61")); g_ratio_sherpa->SetLineColor(TColor::GetColor("#FF6F61"));
            g_ratio_herwig->SetMarkerStyle(34); g_ratio_herwig->SetMarkerColor(TColor::GetColor("#015457")); g_ratio_herwig->SetLineColor(TColor::GetColor("#015457"));


            // Calculate statistical error band for data to be plotted in ratio plot
            TGraphAsymmErrors* g_band = makeStatBand(g_data);
            g_band->SetFillColor(TColor::GetColor("#cccccc"));
            g_band->SetFillStyle(1001);     
            g_band->SetLineColor(TColor::GetColor("#cccccc"));
            g_band->SetMarkerSize(0);

            g_band->GetXaxis()->SetLimits(-10, 10);
            g_band->GetYaxis()->SetTitle("#splitline{Ratio to}{Data}");
            g_band->GetYaxis()->SetRangeUser(0., 2.);
            g_band->GetYaxis()->SetNdivisions(505);
            configurePlot(g_band, (2*ratioHeight - (ratioHeight + 0.1)));


            g_band->Draw("A2");         // axis + filled asymmetric band
            g_ratio_pythia->Draw("P SAME");
            g_ratio_powheg->Draw("P SAME");
            g_ratio_sherpa->Draw("P SAME");
            g_ratio_herwig->Draw("P SAME");

            TLine* line = new TLine(-10, 1, 10, 1);
            line->SetLineStyle(2);
            line->Draw();

            c1->SaveAs("kt_1D_dataMC.pdf");
            delete c1;

        }
    }
    c->SaveAs("kt_1D_dataMC.pdf]");
}




inline void DrawTH1DPlotsWithError(
    // const std::map<std::string, int>& groupColor, 
    const std::vector<std::vector<TH1D*>>& kt_data,
    const std::vector<std::vector<TH1D*>>& kt_pythia,
    const std::vector<std::vector<TH1D*>>& kt_powheg,
    const std::vector<std::vector<TH1D*>>& kt_sherpa,
    const std::vector<std::vector<TH1D*>>& kt_herwig,
    const std::vector<std::vector<TH1D*>>& kt_total,
    const std::vector<std::vector<TH1D*>>& kt_tracking,
    const std::vector<std::vector<TH1D*>>& kt_purw,
    const std::vector<std::vector<TH1D*>>& kt_JES,
    const std::vector<std::vector<TH1D*>>& kt_Modeling
    )
{
    int NpT = kt_pythia.size();
    int NdeltaR = kt_pythia[0].size();

    std::vector<double> kT_bins = {-10,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,4.0,10.};
    std::vector<double> dr_bins = {0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5};
    std::vector<double> pt_bins = {75,118,183,277,415,657,1034,2011,4256};

    TLegend* legUpper = new TLegend(0.7, 0.5, 0.9, 0.9);
    legUpper->SetBorderSize(0);
    legUpper->SetFillStyle(0);


    double ratioHeight = 0.3;

    TCanvas* c = new TCanvas("c_summary", "kT Summary", 800, 800);
    c->SaveAs("kt_1Derror_dataMC.pdf[");

    for (int ipt = 1; ipt < NpT; ++ipt) {
        for (int idr = 0; idr < NdeltaR; ++idr) {


            TGraphAsymmErrors* g_data   = TH1DtoTGraphAsymmErrors(kt_data[ipt][idr]);
            TGraphAsymmErrors* g_pythia = TH1DtoTGraphAsymmErrors(kt_pythia[ipt][idr]);
            TGraphAsymmErrors* g_powheg = TH1DtoTGraphAsymmErrors(kt_powheg[ipt][idr]);
            TGraphAsymmErrors* g_sherpa = TH1DtoTGraphAsymmErrors(kt_sherpa[ipt][idr]);
            TGraphAsymmErrors* g_herwig = TH1DtoTGraphAsymmErrors(kt_herwig[ipt][idr]);

            TGraphAsymmErrors* g_total    = TH1DtoTGraphAsymmErrors(kt_total[ipt][idr]);
            TGraphAsymmErrors* g_tracking = TH1DtoTGraphAsymmErrors(kt_tracking[ipt][idr]);
            TGraphAsymmErrors* g_purw     = TH1DtoTGraphAsymmErrors(kt_purw[ipt][idr]);
            TGraphAsymmErrors* g_modeling = TH1DtoTGraphAsymmErrors(kt_Modeling[ipt][idr]);
            TGraphAsymmErrors* g_jes      = TH1DtoTGraphAsymmErrors(kt_JES[ipt][idr]);
            
            // auto [g_total_up, g_total_down]       = makeSymmetric(g_total, "#cccccc");
            // auto [g_tracking_up, g_tracking_down] = makeSymmetric(g_tracking, "#d62728");
            // auto [g_purw_up, g_purw_down]         = makeSymmetric(g_purw, "#8B4000");
            // auto [g_modeling_up, g_modeling_down] = makeSymmetric(g_modeling, "#FFA62B");
            // auto [g_jes_up, g_jes_down]           = makeSymmetric(g_jes, "#1f77b4");
            
            auto [g_total_up, g_total_down]       = makeSymmetric(kt_total[ipt][idr], "total");
            auto [g_tracking_up, g_tracking_down] = makeSymmetric(kt_tracking[ipt][idr], "tracking");
            auto [g_purw_up, g_purw_down]         = makeSymmetric(kt_purw[ipt][idr], "purw");
            auto [g_modeling_up, g_modeling_down] = makeSymmetric(kt_Modeling[ipt][idr], "modeling");
            auto [g_jes_up, g_jes_down]           = makeSymmetric(kt_JES[ipt][idr], "jes");

            // if (ipt ==1 && idr ==1 ) 
            //     std::cout << " data  " << g_data->GetN() 
            //               <<  " pythia " << g_pythia->GetN() 
            //               << " total " << g_total->GetN()
            //               << " g_total_up " << g_total_up->GetN()
            //               << " modeling " << g_modeling->GetN()
            //               << " g_modeling_up " << g_modeling_up->GetN()
            //               << std::endl;
        
            TCanvas* c1 = new TCanvas(
                Form("c_pt%d_dR%d", ipt, idr),
                "", 800, 800
            );

            // ---------------------------
            // Pads (match Python layout)
            // ---------------------------
            TPad* upperPad = new TPad("pad1","",0.0,2*ratioHeight,1.0,1.0);
            TPad* ratioPad = new TPad("pad2","",0.0,ratioHeight+0.1,1.0,2*ratioHeight);
            TPad* errorPad = new TPad("pad3","",0.0,0.0,1.0,ratioHeight+0.1);

            upperPad->SetBottomMargin(0);
            upperPad->SetLeftMargin(0.15);
            upperPad->SetRightMargin(0.05);

            ratioPad->SetBottomMargin(0);
            ratioPad->SetTopMargin(0);
            ratioPad->SetLeftMargin(0.15);
            ratioPad->SetRightMargin(0.05);

            errorPad->SetBottomMargin(0.35);
            errorPad->SetTopMargin(0);
            errorPad->SetLeftMargin(0.15);
            errorPad->SetRightMargin(0.05);

            upperPad->Draw();
            ratioPad->Draw();
            errorPad->Draw();

            upperPad->cd();

            g_pythia->SetLineColor(TColor::GetColor("#5520c7"));
            g_powheg->SetLineColor(TColor::GetColor("#473868"));
            g_sherpa->SetLineColor(TColor::GetColor("#FF6F61"));
            g_herwig->SetLineColor(TColor::GetColor("#015457"));
            g_data->SetLineColor(TColor::GetColor("#000000"));
            g_data->SetFillColor(TColor::GetColor("#B0B0B0"));
            g_data->SetLineColor(TColor::GetColor("#B0B0B0"));
            
            
            g_pythia->SetMarkerStyle(21);
            g_powheg->SetMarkerStyle(46);
            g_sherpa->SetMarkerStyle(22);
            g_herwig->SetMarkerStyle(34);
            g_data->SetMarkerStyle(20);

            g_pythia->SetMarkerColor(TColor::GetColor("#5520c7"));
            g_powheg->SetMarkerColor(TColor::GetColor("#473868"));
            g_sherpa->SetMarkerColor(TColor::GetColor("#FF6F61"));
            g_herwig->SetMarkerColor(TColor::GetColor("#015457"));
            g_data->SetMarkerColor(TColor::GetColor("#000000"));

            if (ipt == 1 && idr == 0) {
                legUpper->AddEntry(g_data,   "Data",   "fep");
                legUpper->AddEntry(g_pythia, "Pythia", "p");
                legUpper->AddEntry(g_powheg, "Powheg+Pythia", "p");
                legUpper->AddEntry(g_sherpa, "Sherpa (AHADIC)", "p");
                legUpper->AddEntry(g_herwig, "Herwig (AngOrd)", "p");
            }

            g_data->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#sigma}{dln(k_{T})}");
            g_data->GetXaxis()->SetTitle("ln(k_{T} / GeV)");
            g_data->GetXaxis()->SetLimits(-10, 10);
            // g_data->SetMinimum(1e-4);   
            // g_data->SetMaximum(0.3);     
            double errorPadHeight = ratioHeight + 0.1;
            configurePlot(g_data, errorPadHeight);
            g_data->Draw("AP");
            g_data->Draw("E2 SAME");
            g_pythia->Draw("P SAME");
            g_powheg->Draw("P SAME");
            g_sherpa->Draw("P SAME");
            g_herwig->Draw("P SAME");
            legUpper->Draw("SAME");

            std::vector<std::string> labels;
            labels.push_back("#sqrt{s} = 13 TeV, 140 fb^{-1}");
            labels.push_back(Form("%.1f < ln(R/#DeltaR) < %.1f", dr_bins[idr], dr_bins[idr+1]));
            labels.push_back(Form("%.0f < p_{T} / GeV < %.0f", pt_bins[ipt], pt_bins[ipt+1]));
            makeLabels(labels, "Internal", 1.0 - 2*ratioHeight, 0.2);

    

            ratioPad->cd();
            
            TGraphAsymmErrors* g_ratio_pythia = divideGraphs(g_pythia, g_data);
            TGraphAsymmErrors* g_ratio_powheg = divideGraphs(g_powheg, g_data);
            TGraphAsymmErrors* g_ratio_sherpa = divideGraphs(g_sherpa, g_data);
            TGraphAsymmErrors* g_ratio_herwig = divideGraphs(g_herwig, g_data);

            g_ratio_pythia->SetMarkerStyle(21); g_ratio_pythia->SetMarkerColor(TColor::GetColor("#5520c7")); g_ratio_pythia->SetLineColor(TColor::GetColor("#5520c7"));
            g_ratio_powheg->SetMarkerStyle(46); g_ratio_powheg->SetMarkerColor(TColor::GetColor("#473868")); g_ratio_powheg->SetLineColor(TColor::GetColor("#473868"));
            g_ratio_sherpa->SetMarkerStyle(22); g_ratio_sherpa->SetMarkerColor(TColor::GetColor("#FF6F61")); g_ratio_sherpa->SetLineColor(TColor::GetColor("#FF6F61"));
            g_ratio_herwig->SetMarkerStyle(34); g_ratio_herwig->SetMarkerColor(TColor::GetColor("#015457")); g_ratio_herwig->SetLineColor(TColor::GetColor("#015457"));


            // Calculate statistical error band for data to be plotted in ratio plot
            TGraphAsymmErrors* g_band = makeStatBand(g_data);
            g_band->SetFillColor(TColor::GetColor("#cccccc"));
            g_band->SetFillStyle(1001);     
            g_band->SetLineColor(TColor::GetColor("#cccccc"));
            g_band->SetMarkerSize(0);

            g_band->GetXaxis()->SetLimits(-10, 10);
            g_band->GetYaxis()->SetTitle("#splitline{Ratio to}{Data}");
            g_band->GetYaxis()->SetRangeUser(0., 2.);
            g_band->GetYaxis()->SetNdivisions(505);
            configurePlot(g_band, (2*ratioHeight - (ratioHeight + 0.1)));


            g_band->Draw("A2");         // axis + filled asymmetric band
            g_ratio_pythia->Draw("P SAME");
            g_ratio_powheg->Draw("P SAME");
            g_ratio_sherpa->Draw("P SAME");
            g_ratio_herwig->Draw("P SAME");

            TLine* line = new TLine(-10, 1, 10, 1);
            line->SetLineStyle(2);
            line->Draw();

            // now do the error breakdown pad
            errorPad->cd();

            g_total_up->GetYaxis()->SetTitle("#splitline{Relative}{Uncertainty}");
            // tot_up->GetYaxis()->SetRangeUser(-0.2,0.2);
            g_total_up->GetXaxis()->SetTitle("ln(k_{T}/GeV)");
            configurePlot(g_total_up, errorPadHeight);
            
            // tot_up->SetLineColor(kBlack);
            // tot_down->SetLineColor(kBlack);
            
            // tot_up->Draw("axis");
            // TH1D* h_band_err = (TH1D*)kt_total[ipt][idr]->Clone("h_band_err");
            // double maxUnc = 0;

            // for (int b = 1; b < h_band_err->GetNbinsX()+1; ++b) {
            //     double unc = kt_total[ipt][idr]->GetBinContent(b);
            //     h_band_err->SetBinContent(b, 0.0);
            //     h_band_err->SetBinError(b, unc);
            //     maxUnc = std::max(maxUnc, unc);
            // }

            // h_band_err->SetFillColorAlpha(kGray+1, 0.4);
            // h_band_err->SetLineColor(kBlack);

            // h_band_err->GetYaxis()->SetTitle("#splitline{Relative}{Uncertainty}");
            // h_band_err->GetYaxis()->SetRangeUser(-1.2*maxUnc, 1.2*maxUnc);
            // h_band_err->SetMarkerSize(0); 
            // h_band_err->SetMarkerStyle(0); 

    
            std::vector<TH1D*> systUp = {g_total_up, g_tracking_up, g_purw_up, g_jes_up, g_modeling_up};
            std::vector<std::string> errNames = {"Total Syst.", "Tracking", "PURW", "JES+JER", "MC Model"};

            TLegend*legError =  makeLegends(systUp, errNames, 0.4, 0.18, 0.92, 0.6, 5);

            g_total_up->SetMinimum(-2);
            g_total_up->SetMaximum(2);
            // g_total_up->Draw("A");
            g_total_up->Draw("HIST"); g_total_down->Draw("HIST SAME");
            g_tracking_up->Draw("HIST SAME"); g_tracking_down->Draw("HIST SAME");
            g_purw_up->Draw("HIST SAME"); g_purw_down->Draw("HIST SAME");
            g_modeling_up->Draw("HIST SAME"); g_modeling_down->Draw("HIST SAME");
            g_jes_up->Draw("HIST SAME"); g_jes_down->Draw("HIST SAME");
            legError->Draw("SAME");


            c1->SaveAs("kt_1Derror_dataMC.pdf");
            delete c1;
        }
    }

    c->SaveAs("kt_1Derror_dataMC.pdf]");
}




void DrawTH1DPlots( TH1D* herr_total,
                    TH1D* herr_tracking,
                    TH1D* herr_purw,
                    TH1D* herr_JES,
                    TH1D* herr_Modeling, 
                    double ymin, 
                    double ymax
                    )
    {


    TCanvas* c = new TCanvas("", "", 4000, 800);
    TH1D* h_band_err = (TH1D*)herr_total->Clone("h_band");
    double maxUnc = 0;

    for (int b = 1; b < h_band_err->GetNbinsX()+1; ++b) {
        double unc = herr_total->GetBinContent(b);
        h_band_err->SetBinContent(b, 0.0);
        h_band_err->SetBinError(b, unc);
        maxUnc = std::max(maxUnc, unc);
    }

    h_band_err->SetFillColorAlpha(kGray+1, 0.4);
    h_band_err->SetLineColor(kGray+1);

    TLine *line2 = new TLine(h_band_err->GetBinLowEdge(1),0.0,h_band_err->GetBinLowEdge(h_band_err->GetNbinsX()+1),0.0);
    line2->SetLineStyle(2); // line at 0
    line2->SetLineWidth(1); // line at 0

    h_band_err->GetYaxis()->SetNdivisions(505);
    h_band_err->GetYaxis()->SetTitle("");
    h_band_err->GetYaxis()->SetTitle("#splitline{Relative}{uncertainty}");
    h_band_err->GetYaxis()->SetTitleSize(0);
    line2->Draw("same");

    h_band_err->GetYaxis()->SetRangeUser(-1.2*maxUnc, 1.2*maxUnc);
    h_band_err->SetMarkerSize(0); 
    h_band_err->SetMarkerStyle(0); 

    // configurePlot(h_band_err, errorPadHeight);
    // h_band_err->GetXaxis()->SetRangeUser(150,1200);
    h_band_err->GetYaxis()->SetRangeUser(ymin, ymax);
    h_band_err->Draw("HIST E2");


    auto [tot_up, tot_down] = makeSymmetric(herr_total, "total");
    auto [trk_up, trk_down] = makeSymmetric(herr_tracking, "tracking");
    auto [purw_up, purw_down] = makeSymmetric(herr_purw, "purw");
    auto [jes_up, jes_down] = makeSymmetric(herr_JES, "jes");
    auto [modeling_up, modeling_down] = makeSymmetric(herr_Modeling, "modeling");

    trk_up->Draw("hist same");   trk_down->Draw("hist same");
    purw_up->Draw("hist same");  purw_down->Draw("hist same");
    jes_up->Draw("hist same");   jes_down->Draw("hist same");
    modeling_up->Draw("hist same");   modeling_down->Draw("hist same");

    TLegend* leg = new TLegend(0.6, 0.65, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetNColumns(5);
    leg->AddEntry(h_band_err, "Total", "l");
    leg->AddEntry(herr_Modeling, "Modeling", "l");
    leg->AddEntry(herr_purw, "PURW", "l");
    leg->AddEntry(herr_tracking, "Tracking", "l");
    leg->AddEntry(herr_JES, "JES", "l");
    // leg->AddEntry(herr_tile, "TILE", "l");
    
    leg->Draw();

    std::vector<std::string> labels;
    labels.push_back("#sqrt{s} = 13 TeV, 140 fb^{-1}");
    makeLabels(labels, "Simulation",  0.8, 0.2);

    // put pT delimeter
    TLatex latex;
    latex.SetTextAlign(22); 
    latex.SetTextSize(0.03);
    std::vector<double> pt_bins = {75,118,183,277,415,657,1034,2011,4256};
    for (int i = 1; i <= 8; ++i) {
        TLine* l = new TLine(150 * i, ymin , 150 * i, ymax);
        l->SetLineColor(kBlack);
        l->SetLineStyle(3);
        l->Draw();

        double x_center = 150 * (i - 0.5);
        // double y_pos = -1.1 * maxUnc;
        double y_pos = ymin + 0.05 * (ymax - ymin);
        TString label = Form("%.0f < p_{T} / GeV < %.0f", pt_bins[i-1], pt_bins[i]); 
        latex.DrawLatex(x_center, y_pos, label);
    }
    // put dR delimiters
    std::vector<double> dr_bins = {0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5};
    latex.SetTextAngle(90);
    for (int i = 11; i <= 20; ++i) {
        TLine* l = new TLine(15 * i, ymin+1, 15 * i,ymax-1);
        l->SetLineColor(kBlack);
        l->SetLineStyle(3);
        l->Draw();

        double x_center = 15 * (i - 0.5);
        // double y_pos = -0.4 * maxUnc;
        double y_pos = ymin + 0.05 * (ymax - ymin);
        int iterator = i-10;
        TString label = Form("%.1f < ln(R/#Delta R) < %.1f", dr_bins[iterator-1], dr_bins[iterator]); 
        latex.DrawLatex(x_center, y_pos, label);
    }

    


    c->SaveAs("totalUnc_allBins.pdf");
}

#endif

