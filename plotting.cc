#include <fstream>
#include <sstream>
#include <iostream>
#include <map>
#include <vector>
#include <cmath>

#include <AtlasLabels.C>
#include <AtlasStyle.C>
#include <AtlasUtils.C>

#include "TFile.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1.h"
#include "TStyle.h"
#include "TLine.h"

struct SystEntry {
    std::string group;
    std::string up;
    std::string nominal;
    std::string down;
    int groupID;
    std::string type;
    std::string colorHex;
};




TH1D* loadHist(const std::string& base, const std::string& tag) {
    std::string path = base + tag + "/4iter/Udata.root";
    if (TString(tag.c_str()) == "") { path =  "/data/ldelagra/LundPlane/Unfolding/nominal_PrimaryLJP/Final_Files_16-03-2026/4iter/Udata.root"; }
    std::cout << "Loading histogram from: " << path << std::endl;
    TFile* f = TFile::Open(path.c_str());

    TString histName = "gTrue_original"; 
    if (TString(tag.c_str()) == "") {  histName = "UData"; }
    TH1D* h = (TH1D*)f->Get(histName);
    h->SetDirectory(0);
    f->Close();
    return (TH1D*)h->Clone();
}


TH1D* computeRel(TH1D* h_nom, TH1D* h_up, TH1D* h_down, const std::string& type) {

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

std::vector<std::vector<TH1D*>> extract_kT_distributions(TH1D* h_global, TString type = "nominal") {

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
                h_kT->GetYaxis()->SetTitle("Uncertainty");
                

            }

            vec_dR.push_back(h_kT);
        }

        result.push_back(vec_dR);
    }

    return result;
}


void DrawTH1DFinalPlots(const std::vector<std::vector<TH1D*>>& kt_nominal,
                        const std::vector<std::vector<TH1D*>>& kt_total,
                        const std::vector<std::vector<TH1D*>>& kt_tracking,
                        const std::vector<std::vector<TH1D*>>& kt_purw,
                        const std::vector<std::vector<TH1D*>>& kt_JES){
    int NpT = kt_nominal.size();
    int NdeltaR = kt_nominal[0].size();

    std::vector<double> pt_bins = {75,118,183,277,415,657,1034,2011,4256};
    std::vector<double> dr_bins = {0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5};
    
    TCanvas* c = new TCanvas("c_summary", "kT Uncertainties Summary", 1200, 1000);
    c->SaveAs("kT_summary.pdf["); // start multi-page PDF

    for (int ipt = 0; ipt < NpT; ++ipt) {
        for (int idr = 0; idr < NdeltaR; ++idr) {
            
            TCanvas* c1 = new TCanvas( Form("c_pt%d_dR%d", ipt, idr), Form("kT pT%d dR%d", ipt, idr),800, 600);

            TPad* p1 = new TPad("p1", "top", 0, 0.3, 1, 1);
            TPad* p2 = new TPad("p2", "bottom", 0, 0, 1, 0.3);

            p1->SetBottomMargin(0.02);
            p2->SetTopMargin(0.05);
            p2->SetBottomMargin(0.3);

            p1->Draw();
            p2->Draw();
            p1->cd();

            TH1D* h_nom = kt_nominal[ipt][idr];
            h_nom->SetDirectory(0);
            h_nom->SetLineColor(kBlack);
            h_nom->SetLineWidth(2);

            // h_nom->SetTitle(Form("%.0f<p_{T}<%.0f, %.1f<#DeltaR<%.1f",
            //                      pt_bins[ipt], pt_bins[ipt+1],
            //                      dr_bins[idr], dr_bins[idr+1]));

            h_nom->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#it{#sigma}}{d#it{ln(R/#Delta R)}d#it{p}_{T}}");
            h_nom->Draw("hist");

            p2->cd();
            TH1D* h_tot = kt_total[ipt][idr];
            TH1D* h_trk = kt_tracking[ipt][idr];
            TH1D* h_purw = kt_purw[ipt][idr];
            TH1D* h_jes = kt_JES[ipt][idr];
            h_tot->SetDirectory(0);
            h_trk->SetDirectory(0);
            h_purw->SetDirectory(0);
            h_jes->SetDirectory(0);

            h_tot->SetLineColor(kBlack);
            h_trk->SetLineColor(kRed);
            h_purw->SetLineColor(kGreen);
            h_jes->SetLineColor(kBlue);
            
            h_tot->SetLineWidth(2);
            h_trk->SetLineWidth(2);
            h_jes->SetLineWidth(2);
            h_purw->SetLineWidth(2);
            

            h_tot->GetYaxis()->SetTitle("Rel. unc.");
            h_tot->GetXaxis()->SetTitle("k_{T}");

            // remove later once i start adding ATLAS style 
            h_tot->GetYaxis()->SetNdivisions(505);
            h_tot->GetYaxis()->SetTitleSize(0.08);
            h_tot->GetYaxis()->SetLabelSize(0.08);

            h_tot->GetXaxis()->SetTitleSize(0.1);
            h_tot->GetXaxis()->SetLabelSize(0.08);

            h_tot->SetMinimum(0);
            h_tot->SetMaximum(1.0); // adjust if needed

            h_tot->Draw("hist");
            h_trk->Draw("hist same");
            h_jes->Draw("hist same");
            
            TLegend* leg = new TLegend(0.6, 0.65, 0.88, 0.88);
            leg->AddEntry(h_tot, "Total", "l");
            leg->AddEntry(h_purw, "PURW", "l");
            leg->AddEntry(h_trk, "Tracking", "l");
            leg->AddEntry(h_jes, "JES", "l");
            leg->Draw();

            // c1->SaveAs(Form("kT_pt%.0f_%.0f_dR%.1f_%.1f.pdf",
            //                 pt_bins[ipt], pt_bins[ipt+1],
            //                 dr_bins[idr], dr_bins[idr+1]));
            c1->SaveAs("kT_summary.pdf"); // add page to multi-page PDF
            delete c1; 

        }
    }
    c->SaveAs("kT_summary.pdf]"); // end multi-page PDF
}



void configureHist(TH1D* h, int color, int markerStyle = 20, double fillAlpha = 0.0)
{
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    h->SetMarkerStyle(markerStyle);
    h->SetMarkerSize(1.2);
    h->SetLineWidth(2);

    if (fillAlpha > 0)
        h->SetFillColorAlpha(color, fillAlpha);
}

void configureHistFromOther(TH1D* h, TH1D* ref)
{
    h->SetLineColor(ref->GetLineColor());
    h->SetMarkerColor(ref->GetMarkerColor());
    h->SetMarkerStyle(ref->GetMarkerStyle());
    h->SetMarkerSize(ref->GetMarkerSize());
    h->SetLineWidth(ref->GetLineWidth());
}

void configurePlot(TH1D* h, double padHeight)
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

void makeLabels(const std::vector<std::string>& labels,
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

TLegend* makeLegends(const std::vector<TH1*>& hists,
                     const std::vector<std::string>& errNames,
                     double padHeight,
                     double xPos,
                     double yPos,
                     double dX = 0.05,
                     int nColumns = 1)
{
    double Ystart = yPos;

    // double spacing = 0.011 / padHeight;

    // int nRows = std::ceil((double)hists.size() / nColumns);
    // double dY = nRows * spacing / padHeight;
    // dX = 0.97 - xPos;
    double spacing = 0.05; 
    int nRows = std::ceil((double)hists.size() / nColumns);
    double dY = nRows * spacing;

    double yLow = Ystart - dY;
    if (yLow < 0.2) yLow = 0.16;

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

void DrawTH1DPlots(
    TH1D* herr_total,
    TH1D* herr_tracking,
    TH1D* herr_purw,
    TH1D* herr_JES,
    TH1D* herr_Modeling
    // ,TH1D* herr_tile
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
    h_band_err->GetXaxis()->SetRangeUser(150,1200);
    h_band_err->GetYaxis()->SetRangeUser(-5, 5);
    h_band_err->Draw("HIST E2");

    auto makeSymmetric = [&](TH1D* h){
        TH1D* h_up = (TH1D*)h->Clone();
        TH1D* h_down = (TH1D*)h->Clone();

        for (int b = 1; b < h->GetNbinsX()+1; ++b) {
            double val = h->GetBinContent(b);
            h_up->SetBinContent(b, +val);
            h_down->SetBinContent(b, -val);
        }

        return std::make_pair(h_up, h_down);
    };

    auto [tot_up, tot_down] = makeSymmetric(herr_total);
    auto [trk_up, trk_down] = makeSymmetric(herr_tracking);
    auto [purw_up, purw_down] = makeSymmetric(herr_purw);
    auto [jes_up, jes_down] = makeSymmetric(herr_JES);
    auto [modeling_up, modeling_down] = makeSymmetric(herr_Modeling);

    // auto [tile_up, tile_down] = makeSymmetric(herr_tile);
            

    trk_up->SetLineColor(TColor::GetColor("#d62728"));
    trk_down->SetLineColor(TColor::GetColor("#d62728"));
    trk_up->Draw("hist same");
    trk_down->Draw("hist same");

    purw_up->SetLineColor(TColor::GetColor("#8B4000"));
    purw_down->SetLineColor(TColor::GetColor("#8B4000"));
    purw_up->Draw("hist same");
    purw_down->Draw("hist same");

    jes_up->SetLineColor(TColor::GetColor("#1f77b4"));
    jes_down->SetLineColor(TColor::GetColor("#1f77b4"));
    jes_up->Draw("hist same");
    jes_down->Draw("hist same");


    modeling_up->SetLineColor(TColor::GetColor("#FFA62B"));
    modeling_down->SetLineColor(TColor::GetColor("#FFA62B"));
    modeling_up->Draw("hist same");
    modeling_down->Draw("hist same");

    // tile_up->SetLineColor(TColor::GetColor("#188535"));
    // tile_down->SetLineColor(TColor::GetColor("#188535"));
    // tile_up->Draw("hist same");
    // tile_down->Draw("hist same");

    TLegend* leg = new TLegend(0.6, 0.65, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetNColumns(4);
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
    for (int i = 2; i <= 8; ++i) {
        TLine* l = new TLine(150 * i, -5 , 150 * i, 5);
        l->SetLineColor(kBlack);
        l->SetLineStyle(3);
        l->Draw();

        double x_center = 150 * (i - 0.5);
        double y_pos = -1.1 * maxUnc;
        TString label = Form("%.0f < p_{T} / GeV < %.0f", pt_bins[i-1], pt_bins[i]); 
        latex.DrawLatex(x_center, y_pos, label);
    }
    // put dR delimiters
    std::vector<double> dr_bins = {0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5};
    latex.SetTextAngle(90);
    for (int i = 11; i <= 20; ++i) {
        TLine* l = new TLine(15 * i, -3, 15 * i, 3);
        l->SetLineColor(kBlack);
        l->SetLineStyle(3);
        l->Draw();

        double x_center = 15 * (i - 0.5);
        double y_pos = -0.4 * maxUnc;
        int iterator = i-10;
        TString label = Form("%.1f < ln(R/#Delta R) < %.1f", dr_bins[iterator-1], dr_bins[iterator]); 
        latex.DrawLatex(x_center, y_pos, label);
    }

    


    c->SaveAs("totalUnc_allBins.pdf");
}

void DrawTH1DPlotsWithError(
    const std::vector<std::vector<TH1D*>>& kt_pythia,
    const std::vector<std::vector<TH1D*>>& kt_data,
    const std::vector<std::vector<TH1D*>>& kt_powpy,
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

    double ratioHeight = 0.3;

    TCanvas* c = new TCanvas("c_summary", "kT Summary", 800, 800);
    c->SaveAs("kT_summary.pdf[");

    for (int ipt = 1; ipt < NpT; ++ipt) {
        for (int idr = 0; idr < NdeltaR; ++idr) {

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

            TH1D* h_nom = (TH1D*)kt_pythia[ipt][idr]->Clone("h_nom");
            TH1D* h_data = (TH1D*)kt_data[ipt][idr]->Clone("h_data");
            TH1D* h_data_band = (TH1D*)kt_data[ipt][idr]->Clone("h_data_band");
            TH1D* h_powpy = (TH1D*)kt_powpy[ipt][idr]->Clone("h_powpy");
            TH1D* h_sherpa = (TH1D*)kt_sherpa[ipt][idr]->Clone("h_sherpa");
            TH1D* h_herwig = (TH1D*)kt_herwig[ipt][idr]->Clone("h_herwig");
            h_nom->SetDirectory(0);
            h_data->SetDirectory(0);
            h_powpy->SetDirectory(0);
            h_sherpa->SetDirectory(0);
            h_herwig->SetDirectory(0);

            // h_nom->Scale(1./h_nom->Integral());
            // h_data->Scale(1./h_data->Integral());
            // h_powpy->Scale(1./h_powpy->Integral());
            // h_sherpa->Scale(1./h_sherpa->Integral());
            // h_herwig->Scale(1./h_herwig->Integral());

            


            h_nom->GetYaxis()->SetLabelSize(0.08);
            h_nom->GetYaxis()->SetTitleSize(0.05);
            h_nom->GetYaxis()->SetTitleSize(0.05);
            

            h_data->SetFillColor(TColor::GetColor("#cccccc")); 
            h_data->SetFillStyle(1001); 

            h_nom->SetLineColor(TColor::GetColor("#5520c7"));
            h_data->SetLineColor(TColor::GetColor("#cccccc"));
            h_powpy->SetLineColor(TColor::GetColor("#473868"));
            h_sherpa->SetLineColor(TColor::GetColor("#FF6F61"));
            h_herwig->SetLineColor(TColor::GetColor("#015457"));

            h_nom->SetMarkerStyle(21);
            h_data->SetMarkerStyle(20);
            h_powpy->SetMarkerStyle(46);
            h_sherpa->SetMarkerStyle(22);
            h_herwig->SetMarkerStyle(34);

            h_nom->SetMarkerColor(TColor::GetColor("#5520c7"));
            h_data->SetMarkerColor(TColor::GetColor("#000000"));
            h_powpy->SetMarkerColor(TColor::GetColor("#473868"));
            h_sherpa->SetMarkerColor(TColor::GetColor("#FF6F61"));
            h_herwig->SetMarkerColor(TColor::GetColor("#015457"));
            
            h_nom->SetLineStyle(1);
            h_data->SetLineStyle(1);
            h_powpy->SetLineStyle(1);
            h_sherpa->SetLineStyle(1);
            h_herwig->SetLineStyle(1);

            // TLegend* legUpper = new TLegend(0.8, 0.8, 0.88, 0.88);
            TLegend* legUpper = new TLegend(0.8, 0.6, 0.88, 0.9);
            legUpper->SetBorderSize(0);
            legUpper->SetFillStyle(0);
            legUpper->AddEntry(h_data, "Data", "lp");
            legUpper->AddEntry(h_nom, "Pythia", "p");
            legUpper->AddEntry(h_powpy, "Powheg", "p");
            legUpper->AddEntry(h_sherpa, "Sherpa", "p");
            legUpper->AddEntry(h_herwig, "Herwig", "p");
            legUpper->SetNColumns(1);
            legUpper->SetBorderSize(0);
            legUpper->SetLineColor(1);
            legUpper->SetLineStyle(1);
            legUpper->SetLineWidth(1);
            legUpper->SetFillColor(0);
            legUpper->SetFillStyle(0);
            legUpper->SetTextFont(42);
            legUpper->SetTextSize(0.034 / (1.0 - 2*ratioHeight));

            

            configureHist(h_data, kBlack);
            configurePlot(h_data, 1.0 - 2*ratioHeight);


            h_data->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#sigma}{dln(k_{T})}");
            h_data->Draw("p");
            h_nom->Draw("p same");
            h_powpy->Draw("p same");
            h_sherpa->Draw("p same");
            h_herwig->Draw("p same");

            legUpper->Draw();

            // add ATLAS labels 
            std::vector<std::string> labels;
            labels.push_back("#sqrt{s} = 13 TeV, 140 fb^{-1}");
            labels.push_back( Form("%.1f < ln(R/#DeltaR) < %.1f", dr_bins[idr], dr_bins[idr+1]));
            labels.push_back( Form("%.0f < p_{T} / GeV < %.0f", pt_bins[ipt], pt_bins[ipt+1]));
            makeLabels(labels, "Simulation",  1.0 - 2*ratioHeight, 0.2);

            // now do the second pad 
            ratioPad->cd();

            TH1D* h_ratio = (TH1D*)h_data->Clone("h_ratio");
            TH1D* h_ratio_pythia = (TH1D*)h_nom->Clone("h_ratio_pythia");
            TH1D* h_ratio_powpy = (TH1D*)h_powpy->Clone("h_ratio_powpy");
            TH1D* h_ratio_sherpa = (TH1D*)h_sherpa->Clone("h_ratio_sherpa");
            TH1D* h_ratio_herwig = (TH1D*)h_herwig->Clone("h_ratio_herwig");

            for (int b = 1; b < h_ratio->GetNbinsX()+1; b++){

                if (h_data->GetBinContent(b) != 0) {
                    h_ratio->SetBinContent(b, 1);
                    h_ratio_pythia->SetBinContent(b, h_nom->GetBinContent(b) / h_data->GetBinContent(b));
                    h_ratio_powpy->SetBinContent(b, h_powpy->GetBinContent(b) / h_data->GetBinContent(b));
                    h_ratio_sherpa->SetBinContent(b, h_sherpa->GetBinContent(b) / h_data->GetBinContent(b));
                    h_ratio_herwig->SetBinContent(b, h_herwig->GetBinContent(b) / h_data->GetBinContent(b));
                } 
            }

            // TH1D* h_band = (TH1D*)kt_total[ipt][idr]->Clone("h_band");
            TH1D* h_band = (TH1D*)h_data->Clone("h_band");
            // h_band->Divide(h_data);
            h_band->Reset(); 
            for (int b = 1; b < h_data_band->GetNbinsX()+1; b++) {
                double valData = h_data_band->GetBinContent(b);
                double errorData = h_data_band->GetBinError(b);
                h_band->SetBinContent(b, 1.0);
                // std::cout << " bin: " << b  <<  " value : " << valData << " error: " << errorData << std::endl;
                if (valData > 0) { 
                    h_band->SetBinError(b, errorData / valData);
                } else { 
                    h_band->SetBinError(b, 0.0);
                }
            }

            double ratioPadHeight = (2*ratioHeight - (ratioHeight + 0.1));
            configurePlot(h_band, ratioPadHeight);
            h_band->GetYaxis()->SetNdivisions(505);
            h_band->SetFillColor(TColor::GetColor("#cccccc"));
            h_band->SetMarkerSize(0);

            h_band->GetYaxis()->SetTitle("#splitline{Ratio to}{Data}");
            // h_band->GetYaxis()->SetTitle("Ratio to Pythia");
            h_band->GetYaxis()->SetRangeUser(0.5, 1.5);

            h_band->Draw("AXIS");
            // h_band->Draw("e2 ");
            // h_ratio->SetLineWidth(2);
            // h_ratio->Draw("hist");
            h_ratio_pythia->Draw("p same");
            // h_ratio_pythia->Draw(" same");
            h_ratio_powpy->Draw("p same");
            h_ratio_sherpa->Draw("p same");
            h_ratio_herwig->Draw("p same");

            TLine* line = new TLine(
                h_band->GetXaxis()->GetXmin(),1,
                h_band->GetXaxis()->GetXmax(),1
            );
            line->SetLineStyle(2);
            line->Draw();

            // now do the error breakdown pad
            errorPad->cd();

            auto makeSymmetric = [&](TH1D* h){
                TH1D* h_up = (TH1D*)h->Clone();
                TH1D* h_down = (TH1D*)h->Clone();

                for (int b = 1; b < h->GetNbinsX()+1; ++b) {
                    double val = h->GetBinContent(b);
                    h_up->SetBinContent(b, +val);
                    h_down->SetBinContent(b, -val);
                }

                return std::make_pair(h_up, h_down);
            };

            auto [tot_up, tot_down] = makeSymmetric(kt_total[ipt][idr]);
            auto [trk_up, trk_down] = makeSymmetric(kt_tracking[ipt][idr]);
            auto [purw_up, purw_down] = makeSymmetric(kt_purw[ipt][idr]);
            auto [jes_up, jes_down] = makeSymmetric(kt_JES[ipt][idr]);
            auto [modeling_up, modeling_down] = makeSymmetric(kt_Modeling[ipt][idr]); 


            double errorPadHeight = ratioHeight + 0.1;
            tot_up->GetYaxis()->SetTitle("#splitline{Relative}{Uncertainty}");
            tot_up->GetYaxis()->SetRangeUser(-0.2,0.2);
            tot_up->GetXaxis()->SetTitle("ln(k_{T}/GeV)");
            tot_up->SetFillColorAlpha(kGray+1, 0.4);

            
            tot_up->SetLineColor(kBlack);
            tot_down->SetLineColor(kBlack);
            configurePlot(tot_up, errorPadHeight);
            tot_up->Draw("axis");
            TH1D* h_band_err = (TH1D*)kt_total[ipt][idr]->Clone("h_band_err");
            double maxUnc = 0;

            for (int b = 1; b < h_band_err->GetNbinsX()+1; ++b) {
                double unc = kt_total[ipt][idr]->GetBinContent(b);
                h_band_err->SetBinContent(b, 0.0);
                h_band_err->SetBinError(b, unc);
                maxUnc = std::max(maxUnc, unc);
            }

            h_band_err->SetFillColorAlpha(kGray+1, 0.4);
            h_band_err->SetLineColor(kBlack);

            h_band_err->GetYaxis()->SetTitle("#splitline{Relative}{Uncertainty}");
            h_band_err->GetYaxis()->SetRangeUser(-1.2*maxUnc, 1.2*maxUnc);
            h_band_err->SetMarkerSize(0); 
            h_band_err->SetMarkerStyle(0); 

            configurePlot(h_band_err, errorPadHeight);
            h_band_err->Draw("HIST E2");

            trk_up->SetLineColor(TColor::GetColor("#d62728"));
            trk_down->SetLineColor(TColor::GetColor("#d62728"));
            trk_up->Draw("hist same");
            trk_down->Draw("hist same");

            purw_up->SetLineColor( TColor::GetColor("#8B4000"));
            purw_down->SetLineColor( TColor::GetColor("#8B4000"));
            purw_up->Draw("hist same");
            purw_down->Draw("hist same");

            modeling_up->SetLineColor(TColor::GetColor("#FFA62B"));
            modeling_down->SetLineColor(TColor::GetColor("#FFA62B"));
            modeling_up->Draw("hist same");
            modeling_down->Draw("hist same");

            jes_up->SetLineColor(TColor::GetColor("#1f77b4"));
            jes_down->SetLineColor(TColor::GetColor("#1f77b4"));
            jes_up->Draw("hist same");
            jes_down->Draw("hist same");

            // tile_up->SetLineColor(TColor::GetColor("#188535"));
            // tile_down->SetLineColor(TColor::GetColor("#188535"));
            // tile_up->Draw("hist same");
            // tile_down->Draw("hist same");

            std::vector<TH1*> systUp = {tot_up, trk_up, purw_up, jes_up, modeling_up // tile_up
            };
            std::vector<std::string> errNames = {"Total Syst.", "Tracking", "PURW", "JES", "MODELING", 
                // "TILE"
            };
            TLegend*legError =  makeLegends(systUp, errNames, 0.4, 0.18, 0.92, 0.6, 4);
            legError->Draw();


            c1->SaveAs("kT_summary.pdf");
            delete c1;
        }
    }

    c->SaveAs("kT_summary.pdf]");
}

int main(int argc, char* argv[]){

    SetAtlasStyle();
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    gStyle->SetTitleSize(0.04, "XY");
    gStyle->SetLabelSize(0.04, "XY");

    gStyle->SetEndErrorSize(0); 

    if (argc != 2) {
        std::cout << "Use this code like this: ./plotting <error_file.txt>" << std::endl;
        return 0;
    }

    std::string errFileName(argv[1]);

    std::string basePath = "/data/ldelagra/LundPlane/Unfolding/syst_PrimaryLJP_2026/";
    std::string basePath_data = "/data/ldelagra/LundPlane/Unfolding/nominal_PrimaryLJP/Final_Files_16-03-2026/";
    // std::string basePath_data = "/data/ldelagra/LundPlane/Unfolding/nominal_PrimaryLJP/";

    std::string configFile = errFileName;

    std::ifstream infile(configFile);
    std::string line;

    std::vector<SystEntry> entries;


    while (std::getline(infile, line)) {

        if (line.empty() || line[0] == '#') continue;

        std::stringstream ss(line);
        SystEntry e;

        ss >> e.group >> e.up >> e.nominal >> e.down >> e.groupID >> e.type >> e.colorHex;
        entries.push_back(e);
    }

    // Load nominal once 
    TH1D* h_nom_global = loadHist(basePath, "nominal");
    TH1D* h_data_global = loadHist(basePath_data, "");
    TH1D* h_powpy_global = loadHist(basePath, "POWPY"); // This should be POWPY
    TH1D* h_sherpa_global = loadHist(basePath, "SHERPA");
    TH1D* h_herwig_global = loadHist(basePath, "HERWIG");

    TCanvas *coriginal = new TCanvas("", "", 10000, 2000); 
    TPad* upperPad = new TPad("pad1","",0.0,0.3,1.0,1.0);
    TPad* ratioPad = new TPad("pad2","",0.0,0,1.0,0.3);
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
    upperPad->SetLogy();
    h_nom_global->SetDirectory(0);
    h_data_global->SetDirectory(0);
    h_powpy_global->SetDirectory(0);
    h_sherpa_global->SetDirectory(0);
    h_herwig_global->SetDirectory(0);
    
    h_data_global->SetFillColor(TColor::GetColor("#cccccc")); 
    h_data_global->SetFillStyle(1001); 

    h_nom_global->SetLineColor(TColor::GetColor("#5520c7"));
    h_data_global->SetLineColor(TColor::GetColor("#cccccc"));
    h_powpy_global->SetLineColor(TColor::GetColor("#473868"));
    h_sherpa_global->SetLineColor(TColor::GetColor("#FF6F61"));
    h_herwig_global->SetLineColor(TColor::GetColor("#015457"));

    h_nom_global->SetMarkerStyle(21);
    h_data_global->SetMarkerStyle(20);
    h_powpy_global->SetMarkerStyle(46);
    h_sherpa_global->SetMarkerStyle(22);
    h_herwig_global->SetMarkerStyle(34);

    h_nom_global->SetMarkerColor(TColor::GetColor("#5520c7"));
    h_data_global->SetMarkerColor(TColor::GetColor("#000000"));
    h_powpy_global->SetMarkerColor(TColor::GetColor("#473868"));
    h_sherpa_global->SetMarkerColor(TColor::GetColor("#FF6F61"));
    h_herwig_global->SetMarkerColor(TColor::GetColor("#015457"));

    h_nom_global->SetLineStyle(1);
    h_data_global->SetLineStyle(1);
    h_powpy_global->SetLineStyle(1);
    h_sherpa_global->SetLineStyle(1);
    h_herwig_global->SetLineStyle(1);

    TLegend* legUpper = new TLegend(0.8, 0.6, 0.88, 0.9);
    legUpper->SetBorderSize(0);
    legUpper->SetFillStyle(0);
    legUpper->AddEntry(h_data_global, "Data", "lp");
    legUpper->AddEntry(h_nom_global, "Pythia", "p");
    legUpper->AddEntry(h_powpy_global, "Powheg", "p");
    legUpper->AddEntry(h_sherpa_global, "Sherpa", "p");
    legUpper->AddEntry(h_herwig_global, "Herwig", "p");
    legUpper->SetNColumns(1);
    legUpper->SetBorderSize(0);
    legUpper->SetLineColor(1);
    legUpper->SetLineStyle(1);
    legUpper->SetLineWidth(1);
    legUpper->SetFillColor(0);
    legUpper->SetFillStyle(0);
    legUpper->SetTextFont(42);
    // legUpper->SetTextSize(0.034 / (1.0 - 2*ratioHeight));
    h_nom_global->GetXaxis()->SetRangeUser(150, 300);
    h_nom_global->Draw("P");
    h_data_global->Draw("P SAME");
    h_powpy_global->Draw("P SAME");
    h_sherpa_global->Draw("P SAME");
    h_herwig_global->Draw("P SAME"); 
    legUpper->Draw("SAME");

    ratioPad->cd(); 
    
    TH1D* h_ratio = (TH1D*)h_data_global->Clone("h_ratio");
    TH1D* h_ratio_pythia = (TH1D*)h_nom_global->Clone("h_ratio_pythia");
    TH1D* h_ratio_powpy = (TH1D*)h_powpy_global->Clone("h_ratio_powpy");
    TH1D* h_ratio_sherpa = (TH1D*)h_sherpa_global->Clone("h_ratio_sherpa");
    TH1D* h_ratio_herwig = (TH1D*)h_herwig_global->Clone("h_ratio_herwig");

    for (int b = 1; b < h_ratio->GetNbinsX()+1; b++){

        if (h_data_global->GetBinContent(b) != 0) {
            h_ratio->SetBinContent(b, 1);
            h_ratio_pythia->SetBinContent(b, h_nom_global->GetBinContent(b) / h_data_global->GetBinContent(b));
            h_ratio_powpy->SetBinContent(b, h_powpy_global->GetBinContent(b) / h_data_global->GetBinContent(b));
            h_ratio_sherpa->SetBinContent(b, h_sherpa_global->GetBinContent(b) / h_data_global->GetBinContent(b));
            h_ratio_herwig->SetBinContent(b, h_herwig_global->GetBinContent(b) / h_data_global->GetBinContent(b));
        } 
    }
    h_ratio->SetLineWidth(2);
    h_ratio_pythia->GetYaxis()->SetRangeUser(-10, 10);
    h_ratio->Draw("hist");
        h_ratio_pythia->GetXaxis()->SetRangeUser(150, 300);

    h_ratio_pythia->Draw("p");
    // h_ratio_pythia->Draw(" same");
    h_ratio_powpy->Draw("HIST same");
    h_ratio_sherpa->Draw("HIST same");
    h_ratio_herwig->Draw("HIST same");
    coriginal -> SaveAs("./kt_summary_all_pTbin2.pdf");
    std::map<std::string, TH1D*> groupAccum;
    std::map<std::string, int> groupColor;

    for (auto& e : entries) {

        TH1D* h_nom = loadHist(basePath, e.nominal);
        // std::cout << "au moins je suis ici .. "<< std::endl;
        TH1D* h_up  = loadHist(basePath, e.up);
        TH1D* h_down = loadHist(basePath, e.down);

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
    TH1D* h_total = (TH1D*)h_nom_global->Clone("h_total");
    h_total->Reset();

    for (int i = 1; i <= h_total->GetNbinsX(); ++i) {

        double sum2 = 0;

        for (auto& [name, hist] : groupAccum) {
            double v = hist->GetBinContent(i);
            sum2 += v * v;
        }

        h_total->SetBinContent(i, std::sqrt(sum2));
    }

    // Plot 1D ditirbutions for kT in different pT and deltaR bins
    std::cout << " Data NbinsX : " << h_data_global->GetNbinsX() << std::endl;
    std::cout << " Pythia NbinsX : " << h_nom_global->GetNbinsX() << std::endl;
    std::cout << " Sherpa NbinsX : " << h_sherpa_global->GetNbinsX() << std::endl;
    std::cout << " \\\\\\\\\\\\\\\\\\\\\\\\\\\\  " << std::endl;
    std::cout << " Data integral : " << h_data_global->Integral() << std::endl;
    std::cout << " Pythia integral : " <<h_nom_global->Integral() << std::endl;
    std::cout << " Sherpa integral : " <<h_sherpa_global->Integral() << std::endl;

    auto kt_pythia = extract_kT_distributions(h_nom_global, "nominal");
    auto kt_data   = extract_kT_distributions(h_data_global, "DATA");
    auto kt_powpy  = extract_kT_distributions(h_powpy_global, "POWPY");
    auto kt_sherpa = extract_kT_distributions(h_sherpa_global, "SHERPA");
    auto kt_herwig = extract_kT_distributions(h_herwig_global, "HERWIG");

    // auto kt_tile = extract_kT_distributions(groupAccum["TILE"], "tile");
    auto kt_purw     = extract_kT_distributions(groupAccum["PURW"], "purw");
    auto kt_tracking = extract_kT_distributions(groupAccum["TRACKING"], "tracking");
    auto kt_JES      = extract_kT_distributions(groupAccum["JES"], "JES");
    auto kt_Modeling = extract_kT_distributions(groupAccum["MODELING"], "modeling");
    auto kt_total    = extract_kT_distributions(h_total, "total");

    DrawTH1DPlotsWithError(kt_pythia, kt_data, kt_powpy, kt_sherpa, kt_herwig, kt_total, kt_tracking, kt_purw, kt_JES, kt_Modeling);
    // DrawTH1DPlots(h_total, groupAccum["TRACKING"], groupAccum["PURW"], groupAccum["JES"], groupAccum["MODELING"]);
    // TFile* fout = new TFile("syst_from_config.root", "RECREATE");

    // h_total->Write();
    // for (auto& [name, hist] : groupAccum) hist->Write();

    // fout->Close();

    std::cout << "Done." << std::endl;
}