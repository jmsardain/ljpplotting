#ifndef PlottingTools_h
#define PlottingTools_h

#include <iostream>
#include <vector>
#include <cmath>

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TFile.h"
#include "TString.h"
#include "TError.h" // To set gErrorIgnoreLevel

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TF2.h"

#include "Utils.h"
#include "DrawingUtils.h"
#include "AtlasStyle.C"


//**************** Plotting script ****************************//

TString name2Legend(TString name) {
	TString S = name;
	S.ToLower();
	if (S.Contains("udata")) { return "Unfolded data"; }
	if (S.Contains("data"))  { return "Data"; }
	return name;
}


int name2color(TString name) {
	TString S = name;
	S.ToLower();

    if (S.Contains("1 rwiter")) { return kBlue; }
	if (S.Contains("2 rwiter")) { return kRed; }
    if (S.Contains("3 rwiter")) { return kGreen+2; }
	if (S.Contains("4 rwiter"))  { return kMagenta; }
	if (S.Contains("5 rwiter"))  { return kOrange; }

    if (S.Contains("1rwiter")) { return kBlue; }
	if (S.Contains("2rwiter")) { return kRed; }
    if (S.Contains("3rwiter")) { return kGreen+2; }
	if (S.Contains("4rwiter"))  { return kMagenta; }
	if (S.Contains("5rwiter"))  { return kOrange; }

    if (S.Contains("udata")) { return kOrange; }
    if (S.Contains("data")) { return kGreen+2; }

    if (S.Contains("statunc")) { return kBlack; }
    if (S.Contains("_1down")) { return kRed; }
    if (S.Contains("_1up")) { return kRed; }
    if (S.Contains("nominal")) { return kBlue; }

    if (S.Contains("shower"))  { return kTeal+3; }

    if (S.Contains("10iter"))  { return kOrange+1; }
	if (S.Contains("0iter")) { return kAzure+1; }
    if (S.Contains("1iter")) { return kBlue; }
	if (S.Contains("2iter")) { return kRed; }
    if (S.Contains("3iter")) { return kGreen+2; }
	if (S.Contains("4iter"))  { return kMagenta; }
	if (S.Contains("5iter"))  { return kOrange; }
	if (S.Contains("6iter"))  { return kCyan+2; }
    if (S.Contains("7iter"))  { return kPink+1; }
    if (S.Contains("8iter"))  { return kViolet+1; }
    if (S.Contains("9iter"))  { return kGreen+1; }


    if (S.Contains("initial")) { return kBlack; }
    if (S.Contains("final")) { return kRed; }
    
    if (S.Contains("biasids")) { return kAzure+1; }

	return kBlack;
}

void DrawHistos(std::vector<TH1D*> histos, TString outputPath, TString FileName, std::vector<std::pair<double, double>> x_ranges, int iter, double ymin, double ymax, bool logAxis=false, bool sum=false) {
    
	if (histos.empty()) {
        std::cerr << "Error: No histograms provided!" << std::endl;
        return;
    }
    
    TCanvas* can = new TCanvas("can", "", 0., 0., 2400, 1350);
    can->Divide(4, 3);

    Float_t cleft = 0.12, cright = 0.04, cbottom = 0.15, ctop = 0.04;
    can->cd();
    can->SetMargin(cleft, cright, cbottom, ctop);
    
    TLegend* leg = new TLegend(0.2, 0.7, 0.4, 0.9);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetTextSize(0.03);
    
    for (TH1D* histo : histos) {
        histo->SetLineColorAlpha(name2color(histo->GetName()), 0.5);
        histo->SetMarkerColor(name2color(histo->GetName()));
        histo->SetMarkerSize(0.5);
        histo->SetMarkerStyle(20);
        leg->AddEntry(histo, histo->GetName(), "pl");
    }
    
    TLatex *tex = new TLatex;
    tex->SetNDC();
    tex->SetTextSize(0.03);
    tex->SetTextAlign(11);
    tex->SetTextColor(kBlack);
    
    int cw = 1;
    for (size_t i = 0; i < x_ranges.size(); ++i) {
        can->cd(cw);
        double x1 = x_ranges[i].first;
        double x2 = x_ranges[i].second;
        
        if (logAxis) gPad->SetLogy();
        
        for (size_t j = 0; j < histos.size(); ++j) {
            TH1D* h_clone = (TH1D*)histos[j]->Clone();
            h_clone->SetLineWidth(1);
            h_clone->GetXaxis()->SetRangeUser(x1, x2);
            if (ymin != -99. && ymax != -99.) { h_clone->GetYaxis()->SetRangeUser(ymin, ymax); }
            h_clone->Draw(j == 0 ? "oe hist" : "oe hist same");
        }
        
        leg->Draw();
        double ypos = 0.95;
        DrawATLASLabel(0.55, ypos -= 0.05);
        tex->DrawLatex(0.55, ypos -= 0.05, "#sqrt{s} = 13 TeV, 140 fb^{-1}");
        tex->DrawLatex(0.55, ypos -= 0.05, "Anti-#it{k}_{t} #it{R} = 0.4 (PFlow+JES)");
        if (iter!=999) { tex->DrawLatex(0.67, ypos -= 0.05, Form("nIter=%i", iter)); }
        
        cw++;
    }
    
    TString outputFile = logAxis ? outputPath + Form("%s_%diter_log", FileName.Data(), iter) : outputPath + Form("%s_%diter", FileName.Data(), iter);
    std::cout << "Generating file: " << outputFile << std::endl;
    can->SaveAs(outputFile + ".pdf");
    can->SaveAs(outputFile + ".root");

	TCanvas* can2  = new TCanvas("can2", "", 0., 0., 800, 600);

	if (logAxis) gPad->SetLogy(); 

	for (size_t j = 0; j < histos.size(); ++j) {
        TH1D* h_clone = (TH1D*)histos[j]->Clone();
        h_clone->SetLineWidth(2);
        h_clone->GetXaxis()->SetRangeUser(x_ranges[0].first, x_ranges[0].second);
        h_clone->GetXaxis()->SetTitle("k_{T}-DeltaR-p_{T}");
        if (ymin != -99. && ymax != -99.) { h_clone->GetYaxis()->SetRangeUser(ymin, ymax); }
        h_clone->Draw(j == 0 ? "oe hist" : "oe hist same");
    }
    
    leg->Draw();
	double ypos = 0.95;
	DrawATLASLabel(0.55, ypos -= 0.05);
	tex->DrawLatex(0.55, ypos -= 0.05, "#sqrt{s} = 13 TeV, 140 fb^{-1}");
	tex->DrawLatex(0.55, ypos -= 0.05, "Anti-#it{k}_{t} #it{R} = 0.4 (PFlow+JES)");
	if (iter!=999) { tex->DrawLatex(0.67, ypos -= 0.05, Form("nIter=%i", iter)); }

	TString outputFile2 = logAxis ? outputPath + Form("%s_%diter_log_TOT", FileName.Data(), iter) : outputPath + Form("%s_%diter_TOT", FileName.Data(), iter);
    std::cout << "Generating file: " << outputFile2 << std::endl;
    can2->SaveAs(outputFile2 + ".pdf");
    can2->SaveAs(outputFile2 + ".root");
	
    delete can;
	delete can2;
}

void DrawHistosSeparate(std::vector<TH1D*> histos, TString outputPath, TString FileName, std::vector<std::pair<double, double>> x_ranges, int iter, double ymin, double ymax, bool logAxis=false, bool sum=false) {
    
	if (histos.empty()) {
        std::cerr << "Error: No histograms provided!" << std::endl;
        return;
    }

    TString outputFile = logAxis ? outputPath + Form("%s_%diter_log", FileName.Data(), iter) : outputPath + Form("%s_%diter.pdf", FileName.Data(), iter);
    std::cout << "Generating file: " << outputFile << std::endl;
    
    TCanvas* can = new TCanvas("can", "", 0., 0., 800, 600);

    //Float_t cleft = 0.12, cright = 0.04, cbottom = 0.15, ctop = 0.04;
    //can->SetMargin(cleft, cright, cbottom, ctop);
    can->Print(outputFile + "(");
    
    TLegend* leg = new TLegend(0.7, 0.7, 0.92, 0.9);
    //leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetTextSize(0.03);
    
    for (TH1D* histo : histos) {
        histo->SetLineColorAlpha(name2color(histo->GetName()), 0.5);
        histo->SetMarkerColor(name2color(histo->GetName()));
        histo->SetMarkerSize(0.5);
        histo->SetMarkerStyle(20);
        leg->AddEntry(histo, histo->GetName(), "p");
    }
    
    TLatex *tex = new TLatex;
    tex->SetNDC();
    tex->SetTextSize(0.05);
    tex->SetTextAlign(11);
    tex->SetTextColor(kBlack);
    
    for (size_t i = 0; i < x_ranges.size(); ++i) {
        double x1 = x_ranges[i].first;
        double x2 = x_ranges[i].second;
        
        if (logAxis) gPad->SetLogy();
        
        for (size_t j = 0; j < histos.size(); ++j) {
            TH1D* h_clone = (TH1D*)histos[j]->Clone();
            h_clone->SetLineWidth(2);
            //h_clone->SetLineWidth(0); // TEMPORARY
            h_clone->GetXaxis()->SetRangeUser(x1, x2);
            h_clone->GetXaxis()->SetTitle("k_{T}-DeltaR-p_{T}");
            //h_clone->GetYaxis()->SetTitle("data/MC reco");
            h_clone->GetYaxis()->SetTitle(FileName);
            if (ymin != -99. && ymax != -99.) { h_clone->GetYaxis()->SetRangeUser(ymin, ymax); }
            h_clone->Draw(j == 0 ? "o" : "o same");
        }
        
        leg->Draw();
        double ypos = 0.95;
        DrawATLASLabel(0.2, ypos -= 0.05);
        tex->DrawLatex(0.2, ypos -= 0.05, "#sqrt{s} = 13 TeV, 140 fb^{-1}");
        tex->DrawLatex(0.2, ypos -= 0.05, "Anti-#it{k}_{t} #it{R} = 0.4 (PFlow+JES)");
        if (iter!=999) { tex->DrawLatex(0.2, ypos -= 0.05, Form("nIter=%i", iter)); }
        
        can->Print(outputFile);
        can->Clear();
    }

    can->Print(outputFile + ")");
    //can->SaveAs(outputFile + ".pdf");
    //can->SaveAs(outputFile + ".root");

	TCanvas* can2  = new TCanvas("can2", "", 0., 0., 800, 600);

	if (logAxis) gPad->SetLogy(); 

	for (size_t j = 0; j < histos.size(); ++j) {
        TH1D* h_clone = (TH1D*)histos[j]->Clone();
        h_clone->SetLineWidth(2);
        h_clone->GetXaxis()->SetRangeUser(x_ranges[0].first, x_ranges[0].second);
        h_clone->GetXaxis()->SetTitle("k_{T}-DeltaR-p_{T}");
        h_clone->GetYaxis()->SetTitle(FileName);
        if (ymin != -99. && ymax != -99.) { h_clone->GetYaxis()->SetRangeUser(ymin, ymax); }
        h_clone->Draw(j == 0 ? "o" : "o same");
    }
    
    leg->Draw();
	double ypos = 0.95;
	DrawATLASLabel(0.2, ypos -= 0.05);
	tex->DrawLatex(0.2, ypos -= 0.05, "#sqrt{s} = 13 TeV, 140 fb^{-1}");
	tex->DrawLatex(0.2, ypos -= 0.05, "Anti-#it{k}_{t} #it{R} = 0.4 (PFlow+JES)");
	if (iter!=999) { tex->DrawLatex(0.2, ypos -= 0.05, Form("nIter=%i", iter)); }

	TString outputFile2 = logAxis ? outputPath + Form("%s_%diter_log_TOT", FileName.Data(), iter) : outputPath + Form("%s_%diter_TOT", FileName.Data(), iter);
    std::cout << "Generating file: " << outputFile2 << ".pdf" << std::endl;
    can2->SaveAs(outputFile2 + ".pdf");
    can2->SaveAs(outputFile2 + ".root");
	
    delete can;
	delete can2;
}

void DrawHistosWithSubPad(TH1D* h1, TH1D* h2, TH1D* hbottom, TString outputPath, TString FileName, std::pair<double, double> x_range, int niter, bool logAxis=false) {

	Float_t cleft = 0.12, cright = 0.04, cbottom = 0.3, ctop = 0.1;

	TCanvas *can = new TCanvas("", "", 0., 0., 800, 600);
    
    TPad *topPad = new TPad("topPad", "", 0.0, 0.3, 1.0, 1.0);
    topPad->SetMargin(cleft, cright, 0., ctop);
    topPad->Draw();

    TPad *bottomPad = new TPad("bottomPad", "", 0.0, 0.05, 1.0, 0.3);
    bottomPad->SetMargin(cleft, cright, cbottom, 0.0);
    bottomPad->Draw();
    
    topPad->cd();

    h1->SetLineColorAlpha(name2color(h1->GetName()), 0.5);
	h1->SetMarkerSize(0);
	h2->SetLineColorAlpha(name2color(h2->GetName()), 0.5);
	h2->SetMarkerSize(0);

    TLegend *legend = new TLegend(0.35, 0.76, 0.92, 0.83); // Lower left corner (x, y), upper right corner (x, y)
	legend->SetBorderSize(0);
	legend->SetFillColor(0);
	legend->SetTextSize(TEXTSIZE*0.7);
    legend->AddEntry(h1, h1->GetName(), "l");
    legend->AddEntry(h2, h2->GetName(), "l");

    hbottom->SetLineColor(kGray);
	hbottom->SetMarkerColor(kBlack);
	hbottom->SetMarkerSize(0.5);
	hbottom->GetXaxis()->SetTitle("k_{T}-DeltaR-p_{T}");
	hbottom->GetYaxis()->SetTitle(hbottom->GetName());
	hbottom->GetXaxis()->SetLabelColor(1);
    hbottom->GetXaxis()->SetTitleSize(0.15);
    hbottom->GetYaxis()->SetTitleSize(0.15);
    hbottom->GetXaxis()->SetLabelSize(0.1);
    hbottom->GetYaxis()->SetLabelSize(0.1);
    hbottom->GetXaxis()->SetTitleOffset(1.0);
    hbottom->GetYaxis()->SetTitleOffset(0.3);

	TLatex *tex = new TLatex;
	tex->SetNDC();
	tex->SetTextSize(TEXTSIZE);
	tex->SetTextAlign(11);
	tex->SetTextColor(kBlack);

    topPad->cd();

	if (logAxis) gPad->SetLogy();

    h1->GetXaxis()->SetRangeUser(x_range.first, x_range.second);
    h2->GetXaxis()->SetRangeUser(x_range.first, x_range.second);

	h1->Draw("histe");
	h2->Draw("samehiste");

	legend->Draw();

	double ypos = 0.85;
	DrawATLASLabel(0.67, ypos -= 0.05);
	tex->DrawLatex(0.67, ypos -= 0.05, "#sqrt{s} = 13 TeV, 140 fb^{-1}");
	tex->DrawLatex(0.67, ypos -= 0.05, "Anti-#it{k}_{t} #it{R} = 0.4 (PFlow+JES)");
	tex->DrawLatex(0.67, ypos -= 0.05, Form("nIter=%i", niter));

	// Subpad
	bottomPad->cd();
	hbottom->GetXaxis()->SetRangeUser(x_range.first, x_range.second);
	hbottom->Draw("pesame");

	//TLine* line = new TLine(x_range.first, 0.0, x_range.second, 0.0);
    //line->SetLineStyle(1);
	//line->SetLineColorAlpha(kBlue, 0.5);
    //line->Draw();

	// Save figure

    TString outputFileName = logAxis ? outputPath + Form("%s_%diter_log", FileName.Data(), niter) : outputPath + Form("%s_%diter", FileName.Data(), niter);
	std::cout << "Generating file: " << outputFileName << std::endl;
	can->SaveAs(outputFileName + ".pdf");
    can->SaveAs(outputFileName + ".root");

	delete tex;
	delete legend;
	delete topPad;
    delete bottomPad;
	delete can;

}


void TransformToAbsolute(TH1D* h) {
    if (!h) return;

    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        double content = h->GetBinContent(i);
        h->SetBinContent(i, fabs(content));
    }
}

TH1D* GetStatUnc(TString statFileName) {

	TH1D* h_nom = GetHistoCopy<TH1D>(statFileName, "UData_TotStat");

    TH1D* hStatunc = new TH1D("StatUnc", "StatUnc", h_nom->GetNbinsX(), 
                              h_nom->GetXaxis()->GetXmin(), h_nom->GetXaxis()->GetXmax());
							  
    // Loop over bins and set the content and statistical uncertainty
    for (int i = 1; i <= h_nom->GetNbinsX(); ++i) {
        double error = h_nom->GetBinError(i);
        double content = h_nom->GetBinContent(i);
		
		hStatunc->SetBinContent(i, (content != 0) ? TMath::Abs(error / content) : 0);
		hStatunc->SetBinError(i, 0.);
    }

	return hStatunc;
}

TH1D* GetRatioHist(TH1D* h1, TH1D* h2) {

	TH1D *ratio = new TH1D("ratio", "ratio", h2->GetNbinsX(), h2->GetXaxis()->GetXmin(), h2->GetXaxis()->GetXmax());

	for (int i = 1; i <= h2->GetNbinsX(); ++i) {
	    if (h1->GetBinContent(i) != 0. && h2->GetBinContent(i) != 0.) {
	        float b1 = h1->GetBinContent(i);
			float eb1 = h1->GetBinError(i);
			float b2 = h2->GetBinContent(i);
			float eb2 = h2->GetBinError(i);

	        ratio->SetBinContent(i, b2/b1);
	        ratio->SetBinError(i, TMath::Power( TMath::Power(eb1/b1, 2) + TMath::Power(eb2/b2, 2), 0.5) );
	    } else { 
			ratio->SetBinContent(i, 0.);
	        ratio->SetBinError(i, 0.);
		}
	}

    return ratio;
}

void DrawHistosWithRatioSubPad(std::vector<TH1D*> histos,
                              std::vector<TH1D*> histoRatioAll,
                              TString outputPath,
                              TString plotNameBase,
                              std::vector<std::pair<double, double>> x_ranges,
                              std::vector<int> iterations,
                              double ymin, double ymax,
                              bool logY = true,
                              TString yAxisLabel = "",
                              TString ratioYAxisLabel = "Iter / Iter0")
{


    // Check input
    if(histos.size() != iterations.size()){
        printf("[DrawHistosWithRatioSubPad] Error: histos.size() != iterations.size()\n");
        std::cout << histos.size() << std::endl;
        std::cout << histoRatioAll.size() << std::endl;
        std::cout << iterations.size() << std::endl;
        return;
    }
    // Define colors (same order as iterations)
    std::vector<int> colors = {kBlack, kRed+1, kBlue+1, kGreen+2, kMagenta+2, kOrange+7, kCyan+2, kPink+9, kViolet+9, kAzure+6};
    // Extend colors if needed
    while((int)colors.size() < (int)iterations.size()){
        colors.push_back(kGray + (colors.size() % 9));
    }

    for (size_t xidx=0; xidx<x_ranges.size(); ++xidx){
        auto [xmin, xmax] = x_ranges[xidx];
        TString canvasName = Form("%s_xrange%d_%d", plotNameBase.Data(), (int)xmin, (int)xmax);
        TCanvas* c = new TCanvas(canvasName, canvasName, 800, 800);
        c->Clear();

        // Define pads
        TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
        TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
        pad1->SetBottomMargin(0.02);
        pad2->SetTopMargin(0.03);
        pad2->SetBottomMargin(0.3);
        pad1->SetLogy(logY);
        //pad2->SetLogy(logY);
        pad1->Draw();
        pad2->Draw();

        // Prepare histograms for this x_range: clone and restrict range
        std::vector<TH1D*> h_clones;
        for (size_t i=0; i<histos.size(); ++i){
            TH1D* h = (TH1D*)histos[i]->Clone(Form("%s_clone_iter%d_x%d", histos[i]->GetName(), iterations[i], (int)xmin));
            // Restrict x-axis to xmin,xmax
            //int binmin = h->FindBin(xmin);
            //int binmax = h->FindBin(xmax);
            h->GetXaxis()->SetRange(xmin, xmax);
            h_clones.push_back(h);
        }

        // Draw upper pad
        pad1->cd();

        // Prepare frame histogram with axis limits and labels (use first histo as template)
        TH1D* hframe = (TH1D*)h_clones[0]->Clone("frame");
        hframe->SetTitle("");
        hframe->SetMinimum(ymin);
        hframe->SetMaximum(ymax);
        //hframe->GetXaxis()->SetTitle("p_{T} [GeV]");
        hframe->GetYaxis()->SetTitle(yAxisLabel);
        hframe->GetXaxis()->SetRangeUser(xmin, xmax);
        hframe->Draw("axis");

        // Draw histograms with lines
        for (size_t i=0; i<h_clones.size(); ++i){
            h_clones[i]->SetLineWidth(2);
            h_clones[i]->SetLineColor(colors[i]);
            if (i==4) {
                h_clones[i]->SetLineWidth(4);
            };
            h_clones[i]->SetMarkerColor(colors[i]);
            h_clones[i]->SetMarkerStyle(20+i); // different marker styles optionally
            h_clones[i]->Draw("hist same");
        }

        // Legend
        TLegend *leg = new TLegend(0.6,0.62,0.98,0.92);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.025);
        leg->SetNColumns(3);
        for (size_t i=0; i<iterations.size(); ++i){
            leg->AddEntry(h_clones[i], Form("Iter %d", iterations[i]), "l");
        }
        leg->Draw();

        // ATLAS labels (adjust coordinates as needed)
        TLatex *tex = new TLatex;
        tex->SetNDC();
        tex->SetTextSize(0.04);
        tex->SetTextAlign(11);
        tex->SetTextColor(kBlack);
	    double ypos = 0.95;
	    DrawATLASLabel(0.22, ypos -= 0.05);
	    tex->DrawLatex(0.22, ypos -= 0.05, "#sqrt{s} = 13 TeV, 140 fb^{-1}");
	    tex->DrawLatex(0.22, ypos -= 0.05, "Anti-#it{k}_{t} #it{R} = 0.4 (PFlow+JES)");

        // Draw lower pad (ratios)
        pad2->cd();

        // Draw axis frame for ratios first
        histoRatioAll[0]->Draw("axis");
        for (size_t i=0; i<histoRatioAll.size(); ++i){
            histoRatioAll[i]->SetLineColor(colors[i]);
            histoRatioAll[i]->SetLineWidth(2);
            if (i==4) {
                histoRatioAll[i]->SetLineWidth(4);
            };
            histoRatioAll[i]->SetMarkerColor(colors[i]);
            histoRatioAll[i]->GetXaxis()->SetRangeUser(xmin, xmax);
            histoRatioAll[i]->SetMarkerStyle(20+i); // different marker styles optionally
	        histoRatioAll[i]->GetXaxis()->SetTitle("k_{T}-DeltaR-p_{T}");
	        histoRatioAll[i]->GetYaxis()->SetTitle(ratioYAxisLabel);
	        histoRatioAll[i]->GetXaxis()->SetLabelColor(1);
            histoRatioAll[i]->GetXaxis()->SetTitleSize(0.15);
            histoRatioAll[i]->GetYaxis()->SetTitleSize(0.15);
            histoRatioAll[i]->GetXaxis()->SetLabelSize(0.1);
            histoRatioAll[i]->GetYaxis()->SetLabelSize(0.1);
            histoRatioAll[i]->GetXaxis()->SetTitleOffset(1.0);
            histoRatioAll[i]->GetYaxis()->SetTitleOffset(0.3);
            histoRatioAll[i]->Draw("hist same");
            if(i==0) histoRatioAll[i]->Draw("hist");
        }

        TLine* line = new TLine(xmin,1.,xmax,1.);
        line->SetLineStyle(2);
        line->Draw("same");

        c->cd();
        c->Update();

        // Save canvas
        TString outputFile = outputPath + "/" + canvasName;
        std::cout << "Generating file: " << outputFile << std::endl;
        c->SaveAs(outputFile + ".pdf");
        c->SaveAs(outputFile + ".root");

        // Cleanup clones and pads
        for(auto h : h_clones) delete h;
        delete leg;
        delete line;
        delete pad1;
        delete pad2;
        delete c;
    }
    
    for(auto h : histoRatioAll) delete h;
}

void Draw2DLund(TCanvas*can, TH2D* h, TString outputDir, int ptBinIndex, TString ptLabel, TString hLabel) {
    // TCanvas* can = new TCanvas("", "", 1200, 1200);
    can->cd();
    can->Clear();
    TPad* pad = new TPad("pad", "pad", 0.05, 0.05, 0.95, 0.95);
    pad->SetMargin(0.2, 0.2, 0.15, 0.28);
    pad->Draw();
    pad->cd();

    h->GetYaxis()->SetTitle("ln(k_{T})");
    h->GetYaxis()->SetTitleSize(0.05);
    // h->GetYaxis()->SetRangeUser(-3.,3.);
    h->GetXaxis()->SetTitle("ln(R/#DeltaR)");
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetTitleOffset(1.2);
    h->GetYaxis()->SetTitleOffset(1.2);
    gStyle->SetPaintTextFormat(".2f"); 
    h->Draw("COLZ0 TEXT");

    TLatex* tex = new TLatex;
    tex->SetNDC();
    tex->SetTextSize(0.04);
    tex->SetTextColor(kBlack);

    double ypos = 0.95;
    DrawATLASLabel(0.22, ypos -= 0.04);
    tex->DrawLatex(0.22, ypos -= 0.04, "#sqrt{s} = 13 TeV,  140 fb^{-1}");
    //tex->DrawLatex(0.22, ypos -= 0.04, "#sqrt{s} = 13 TeV,  Pythia, truth");
    tex->DrawLatex(0.22, ypos -= 0.04, "Anti-#it{k}_{t} #it{R} = 0.4");
    tex->DrawLatex(0.22, ypos -= 0.04, hLabel);
    tex->DrawLatex(0.22, ypos -= 0.04, ptLabel);

    std::cout << "Saving " << outputDir << ".pdf" << std::endl;

    can->SaveAs(outputDir + ".pdf");
    // can->SaveAs(outputDir + ".root");

    delete tex;
    // delete can;
    // return can;
}


#endif // #ifdef Utils_h
