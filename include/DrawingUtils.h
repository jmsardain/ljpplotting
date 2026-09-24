#include <TCanvas.h>
#include <TLatex.h>
#include <TMarker.h>
#include <TLine.h>
#include <TGraphErrors.h>
#include <TGraph2DErrors.h>
#include <TLegend.h>
#include <TPaletteAxis.h>

#include "AtlasLabels.C"

TCanvas *_can = new TCanvas;
TLatex *_tex = new TLatex;
double TEXTSIZE = 0.05;

std::vector<int> cols = {kBlack, kRed, kBlue, kViolet + 1, kGreen + 3, kOrange + 2, kCyan + 1, kMagenta + 2, kYellow + 2, kBlue - 2, kCyan + 3, kPink + 3, kMagenta + 2, kBlue - 10};
std::vector<int> mark = {20, 21, 22, 23, 26, 25, 24, 27, 28, 29, 30, 20, 2, 3, 5, 31, 32, 33, 34};
std::vector<int> mark2 = {24, 25, 26, 27, 28, 29, 30, 2, 3, 5, 31, 32, 33, 34, 20, 21, 22, 23, 24};

int _jetR;
TString _calibration, _methodDesc, _name;

void FormatHisto(TH1 *h, int mstyle, int col, double msize);
void FormatHisto(TH1 *h, int style, double msize = 0.8); 
void FormatGraph(TGraphErrors *h, int mstyle, int col, double msize);
void FormatGraph(TGraphErrors *h, int style, double msize = 0.8); 

void DrawText( TString txt, int col = kBlack, double x = 0.145, double y = 0.88);
void DrawTextR(TString txt, int col = kBlack, double x = 0.92,  double y = 0.8);
void DrawLabel(TString txt, int style, double x, double y);
void DrawLabel(TString txt, double x, double y, int mstyle, int col, double msize);


void FormatHisto(TH1 *h, int style, double msize) {
	if      (style == 1)  { FormatHisto(h, 20, kBlack,       msize + 0.5); }
	else if (style == 2)  { FormatHisto(h, 22, kRed,         msize + 0.5); }
	else if (style == 3)  { FormatHisto(h, 23, kBlue,        msize + 0.5); }
	else if (style == 4)  { FormatHisto(h, 24, kOrange + 7,  msize + 0.5); }
	else if (style == 5)  { FormatHisto(h, 25, kGreen + 2,   msize + 0.5); }
	else if (style == 6)  { FormatHisto(h, 26, kGray + 2,    msize + 0.5); }
	else if (style == 7)  { FormatHisto(h, 27, kBlue - 1,    msize + 0.5); }
	else if (style == 8)  { FormatHisto(h, 28, kCyan + 1,    msize + 0.5); }
	else if (style == 9)  { FormatHisto(h, 29, kMagenta + 2, msize + 0.5); }
	else if (style == 10) { FormatHisto(h, 30, kPink + 1,    msize + 0.5); }

	h->SetXTitle("");
	h->SetYTitle("#it{R}_{rel}");
	h->GetXaxis()->SetTitleOffset(0.95);
	h->GetYaxis()->SetTitleOffset(1.1);
	h->GetYaxis()->SetNdivisions(505);
}

void FormatHisto(TH1 *h, int mstyle, int col, double msize) {
	h->SetMarkerStyle(mstyle);
	h->SetMarkerSize(msize);
	h->SetMarkerColor(col);
	h->SetLineColor(col);
}

void FormatGraph(TGraphErrors *h, int style, double msize) {

	if      (style == 1)  { FormatGraph(h, 20, kBlack,      msize + 0.5); }
	else if (style == 2)  { FormatGraph(h, 23, kRed,        msize + 0.5); }
	else if (style == 3)  { FormatGraph(h, 22, kBlue,       msize + 0.5); }
	else if (style == 4)  { FormatGraph(h, 26, kOrange + 7, msize + 0.5); }
	else if (style == 5)  { FormatGraph(h, 33, kGreen + 2,  msize + 0.5); }
	else if (style == 6)  { FormatGraph(h, 28, kGray + 2,   msize + 0.5); }
	else if (style == 10) { FormatGraph(h, 24, kBlack,      msize + 0.5); }
	else if (style == 11) { FormatGraph(h, 24, kBlack,      msize + 0.5); }
	else if (style == 12) { FormatGraph(h, 25, kRed,        msize + 0.5); }
	else if (style == 13) { FormatGraph(h, 27, kBlue,       msize + 0.5); }

}

void FormatGraph(TGraphErrors *h, int mstyle, int col, double msize) {
	h->SetMarkerStyle(mstyle);
	h->SetMarkerSize(msize);
	h->SetMarkerColor(col);
	h->SetLineColor(col);
	h->SetLineWidth(2);
}

double SUBFIGURE_MARGIN = 0.06;
TPad* DrawSubPad(TCanvas* can, Float_t left = 0.12, Float_t right = 0.04, Float_t bottom = 0.125, double FIGURE2_RATIO = 0.36) {
	can->SetBottomMargin(FIGURE2_RATIO);
	// create new pad, fullsize to have equal font-sizes in both plots
	TPad *p = new TPad("p_test", "", 0, 0, 1, 1.0 - SUBFIGURE_MARGIN, 0, 0, 0); 
	p->SetFillStyle(0);
	p->SetMargin(left, right, bottom, 1.0 - FIGURE2_RATIO);
	p->Draw();
	return p;
}

void DrawATLASLabel(double x = 0.15, double y = 0.875) {
	ATLASLabel(x, y, (char*)"Work in progress", kBlack); //Internal
}

void PrintInfo() {
	DrawText(Form("Anti-#it{k}_{t} #it{R} = 0.%d ", _jetR) + _calibration + ", " + _methodDesc, kBlack, 0.12, 0.87);
	DrawTextR(_name, kBlack, 0.92, 0.97);
}

void PrintEtaBinInfo(double ptmin, double ptmax, TString Calibration = _calibration) {
	DrawText(Form("Anti-#it{k}_{t} #it{R} = 0.%d ", _jetR) + Calibration, kBlack, 0.92, 0.77);
	DrawTextR(Form("%.0f #leq #it{p}_{T}^{avg} < %.0f GeV", ptmin, ptmax), 1, 0.92, 0.72);
	//DrawTextR(_methodDesc, kBlack, 0.92, 0.875);
	DrawTextR(_name, kBlack, 0.92, 0.75);
}

void PrintPtBinInfo(double etamin, double etamax, TString Calibration = _calibration) {
	DrawText(Form("Anti-#it{k}_{t} #it{R} = 0.%d ", _jetR) + Calibration, kBlack, 0.92, 0.77);
	DrawTextR(Form("%.1f #leq #it{#eta}_{det} < %.1f", etamin, etamax), 1, 0.92, 0.72);
	//DrawTextR(_methodDesc, kBlack, 0.92, 0.875);
	DrawTextR(_name, kBlack, 0.92, 0.75);
}

void DrawText(TString txt, int col, double x, double y) {
	_tex->SetTextAlign(11); _tex->SetTextColor(col); _tex->DrawLatex(x, y, txt);
}

void DrawTextR(TString txt, int col, double x, double y) {
	_tex->SetTextAlign(31); _tex->SetTextColor(col); _tex->DrawLatex(x, y, txt);
	_tex->SetTextAlign(11); _tex->SetTextColor(kBlack);
}

void DrawGuideLineSyst(double min, double max) {
	static TLine *line = new TLine();
	line->SetLineStyle(2);
	line->SetLineWidth(1); line->SetLineColor(kRed);
	line->DrawLine(min, 0.05, max, 0.05);
	line->SetLineColor(kBlue);
	line->DrawLine(min, 0.02, max, 0.02);
}

void DrawGuideLines(double xmin, double xmax, double y1 = 0.97, double y2 = 1.03) {
	static TLine *line = new TLine();
	line->SetLineStyle(2);
	line->SetLineWidth(1);
	line->SetLineColor(kGray + 1); 
	line->SetLineColor(kBlack);
	line->DrawLine(-0.8, y1, -0.8, y2);
	line->DrawLine(0.8, y1, 0.8, y2);
	line->SetLineWidth(1); line->SetLineColor(kBlack);
	line->DrawLine(xmin, 1.0, xmax, 1.0);
	line->SetLineWidth(1);
	line->SetLineColor(kBlue);
	line->DrawLine(xmin, y1, xmax, y1);
	line->DrawLine(xmin, y2, xmax, y2);
}

void DrawCenterLine(double xmin, double xmax, double center = 1.0) {
	static TLine *line = new TLine();
	line->SetLineStyle(3); line->SetLineWidth(1);
	line->SetLineColor(kBlack);
	line->DrawLine(xmin, center, xmax, center);
}

void DrawLine(double xmin, double ymin, double xmax, double ymax) {
	static TLine *line = new TLine();
	line->SetLineStyle(3); line->SetLineWidth(1);
	line->SetLineColor(kBlack);
	line->DrawLine(xmin, ymin, xmax, ymax);
}

void DrawGuideLinesVsEta(double etaMin = -4.5, double etaMax = 4.5) {
	DrawGuideLines(etaMin, 0.97, etaMax, 1.03);
}

void DrawGuideLinesVsPt() {
	DrawGuideLines(20, 0.97, 2500, 1.03);
}

void DrawLabel(TString txt, double x, double y, int mstyle, int col, double msize) {
	TMarker *m = new TMarker(x, y, mstyle);
	m->SetNDC();
	m->SetMarkerSize(msize);
	m->SetMarkerColor(col);
	TLine *l = new TLine();
	l->SetLineWidth(2);
	l->SetLineColor(col);
	l->DrawLineNDC(x - 0.02, y, x - 0.005, y);
	l->DrawLineNDC(x + 0.005, y, x + 0.02, y);
	m->Draw();
	//_tex->SetTextSize(0.04);
	_tex->SetTextSize(TEXTSIZE);
	_tex->SetTextAlign(12);
	_tex->SetTextColor(kBlack);
	_tex->DrawLatex(x + 0.04, y, txt);
	_tex->SetTextSize(TEXTSIZE);
}

void DrawLabel(TString txt, int style, double x, double y) {
	double size = 1.2;
	if      (style == 1)  { DrawLabel(txt, x, y, 20, kBlack,      size + 0.5); }
	else if (style == 2)  { DrawLabel(txt, x, y, 22, kRed,        size + 0.5); }
	else if (style == 3)  { DrawLabel(txt, x, y, 23, kBlue,       size + 0.5); }
	else if (style == 4)  { DrawLabel(txt, x, y, 26, kOrange + 7, size + 0.5); }
	else if (style == 5)  { DrawLabel(txt, x, y, 33, kGreen + 2,  size + 0.5); }
	else if (style == 6)  { DrawLabel(txt, x, y, 28, kGray + 2,   size + 0.5); }
	else if (style == 10) { DrawLabel(txt, x, y, 24, kBlack,      size + 0.5); }
	else if (style == 11) { DrawLabel(txt, x, y, 24, kBlack,      size + 0.5); }
	else if (style == 12) { DrawLabel(txt, x, y, 25, kRed,        size + 0.5); }
	else if (style == 13) { DrawLabel(txt, x, y, 27, kBlue,       size + 0.5); }
}

void PrintIntegratedLumi(TString Lumi) {
	//DrawText(Form("#int #it{L dt} = %s fb^{-1}", Lumi.Data()), kBlack, 0.15, 0.8); //match public plots
	DrawText(Form("#sqrt{s} = 13 TeV, %s fb^{-1}", Lumi.Data()), kBlack, 0.145, 0.82); //match public plots
}

void PrintDataMC(TString dataDesc, TString mc1, TString mc2, bool drawSystMC = true) {
	//DrawLabel("Data " + dataDesc + ", #sqrt{s} = 13 TeV", 1, 0.38, 0.5); //match public plots
	DrawLabel("Data " + dataDesc, 1, 0.65, 0.87); //match public plots
	if (drawSystMC) {
		DrawLabel(mc1, 2, 0.65, 0.82);
		DrawLabel(mc2, 3, 0.65, 0.77);
	} else {
		DrawLabel(mc1, 2, 0.65, 0.82);
	}
}

void DrawPtEtaPoints(TGraph2DErrors *graph, VecD ptBins, VecD etaBins, double freezeCalibEta = 4.5, double minPt = 20.0,  double maxPt = 1500.0) {

	TLine *line = new TLine();
	line->SetLineColor(kGray);
	TGraphErrors *bin_xy = new TGraphErrors(), *avg_xy = new TGraphErrors();
	static TH2F *h = new TH2F("a", "", 1, 17.5, 2000, 1, etaBins[0], etaBins[etaBins.size() - 1]);
	h->GetXaxis()->SetMoreLogLabels();
	h->GetYaxis()->SetTitleOffset(0.8);
	h->SetXTitle("p_{T}^{avg}");
	h->SetYTitle("#eta_{det}");
	h->SetXTitle("p_{T}^{probe} [GeV]");
	h->Draw();

	for (unsigned int ieta = 0 ; ieta < etaBins.size() - 1 ; ++ieta) {
		double etal = etaBins[ieta];
		double eta = 0.5 * (etal + etaBins[ieta + 1]);
		line->DrawLine(25, etal, 1500, etal);
		if (ieta == etaBins.size() - 2) { line->DrawLine(minPt, etaBins[ieta + 1], 1500, etaBins[ieta + 1]); }
		for (unsigned int ipt = 0 ; ipt < ptBins.size() - 1 ; ++ipt) {
			double ptl = ptBins[ipt];
			double pt  = 0.5 * (ptl + ptBins[ipt + 1]);
			line->DrawLine(ptl, etaBins[0], ptl, etaBins[etaBins.size() - 1]);
			if (ipt == ptBins.size() - 2) { line->DrawLine(ptBins[ipt + 1], etaBins[0], ptBins[ipt + 1], 4.5); }
			int n = bin_xy->GetN();
			bin_xy->SetPoint(n, pt, eta);
		}
	}
	for (int i = 0 ; i < graph->GetN(); ++i) {
		avg_xy->SetPoint(i, graph->GetX()[i], graph->GetY()[i]);
		avg_xy->SetPointError(i, graph->GetEX()[i], graph->GetEY()[i]);
	}
	avg_xy->SetMarkerSize(0.4);
	avg_xy->Draw("P");

	float minEta = -freezeCalibEta;
	if (etaBins[0] >= 0) { minEta = 0; }

	line->SetLineColor(kRed + 1);
	line->SetLineWidth(2);

	line->DrawLine(minPt, minEta, minPt, freezeCalibEta);
	line->DrawLine(minPt, minEta, maxPt, minEta);
	line->DrawLine(maxPt, minEta, maxPt, freezeCalibEta);
	line->DrawLine(minPt, freezeCalibEta, maxPt, freezeCalibEta);
	h->Draw("sameaxis");
}

TString toLegend_method(TString method) {
	if (method == "MM") {
		return "matrix method";
	} else if (method == "SM") {
		return "central reference method";
	} else {
		return method;
	}
}

TString toLegend_jetAlgo(TString jetAlgo) {
	if (jetAlgo == "AntiKt4Truth") {
		return "truth";
	} else if (jetAlgo == "AntiKt4EMPFlow") {
		return "reco";
	} else {
		return jetAlgo;
	}
}

TString toLegend_generator(TString generator) {
	if (generator == "PowhegPythia8") {
		return "Powheg+Pythia8";
	} else if (generator == "PowhegHerwig7Angular") {
		return "Powheg+Herwig7";
	} else {
		return generator;
	}
}

TString toLegend_systName(TString systName) {
	if (systName == "j3") {
		return "#it{p}_{T}^{3} / #it{p}_{T}^{avg}";
	} else if (systName == "dphi") {
		return "#Delta#phi^{1,2}";
	} else {
		return systName;
	}
}
