/*
 *  Utils.h -- a collection of helper functions for
 *             HEP analysis
 */

#ifndef Utils_h
#define Utils_h

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <bitset>
#include <time.h>
#include <algorithm>
#include <sys/stat.h>

#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TChain.h"
#include "TFile.h"
#include "TEnv.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TROOT.h"
#include "TApplication.h"
#include "TError.h" // To set gErrorIgnoreLevel

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TF2.h"

#include <TGraphErrors.h>
#include <TGraph2DErrors.h>
#include "TProfile2D.h"
#include "TLorentzVector.h"

typedef std::vector<double> VecD;
typedef std::vector<float> VecF;
typedef std::vector<int> VecI;

///////////////////////////////////////////////////////////////////////////////////////

void error(TString msg) {
	std::cout << "ERROR:" << std::endl << std::endl << "  " << msg << std::endl << std::endl;
	abort();
}

//******************************************************************
//**************** vectorization methods/functions
//******************************************************************


inline TH1D* loadHist(const std::string& path, const std::string& tag, TString histName) {

    std::string fullPath = path + tag + "/4iter/Udata.root";
    if (histName == "UData_TotStat") { 
        fullPath = path + tag + "/4iter/StatUnc.root";
    } else if (histName == "biasIDS"){
        fullPath = path + "/4iter/Bias_results.root";
    }
    std::cout << "Loading histogram " << histName << " from: " << path << "/" << tag << std::endl;
    TFile* f = TFile::Open(fullPath.c_str());
    TH1D* h = (TH1D*)f->Get(histName);
    h->SetDirectory(0);
    f->Close();
    return (TH1D*)h->Clone();
}

std::vector<TString> Vectorize(TString str, TString sep = " ") {
	std::vector<TString> result; TObjArray *strings = str.Tokenize(sep.Data());
	if (strings->GetEntries() == 0) { return result; }
	TIter istr(strings);
	while (TSubString* os = (TSubString*)istr()) {
		if (os->String()[0] != '#') { result.push_back(os->String()); }
		else { break; }
	}
	return result;
}

VecD VectorizeD(TString str) {
	VecD result; std::vector<TString> vecS = Vectorize(str);
	for (unsigned int i = 0 ; i < vecS.size() ; ++i)
		result.push_back(atof(vecS[i]));
	return result;
}

VecD MakeUniformVector(int N, double min, double max) {
	VecD vec;
	double dx = (max - min) / N;
	for (int i = 0 ; i <= N ; ++i) { vec.push_back(min + i * dx); }
	return vec;
}

VecD MakeLogVector(int N, double min, double max) {
	VecD vec;
	double dx = (log(max) - log(min)) / N;
	for (int i = 0 ; i <= N ; ++i) { vec.push_back(exp(log(min) + i * dx)); }
	return vec;
}

//******************************************************************
//**************** functions to obtain histos ect.
//******************************************************************

TFile *Open(TString fn) { TFile *f = TFile::Open(fn); if (f == nullptr) { error("Can't open " + fn);} return f; }

TObject* Get(TFile *f, TString jetAlgo, TString objName) {
	f->cd(jetAlgo);
	TObject *obj = gDirectory->Get(objName);
	if (obj == nullptr) { error("Can't access " + objName + " in " + f->GetName()); }
	return obj;
}

template<typename THnD>
THnD* GetHisto(TFile *f, TString jetAlgo, TString hname) { return (THnD*)Get(f, jetAlgo, hname); }

template<typename THnD>
THnD* GetHisto(TFile *f, TString hn) { if (f->Get(hn) == nullptr) { error("Cannot access: " + hn + " in file " + f->GetName()); } return (THnD*)f->Get(hn); }

template<typename THnD>
THnD* GetHistoCopy(TString file_name, TString histo_name) {
	TFile *file = Open(file_name);
	THnD *histo = GetHisto<THnD>(file, histo_name);
	histo->SetDirectory(0); // Keep a copy of the histogram and allow for the file to be closed
	file->Close();
	return histo;
}

template<typename TFn>
TFn* GetFunctionCopy(TString file_name, TString func_name) {
	TFile *file = Open(file_name);
	TFn *func = GetHisto<TFn>(file, func_name);
	file->Close();
	return func;
}

//******************************************************************
//**************** functions to retrive calibration input points
//******************************************************************

TGraph2DErrors *GetCalibrationPoints(TGraph2DErrors *d, TGraph2DErrors *mcPrediction) {
	if (d->GetN() != mcPrediction->GetN()) {
		std::cout << "d->GetN() = " << d->GetN() << " and mcPrediction->GetN() = " << mcPrediction->GetN() << std::endl;
		error("The 'd' TGraph2D doesn't have the same number of points as the 'mcPrediction' TGraph2D.");
	}
	TGraph2DErrors *calib = new TGraph2DErrors();
	for (int i = 0 ; i < d->GetN() ; ++i) {
		double x = d->GetX()[i], y = d->GetY()[i];
		double ex = d->GetEX()[i], ey = d->GetEY()[i];
		double R_d = d->GetZ()[i], dR_d = d->GetEZ()[i];
		double R_mc = mcPrediction->GetZ()[i], dR_mc = mcPrediction->GetEZ()[i];
		double c = R_mc / R_d;
		if (R_d == 0) { std::cout << "Warning in GetCalibrationPoints: R_d == 0, err_c will be nan." << std::endl; }
		if (R_mc == 0) { std::cout << "Warning in GetCalibrationPoints: R_mc == 0, err_c will be nan." << std::endl; }
		double err_c = c * sqrt(pow(dR_d / R_d, 2) + pow(dR_mc / R_mc, 2));
		calib->SetPoint(i, x, y, c);
		calib->SetPointError(i, ex, ey, err_c);
	}
	return calib;
}

bool numberOfPointsEquality(TGraph2DErrors *g0, TGraph2DErrors *g, TString gn) {
	if (g != nullptr && g0->GetN() != g->GetN()) {
		std::cout << "Npoints g1 = " << g0->GetN() << ", g" << gn << " = " << g->GetN() << std::endl;
		return false;
	}
	return true;
}

bool greaterMaxErr(TGraph2DErrors *g, int i, double maxErr) {
	return g != nullptr && g->GetEZ()[i] > maxErr;
}

bool highPt(TGraph2DErrors *g, int i, double high_pT) {
	return g != nullptr && g->GetX()[i] > high_pT;
}

bool negativeZ(TGraph2DErrors *g, int i) {
	return g != nullptr && g->GetZ()[i] <= 0;
}

TGraph2DErrors* TrimGraph(TGraph2DErrors *g1, TGraph2DErrors *g2 = nullptr, TGraph2DErrors *g3 = nullptr, TGraph2DErrors *g4 = nullptr) {
	if (!numberOfPointsEquality(g1, g2, "2") || !numberOfPointsEquality(g1, g3, "3") || !numberOfPointsEquality(g1, g4, "4")) {
		error("Trimming issue!");
	}
	bool verbose = true;
	double high_pT = 100.0;
	double maxErr = 0.3;
	if (verbose) { std::cout << "TrimGraph: Maximum error for a point to be kept above " << high_pT << " GeV: " << maxErr << std::endl; }
	int nn = 0;
	TGraph2DErrors *g = new TGraph2DErrors();
	int zInf0 = 0, ezSupErr = 0;
	for (int i = 0 ; i < g1->GetN() ; ++i) {
		if (negativeZ(g1, i) || negativeZ(g2, i) || negativeZ(g3, i) || negativeZ(g4, i)) {
			++zInf0; continue;
		}
		bool high_pT = (highPt(g1, i, high_pT) || highPt(g2, i, high_pT) || highPt(g3, i, high_pT) || highPt(g4, i, high_pT));
		if (high_pT && (greaterMaxErr(g1, i, maxErr) || greaterMaxErr(g2, i, maxErr) || greaterMaxErr(g3, i, maxErr) || greaterMaxErr(g4, i, maxErr))) {
			++ezSupErr; continue;
		}
		g->SetPoint(nn, g1->GetX()[i], g1->GetY()[i], g1->GetZ()[i]);
		g->SetPointError(nn, g1->GetEX()[i], g1->GetEY()[i], g1->GetEZ()[i]);
		++nn;
	}
	if (verbose && (zInf0 > 0 || ezSupErr > 0)) {
		std::cout << "TrimGraph: Number of points for which z <= 0:                        " << zInf0 << std::endl;
		std::cout << "TrimGraph: Number of points for which pT > " << high_pT << " and error on z > " << maxErr << ": " << ezSupErr << std::endl;
	}
	return g;
}

TH2D* GetCalibration2DHist(TFile *f, TString jetAlgo, TString method, bool numInv = true, double pTmax = 0) {

	TH2D *avg_pt  = GetHisto<TH2D>(f, jetAlgo + "_AvgPtAvg_vs_EtaPtAvg_" + method);
	TH2D *avg_eta = GetHisto<TH2D>(f, jetAlgo + "_AvgEta_vs_EtaPtAvg_" + method);

	TH2D *calibGraph = (TH2D*)avg_pt->Clone("pTvsEtaDet_contentisR");

	int Npt = avg_pt->GetNbinsY();
	int Neta = avg_pt->GetNbinsX();

	for (int ipt = 1 ; ipt <= Npt ; ++ipt) {

		double ptmin = avg_pt->GetYaxis()->GetBinLowEdge(ipt);
		double ptmax = avg_pt->GetYaxis()->GetBinLowEdge(ipt + 1);

		if (pTmax > 0 && pTmax <= ptmin) { continue; }

		TString Rhname = jetAlgo + "/" + jetAlgo + "_" + Form("pt%.0fto%.0f_R_", ptmin, ptmax) + method;
		TH1D *hR = GetHisto<TH1D>(f, Rhname);

		// cut the Eta bin edges of SM distribution to same one of MM
		int ieta_min = 1, ieta_max = Neta;
		if (method == "SM") {
			TH1D *hR_MM;
			TString Rhname_MM = jetAlgo + "/" + jetAlgo + "_" + Form("pt%.0fto%.0f_R_", ptmin, ptmax) + "MM";
			hR_MM = GetHisto<TH1D>(f, Rhname_MM);
			//get the first non negative R
			while (hR_MM->GetBinContent(ieta_min) <= 0 && ieta_min < Neta) {
				++ieta_min;
			}
			//get the last non negative R
			while (hR_MM->GetBinContent(ieta_max) > 0 && ieta_max > 1) {
				--ieta_max;
			}
			delete hR_MM;
		}

		if (hR->GetNbinsX() != Neta) { error("Binning problem"); }
		for (int ieta = 1 ; ieta <= Neta ; ++ieta) {

			double dR      = hR->GetBinError(ieta);
			double pt      = avg_pt->GetBinContent(ieta, ipt);
			double err_pt  = avg_pt->GetBinError(ieta, ipt);
			double eta     = avg_eta->GetBinContent(ieta, ipt);

			TAxis* xaxis = calibGraph->GetXaxis(); TAxis* yaxis = calibGraph->GetYaxis();
			Int_t binx = xaxis->FindBin(eta);
			Int_t biny = yaxis->FindBin(pt);

			double R = hR->GetBinContent(ieta);

			// cut out of range bins
			if (ieta < ieta_min || ieta > ieta_max) { R = -99; }

			if (R < 0) {
				calibGraph->SetBinContent(binx, biny, 0);
				continue;
			}

			//double err_eta = avg_eta->GetBinError(ieta, ipt); // Unused variable
			//printf("Pt,Eta bins: (%2d,%2d), (%6.2f +/- %5.2f)%%\n", ipt, ieta, R * 100 - 100, dR * 100);
			//printf("%.2f < %.2f\n", avg_pt->GetXaxis()->GetBinLowEdge(ieta), eta);
			//printf("%.2f < %.2f\n", avg_pt->GetYaxis()->GetBinLowEdge(ipt), pt);
			if (numInv) {
				// pTavg = (pTprobe + pTref) / 2  and  R = pTprobe / pTref   => pTprobe = 2R / (R + 1) pTavg
				pt *= 2.0 * R / (R + 1); // pTavg -> pTprobe
				// from error propagation:
				err_pt = sqrt(pow(err_pt * 2.0 * R / (R + 1), 2) + pow(2.0 * pt / (R + 1) / (R + 1) * dR, 2));
			}
			//printf("( ptmin : ptmax :n : pt : eta : R )  = ( %.0f : %.0f : %d : %.0f : %.1f : %.2f )\n", ptmin, ptmax, n, pt, eta, R);

			calibGraph->SetBinContent(binx, biny, R);
			calibGraph->SetBinError(binx, biny, dR);
		}
		delete hR;
	}
	delete avg_pt;
	delete avg_eta;
	return calibGraph;
}

TGraph2DErrors* GetCalibration2DGraph(TFile *f, TString jetAlgo, TString method, bool numInv, double pTmax) {

	TH2D *avg_pt  = GetHisto<TH2D>(f, jetAlgo + "_AvgPtAvg_vs_EtaPtAvg_" + method);
	TH2D *avg_eta = GetHisto<TH2D>(f, jetAlgo + "_AvgEta_vs_EtaPtAvg_" + method);

	TGraph2DErrors *calibGraph = new TGraph2DErrors(); 

	int Npt  = avg_pt->GetNbinsY();
	int Neta = avg_pt->GetNbinsX();

	for (int ipt = 1 ; ipt <= Npt ; ++ipt) {

		double ptmin = avg_pt->GetYaxis()->GetBinLowEdge(ipt);
		double ptmax = avg_pt->GetYaxis()->GetBinLowEdge(ipt + 1);

		if (pTmax > 0 && pTmax <= ptmin) { continue; }

		TString Rhname = jetAlgo + "/" + jetAlgo + "_" + Form("pt%.0fto%.0f_R_", ptmin, ptmax) + method;
		TH1D *hR = GetHisto<TH1D>(f, Rhname);

		// cut the Eta bin edges of SM distribution to same one of MM
		int ieta_min = 1, ieta_max = Neta;
		if (method == "SM") {
			TH1D *hR_MM;
			TString Rhname_MM = jetAlgo + "/" + jetAlgo + "_" + Form("pt%.0fto%.0f_R_", ptmin, ptmax) + "MM";
			hR_MM = GetHisto<TH1D>(f, Rhname_MM);
			//get the first non negative R
			while (hR_MM->GetBinContent(ieta_min) <= 0 && ieta_min < Neta) {
				++ieta_min;
			}
			//get the last non negative R
			while (hR_MM->GetBinContent(ieta_max) > 0 && ieta_max > 1) {
				--ieta_max;
			}
			delete hR_MM;
		}

		if (hR->GetNbinsX() != Neta) { error("Binning problem"); }
		for (int ieta = 1 ; ieta <= Neta ; ++ieta) {
			int n = calibGraph->GetN(); 
			double R = hR->GetBinContent(ieta);

			// cut out of range bins
			if (ieta < ieta_min || ieta > ieta_max) { R = -99; }

			if (R < 0) {
				calibGraph->SetPoint(n, -99, -99, -99);
				continue;
			}

			double dR      = hR->GetBinError(ieta);
			double pt      = avg_pt->GetBinContent(ieta, ipt);
			double err_pt  = avg_pt->GetBinError(ieta, ipt);
			double eta     = avg_eta->GetBinContent(ieta, ipt);
			double err_eta = avg_eta->GetBinError(ieta, ipt);

			//printf("Pt,Eta bins: (%2d,%2d), (%6.2f +/- %5.2f)%%\n", ipt, ieta, R * 100 - 100, dR * 100);
			//printf("%.2f < %.2f\n", avg_pt->GetXaxis()->GetBinLowEdge(ieta), eta);
			//printf("%.2f < %.2f\n", avg_pt->GetYaxis()->GetBinLowEdge(ipt), pt);
			if (numInv) {
				// pTavg = (pTprobe + pTref) / 2  and  R = pTprobe / pTref   => pTprobe = 2R / (R + 1) pTavg
				pt *= 2.0 * R / (R + 1); // pTavg -> pTprobe
				// from error propagation:
				err_pt = sqrt(pow(err_pt * 2.0 * R / (R + 1), 2) + pow(2.0 * pt / (R + 1) / (R + 1) * dR, 2));
			}
			//printf("( ptmin : ptmax :n : pt : eta : R )  = ( %.0f : %.0f : %d : %.0f : %.1f : %.2f )\n", ptmin, ptmax, n, pt, eta, R);

			calibGraph->SetPoint(n, pt, eta, R);
			calibGraph->SetPointError(n, err_pt, err_eta, dR);
		}
		delete hR;
	}
	delete avg_pt;
	delete avg_eta;
	return calibGraph;
}

TGraph2DErrors* GetBinWidth2DGraph(TFile *f, TString jetAlgo, TString method, double pTmax = 0) {

	TGraph2DErrors *calibGraph = new TGraph2DErrors();
	TH2D *avg_pt  = GetHisto<TH2D>(f, jetAlgo + "_AvgPtAvg_vs_EtaPtAvg_" + method);
	TH2D *avg_eta = GetHisto<TH2D>(f, jetAlgo + "_AvgEta_vs_EtaPtAvg_" + method);

	f->cd(jetAlgo);
	TList *l = gDirectory->GetListOfKeys();

	int Npt = avg_pt->GetNbinsY(), Neta=avg_pt->GetNbinsX();

	for (int ipt = 1 ; ipt <= Npt ; ++ipt) {

		double ptmin = avg_pt->GetYaxis()->GetBinLowEdge(ipt);
		double ptmax = avg_pt->GetYaxis()->GetBinLowEdge(ipt + 1);

		if (pTmax > 0 && pTmax <= ptmin) { continue; }

		TString Rhname = jetAlgo + "/" + jetAlgo + "_" + Form("pt%.0fto%.0f_R_", ptmin, ptmax) + method;
		TH1D *hR = GetHisto<TH1D>(f, Rhname);

		// cut the Eta bin edges of SM distribution to same one of MM
		int ieta_min = 1, ieta_max = Neta;
		TH1D *hR_MM;

		if (method == "SM") {
			TString Rhname_MM = jetAlgo + "/" + jetAlgo + "_" + Form("pt%.0fto%.0f_R_", ptmin, ptmax) + "MM";
			hR_MM = GetHisto<TH1D>(f, Rhname_MM);
			//get the first non negative R
			for (int ieta = 1 ; ieta <= Neta ; ++ieta) {
				ieta_min = ieta;
				if (hR_MM->GetBinContent(ieta) > 0) { break; }
			}
			//get the last non negative R
			for (int ieta = Neta ; ieta >= 1 ; --ieta) {
				ieta_max = ieta;
				if (hR_MM->GetBinContent(ieta) > 0) { break; }
			}
			delete hR_MM;
		}

		TGraphErrors *gR = nullptr;
		int gCount = 0;
		if (method == "MM") {
			TString gRname;
			int lCount = 0;
			while (true) {
				TObject *obj = l->At(lCount);
				if (obj == nullptr) {
					std::cout << "error: cant find object in list of keys" << std::endl;
					abort();
				}
				gRname = obj->GetName();
				if ((gRname.Contains(jetAlgo + "_HLT") || gRname.Contains(jetAlgo + "_MC")) && gRname.Contains(Form("_R_vs_eta_pt%.0fto%.0f", ptmin, ptmax))) { break; }
				++lCount;
			}
			gR = (TGraphErrors*)Get(f, jetAlgo, gRname);
		} else if (method != "SM") {
			std::cout << "unknown method, aborting" << std::endl;
			abort();
		}

		for (int ieta = 1 ; ieta <= Neta ; ++ieta) {
			int n = calibGraph->GetN();
			double R  = hR->GetBinContent(ieta);
			double dR = hR->GetBinError(ieta);

			// cut out of range bins
			if (ieta < ieta_min || ieta > ieta_max) { R = -99; }

			if (R < 0) {
				calibGraph->SetPoint(n, -99, -99, -99);
				continue;
			}
			double Y, EY;
			if (method == "MM") {
				Y = gR->GetX()[gCount];
				EY = gR->GetEX()[gCount];
				++gCount;
			} else {
				Y = hR->GetXaxis()->GetBinCenter(ieta);
				EY = 0.5 * hR->GetXaxis()->GetBinWidth(ieta);
			}
			calibGraph->SetPoint(n, 0.5 * (ptmin + ptmax), Y, R);
			calibGraph->SetPointError(n, 0.5 * (ptmax - ptmin), EY, dR);
		}
		delete hR;
		if (method == "MM") { delete gR; }
	}
	delete avg_pt;
	delete avg_eta;
	return calibGraph;
}

void PrintPoint(TGraph2DErrors *g, int i) {
	printf("Point %3d (x,y,z) = (%6.1f,%6.2f,%7.4f ) +/- (%6.1f,%6.2f,%7.4f)\n",
			i, g->GetX()[i], g->GetY()[i], g->GetZ()[i], g->GetEX()[i], g->GetEY()[i], g->GetEZ()[i]);
}

TGraph2DErrors *GetAvg(TGraph2DErrors *g1, TGraph2DErrors *g2) {
	if (g1->GetN() != g2->GetN()) { error("GetAvg error."); }
	TGraph2DErrors *g = new TGraph2DErrors();
	for (int i = 0 ; i < g1->GetN() ; ++i) {
		double x1 = g1->GetX()[i], y1 = g1->GetY()[i];
		double ex1 = g1->GetEX()[i], ey1 = g1->GetEY()[i];
		double z1 = g1->GetZ()[i], z2 = g2->GetZ()[i];
		double e1 = g1->GetEZ()[i], e2 = g2->GetEZ()[i];
		g->SetPoint(i, x1, y1, (z1 + z2) / 2);
		g->SetPointError(i, ex1, ey1, 0.5 * sqrt(e1 * e1 + e2 * e2));
	}
	return g;
}

TGraph2DErrors *GetWeightedAvg(TGraph2DErrors *g1, TGraph2DErrors *g2) {
	if (g1->GetN() != g2->GetN()) { error("GetAvg error."); }
	TGraph2DErrors *g = new TGraph2DErrors();
	for (int i = 0 ; i < g1->GetN() ; ++i) {
		double x1 = g1->GetX()[i], y1 = g1->GetY()[i];
		double ex1 = g1->GetEX()[i], ey1 = g1->GetEY()[i];
		double z1 = g1->GetZ()[i], z2 = g2->GetZ()[i];
		double e1 = g1->GetEZ()[i], e2 = g2->GetEZ()[i];
		g->SetPoint(i, x1, y1, (z1 * e1 + z2 * e2) / (e1 + e2));
		g->SetPointError(i, ex1, ey1, 0.5 * sqrt(e1 * e1 + e2 * e2));
	}
	return g;
}

//******************************************************************
//**************** functions to retrive calibration points histos and values
//******************************************************************

void PrintVecI(VecI vecToPrint, TString startMessage = "Printing vector of integers: ") {
	std::cout << startMessage;
	for (unsigned int i = 0 ; i < vecToPrint.size() ; ++i) {
		std::cout << vecToPrint[i] << ",";
	}
	std::cout << std::endl;
}

VecI GetRvsPtPoints(TGraph2DErrors *g2d, double etamin, double etamax, double maxPt = 1500) {
	VecI res;
	for (int pi = 0 ; pi < g2d->GetN() ; ++pi) {
		double eta = g2d->GetY()[pi];
		if (g2d->GetX()[pi] > maxPt) { continue; }
		if (etamin < eta && eta < etamax) { res.push_back(pi); }
	}
	return res;
}

VecI GetRvsEtaPoints(TGraph2DErrors *g2d, double ptmin, double ptmax) {
	VecI res;
	for (int pi = 0 ; pi < g2d->GetN() ; ++pi) {
		double pt = g2d->GetX()[pi], R = g2d->GetZ()[pi];
		double sf = 2.0 * R / (R + 1);
		if (ptmin < pt / sf && pt / sf < ptmax) { res.push_back(pi); }
	}
	return res;
}
/*
   TGraphErrors* GetEtaPoints(TGraph2DErrors* g2d, double pT) {
   TGraphErrors *res = new TGraphErrors(); double sumpt = 0;
//get bin range
int xRange = -1;
for (int i = 0 ; i < g2d->GetNpx() ; ++i) {
if (g2d->GetX()[i] < pT) { xRange = i; }
}

TAxis *xaxis = g2d->GetXaxis();
Int_t binx = xaxis->FindBin(x);
Int_t biny = yaxis->FindBin(y);
}*/
TGraphErrors *GetRvsEtaGraph(TGraph2DErrors *g2d, VecI points) {
	TGraphErrors *res = new TGraphErrors();
	//double sumpt = 0; // Unused variable
	for (unsigned int pi = 0 ; pi < points.size() ; ++pi) {
		int n = res->GetN();
		int index = points[pi];
		res->SetPoint(n, g2d->GetY()[index], g2d->GetZ()[index]);
		res->SetPointError(n, g2d->GetEY()[index], g2d->GetEZ()[index]);
	}
	return res;
}

// g2d is the TGraph2DErrors from which you want to extract points. This is not necessarily a calibration (can contain uncertainties for instance), sor the scalefactor has to be obtained an other way.
// g2d_calib is the TGraph2DErrors that contains calibration points, used to obtain the scalefactor.
// g2d and g2d_calib do not necessarily have the same number of points.
// There should be maximum one point per bin in each TGraph2DErrors. This is not necessarily the case if the scalefactor is big.
TGraphErrors *GetRvsEtaGraph(TGraph2DErrors *g2d, TGraph2DErrors *g2d_calib, VecD etaBins, double ptmin, double ptmax) {
	TGraphErrors *res = new TGraphErrors();
	for (unsigned int ieta = 0 ; ieta < etaBins.size() - 1 ; ++ieta) {
		double etamin = etaBins[ieta];
		double etamax = etaBins[ieta + 1];
		double sf_calib = -99;
		bool   calibPointFound = false;
		// Find out if there is a calibration value in this (pT, eta) bin, and retrieve the calibration scalefactor associated to it (that is used to move from pTavg to pT)
		for (int pi = 0 ; pi < g2d_calib->GetN() ; ++pi) {
			double pt  = g2d_calib->GetX()[pi];
			double eta = g2d_calib->GetY()[pi];
			double R   = g2d_calib->GetZ()[pi];
			double sf  = 2.0 * R / (R + 1);
			if (ptmin < pt / sf && pt / sf < ptmax && etamin < eta && eta < etamax) {
				sf_calib = sf;
				if (calibPointFound) {
					std::cout << "Warning: more than one calibration points have been found in this bin (ptmin, ptmax, etamin, etamax) = " << ptmin << ", " << ptmax << ", " << etamin << ", " << etamax << std::endl;
				}
				calibPointFound = true;
			}
		}
		if (calibPointFound) {
			bool pointFound = false;
			// Find out if there is a point in g2d that falls into this (pT / sf_calib, eta) bin, and store
			for (int pi = 0 ; pi < g2d->GetN() ; ++pi) {
				double pt  = g2d->GetX()[pi];
				double eta = g2d->GetY()[pi];
				if (ptmin < pt / sf_calib && pt / sf_calib < ptmax && etamin < eta && eta < etamax) {
					int n = res->GetN();
					res->SetPoint(n, g2d->GetY()[pi], g2d->GetZ()[pi]);
					res->SetPointError(n, g2d->GetEY()[pi], g2d->GetEZ()[pi]);
					if (pointFound) {
						std::cout << "Warning: more than one points have been found in this bin (ptmin, ptmax, etamin, etamax) = " << ptmin << ", " << ptmax << ", " << etamin << ", " << etamax << std::endl;
					}
					pointFound = true;
				}
			}
		}
	}
	return res;
}

TGraphErrors *GetRvsEtaGraph(TGraph2DErrors *g2d, double ptmin, double ptmax) {
	// fix the pt range
	TGraphErrors *res = new TGraphErrors(); double sumpt = 0;
	for (int pi = 0 ; pi < g2d->GetN() ; ++pi) {
		double pt = g2d->GetX()[pi];
		double R = g2d->GetZ()[pi];
		double sf = 2.0 * R / (R + 1);
		if (ptmin < pt / sf && pt / sf < ptmax) {
			int n = res->GetN(); sumpt += pt;
			res->SetPoint(n, g2d->GetY()[pi], g2d->GetZ()[pi]);
			res->SetPointError(n, g2d->GetEY()[pi], g2d->GetEZ()[pi]);
		}
	}
	//int n = res->GetN();
	//  avg_pt = n > 0 ? sumpt / n : 0;
	return res;
}

TGraphErrors *GetRvsPtGraph(TGraph2DErrors *g2d, double etamin, double etamax) {
	TGraphErrors *res = new TGraphErrors(); double sumeta = 0;
	for (int pi = 0 ; pi < g2d->GetN() ; ++pi) {
		double eta = g2d->GetY()[pi];
		if (etamin < eta && eta < etamax) {
			int n = res->GetN(); sumeta += eta;
			res->SetPoint(n, g2d->GetX()[pi], g2d->GetZ()[pi]);
			res->SetPointError(n, g2d->GetEX()[pi], g2d->GetEZ()[pi]);
		}
	}
	return res;
}

TGraphErrors *GetRvsPtGraph(TGraph2DErrors *g2d, VecI points) {
	TGraphErrors *res = new TGraphErrors();
	for (unsigned int pi = 0 ; pi < points.size() ; ++pi) {
		int n = res->GetN();
		int index = points[pi];
		res->SetPoint(n, g2d->GetX()[index], g2d->GetZ()[index]);
		res->SetPointError(n, g2d->GetEX()[index], g2d->GetEZ()[index]);
	}
	return res;
}

// g2d is the TGraph2DErrors from which you want to extract points. This is not necessarily a calibration (can contain uncertainties for instance), sor the scalefactor has to be obtained an other way.
// g2d_calib is the TGraph2DErrors that contains calibration points, used to obtain the scalefactor.
// g2d and g2d_calib do not necessarily have the same number of points.
// There should be maximum one point per bin in each TGraph2DErrors. This is not necessarily the case if the scalefactor is big.
TGraphErrors *GetRvsPtGraph(TGraph2DErrors *g2d, TGraph2DErrors *g2d_calib, VecD ptBins, double etamin, double etamax, double pT_freeze = 2000) {
	TGraphErrors *res = new TGraphErrors();
	for (unsigned int ipt = 0 ; ipt < ptBins.size() - 1 ; ++ipt) {
		double ptmin = ptBins[ipt];
		double ptmax = ptBins[ipt + 1];
		double sf_calib = -99;
		bool   calibPointFound = false;
		// Find out if there is a calibration value in this (pT, eta) bin, and retrieve the calibration scalefactor associated to it (that is used to move from pTavg to pT)
		for (int pi = 0 ; pi < g2d_calib->GetN() ; ++pi) {
			double pt  = g2d_calib->GetX()[pi];
			double eta = g2d_calib->GetY()[pi];
			double R   = g2d_calib->GetZ()[pi];
			double sf  = 2.0 * R / (R + 1);
			if (ptmin < pt / sf && pt / sf < std::min(ptmax, pT_freeze) && etamin < eta && eta < etamax) {
				sf_calib = sf;
				if (calibPointFound) {
					std::cout << "Warning: more than one calibration points have been found in this bin (ptmin, ptmax, etamin, etamax) = " << ptmin << ", " << ptmax << ", " << etamin << ", " << etamax << std::endl;
				}
				calibPointFound = true;
			}
		}
		if (calibPointFound) {
			bool pointFound = false;
			// Find out if there is a point in g2d that falls into this (pT / sf_calib, eta) bin, and store
			for (int pi = 0 ; pi < g2d->GetN() ; ++pi) {
				double pt  = g2d->GetX()[pi];
				double eta = g2d->GetY()[pi];
				if (ptmin < pt / sf_calib && pt / sf_calib < std::min(ptmax, pT_freeze) && etamin < eta && eta < etamax) {
					int n = res->GetN();
					res->SetPoint(n, g2d->GetX()[pi], g2d->GetZ()[pi]);
					res->SetPointError(n, g2d->GetEX()[pi], g2d->GetEZ()[pi]);
					if (pointFound) {
						std::cout << "Warning: more than one points have been found in this bin (ptmin, ptmax, etamin, etamax) = " << ptmin << ", " << ptmax << ", " << etamin << ", " << etamax << std::endl;
					}
					pointFound = true;
				}
			}
		}
	}
	return res;
}

double GetAvgEta(TGraph2DErrors *g2d, VecI points) {
	double sumw = 0, sumwx = 0;
	for (unsigned int pi = 0 ; pi < points.size() ; ++pi) {
		int i = points[pi];
		double e = g2d->GetEY()[i];
		if (e <= 0) { continue; }
		double w = 1.0 / e / e;
		sumw += w;
		sumwx += w * g2d->GetY()[i];
	}
	return sumw == 0 ? 0 : sumwx / sumw;
}

double GetAvgPt(TGraph2DErrors *g2d, VecI points) {
	double sumw = 0, sumwx = 0;
	for (unsigned int pi = 0 ; pi < points.size() ; ++pi) {
		int i = points[pi];
		double e = g2d->GetEX()[i];
		if (e <= 0) { continue; }
		double w = 1.0 / e / e;
		sumw += w;
		sumwx += w * g2d->GetX()[i];
	}
	return sumw == 0 ? 0 : sumwx / sumw;
}

TH1D *GetSmoothVsEta(TH2D *h2d, double avg, VecD bins, double maxEta = 4.5) {
	static int si = 0;
	TH1D *h = new TH1D(Form("vseta%d", ++si), "", bins.size() - 1, &bins[0]);
	for (int i = 1 ; i <= h->GetNbinsX() ; ++i) {
		double eta = h->GetBinCenter(i);
		if (fabs(eta) > maxEta) { continue; }
		h->SetBinContent(i, h2d->Interpolate(avg, eta));
		double eta_edge = h->GetXaxis()->GetBinLowEdge(i);
		if (eta_edge <= 0 && eta_edge >= -0.8) { h->SetBinContent(i, 0.0); }
		if (eta_edge >= 0 && eta_edge <= 0.7)  { h->SetBinContent(i, 0.0); }
	}
	return h;
}

TH1D *GetSmoothVsEta(TH2D *h2d, double avg, double min, double max) {
	static int si = 0;
	TH1D *h = new TH1D(Form("vseta%d", ++si), "", int(5 * 10 * (max - min)), min, max); 
	for (int i = 1 ; i <= h->GetNbinsX() ; ++i) {
		double eta = h->GetBinCenter(i);
		h->SetBinContent(i, h2d->Interpolate(avg, eta));
	}
	return h;
}

TH1D *GetSmoothVsPt(TH2D *h2d, double avg, double min, double max) {
	static int si = 0;
	VecD bins = MakeLogVector(200, min, max);
	TH1D *h = new TH1D(Form("vspt%d", ++si), "", 200, &bins[0]);
	for (int i = 1 ; i <= h->GetNbinsX() ; ++i) {
		double pt = h->GetBinCenter(i);
		h->SetBinContent(i, h2d->Interpolate(pt, avg));
	}
	return h;
}

TH1D *ProjectVsEta(TH2D *h, double avgPt, double maxPt = -1) {
	int ptBin = h->GetXaxis()->FindBin(avgPt);
	if (maxPt > 0 && avgPt > maxPt) { ptBin = h->GetXaxis()->FindBin(maxPt); }
	TH1D *p = h->ProjectionY(Form("%s_%d", h->GetName(), ptBin), ptBin, ptBin);
	return p;
}


TH1D *ProjectVsPt(TH2D *h, double eta) {
	int etaBin = h->GetYaxis()->FindBin(eta);
	TH1D *p = h->ProjectionX(Form("%s_pt%d", h->GetName(), etaBin), etaBin, etaBin);
	p->SetLineColor(kViolet + 3); p->SetLineWidth(3);
	return p;
}

//************** errors etc.

void RemoveRefRegionPoints(TGraphErrors *g) {
	for (int i = 0 ; i < g->GetN() ; ++i) {
		if (g->GetX()[i] >= -0.8 && g->GetX()[i] <= 0.8) { g->SetPoint(i, g->GetX()[i], 0); g->SetPointError(i, 0); }
	}
}

TH1D* AddError(TH1D *h, TH1D *herr) {
	static int a = 0; TH1D *hh = (TH1D*)herr->Clone(Form("cerr%d", a++));
	for (int i = 1 ; i <= h->GetNbinsX() ; ++i) {
		hh->SetBinContent(i, h->GetBinContent(i));
		hh->SetBinError(i, herr->GetBinContent(i)); 
	}
	return hh;
}

TH1D *AddInQuad(TH1D *h1, TH1D *h2, TString hname = "clone") {
	TH1D *total = (TH1D*)h1->Clone(hname);
	total->Reset();
	for (int i = 1 ; i <= total->GetNbinsX() ; ++i) {

		double x1 = h1->GetBinContent(i);
		double x2 = h2->GetBinContent(i);
		double sumx2 = x1 * x1 + x2 * x2;
		total->SetBinContent(i, sqrt(sumx2));
	}
	return total;
}

TH2D *AddInQuad(std::vector<TH2D*> hvec, TString hname) {
	TH2D *total = (TH2D*)hvec[0]->Clone(hname);
	total->Reset();
	for (int ipt = 0 ; ipt <= total->GetNbinsX() + 1 ; ++ipt) {
		for (int ieta = 0 ; ieta <= total->GetNbinsY() + 1 ; ++ieta) {
			double z = 0, sumz2 = 0;
			for (unsigned int hi = 0 ; hi < hvec.size() ; ++hi) {
				z = hvec[hi]->GetBinContent(ipt, ieta);
				sumz2 += z * z;
			}
			total->SetBinContent(ipt, ieta, sqrt(sumz2));
		}
	}
	return total;
}

TH1D* AddError(TH1D *h, TH1D *herr, int color) {
	static int a = 0;
	TH1D *hh = (TH1D*)herr->Clone(Form("cerr%d", a++));
	for (int i = 1 ; i <= h->GetNbinsX() ; ++i) {

		double binctr = h->GetXaxis()->GetBinCenter(i);
		int errbin = herr->GetXaxis()->FindBin(binctr);
		hh->SetBinContent(i, h->GetBinContent(i));
		hh->SetBinError(i, herr->GetBinContent(errbin)); 
	}
	hh->SetFillColor(color);
	hh->SetMarkerStyle(1);
	return hh;
}

TH1D* AddErrorToUnity(TH1D *herr, int color) {
	static int a = 0;
	TH1D *hh = (TH1D*)herr->Clone(Form("cerr%d", a++));
	for (int i = 1 ; i <= herr->GetNbinsX() ; ++i) {
		hh->SetBinContent(i, 1.0);
		hh->SetBinError(i, herr->GetBinContent(i));
	}

	hh->SetFillColor(color);
	hh->SetMarkerStyle(1);
	return hh;
}

TH2D *GetCalibRatio(TH2D *h1, TH2D *h2, TString hname = "", bool applyAbs = false, TH2D *h3 = nullptr) {
	TH2D *hratio = (TH2D*)h1->Clone(hname);
	hratio->Reset();
	for (int ax = 1 ; ax <= h1->GetNbinsX() ; ++ax) {
		for (int ay = 1 ; ay <= h1->GetNbinsY() ; ++ay) {
			double Sys1 = -999, Sys2 = -999;
			Sys1 = h2->GetBinContent(ax, ay) / h1->GetBinContent(ax, ay);
			if (h3 != nullptr) { Sys2 = h3->GetBinContent(ax, ay) / h1->GetBinContent(ax, ay); }
			double sys_max = (Sys1 < Sys2) ? Sys2 : Sys1;
			double SYST = (applyAbs ? fabs(1 - sys_max) : 1 - sys_max);
			hratio->SetBinContent(ax, ay, SYST);
		}
	}
	return hratio;
}

TH1D* GetRelError(TString f_name, TString h_name, TString h_newName, double scale = 1.) {
	TH1D* h = GetHistoCopy<TH1D>(f_name, h_name);
	TH1D* hRelErr = (TH1D*)h->Clone(h_newName);
	hRelErr->Reset();
	for (int ix = 0 ; ix <= h->GetNbinsX() + 1 ; ++ix) {
		double err = 0;
		if (h->GetBinContent(ix) != 0) {
			err = h->GetBinError(ix) / h->GetBinContent(ix) * scale;
		}
		hRelErr->SetBinContent(ix, err);
	}
	delete h;

	return hRelErr;
}

TH1D* MakeAbsolute(TH1D* h, double scale = 1.) {
	TH1D* hAbs = (TH1D*)h->Clone((TString)h->GetName() + "_abs");
	for (int ix = 0 ; ix <= h->GetNbinsX() + 1 ; ++ix) {
		hAbs->SetBinContent(ix, std::abs(h->GetBinContent(ix)) * scale);
		hAbs->SetBinError(  ix,          h->GetBinError(  ix)  * scale);
	}
	return hAbs;
}

//**************** Retrieve pTmin and pTmax *****************************************//

// Get pT min and pT max for one value of iter, to not plot the systematics outside of the range in which there are data points
void getPtMinMax(std::vector<float> &pTmin, std::vector<float> &pTmax, const int ny, TString path, int iter) {
	for (int iy = 0 ; iy < ny ; ++iy) {
		std::ifstream file;
		file.open((path + Form("%diter/ybin_%d/data_ptminmax.txt", iter, iy)).Data());
		std::string line;
		if (file.is_open()) {
			while (file) {
				std::getline(file, line);
				if (static_cast<TString>(line).Contains("pTmin")) {
					line.erase(0, 7);
					pTmin.push_back(std::stof(line));
				}
				if (static_cast<TString>(line).Contains("pTmax")) {
					line.erase(0, 7);
					pTmax.push_back(std::stof(line));
				}
			}
		}
		file.close();
	}

	bool verbose = false;
	if (verbose) {
		for (int iy = 0 ; iy < ny ; ++iy) {
			std::cout << "iter = " << iter << std::endl;
			std::cout << "iy = " << iy << ", pTmin = " << pTmin[iy] << ", pTmax = " << pTmax[iy] << std::endl;
		}
	}
}

// Get pT min and pT max for all values in iterALL, to not plot the systematics outside of the range in which there are data points
void getPtMinMax(std::vector<float> &pTmin, std::vector<float> &pTmax, const unsigned int ny, TString path, std::vector<int> iterALL) {
	for (int iter : iterALL) {
		// Retrieve min/max of all rapidity bins for one value of iter
		std::vector<float> pTmin_oneIter = {};
		std::vector<float> pTmax_oneIter = {};
		getPtMinMax(pTmin_oneIter, pTmax_oneIter, ny, path, iter);

		// Get the min/max of all rapidity bins, for all the iterations retrieved so far
		for (unsigned int iy = 0 ; iy < ny ; ++iy) {
			if (iy < pTmin.size()) {
				pTmin[iy] = std::min(pTmin[iy], pTmin_oneIter[iy]);
			} else {
				pTmin.push_back(pTmin_oneIter[iy]);
			}
			if (iy < pTmax.size()) {
				pTmax[iy] = std::max(pTmax[iy], pTmax_oneIter[iy]);
			} else {
				pTmax.push_back(pTmax_oneIter[iy]);
			}
		}
	}

	bool verbose = false;
	if (verbose) {
		for (unsigned int iy = 0 ; iy < ny ; ++iy) {
			std::cout << "All iterations" << std::endl;
			std::cout << "iy = " << iy << ", pTmin = " << pTmin[iy] << ", pTmax = " << pTmax[iy] << std::endl;
		}
	}
}

// Get minimum of the second elements of pairs in a map
template<typename t1, typename t2>
t2 minSecondElementFromMap(std::map<t1, t2> m) {
	return (*std::min_element(m.begin(), m.end(), [](const auto& l, const auto& r) { return l.second < r.second; })).second;
}

// Get maximum of the second elements of pairs in a map
template<typename t1, typename t2>
t2 maxSecondElementFromMap(std::map<t1, t2> m) {
	return (*std::max_element(m.begin(), m.end(), [](const auto& l, const auto& r) { return l.second < r.second; })).second;
}

// Get the index of the x bin
template<typename T>
std::vector<T> getIndices(std::vector<TH1D*> hALL, std::vector<float> x, float shift = 0.) {
	std::vector<T> indices = {};
	for (unsigned int iy = 0 ; iy < hALL.size() ; ++iy) {
		indices.push_back(static_cast<T>(hALL[iy]->GetXaxis()->FindBin(x[iy] + 1e-10)) + static_cast<T>(shift));
	}

	return indices;
}

// Get the indices of the bins corresponding to xmin and xmax
void getMinMaxIndices(std::vector<TH1D*> hALL, std::vector<int> &xBinMin, std::vector<int> &xBinMax, std::vector<float> xmin, std::vector<float> xmax) {
	if (xBinMin.size() != 0 || xBinMax.size() != 0) {
		std::cout << "Warning: xBinMin or xBinMax is not empty." << std::endl;
	}
	xBinMin.clear();
	xBinMax.clear();
	for (unsigned int iy = 0 ; iy < hALL.size() ; ++iy) {
		xBinMin.push_back(hALL[iy]->GetXaxis()->FindBin(xmin[iy] + 1e-10));
		xBinMax.push_back(hALL[iy]->GetXaxis()->FindBin(xmax[iy] + 1e-10) - 1);
	}
}

// Get the number of bins between xmin and xmax
std::vector<int> getNBins(std::vector<TH1D*> hALL, std::vector<float> xmin, std::vector<float> xmax) {
	std::vector<int> xBinMin = {};
	std::vector<int> xBinMax = {};
	getMinMaxIndices(hALL, xBinMin, xBinMax, xmin, xmax);

	std::vector<int> nBins = {};
	for (unsigned int iy = 0 ; iy < hALL.size() ; ++iy) {
		nBins.push_back(xBinMax[iy] - xBinMin[iy] + 1);
	}

	return nBins;
}

//**************** Produce histograms with other histograms *************************//

// Combine all the histograms in hALL, excluding the bins below xmin and above xmax
TH1D* combineHistos(std::vector<TH1D*> hALL, std::vector<float> xmin, std::vector<float> xmax) {
	std::vector<int> nBins = getNBins(hALL, xmin, xmax);

	int nbinsx = 0;
	for (unsigned int iy = 0 ; iy < hALL.size() ; ++iy) {
		nbinsx += nBins[iy];
	}

	TH1D* hCombined = new TH1D("", "", nbinsx, 0., (Double_t)nbinsx);

	std::vector<int> xBinMin = {};
	std::vector<int> xBinMax = {};
	getMinMaxIndices(hALL, xBinMin, xBinMax, xmin, xmax);

	int i = 1;
	for (unsigned int iy = 0 ; iy < hALL.size() ; ++iy) {
		for (int ibin = xBinMin[iy] ; ibin <= xBinMax[iy] ; ++ibin) {
			hCombined->SetBinContent(i, hALL[iy]->GetBinContent(ibin));
			hCombined->SetBinError(i, hALL[iy]->GetBinError(ibin));
			++i;
		}
	}

	return hCombined;
}

// Returns h1 - h2, set its name
template<typename THnD>
THnD* subtractTH(TString name, THnD* h1, THnD* h2) {
	THnD* h3 = static_cast<THnD*>(h1->Clone(name));
	h3->Add(h2, -1.);
	return h3;
}

// Returns h1 * h2, set its name
TH1D* multiplyTH1D(TString name, TH1D* h1, TH1D* h2) {
	TH1D* h3 = static_cast<TH1D*>(h1->Clone(name));
	h3->Multiply(h2);
	return h3;
}

// Returns h1 / h2, set its name
TH1D* divideTH1D(TString name, TH1D* h1, TH1D* h2) {
	TH1D* h3 = static_cast<TH1D*>(h1->Clone(name));
	h3->Divide(h2);
	return h3;
}

// Normalize a TH1D histo
TH1D* normalizeTH1D(TString name, TH1D* hi) {
	TH1D* ho = static_cast<TH1D*>(hi->Clone(name));
	float total = 0;
	for (int ix = 1 ; ix <= hi->GetNbinsX() ; ++ix) {
		total += hi->GetBinContent(ix);
	}
	if (total != 0) {
		ho->Scale(1 / total);
	} else {
		ho->Scale(0.);
	}
	return ho;
}

// Clone input histogram and rescale it
TH1D* cloneScale(TString name, TH1D* hi, double factor) {
	TH1D* ho = static_cast<TH1D*>(hi->Clone(name));
	ho->Scale(factor);
	return ho;
}

// Performs the operation (hi + value) * factor
TH1D* addDoubleAndScale(TString name, TH1D* hi, double value, double factor) {
	TH1D* ho = static_cast<TH1D*>(hi->Clone(name));
	ho->Reset();
	for (int i = 0 ; i <= ho->GetNbinsX() ; ++i) {
		ho->SetBinContent(i, hi->GetBinContent(i) + value);
		ho->SetBinError(  i, hi->GetBinError(  i));
	}
	ho->Scale(factor);
	return ho;
}

// Normalize each line of a TH2D histo
TH2D* normalizeLinesTH2D(TString name, TH2D* hi) {
	TH2D* ho = static_cast<TH2D*>(hi->Clone(name));
	for (int iy = 1 ; iy <= hi->GetNbinsY() ; ++iy) {
		float total = 0;
		for (int ix = 1 ; ix <= hi->GetNbinsX() ; ++ix) {
			total += hi->GetBinContent(ix, iy);
		}
		for (int ix = 1 ; ix <= hi->GetNbinsX() ; ++ix) {
			if (total != 0) {
				ho->SetBinContent(ix, iy, hi->GetBinContent(ix, iy) / total);
				ho->SetBinError(  ix, iy, hi->GetBinError(  ix, iy) / total);
			} else {
				ho->SetBinContent(ix, iy, 0.);
				ho->SetBinError(  ix, iy, 0.);
			}
		}
	}
	return ho;
}

// Get TF1 from the TF2 parametrized by (bin number, pT)
TF1* TF1FromTF2(TF2* f2, float pTmin, float pTmax) {
	TString formula = f2->GetExpFormula();
	formula.Remove(0, formula.Index('(', 1, 1, TString::ECaseCompare::kExact) + 1);
	formula.Chop();
	formula.ReplaceAll('y', 'x');
	TString name = static_cast<TString>(f2->GetName()) + "_f1";
	TF1* f1 = new TF1(name, formula, pTmin, pTmax);
	for (int i = 0 ; i < f2->GetNpar() ; ++i) {
		f1->SetParameter(i, f2->GetParameter(i));
		f1->SetParError( i, f2->GetParError(i));
	}
	//std::cout << f1->GetExpFormula() << std::endl;
	//std::cout << f1->GetExpFormula("P") << std::endl;
	return f1;
}

// Set error of all bins to 0
void setErrorsTo0(TH1D* h) {
	for (int i = 0 ; i <= h->GetNbinsX() + 1 ; ++i) {
		h->SetBinError(i, 0);
	}
}

//**************** cool method to keep track of runtimes ****************************//

TString getTime() {
	time_t aclock;
	::time(&aclock);
	return TString(asctime(localtime(&aclock)));
}

void PrintTime() {
	static bool first = true;
	static time_t start;
	if (first) { first = false; ::time(&start); }
	time_t aclock; ::time(&aclock);
	char tbuf[25]; ::strncpy(tbuf, asctime(localtime(&aclock)), 24);
	tbuf[24] = 0;
	std::cout <<  "Current time: " << tbuf << " ( " << ::difftime(aclock, start) << " s elapsed )" << std::endl;
}

#endif // #ifdef Utils_h
