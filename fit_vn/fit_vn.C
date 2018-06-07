//----------------------------------------------------
// Take published v2(pT) for pizeros (PHENIX PPG129)
// and fit the published points with a smooth curve
// for each centrality.
//
// JDOK
// 06-06-18
//----------------------------------------------------

#include <iostream>

using namespace std;

//---------------------------------
// Variables
//---------------------------------

//Published data points
TGraphErrors *g_v2_0_20;
TGraphErrors *g_v2_20_40;
TGraphErrors *g_v2_40_60;

//Fit functions
TF1 *f_v2_0_20;
TF1 *f_v2_20_40;
TF1 *f_v2_40_60;

//---------------------------------
// Functions
//---------------------------------

/*
 * Take published data points and construct TGraphError objects
 */
void constructTGraphs()
{
	//v2(pT) for 0-20% centrality
	float pT_0_20[12] = {1.346590909,
	                     1.851136364,
	                     2.354545455,
	                     2.825,
	                     3.325568182,
	                     3.822159091,
	                     4.289772727,
	                     4.788636364,
	                     5.535227273,
	                     6.5,
	                     7.469886364,
	                     8.498295455
	                    };

	float v2_0_20[12] = {0.069090909,
	                     0.083636364,
	                     0.094545455,
	                     0.1,
	                     0.101818182,
	                     0.090909091,
	                     0.087272727,
	                     0.083636364,
	                     0.072727273,
	                     0.06,
	                     0.063636364,
	                     0.054545455
	                    };

	float v2_ex[12] = {0.0};
	float v2_ey[12] = {0.0};

	g_v2_0_20 = new TGraphErrors(12, pT_0_20, v2_0_20, v2_ex, v2_ey);
}


/*
 * Fit v2(pT) for the 0-20% centrality
 */
void fit_v2_0_20()
{
	f_v2_0_20 = new TF1("f_v2_0_20", "pol5", 0, 10);
	g_v2_0_20->Fit("f_v2_0_20", "Q0R");
}


/*
 * Draw published points and fits
 */
void draw()
{
	g_v2_0_20->SetMarkerStyle(20);
	g_v2_0_20->SetMarkerColor(kBlue);

	TCanvas *c_v2_0_20 = new TCanvas("c_v2_0_20", "v_{2} 0-20", 500, 500);
	g_v2_0_20->SetTitle("");
	g_v2_0_20->GetYaxis()->SetRangeUser(0.0, 0.4);
	g_v2_0_20->GetYaxis()->SetTitle("v_{2}");
	g_v2_0_20->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	g_v2_0_20->Draw("AP");
	f_v2_0_20->Draw("same");
}

void fit_vn()
{
	constructTGraphs();
	fit_v2_0_20();
	draw();

	//Save functions
	TFile *fout = new TFile("f_v2_param.root", "RECREATE");
	f_v2_0_20->Write();
	fout->Close();
}