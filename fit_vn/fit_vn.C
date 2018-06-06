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
}

void fit_vn()
{

}