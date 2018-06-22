//------------------------------------------------------------
// This code solves a differential equation numerically to 
// create a mapping from a uniform probability distribution
// to P(phi) = (1/2pi) (1 + 2*v_2*Cos(2*phi)).
// 
// A mapping function is computed within [0, 2pi] for a wide
// range of v2 values, within [0.02, 0.3]
//------------------------------------------------------------

#include <iostream>

using namespace std;

//-----------------------------------
// Variables
//-----------------------------------

//Value of v2 for which we want to compute the mapping
float v2;

//Value of phi for which we want to compute the mapping
float phi;

//Range of v2 values to compute mapping
const float V2LO = 0.02;
const float V2HI = 0.3;

//Step size in v2
float v2_step_size = 0.01;

//Number of steps 
int n_v2_steps = (V2HI - V2LO)/v2_step_size;

//Mapping from phi to phi' for a range of v2 values
std::vector<TH1F*> h_map;

//Histogram to document which v2 value each entry in h_map corresponds to
TH1F *h_index;

//Number of phi values for every map
const int NPHIVALS = 100;

//Step size for secant method
const float STEPSIZE = 1E-3;

//-----------------------------------
// Functions
//-----------------------------------

float f(float x)
{
	//This transcendental equation comes from solving the differential equation described in documentation.pptx
	float f = x + v2*TMath::Sin(2*x) - phi - 0.5;
	return f;
}


float secant(float x1, float x2, float E)
{
	//Code from https://www.geeksforgeeks.org/program-to-find-root-of-an-equations-using-secant-method/

	float n = 0, xm, x0, c;
	if (f(x1) * f(x2) < 0) {
		do {
			// calculate the intermediate value
			x0 = (x1 * f(x2) - x2 * f(x1)) / (f(x2) - f(x1));

			// check if x0 is root of equation or not
			c = f(x1) * f(x0);

			// update the value of interval
			x1 = x2;
			x2 = x0;

			// update number of iteration
			n++;

			// if x0 is the root of equation then break the loop
			if (c == 0)
				break;
			xm = (x1 * f(x2) - x2 * f(x1)) / (f(x2) - f(x1));
		} while (fabs(xm - x0) >= E); // repeat the loop
		// until the convergence

		return x0;
		//cout << "Root of the given equation=" << x0 << endl;
		//cout << "No. of iterations = " << n << endl;
	}
	else{
		return -9999;
		//cout << "Can not find a root in the given inteval";
	}
}

void compute_flow_mapping()
{
	//Initialize index, which tells you what v2 values a given mapping corresponds to
	h_index = new TH1F("h_index", "", n_v2_steps, -0.5, n_v2_steps - 0.5);

	//Loop over the range of desired v2 values
	for(int i=0; i<n_v2_steps; i++)
	{
		//Initialize mapping histogram to be filled for the current value of v2
		TH1F *h_map_aux = new TH1F(Form("h_map_%i", i), ";#phi;#phi '", NPHIVALS, 0, 2*TMath::Pi());

		//Current value of v2
		v2 = V2LO + v2_step_size*i;
		h_index->SetBinContent(i+1, v2);

		cout << "Processing Step " << i << " for v2 = " << v2 << endl;

		//Loop over all values of phi to populate the mapping
		for(int j=0; j<NPHIVALS; j++)
		{
			phi = (float) j*2*TMath::Pi()/NPHIVALS;

			//Find root of transcendental equation f(x)
			float phi_prime = secant(0, 20, STEPSIZE);
			h_map_aux->SetBinContent(j+1, phi_prime);
 		}

		//Save mapping
		h_map.push_back(h_map_aux);
	}	

	//Write mapping out to file
	TFile *fout = new TFile("flow_maps.root", "RECREATE");
	for(int i=0; i<h_map.size(); i++)
	{
		h_map[i]->Write();
	}

	h_index->Write();
	fout->Close();
}