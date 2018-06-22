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



//Mapping from phi to phi' for a range of v2 values
TH1F *h_map[]

//-----------------------------------
// Functions
//-----------------------------------

float f(float x)
{
	// we are taking equation as x^3+x-1
	float f = x + v2*TMath::Sin(2*x) - phi - 0.5;
	return f;
}


void secant(float x1, float x2, float E)
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

		cout << "Root of the given equation=" << x0 << endl;
		cout << "No. of iterations = " << n << endl;
	} else
		cout << "Can not find a root in the given inteval";
}

void compute_flow_mapping()
{
	v2 = 0.5;
	phi = TMath::Pi();

	// initializing the values
	float x1 = 0, x2 = 2*TMath::Pi(), E = 0.0001;
	secant(x1, x2, E);
}