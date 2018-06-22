//-------------------------------------------
// Read in a HepMC file with HIJING event
// output, and impose azimuthal modulation
// on individual tracks based on their
// transverse momentum and event centrality,
// writing out a new HepMC file
//
// JDOK
// 06-06-18
//-------------------------------------------

#include <iostream>
#include <fstream>

using namespace std;

//---------------------------------
// Variables
//---------------------------------

struct particle {
	int barcode;
	int pdg_id;
	float px;
	float py;
	float pz;
	float energy;
	float mass;
	int status;
};

//Parameterization of v2(pT) for each centrality
//For now, just use 0-20% central events
TF1 *f_v2_0_20;

//Mapping from uniform distribution to azimuthal flow modulation
TH1F *h_map;

//File to write out
ofstream myfile;

//Reaction plane angle and impact parameter for the event at hand
float current_psi;
float current_b;

TH1F *h_phi = new TH1F("h_phi", "", 100, 0, 2 * TMath::Pi());
TH1F *h_phi_prime = new TH1F("h_phi_prime", "", 100, 0, 2 * TMath::Pi());

TProfile *h_v2 = new TProfile("h_v2", ";p_{T};v_{2}", 50, 0, 5, -10, 10);
TProfile *h_v2_modif = new TProfile("h_v2_modif", ";p_{T};v_{2}", 50, 0, 5, -10, 10);

//---------------------------------
// Functions
//---------------------------------

/*
 *
 */
float f(float x, float phi, float v2)
{
	//This transcendental equation comes from solving the differential equation described in documentation.pptx
	float f = x + v2 * TMath::Sin(2 * x) - phi - 0.5;
	return f;
}


/*
 *
 */
float secant(float x1, float x2, float E, float phi, float v2)
{
	//Code from https://www.geeksforgeeks.org/program-to-find-root-of-an-equations-using-secant-method/

	float n = 0, xm, x0, c;
	if (f(x1, phi, v2) * f(x2, phi, v2) < 0) {
		do {
			// calculate the intermediate value
			x0 = (x1 * f(x2, phi, v2) - x2 * f(x1, phi, v2)) / (f(x2, phi, v2) - f(x1, phi, v2));

			// check if x0 is root of equation or not
			c = f(x1, phi, v2) * f(x0, phi, v2);

			// update the value of interval
			x1 = x2;
			x2 = x0;

			// update number of iteration
			n++;

			// if x0 is the root of equation then break the loop
			if (c == 0)
				break;
			xm = (x1 * f(x2, phi, v2) - x2 * f(x1, phi, v2)) / (f(x2, phi, v2) - f(x1, phi, v2));
		} while (fabs(xm - x0) >= E); // repeat the loop
		// until the convergence

		return x0;
	}
	else {
		return -9999;
	}
}


/*
 * For a given particle, take the event impact parameter,
 * and the event plane angle and add flow modulations based on
 * particle pT and event centrality
 */
particle processParticle(particle p)
{
	float pT = TMath::Sqrt(p.px * p.px + p.py * p.py);
	float phi = TMath::ATan2(p.py, p.px);

	//The range of ATan2 is [-pi, pi], so add pi to distribute the particles within [0, 2pi]
	phi = phi + TMath::Pi();
	h_phi->Fill(phi);

	//Eventually, we will select the mapping based on the v2 corresponding to the particle pT
	//float v2 = f_v2_0_20->Eval(pT);
	float v2 = 0.5;

	//Apply mapping
	float phi_prime = secant(0, 20, 1E-3, phi, v2);

	//Wrap phi_prime around
	if (phi_prime > 2 * TMath::Pi())
	{
		phi_prime = phi_prime - 2 * TMath::Pi();
	}
	else if (phi_prime < 0)
	{
		phi_prime = phi_prime + 2 * TMath::Pi();
	}

	h_phi_prime->Fill(phi_prime);

	//p.px = pT * TMath::Cos(phi_prime + current_psi);
	//p.py = pT * TMath::Sin(phi_prime + current_psi);

	p.px = pT * TMath::Cos(phi_prime);
	p.py = pT * TMath::Sin(phi_prime);

	h_v2->Fill(pT, TMath::Cos(2 * phi));
	h_v2_modif->Fill(pT, TMath::Cos(2 * phi_prime));

	return p;
}


/*
 * Take a HepMC file, from HIJING output, and parse it line by line.
 * For every event, the final-state particles are stored in a vector and processed to add flow modulations
 */
void read_hepmc()
{
	//Load mapping from a uniform azimuthal distribution to a flow modulation
	TFile *f_mapping = new TFile("flow_maps.root");
	h_map = (TH1F*) f_mapping->Get("h_map_28");

	//Load vn parameterizations
	TFile *f_param = new TFile("fit_vn/f_v2_param.root");
	f_v2_0_20 = (TF1*) f_param->Get("f_v2_0_20");

	//Load in HepMC file and parse line by line
	ifstream infile("sHijing.dat");
	string line;

	//Open file to write out modified HepMC file
	myfile.open ("sHijing.modif.dat");

	int evtnumber = 0;
	float b, psi;

	while (getline(infile, line))
	{
		//If this is not a particle line, write it out to the modified file
		if (line.c_str()[0] != 'P')	myfile << line << endl;

		//Have we found a new event?
		if (line.c_str()[0] == 'E')
		{
			evtnumber++;

			//Look two lines down for impact parameter and event plane information
			getline(infile, line);
			myfile << line << endl;
			getline(infile, line);
			myfile << line << endl;

			string junk1;
			float n_hardscat;
			float n_proj_part;
			float n_targ_part;
			float n_coll;
			float n_spect_n;
			float n_spect_p;
			float n_wound1;
			float n_wound2;
			float n_wound3;
			float impact_par;
			float ep_angle;
			float eccen;
			float cross_sec;

			stringstream s(line);
			if (!(s >> junk1 >> n_hardscat >> n_proj_part >> n_targ_part >> n_coll >> n_spect_n >> n_spect_p >> n_wound1 >> n_wound2 >> n_wound3 >> impact_par >> ep_angle >> eccen >> cross_sec)) break;

			current_b = impact_par;
			current_psi = ep_angle;
		}

		//Have we found a particle?
		if (line.c_str()[0] == 'P')
		{
			string junk1;
			int barcode;
			int pid;
			float px, py, pz;
			float energy;
			float mass;
			int status;
			float pol_theta;
			float pol_phi;
			int barcode_vtx;
			int num_entries;

			stringstream s(line);
			if (!(s >> junk1 >> barcode >> pid >> px >> py >> pz >> energy >> mass >> status >> pol_theta >> pol_phi >> barcode_vtx >> num_entries)) break;

			//Store it in the vector containing the particles of the event at hand
			particle p;
			p.px = px;
			p.py = py;
			p.pz = pz;
			p.barcode = barcode;
			p.pdg_id = pid;
			p.energy = energy;
			p.mass = mass;
			p.status = status;

			particle p_modif = processParticle(p);
			myfile << "P " << p_modif.barcode << " " << p_modif.pdg_id << " " << p_modif.px << " " << p_modif.py << " " << p_modif.pz << " " << p_modif.energy << " " << p_modif.mass << " " << " " << p_modif.status << " 0 0 0 0" << endl;
		}
	}

	TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
	h_phi->Draw();

	TCanvas *c2 = new TCanvas("c2", "c2", 600, 600);
	h_phi_prime->Draw();

	TCanvas *c3 = new TCanvas("c3", "c3", 600, 600);
	h_v2_modif->Draw();
}