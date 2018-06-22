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
	float px;
	float py;
};

//Particles in a given event
std::vector<particle> event_particles;

//Parameterization of v2(pT) for each centrality
//For now, just use 0-20% central events
TF1 *f_v2_0_20;

//Mapping from uniform distribution to azimuthal flow modulation
TH1F *h_map;

//---------------------------------
// Functions
//---------------------------------

/*
 * For a given event, take the particle list, the impact parameter,
 * and the event plane angle and add flow modulations based on
 * particle pT and event centrality
 */
void processEvent(float b, float psi)
{
	/*
	cout << b << "   " << psi << endl;
	//Azimuthal modulation
	//TODO: Implement centrality determination; for now, assume everything is central
	TF1 *f_modulation = new TF1("f_modulation", "1 + 2*[0]*TMath::Cos(2*x)", 0, 2 * TMath::Pi());

	//Loop over event particles
	for (int i = 0; i < event_particles.size(); i++)
	{
		particle p = event_particles[i];
		float pT = TMath::Sqrt(p.px * p.px + p.py * p.py);
		float v2 = f_v2_0_20->Eval(pT);
		f_modulation->SetParameter(0, v2);

		//Randomize orientation of momentum in xy plane
		float phi = f_modulation->GetRandom();
		p.px = pT *TMath::Cos(phi + psi);
		p.py = pT *TMath::Sin(phi + psi);
	}

	event_particles.clear();
	*/
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

	int evtnumber = 0;
	float b, psi;

	while (getline(infile, line))
	{
		//Have we found a new event?
		if (line.c_str()[0] == 'E')
		{
			evtnumber++;

			//If this is not the first event, process the particles in the previous event
			if (evtnumber > 1)
			{
				processEvent(b, psi);
			}

			//Look two lines down for impact parameter and event plane information
			getline(infile, line);
			getline(infile, line);

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

			b = impact_par;
			psi = ep_angle;
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
			event_particles.push_back(p);
		}

		//If we reach the end of the file, process the last event
		if (line == "HepMC::IO_GenEvent-END_EVENT_LISTING")
		{
			processEvent(b, psi);
		}
	}
}