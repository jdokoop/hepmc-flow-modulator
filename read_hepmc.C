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

//---------------------------------
// Functions
//---------------------------------

void processEvent(float b, float psi)
{
	cout << b << "  " << psi <<  "   " << event_particles.size() <<endl;    
	event_particles.clear();
}


void read_hepmc()
{
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
			if(evtnumber > 1)
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
			if(!(s >> junk1 >> barcode >> pid >> px >> py >> pz >> energy >> mass >> status >> pol_theta >> pol_phi >> barcode_vtx >> num_entries)) break;

			particle p;
			p.px = px;
			p.py = py;
			event_particles.push_back(p);
		}

		//If we reach the end of the file, process the last event
		if(line == "HepMC::IO_GenEvent-END_EVENT_LISTING")
		{
			processEvent(b, psi);
		}
	}
}