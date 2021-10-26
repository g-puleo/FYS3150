#include<iostream>
#include<armadillo>
#include<fstream>
#include<cmath>
#include<string>
#include<unistd.h>
#include "Particle.hpp"
#include "PenningTrap.hpp"
using namespace std;
using namespace arma;
int main (int argc, char* argv[]) {
//initializing parameters of the problem
//string algoname = string(argv[0]);
 double m = 40.078; //(in amu) //mass of each particle
 double q = 1.0;//(in units of elementary charge) charge of each particle
 double B0 = 9.65e1; //(in these units: 1 Tesla = 96.4852558) magnetic field in the penning trap
 double V0 = 241213.1395; //(in these units: 1 Volt = 9.64852558 x 10^7)
 double d = 5.0e2; //in micrometres
 //set this to true if you want to also output files containing analytical solution
 //bool analytical_prediction = false;

//choose a total time (in microsecs)
double total_time = 500;
//choose a stepsize (in microsecs)
double h = 0.1;

//initialize penning trap
PenningTrap penn_1 = PenningTrap(B0, V0, d);

//depending on input arguments, set algorithm type, coulomb intractions, and number of particles
string algo = "RK4";
int N = 1;
if (argc>1)
{
      algo = string(argv[1]);
      if (algo!="FE" && algo!="RK4")
      {
        cout << "First input argument (algorithm type) must be either FE or RK4.\n"<<"Selecting RK4 by default." << endl;
        algo = "RK4";
      }
}
if (argc>2)
{
    string num = string(argv[2]);
    N = stoi(num); //converts string to int
}
if (argc>3)
{
   string c_int = string(argv[3]);
   if (c_int == "off")
      penn_1.coul_int_ = false;
   else if (c_int == "on")
      ;//by default coul_int_ is on. This line does nothing
   else
     cout << "Third input argument (coulomb interactions) must be either on or off.\nSelecting on by default."<<endl;
}
if (argc>4)
{
  string osc = string(argv[4]);
  if (osc == "on")
  	penn_1.oscil_E_field_ = true;
  else if ( osc == "off")
  	;
  else
   	cout << "Fourth input argument (oscillations in the electric field) must be either on or off.\n Selecting off by default."<<endl;
}
if (argc>5)
{
  cout << "The program does not support more than 4 arguments. The 5th and the the next will be ignored." <<endl;
}

arma_rng::set_seed_random(); //seeding random number generator.
//SET INITIAL CONDITIONS FOR EACH PARTICLE and ADD PARTICLES TO THE PENNING TRAP
mat rstart = mat(3, N, fill::zeros); //matrices which will store initial conditions
mat vstart = mat(3, N, fill::zeros);
for (int j=0; j<N; j++)
{
  //randomly set initial conditions for particle

   vec rs ;
   vec vs;
   rs = vec(3, fill::randu)*0.1*d;
   vs = vec(3, fill::randu)*0.1*d;
   rstart.col(j) = rs;
   vstart.col(j) = vs;

   Particle p = Particle(q, m, rstart.col(j), vstart.col(j)); //initializing one particle
   penn_1.add_particle( p ); //adding a particle to the PENNING TRAP
}
penn_1.trapInfo(algo); //display info of Penning Trap before simulation.

//initiating vectors for frequencies (omega_v) and amplitudes (amps) we are simulating
vec amps = {0.1, 0.4, 0.7};
int N_frequencies = 101;
vec omega_v = linspace(0.35, 0.55, N_frequencies);
//vec omega_v = {2*0.2194-0.001,2*0.2194, 2*0.2194+0.001};

for (int a=0; a<3; a++) {
	double amplitude = amps(a);
	penn_1.amplitude_ = amplitude;
	string ampfile = "PD_" + algo + "_AMP_" + to_string(a) + ".csv";
	ofstream ofile;
	ofile.open(ampfile);
	ofile << "Frequency;N" << endl;
	ofile.close();
	for (int b=0; b<N_frequencies; b++) {
		double frequency = omega_v(b);
		penn_1.frequency_ = frequency;
  		double elapsed_time = 0; //keeps track of the elapsed time
  		int counter = 0; //this is to count how many times the loops executes (want to overwrite old files when counter==0)

  		///////////////////// THIS IS THE PART WHERE THE EQUATION IS SOLVED NUMERICALLY ///////////////////
  		//now evolve the penning trap using whatever method selected, until the elapsed time reaches the total time set
  		//int counter = 0; Before evolving the system, I write current position, velocity and elapsed time to different files
  		//(one file for each particle: e.g. particle 1 data, calculated with RK4 will be "PenningData_RK4_1.csv")
  		do {
  			// update elapsed time in Penning trap
  			penn_1.elapsed_time_ = elapsed_time;

  			/*
    			//for loop on particles, write on output file(s) the current particle data.
    			for (int j=0; j<5; j++){

      				string filename = "PenningData_"  + algo + "_" + to_string(j+1)+ "_f_" + to_string(b)+ "_AMP_"+to_string(a)+ "off.csv";
      				ofstream output_file;
      				if (counter == 0)
      				{	//comparison between int and 0 is more reliable than comparison between double and 0.0
        				output_file.open(filename, ios::trunc);//this deletes the already existing file with the same name
        				output_file << "Time;pos_x;pos_y;pos_z;vel_x;vel_y;vel_z;" << endl;
      				}
      				else
        				output_file.open(filename, ios::app);//appends data to already existing file
      				//write the elapsed time
      				output_file << elapsed_time << ";";
      				for (int k=0; k<3; k++)
       			{
         				output_file << penn_1.Particles_[j].r_(k) << ";";
       			}
      				for(int k=0; k<3; k++)
       			{
         				output_file << penn_1.Particles_[j].v_(k) << ";";
       			}
      				output_file << endl;
      				output_file.close();


     			} //end of for loop to write on files for each particle
			*/
  			//now evolve the system using the chosen algorithm
    			if (algo=="FE")
      			{
        			penn_1.evolve_fw_Euler(h);
      			}
    			else
      			{
        			penn_1.evolve_RK4(h);
      			}
      			elapsed_time+=h;
      			counter++;
    		} while (elapsed_time<total_time); //end of "do while" loop (each iteration correspond to one time step)

    		// compute particles left in trap
		int particles_left;
		particles_left = penn_1.countParticles();
		// write number of particles to file
		ofstream ofile;
		if (counter ==0)
		{
			ofile.open(ampfile, ios::trunc);
		}
		else
			ofile.open(ampfile, ios::app);
		ofile << frequency << ";" << particles_left << endl;
		ofile.close();
 		//before starting simulation for another value of h, reset the penning trap
 		penn_1.reset_PenningTrap(rstart, vstart);
	}
}

/* Uncomment this to loop over same frequencies, but now with Coulomb interactions between the particles.
penn_1.coul_int_ = true;
for (int a=0; a<1; a++) {
	double amplitude = amps(a);
	penn_1.amplitude_ = amplitude;
	string ampfile = "PD_" + algo + "_AMP_" + to_string(a) + "on.csv";
	ofstream ofile;
	ofile.open(ampfile);
	ofile << "Frequency;N" << endl;
	ofile.close();
	for (int b=0; b<N_frequencies; b++) {
		double frequency = omega_v(b);
		penn_1.frequency_ = frequency;
  		double elapsed_time = 0; //keeps track of the elapsed time
  		int counter = 0; //this is to count how many times the loops executes (want to overwrite old files when counter==0)

  		do {
  			// update elapsed time in Penning trap
  			penn_1.elapsed_time_ = elapsed_time;

  			te on output file(s) the current particle data.
    			for (int j=0; j<5; j++){
    			//for loop on particles, write on output file(s) the current particle data.
    			for (int j=0; j<5; j++){

      				string filename = "PenningData_"  + algo + "_" + to_string(j+1)+ "_f_" + to_string(b)+ "_AMP_"+to_string(a)+ "on.csv";
      				ofstream output_file;
      				if (counter == 0)
      				{	//comparison between int and 0 is more reliable than comparison between double and 0.0
        				output_file.open(filename, ios::trunc);//this deletes the already existing file with the same name
        				output_file << "Time;pos_x;pos_y;pos_z;vel_x;vel_y;vel_z;" << endl;
      				}
      				else
        				output_file.open(filename, ios::app);//appends data to already existing file
      				//write the elapsed time
      				output_file << elapsed_time << ";";
      				for (int k=0; k<3; k++)
       			{
         				output_file << penn_1.Particles_[j].r_(k) << ";";
       			}
      				for(int k=0; k<3; k++)
       			{
         				output_file << penn_1.Particles_[j].v_(k) << ";";
       			}
      				output_file << endl;
      				output_file.close();


     			} //end of for loop to write on files for each particle

  			//now evolve the system using the chosen algorithm
    			if (algo=="FE")
      			{
        			penn_1.evolve_fw_Euler(h);
      			}
    			else
      			{
        			penn_1.evolve_RK4(h);
      			}
      			elapsed_time+=h;
      			counter++;
    		} while (elapsed_time<total_time); //end of "do while" loop (each iteration correspond to one time step)

    		// compute particles left in trap
		int particles_left;
		particles_left = penn_1.countParticles();
		// write number of particles to file
		ofstream ofile;
		if (counter ==0)
		{
			ofile.open(ampfile, ios::trunc);
		}
		else
			ofile.open(ampfile, ios::app);
		ofile << frequency << ";" << particles_left << endl;
		ofile.close();
 		//before starting simulation for another value of h, reset the penning trap
 		penn_1.reset_PenningTrap(rstart, vstart);
	}
}
*/

  return 0;//THIS IS THE END OF MAIN
}
