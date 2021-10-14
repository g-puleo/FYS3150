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
 double V0 = 9.65e8; //(in these units: 1 Volt = 9.64852558 x 10^7) p
 double d = 1.0e4; //in micrometres

//choose a total time (in microsecs)
double total_time = 100.0;
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
    N = stoi(num);
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
      //for now, 0=FWD euler; 1=RK4  . I could  update the code and take this as command line input
}
if (argc>4)
{
  cout << "The program does not support more than 3 arguments. The 4th and the the next will be ignored." <<endl;
}


//SET INITIAL CONDITIONS FOR EACH PARTICLE and ADD PARTICLES TO THE PENNING TRAP
for (int j=0; j<N; j++)
{
  arma_rng::set_seed_random();
  //randomly set initial conditions for particle
    vec rstart = vec(3, fill::randn)*0.1*d;
    vec vstart = vec(3, fill::randn)*0.1*d;
    Particle p = Particle(q, m, rstart, vstart); //initializing one particle
    penn_1.add_particle( p ); //adding a particle to the PENNING TRAP
}
penn_1.trapInfo( algo); //display info of Penning Trap before simulation.

// vec rstart = {0,0,0};
// vec vstart = {0.0, 0.0, 0.001*d};
double elapsed_time = 0; //keeps track of the elapsed time
int counter = 0; //this is to count how many times the loops executes (want to overwrite old files when counter==0)


///////////////////// THIS IS THE PART WHERE THE EQUATION IS SOLVED NUMERICALLY ///////////////////
//now evolve the penning trap using whatever method selected, until the elapsed time reaches the total time set
//int counter = 0; Before evolving the system, I write current position, velocity and elapsed time to different files
//(one file for each particle: e.g. particle 1 data, calculated with RK4 will be "PenningData_RK4_1.csv")
do {
  //for loop on particles, write on output file(s) the current particle data.
  for (int j=0; j<N; j++){
    string filename = "PenningData_"  + algo + "_" + to_string(j+1) +  ".csv";
    ofstream output_file;
    if (counter == 0) //comparison between int and 0 is more reliable than comparison between double and 0.0
      output_file.open(filename, ios::trunc);//this deletes the already existing file with the same name
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
       output_file << penn_1.Particles_[j].v_(k) << ";" ;
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
  } while (elapsed_time<total_time);

  return 0;//THIS IS THE END OF MAIN
}
