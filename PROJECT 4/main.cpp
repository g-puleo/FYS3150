#include<iostream>
#include<string>
#include<fstream>
#include<iomanip>
#include<cmath>
#include<chrono>
#include "omp.h"  // OpenMP header
#include "Lattice.hpp"

using namespace std;
using namespace arma;

int main (int argc, const char* argv[])
{
  //first input argument is temperature
  if (argc!=3)  //change this to argc!=4 if you want also to take temperature as input
  {
    cout << "Bad usage: you should input 3 arguments" << endl;
    cout << "proper use: " << string(argv[0]) << "<Ncycles> <L>" << endl;
    cout << "<Ncycle> = number of montecarlo cycles" << endl;
    cout << "<L> = dimension of the lattice (LxL)" << endl;

    //cout << "<T> = temperature" << endl;
    exit(1);
  }
  //converting command line input to variables
  const int N_MC_cycles = atoi(argv[1]);
  const int burn_in_time = 25000;
  int L = atoi(argv[2]);
  //double T = atof(argv[3]);      //uncomment if temperature is input

  double T_start = 2.25;
  double T_end = 2.35;
  int N_T =50;
  vec T_array = vec (N_T+1);
  for (int i=0; i<N_T+1; i++) {
  	T_array(i) = T_start + (T_end-T_start)/N_T*i;
  }
  //cout << T_array << endl;
  //the following array contains the possible energy changes occurring in one spin flip.
  int N = L*L; //total number of spins
  bool order = false; //if true, starts from an ordered grid (max energy, improbable state)


  //omp_set_dynamic(0);     // Explicitly disable dynamic teams
  //omp_set_num_threads(2); // Use 4 threads for all consecutive parallel regions
  // Parallelize from here
  string ordlabel; //again, editing file name
  if (order)
    ordlabel = "O"; //ordered configuration
  else
    ordlabel = "D"; //disordered configuration
  //#pragma omp parallel
  //start parallel region
  #pragma omp parallel
  {
  //each thread must have its own i and T
  int i;
  double T;
  const int thread_num = omp_get_thread_num(); //get thread number
  #pragma omp for
  for (i=0; i<N_T+1; i++)
  {
  //	cout << "threadnum " + to_string(thread_num)+ "index "+to_string(i) << endl;

	  T = T_array(i);
    const double p_ratio[] = {exp(8.0/T), exp(4.0/T), 1.0, exp(-4.0/T), exp(-8.0/T)}; //units of J

  	Lattice l1 = Lattice(L, order, time(NULL)-191121*thread_num);
     //initialize one lattice, seed depends on thread number.
  	//cout << "Initial state of grid : \n" << l1.spin_grid_(span(0,L-1), span(0,L-1)) << endl;
	//initialize variables to store various quantities
  	double e_avg = 0.0;
  	double m_avg = 0.0;
  	double E2_avg = 0.0;
  	double M2_avg = 0.0;
  	double Cv = 0.0;
  	double chi = 0.0;
  	//temp is an array containing the temperature, used to ouptut file with proper name
  	char temp[10] ;
  	sprintf( temp , "%.3f", T); //sprintf lets you adjust the number of digits, better than to_string()

	//finally create the string containing name of output file
  	string filename = "LatticeData_" + to_string(L) + "_" + string(temp) + "_" + ordlabel + ".csv";
  	ofstream outputf;
  	//cout << "Writing output to " << filename << endl;
  	outputf.open(filename);
  	//outputf << "cycles,average energy,average magnetization,specific heat,susceptibility,energy,magnetization" << endl;
    if(!i) //only print header when i==0
    outputf << "T,average energy,average magnetization,specific heat,susceptibility" << endl;

  	//float t = 0.0;

    //DISCARDING BURNIN  TIME
    for (int j=0; j<burn_in_time; j++)
    {
        l1.evolve_Metropolis(p_ratio);
    }
    cout << "thread " + to_string(thread_num) + "starts simulation for T= "  + to_string(T) << endl;
    string current_status;
  	for (int j=0; j<N_MC_cycles; j++)
    {
         if (!(j%100000) )
         {
           current_status = "thread " + to_string(thread_num) + ": j=" + to_string(j);
           cout << current_status << endl;
         }
  	 /* auto clock_start = chrono::system_clock::now();*/
    /* cout << l1.spin_grid_ << endl;
     cout << l1.E_ << endl;
     cout << l1.tot_energy() << endl;*/
          l1.evolve_Metropolis(p_ratio);//this corresponds to one MC cycle (N=L^2 attempted flips)
        /*  auto clock_now = chrono::system_clock::now();
          float currentTime = float(chrono::duration_cast <chrono::microseconds> (clock_now - clock_start).count());
          t += currentTime;*/
          //cout << "Montecarlo cycle: " << j<< "\nGrid state: \n" << l1.spin_grid_(span(0,1), span(0,1)) << "\n\n"  << endl;
          //update expected values of E/spin and M/spin, their squares ,specific heat Cv ,susceptibility chi
          e_avg = (double) ( e_avg*(j) + ((double) l1.E_/N)  )  /   (double)  (j+1) ;
          m_avg = (double) ( m_avg*(j) + ((double) abs(l1.M_)/N)  )  /  (double) (j+1);

          E2_avg =  (double)(E2_avg*(j) + pow(((double) l1.E_), 2.0) )  /   (double) (j+1);
          M2_avg =  (double)(M2_avg*(j) + pow(((double) l1.M_), 2.0) )  /  (double)(j+1);

          Cv = ( E2_avg - N*N*e_avg*e_avg ) / (N*T*T);
          chi = (M2_avg - N*N*m_avg*m_avg) / (N*T);
          //output to file
      /*   outputf << setprecision(8) << T << "," << e_avg << "," << m_avg << ","
              << Cv << "," << chi << "," << l1.E_ << "," << l1.M_ << endl; */
   	}
   outputf << setprecision(8) << T << ","<< e_avg << "," << m_avg << ","
            << Cv << "," << chi <<  endl;

  	//cout << "Elapsed Time: " << t /1000000 << " S \n";
   } // End parallelized for loop
 } //end of parallel region
   return 0;
}
