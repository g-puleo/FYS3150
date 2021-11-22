#include<iostream>
#include<string>
#include<fstream>
#include<iomanip>
#include<cmath>
#include "omp.h"  // OpenMP header
#include "Lattice.hpp"

using namespace std;
using namespace arma;

int main (int argc, const char* argv[])
{

  if (argc!=3)  //we force the user to input exactly TWO command line arguments.
  //they are the total number of MC cycles and the dimension of the lattice.
  {
    cout << "Bad usage: you should input 2 arguments" << endl;
    cout << "proper use: " << string(argv[0]) << " <Ncycles> <L>" << endl;
    cout << "<Ncycle> = number of Monte Carlo cycles" << endl;
    cout << "<L> = dimension of the lattice (LxL)" << endl;

    //cout << "<T> = temperature" << endl;
    exit(1);
  }
  //converting command line input to variables
  const int N_MC_cycles = atoi(argv[1]);
  const int burn_in_time = 25000;
  int L = atoi(argv[2]);
  //double T = atof(argv[3]);  //We could take T as a command line input.
  
  //T_start and T_end are the limits of the temperature interval that is going
  //to be scanned.
  double T_start = 2.25;
  double T_end = 2.35;
  int N_T =3; //number of temperatures scanned, including
  //the interval limits. Note that the total number of scans will be N_T+1
  vec T_array = vec (N_T+1); //set up an array containing the temperatures we want to scan
  for (int i=0; i<N_T+1; i++) {
  	T_array(i) = T_start + (T_end-T_start)/N_T*i;
  }

  //the following array contains the possible energy changes occurring in one spin flip.
  int N = L*L; //total number of spins
  bool order = false; //if true, starts from an ordered grid (max energy, improbable state)

  // Parallelize from here
  string ordlabel; //again, editing file name
  if (order)
    ordlabel = "O"; //ordered configuration
  else
    ordlabel = "D"; //disordered configuration

  //start parallel region
  #pragma omp parallel
  {
  //each thread must have its own private i and T
  int i;
  double T;
  const int thread_num = omp_get_thread_num(); //get thread number
  #pragma omp for
  for (i=0; i<N_T+1; i++)
  {

	  T = T_array(i);
    //initialize array of ratios of probabilities : p(s')/p(s_i).
    //this depends on temperature and is passed to the Lattice constructor.
    const double p_ratio[] = {exp(8.0/T), exp(4.0/T), 1.0, exp(-4.0/T), exp(-8.0/T)}; //units of J
    //initialize one lattice, seed depends on thread number.
  	Lattice l1 = Lattice(L, order, time(NULL)-191121*thread_num);

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
         if (!(j%100000) )//prints some output every 100 000 MC cycles
         {
           current_status = "thread " + to_string(thread_num) + ": j=" + to_string(j);
           cout << current_status << endl;
         }

          l1.evolve_Metropolis(p_ratio);//this corresponds to one MC cycle (N=L^2 attempted flips)

          //update expected values of E/spin and M/spin, their squares ,specific heat Cv ,susceptibility chi
          e_avg = (double) ( e_avg*(j) + ((double) l1.E_/N)  )  /   (double)  (j+1) ;
          m_avg = (double) ( m_avg*(j) + ((double) abs(l1.M_)/N)  )  /  (double) (j+1);

          E2_avg =  (double)(E2_avg*(j) + pow(((double) l1.E_), 2.0) )  /   (double) (j+1);
          M2_avg =  (double)(M2_avg*(j) + pow(((double) l1.M_), 2.0) )  /  (double)(j+1);

          Cv = ( E2_avg - N*N*e_avg*e_avg ) / (N*T*T);
          chi = (M2_avg - N*N*m_avg*m_avg) / (N*T);

    ///NOTE WELL: UNCOMMENT THE FOLLOWING OUTPUT STATEMENT IF YOU WANT ESTIMATES AS FUNCTION OF THE NUMBER OF
    ///MC CYCLES. IF THIS NUMBER IS LARGE, THIS COULD PRINT OUT A PRETTY HEAVY FILE
      /*   outputf << setprecision(8) << T << "," << e_avg << "," << m_avg << ","
              << Cv << "," << chi << "," << l1.E_ << "," << l1.M_ << endl; */
   	}
  ///THIS OUTPUT STATEMENT ONLY EXECUTES ONE TIME AT THE END OF THE SIMULATION, SO IT PRINTS TO A FILE
  ///ONLY THE LAST ESTIMATES FOUND.
   outputf << setprecision(8) << T << ","<< e_avg << "," << m_avg << ","
            << Cv << "," << chi <<  endl;

   } // End parallelized for loop
 } //end of parallel region
   return 0;
}
