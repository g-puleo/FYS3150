#include<iostream>
#include<string>
#include<fstream>
#include<iomanip>
#include<cmath>
#include "Lattice.hpp"

using namespace std;
using namespace arma;

int main (int argc, const char* argv[])
{
  //first input argument is temperature
  if (argc!=4)
  {
    cout << "Bad usage: you should input 3 arguments" << endl;
    cout << "proper use: " << string(argv[0]) << "<Ncycles> <L> <T> " << endl;
    cout << "<L> = dimension of the lattice (LxL)" << endl;
    cout << "<T> = temperature (units of J/kB)" << endl;
    cout << "<Ncycle> = number of montecarlo cycles" << endl;
    exit(1);
  }
  //converting command line input to variables
  int N_MC_cycles = atoi(argv[1]);
  int L = atoi(argv[2]);
  double T = atof(argv[3]);
  //the following array contains the possible energy changes occurring in one spin flip.
  const double p_ratio[] = {exp(8.0/T), exp(4.0/T), 1.0, exp(-4.0/T), exp(-8.0/T)}; //units of J
//  cout<< "array containing p_old/p_new , called p_ratio[]:\n"<< p_ratio << "\n" << endl;
// for(int jj=0; jj<5; jj++){
//   cout << p_ratio[jj] << endl;
// }
  int N = L*L; //total number of spins
  bool order = false; //if true, starts from an ordered grid (max energy, improbable state)
  Lattice l1 = Lattice(L, order); //initialize one lattice
  cout << "Initial state of grid : \n" << l1.spin_grid_(span(0,L-1), span(0,L-1)) << endl;
//initialize variables to store various quantities
  double e_avg = 0.0;
  double m_avg = 0.0;
  double E2_avg = 0.0;
  double M2_avg = 0.0;
  double Cv = 0.0;
  double chi = 0.0;
  //temp is an array containing the temperature, used to ouptut file with proper name
  char temp[10] ;
  sprintf( temp , "%.1f", T); //sprintf lets you adjust the number of digits, better than to_string()
  string ordlabel; //again, editing file name
  if (order)
    ordlabel = "O";
  else
    ordlabel = "D";
//finally create the string containing name of output file
  string filename = "LatticeData_" + to_string(L) + "_" + string(temp) + "_" + ordlabel + ".csv";
  ofstream outputf;
  cout << "Writing output to " << filename << endl;
  outputf.open(filename);
  outputf << "cycles,average energy, average magnetization, specific heat, susceptibility, energy,magnetization" << endl;
  for (int j=0; j<N_MC_cycles; j++){
          l1.evolve_Metropolis(p_ratio);//this corresponds to one MC cycle (N=L^2 attempted flips)
          //cout << "Montecarlo cycle: " << j<< "\nGrid state: \n" << l1.spin_grid_(span(0,1), span(0,1)) << "\n\n"  << endl;
          //update expected values of E/spin and M/spin, their squares ,specific heat Cv ,susceptibility chi
          e_avg = (double) ( e_avg*(j) + ((double) l1.E_/N)  )  /   (double)  (j+1) ;
          m_avg = (double) ( m_avg*(j) + ((double) abs(l1.M_)/N)  )  /  (double) (j+1);

          E2_avg =  (double)(E2_avg*(j) + pow(((double) l1.E_), 2.0) )  /   (double) (j+1);
          M2_avg =  (double)(M2_avg*(j) + pow(((double) l1.M_), 2.0) )  /  (double)(j+1);

          Cv = ( E2_avg - N*N*e_avg*e_avg ) / (N*T*T);
          chi = (M2_avg - N*N*m_avg*m_avg) / (N*T);
          //output to file
          outputf << setprecision(8) << j+1 << "," << e_avg << "," << m_avg << ","
                  << Cv << "," << chi << "," << l1.E_ << "," << l1.M_ << endl;
   }
}
