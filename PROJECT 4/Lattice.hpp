#ifndef __Lattice_hpp__
#define __Lattice_hpp__
#include<armadillo>
class Lattice
{

public:

  arma::Mat<short int> spin_grid_;//array containing all of the spins
  int L_; //  number of spins
  int E_; //total energy
  int M_; //total magnetization
  Lattice(int L_in, bool ord); //constructor

  //funcitons which computes for you the total energy and the magnetization
  double tot_energy( void );
  double tot_magnetization( void );
  void evolve_Metropolis( const double* );
};

#endif
