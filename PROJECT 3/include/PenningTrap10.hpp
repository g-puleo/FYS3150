#ifndef __PenningTrap_hpp__
#define __PenningTrap_hpp__
#include<vector>
#include<iostream>
#include<armadillo>
#include<string>
#include"Particle.hpp"
#include<cmath>
class PenningTrap {

public:
  double B0_ ; //magnetic field
  double V0_ ; //potential
  double d_ ;//characteristic dimension
  double V0d2_; //store also the square of d, it avoids some FLOPs
  int N_; //number of particles
  bool coul_int_; //swithces on/off the coulomb interactions
  bool oscil_E_field_; //switches on/off the oscillations in the electric field
  double amplitude_; // amplitude of oscillations
  double frequency_; // frequency of oscillations
  double elapsed_time_; //time for oscillations
  std::vector<Particle> Particles_; //vector containing all of the particles

  //constructor
  PenningTrap(double B0_in, double V0_in, double d_in);

  //member functions
  //adds a particle to the trap
  void add_particle ( Particle p_in ) ;

  //evaluates the E field at one point
  arma::vec ext_E_field( arma::vec r);

  //evaluates the B field at one point
  arma::vec ext_B_field( arma::vec r);

  //force from particle i on particle j
  arma::vec coulomb_interaction( int i, int j);

  //force on particle i due to the external fields
  arma::vec total_force_extfields( int i);

  //force on particle i due to the coulomb interaction with other particles
  arma::vec total_force_particles( int i);

  //total force on particle i
  arma::vec total_force( int i);

  //creates 3xN matrix storing all of the forces (i-th colum = force on particle i)
  arma::mat store_accelerations( void );
  //create 3XN matrix storing all of the positions
  arma::mat store_positions( void );

  //create 3xN matrix storing all of the velocities
  arma::mat store_velocities ( void );
  //reset penning trap given positions and velocities
  void reset_PenningTrap( arma::mat r, arma::mat v ) ;

  //useful for taking intermediate step in RK4. Incidentally, almost implements the whole FE algorithm.
  void evolve_particles( double dt, arma::mat k_a, arma::mat k_b);

  //Evolve one time step using Runge-Kutta 4
  void evolve_RK4 (double dt);

  //Evolve one time step using forward euler
  void evolve_fw_Euler(double dt);

  //Print useful info (I wrote this just to practice a lil bit)
  void particlesInfo ( void );
  void trapInfo ( std::string );
  int countParticles( void ); 
};

#endif
