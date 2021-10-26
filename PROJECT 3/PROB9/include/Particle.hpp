//this is the particle class header file
#ifndef __Particle_hpp__
#define __Particle_hpp__
#include<armadillo>
class Particle {

public:
//stuff that can be updated by external functions
arma::vec r_ ; //position
arma::vec v_ ; //velocity
double q_ ; //it is the charge of the particle
double m_ ; //it is the mass
//constructor
Particle(double q_in , double m_in, arma::vec pos, arma::vec vel);

//method that returns the position
arma::vec position();

//method that returns the velocity
arma::vec velocity();

};
#endif
