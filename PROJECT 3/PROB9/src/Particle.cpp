#include "Particle.hpp"

//Constructor
Particle::Particle(double q_in, double m_in, arma::vec pos , arma::vec vel)
{
  m_ = m_in;
  q_ = q_in;
  r_ = pos;
  v_ = vel;
}
//to return the position
arma::vec Particle::position() {
  return r_;
}
//to return the velocity
arma::vec Particle::velocity() {
  return v_;
}
