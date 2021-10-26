#include "PenningTrap9.hpp"
static const double K_e = 1.38935333e5;
//constructor
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in)
{
  B0_ = B0_in;
  V0_ = V0_in;
  d_ = d_in;
  V0d2_ = V0_in/(d_in*d_in);
  N_ = 0;
  coul_int_ = true; //by default , turns on coulomb interactions
}

 //adds a particle at the end of the "Particles_" vector
void PenningTrap::add_particle(Particle p_in)
{
  Particles_.push_back(p_in);
  N_ = Particles_.size();
  return;
}

//evaluates E at x=r(1), y=r(2), z=r(3)
arma::vec PenningTrap::ext_E_field( arma::vec r)
{
  arma::vec f = {r(0), r(1), -2.0*r(2)};
  return V0d2_*f ;
}

//evaluates B at x=r(1).... z=r(3)
arma::vec PenningTrap::ext_B_field ( arma::vec r)
{
  arma::vec f = {0.0, 0.0, 1.0};
  return B0_*f;
}
//evaluates coulomb interaction ACTING ON particle i DUE TO particle j
arma::vec PenningTrap::coulomb_interaction( int i, int j)
{
  //std::cout << "For now it is just empty. Returns (0, 0 , 0)";
  if (coul_int_)
  {
   arma::vec r = Particles_[i].r_ - Particles_[j].r_;

   return  (K_e*Particles_[i].q_*Particles_[j].q_/pow(norm(r), 3.0) )*r;
  }
  else
  {
    arma::vec force = {0, 0 ,0};
    return force;
  }
}
//evaluates total force due to ext fields on particle i
arma::vec PenningTrap::total_force_extfields(int i)
{
  arma::vec F = {0, 0 , 0}; //initialize zeros vector to be returned
  //initialize variables containing particle data
  double q = Particles_[i].q_;

  //std::cout << "Calculating external fields in funciton\ntotal_force_extfields(). charge of the particle:" << q << std::endl;
  arma::vec r = Particles_[i].position();
  arma::vec v = Particles_[i].velocity();

  //calculate fields E and B at position r
  arma::vec Br = ext_B_field(r);
  arma::vec Er = ext_E_field(r);

  //std::cout << "Z compontent of electric field : Ez = " << Er(2) << std::endl;
  //now ready to calculate force (in components, using Lorentz Force law)
  F(0) = q*(Er(0) + v(1)*Br(2) - v(2)*Br(1));
  F(1) = q*(Er(1) - v(0)*Br(2) + v(2)*Br(0));
  F(2) = q*(Er(2) + v(0)*Br(1) - v(1)*Br(0));

  return F;
}
//calculates the total force on particle i due to ALL of the other particles
arma::vec PenningTrap::total_force_particles(int i)
{
  arma::vec F_out = arma::vec(3,arma::fill::zeros); //initialize empty vector to be returned
  //
  for (int j=0; j<N_; j++)//adds to F all of the contributions from other particles
    {
      if (j==i) //particle i does not act a force on itself
        continue;
      else
      {
        F_out = F_out + coulomb_interaction(i,j);//adds the coulomb interaction to the total contribution
      }
    }
  return F_out;
}

//evaluates total force acting on particle i (due to both particles and field)
arma::vec PenningTrap::total_force(int i)
{
  arma::vec F_out = arma::vec(3, arma::fill::zeros);
  F_out = total_force_extfields(i) + total_force_particles(i);
  return F_out;
}

//creates a 3xN matrix whose i-th column is the total force acting on particle i
arma::mat PenningTrap::store_accelerations ( void )
{
  //create a 3xN matrix, the i-th column will represent the total force acting on particle i
  arma::mat forces = arma::mat(3,N_, arma::fill::zeros);
//first store all of the accelerations before evolving the system
  for (int j=0; j<N_; j++)
  {
    forces.col(j) = total_force(j)/Particles_[j].m_;
  }
  return forces;
}
arma::mat PenningTrap::store_positions(void )
{
  int N_ = Particles_.size();
  arma::mat pos = arma::mat(3,N_, arma::fill::zeros);
  for (int j=0; j<N_; j++)
  {
    pos.col(j) = Particles_[j].position();
  }
  return pos;
}
arma::mat PenningTrap::store_velocities( void )
{
  arma::mat vel = arma::mat(3,N_, arma::fill::zeros);
  for (int j=0; j<N_; j++)
  {
    vel.col(j) = Particles_[j].velocity();
  }
  return vel;
}
//evolve the system using forward euler
void PenningTrap::evolve_fw_Euler(double dt) {
  arma::mat acc = arma::mat(3,N_, arma::fill::zeros);
  acc = store_accelerations(); //have a 3xN matrix with all of the current accelerations
  //now evolve the system, using a for loop on the particles
  for (int i=0; i<N_; i++)
  {
    Particle Pi = Particles_[i]; //for user readability
    Particles_[i].r_ = Pi.r_ + dt*Pi.v_;
    Particles_[i].v_ = Pi.v_ + dt*acc.col(i);
  }
  //std::cout<<"this is just a try" << std::endl;
  return;
}

//resets the penning trap according
void PenningTrap::reset_PenningTrap( arma::mat r, arma::mat v)
{
  for (int j = 0 ; j<N_; j++)
  {
    Particles_[j].r_ = r.col(j);
    Particles_[j].v_ = v.col(j);
  }
}
//evoles the particles linearly in time
void PenningTrap::evolve_particles( double dt, arma::mat ka, arma::mat kb)
{
  for (int i = 0; i<N_; i++)
  {
    Particles_[i].v_ = Particles_[i].v_ + dt*kb.col(i);
    Particles_[i].r_ = Particles_[i].r_ + dt*ka.col(i);
  }
  return;
}
//evolve using RUNGE KUTTA 4
void PenningTrap::evolve_RK4 ( double dt)
{
  double h = dt/2;

  arma::mat vel_0 = arma::mat(store_velocities());//3xN matrix with all of the velocities

  arma::mat pos_0 = arma::mat(store_positions());//3xN matrix with all positions
//RK4 algorithm requires to evaluate k1, k2, k3, k4. Here,  each k_n  for each particle j, is
//represented by two column vectors as follows:
//k(1:3, j, n) = column vector of 3 entries (used to update position)
//k(4:6, j, n) =  column vector of 3 entries (used to update velocity)
//...and so on
  arma::cube k = arma::cube(6,N_, 4, arma::fill::zeros);
  for (int j =0; j<4; j++)//FOR LOOP ON K1 K2 K3 K4.
  {
 //stores velocities and accelerations of ALL of the particle
  k.slice(j).rows(0,2) = store_velocities();
  k.slice(j).rows(3,5) = store_accelerations();
  reset_PenningTrap(pos_0, vel_0); //at the first iteration, this does nothing.
                                  //it is important for the next iterations
  if (j==3) //when I have k4, I don't need to evolve the system according to k4.
    break;
  if (j==2) //when I have just stored k3, I want to evolve the system according to k3 and dt, not dt/2
    h=dt;
  evolve_particles(h, k.slice(j).rows(0,2) , k.slice(j).rows(3,5));
} //end of for loop for calculating different values of k

  //THIS IS THE ACTUAL EVOLUTION OF THE system
  arma::mat k_update = (1.0/6.0)*(k.slice(0)+2*k.slice(1)+2*k.slice(2) + k.slice(3));
  evolve_particles(dt, k_update.rows(0,2), k_update.rows(3,5));
}

//print some useful info (data of the particles contained)
void PenningTrap::particlesInfo (void )
{
  std::cout<< "The penning trap contains " << N_ << " particles.\n" <<std::endl;
  for ( int i=0; i<N_; i++){
    std::cout << "Particle " << i << ":\n" << std::endl;
    std::cout << "position:\n" << Particles_[i].position() <<std::endl;
    std::cout << "velocity:\n" << Particles_[i].velocity() <<std::endl;
    // std::cout << "charge:\n" << Particles_[i].q_ <<std::endl;
    // std::cout << "mass:\n" << Particles_[i].m_ << "\n\n" << std::endl;
  }
  return;
}

//just to practice a lil bit I write these functions
void PenningTrap::trapInfo ( std::string algo )
{

  std::cout << "Running program with the following parameters:\n" << std::endl;
  std::cout << "Algorithm: " << algo << std::endl;
  if (coul_int_)
    std::cout << "Coulomb interactions: on" << std::endl;
  else
    std::cout << "Coulomb interactions: off" <<std::endl;
  std::cout << "Number of particles: " << N_ << std::endl;
  std::cout << "B0 = " << B0_ << std::endl;
  std::cout << "V0 = " << V0_ << std::endl;
  std::cout << "d = " << d_ << std::endl;
  return;
}
