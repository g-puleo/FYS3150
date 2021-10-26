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
 double V0 = 9.65e8; //(in these units: 1 Volt = 9.64852558 x 10^7)
 double d = 1.0e4; //in micrometres
 //set this to true if you want to also output files containing analytical solution
 bool analytical_prediction = false;

//choose a total time (in microsecs)
double total_time = 20;
//choose a stepsize (in microsecs)
//this is a arma::vec bcause the code is designed to repeat simulations for possibly more than 1 value of h
vec hvalues = { 0.001 };//{0.11, 0.05, 0.023, 0.01, 0.0045};

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
  cout << "The program does not support more than 3 arguments. The 4th and the the next will be ignored." <<endl;
}

arma_rng::set_seed_random(); //seeding random number generator.
//SET INITIAL CONDITIONS FOR EACH PARTICLE and ADD PARTICLES TO THE PENNING TRAP
mat rstart = mat(3, N, fill::zeros); //matrices which will store initial conditions
mat vstart = mat(3, N, fill::zeros);
for (int j=0; j<N; j++)
{
  //randomly set initial conditions for particle

    //CASE OF 2 PARTICLES. COMMENT OUT THE TWO IFs IN ORDER TO HAVE STANDARD RUN
   //(every particle initialized at random)
   vec rs ;
   vec vs;
      if(j==0)
      { //set initial condition of first particle
        rs = {-9.773, -5.226, 10.226};
        vs = {0.0, 0.0 ,0.0};
      }
      else if (j==1)
      {//initial condition of second particle
        rs = {-4, -11, 16};
        vs= {0.0, 0.0, 0.0};
      }
      else
      { //keep only this to have standard run
        rs = vec(3, fill::randn)*0.1*d;
        vs = vec(3, fill::randn)*0.1*d;
      }
    rstart.col(j) = rs;
    vstart.col(j) = vs;

    Particle p = Particle(q, m, rstart.col(j), vstart.col(j)); //initializing one particle
    penn_1.add_particle( p ); //adding a particle to the PENNING TRAP
}
penn_1.trapInfo(algo); //display info of Penning Trap before simulation.

string labels = "ABCDE"; //different file names will be assigned to different values of h

//comment out this for loop to execute the code for just one value of h
for(int jj=0; jj<hvalues.size(); jj++)
{
  double h = hvalues(jj); //set value of h of current simulation
  double elapsed_time = 0; //keeps track of the elapsed time
  int counter = 0; //this is to count how many times the loops executes (want to overwrite old files when counter==0)

  ///////////////////// THIS IS THE PART WHERE THE EQUATION IS SOLVED NUMERICALLY ///////////////////
  //now evolve the penning trap using whatever method selected, until the elapsed time reaches the total time set
  //int counter = 0; Before evolving the system, I write current position, velocity and elapsed time to different files
  //(one file for each particle: e.g. particle 1 data, calculated with RK4 will be "PenningData_RK4_1.csv")
  do {
    //for loop on particles, write on output file(s) the current particle data.
    for (int j=0; j<N; j++){

      string filename = "PenningData_"  + algo + "_" + labels.substr(jj,1) + "_" + to_string(j+1)+ ".csv";
      ofstream output_file;
      if (counter == 0)
      { //comparison between int and 0 is more reliable than comparison between double and 0.0
        output_file.open(filename, ios::trunc);//this deletes the already existing file with the same name
      //  output_file << "First run" << endl;
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


      //if the coulomb interactions are switched off, also output files
      //containing the analytical solutions.
      //remember that this solution is correct under suitable initial conditions
      //the flag analytical_prediction allows to turn this on and off.
       if ( !penn_1.coul_int_ && analytical_prediction)
       {
         //calculating the analytical solution
         double omega_0 = q*B0/m ;
         double omega_z = sqrt(2.0*q*V0/(m*d*d));

         double omega_plus = 0.5*(omega_0 + sqrt(omega_0*omega_0 - 2.0*omega_z*omega_z) ) ;
         double omega_minus = 0.5*(omega_0 - sqrt(omega_0*omega_0 - 2.0*omega_z*omega_z) ) ;

         double A_plus = (vstart.col(j)[1] + omega_minus*rstart.col(j)[0])/(omega_minus - omega_plus);
         double A_minus = -(vstart.col(j)[1] + omega_plus*rstart.col(j)[0])/(omega_minus - omega_plus);

         complex<double> fxy = A_plus*exp(-1i*omega_plus*elapsed_time) + A_minus*exp(-1i*omega_minus*elapsed_time) ;
         //take imaginary and real part which correspond to motion along y and x
         vec r_pred = vec(3, fill::zeros);
         r_pred(0) = real(fxy);
         r_pred(1) = imag(fxy);
         //calculate also motion along z-axis
         r_pred(2) = rstart.col(j)[2]*cos(omega_z*elapsed_time);

         ////print analytical prediction to a different file called "r_ana_x.csv" for particle x

         //first I open the file in appropriate mode
         if (counter == 0) //comparison between int and 0 is more reliable than comparison between double and 0.0
           output_file.open("r_ana_"+labels.substr(jj,1)+"_" + to_string(j+1)+".csv", ios::trunc); //deletes existing file
         else
           output_file.open("r_ana_"+labels.substr(jj,1)+"_" + to_string(j+1)+".csv", ios::app); //appends to already existing file


         //then start writing the output
         output_file << elapsed_time << ";"; //store elapsed time

         for (int k=0; k<3; k++)
             {
               output_file << r_pred(k) << ";"; //store analytical prediction
             }
         output_file << norm(r_pred-penn_1.Particles_[j].r_)/norm(r_pred) << ";" ; //also store relative error
         output_file << norm(r_pred-penn_1.Particles_[j].r_) << endl ; //store absolute error
         output_file.close();
        }
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

 //before starting simulation for another value of h, reset the penning trap
 penn_1.reset_PenningTrap(rstart, vstart);

}//end of for loop over different values of h

  return 0;//THIS IS THE END OF MAIN
}
