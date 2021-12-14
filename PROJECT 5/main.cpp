#include<armadillo>
#include<iostream>
#include<string>
#include<fstream>
#include<vector>
#include"Box.hpp"
#include"utils.hpp"
using namespace arma;
using namespace std;

int main(int argc, const char* argv[])

{

	if (argc!=3) {
		cout << "ERROR. PROPER USE: " << argv[0] << " <init_cond> <potential>" << endl;
		cout << "<init_conds> must be a .csv file containing the initial conditions." << endl;
		cout << "<potential> must be a .dat file containing the shape of the potential." << endl;
		exit(1);
	}
	int M = 201; //size of box

	//argv[1] is the name of the .csv file containing gaussian wave packet (gwp)and other parameters
	string init_conditions = argv[1];
	vector<string> initial_conditions = read_file(init_conditions, ','); //import parameters
	//store various parameters as variables
	double V0 = stod(initial_conditions[0]);
	double time_step = stod(initial_conditions[1]);
	double T = stod(initial_conditions[2]);
	double x_c = stod(initial_conditions[3]);
	double y_c = stod(initial_conditions[4]);
	double s_x = stod(initial_conditions[5]);
	double s_y = stod(initial_conditions[6]);
	double mom_x = stod(initial_conditions[7]);
	double mom_y = stod(initial_conditions[8]);
	//initialize structure containing gwp parameters
  GaussParams gpars = {x_c,y_c, s_x, s_y, mom_x, mom_y};
	Box b1 = Box(M, gpars); //create box
	// the value of r is necessary to initialize A and B
	complex<double> r = time_step/(2*b1.h_*b1.h_)*1.0i;

  //initialize the potential :
  mat V;
	//argv[2] is the name of the .dat file containing the potential matrix
	string potential = argv[2];
	V.load(potential, raw_ascii);
	V = V*V0; //scale the potential using the input parameter V0

	//prepare matrices A and B, first filling vectors a and b.
	b1.fill_a_and_b(r, time_step, V);
	b1.prepare_matrices(r);

	int Ntot  = ceil(T/time_step); //initialize total number of iterations
	//just to inform the user:
	cout << "The system will be evolved for total number of " << Ntot << " time steps." << endl;

	//Now we initialize cubes that will contain the data we want to import in python

	cube mod_U_sq = cube(M, M, Ntot, fill::zeros); //squared modulus of U
	cube real_wf = cube(M, M, Ntot, fill::zeros);  //real part of U
	cube imag_wf = cube(M, M, Ntot, fill::zeros);  //imaginary part of U

	double e_time = 0.0; //keeps track of the elapsed time
	//this for loop is the actual time evolution
 	for(int ii=0; ii<Ntot ; ii++){

		if(!(ii%(Ntot/10)) ) {
			//inform user every time the program has run ~1/10 of the simulation
			 cout << "Time step n. " << ii <<  endl;
		}
		//store the current state of the absolute value of U, and also its real part and its imaginary part
		mod_U_sq.slice(ii)=real(b1.U_%conj(b1.U_));
		real_wf.slice(ii)=real(b1.U_);
		imag_wf.slice(ii)=imag(b1.U_);
	  b1.advance(); //advance the system!!
		e_time+=time_step; //update elapsed time
	}

	//finally, save results to output files.
	mod_U_sq.save("born_prob.bin", arma_binary);
	real_wf.save("real_wf.bin", arma_binary);
	imag_wf.save("imag_wf.bin", arma_binary);

	return 0;
}
