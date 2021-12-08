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
	int M = 201; //size of box
	string init_conditions = argv[1];
	vector<string> initial_conditions = read_file(init_conditions, ','); 
	string potential = argv[2]; 
	double V0 = stod(initial_conditions[0]); 
	double time_step = stod(initial_conditions[1]);
	double T = stod(initial_conditions[2]);
	double x_c = stod(initial_conditions[3]);
	double y_c = stod(initial_conditions[4]);
	double s_x = stod(initial_conditions[5]);
	double s_y = stod(initial_conditions[6]);
	double mom_x = stod(initial_conditions[7]);
	double mom_y = stod(initial_conditions[8]);
  	GaussParams gpars = {x_c,y_c, s_x, s_y, mom_x, mom_y};
	Box b1 = Box(M, gpars); //create box
	complex<double> r = time_step/(2*b1.h_*b1.h_)*1.0i;
  	
  	//initialize the potential :
  	mat V; 
	V.load(potential, raw_ascii);  
	V = V*V0; 

	//prepare matrices A and B
	b1.fill_a_and_b(r, time_step, V);
	b1.prepare_matrices(r);

	int Ntot  = ceil(T/time_step); //initialize total number of iterations
	cout << "total n. iterations" << Ntot << endl;
  	cube Ufinal = cube(M, M, Ntot, fill::zeros);
  	cube real_wf = cube(M, M, Ntot, fill::zeros); 
  	cube imag_wf = cube(M, M, Ntot, fill::zeros); 
	//std::cout <<  "DEBUG: " << __LINE__ << std::endl;
	//cout << "entering for loop" << endl;
	double e_time = 0.0;
	mat probability;
 	for(int ii=0; ii<Ntot ; ii++){
	//	cout << "ii = " << ii << endl ;
		probability = real(b1.U_%conj(b1.U_));
	//	cout << "advanced" << endl ;
		Ufinal.slice(ii)=probability;
		real_wf.slice(ii)=real(b1.U_); 
		imag_wf.slice(ii)=imag(b1.U_); 
	  	b1.advance();
		e_time+=time_step;
	}

	Ufinal.save("born_prob.bin", arma_binary);
	real_wf.save("real_wf.bin", arma_binary); 
	imag_wf.save("imag_wf.bin", arma_binary); 

	return 0;
}
