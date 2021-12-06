#include<armadillo>
#include<iostream>
#include<string>
#include"Box.hpp"
#include"utils.hpp"
using namespace arma;
using namespace std;

int main(int argc, const char* argv[])

{
	int M = 201; //size of box
	double time_step = 2.5e-5;
	double T = 0.008;

	double x_c = 0.25;
	double y_c = 0.5;
	double mom_x = 200.0;
	double mom_y = 0.0;
	double s_x = 0.05;
	double s_y = 0.05;
  GaussParams gpars = {x_c,y_c, s_x, s_y, mom_x, mom_y};
	Box b1 = Box(M, gpars); //create box
	complex<double> r = time_step/(2*b1.h_*b1.h_)*1.0i;
  //initialize the potential :

	sp_mat V = sp_mat(M-2,M-2);

	//prepare matrices A and B
	b1.fill_a_and_b(r, time_step, V);
	b1.prepare_matrices(r);

	int Ntot  = ceil(T/time_step); //initialize total number of iterations
	cout << "total n. iterations" << Ntot << endl;
  cube Ufinal = cube(M, M, Ntot, fill::zeros);
//std::cout <<  "DEBUG: " << __LINE__ << std::endl;
//cout << "entering for loop" << endl;
	double e_time = 0.0;
	mat probability;
 	for(int ii=0; ii<Ntot ; ii++){
	//	cout << "ii = " << ii << endl ;
		probability = real(b1.U_%conj(b1.U_));
	//	cout << "advanced" << endl ;
		Ufinal.slice(ii)=probability;
	  b1.advance();
		e_time+=time_step;
	}

Ufinal.save("schrod.bin", arma_binary);

	return 0;
}
