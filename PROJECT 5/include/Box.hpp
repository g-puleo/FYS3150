#ifndef __Box_hpp__
#define __Box_hpp__
#include<armadillo>
#include<cmath>
#include<complex>
#include"utils.hpp"

class Box
{

public:
//	arma::mat box_;
	int M_;   // dimension of the grid
	double h_;  // step size of the grid

	//these are useful matrices and vector to run the Crank Nicholson algorithm
	arma::sp_cx_mat A_;
	arma::sp_cx_mat B_;
	arma::cx_mat U_; //contains state of the system
	arma::cx_vec a_; //a and b are handy vectors during the CN algorithm
	arma::cx_vec b_;
	arma::cx_vec b_tmp_;
	arma::cx_mat u_vec_new_ ;
	//constructor
	Box(int M, GaussParams gp);
	//this will be the time evolution function
  	void prepare_matrices(std::complex<double> r);
  	void check_norm(arma::cx_mat u) ;
	void fill_a_and_b(std::complex<double> r, double time_step, arma::mat V);
	void advance(void);
};

#endif
