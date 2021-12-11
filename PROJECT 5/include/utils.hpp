#ifndef __utils_hpp__
#define __utils_hpp__
#include<complex>
#include<armadillo>
#include<string>
#include<vector>
#include<iostream> 
#include<sstream>
#include<fstream>

//this struct  type is designed to contain parameters of a gaussian wave packet
struct gaussParams {
double xc;
double yc;
double sigma_x;
double sigma_y;
double px;
double py;
};

typedef struct gaussParams GaussParams ;

//this function calculates the complex value of the gaussian wave packet
arma::cx_double gaussian_wave_packet( double x, double y, GaussParams pars );

//this function normalizes a complex matrix
void normalizer( arma::cx_mat& wavefunction );

//convert matrix subscripts to linear indices
int sub2ind(int ii, int jj, int size);

//returns a vector with elements of the second line of a file
std::vector<std::string> read_file(const std::string filename, char delim); 

//splits a string using a delimiter
std::vector<std::string> split(const std::string &s, char delim); 


#endif
