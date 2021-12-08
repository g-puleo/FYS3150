#include"utils.hpp"
//this function calculates the complex value of the gaussian wave packet at the point
//x, y
const arma::cx_double i(0.0,1.0);
arma::cx_double gaussian_wave_packet( double x, double y, GaussParams pars )
    {

      double x1 = (x-pars.xc)*(x-pars.xc)/(2.0*pars.sigma_x*pars.sigma_x);
      double x2 = (y-pars.yc)*(y-pars.yc)/(2.0*pars.sigma_y*pars.sigma_y);
      double p1 = (pars.px*(x-pars.xc));
      double p2 = (pars.py*(y-pars.yc));
      return std::exp(-x1-x2+(p1+p2)*i);

    }


//this function normalizes the matrix u so that the sum of the absolute values (squared) of
//all of its elements equals 1.
void normalizer( arma::cx_mat& u )
  {
    double s = 0;
    int M = u.n_rows;
    std::cout << "normalizing wavefunction" << std::endl;
    arma::cx_mat conj_u = arma::conj(u);
    arma::cx_mat u2 = u % conj_u;
    for (int ii=0; ii<M; ii++)
    {
      for (int jj=0; jj<M; jj++)
     {
        s = s + real(u2(ii,jj));
      }
    }
    u = u/sqrt(s);
  }

//translate matrix indices (i,j) to vector subscript k, given the size of a matrix
int sub2ind(int ii, int jj, int size)
  {
  	int k = (ii) + (jj)*(size);
  	return k;
  }


std::vector<std::string> split (const std::string &s, char delim) {
    std::vector<std::string> result;
    std::stringstream ss (s);
    std::string item;

    while (getline (ss, item, delim)) {
        result.push_back (item);
    }

    return result;
}
  
//reading file with two lines
std::vector<std::string> read_file(const std::string filename, char delim) {
	std::fstream initial_conditions; 
	initial_conditions.open(filename, std::ios::in);
	std::string input; 
	std::string values; 
	if (!initial_conditions) {
		std::cout << "No such file";
	}
	else {
		getline(initial_conditions, input); 
		getline(initial_conditions, values); 
	}
	std::vector<std::string> vals = split(values, delim); 
	return vals; 

}



