#include"Box.hpp"
// Constructor of box
Box::Box(int M, GaussParams gp)
{
		M_ = M; //dimension of box including boundary
		h_ = 1.0/(M_-1);//space step
//std::cout <<  "DEBUG: " << __LINE__ << std::endl;
	  A_ = arma::sp_cx_mat((M_-2)*(M_-2), (M_-2)*(M_-2));
		//std::cout <<  "DEBUG: " error: arma::memory::acquire(): out of memoryerror: arma::memory::acquire(): out of memoryerror: arma::memory::acquire(): out of memoryerror: arma::memory::acquire(): out of memoryerror: arma::memory::acquire(): out of memory<< __LINE__ << std::endl;
		//invA_ = arma::sp_cx_mat((M_-2)*(M_-2), (M_-2)*(M_-2));
		//std::cout <<  "DEBUG: " << __LINE__ << std::endl;

		B_ = arma::sp_cx_mat((M_-2)*(M_-2), (M_-2)*(M_-2));
		std::cout <<  "DEBUG: " << __LINE__ << std::endl;
		U_ = arma::cx_mat(M_, M_, arma::fill::zeros);//matrix containing the wave function
		std::cout <<  "DEBUG: " << __LINE__ << std::endl;
		//initialize matrix using gaussian wave packet:

		for (int ii = 1; ii<M_-1; ii++)
			{
				for (int jj = 1; jj<M_-1; jj++)
				{
					//std::cout <<  "DEBUG: " << __LINE__ << std::endl;
					U_(ii, jj) = gaussian_wave_packet(ii*h_, jj*h_, gp);
				}
			}
    normalizer(U_);
		a_ = arma::cx_vec((M_-2)*(M_-2));
		b_ = arma::cx_vec((M_-2)*(M_-2));
	  //Uinternal_ = arma::cx_mat(M_-2, M_-2);
	  b_tmp_ = arma::cx_vec((M_-2)*(M_-2) );
		u_vec_new_ = arma::cx_mat((M_-2)*(M_-2),1);
}

void Box::check_norm(arma::cx_mat u) {
	double s=0;
	arma::cx_mat u2 = u % arma::conj(u); //elementwise multiplication
	for (int i=1; i<M_-1; i++) {
		for (int j=1; j<M_-1; j++) {
			s = s + real(u2(i,j));
		}
	}
	double deviation = abs(s-1.0);
	if (deviation>0.01) {
		std::cout << "Not normalized" << std::endl;
	}
}

//this function prepares matrices A and B to run the Crank Nicholson algorithm
void Box::prepare_matrices(std::complex<double> r)
{
	int size_sub = M_-2;

	for (int i1=0; i1<size_sub; i1++) //loop over diagonal submatrices
	{
		// Filling diagonal submatrices
		for (int i2=0; i2<size_sub; i2++) //this for loop fills the diagonal elements
						{
							//std::cout <<  "DEBUG: " << __LINE__ << std::endl;
							A_(i1*size_sub + i2, i1*size_sub + i2) = a_(i1*size_sub + i2);
							//std::cout <<  "DEBUG: " << __LINE__ << std::endl;
							B_(i1*size_sub + i2, i1*size_sub + i2) = b_(i1*size_sub + i2);
						}
	  for (int i2=0; i2<size_sub-1; i2++) //fill sub and super diagonal elements of diagonal submatrix
						{
							//std::cout <<  "DEBUG: " << __LINE__ << std::endl;
							A_(i1*size_sub + i2+1, i1*size_sub + i2) = -r;
							//std::cout <<  "DEBUG: " << __LINE__ << std::endl;
							A_(i1*size_sub + i2, i1*size_sub + i2+1) = -r;
							//std::cout <<  "DEBUG: " << __LINE__ << std::endl;
							B_(i1*size_sub + i2+1, i1*size_sub + i2) = r;
							//std::cout <<  "DEBUG: " << __LINE__ << std::endl;
							B_(i1*size_sub + i2, i1*size_sub + i2+1) = r;
						}
	}

	for (int i1=0; i1<size_sub-1; i1++) //loop over"Sub- and superdiagonal submatrices"
		{
			for (int i2=0; i2<size_sub; i2++)
			 			{ //filling diagonal elements of sub and superdiagonal submatrices
							//std::cout <<  "DEBUG: " << __LINE__ << std::endl;
							A_((i1+1)*size_sub + i2, i1*size_sub + i2) = -r;
							//std::cout <<  "DEBUG: " << __LINE__ << std::endl;
							A_(i1*size_sub + i2, (i1+1)*size_sub + i2) = -r;
							//std::cout <<  "DEBUG: " << __LINE__ << std::endl;
							B_((i1+1)*size_sub + i2, i1*size_sub + i2) = r;
							//std::cout <<  "DEBUG: " << __LINE__ << std::endl;
							B_(i1*size_sub + i2, (i1+1)*size_sub + i2) = r;

						}

		}

		//invA_ = arma::inv(A_);
}
void Box::fill_a_and_b(std::complex<double> r, double time_step, arma::sp_mat V)
{
 arma::cx_double z1 = time_step*0.5j;
 //std::cout <<  "DEBUG: " << __LINE__ << std::endl;
	for (int ii=0; ii<(M_-2); ii++)
	{
		for (int jj=0; jj<M_-2; jj++)
		{
		//std::cout <<  "DEBUG: " << __LINE__ << std::endl;
			a_(sub2ind(ii,jj, M_-2)) = 1.0 + 4.0*r + z1*  (   (double)V(ii,jj) ) ;
			b_(sub2ind(ii,jj, M_-2)) = 1.0 - 4.0*r - z1*  (   (double)V(ii,jj) ) ;
	 	}
	}
}// 



void Box::advance() {
//	std::cout <<  "DEBUG: " << __LINE__ << std::endl;
  //Uinternal_ = ;
//	std::cout <<  "DEBUG: " << __LINE__ << std::endl;
	b_tmp_ = B_*(U_(arma::span(1,M_-2), arma::span(1,M_-2)).as_col());
//	std::cout <<  "DEBUG: " << __LINE__ << std::endl
	//this is the correct line!!!!
	arma::spsolve(u_vec_new_, A_, b_tmp_);
//  u_vec_new_ = invA_ * b_tmp_;
//	std::cout << u_vec_new.size() << std::endl;
//	std::cout <<  "DEBUG: " << __LINE__ << std::endl;
  u_vec_new_.reshape(M_-2,M_-2);
	U_(arma::span(1,M_-2), arma::span(1,M_-2)) = u_vec_new_;
	//check_norm(u_new);
}
