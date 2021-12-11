#include"Box.hpp"
// Constructor of box
//inputs: an integer M (dimension of the box) and a struct GaussParams which contains the
//Gaussian wave packet parameters that we want to set.
Box::Box(int M, GaussParams gp)
{
			M_ = M; //dimension of box including boundary
			h_ = 1.0/(M_-1);//space step
			//these are the matrices A and B which are necessary for the Crank Nicholson
			A_ = arma::sp_cx_mat((M_-2)*(M_-2), (M_-2)*(M_-2));
			B_ = arma::sp_cx_mat((M_-2)*(M_-2), (M_-2)*(M_-2));

			//matrix containing the complex wave function
			U_ = arma::cx_mat(M_, M_, arma::fill::zeros);
			//NOTE: THE BOUNDARY CONDITION IS ENSURED BECAUSE WE INITIALIZE THE FUNCTION U TO
			//ALL ZEROS . When initializing the gaussian wave packet, we do not loop on the boundary points.
			//Also, we will never update them.
			
			//initialize wave function matrix using gaussian wave packet:
			for (int ii = 1; ii<M_-1; ii++)
			{
					for (int jj = 1; jj<M_-1; jj++)
					{
						//this function is defined in "utils.cpp"
						U_(ii, jj) = gaussian_wave_packet(ii*h_, jj*h_, gp);
					}
			}
			//this function normalizes the wave packet. Defined in "utils.cpp"
			normalizer(U_);
			//these are the vectors a and b which will contain the diagonal elements of
			//matrices A and B.
			a_ = arma::cx_vec((M_-2)*(M_-2));
			b_ = arma::cx_vec((M_-2)*(M_-2));

			//these two variables come in handy when it's necessary to solve the system Au^{n+1}= Bu^n
			//Having them here avoids re-allocating memory every time the system evolves.
			b_tmp_ = arma::cx_vec((M_-2)*(M_-2) );
			u_vec_new_ = arma::cx_mat((M_-2)*(M_-2),1);
}

//this function checks that the total probability is "not too far" from 1.
//Might be a useful debugging tool.
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
//Before running this, it is necessary to have vectors a and b properly initialized.
//The initialization of a and b is done by the function fill_a_and_b (see below).
void Box::prepare_matrices(std::complex<double> r)
{
	int size_sub = M_-2; //size of the interesting part of the wave function matrix (excluding boundary)

	for (int i1=0; i1<size_sub; i1++) //loop over diagonal submatrices
	{
		// Filling diagonal submatrices
		for (int i2=0; i2<size_sub; i2++) //this for loop fills the diagonal elements
						{													//of each diagonal submatrix D_m^A and D_m^B (see notation in the report)
							A_(i1*size_sub + i2, i1*size_sub + i2) = a_(i1*size_sub + i2);
							B_(i1*size_sub + i2, i1*size_sub + i2) = b_(i1*size_sub + i2);
						}
	  for (int i2=0; i2<size_sub-1; i2++) //fills sub and super diagonal elements of diagonal submatrix
						{
							A_(i1*size_sub + i2+1, i1*size_sub + i2) = -r;
							A_(i1*size_sub + i2, i1*size_sub + i2+1) = -r;
							B_(i1*size_sub + i2+1, i1*size_sub + i2) = r;
							B_(i1*size_sub + i2, i1*size_sub + i2+1) = r;
						}
	}

	for (int i1=0; i1<size_sub-1; i1++) //loop over "sub-" and "super-" diagonal submatrices
		{
			for (int i2=0; i2<size_sub; i2++)
			 			{ //filling diagonal elements of each sub and superdiagonal submatrix

							A_((i1+1)*size_sub + i2, i1*size_sub + i2) = -r;
							A_(i1*size_sub + i2, (i1+1)*size_sub + i2) = -r;
							B_((i1+1)*size_sub + i2, i1*size_sub + i2) = r;
							B_(i1*size_sub + i2, (i1+1)*size_sub + i2) = r;

						}

		}

}//end of function "prepare_matrices"

//this function fills the vectors a and b, which contain the diagonal elements of the matrices A and B
void Box::fill_a_and_b(std::complex<double> r, double time_step, arma::mat V)
{
 arma::cx_double z1 = time_step*0.5j; //doing this before entering the for loop spares some 	FLOPs
	for (int ii=0; ii<(M_-2); ii++)
	{
		for (int jj=0; jj<M_-2; jj++)
		{
			//sub2ind properly converts matrix subscripts to linear indices. See utils.cpp for the definition.
			a_(sub2ind(ii,jj, M_-2)) = 1.0 + 4.0*r + z1*  (   (double)V(ii,jj) ) ;
			b_(sub2ind(ii,jj, M_-2)) = 1.0 - 4.0*r - z1*  (   (double)V(ii,jj) ) ;
	 	}
	}
}


//This function updates the state of the system, by solving the linear system A u^{n+1} = B u^n
void Box::advance() {
	//The built-in method .as_col is really useful because it converts our matrix U_ to a column vector
	//in the exact way we want.
	b_tmp_ = B_*(U_(arma::span(1,M_-2), arma::span(1,M_-2)).as_col());
	//Solve a sparse system of equations using library solver spsolve()
	arma::spsolve(u_vec_new_, A_, b_tmp_);
	//convert the resulting column vector back to matrix form
  u_vec_new_.reshape(M_-2,M_-2);
  //update current state of the system
	U_(arma::span(1,M_-2), arma::span(1,M_-2)) = u_vec_new_;
}
