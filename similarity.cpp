#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
double cosine_similarity(arma::vec A, arma::vec B)
{
	double theta = 0;
	theta = dot(A, B)/(norm(A)*norm(B));
	return theta;
}
/*This function checks the row A against all rows of matrix X.
  Matrix X should be subsetted s.t all obs belong to same clust.
  Should also be that the user has not visited the place before*/
// [[Rcpp::export]]
int best_cosine_row(arma::mat X, arma::vec A, int n)
{
	int row = 0;
	int it_end = n;
	double sim = 1;
	for(int i = 0; i != it_end; i++)
	{
		double current_sim = cosine_similarity(X.row(i), A);
		if(current_sim < sim) {
			sim = current_sim;
			row = i;
		}
	}
	return(row);
}
