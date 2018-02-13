#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

/*calculate the support (# non skips) of rest or user*/
// [[Rcpp::export]]
IntegerVector support(CharacterVector uniques, CharacterVector dat_vec)
{
	int n = uniques.size();
	int n_dat = dat_vec.size();
	IntegerVector vec(n);

	for(int i = 0; i < n; i++){
		for(int j = 0; j < n_dat; j++){
			if(uniques(i) == dat_vec(j)){
				vec(i) += 1;
			}
		}
	}
	return(vec);
}

/*returns the residuals from r_ui = theta_u*x_ui + resid.
  This function is written assuming that no user kan rate
  a restaurant twice*/ 
// [[Rcpp::export]]
NumericVector resid_rating(IntegerVector score, NumericVector x_u, double mean_x = 0)
{
	int n = score.size(); 
	double sum_x_u = 0;
	double theta = 0;
	NumericVector res(n);

	for(int i = 0; i != n; ++i){
		sum_x_u += pow(x_u(i), 2);
	}

	for(int j = 0; j != n; ++j){
		theta += (score(j)*x_u(j))/sum_x_u;
	}

	for(int i = 0; i != n; ++i){
		res(i) = score(i) - theta*(x_u(i)-mean_x);
	}

	return(res);
}
