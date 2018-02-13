#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

double cosine_similarity(arma::vec A, arma::vec B)
{
	double theta = 0;
	theta = dot(A, B)/(norm(A)*norm(B));
	return theta;
}

// [[Rcpp::export]]
NumericMatrix similar_users_rcpp(NumericMatrix X, NumericVector y, int n)
{
	NumericMatrix res(n, 2);
	for(int i = 0; i != n; i++)
	{
		res(i, 1) = cosine_similarity(y, X.row(i));
		res(i, 0) = i;
	}

	return(res);
}


/*using functions written by hand*/
/*This cant be written as cppfunction in the R-script. Could write dot() and norm() by hand.
  Default for norm is p = 2*/
double dot_man(NumericVector A, NumericVector B)
{
	int n = A.size();
	double ret = 0;

	for(int i = 0; i != n; i++){
		ret += A(i) * B(i);
	}

	return(ret);
}
/*should be sum of absolute values, but values are 0,1 so doesn't matter*/
double norm_man(NumericVector A, double p = 2)
{
	int n = A.size();
	double ret = 0;

	for(int i = 0; i != n; i++){
		ret =+ pow(A(i), 2);
	}

	return(pow(ret, 1/p));
}

double cosine_similarity_man(NumericVector A, NumericVector B)
{
	double theta = 0;
	theta = dot_man(A, B)/(norm_man(A)*norm_man(B));
	return theta;
}

// [[Rcpp::export]]
NumericMatrix similar_users_rcpp_man(NumericMatrix X, NumericVector y, int n)
{
	NumericMatrix res(n, 2);
	for(int i = 0; i != n; i++)
	{
		res(i, 1) = cosine_similarity_man(y, X.row(i));
		res(i, 0) = i;
	}

	return(res);
}


