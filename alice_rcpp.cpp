#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

/*guess we need a function to extract submatricies.
  For the delta_Uik & delta_Vjk*/
arma::vec sub_Mat(arma::mat A, arma::vec IDX, int j,
	bool user = true)
{
	int n = IDX.size();
	arma::vec return_vec(n);

	if(user == true)
	{
	arma::vec x = A.col(j);

	for(int i = 0; i != n; i++){
		/*everything starts one position back*/
		return_vec(i) = x(IDX(i)-1);
		}
	} else if(user != true)
	{
		arma::rowvec x = A.row(j);

		for(int i = 0; i != (n-1); i++){
			return_vec(i) = x(IDX(i)-1);
		}

	}

	return(return_vec);
}

/* Need to deal with updating U by helper function
   U[userIDX[[j]], f] <- U[userIDX[[j]], f] + delta_Uik*/
arma::mat reader(arma::mat U, arma::vec idx, int f, 
	arma::vec U_sub, arma::vec delta_Uik)
{
	int n = idx.size();
	for(int i = 0; i != n; i++){
		U(idx(i) - 1, f) = U_sub(i) + delta_Uik(i);
	}

	return(U);
}

/*NA remover for following line of code
  new_error <- sqrt(sum(abs(x - p)^2, na.rm = TRUE)/length(x))
  this function should sum square of absolute values of all non NA's
  of a matrix. Take a matrix which is the difference of two!*/

double error_update(arma::mat X, arma::mat p)
{
	int row = X.n_rows;
	int col = X.n_cols;

	arma::mat diff(row, col);
	diff = X - p;

	double ret = 0;

	for(int i = 0; i != row; i++)
	{
		for(int j = 0; j != col; j++) 
		{
			if(!R_IsNA(X(i, j))) {
				ret += pow(X(i, j) - p(i, j), 2);
			} else {
				continue;
			}
		}
	}
	return(pow(ret / (col*row), 0.5));
}


// [[Rcpp::export]]
List alice_fast(arma::mat X, arma::mat U, arma::mat V, 
	List itemIDX, List userIDX, int k, double lambda, double gamma,
	double min_improvement, int max_epochs, int min_epochs)
{
	double errors_length = 1.0e+5;
	arma::vec errors(errors_length);

	/*go through all features*/
	for(int f = 0; f != k; f++)
	{
		int col = X.n_cols;
		int row = X.n_rows;
		/*convergence check. Can use arma::math::inf() instead
		but gives a cryptic error*/
		double last_error =  1.0e+10; 
		double delta_error =  1.0e+10;
		int epoch = 0 ;

		arma::mat p(row, col);
		p = U * arma::trans(V);
		arma::mat error(row, col);

		while(epoch < max_epochs && delta_error > min_improvement)
		{
			error = X - p;
			arma::mat temp_U = U;

			/*user features*/
			for(int j = 0; j != col; j++)
			{
				arma::vec u = as<NumericVector>(userIDX[j]);
				int n = u.size();
				arma::vec delta_Uik(n);
				arma::vec error_temp(n);
				error_temp = sub_Mat(error, u, j);
				arma::vec U_j(n);
				U_j = sub_Mat(U, u, f);


				delta_Uik = lambda * (error_temp * V(j, f) - gamma * U_j);
				U = reader(U, u, f, U_j, delta_Uik);

			}
			/*item features*/
			for(int i = 0; i != row; i++)
			{
				arma::vec v = as<NumericVector>(itemIDX[i]);
				int n = v.size();
				arma::vec delta_Vik(n);
				arma::vec error_temp(n);
				error_temp = sub_Mat(error, v, i, false);
				arma::vec V_i(n);
				V_i = sub_Mat(V, v, f);

				delta_Vik = lambda * (error_temp * temp_U(i, f) - gamma * V_i);
				V = reader(V, v, f, V_i, delta_Vik);

			}
			/*update error*/
			p = U * arma::trans(V);
			double new_error = error_update(X, p);
			/*absolute value, to lazy to google*/
			delta_error = pow(pow(last_error - new_error, 2), 0.5);
			epoch += 1;
			errors(epoch) = new_error;

			last_error = new_error;
		
		}
	}
	return List::create(Named("U") = U,
						Named("V") = V,
						Named("errors") = errors);

}




