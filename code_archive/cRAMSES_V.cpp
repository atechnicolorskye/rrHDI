#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::plugins(cpp11)]]

// JLK: below function is taken from
// http://gallery.rcpp.org/articles/fast-linear-model-with-armadillo/
// [[Rcpp::export]]
Rcpp::List fastLm(const arma::vec & y, const arma::mat & X) {
  
  int n = X.n_rows, k = X.n_cols;
  
  arma::colvec coef = arma::solve(X, y); 
  arma::colvec resid = y - X*coef; 
  
  double sig2 = arma::as_scalar(arma::trans(resid)*resid/(n-k));
  arma::colvec stderrest = 
    arma::sqrt(sig2 * arma::diagvec( arma::inv(arma::trans(X)*X)) );
  
  return Rcpp::List::create(Rcpp::Named("coefficients") = coef,
                            Rcpp::Named("stderr")       = stderrest);
}


// [[Rcpp::export]]
arma::vec OLS_c(arma::vec y,
                arma::mat X) {
  
  arma::vec coef = arma::solve(X, y); 
  return coef;
}

// [[Rcpp::export]]
double OLS_c_1d(arma::vec y,
                arma::vec x) {
  
  arma::vec coef = arma::solve(x, y); 
  return coef[0];
}

// [[Rcpp::export]]
arma::mat restricted_OLS_c(arma::vec y,
                           arma::mat X,
                           arma::vec bhat,
                           arma::mat Q,
                           double c) {
  arma::mat Sx = inv(trans(X) * X);
  
  return bhat - Sx * Q * inv(trans(Q) * Sx * Q) * (trans(Q) * bhat - c);
}


// [[Rcpp::export]]
double Tn_c(arma::vec eps,
            arma::mat X,
            arma::vec lam) {
  
  arma::vec coef = fastLm(eps, X+0)[0]; // 0th index: coefficients, 1st index: stderr
  return dot(lam, coef);
}


// [[Rcpp::export]]
arma::vec g_c(Rcpp::List cluster_eps,
              bool use_perm,
              bool use_sign) { 
  // inputted partition will be a list of list, e.g. {{1, 2}, {3, 4, 5}, {6, 7}}
  // this is partitioned "er"
  int n_cluster = cluster_eps.size();
  // initilizing list with size n_cluster
  Rcpp::List permuted_cluster(n_cluster);
  
  // initializing vector that will be returned
  arma::vec out; //
  
  // IF there is only 1 cluster, just shuffle
  if (n_cluster == 1) {
    // Rprintf("entering here");
    arma::vec temp = cluster_eps[0];
    out = temp;
    if (use_perm) {
      out = arma::shuffle(temp);   
    }
  }
  
  // IF more than 1 cluster, permute and randomly change sign on cluster level
  else {
    // 1. samples random signs (as many clusters)
    arma::vec rand_int = arma::randi<arma::vec>(n_cluster, arma::distr_param(0, 1)); // samples 0 and 1 randomly
    rand_int.replace(0, -1); // replace 0 with -1
    // rand_int.print();
    
    // 2. permutes the elements within the specified clusters
    // and multiply by randomly generated sign on cluster level
    for (int i=0; i<n_cluster; ++i) {
      // 1. grab i^th cluster
      arma::vec temp = cluster_eps[i]; 
      
      // 2. randomly shuffle i^th cluster, multiply random sign, and append to out ("er")
      if (use_perm && use_sign) {
        out = arma::join_cols(out, arma::shuffle(temp) * rand_int[i]);    
      } 
      else if (use_perm && !use_sign) {
        out = arma::join_cols(out, arma::shuffle(temp));    
      }
      else if (!use_perm && use_sign) {
        out = arma::join_cols(out, temp * rand_int[i]);    
      }
    }
    // print(cluster);
    // out.print();
  }
  return out;
}



// [[Rcpp::export]]
Rcpp::List r_test_c(arma::vec y,
                    arma::mat X,
                    arma::vec lam,
                    double lam0,
                    Rcpp::List cluster_eps_r,
                    bool use_perm,
                    bool use_sign,
                    int num_R) {
  if (y.size() != X.n_rows) {
    Rcpp::Rcout << "Error: length of y and nrow of X not matching, returning Inf" << std::endl;
    return std::numeric_limits<double>::infinity();
  }
  
  if (lam.size() != X.n_cols) {
    Rcpp::Rcout << "Error: length of lam and ncol of X not matching, returning Inf" << std::endl;
    return std::numeric_limits<double>::infinity();
  }
  
  arma::vec bhat = fastLm(y, X+0)[0]; // 0th index: coefficients, 1st index: stderr
  double tobs = dot(lam, bhat) - lam0;
  
  // JLK: calculating "er" in example instead
  // arma::mat Q = lam;
  // arma::mat bhat_r = restricted_OLS_c(y, X, bhat, Q, lam0);
  // arma::mat er = y - X * bhat_r;
  // create clustered version of er = "restricted residuals".
  // Rcpp::List cluster_er = Rcpp::List::create();
  // for(int j=0; j < clustering.size(); ++j) {
  //   arma::vec ind = clustering[j]; // indexes that correspond to cluster j (e.g., 11, 12, 13...)
  //   arma::vec temp = arma::zeros(ind.size());
  //   for(int k = 0; k < temp.size(); ++k) {
  //     temp[k] = er[ind[k]]; //populate with residuals.
  //   }
  //   cluster_er.push_back(temp);
  // }
  
  arma::vec tvals = arma::zeros(num_R + 1);
  for (int i=0; i<num_R; ++i) {
    // JLK: given partitioned er (list of list), performs operations and outputs new er (vector)
    arma::vec er_new = g_c(cluster_eps_r, use_perm, use_sign);
    // Rcpp::print(er_new_partitioned);
    tvals(i) = Tn_c(er_new, X, lam);
  }
  tvals(num_R) = tobs;
  // Return test statistic values.
  return Rcpp::List::create(Rcpp::Named("tobs") = tobs,
                            Rcpp::Named("tvals") = tvals);
  
  // arma::vec pvalavgs = arma::zeros(2);
  // pvalavgs(0) = mean(bigger);
  // pvalavgs(1) = mean(smaller);
  // double pval = pvalavgs.min();
  // // 
  // return(pval);
}

// [[Rcpp::export]]
Rcpp::List r_test_exact(arma::vec y,
                   arma::vec x,
                   double beta0, double beta1,
                   int num_R) {
  // y = b1 * x + e  with Ho: b1=lam0
  // Performs inference based on exchangeability.
  //
  if (y.size() != x.size()) {
    Rcpp::Rcout << "Error: length of y and nrow of x not matching, returning Inf" << std::endl;
    return std::numeric_limits<double>::infinity();
  }
  
  arma::vec e = y - beta0 -beta1 * x; // residuals
  double tobs = OLS_c_1d(e, x); // test statistic.
    
  arma::vec tvals = arma::zeros(num_R + 1);
  for (int i=0; i<num_R; ++i) {
    arma::vec enew = arma::shuffle(e);  
    tvals(i) = OLS_c_1d(enew, x);
  }
  tvals(num_R) = tobs;
  // Return test statistic values.
  return Rcpp::List::create(Rcpp::Named("tobs") = tobs,
                            Rcpp::Named("tvals") = tvals);
  
}

// [[Rcpp::export]]
Rcpp::List r_test_composite_c(arma::vec y,
                            arma::vec x,
                            double beta1_H0,
                            arma::vec beta0_seq,
                            int num_R) {
  // y = b1 * x + e  with Ho: b1=lam0
  // Performs inference based on exchangeability.
  //
  if (y.size() != x.size()) {
    Rcpp::Rcout << "Error: length of y and nrow of x not matching, returning Inf" << std::endl;
    return std::numeric_limits<double>::infinity();
  }
  
  arma::vec y2 = y - arma::mean(y) - beta1_H0 * (x-arma::mean(x));

  int n = y2.size();
  arma::vec pvals = arma::zeros(beta0_seq.size());
  double tobs = OLS_c_1d(y2, x); // observed test statistic.
  
  for(int i_pval=0; i_pval<beta0_seq.size(); ++i_pval) {
    // calculate mean(eps) 
    double mean_eps = arma::mean(y) - beta0_seq(i_pval) - beta1_H0* arma::mean(x);
    // New test for specific beta0.
    //
    arma::vec tvals = arma::zeros(num_R);
    double bigger = 0;
    //
    for (int j=0; j<num_R; ++j) {
      arma::vec signs = arma::randi<arma::vec>(n, arma::distr_param(0, 1)); // samples 0 and 1 randomly
      signs.replace(0, -1); // replace 0 with -1
      // new dataset
      arma::vec y2_new = arma::zeros(n);
      for(int k=0; k<n; ++k) {
        y2_new(k) = signs(k) * y2(k) + signs(k) * mean_eps; //signs * mean_eps;
      }
      
      tvals(j) = OLS_c_1d(y2_new, x);
      if(tvals(j) >= tobs) {
        bigger = bigger+1;
      }
    }
    double m = bigger/num_R;
    
    pvals(i_pval) = std::min(m, 1-m);
    
  }
  
  // Return test statistic values.
  return Rcpp::List::create(Rcpp::Named("pvals") = pvals);
}


