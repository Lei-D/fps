//
// irlb.h
// 
// Copyright  All rights reserved
//

#ifndef __IRLB_H
#define __IRLB_H

#include <RcppArmadillo.h>

Rcpp::List IRLB(arma::mat& X,
          int nu,
          int work,
          int maxit=1000,
          double tol=1e-5,
          double eps=1e-9,
          double svtol=1e-5);
  
#endif
