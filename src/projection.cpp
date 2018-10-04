//
// projection.cpp
// 
// Created by Vincent Q. Vu on 2014-07-03
// Copyright 2014 Vincent Q. Vu. All rights reserved
// 
#include "projection.h"
#include "simplex.h"
#include "IRLB.h"
#include "Rcpp.h"

using namespace Rcpp;
using namespace arma;

void IrlbaProjection::operator()(mat& x) const {
  
  int rank;
  int active = x.n_cols + 1;
  int n = std::min(ncomp, active);  // data type inside std::min must be the same
  Rcpp::List decomp(5);

  decomp = IRLB(x, std::max(n,1), std::max(n,1)+7);
  // cout << "ncomp";
  // cout << std::max(n,1);
  
  arma::vec s = decomp["d"];
  arma::mat u = decomp["u"];
  arma::mat v = decomp["v"];
  rank = simplex(s, d);
  
  // cout << "eigenvalues";
  // cout << s.subvec(0, rank-1);
  // cout << "rank";
  // cout << rank;
  
  // Reconstruct
  x = (
      u.cols(0, rank - 1) *
      diagmat(s.subvec(0, rank - 1)) *
      v.cols(0, rank - 1).t()
  );
  
  return;
}

void FantopeProjection::operator()(mat& x) const {

  int rank;
  vec eigval;
  mat eigvec;

  eig_sym(eigval, eigvec, x);
  rank = simplex(eigval, d);

  // Reconstruct
  x = (
    eigvec.cols(eigvec.n_cols - rank, eigvec.n_cols - 1) *
      diagmat(eigval.subvec(eigval.n_elem - rank, eigval.n_elem - 1)) *
      eigvec.cols(eigvec.n_cols - rank, eigvec.n_cols - 1).t()
  );
  
  return;
}

void SingularValueProjection::operator()(mat& x) const {  

  int rank;
  vec s;
  mat u, v;

  svd(u, s, v, x);
  rank = simplex(s, d, true);

  // Reconstruct
  x = (
    u.cols(0, rank - 1) * 
    diagmat(s.subvec(0, rank - 1)) *
    v.cols(0, rank - 1).t()
  );

  return;
}
