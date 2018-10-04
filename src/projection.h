//
// projection.h
// 
// Copyright 2014 Vincent Q. Vu. All rights reserved
//

#ifndef __PROJECTION_H
#define __PROJECTION_H

#include <RcppArmadillo.h>

struct IrlbaProjection {
  
  IrlbaProjection(double d, int ncomp) : d(d), ncomp(ncomp) {}
  void operator()(arma::mat& x) const;
  
private:
  double d;
  int ncomp;
};

struct FantopeProjection {

  FantopeProjection(double d) : d(d) {}
  void operator()(arma::mat& x) const;

private:
  double d;
};

struct SingularValueProjection {

  SingularValueProjection(double d) : d(d) {}
  void operator()(arma::mat& x) const;

private:
  double d;
};

#endif
