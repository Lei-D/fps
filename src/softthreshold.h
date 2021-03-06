//
// softthreshold.h
// 
// Copyright 2014 Vincent Q. Vu. All rights reserved
//

#ifndef __SOFTTHRESHOLD_H
#define __SOFTTHRESHOLD_H

#include <algorithm>

// Proximal operator for \lambda |x|_1
struct SoftThresholdOp
{
  SoftThresholdOp(const double& z) : z(z) {}
  const double operator() (const double& x) const {
    return ((x > 0) - (x < 0)) * std::max(0.0, std::abs(x) - z);
  }

private:
    const double z;
};

// Proximal operator for \lambda (\alpha |x|_1 + 0.5 (1-\alpha) |x|_2^2)
struct ElasticSoftThresholdOp
{
  ElasticSoftThresholdOp(const double& z, const double& alpha) :
    z1(z * alpha), z2(1.0 / (1.0 + 0.5 * (1.0 - alpha) * z)) {}
  const double operator() (const double& x) const {
    return ((x > 0) - (x < 0)) * z2 * std::max(0.0, std::abs(x) - z1);
  }

private:
  const double z1, z2;
};

struct EntrywiseSoftThreshold
{
  EntrywiseSoftThreshold(const double& lambda) : lambda(lambda) {}
  void operator()(arma::mat& x, const double& z) const {
    x.transform( SoftThresholdOp(z * lambda) );
  }
private:
  const double lambda;
};

#endif
