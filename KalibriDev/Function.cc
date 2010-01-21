#include "Function.h"


double Function::parCov(int i, int j) const {
  assert( i >= 0 && i < nPars() );
  assert( j >= 0 && j < nPars() );

  int idx = (i*i + i)/2 + j;
  return cov()[covIdx_[idx]];
}
