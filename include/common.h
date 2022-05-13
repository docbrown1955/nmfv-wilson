#ifndef CSL_LIB_C9_NMFV_COMMON_H_INCLUDED
#define CSL_LIB_C9_NMFV_COMMON_H_INCLUDED
#include <complex>
#include <quadmath.h>
#define CSL_LT_DISABLE_ITERATOR
#include "librarytensor.h"

namespace c9_nmfv {

using real_t = __float128;
using complex_t = __complex128;
static constexpr complex_t _i_{0, 1};

void setMu(const double mu);
void setLambda2(const double lambda2);
void setUVDiv(const double x);

} // End of namespace c9_nmfv

#endif
