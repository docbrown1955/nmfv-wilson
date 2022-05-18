// This file is part of NMFV_WILSONS.
//
// NMFV_WILSONS is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// NMFV_WILSONS is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with NMFV_WILSONS. If not, see <https://www.gnu.org/licenses/>.

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
