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

#include <quadmath.h>
#include "clooptools.h"
#include "marty/looptools_init.h"
#include "marty/looptools_quad_extension.h"
#include <cmath>
#include "marty/looptools_interface.h"
#include "m_Z.h"
#include "common.h"

#include "params.h"
#include "libcomplexop.h"

namespace c9_nmfv {

complex_t m_Z(
        param_t const &param
        )
{
    clearcache();
    const real_t M_W = param.M_W;
    const real_t theta_W = param.theta_W;
    const complex_t IT_0000 = tanq(theta_W);
    const complex_t IT_0001 = cpowq(IT_0000, 2);
    const complex_t IT_0002 = cpowq(1 + IT_0001, 0.5);
    const complex_t IT_0003 = M_W*IT_0002;
    return IT_0003;
}
} // End of namespace c9_nmfv
