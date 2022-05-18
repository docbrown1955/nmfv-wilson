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

#include "ckm.h"

using namespace c9_nmfv;

void fillCKM(c9_nmfv::param_t &params)
{
    params.V_ud = real_t{0.97434};
    params.V_us = real_t{0.22506};
    auto V_ub = 0.00117;// + i*-0.00338;
    params.V_ub_mod = cabsq(V_ub);
    params.delta_wolf = -cargq(V_ub);
    params.V_cd = complex_t{-0.22492 , -0.00014};
    params.V_cs = complex_t{0.97351 , -0.00003};
    params.V_cb = real_t{0.04108};
    params.V_td = complex_t{0.00810 , -0.00329};
    params.V_ts = complex_t{-0.04029 , -0.00076};
    params.V_tb = real_t{0.99915};
}

