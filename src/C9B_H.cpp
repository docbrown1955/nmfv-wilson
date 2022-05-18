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
#include "C9B_H.h"
#include "common.h"

#include "params.h"
#include "libcomplexop.h"

namespace c9_nmfv {

complex_t C9B_H(
        param_t const &param
        )
{
    clearcache();
    const real_t M_W = param.M_W;
    const real_t m_b = param.m_b;
    const real_t m_c = param.m_c;
    const real_t m_s = param.m_s;
    const real_t m_t = param.m_t;
    const real_t m_u = param.m_u;
    const real_t V_cb = param.V_cb;
    const real_t V_tb = param.V_tb;
    const real_t V_us = param.V_us;
    const real_t beta = param.beta;
    const real_t e_em = param.e_em;
    const real_t m_Hp = param.m_Hp;
    const real_t m_mu = param.m_mu;
    const real_t theta_W = param.theta_W;
    const real_t V_ub_mod = param.V_ub_mod;
    const real_t delta_wolf = param.delta_wolf;
    const complex_t V_cs = param.V_cs;
    const complex_t V_ts = param.V_ts;
    const complex_t IT_0000 = powq(M_W, 2);
    const complex_t IT_0001 = powq(V_tb, -1);
    const complex_t IT_0002 = cpowq(conjq(V_ts), -1);
    const complex_t IT_0003 = powq(e_em, -4);
    const complex_t IT_0004 = 9.86960440108936*IT_0000*IT_0001*IT_0002*IT_0003;
    const complex_t IT_0005 = powq(M_W, -1);
    const complex_t IT_0006 = cexpq((complex_t{0, -1})*delta_wolf);
    const complex_t IT_0007 = cosq(beta);
    const complex_t IT_0008 = sinq(beta);
    const complex_t IT_0009 = cpowq(IT_0008, -1);
    const complex_t IT_0010 = sinq(theta_W);
    const complex_t IT_0011 = cpowq(IT_0010, -1);
    const complex_t IT_0012 = (complex_t{0, 1.4142135623731})*m_u*e_em*IT_0005
      *IT_0006*IT_0007*IT_0009*IT_0011*V_ub_mod;
    const complex_t IT_0013 = 0.5*IT_0012;
    const complex_t IT_0014 = (complex_t{0, 1.4142135623731})*m_u*V_us*e_em
      *IT_0005*IT_0007*IT_0009*IT_0011;
    const complex_t IT_0015 = 0.5*IT_0014;
    const complex_t IT_0016 = powq(m_u, 2);
    const complex_t IT_0017 = powq(m_Hp, 2);
    const complex_t IT_0018 = mty::lt::D0iC(12, 0, 0, 0, 0, 0, 0, 0, IT_0016,
       IT_0017, IT_0017, mty::lt::reg_int);
    const complex_t IT_0019 = IT_0013*IT_0015*IT_0018;
    const complex_t IT_0020 = cpowq(IT_0007, -1);
    const complex_t IT_0021 = (complex_t{0, 1.4142135623731})*e_em*m_mu
      *IT_0005*IT_0008*IT_0011*IT_0020;
    const complex_t IT_0022 = 0.5*IT_0021;
    const complex_t IT_0023 = cpowq(IT_0022, 2);
    const complex_t IT_0024 = (complex_t{0, 0.101321183642338})*IT_0023;
    const complex_t IT_0025 = IT_0019*IT_0024;
    const complex_t IT_0026 = (complex_t{0, 1.4142135623731})*m_c*V_cb*e_em
      *IT_0005*IT_0007*IT_0009*IT_0011;
    const complex_t IT_0027 = 0.5*IT_0026;
    const complex_t IT_0028 = (complex_t{0, 1.4142135623731})*m_c*conjq(V_cs)
      *e_em*IT_0005*IT_0007*IT_0009*IT_0011;
    const complex_t IT_0029 = 0.5*IT_0028;
    const complex_t IT_0030 = powq(m_c, 2);
    const complex_t IT_0031 = mty::lt::D0iC(12, 0, 0, 0, 0, 0, 0, 0, IT_0030,
       IT_0017, IT_0017, mty::lt::reg_int);
    const complex_t IT_0032 = IT_0027*IT_0029*IT_0031;
    const complex_t IT_0033 = IT_0024*IT_0032;
    const complex_t IT_0034 = (complex_t{0, 1.4142135623731})*m_t*V_tb*e_em
      *IT_0005*IT_0007*IT_0009*IT_0011;
    const complex_t IT_0035 = 0.5*IT_0034;
    const complex_t IT_0036 = (complex_t{0, 1.4142135623731})*m_t*conjq(V_ts)
      *e_em*IT_0005*IT_0007*IT_0009*IT_0011;
    const complex_t IT_0037 = 0.5*IT_0036;
    const complex_t IT_0038 = powq(m_t, 2);
    const complex_t IT_0039 = mty::lt::D0iC(12, 0, 0, 0, 0, 0, 0, 0, IT_0038,
       IT_0017, IT_0017, mty::lt::reg_int);
    const complex_t IT_0040 = IT_0035*IT_0037*IT_0039;
    const complex_t IT_0041 = IT_0024*IT_0040;
    const complex_t IT_0042 = (complex_t{0, 1.4142135623731})*m_b*V_cb*e_em
      *IT_0005*IT_0008*IT_0011*IT_0020;
    const complex_t IT_0043 = 0.5*IT_0042;
    const complex_t IT_0044 = (complex_t{0, 1.4142135623731})*m_s*conjq(V_cs)
      *e_em*IT_0005*IT_0008*IT_0011*IT_0020;
    const complex_t IT_0045 = 0.5*IT_0044;
    const complex_t IT_0046 = IT_0031*IT_0043*IT_0045;
    const complex_t IT_0047 = IT_0024*IT_0046;
    const complex_t IT_0048 = (complex_t{0, 1.4142135623731})*m_b*e_em*IT_0005
      *IT_0006*IT_0008*IT_0011*IT_0020*V_ub_mod;
    const complex_t IT_0049 = 0.5*IT_0048;
    const complex_t IT_0050 = (complex_t{0, 1.4142135623731})*m_s*V_us*e_em
      *IT_0005*IT_0008*IT_0011*IT_0020;
    const complex_t IT_0051 = 0.5*IT_0050;
    const complex_t IT_0052 = IT_0018*IT_0049*IT_0051;
    const complex_t IT_0053 = IT_0024*IT_0052;
    const complex_t IT_0054 = (complex_t{0, 1.4142135623731})*m_b*V_tb*e_em
      *IT_0005*IT_0008*IT_0011*IT_0020;
    const complex_t IT_0055 = 0.5*IT_0054;
    const complex_t IT_0056 = (complex_t{0, 1.4142135623731})*m_s*conjq(V_ts)
      *e_em*IT_0005*IT_0008*IT_0011*IT_0020;
    const complex_t IT_0057 = 0.5*IT_0056;
    const complex_t IT_0058 = IT_0039*IT_0055*IT_0057;
    const complex_t IT_0059 = IT_0024*IT_0058;
    const complex_t IT_0060 = cpowq(IT_0010, 2);
    const complex_t IT_0061 = (IT_0025 + IT_0033 + IT_0041 + IT_0047 + IT_0053
       + IT_0059)*IT_0060;
    const complex_t IT_0062 = IT_0004*IT_0061;
    const complex_t IT_0063 = (IT_0025 + IT_0033 + IT_0041 + -IT_0047 + 
      -IT_0053 + -IT_0059)*IT_0060;
    const complex_t IT_0064 = IT_0004*IT_0063;
    const complex_t IT_0065 = (-0.5)*IT_0064;
    return (complex_t{0, (-0.5)})*IT_0062 + (complex_t{0, 1})*IT_0065;
}
} // End of namespace c9_nmfv
