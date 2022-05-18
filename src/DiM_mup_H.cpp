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
#include "DiM_mup_H.h"
#include "common.h"

#include "params.h"
#include "libcomplexop.h"

namespace c9_nmfv {

complex_t DiM_mup_H(
        param_t const &param
        )
{
    clearcache();
    const real_t M_W = param.M_W;
    const real_t beta = param.beta;
    const real_t e_em = param.e_em;
    const real_t m_Hp = param.m_Hp;
    const real_t m_mu = param.m_mu;
    const real_t s_12 = param.s_12;
    const real_t theta_W = param.theta_W;
    const complex_t IT_0000 = powq(e_em, -1);
    const complex_t IT_0001 = powq(M_W, -1);
    const complex_t IT_0002 = cosq(beta);
    const complex_t IT_0003 = cpowq(IT_0002, -1);
    const complex_t IT_0004 = sinq(beta);
    const complex_t IT_0005 = sinq(theta_W);
    const complex_t IT_0006 = cpowq(IT_0005, -1);
    const complex_t IT_0007 = (complex_t{0, 1.4142135623731})*e_em*m_mu
      *IT_0001*IT_0003*IT_0004*IT_0006;
    const complex_t IT_0008 = 0.5*IT_0007;
    const complex_t IT_0009 = cpowq(IT_0008, 2);
    const complex_t IT_0010 = 0.101321183642338*IT_0009;
    const complex_t IT_0011 = cosq(theta_W);
    const complex_t IT_0012 = cpowq(IT_0011, -1);
    const complex_t IT_0013 = tanq(theta_W);
    const complex_t IT_0014 = cpowq(IT_0013, 2);
    const complex_t IT_0015 = cpowq(1 + IT_0014, (-0.5));
    const complex_t IT_0016 = (complex_t{0, 1})*e_em*IT_0012*IT_0015;
    const complex_t IT_0017 = powq(m_mu, 2);
    const complex_t IT_0018 = powq(m_Hp, 2);
    const complex_t IT_0019 = mty::lt::C0iC(3, IT_0017, 2*IT_0017 + (-2)*s_12,
       IT_0017, 0, IT_0018, IT_0018, mty::lt::reg_int);
    const complex_t IT_0020 = IT_0016*IT_0019;
    const complex_t IT_0021 = m_mu*IT_0020;
    const complex_t IT_0022 = 2*IT_0016;
    const complex_t IT_0023 = mty::lt::C0iC(12, IT_0017, 2*IT_0017 + (-2)
      *s_12, IT_0017, 0, IT_0018, IT_0018, mty::lt::reg_int);
    const complex_t IT_0024 = IT_0022*IT_0023;
    const complex_t IT_0025 = m_mu*IT_0024;
    const complex_t IT_0026 = mty::lt::C0iC(6, IT_0017, 2*IT_0017 + (-2)*s_12,
       IT_0017, 0, IT_0018, IT_0018, mty::lt::reg_int);
    const complex_t IT_0027 = IT_0016*IT_0026;
    const complex_t IT_0028 = m_mu*IT_0027;
    const complex_t IT_0029 = mty::lt::C0iC(15, IT_0017, 2*IT_0017 + (-2)
      *s_12, IT_0017, 0, IT_0018, IT_0018, mty::lt::reg_int);
    const complex_t IT_0030 = IT_0022*IT_0029;
    const complex_t IT_0031 = m_mu*IT_0030;
    const complex_t IT_0032 = IT_0019*IT_0022;
    const complex_t IT_0033 = m_mu*IT_0032;
    const complex_t IT_0034 = IT_0022*IT_0026;
    const complex_t IT_0035 = m_mu*IT_0034;
    const complex_t IT_0036 = mty::lt::C0iC(18, IT_0017, 2*IT_0017 + (-2)
      *s_12, IT_0017, 0, IT_0018, IT_0018, mty::lt::reg_int);
    const complex_t IT_0037 = IT_0022*IT_0036;
    const complex_t IT_0038 = m_mu*IT_0037;
    const complex_t IT_0039 = (complex_t{0, (-0.03125)})*m_mu*IT_0000*IT_0010*
      (IT_0021 + IT_0025 + IT_0028 + IT_0031) + (complex_t{0, (-0.03125)})*m_mu
      *IT_0000*IT_0010*(IT_0021 + IT_0025 + -IT_0028 + -IT_0031) + (complex_t{0,
       0.03125})*m_mu*IT_0000*IT_0010*(IT_0021 + -IT_0028 + -IT_0031 + -IT_0033 
      + IT_0035 + IT_0038) + (complex_t{0, 0.03125})*m_mu*IT_0000*IT_0010*
      (IT_0021 + IT_0028 + -IT_0031 + -IT_0033 + -IT_0035 + -IT_0038);
    return IT_0039;
}
} // End of namespace c9_nmfv
