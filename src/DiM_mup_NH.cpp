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
#include "DiM_mup_NH.h"
#include "common.h"

#include "params.h"
#include "libcomplexop.h"

namespace c9_nmfv {

complex_t DiM_mup_NH(
        param_t const &param
        )
{
    clearcache();
    const real_t M_W = param.M_W;
    const real_t beta = param.beta;
    const real_t e_em = param.e_em;
    const real_t m_A0 = param.m_A0;
    const real_t m_H0 = param.m_H0;
    const real_t m_mu = param.m_mu;
    const real_t s_12 = param.s_12;
    const real_t alpha = param.alpha;
    const real_t theta_W = param.theta_W;
    const complex_t IT_0000 = powq(e_em, -1);
    const complex_t IT_0001 = m_mu*IT_0000;
    const complex_t IT_0002 = cosq(theta_W);
    const complex_t IT_0003 = cpowq(IT_0002, -1);
    const complex_t IT_0004 = tanq(theta_W);
    const complex_t IT_0005 = cpowq(IT_0004, 2);
    const complex_t IT_0006 = cpowq(1 + IT_0005, (-0.5));
    const complex_t IT_0007 = (complex_t{0, 1})*e_em*IT_0003*IT_0006;
    const complex_t IT_0008 = -IT_0007;
    const complex_t IT_0009 = powq(M_W, -1);
    const complex_t IT_0010 = cosq(beta);
    const complex_t IT_0011 = cpowq(IT_0010, -1);
    const complex_t IT_0012 = sinq(beta);
    const complex_t IT_0013 = sinq(theta_W);
    const complex_t IT_0014 = cpowq(IT_0013, -1);
    const complex_t IT_0015 = e_em*m_mu*IT_0009*IT_0011*IT_0012*IT_0014;
    const complex_t IT_0016 = (-0.5)*IT_0015;
    const complex_t IT_0017 = cpowq(IT_0016, 2);
    const complex_t IT_0018 = 0.101321183642338*IT_0017;
    const complex_t IT_0019 = IT_0008*IT_0018;
    const complex_t IT_0020 = powq(m_mu, 2);
    const complex_t IT_0021 = powq(m_A0, 2);
    const complex_t IT_0022 = mty::lt::C0iC(0, IT_0020, IT_0020, 2*IT_0020 + (
      -2)*s_12, IT_0020, IT_0021, IT_0020, mty::lt::reg_int);
    const complex_t IT_0023 = m_mu*IT_0022;
    const complex_t IT_0024 = mty::lt::C0iC(6, IT_0020, IT_0020, 2*IT_0020 + (
      -2)*s_12, IT_0020, IT_0021, IT_0020, mty::lt::reg_int);
    const complex_t IT_0025 = m_mu*IT_0024;
    const complex_t IT_0026 = mty::lt::C0iC(12, IT_0020, IT_0020, 2*IT_0020 + 
      (-2)*s_12, IT_0020, IT_0021, IT_0020, mty::lt::reg_int);
    const complex_t IT_0027 = m_mu*IT_0026;
    const complex_t IT_0028 = mty::lt::C0iC(15, IT_0020, IT_0020, 2*IT_0020 + 
      (-2)*s_12, IT_0020, IT_0021, IT_0020, mty::lt::reg_int);
    const complex_t IT_0029 = m_mu*IT_0028;
    const complex_t IT_0030 = IT_0023 + IT_0025 + IT_0027 + IT_0029;
    const complex_t IT_0031 = mty::lt::C0iC(3, IT_0020, IT_0020, 2*IT_0020 + (
      -2)*s_12, IT_0020, IT_0021, IT_0020, mty::lt::reg_int);
    const complex_t IT_0032 = m_mu*IT_0031;
    const complex_t IT_0033 = 2*IT_0032;
    const complex_t IT_0034 = IT_0030 + IT_0033;
    const complex_t IT_0035 = IT_0019*IT_0034;
    const complex_t IT_0036 = cosq(alpha);
    const complex_t IT_0037 = (complex_t{0, 1})*e_em*m_mu*IT_0009*IT_0011
      *IT_0014*IT_0036;
    const complex_t IT_0038 = (-0.5)*IT_0037;
    const complex_t IT_0039 = cpowq(IT_0038, 2);
    const complex_t IT_0040 = 0.101321183642338*IT_0039;
    const complex_t IT_0041 = IT_0008*IT_0040;
    const complex_t IT_0042 = powq(m_H0, 2);
    const complex_t IT_0043 = mty::lt::C0iC(0, IT_0020, IT_0020, 2*IT_0020 + (
      -2)*s_12, IT_0020, IT_0042, IT_0020, mty::lt::reg_int);
    const complex_t IT_0044 = m_mu*IT_0043;
    const complex_t IT_0045 = mty::lt::C0iC(6, IT_0020, IT_0020, 2*IT_0020 + (
      -2)*s_12, IT_0020, IT_0042, IT_0020, mty::lt::reg_int);
    const complex_t IT_0046 = m_mu*IT_0045;
    const complex_t IT_0047 = IT_0044 + IT_0046;
    const complex_t IT_0048 = mty::lt::C0iC(12, IT_0020, IT_0020, 2*IT_0020 + 
      (-2)*s_12, IT_0020, IT_0042, IT_0020, mty::lt::reg_int);
    const complex_t IT_0049 = m_mu*IT_0048;
    const complex_t IT_0050 = mty::lt::C0iC(15, IT_0020, IT_0020, 2*IT_0020 + 
      (-2)*s_12, IT_0020, IT_0042, IT_0020, mty::lt::reg_int);
    const complex_t IT_0051 = m_mu*IT_0050;
    const complex_t IT_0052 = -IT_0049 + -IT_0051;
    const complex_t IT_0053 = IT_0047 + IT_0052;
    const complex_t IT_0054 = IT_0041*IT_0053;
    const complex_t IT_0055 = IT_0035 + IT_0054;
    const complex_t IT_0056 = IT_0001*IT_0055;
    const complex_t IT_0057 = IT_0025 + IT_0029;
    const complex_t IT_0058 = IT_0019*IT_0057;
    const complex_t IT_0059 = -IT_0051;
    const complex_t IT_0060 = IT_0046 + IT_0059;
    const complex_t IT_0061 = IT_0041*IT_0060;
    const complex_t IT_0062 = IT_0058 + IT_0061;
    const complex_t IT_0063 = IT_0001*IT_0062;
    return (complex_t{0, 0.125})*IT_0056 + (complex_t{0, (-0.125)})*IT_0063;
}
} // End of namespace c9_nmfv
