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
#include "C9A_H.h"
#include "common.h"

#include "params.h"
#include "libcomplexop.h"

namespace c9_nmfv {

complex_t C9A_H(
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
    const real_t s_12 = param.s_12;
    const real_t s_34 = param.s_34;
    const real_t theta_W = param.theta_W;
    const real_t V_ub_mod = param.V_ub_mod;
    const real_t reg_prop = param.reg_prop;
    const real_t delta_wolf = param.delta_wolf;
    const complex_t V_cs = param.V_cs;
    const complex_t V_ts = param.V_ts;
    const complex_t IT_0000 = powq(m_b + -m_s, -1);
    const complex_t IT_0001 = cosq(theta_W);
    const complex_t IT_0002 = cpowq(IT_0001, -1);
    const complex_t IT_0003 = tanq(theta_W);
    const complex_t IT_0004 = cpowq(IT_0003, 2);
    const complex_t IT_0005 = cpowq(1 + IT_0004, (-0.5));
    const complex_t IT_0006 = (complex_t{0, 1})*e_em*IT_0002*IT_0005;
    const complex_t IT_0007 = -IT_0006;
    const complex_t IT_0008 = powq(m_b, 2);
    const complex_t IT_0009 = powq(m_s, 2);
    const complex_t IT_0010 = powq(m_mu, 2);
    const complex_t IT_0011 = cpowq(s_34 + IT_0010 + 0.5*reg_prop, -1);
    const complex_t IT_0012 = powq(M_W, 2);
    const complex_t IT_0013 = powq(V_tb, -1);
    const complex_t IT_0014 = cpowq(conjq(V_ts), -1);
    const complex_t IT_0015 = powq(e_em, -4);
    const complex_t IT_0016 = 9.86960440108936*IT_0012*IT_0013*IT_0014*IT_0015;
    const complex_t IT_0017 = powq(M_W, -1);
    const complex_t IT_0018 = cosq(beta);
    const complex_t IT_0019 = sinq(beta);
    const complex_t IT_0020 = cpowq(IT_0019, -1);
    const complex_t IT_0021 = sinq(theta_W);
    const complex_t IT_0022 = cpowq(IT_0021, -1);
    const complex_t IT_0023 = (complex_t{0, 1.4142135623731})*m_t*V_tb*e_em
      *IT_0017*IT_0018*IT_0020*IT_0022;
    const complex_t IT_0024 = 0.5*IT_0023;
    const complex_t IT_0025 = (complex_t{0, 1.4142135623731})*m_t*conjq(V_ts)
      *e_em*IT_0017*IT_0018*IT_0020*IT_0022;
    const complex_t IT_0026 = 0.5*IT_0025;
    const complex_t IT_0027 = IT_0024*IT_0026;
    const complex_t IT_0028 = 0.101321183642338*IT_0027;
    const complex_t IT_0029 = powq(m_t, 2);
    const complex_t IT_0030 = powq(m_Hp, 2);
    const complex_t IT_0031 = mty::lt::C0iC(3, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0029, IT_0030, IT_0030, mty::lt::reg_int);
    const complex_t IT_0032 = IT_0006*IT_0031;
    const complex_t IT_0033 = m_b*IT_0032;
    const complex_t IT_0034 = mty::lt::C0iC(6, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0029, IT_0030, IT_0030, mty::lt::reg_int);
    const complex_t IT_0035 = IT_0006*IT_0034;
    const complex_t IT_0036 = m_s*IT_0035;
    const complex_t IT_0037 = 2*IT_0006;
    const complex_t IT_0038 = mty::lt::C0iC(15, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0029, IT_0030, IT_0030, mty::lt::reg_int);
    const complex_t IT_0039 = IT_0037*IT_0038;
    const complex_t IT_0040 = m_s*IT_0039;
    const complex_t IT_0041 = mty::lt::C0iC(12, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0029, IT_0030, IT_0030, mty::lt::reg_int);
    const complex_t IT_0042 = IT_0037*IT_0041;
    const complex_t IT_0043 = m_b*IT_0042;
    const complex_t IT_0044 = IT_0033 + IT_0036 + IT_0040 + IT_0043;
    const complex_t IT_0045 = IT_0028*IT_0044;
    const complex_t IT_0046 = cexpq((complex_t{0, -1})*delta_wolf);
    const complex_t IT_0047 = cpowq(IT_0018, -1);
    const complex_t IT_0048 = (complex_t{0, 1.4142135623731})*m_b*e_em*IT_0017
      *IT_0019*IT_0022*IT_0046*IT_0047*V_ub_mod;
    const complex_t IT_0049 = 0.5*IT_0048;
    const complex_t IT_0050 = (complex_t{0, 1.4142135623731})*m_s*V_us*e_em
      *IT_0017*IT_0019*IT_0022*IT_0047;
    const complex_t IT_0051 = 0.5*IT_0050;
    const complex_t IT_0052 = IT_0049*IT_0051;
    const complex_t IT_0053 = 0.101321183642338*IT_0052;
    const complex_t IT_0054 = powq(m_u, 2);
    const complex_t IT_0055 = mty::lt::C0iC(3, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0054, IT_0030, IT_0030, mty::lt::reg_int);
    const complex_t IT_0056 = IT_0006*IT_0055;
    const complex_t IT_0057 = m_b*IT_0056;
    const complex_t IT_0058 = mty::lt::C0iC(6, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0054, IT_0030, IT_0030, mty::lt::reg_int);
    const complex_t IT_0059 = IT_0006*IT_0058;
    const complex_t IT_0060 = m_s*IT_0059;
    const complex_t IT_0061 = mty::lt::C0iC(15, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0054, IT_0030, IT_0030, mty::lt::reg_int);
    const complex_t IT_0062 = IT_0037*IT_0061;
    const complex_t IT_0063 = m_s*IT_0062;
    const complex_t IT_0064 = mty::lt::C0iC(12, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0054, IT_0030, IT_0030, mty::lt::reg_int);
    const complex_t IT_0065 = IT_0037*IT_0064;
    const complex_t IT_0066 = m_b*IT_0065;
    const complex_t IT_0067 = IT_0057 + IT_0060 + IT_0063 + IT_0066;
    const complex_t IT_0068 = IT_0053*IT_0067;
    const complex_t IT_0069 = 0.101321183642338*m_c;
    const complex_t IT_0070 = 0.666666666666667*IT_0006;
    const complex_t IT_0071 = (complex_t{0, 1.4142135623731})*m_c*V_cb*e_em
      *IT_0017*IT_0018*IT_0020*IT_0022;
    const complex_t IT_0072 = 0.5*IT_0071;
    const complex_t IT_0073 = (complex_t{0, 1.4142135623731})*m_s*conjq(V_cs)
      *e_em*IT_0017*IT_0019*IT_0022*IT_0047;
    const complex_t IT_0074 = 0.5*IT_0073;
    const complex_t IT_0075 = IT_0070*IT_0072*IT_0074;
    const complex_t IT_0076 = IT_0069*IT_0075;
    const complex_t IT_0077 = powq(m_c, 2);
    const complex_t IT_0078 = mty::lt::C0iC(0, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0077, IT_0077, IT_0030, mty::lt::reg_int);
    const complex_t IT_0079 = mty::lt::C0iC(3, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0077, IT_0077, IT_0030, mty::lt::reg_int);
    const complex_t IT_0080 = mty::lt::C0iC(6, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0077, IT_0077, IT_0030, mty::lt::reg_int);
    const complex_t IT_0081 = IT_0078 + IT_0079 + IT_0080;
    const complex_t IT_0082 = IT_0076*IT_0081;
    const complex_t IT_0083 = (complex_t{0, 1.4142135623731})*m_c*conjq(V_cs)
      *e_em*IT_0017*IT_0018*IT_0020*IT_0022;
    const complex_t IT_0084 = 0.5*IT_0083;
    const complex_t IT_0085 = IT_0070*IT_0072*IT_0084;
    const complex_t IT_0086 = 0.101321183642338*IT_0085;
    const complex_t IT_0087 = m_s*IT_0079;
    const complex_t IT_0088 = mty::lt::C0iC(12, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0077, IT_0077, IT_0030, mty::lt::reg_int);
    const complex_t IT_0089 = m_s*IT_0088;
    const complex_t IT_0090 = mty::lt::C0iC(15, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0077, IT_0077, IT_0030, mty::lt::reg_int);
    const complex_t IT_0091 = m_s*IT_0090;
    const complex_t IT_0092 = IT_0087 + IT_0089 + IT_0091;
    const complex_t IT_0093 = m_b*IT_0088;
    const complex_t IT_0094 = m_b*IT_0090;
    const complex_t IT_0095 = mty::lt::C0iC(18, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0077, IT_0077, IT_0030, mty::lt::reg_int);
    const complex_t IT_0096 = m_b*IT_0095;
    const complex_t IT_0097 = m_b*IT_0079;
    const complex_t IT_0098 = m_b*IT_0080;
    const complex_t IT_0099 = -IT_0093 + (-2)*IT_0094 + -IT_0096 + -IT_0097 + 
      -IT_0098;
    const complex_t IT_0100 = IT_0092 + IT_0099;
    const complex_t IT_0101 = IT_0086*IT_0100;
    const complex_t IT_0102 = 0.101321183642338*m_t;
    const complex_t IT_0103 = (complex_t{0, 1.4142135623731})*m_s*conjq(V_ts)
      *e_em*IT_0017*IT_0019*IT_0022*IT_0047;
    const complex_t IT_0104 = 0.5*IT_0103;
    const complex_t IT_0105 = IT_0024*IT_0070*IT_0104;
    const complex_t IT_0106 = IT_0102*IT_0105;
    const complex_t IT_0107 = mty::lt::C0iC(0, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0029, IT_0029, IT_0030, mty::lt::reg_int);
    const complex_t IT_0108 = mty::lt::C0iC(3, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0029, IT_0029, IT_0030, mty::lt::reg_int);
    const complex_t IT_0109 = mty::lt::C0iC(6, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0029, IT_0029, IT_0030, mty::lt::reg_int);
    const complex_t IT_0110 = IT_0107 + IT_0108 + IT_0109;
    const complex_t IT_0111 = IT_0106*IT_0110;
    const complex_t IT_0112 = (complex_t{0, 1.4142135623731})*m_b*V_tb*e_em
      *IT_0017*IT_0019*IT_0022*IT_0047;
    const complex_t IT_0113 = 0.5*IT_0112;
    const complex_t IT_0114 = IT_0070*IT_0104*IT_0113;
    const complex_t IT_0115 = 0.101321183642338*IT_0114;
    const complex_t IT_0116 = m_s*IT_0108;
    const complex_t IT_0117 = mty::lt::C0iC(12, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0029, IT_0029, IT_0030, mty::lt::reg_int);
    const complex_t IT_0118 = m_s*IT_0117;
    const complex_t IT_0119 = mty::lt::C0iC(15, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0029, IT_0029, IT_0030, mty::lt::reg_int);
    const complex_t IT_0120 = m_s*IT_0119;
    const complex_t IT_0121 = IT_0116 + IT_0118 + IT_0120;
    const complex_t IT_0122 = m_b*IT_0119;
    const complex_t IT_0123 = m_b*IT_0109;
    const complex_t IT_0124 = m_b*IT_0117;
    const complex_t IT_0125 = mty::lt::C0iC(18, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0029, IT_0029, IT_0030, mty::lt::reg_int);
    const complex_t IT_0126 = m_b*IT_0125;
    const complex_t IT_0127 = m_b*IT_0108;
    const complex_t IT_0128 = (-2)*IT_0122 + -IT_0123 + -IT_0124 + -IT_0126 + 
      -IT_0127;
    const complex_t IT_0129 = IT_0121 + IT_0128;
    const complex_t IT_0130 = IT_0115*IT_0129;
    const complex_t IT_0131 = 0.101321183642338*m_u;
    const complex_t IT_0132 = (complex_t{0, 1.4142135623731})*m_u*e_em*IT_0017
      *IT_0018*IT_0020*IT_0022*IT_0046*V_ub_mod;
    const complex_t IT_0133 = 0.5*IT_0132;
    const complex_t IT_0134 = IT_0051*IT_0070*IT_0133;
    const complex_t IT_0135 = IT_0131*IT_0134;
    const complex_t IT_0136 = mty::lt::C0iC(0, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0054, IT_0054, IT_0030, mty::lt::reg_int);
    const complex_t IT_0137 = mty::lt::C0iC(3, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0054, IT_0054, IT_0030, mty::lt::reg_int);
    const complex_t IT_0138 = mty::lt::C0iC(6, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0054, IT_0054, IT_0030, mty::lt::reg_int);
    const complex_t IT_0139 = IT_0136 + IT_0137 + IT_0138;
    const complex_t IT_0140 = IT_0135*IT_0139;
    const complex_t IT_0141 = IT_0049*IT_0051*IT_0070;
    const complex_t IT_0142 = 0.101321183642338*IT_0141;
    const complex_t IT_0143 = m_s*IT_0137;
    const complex_t IT_0144 = mty::lt::C0iC(15, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0054, IT_0054, IT_0030, mty::lt::reg_int);
    const complex_t IT_0145 = m_s*IT_0144;
    const complex_t IT_0146 = mty::lt::C0iC(12, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0054, IT_0054, IT_0030, mty::lt::reg_int);
    const complex_t IT_0147 = m_s*IT_0146;
    const complex_t IT_0148 = IT_0143 + IT_0145 + IT_0147;
    const complex_t IT_0149 = mty::lt::C0iC(18, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0054, IT_0054, IT_0030, mty::lt::reg_int);
    const complex_t IT_0150 = m_b*IT_0149;
    const complex_t IT_0151 = m_b*IT_0137;
    const complex_t IT_0152 = m_b*IT_0138;
    const complex_t IT_0153 = m_b*IT_0144;
    const complex_t IT_0154 = m_b*IT_0146;
    const complex_t IT_0155 = -IT_0150 + -IT_0151 + -IT_0152 + (-2)*IT_0153 + 
      -IT_0154;
    const complex_t IT_0156 = IT_0148 + IT_0155;
    const complex_t IT_0157 = IT_0142*IT_0156;
    const complex_t IT_0158 = (complex_t{0, 1.4142135623731})*m_b*V_cb*e_em
      *IT_0017*IT_0019*IT_0022*IT_0047;
    const complex_t IT_0159 = 0.5*IT_0158;
    const complex_t IT_0160 = IT_0070*IT_0084*IT_0159;
    const complex_t IT_0161 = IT_0069*IT_0160;
    const complex_t IT_0162 = IT_0081*IT_0161;
    const complex_t IT_0163 = IT_0072*IT_0074;
    const complex_t IT_0164 = IT_0069*IT_0163;
    const complex_t IT_0165 = mty::lt::C0iC(0, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0077, IT_0030, IT_0030, mty::lt::reg_int);
    const complex_t IT_0166 = IT_0006*IT_0165;
    const complex_t IT_0167 = mty::lt::C0iC(3, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0077, IT_0030, IT_0030, mty::lt::reg_int);
    const complex_t IT_0168 = IT_0037*IT_0167;
    const complex_t IT_0169 = IT_0166 + IT_0168;
    const complex_t IT_0170 = IT_0164*IT_0169;
    const complex_t IT_0171 = IT_0024*IT_0104;
    const complex_t IT_0172 = IT_0102*IT_0171;
    const complex_t IT_0173 = mty::lt::C0iC(0, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0029, IT_0030, IT_0030, mty::lt::reg_int);
    const complex_t IT_0174 = IT_0006*IT_0173;
    const complex_t IT_0175 = IT_0031*IT_0037;
    const complex_t IT_0176 = IT_0174 + IT_0175;
    const complex_t IT_0177 = IT_0172*IT_0176;
    const complex_t IT_0178 = IT_0026*IT_0070*IT_0113;
    const complex_t IT_0179 = IT_0102*IT_0178;
    const complex_t IT_0180 = IT_0110*IT_0179;
    const complex_t IT_0181 = IT_0051*IT_0133;
    const complex_t IT_0182 = IT_0131*IT_0181;
    const complex_t IT_0183 = mty::lt::C0iC(0, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0054, IT_0030, IT_0030, mty::lt::reg_int);
    const complex_t IT_0184 = IT_0006*IT_0183;
    const complex_t IT_0185 = IT_0037*IT_0055;
    const complex_t IT_0186 = IT_0184 + IT_0185;
    const complex_t IT_0187 = IT_0182*IT_0186;
    const complex_t IT_0188 = (complex_t{0, 1.4142135623731})*m_u*V_us*e_em
      *IT_0017*IT_0018*IT_0020*IT_0022;
    const complex_t IT_0189 = 0.5*IT_0188;
    const complex_t IT_0190 = IT_0133*IT_0189;
    const complex_t IT_0191 = 0.101321183642338*IT_0190;
    const complex_t IT_0192 = IT_0067*IT_0191;
    const complex_t IT_0193 = IT_0024*IT_0026*IT_0070;
    const complex_t IT_0194 = 0.101321183642338*IT_0193;
    const complex_t IT_0195 = IT_0129*IT_0194;
    const complex_t IT_0196 = IT_0026*IT_0113;
    const complex_t IT_0197 = IT_0102*IT_0196;
    const complex_t IT_0198 = IT_0176*IT_0197;
    const complex_t IT_0199 = IT_0070*IT_0074*IT_0159;
    const complex_t IT_0200 = 0.101321183642338*IT_0199;
    const complex_t IT_0201 = IT_0100*IT_0200;
    const complex_t IT_0202 = mty::lt::C0iC(6, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0077, IT_0030, IT_0030, mty::lt::reg_int);
    const complex_t IT_0203 = IT_0006*IT_0202;
    const complex_t IT_0204 = m_s*IT_0203;
    const complex_t IT_0205 = mty::lt::C0iC(15, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0077, IT_0030, IT_0030, mty::lt::reg_int);
    const complex_t IT_0206 = IT_0037*IT_0205;
    const complex_t IT_0207 = m_s*IT_0206;
    const complex_t IT_0208 = IT_0006*IT_0167;
    const complex_t IT_0209 = m_b*IT_0208;
    const complex_t IT_0210 = mty::lt::C0iC(12, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0077, IT_0030, IT_0030, mty::lt::reg_int);
    const complex_t IT_0211 = IT_0037*IT_0210;
    const complex_t IT_0212 = m_b*IT_0211;
    const complex_t IT_0213 = IT_0204 + IT_0207 + IT_0209 + IT_0212;
    const complex_t IT_0214 = IT_0074*IT_0159;
    const complex_t IT_0215 = 0.101321183642338*IT_0214;
    const complex_t IT_0216 = IT_0213*IT_0215;
    const complex_t IT_0217 = IT_0084*IT_0159;
    const complex_t IT_0218 = IT_0069*IT_0217;
    const complex_t IT_0219 = IT_0169*IT_0218;
    const complex_t IT_0220 = IT_0049*IT_0189;
    const complex_t IT_0221 = IT_0131*IT_0220;
    const complex_t IT_0222 = IT_0186*IT_0221;
    const complex_t IT_0223 = IT_0104*IT_0113;
    const complex_t IT_0224 = 0.101321183642338*IT_0223;
    const complex_t IT_0225 = IT_0044*IT_0224;
    const complex_t IT_0226 = IT_0070*IT_0133*IT_0189;
    const complex_t IT_0227 = 0.101321183642338*IT_0226;
    const complex_t IT_0228 = IT_0156*IT_0227;
    const complex_t IT_0229 = IT_0072*IT_0084;
    const complex_t IT_0230 = 0.101321183642338*IT_0229;
    const complex_t IT_0231 = IT_0213*IT_0230;
    const complex_t IT_0232 = IT_0049*IT_0070*IT_0189;
    const complex_t IT_0233 = IT_0131*IT_0232;
    const complex_t IT_0234 = IT_0139*IT_0233;
    const complex_t IT_0235 = IT_0045 + IT_0068 + (-2)*IT_0082 + (-2)*IT_0101 
      + (-2)*IT_0111 + (-2)*IT_0130 + (-2)*IT_0140 + (-2)*IT_0157 + (-2)*IT_0162
       + -IT_0170 + -IT_0177 + (-2)*IT_0180 + -IT_0187 + IT_0192 + (-2)*IT_0195 
      + -IT_0198 + (-2)*IT_0201 + IT_0216 + -IT_0219 + -IT_0222 + IT_0225 + (-2)
      *IT_0228 + IT_0231 + (-2)*IT_0234;
    const complex_t IT_0236 = cpowq(IT_0021, 2);
    const complex_t IT_0237 = IT_0235*IT_0236;
    const complex_t IT_0238 = -IT_0237;
    const complex_t IT_0239 = IT_0016*IT_0238;
    const complex_t IT_0240 = -IT_0239;
    const complex_t IT_0241 = IT_0108*IT_0179;
    const complex_t IT_0242 = IT_0137*IT_0233;
    const complex_t IT_0243 = IT_0037*IT_0165;
    const complex_t IT_0244 = IT_0037*IT_0202;
    const complex_t IT_0245 = -IT_0243 + -IT_0244;
    const complex_t IT_0246 = IT_0166 + IT_0245;
    const complex_t IT_0247 = IT_0218*IT_0246;
    const complex_t IT_0248 = IT_0087 + IT_0089;
    const complex_t IT_0249 = -IT_0093 + -IT_0094 + -IT_0097;
    const complex_t IT_0250 = IT_0248 + IT_0249;
    const complex_t IT_0251 = IT_0200*IT_0250;
    const complex_t IT_0252 = IT_0037*IT_0183;
    const complex_t IT_0253 = IT_0037*IT_0058;
    const complex_t IT_0254 = -IT_0252 + -IT_0253;
    const complex_t IT_0255 = IT_0184 + IT_0254;
    const complex_t IT_0256 = IT_0221*IT_0255;
    const complex_t IT_0257 = IT_0076*IT_0079;
    const complex_t IT_0258 = IT_0086*IT_0250;
    const complex_t IT_0259 = IT_0057 + IT_0060;
    const complex_t IT_0260 = m_b*IT_0185;
    const complex_t IT_0261 = m_b*IT_0062;
    const complex_t IT_0262 = m_s*IT_0253;
    const complex_t IT_0263 = mty::lt::C0iC(18, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0054, IT_0030, IT_0030, mty::lt::reg_int);
    const complex_t IT_0264 = IT_0037*IT_0263;
    const complex_t IT_0265 = m_s*IT_0264;
    const complex_t IT_0266 = -IT_0260 + -IT_0261 + -IT_0262 + -IT_0265;
    const complex_t IT_0267 = IT_0259 + IT_0266;
    const complex_t IT_0268 = IT_0191*IT_0267;
    const complex_t IT_0269 = IT_0164*IT_0246;
    const complex_t IT_0270 = IT_0037*IT_0173;
    const complex_t IT_0271 = IT_0034*IT_0037;
    const complex_t IT_0272 = -IT_0270 + -IT_0271;
    const complex_t IT_0273 = IT_0174 + IT_0272;
    const complex_t IT_0274 = IT_0172*IT_0273;
    const complex_t IT_0275 = IT_0033 + IT_0036;
    const complex_t IT_0276 = m_b*IT_0175;
    const complex_t IT_0277 = m_b*IT_0039;
    const complex_t IT_0278 = m_s*IT_0271;
    const complex_t IT_0279 = mty::lt::C0iC(18, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0029, IT_0030, IT_0030, mty::lt::reg_int);
    const complex_t IT_0280 = IT_0037*IT_0279;
    const complex_t IT_0281 = m_s*IT_0280;
    const complex_t IT_0282 = -IT_0276 + -IT_0277 + -IT_0278 + -IT_0281;
    const complex_t IT_0283 = IT_0275 + IT_0282;
    const complex_t IT_0284 = IT_0028*IT_0283;
    const complex_t IT_0285 = IT_0204 + IT_0209;
    const complex_t IT_0286 = m_b*IT_0168;
    const complex_t IT_0287 = m_b*IT_0206;
    const complex_t IT_0288 = m_s*IT_0244;
    const complex_t IT_0289 = mty::lt::C0iC(18, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0077, IT_0030, IT_0030, mty::lt::reg_int);
    const complex_t IT_0290 = IT_0037*IT_0289;
    const complex_t IT_0291 = m_s*IT_0290;
    const complex_t IT_0292 = -IT_0286 + -IT_0287 + -IT_0288 + -IT_0291;
    const complex_t IT_0293 = IT_0285 + IT_0292;
    const complex_t IT_0294 = IT_0215*IT_0293;
    const complex_t IT_0295 = IT_0106*IT_0108;
    const complex_t IT_0296 = IT_0116 + IT_0118;
    const complex_t IT_0297 = -IT_0122 + -IT_0124 + -IT_0127;
    const complex_t IT_0298 = IT_0296 + IT_0297;
    const complex_t IT_0299 = IT_0194*IT_0298;
    const complex_t IT_0300 = IT_0135*IT_0137;
    const complex_t IT_0301 = IT_0143 + IT_0147;
    const complex_t IT_0302 = -IT_0151 + -IT_0153 + -IT_0154;
    const complex_t IT_0303 = IT_0301 + IT_0302;
    const complex_t IT_0304 = IT_0142*IT_0303;
    const complex_t IT_0305 = IT_0224*IT_0283;
    const complex_t IT_0306 = IT_0197*IT_0273;
    const complex_t IT_0307 = IT_0053*IT_0267;
    const complex_t IT_0308 = IT_0079*IT_0161;
    const complex_t IT_0309 = IT_0115*IT_0298;
    const complex_t IT_0310 = IT_0227*IT_0303;
    const complex_t IT_0311 = IT_0182*IT_0255;
    const complex_t IT_0312 = IT_0230*IT_0293;
    const complex_t IT_0313 = IT_0241 + IT_0242 + 0.5*IT_0247 + IT_0251 + 0.5
      *IT_0256 + IT_0257 + IT_0258 + (-0.5)*IT_0268 + 0.5*IT_0269 + 0.5*IT_0274 
      + (-0.5)*IT_0284 + (-0.5)*IT_0294 + IT_0295 + IT_0299 + IT_0300 + IT_0304 
      + (-0.5)*IT_0305 + 0.5*IT_0306 + (-0.5)*IT_0307 + IT_0308 + IT_0309 +
       IT_0310 + 0.5*IT_0311 + (-0.5)*IT_0312;
    const complex_t IT_0314 = IT_0236*IT_0313;
    const complex_t IT_0315 = 2*IT_0314;
    const complex_t IT_0316 = IT_0016*IT_0315;
    const complex_t IT_0317 = powq(m_b + m_s, -1);
    const complex_t IT_0318 = IT_0087 + IT_0089 + IT_0091 + IT_0093 + IT_0096 
      + IT_0097 + IT_0098;
    const complex_t IT_0319 = 2*IT_0094;
    const complex_t IT_0320 = IT_0318 + IT_0319;
    const complex_t IT_0321 = IT_0200*IT_0320;
    const complex_t IT_0322 = IT_0116 + IT_0118 + IT_0120 + IT_0123 + IT_0124 
      + IT_0126 + IT_0127;
    const complex_t IT_0323 = 2*IT_0122;
    const complex_t IT_0324 = IT_0322 + IT_0323;
    const complex_t IT_0325 = IT_0194*IT_0324;
    const complex_t IT_0326 = IT_0115*IT_0324;
    const complex_t IT_0327 = IT_0143 + IT_0145 + IT_0147 + IT_0150 + IT_0151 
      + IT_0152 + IT_0154;
    const complex_t IT_0328 = 2*IT_0153;
    const complex_t IT_0329 = IT_0327 + IT_0328;
    const complex_t IT_0330 = IT_0142*IT_0329;
    const complex_t IT_0331 = IT_0204 + IT_0207;
    const complex_t IT_0332 = -IT_0209 + -IT_0212;
    const complex_t IT_0333 = IT_0331 + IT_0332;
    const complex_t IT_0334 = IT_0230*IT_0333;
    const complex_t IT_0335 = IT_0036 + IT_0040;
    const complex_t IT_0336 = -IT_0033 + -IT_0043;
    const complex_t IT_0337 = IT_0335 + IT_0336;
    const complex_t IT_0338 = IT_0028*IT_0337;
    const complex_t IT_0339 = IT_0215*IT_0333;
    const complex_t IT_0340 = IT_0224*IT_0337;
    const complex_t IT_0341 = IT_0086*IT_0320;
    const complex_t IT_0342 = IT_0060 + IT_0063;
    const complex_t IT_0343 = -IT_0057 + -IT_0066;
    const complex_t IT_0344 = IT_0342 + IT_0343;
    const complex_t IT_0345 = IT_0191*IT_0344;
    const complex_t IT_0346 = IT_0053*IT_0344;
    const complex_t IT_0347 = IT_0227*IT_0329;
    const complex_t IT_0348 = IT_0082 + IT_0111 + IT_0140 + -IT_0162 + 0.5
      *IT_0170 + 0.5*IT_0177 + -IT_0180 + 0.5*IT_0187 + (-0.5)*IT_0198 + (-0.5)
      *IT_0219 + (-0.5)*IT_0222 + -IT_0234 + -IT_0321 + IT_0325 + -IT_0326 + 
      -IT_0330 + (-0.5)*IT_0334 + (-0.5)*IT_0338 + 0.5*IT_0339 + 0.5*IT_0340 +
       IT_0341 + (-0.5)*IT_0345 + 0.5*IT_0346 + IT_0347;
    const complex_t IT_0349 = IT_0236*IT_0348;
    const complex_t IT_0350 = (-2)*IT_0349;
    const complex_t IT_0351 = (-0.5)*IT_0350;
    const complex_t IT_0352 = 2*IT_0351;
    const complex_t IT_0353 = IT_0016*IT_0352;
    const complex_t IT_0354 = IT_0116 + IT_0118 + IT_0122 + IT_0124 + IT_0127;
    const complex_t IT_0355 = IT_0194*IT_0354;
    const complex_t IT_0356 = IT_0143 + IT_0147 + IT_0151 + IT_0153 + IT_0154;
    const complex_t IT_0357 = IT_0227*IT_0356;
    const complex_t IT_0358 = IT_0142*IT_0356;
    const complex_t IT_0359 = -IT_0209 + -IT_0288 + -IT_0291;
    const complex_t IT_0360 = IT_0204 + IT_0286 + IT_0287;
    const complex_t IT_0361 = IT_0359 + IT_0360;
    const complex_t IT_0362 = IT_0230*IT_0361;
    const complex_t IT_0363 = IT_0036 + IT_0276 + IT_0277;
    const complex_t IT_0364 = -IT_0033 + -IT_0278 + -IT_0281;
    const complex_t IT_0365 = IT_0363 + IT_0364;
    const complex_t IT_0366 = IT_0224*IT_0365;
    const complex_t IT_0367 = IT_0087 + IT_0089 + IT_0093 + IT_0094 + IT_0097;
    const complex_t IT_0368 = IT_0086*IT_0367;
    const complex_t IT_0369 = IT_0060 + IT_0260 + IT_0261;
    const complex_t IT_0370 = -IT_0057 + -IT_0262 + -IT_0265;
    const complex_t IT_0371 = IT_0369 + IT_0370;
    const complex_t IT_0372 = IT_0191*IT_0371;
    const complex_t IT_0373 = IT_0028*IT_0365;
    const complex_t IT_0374 = IT_0215*IT_0361;
    const complex_t IT_0375 = IT_0053*IT_0371;
    const complex_t IT_0376 = IT_0200*IT_0367;
    const complex_t IT_0377 = IT_0115*IT_0354;
    const complex_t IT_0378 = IT_0241 + IT_0242 + 0.5*IT_0247 + 0.5*IT_0256 + 
      -IT_0257 + (-0.5)*IT_0269 + (-0.5)*IT_0274 + -IT_0295 + -IT_0300 + 0.5
      *IT_0306 + IT_0308 + (-0.5)*IT_0311 + -IT_0355 + -IT_0357 + IT_0358 + 0.5
      *IT_0362 + (-0.5)*IT_0366 + -IT_0368 + 0.5*IT_0372 + 0.5*IT_0373 + (-0.5)
      *IT_0374 + (-0.5)*IT_0375 + IT_0376 + IT_0377;
    const complex_t IT_0379 = IT_0236*IT_0378;
    const complex_t IT_0380 = -IT_0379;
    const complex_t IT_0381 = 2*IT_0380;
    const complex_t IT_0382 = -IT_0381;
    const complex_t IT_0383 = 0.5*IT_0382;
    const complex_t IT_0384 = -IT_0383;
    const complex_t IT_0385 = 2*IT_0384;
    const complex_t IT_0386 = IT_0016*IT_0385;
    const complex_t IT_0387 = -IT_0386;
    const complex_t IT_0388 = (-0.5)*IT_0000*IT_0007*IT_0011*(IT_0240 + 
      -IT_0316)*(s_12 + -1./2*IT_0008 + -1./2*IT_0009 + -1./2*reg_prop) + (-0.5)
      *IT_0007*IT_0011*IT_0317*(IT_0353 + -IT_0387)*(s_12 + -1./2*IT_0008 + -1.
      /2*IT_0009 + -1./2*reg_prop);
    return IT_0388;
}
} // End of namespace c9_nmfv
