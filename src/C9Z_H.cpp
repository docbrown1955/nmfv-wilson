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
#include "C9Z_H.h"
#include "common.h"

#include "params.h"
#include "libcomplexop.h"

namespace c9_nmfv {

complex_t C9Z_H(
        param_t const &param
        )
{
    clearcache();
    const real_t M_W = param.M_W;
    const real_t M_Z = param.M_Z;
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
    const real_t s_34 = param.s_34;
    const real_t Finite = param.Finite;
    const real_t theta_W = param.theta_W;
    const real_t V_ub_mod = param.V_ub_mod;
    const real_t reg_prop = param.reg_prop;
    const real_t delta_wolf = param.delta_wolf;
    const complex_t V_cs = param.V_cs;
    const complex_t V_ts = param.V_ts;
    const complex_t IT_0000 = powq(M_W, 2);
    const complex_t IT_0001 = powq(V_tb, -1);
    const complex_t IT_0002 = cpowq(conjq(V_ts), -1);
    const complex_t IT_0003 = powq(e_em, -4);
    const complex_t IT_0004 = 9.86960440108936*IT_0000*IT_0001*IT_0002*IT_0003;
    const complex_t IT_0005 = powq(M_Z, 2);
    const complex_t IT_0006 = powq(m_mu, 2);
    const complex_t IT_0007 = cpowq((-2)*s_34 + IT_0005 + (-2)*IT_0006 + 
      -reg_prop, -1);
    const complex_t IT_0008 = cosq(theta_W);
    const complex_t IT_0009 = cpowq(IT_0008, -2);
    const complex_t IT_0010 = sinq(theta_W);
    const complex_t IT_0011 = tanq(theta_W);
    const complex_t IT_0012 = cpowq(IT_0011, 2);
    const complex_t IT_0013 = cpowq(1 + IT_0012, (-0.5));
    const complex_t IT_0014 = (complex_t{0, 1})*e_em*IT_0009*IT_0010*IT_0013;
    const complex_t IT_0015 = powq(m_b, 2);
    const complex_t IT_0016 = powq(m_s, 2);
    const complex_t IT_0017 = cpowq(IT_0015 + -IT_0016 + reg_prop, -1);
    const complex_t IT_0018 = IT_0009*IT_0010*IT_0013;
    const complex_t IT_0019 = e_em*IT_0018;
    const complex_t IT_0020 = cpowq(IT_0010, -1);
    const complex_t IT_0021 = IT_0013*IT_0020;
    const complex_t IT_0022 = e_em*IT_0021;
    const complex_t IT_0023 = (complex_t{0, 1})*(IT_0019 + 3*IT_0022);
    const complex_t IT_0024 = (-0.166666666666667)*IT_0023;
    const complex_t IT_0025 = powq(M_W, -1);
    const complex_t IT_0026 = cexpq((complex_t{0, -1})*delta_wolf);
    const complex_t IT_0027 = cosq(beta);
    const complex_t IT_0028 = sinq(beta);
    const complex_t IT_0029 = cpowq(IT_0028, -1);
    const complex_t IT_0030 = (complex_t{0, 1.4142135623731})*m_u*e_em*IT_0020
      *IT_0025*IT_0026*IT_0027*IT_0029*V_ub_mod;
    const complex_t IT_0031 = 0.5*IT_0030;
    const complex_t IT_0032 = (complex_t{0, 1.4142135623731})*m_u*V_us*e_em
      *IT_0020*IT_0025*IT_0027*IT_0029;
    const complex_t IT_0033 = 0.5*IT_0032;
    const complex_t IT_0034 = powq(m_u, 2);
    const complex_t IT_0035 = powq(m_Hp, 2);
    const complex_t IT_0036 = mty::lt::B0iC(3, IT_0015, IT_0034, IT_0035,
       mty::lt::reg_int);
    const complex_t IT_0037 = IT_0015*IT_0036;
    const complex_t IT_0038 = IT_0024*IT_0031*IT_0033*IT_0037;
    const complex_t IT_0039 = 0.101321183642338*IT_0017*IT_0038;
    const complex_t IT_0040 = IT_0014*IT_0039;
    const complex_t IT_0041 = (complex_t{0, 1})*(IT_0019 + -IT_0022);
    const complex_t IT_0042 = 0.5*IT_0041;
    const complex_t IT_0043 = IT_0039*IT_0042;
    const complex_t IT_0044 = 0.101321183642338*m_c;
    const complex_t IT_0045 = (complex_t{0, 1.4142135623731})*m_c*conjq(V_cs)
      *e_em*IT_0020*IT_0025*IT_0027*IT_0029;
    const complex_t IT_0046 = 0.5*IT_0045;
    const complex_t IT_0047 = cpowq(IT_0027, -1);
    const complex_t IT_0048 = (complex_t{0, 1.4142135623731})*m_b*V_cb*e_em
      *IT_0020*IT_0025*IT_0028*IT_0047;
    const complex_t IT_0049 = 0.5*IT_0048;
    const complex_t IT_0050 = powq(m_c, 2);
    const complex_t IT_0051 = mty::lt::B0iC(0, 0, IT_0050, IT_0035,
       mty::lt::reg_int);
    const complex_t IT_0052 = m_b*IT_0051;
    const complex_t IT_0053 = IT_0024*IT_0046*IT_0049*IT_0052;
    const complex_t IT_0054 = IT_0017*IT_0044*IT_0053;
    const complex_t IT_0055 = IT_0014*IT_0054;
    const complex_t IT_0056 = 0.101321183642338*m_s*m_t;
    const complex_t IT_0057 = (complex_t{0, 1.4142135623731})*m_t*V_tb*e_em
      *IT_0020*IT_0025*IT_0027*IT_0029;
    const complex_t IT_0058 = 0.5*IT_0057;
    const complex_t IT_0059 = (complex_t{0, 1.4142135623731})*m_s*conjq(V_ts)
      *e_em*IT_0020*IT_0025*IT_0028*IT_0047;
    const complex_t IT_0060 = 0.5*IT_0059;
    const complex_t IT_0061 = powq(m_t, 2);
    const complex_t IT_0062 = mty::lt::B0iC(0, 0, IT_0061, IT_0035,
       mty::lt::reg_int);
    const complex_t IT_0063 = IT_0024*IT_0058*IT_0060*IT_0062;
    const complex_t IT_0064 = IT_0017*IT_0056*IT_0063;
    const complex_t IT_0065 = IT_0042*IT_0064;
    const complex_t IT_0066 = (complex_t{0, 1.4142135623731})*m_t*conjq(V_ts)
      *e_em*IT_0020*IT_0025*IT_0027*IT_0029;
    const complex_t IT_0067 = 0.5*IT_0066;
    const complex_t IT_0068 = mty::lt::B0iC(3, IT_0015, IT_0061, IT_0035,
       mty::lt::reg_int);
    const complex_t IT_0069 = IT_0015*IT_0068;
    const complex_t IT_0070 = IT_0024*IT_0058*IT_0067*IT_0069;
    const complex_t IT_0071 = 0.101321183642338*IT_0017*IT_0070;
    const complex_t IT_0072 = IT_0014*IT_0071;
    const complex_t IT_0073 = IT_0042*IT_0071;
    const complex_t IT_0074 = 0.101321183642338*m_s;
    const complex_t IT_0075 = (complex_t{0, 1.4142135623731})*m_s*conjq(V_cs)
      *e_em*IT_0020*IT_0025*IT_0028*IT_0047;
    const complex_t IT_0076 = 0.5*IT_0075;
    const complex_t IT_0077 = mty::lt::B0iC(3, IT_0015, IT_0050, IT_0035,
       mty::lt::reg_int);
    const complex_t IT_0078 = m_b*IT_0077;
    const complex_t IT_0079 = IT_0024*IT_0049*IT_0076*IT_0078;
    const complex_t IT_0080 = IT_0017*IT_0074*IT_0079;
    const complex_t IT_0081 = IT_0014*IT_0080;
    const complex_t IT_0082 = IT_0042*IT_0080;
    const complex_t IT_0083 = (complex_t{0, 1.4142135623731})*m_b*e_em*IT_0020
      *IT_0025*IT_0026*IT_0028*IT_0047*V_ub_mod;
    const complex_t IT_0084 = 0.5*IT_0083;
    const complex_t IT_0085 = (complex_t{0, 1.4142135623731})*m_s*V_us*e_em
      *IT_0020*IT_0025*IT_0028*IT_0047;
    const complex_t IT_0086 = 0.5*IT_0085;
    const complex_t IT_0087 = m_b*IT_0036;
    const complex_t IT_0088 = IT_0024*IT_0084*IT_0086*IT_0087;
    const complex_t IT_0089 = IT_0017*IT_0074*IT_0088;
    const complex_t IT_0090 = IT_0014*IT_0089;
    const complex_t IT_0091 = (complex_t{0, 1.4142135623731})*m_b*V_tb*e_em
      *IT_0020*IT_0025*IT_0028*IT_0047;
    const complex_t IT_0092 = 0.5*IT_0091;
    const complex_t IT_0093 = m_b*IT_0068;
    const complex_t IT_0094 = IT_0024*IT_0060*IT_0092*IT_0093;
    const complex_t IT_0095 = IT_0017*IT_0074*IT_0094;
    const complex_t IT_0096 = IT_0014*IT_0095;
    const complex_t IT_0097 = 0.101321183642338*IT_0061;
    const complex_t IT_0098 = (complex_t{0, 1})*(IT_0019 + (-3)*IT_0022);
    const complex_t IT_0099 = (-0.166666666666667)*IT_0098;
    const complex_t IT_0100 = mty::lt::C0iC(0, 0, 0, 0, IT_0061, IT_0061,
       IT_0035, mty::lt::reg_int);
    const complex_t IT_0101 = IT_0058*IT_0067*IT_0099*IT_0100;
    const complex_t IT_0102 = IT_0097*IT_0101;
    const complex_t IT_0103 = IT_0014*IT_0102;
    const complex_t IT_0104 = cpowq(IT_0015 + -IT_0016 + -reg_prop, -1);
    const complex_t IT_0105 = 0.101321183642338*m_b*m_t;
    const complex_t IT_0106 = IT_0024*IT_0062*IT_0067*IT_0092;
    const complex_t IT_0107 = IT_0104*IT_0105*IT_0106;
    const complex_t IT_0108 = IT_0014*IT_0107;
    const complex_t IT_0109 = IT_0042*IT_0107;
    const complex_t IT_0110 = 0.333333333333333*IT_0014;
    const complex_t IT_0111 = (complex_t{0, 1.4142135623731})*m_c*V_cb*e_em
      *IT_0020*IT_0025*IT_0027*IT_0029;
    const complex_t IT_0112 = 0.5*IT_0111;
    const complex_t IT_0113 = IT_0046*IT_0078*IT_0110*IT_0112;
    const complex_t IT_0114 = IT_0017*IT_0074*IT_0113;
    const complex_t IT_0115 = IT_0014*IT_0114;
    const complex_t IT_0116 = IT_0042*IT_0114;
    const complex_t IT_0117 = IT_0058*IT_0060*IT_0062*IT_0110;
    const complex_t IT_0118 = IT_0104*IT_0105*IT_0117;
    const complex_t IT_0119 = IT_0014*IT_0118;
    const complex_t IT_0120 = IT_0042*IT_0118;
    const complex_t IT_0121 = 0.101321183642338*m_u;
    const complex_t IT_0122 = mty::lt::B0iC(0, 0, IT_0034, IT_0035,
       mty::lt::reg_int);
    const complex_t IT_0123 = m_s*IT_0122;
    const complex_t IT_0124 = IT_0033*IT_0084*IT_0110*IT_0123;
    const complex_t IT_0125 = IT_0104*IT_0121*IT_0124;
    const complex_t IT_0126 = IT_0014*IT_0125;
    const complex_t IT_0127 = (-0.666666666666667)*IT_0014;
    const complex_t IT_0128 = IT_0060*IT_0092*IT_0100*IT_0127;
    const complex_t IT_0129 = IT_0097*IT_0128;
    const complex_t IT_0130 = IT_0042*IT_0129;
    const complex_t IT_0131 = mty::lt::C0iC(9, 0, 0, 0, IT_0050, IT_0050,
       IT_0035, mty::lt::reg_int);
    const complex_t IT_0132 = (-4)*IT_0131;
    const complex_t IT_0133 = Finite + IT_0132;
    const complex_t IT_0134 = IT_0046*IT_0112*IT_0127*IT_0133;
    const complex_t IT_0135 = 0.101321183642338*IT_0134;
    const complex_t IT_0136 = IT_0042*IT_0135;
    const complex_t IT_0137 = m_s*IT_0051;
    const complex_t IT_0138 = IT_0046*IT_0049*IT_0110*IT_0137;
    const complex_t IT_0139 = IT_0044*IT_0104*IT_0138;
    const complex_t IT_0140 = IT_0042*IT_0139;
    const complex_t IT_0141 = IT_0042*IT_0125;
    const complex_t IT_0142 = 0.101321183642338*m_c*m_s;
    const complex_t IT_0143 = IT_0024*IT_0051*IT_0076*IT_0112;
    const complex_t IT_0144 = IT_0017*IT_0142*IT_0143;
    const complex_t IT_0145 = IT_0042*IT_0144;
    const complex_t IT_0146 = m_b*IT_0122;
    const complex_t IT_0147 = IT_0024*IT_0033*IT_0084*IT_0146;
    const complex_t IT_0148 = IT_0017*IT_0121*IT_0147;
    const complex_t IT_0149 = IT_0014*IT_0148;
    const complex_t IT_0150 = 0.101321183642338*IT_0034;
    const complex_t IT_0151 = mty::lt::C0iC(0, 0, 0, 0, IT_0034, IT_0034,
       IT_0035, mty::lt::reg_int);
    const complex_t IT_0152 = IT_0031*IT_0033*IT_0099*IT_0151;
    const complex_t IT_0153 = IT_0150*IT_0152;
    const complex_t IT_0154 = IT_0042*IT_0153;
    const complex_t IT_0155 = IT_0024*IT_0076*IT_0112*IT_0137;
    const complex_t IT_0156 = IT_0044*IT_0104*IT_0155;
    const complex_t IT_0157 = IT_0042*IT_0156;
    const complex_t IT_0158 = 0.101321183642338*m_b;
    const complex_t IT_0159 = mty::lt::B0iC(3, IT_0016, IT_0061, IT_0035,
       mty::lt::reg_int);
    const complex_t IT_0160 = m_s*IT_0159;
    const complex_t IT_0161 = IT_0024*IT_0060*IT_0092*IT_0160;
    const complex_t IT_0162 = IT_0104*IT_0158*IT_0161;
    const complex_t IT_0163 = IT_0014*IT_0162;
    const complex_t IT_0164 = mty::lt::B0iC(3, IT_0016, IT_0034, IT_0035,
       mty::lt::reg_int);
    const complex_t IT_0165 = m_s*IT_0164;
    const complex_t IT_0166 = IT_0031*IT_0033*IT_0110*IT_0165;
    const complex_t IT_0167 = IT_0104*IT_0158*IT_0166;
    const complex_t IT_0168 = IT_0042*IT_0167;
    const complex_t IT_0169 = IT_0062*IT_0067*IT_0092*IT_0110;
    const complex_t IT_0170 = IT_0017*IT_0056*IT_0169;
    const complex_t IT_0171 = IT_0042*IT_0170;
    const complex_t IT_0172 = IT_0042*IT_0148;
    const complex_t IT_0173 = mty::lt::C0iC(9, 0, 0, 0, IT_0061, IT_0061,
       IT_0035, mty::lt::reg_int);
    const complex_t IT_0174 = (-4)*IT_0173;
    const complex_t IT_0175 = Finite + IT_0174;
    const complex_t IT_0176 = IT_0060*IT_0092*IT_0099*IT_0175;
    const complex_t IT_0177 = 0.101321183642338*IT_0176;
    const complex_t IT_0178 = IT_0014*IT_0177;
    const complex_t IT_0179 = IT_0042*IT_0177;
    const complex_t IT_0180 = IT_0014*IT_0153;
    const complex_t IT_0181 = (complex_t{0, 1})*e_em*IT_0013*(IT_0009*IT_0010 
      + -IT_0020);
    const complex_t IT_0182 = mty::lt::C0iC(9, 0, 0, 0, IT_0050, IT_0035,
       IT_0035, mty::lt::reg_int);
    const complex_t IT_0183 = IT_0181*IT_0182;
    const complex_t IT_0184 = IT_0046*IT_0112*IT_0183;
    const complex_t IT_0185 = 0.101321183642338*IT_0184;
    const complex_t IT_0186 = IT_0014*IT_0185;
    const complex_t IT_0187 = IT_0042*IT_0185;
    const complex_t IT_0188 = mty::lt::C0iC(9, 0, 0, 0, IT_0061, IT_0035,
       IT_0035, mty::lt::reg_int);
    const complex_t IT_0189 = IT_0181*IT_0188;
    const complex_t IT_0190 = IT_0058*IT_0067*IT_0189;
    const complex_t IT_0191 = 0.101321183642338*IT_0190;
    const complex_t IT_0192 = IT_0014*IT_0191;
    const complex_t IT_0193 = IT_0042*IT_0191;
    const complex_t IT_0194 = IT_0049*IT_0076*IT_0183;
    const complex_t IT_0195 = 0.101321183642338*IT_0194;
    const complex_t IT_0196 = IT_0014*IT_0195;
    const complex_t IT_0197 = IT_0042*IT_0195;
    const complex_t IT_0198 = 0.101321183642338*m_s*m_u;
    const complex_t IT_0199 = IT_0024*IT_0031*IT_0086*IT_0122;
    const complex_t IT_0200 = IT_0017*IT_0198*IT_0199;
    const complex_t IT_0201 = IT_0014*IT_0200;
    const complex_t IT_0202 = IT_0042*IT_0200;
    const complex_t IT_0203 = IT_0014*IT_0144;
    const complex_t IT_0204 = IT_0015*IT_0077;
    const complex_t IT_0205 = IT_0024*IT_0046*IT_0112*IT_0204;
    const complex_t IT_0206 = 0.101321183642338*IT_0017*IT_0205;
    const complex_t IT_0207 = IT_0014*IT_0206;
    const complex_t IT_0208 = IT_0042*IT_0206;
    const complex_t IT_0209 = IT_0014*IT_0064;
    const complex_t IT_0210 = IT_0042*IT_0054;
    const complex_t IT_0211 = IT_0042*IT_0089;
    const complex_t IT_0212 = IT_0042*IT_0095;
    const complex_t IT_0213 = 0.101321183642338*m_t;
    const complex_t IT_0214 = m_b*IT_0062;
    const complex_t IT_0215 = IT_0024*IT_0067*IT_0092*IT_0214;
    const complex_t IT_0216 = IT_0017*IT_0213*IT_0215;
    const complex_t IT_0217 = IT_0014*IT_0216;
    const complex_t IT_0218 = IT_0042*IT_0216;
    const complex_t IT_0219 = IT_0042*IT_0102;
    const complex_t IT_0220 = mty::lt::C0iC(9, 0, 0, 0, IT_0034, IT_0034,
       IT_0035, mty::lt::reg_int);
    const complex_t IT_0221 = (-4)*IT_0220;
    const complex_t IT_0222 = Finite + IT_0221;
    const complex_t IT_0223 = IT_0084*IT_0086*IT_0099*IT_0222;
    const complex_t IT_0224 = 0.101321183642338*IT_0223;
    const complex_t IT_0225 = IT_0014*IT_0224;
    const complex_t IT_0226 = IT_0042*IT_0224;
    const complex_t IT_0227 = IT_0024*IT_0031*IT_0086*IT_0123;
    const complex_t IT_0228 = IT_0104*IT_0121*IT_0227;
    const complex_t IT_0229 = IT_0014*IT_0228;
    const complex_t IT_0230 = IT_0042*IT_0228;
    const complex_t IT_0231 = IT_0016*IT_0164;
    const complex_t IT_0232 = IT_0024*IT_0031*IT_0033*IT_0231;
    const complex_t IT_0233 = 0.101321183642338*IT_0104*IT_0232;
    const complex_t IT_0234 = IT_0014*IT_0233;
    const complex_t IT_0235 = IT_0042*IT_0233;
    const complex_t IT_0236 = IT_0014*IT_0156;
    const complex_t IT_0237 = mty::lt::B0iC(3, IT_0016, IT_0050, IT_0035,
       mty::lt::reg_int);
    const complex_t IT_0238 = IT_0016*IT_0237;
    const complex_t IT_0239 = IT_0024*IT_0046*IT_0112*IT_0238;
    const complex_t IT_0240 = 0.101321183642338*IT_0104*IT_0239;
    const complex_t IT_0241 = IT_0014*IT_0240;
    const complex_t IT_0242 = IT_0042*IT_0240;
    const complex_t IT_0243 = m_s*IT_0062;
    const complex_t IT_0244 = IT_0024*IT_0058*IT_0060*IT_0243;
    const complex_t IT_0245 = IT_0104*IT_0213*IT_0244;
    const complex_t IT_0246 = IT_0014*IT_0245;
    const complex_t IT_0247 = IT_0042*IT_0245;
    const complex_t IT_0248 = IT_0016*IT_0159;
    const complex_t IT_0249 = IT_0024*IT_0058*IT_0067*IT_0248;
    const complex_t IT_0250 = 0.101321183642338*IT_0104*IT_0249;
    const complex_t IT_0251 = IT_0014*IT_0250;
    const complex_t IT_0252 = IT_0042*IT_0250;
    const complex_t IT_0253 = m_s*IT_0237;
    const complex_t IT_0254 = IT_0024*IT_0049*IT_0076*IT_0253;
    const complex_t IT_0255 = IT_0104*IT_0158*IT_0254;
    const complex_t IT_0256 = IT_0014*IT_0255;
    const complex_t IT_0257 = IT_0042*IT_0255;
    const complex_t IT_0258 = 0.101321183642338*m_b*m_c;
    const complex_t IT_0259 = IT_0024*IT_0046*IT_0049*IT_0051;
    const complex_t IT_0260 = IT_0104*IT_0258*IT_0259;
    const complex_t IT_0261 = IT_0014*IT_0260;
    const complex_t IT_0262 = IT_0042*IT_0260;
    const complex_t IT_0263 = IT_0024*IT_0084*IT_0086*IT_0165;
    const complex_t IT_0264 = IT_0104*IT_0158*IT_0263;
    const complex_t IT_0265 = IT_0014*IT_0264;
    const complex_t IT_0266 = IT_0042*IT_0264;
    const complex_t IT_0267 = 0.101321183642338*m_b*m_u;
    const complex_t IT_0268 = IT_0024*IT_0033*IT_0084*IT_0122;
    const complex_t IT_0269 = IT_0104*IT_0267*IT_0268;
    const complex_t IT_0270 = IT_0014*IT_0269;
    const complex_t IT_0271 = IT_0042*IT_0269;
    const complex_t IT_0272 = IT_0042*IT_0162;
    const complex_t IT_0273 = 0.101321183642338*IT_0050;
    const complex_t IT_0274 = mty::lt::C0iC(0, 0, 0, 0, IT_0050, IT_0050,
       IT_0035, mty::lt::reg_int);
    const complex_t IT_0275 = IT_0046*IT_0099*IT_0112*IT_0274;
    const complex_t IT_0276 = IT_0273*IT_0275;
    const complex_t IT_0277 = IT_0014*IT_0276;
    const complex_t IT_0278 = IT_0042*IT_0276;
    const complex_t IT_0279 = IT_0049*IT_0076*IT_0099*IT_0133;
    const complex_t IT_0280 = 0.101321183642338*IT_0279;
    const complex_t IT_0281 = IT_0014*IT_0280;
    const complex_t IT_0282 = IT_0042*IT_0280;
    const complex_t IT_0283 = mty::lt::C0iC(9, 0, 0, 0, IT_0034, IT_0035,
       IT_0035, mty::lt::reg_int);
    const complex_t IT_0284 = IT_0181*IT_0283;
    const complex_t IT_0285 = IT_0031*IT_0033*IT_0284;
    const complex_t IT_0286 = 0.101321183642338*IT_0285;
    const complex_t IT_0287 = IT_0014*IT_0286;
    const complex_t IT_0288 = IT_0042*IT_0286;
    const complex_t IT_0289 = IT_0084*IT_0086*IT_0284;
    const complex_t IT_0290 = 0.101321183642338*IT_0289;
    const complex_t IT_0291 = IT_0014*IT_0290;
    const complex_t IT_0292 = IT_0042*IT_0290;
    const complex_t IT_0293 = IT_0060*IT_0092*IT_0189;
    const complex_t IT_0294 = 0.101321183642338*IT_0293;
    const complex_t IT_0295 = IT_0014*IT_0294;
    const complex_t IT_0296 = IT_0042*IT_0294;
    const complex_t IT_0297 = IT_0031*IT_0086*IT_0110*IT_0122;
    const complex_t IT_0298 = IT_0104*IT_0267*IT_0297;
    const complex_t IT_0299 = IT_0014*IT_0298;
    const complex_t IT_0300 = IT_0042*IT_0298;
    const complex_t IT_0301 = IT_0014*IT_0167;
    const complex_t IT_0302 = IT_0031*IT_0086*IT_0110*IT_0146;
    const complex_t IT_0303 = IT_0017*IT_0121*IT_0302;
    const complex_t IT_0304 = IT_0014*IT_0303;
    const complex_t IT_0305 = IT_0042*IT_0303;
    const complex_t IT_0306 = IT_0031*IT_0033*IT_0087*IT_0110;
    const complex_t IT_0307 = IT_0017*IT_0074*IT_0306;
    const complex_t IT_0308 = IT_0014*IT_0307;
    const complex_t IT_0309 = IT_0042*IT_0307;
    const complex_t IT_0310 = IT_0031*IT_0033*IT_0127*IT_0222;
    const complex_t IT_0311 = 0.101321183642338*IT_0310;
    const complex_t IT_0312 = IT_0014*IT_0311;
    const complex_t IT_0313 = IT_0042*IT_0311;
    const complex_t IT_0314 = IT_0051*IT_0076*IT_0110*IT_0112;
    const complex_t IT_0315 = IT_0104*IT_0258*IT_0314;
    const complex_t IT_0316 = IT_0014*IT_0315;
    const complex_t IT_0317 = IT_0042*IT_0315;
    const complex_t IT_0318 = IT_0046*IT_0110*IT_0112*IT_0253;
    const complex_t IT_0319 = IT_0104*IT_0158*IT_0318;
    const complex_t IT_0320 = IT_0014*IT_0319;
    const complex_t IT_0321 = IT_0042*IT_0319;
    const complex_t IT_0322 = IT_0014*IT_0135;
    const complex_t IT_0323 = IT_0052*IT_0076*IT_0110*IT_0112;
    const complex_t IT_0324 = IT_0017*IT_0044*IT_0323;
    const complex_t IT_0325 = IT_0014*IT_0324;
    const complex_t IT_0326 = IT_0042*IT_0324;
    const complex_t IT_0327 = IT_0058*IT_0067*IT_0110*IT_0160;
    const complex_t IT_0328 = IT_0104*IT_0158*IT_0327;
    const complex_t IT_0329 = IT_0014*IT_0328;
    const complex_t IT_0330 = IT_0042*IT_0328;
    const complex_t IT_0331 = IT_0058*IT_0060*IT_0110*IT_0214;
    const complex_t IT_0332 = IT_0017*IT_0213*IT_0331;
    const complex_t IT_0333 = IT_0014*IT_0332;
    const complex_t IT_0334 = IT_0042*IT_0332;
    const complex_t IT_0335 = IT_0058*IT_0067*IT_0093*IT_0110;
    const complex_t IT_0336 = IT_0017*IT_0074*IT_0335;
    const complex_t IT_0337 = IT_0014*IT_0336;
    const complex_t IT_0338 = IT_0042*IT_0336;
    const complex_t IT_0339 = IT_0058*IT_0067*IT_0127*IT_0175;
    const complex_t IT_0340 = 0.101321183642338*IT_0339;
    const complex_t IT_0341 = IT_0014*IT_0340;
    const complex_t IT_0342 = IT_0042*IT_0340;
    const complex_t IT_0343 = IT_0049*IT_0076*IT_0110*IT_0238;
    const complex_t IT_0344 = 0.101321183642338*IT_0104*IT_0343;
    const complex_t IT_0345 = IT_0014*IT_0344;
    const complex_t IT_0346 = IT_0042*IT_0344;
    const complex_t IT_0347 = IT_0014*IT_0139;
    const complex_t IT_0348 = IT_0049*IT_0076*IT_0127*IT_0274;
    const complex_t IT_0349 = IT_0273*IT_0348;
    const complex_t IT_0350 = IT_0014*IT_0349;
    const complex_t IT_0351 = IT_0042*IT_0349;
    const complex_t IT_0352 = IT_0049*IT_0076*IT_0110*IT_0204;
    const complex_t IT_0353 = 0.101321183642338*IT_0017*IT_0352;
    const complex_t IT_0354 = IT_0014*IT_0353;
    const complex_t IT_0355 = IT_0042*IT_0353;
    const complex_t IT_0356 = IT_0046*IT_0049*IT_0051*IT_0110;
    const complex_t IT_0357 = IT_0017*IT_0142*IT_0356;
    const complex_t IT_0358 = IT_0014*IT_0357;
    const complex_t IT_0359 = IT_0042*IT_0357;
    const complex_t IT_0360 = IT_0084*IT_0086*IT_0110*IT_0231;
    const complex_t IT_0361 = 0.101321183642338*IT_0104*IT_0360;
    const complex_t IT_0362 = IT_0014*IT_0361;
    const complex_t IT_0363 = IT_0042*IT_0361;
    const complex_t IT_0364 = IT_0060*IT_0092*IT_0110*IT_0248;
    const complex_t IT_0365 = 0.101321183642338*IT_0104*IT_0364;
    const complex_t IT_0366 = IT_0014*IT_0365;
    const complex_t IT_0367 = IT_0042*IT_0365;
    const complex_t IT_0368 = IT_0067*IT_0092*IT_0110*IT_0243;
    const complex_t IT_0369 = IT_0104*IT_0213*IT_0368;
    const complex_t IT_0370 = IT_0014*IT_0369;
    const complex_t IT_0371 = IT_0042*IT_0369;
    const complex_t IT_0372 = IT_0037*IT_0084*IT_0086*IT_0110;
    const complex_t IT_0373 = 0.101321183642338*IT_0017*IT_0372;
    const complex_t IT_0374 = IT_0014*IT_0373;
    const complex_t IT_0375 = IT_0042*IT_0373;
    const complex_t IT_0376 = IT_0033*IT_0084*IT_0110*IT_0122;
    const complex_t IT_0377 = IT_0017*IT_0198*IT_0376;
    const complex_t IT_0378 = IT_0014*IT_0377;
    const complex_t IT_0379 = IT_0042*IT_0377;
    const complex_t IT_0380 = IT_0084*IT_0086*IT_0127*IT_0151;
    const complex_t IT_0381 = IT_0150*IT_0380;
    const complex_t IT_0382 = IT_0014*IT_0381;
    const complex_t IT_0383 = IT_0042*IT_0381;
    const complex_t IT_0384 = IT_0060*IT_0069*IT_0092*IT_0110;
    const complex_t IT_0385 = 0.101321183642338*IT_0017*IT_0384;
    const complex_t IT_0386 = IT_0014*IT_0385;
    const complex_t IT_0387 = IT_0042*IT_0385;
    const complex_t IT_0388 = IT_0014*IT_0170;
    const complex_t IT_0389 = IT_0014*IT_0129;
    const complex_t IT_0390 = IT_0040 + IT_0043 + -IT_0055 + -IT_0065 +
       IT_0072 + IT_0073 + IT_0081 + IT_0082 + IT_0090 + IT_0096 + -IT_0103 +
       IT_0108 + IT_0109 + IT_0115 + IT_0116 + IT_0119 + IT_0120 + IT_0126 + 
      -IT_0130 + (-0.5)*IT_0136 + IT_0140 + IT_0141 + -IT_0145 + -IT_0149 + 
      -IT_0154 + IT_0157 + -IT_0163 + -IT_0168 + -IT_0171 + -IT_0172 + (-0.5)
      *IT_0178 + (-0.5)*IT_0179 + -IT_0180 + IT_0186 + IT_0187 + IT_0192 +
       IT_0193 + IT_0196 + IT_0197 + -IT_0201 + -IT_0202 + -IT_0203 + IT_0207 +
       IT_0208 + -IT_0209 + -IT_0210 + IT_0211 + IT_0212 + -IT_0217 + -IT_0218 +
       -IT_0219 + (-0.5)*IT_0225 + (-0.5)*IT_0226 + IT_0229 + IT_0230 + -IT_0234
       + -IT_0235 + IT_0236 + -IT_0241 + -IT_0242 + IT_0246 + IT_0247 + -IT_0251
       + -IT_0252 + -IT_0256 + -IT_0257 + IT_0261 + IT_0262 + -IT_0265 + 
      -IT_0266 + IT_0270 + IT_0271 + -IT_0272 + -IT_0277 + -IT_0278 + (-0.5)
      *IT_0281 + (-0.5)*IT_0282 + IT_0287 + IT_0288 + IT_0291 + IT_0292 +
       IT_0295 + IT_0296 + IT_0299 + IT_0300 + -IT_0301 + -IT_0304 + -IT_0305 +
       IT_0308 + IT_0309 + (-0.5)*IT_0312 + (-0.5)*IT_0313 + IT_0316 + IT_0317 +
       -IT_0320 + -IT_0321 + (-0.5)*IT_0322 + -IT_0325 + -IT_0326 + -IT_0329 + 
      -IT_0330 + -IT_0333 + -IT_0334 + IT_0337 + IT_0338 + (-0.5)*IT_0341 + (
      -0.5)*IT_0342 + -IT_0345 + -IT_0346 + IT_0347 + -IT_0350 + -IT_0351 +
       IT_0354 + IT_0355 + -IT_0358 + -IT_0359 + -IT_0362 + -IT_0363 + -IT_0366 
      + -IT_0367 + IT_0370 + IT_0371 + IT_0374 + IT_0375 + -IT_0378 + -IT_0379 +
       -IT_0382 + -IT_0383 + IT_0386 + IT_0387 + -IT_0388 + -IT_0389;
    const complex_t IT_0391 = cpowq(IT_0010, 2);
    const complex_t IT_0392 = IT_0007*IT_0390*IT_0391;
    const complex_t IT_0393 = -IT_0392;
    const complex_t IT_0394 = IT_0004*IT_0393;
    const complex_t IT_0395 = (-0.5)*IT_0394;
    const complex_t IT_0396 = IT_0040 + IT_0043 + -IT_0055 + -IT_0065 +
       IT_0072 + IT_0073 + IT_0081 + IT_0082 + IT_0090 + IT_0096 + -IT_0103 +
       IT_0108 + IT_0109 + -IT_0115 + -IT_0116 + -IT_0119 + -IT_0120 + -IT_0126 
      + IT_0130 + (-0.5)*IT_0136 + -IT_0140 + -IT_0141 + -IT_0145 + -IT_0149 + 
      -IT_0154 + IT_0157 + -IT_0163 + IT_0168 + IT_0171 + -IT_0172 + 0.5*IT_0178
       + 0.5*IT_0179 + -IT_0180 + IT_0186 + IT_0187 + IT_0192 + IT_0193 + 
      -IT_0196 + -IT_0197 + -IT_0201 + -IT_0202 + -IT_0203 + IT_0207 + IT_0208 +
       -IT_0209 + -IT_0210 + IT_0211 + IT_0212 + -IT_0217 + -IT_0218 + -IT_0219 
      + 0.5*IT_0225 + 0.5*IT_0226 + IT_0229 + IT_0230 + -IT_0234 + -IT_0235 +
       IT_0236 + -IT_0241 + -IT_0242 + IT_0246 + IT_0247 + -IT_0251 + -IT_0252 +
       -IT_0256 + -IT_0257 + IT_0261 + IT_0262 + -IT_0265 + -IT_0266 + IT_0270 +
       IT_0271 + -IT_0272 + -IT_0277 + -IT_0278 + 0.5*IT_0281 + 0.5*IT_0282 +
       IT_0287 + IT_0288 + -IT_0291 + -IT_0292 + -IT_0295 + -IT_0296 + -IT_0299 
      + -IT_0300 + IT_0301 + IT_0304 + IT_0305 + -IT_0308 + -IT_0309 + (-0.5)
      *IT_0312 + (-0.5)*IT_0313 + -IT_0316 + -IT_0317 + IT_0320 + IT_0321 + (
      -0.5)*IT_0322 + IT_0325 + IT_0326 + IT_0329 + IT_0330 + IT_0333 + IT_0334 
      + -IT_0337 + -IT_0338 + (-0.5)*IT_0341 + (-0.5)*IT_0342 + IT_0345 +
       IT_0346 + -IT_0347 + IT_0350 + IT_0351 + -IT_0354 + -IT_0355 + IT_0358 +
       IT_0359 + IT_0362 + IT_0363 + IT_0366 + IT_0367 + -IT_0370 + -IT_0371 + 
      -IT_0374 + -IT_0375 + IT_0378 + IT_0379 + IT_0382 + IT_0383 + -IT_0386 + 
      -IT_0387 + IT_0388 + IT_0389;
    const complex_t IT_0397 = IT_0007*IT_0391*IT_0396;
    const complex_t IT_0398 = IT_0004*IT_0397;
    const complex_t IT_0399 = (-0.5)*IT_0398;
    return -IT_0395 + IT_0399;
}
} // End of namespace c9_nmfv
