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
#include "C9Z_G.h"
#include "common.h"

#include "params.h"
#include "libcomplexop.h"

namespace c9_nmfv {

complex_t C9Z_G(
        param_t const &param
        )
{
    clearcache();
    const real_t M_W = param.M_W;
    const real_t g_s = param.g_s;
    const real_t M_Z = param.M_Z;
    const real_t m_b = param.m_b;
    const real_t m_s = param.m_s;
    const real_t V_tb = param.V_tb;
    const real_t e_em = param.e_em;
    const real_t m_mu = param.m_mu;
    const real_t m_sG = param.m_sG;
    const real_t s_34 = param.s_34;
    const real_t m_sb_L = param.m_sb_L;
    const real_t m_sb_R = param.m_sb_R;
    const real_t m_sd_L = param.m_sd_L;
    const real_t m_sd_R = param.m_sd_R;
    const real_t m_ss_L = param.m_ss_L;
    const real_t m_ss_R = param.m_ss_R;
    const real_t theta_W = param.theta_W;
    const real_t reg_prop = param.reg_prop;
    const complex_t V_ts = param.V_ts;
    const complex_t U_sd_00 = param.U_sd_00;
    const complex_t U_sd_01 = param.U_sd_01;
    const complex_t U_sd_02 = param.U_sd_02;
    const complex_t U_sd_03 = param.U_sd_03;
    const complex_t U_sd_04 = param.U_sd_04;
    const complex_t U_sd_05 = param.U_sd_05;
    const complex_t U_sd_10 = param.U_sd_10;
    const complex_t U_sd_11 = param.U_sd_11;
    const complex_t U_sd_12 = param.U_sd_12;
    const complex_t U_sd_13 = param.U_sd_13;
    const complex_t U_sd_14 = param.U_sd_14;
    const complex_t U_sd_15 = param.U_sd_15;
    const complex_t U_sd_20 = param.U_sd_20;
    const complex_t U_sd_21 = param.U_sd_21;
    const complex_t U_sd_22 = param.U_sd_22;
    const complex_t U_sd_23 = param.U_sd_23;
    const complex_t U_sd_24 = param.U_sd_24;
    const complex_t U_sd_25 = param.U_sd_25;
    const complex_t U_sd_30 = param.U_sd_30;
    const complex_t U_sd_31 = param.U_sd_31;
    const complex_t U_sd_32 = param.U_sd_32;
    const complex_t U_sd_33 = param.U_sd_33;
    const complex_t U_sd_34 = param.U_sd_34;
    const complex_t U_sd_35 = param.U_sd_35;
    const complex_t U_sd_40 = param.U_sd_40;
    const complex_t U_sd_41 = param.U_sd_41;
    const complex_t U_sd_42 = param.U_sd_42;
    const complex_t U_sd_43 = param.U_sd_43;
    const complex_t U_sd_44 = param.U_sd_44;
    const complex_t U_sd_45 = param.U_sd_45;
    const complex_t U_sd_50 = param.U_sd_50;
    const complex_t U_sd_51 = param.U_sd_51;
    const complex_t U_sd_52 = param.U_sd_52;
    const complex_t U_sd_53 = param.U_sd_53;
    const complex_t U_sd_54 = param.U_sd_54;
    const complex_t U_sd_55 = param.U_sd_55;
    const complex_t IT_0000 = powq(M_W, 2);
    const complex_t IT_0001 = powq(V_tb, -1);
    const complex_t IT_0002 = cpowq(conjq(V_ts), -1);
    const complex_t IT_0003 = powq(e_em, -4);
    const complex_t IT_0004 = 9.86960440108936*IT_0000*IT_0001*IT_0002*IT_0003;
    const complex_t IT_0005 = cosq(theta_W);
    const complex_t IT_0006 = cpowq(IT_0005, -2);
    const complex_t IT_0007 = sinq(theta_W);
    const complex_t IT_0008 = tanq(theta_W);
    const complex_t IT_0009 = cpowq(IT_0008, 2);
    const complex_t IT_0010 = cpowq(1 + IT_0009, (-0.5));
    const complex_t IT_0011 = (complex_t{0, 1})*e_em*IT_0006*IT_0007*IT_0010;
    const complex_t IT_0012 = (complex_t{0, 1.4142135623731})*g_s*conjq
      (U_sd_20);
    const complex_t IT_0013 = (complex_t{0, 1.4142135623731})*g_s*U_sd_10;
    const complex_t IT_0014 = U_sd_20*conjq(U_sd_20);
    const complex_t IT_0015 = U_sd_10*conjq(U_sd_10);
    const complex_t IT_0016 = U_sd_00*conjq(U_sd_00);
    const complex_t IT_0017 = IT_0014 + IT_0015 + IT_0016;
    const complex_t IT_0018 = cpowq(IT_0007, -1);
    const complex_t IT_0019 = (complex_t{0, 1})*e_em*IT_0010*((-1.5)*IT_0017
      *IT_0018 + IT_0006*IT_0007*((-0.5)*IT_0017 + U_sd_30*conjq(U_sd_30) +
       U_sd_40*conjq(U_sd_40) + U_sd_50*conjq(U_sd_50)));
    const complex_t IT_0020 = (-0.666666666666667)*IT_0019;
    const complex_t IT_0021 = powq(m_sG, 2);
    const complex_t IT_0022 = powq(m_sd_L, 2);
    const complex_t IT_0023 = mty::lt::C0iC(9, 0, 0, 0, IT_0021, IT_0022,
       IT_0022, mty::lt::reg_int);
    const complex_t IT_0024 = IT_0020*IT_0023;
    const complex_t IT_0025 = IT_0012*IT_0013*IT_0024;
    const complex_t IT_0026 = 0.101321183642338*IT_0025;
    const complex_t IT_0027 = IT_0011*IT_0026;
    const complex_t IT_0028 = IT_0006*IT_0007*IT_0010;
    const complex_t IT_0029 = e_em*IT_0028;
    const complex_t IT_0030 = IT_0010*IT_0018;
    const complex_t IT_0031 = e_em*IT_0030;
    const complex_t IT_0032 = (complex_t{0, 1})*(IT_0029 + -IT_0031);
    const complex_t IT_0033 = 0.5*IT_0032;
    const complex_t IT_0034 = IT_0026*IT_0033;
    const complex_t IT_0035 = (complex_t{0, 1.4142135623731})*g_s*conjq
      (U_sd_50);
    const complex_t IT_0036 = (complex_t{0, 1.4142135623731})*g_s*U_sd_40;
    const complex_t IT_0037 = IT_0024*IT_0035*IT_0036;
    const complex_t IT_0038 = 0.101321183642338*IT_0037;
    const complex_t IT_0039 = IT_0011*IT_0038;
    const complex_t IT_0040 = IT_0033*IT_0038;
    const complex_t IT_0041 = powq(m_b, 2);
    const complex_t IT_0042 = powq(m_s, 2);
    const complex_t IT_0043 = cpowq(IT_0041 + -IT_0042 + reg_prop, -1);
    const complex_t IT_0044 = 0.101321183642338*m_s;
    const complex_t IT_0045 = (complex_t{0, 1})*(IT_0029 + 3*IT_0031);
    const complex_t IT_0046 = (-0.166666666666667)*IT_0045;
    const complex_t IT_0047 = mty::lt::B0iC(0, 0, IT_0021, IT_0022,
       mty::lt::reg_int);
    const complex_t IT_0048 = m_sG*IT_0047;
    const complex_t IT_0049 = IT_0012*IT_0036*IT_0046*IT_0048;
    const complex_t IT_0050 = IT_0043*IT_0044*IT_0049;
    const complex_t IT_0051 = IT_0011*IT_0050;
    const complex_t IT_0052 = IT_0033*IT_0050;
    const complex_t IT_0053 = mty::lt::B0iC(3, IT_0041, IT_0021, IT_0022,
       mty::lt::reg_int);
    const complex_t IT_0054 = m_b*IT_0053;
    const complex_t IT_0055 = IT_0035*IT_0036*IT_0046*IT_0054;
    const complex_t IT_0056 = IT_0043*IT_0044*IT_0055;
    const complex_t IT_0057 = IT_0011*IT_0056;
    const complex_t IT_0058 = IT_0033*IT_0056;
    const complex_t IT_0059 = (complex_t{0, 1.4142135623731})*g_s*conjq
      (U_sd_21);
    const complex_t IT_0060 = (complex_t{0, 1.4142135623731})*g_s*U_sd_41;
    const complex_t IT_0061 = powq(m_ss_L, 2);
    const complex_t IT_0062 = mty::lt::B0iC(0, 0, IT_0021, IT_0061,
       mty::lt::reg_int);
    const complex_t IT_0063 = m_sG*IT_0062;
    const complex_t IT_0064 = IT_0046*IT_0059*IT_0060*IT_0063;
    const complex_t IT_0065 = IT_0043*IT_0044*IT_0064;
    const complex_t IT_0066 = IT_0011*IT_0065;
    const complex_t IT_0067 = IT_0033*IT_0065;
    const complex_t IT_0068 = (complex_t{0, 1.4142135623731})*g_s*U_sd_11;
    const complex_t IT_0069 = mty::lt::B0iC(3, IT_0041, IT_0021, IT_0061,
       mty::lt::reg_int);
    const complex_t IT_0070 = IT_0041*IT_0069;
    const complex_t IT_0071 = IT_0046*IT_0059*IT_0068*IT_0070;
    const complex_t IT_0072 = 0.101321183642338*IT_0043*IT_0071;
    const complex_t IT_0073 = IT_0011*IT_0072;
    const complex_t IT_0074 = IT_0033*IT_0072;
    const complex_t IT_0075 = (complex_t{0, 1.4142135623731})*g_s*conjq
      (U_sd_51);
    const complex_t IT_0076 = m_b*IT_0069;
    const complex_t IT_0077 = IT_0046*IT_0060*IT_0075*IT_0076;
    const complex_t IT_0078 = IT_0043*IT_0044*IT_0077;
    const complex_t IT_0079 = IT_0011*IT_0078;
    const complex_t IT_0080 = IT_0033*IT_0078;
    const complex_t IT_0081 = (complex_t{0, 1.4142135623731})*g_s*conjq
      (U_sd_22);
    const complex_t IT_0082 = (complex_t{0, 1.4142135623731})*g_s*U_sd_42;
    const complex_t IT_0083 = powq(m_sb_L, 2);
    const complex_t IT_0084 = mty::lt::B0iC(0, 0, IT_0021, IT_0083,
       mty::lt::reg_int);
    const complex_t IT_0085 = m_sG*IT_0084;
    const complex_t IT_0086 = IT_0046*IT_0081*IT_0082*IT_0085;
    const complex_t IT_0087 = IT_0043*IT_0044*IT_0086;
    const complex_t IT_0088 = IT_0011*IT_0087;
    const complex_t IT_0089 = IT_0033*IT_0087;
    const complex_t IT_0090 = (complex_t{0, 1.4142135623731})*g_s*conjq
      (U_sd_53);
    const complex_t IT_0091 = (complex_t{0, 1.4142135623731})*g_s*U_sd_43;
    const complex_t IT_0092 = powq(m_sd_R, 2);
    const complex_t IT_0093 = mty::lt::B0iC(3, IT_0041, IT_0021, IT_0092,
       mty::lt::reg_int);
    const complex_t IT_0094 = m_b*IT_0093;
    const complex_t IT_0095 = IT_0046*IT_0090*IT_0091*IT_0094;
    const complex_t IT_0096 = IT_0043*IT_0044*IT_0095;
    const complex_t IT_0097 = IT_0011*IT_0096;
    const complex_t IT_0098 = IT_0033*IT_0096;
    const complex_t IT_0099 = (complex_t{0, 1.4142135623731})*g_s*conjq
      (U_sd_24);
    const complex_t IT_0100 = (complex_t{0, 1.4142135623731})*g_s*U_sd_44;
    const complex_t IT_0101 = powq(m_ss_R, 2);
    const complex_t IT_0102 = mty::lt::B0iC(0, 0, IT_0021, IT_0101,
       mty::lt::reg_int);
    const complex_t IT_0103 = m_sG*IT_0102;
    const complex_t IT_0104 = IT_0046*IT_0099*IT_0100*IT_0103;
    const complex_t IT_0105 = IT_0043*IT_0044*IT_0104;
    const complex_t IT_0106 = IT_0011*IT_0105;
    const complex_t IT_0107 = IT_0033*IT_0105;
    const complex_t IT_0108 = cpowq(IT_0041 + -IT_0042 + -reg_prop, -1);
    const complex_t IT_0109 = 0.101321183642338*m_b;
    const complex_t IT_0110 = (complex_t{0, 1.4142135623731})*g_s*conjq
      (U_sd_54);
    const complex_t IT_0111 = (complex_t{0, 1.4142135623731})*g_s*U_sd_14;
    const complex_t IT_0112 = IT_0046*IT_0103*IT_0110*IT_0111;
    const complex_t IT_0113 = IT_0108*IT_0109*IT_0112;
    const complex_t IT_0114 = IT_0033*IT_0113;
    const complex_t IT_0115 = (complex_t{0, 1.4142135623731})*g_s*conjq
      (U_sd_55);
    const complex_t IT_0116 = (complex_t{0, 1.4142135623731})*g_s*U_sd_15;
    const complex_t IT_0117 = powq(m_sb_R, 2);
    const complex_t IT_0118 = mty::lt::B0iC(0, 0, IT_0021, IT_0117,
       mty::lt::reg_int);
    const complex_t IT_0119 = m_sG*IT_0118;
    const complex_t IT_0120 = IT_0046*IT_0115*IT_0116*IT_0119;
    const complex_t IT_0121 = IT_0108*IT_0109*IT_0120;
    const complex_t IT_0122 = IT_0011*IT_0121;
    const complex_t IT_0123 = IT_0033*IT_0121;
    const complex_t IT_0124 = 0.333333333333333*IT_0011;
    const complex_t IT_0125 = IT_0115*IT_0116*IT_0119*IT_0124;
    const complex_t IT_0126 = IT_0043*IT_0044*IT_0125;
    const complex_t IT_0127 = IT_0033*IT_0126;
    const complex_t IT_0128 = (complex_t{0, 1.4142135623731})*g_s*conjq
      (U_sd_23);
    const complex_t IT_0129 = (complex_t{0, 1.4142135623731})*g_s*U_sd_13;
    const complex_t IT_0130 = U_sd_23*conjq(U_sd_23);
    const complex_t IT_0131 = U_sd_13*conjq(U_sd_13);
    const complex_t IT_0132 = U_sd_03*conjq(U_sd_03);
    const complex_t IT_0133 = IT_0130 + IT_0131 + IT_0132;
    const complex_t IT_0134 = (complex_t{0, 1})*e_em*IT_0010*((-1.5)*IT_0018
      *IT_0133 + IT_0006*IT_0007*((-0.5)*IT_0133 + U_sd_33*conjq(U_sd_33) +
       U_sd_43*conjq(U_sd_43) + U_sd_53*conjq(U_sd_53)));
    const complex_t IT_0135 = (-0.666666666666667)*IT_0134;
    const complex_t IT_0136 = mty::lt::C0iC(9, 0, 0, 0, IT_0021, IT_0092,
       IT_0092, mty::lt::reg_int);
    const complex_t IT_0137 = IT_0135*IT_0136;
    const complex_t IT_0138 = IT_0128*IT_0129*IT_0137;
    const complex_t IT_0139 = 0.101321183642338*IT_0138;
    const complex_t IT_0140 = IT_0033*IT_0139;
    const complex_t IT_0141 = (complex_t{0, 1.4142135623731})*g_s*U_sd_12;
    const complex_t IT_0142 = U_sd_22*conjq(U_sd_22);
    const complex_t IT_0143 = U_sd_12*conjq(U_sd_12);
    const complex_t IT_0144 = U_sd_02*conjq(U_sd_02);
    const complex_t IT_0145 = IT_0142 + IT_0143 + IT_0144;
    const complex_t IT_0146 = (complex_t{0, 1})*e_em*IT_0010*((-1.5)*IT_0018
      *IT_0145 + IT_0006*IT_0007*((-0.5)*IT_0145 + U_sd_32*conjq(U_sd_32) +
       U_sd_42*conjq(U_sd_42) + U_sd_52*conjq(U_sd_52)));
    const complex_t IT_0147 = (-0.666666666666667)*IT_0146;
    const complex_t IT_0148 = mty::lt::C0iC(9, 0, 0, 0, IT_0021, IT_0083,
       IT_0083, mty::lt::reg_int);
    const complex_t IT_0149 = IT_0147*IT_0148;
    const complex_t IT_0150 = IT_0081*IT_0141*IT_0149;
    const complex_t IT_0151 = 0.101321183642338*IT_0150;
    const complex_t IT_0152 = IT_0011*IT_0151;
    const complex_t IT_0153 = IT_0033*IT_0151;
    const complex_t IT_0154 = (complex_t{0, 1.4142135623731})*g_s*conjq
      (U_sd_52);
    const complex_t IT_0155 = IT_0082*IT_0149*IT_0154;
    const complex_t IT_0156 = 0.101321183642338*IT_0155;
    const complex_t IT_0157 = IT_0011*IT_0156;
    const complex_t IT_0158 = IT_0033*IT_0156;
    const complex_t IT_0159 = U_sd_21*conjq(U_sd_21);
    const complex_t IT_0160 = U_sd_11*conjq(U_sd_11);
    const complex_t IT_0161 = U_sd_01*conjq(U_sd_01);
    const complex_t IT_0162 = IT_0159 + IT_0160 + IT_0161;
    const complex_t IT_0163 = (complex_t{0, 1})*e_em*IT_0010*((-1.5)*IT_0018
      *IT_0162 + IT_0006*IT_0007*((-0.5)*IT_0162 + U_sd_31*conjq(U_sd_31) +
       U_sd_41*conjq(U_sd_41) + U_sd_51*conjq(U_sd_51)));
    const complex_t IT_0164 = (-0.666666666666667)*IT_0163;
    const complex_t IT_0165 = mty::lt::C0iC(9, 0, 0, 0, IT_0021, IT_0061,
       IT_0061, mty::lt::reg_int);
    const complex_t IT_0166 = IT_0164*IT_0165;
    const complex_t IT_0167 = IT_0059*IT_0068*IT_0166;
    const complex_t IT_0168 = 0.101321183642338*IT_0167;
    const complex_t IT_0169 = IT_0011*IT_0168;
    const complex_t IT_0170 = (complex_t{0, 1.4142135623731})*g_s*conjq
      (U_sd_25);
    const complex_t IT_0171 = (complex_t{0, 1.4142135623731})*g_s*U_sd_45;
    const complex_t IT_0172 = m_b*m_sG;
    const complex_t IT_0173 = IT_0118*IT_0172;
    const complex_t IT_0174 = IT_0124*IT_0170*IT_0171*IT_0173;
    const complex_t IT_0175 = 0.101321183642338*IT_0043*IT_0174;
    const complex_t IT_0176 = IT_0011*IT_0175;
    const complex_t IT_0177 = mty::lt::B0iC(0, 0, IT_0021, IT_0092,
       mty::lt::reg_int);
    const complex_t IT_0178 = m_sG*IT_0177;
    const complex_t IT_0179 = IT_0046*IT_0091*IT_0128*IT_0178;
    const complex_t IT_0180 = IT_0043*IT_0044*IT_0179;
    const complex_t IT_0181 = IT_0011*IT_0180;
    const complex_t IT_0182 = IT_0046*IT_0119*IT_0170*IT_0171;
    const complex_t IT_0183 = IT_0043*IT_0044*IT_0182;
    const complex_t IT_0184 = IT_0011*IT_0183;
    const complex_t IT_0185 = mty::lt::B0iC(3, IT_0041, IT_0021, IT_0117,
       mty::lt::reg_int);
    const complex_t IT_0186 = IT_0041*IT_0185;
    const complex_t IT_0187 = IT_0046*IT_0116*IT_0170*IT_0186;
    const complex_t IT_0188 = 0.101321183642338*IT_0043*IT_0187;
    const complex_t IT_0189 = IT_0033*IT_0188;
    const complex_t IT_0190 = m_b*IT_0185;
    const complex_t IT_0191 = IT_0046*IT_0115*IT_0171*IT_0190;
    const complex_t IT_0192 = IT_0043*IT_0044*IT_0191;
    const complex_t IT_0193 = IT_0011*IT_0192;
    const complex_t IT_0194 = IT_0033*IT_0192;
    const complex_t IT_0195 = IT_0047*IT_0172;
    const complex_t IT_0196 = IT_0013*IT_0035*IT_0046*IT_0195;
    const complex_t IT_0197 = 0.101321183642338*IT_0043*IT_0196;
    const complex_t IT_0198 = IT_0011*IT_0197;
    const complex_t IT_0199 = IT_0033*IT_0197;
    const complex_t IT_0200 = m_s*m_sG;
    const complex_t IT_0201 = IT_0047*IT_0200;
    const complex_t IT_0202 = IT_0012*IT_0036*IT_0046*IT_0201;
    const complex_t IT_0203 = 0.101321183642338*IT_0108*IT_0202;
    const complex_t IT_0204 = IT_0011*IT_0203;
    const complex_t IT_0205 = IT_0033*IT_0203;
    const complex_t IT_0206 = mty::lt::B0iC(3, IT_0042, IT_0021, IT_0022,
       mty::lt::reg_int);
    const complex_t IT_0207 = m_s*IT_0206;
    const complex_t IT_0208 = IT_0035*IT_0036*IT_0046*IT_0207;
    const complex_t IT_0209 = IT_0108*IT_0109*IT_0208;
    const complex_t IT_0210 = IT_0011*IT_0209;
    const complex_t IT_0211 = IT_0033*IT_0209;
    const complex_t IT_0212 = IT_0062*IT_0200;
    const complex_t IT_0213 = IT_0046*IT_0059*IT_0060*IT_0212;
    const complex_t IT_0214 = 0.101321183642338*IT_0108*IT_0213;
    const complex_t IT_0215 = IT_0011*IT_0214;
    const complex_t IT_0216 = IT_0062*IT_0172;
    const complex_t IT_0217 = IT_0059*IT_0060*IT_0124*IT_0216;
    const complex_t IT_0218 = 0.101321183642338*IT_0043*IT_0217;
    const complex_t IT_0219 = IT_0033*IT_0218;
    const complex_t IT_0220 = IT_0116*IT_0124*IT_0170*IT_0190;
    const complex_t IT_0221 = IT_0043*IT_0044*IT_0220;
    const complex_t IT_0222 = IT_0011*IT_0221;
    const complex_t IT_0223 = U_sd_20*conjq(U_sd_21);
    const complex_t IT_0224 = U_sd_10*conjq(U_sd_11);
    const complex_t IT_0225 = U_sd_00*conjq(U_sd_01);
    const complex_t IT_0226 = IT_0223 + IT_0224 + IT_0225;
    const complex_t IT_0227 = (complex_t{0, 1})*e_em*IT_0010*((-1.5)*IT_0018
      *IT_0226 + IT_0006*IT_0007*((-0.5)*IT_0226 + U_sd_30*conjq(U_sd_31) +
       U_sd_40*conjq(U_sd_41) + U_sd_50*conjq(U_sd_51)));
    const complex_t IT_0228 = (-0.666666666666667)*IT_0227;
    const complex_t IT_0229 = mty::lt::C0iC(9, 0, 0, 0, IT_0021, IT_0022,
       IT_0061, mty::lt::reg_int);
    const complex_t IT_0230 = IT_0228*IT_0229;
    const complex_t IT_0231 = IT_0035*IT_0060*IT_0230;
    const complex_t IT_0232 = 0.101321183642338*IT_0231;
    const complex_t IT_0233 = IT_0033*IT_0232;
    const complex_t IT_0234 = mty::lt::B0iC(3, IT_0041, IT_0021, IT_0101,
       mty::lt::reg_int);
    const complex_t IT_0235 = m_b*IT_0234;
    const complex_t IT_0236 = IT_0046*IT_0100*IT_0110*IT_0235;
    const complex_t IT_0237 = IT_0043*IT_0044*IT_0236;
    const complex_t IT_0238 = IT_0011*IT_0237;
    const complex_t IT_0239 = IT_0046*IT_0068*IT_0075*IT_0216;
    const complex_t IT_0240 = 0.101321183642338*IT_0043*IT_0239;
    const complex_t IT_0241 = IT_0011*IT_0240;
    const complex_t IT_0242 = IT_0033*IT_0240;
    const complex_t IT_0243 = IT_0084*IT_0172;
    const complex_t IT_0244 = IT_0046*IT_0141*IT_0154*IT_0243;
    const complex_t IT_0245 = 0.101321183642338*IT_0043*IT_0244;
    const complex_t IT_0246 = IT_0011*IT_0245;
    const complex_t IT_0247 = IT_0042*IT_0206;
    const complex_t IT_0248 = IT_0012*IT_0013*IT_0046*IT_0247;
    const complex_t IT_0249 = 0.101321183642338*IT_0108*IT_0248;
    const complex_t IT_0250 = IT_0011*IT_0249;
    const complex_t IT_0251 = IT_0033*IT_0249;
    const complex_t IT_0252 = IT_0013*IT_0035*IT_0046*IT_0048;
    const complex_t IT_0253 = IT_0108*IT_0109*IT_0252;
    const complex_t IT_0254 = IT_0011*IT_0253;
    const complex_t IT_0255 = IT_0102*IT_0172;
    const complex_t IT_0256 = IT_0099*IT_0100*IT_0124*IT_0255;
    const complex_t IT_0257 = 0.101321183642338*IT_0043*IT_0256;
    const complex_t IT_0258 = IT_0033*IT_0257;
    const complex_t IT_0259 = IT_0115*IT_0124*IT_0171*IT_0186;
    const complex_t IT_0260 = 0.101321183642338*IT_0043*IT_0259;
    const complex_t IT_0261 = IT_0033*IT_0260;
    const complex_t IT_0262 = conjq(U_sd_21)*U_sd_25;
    const complex_t IT_0263 = conjq(U_sd_11)*U_sd_15;
    const complex_t IT_0264 = conjq(U_sd_01)*U_sd_05;
    const complex_t IT_0265 = IT_0262 + IT_0263 + IT_0264;
    const complex_t IT_0266 = (complex_t{0, 1})*e_em*IT_0010*((-1.5)*IT_0018
      *IT_0265 + IT_0006*IT_0007*((-0.5)*IT_0265 + conjq(U_sd_31)*U_sd_35 +
       conjq(U_sd_41)*U_sd_45 + conjq(U_sd_51)*U_sd_55));
    const complex_t IT_0267 = (-0.666666666666667)*IT_0266;
    const complex_t IT_0268 = mty::lt::C0iC(9, 0, 0, 0, IT_0021, IT_0117,
       IT_0061, mty::lt::reg_int);
    const complex_t IT_0269 = IT_0267*IT_0268;
    const complex_t IT_0270 = IT_0060*IT_0115*IT_0269;
    const complex_t IT_0271 = 0.101321183642338*IT_0270;
    const complex_t IT_0272 = IT_0011*IT_0271;
    const complex_t IT_0273 = IT_0033*IT_0271;
    const complex_t IT_0274 = IT_0068*IT_0170*IT_0269;
    const complex_t IT_0275 = 0.101321183642338*IT_0274;
    const complex_t IT_0276 = IT_0011*IT_0275;
    const complex_t IT_0277 = IT_0033*IT_0275;
    const complex_t IT_0278 = U_sd_21*conjq(U_sd_24);
    const complex_t IT_0279 = U_sd_11*conjq(U_sd_14);
    const complex_t IT_0280 = U_sd_01*conjq(U_sd_04);
    const complex_t IT_0281 = IT_0278 + IT_0279 + IT_0280;
    const complex_t IT_0282 = (complex_t{0, 1})*e_em*IT_0010*((-1.5)*IT_0018
      *IT_0281 + IT_0006*IT_0007*((-0.5)*IT_0281 + U_sd_31*conjq(U_sd_34) +
       U_sd_41*conjq(U_sd_44) + U_sd_51*conjq(U_sd_54)));
    const complex_t IT_0283 = (-0.666666666666667)*IT_0282;
    const complex_t IT_0284 = mty::lt::C0iC(9, 0, 0, 0, IT_0021, IT_0061,
       IT_0101, mty::lt::reg_int);
    const complex_t IT_0285 = IT_0283*IT_0284;
    const complex_t IT_0286 = IT_0059*IT_0111*IT_0285;
    const complex_t IT_0287 = 0.101321183642338*IT_0286;
    const complex_t IT_0288 = IT_0011*IT_0287;
    const complex_t IT_0289 = IT_0033*IT_0287;
    const complex_t IT_0290 = IT_0075*IT_0100*IT_0285;
    const complex_t IT_0291 = 0.101321183642338*IT_0290;
    const complex_t IT_0292 = IT_0011*IT_0291;
    const complex_t IT_0293 = IT_0033*IT_0291;
    const complex_t IT_0294 = conjq(U_sd_21)*U_sd_22;
    const complex_t IT_0295 = conjq(U_sd_11)*U_sd_12;
    const complex_t IT_0296 = conjq(U_sd_01)*U_sd_02;
    const complex_t IT_0297 = IT_0294 + IT_0295 + IT_0296;
    const complex_t IT_0298 = (complex_t{0, 1})*e_em*IT_0010*((-1.5)*IT_0018
      *IT_0297 + IT_0006*IT_0007*((-0.5)*IT_0297 + conjq(U_sd_31)*U_sd_32 +
       conjq(U_sd_41)*U_sd_42 + conjq(U_sd_51)*U_sd_52));
    const complex_t IT_0299 = (-0.666666666666667)*IT_0298;
    const complex_t IT_0300 = mty::lt::C0iC(9, 0, 0, 0, IT_0021, IT_0083,
       IT_0061, mty::lt::reg_int);
    const complex_t IT_0301 = IT_0299*IT_0300;
    const complex_t IT_0302 = IT_0060*IT_0154*IT_0301;
    const complex_t IT_0303 = 0.101321183642338*IT_0302;
    const complex_t IT_0304 = IT_0011*IT_0303;
    const complex_t IT_0305 = IT_0033*IT_0303;
    const complex_t IT_0306 = IT_0068*IT_0081*IT_0301;
    const complex_t IT_0307 = 0.101321183642338*IT_0306;
    const complex_t IT_0308 = IT_0011*IT_0307;
    const complex_t IT_0309 = IT_0033*IT_0307;
    const complex_t IT_0310 = U_sd_22*conjq(U_sd_24);
    const complex_t IT_0311 = U_sd_12*conjq(U_sd_14);
    const complex_t IT_0312 = U_sd_02*conjq(U_sd_04);
    const complex_t IT_0313 = IT_0310 + IT_0311 + IT_0312;
    const complex_t IT_0314 = (complex_t{0, 1})*e_em*IT_0010*((-1.5)*IT_0018
      *IT_0313 + IT_0006*IT_0007*((-0.5)*IT_0313 + U_sd_32*conjq(U_sd_34) +
       U_sd_42*conjq(U_sd_44) + U_sd_52*conjq(U_sd_54)));
    const complex_t IT_0315 = (-0.666666666666667)*IT_0314;
    const complex_t IT_0316 = mty::lt::C0iC(9, 0, 0, 0, IT_0021, IT_0083,
       IT_0101, mty::lt::reg_int);
    const complex_t IT_0317 = IT_0315*IT_0316;
    const complex_t IT_0318 = IT_0081*IT_0111*IT_0317;
    const complex_t IT_0319 = 0.101321183642338*IT_0318;
    const complex_t IT_0320 = IT_0011*IT_0319;
    const complex_t IT_0321 = IT_0033*IT_0319;
    const complex_t IT_0322 = IT_0100*IT_0154*IT_0317;
    const complex_t IT_0323 = 0.101321183642338*IT_0322;
    const complex_t IT_0324 = IT_0011*IT_0323;
    const complex_t IT_0325 = IT_0033*IT_0323;
    const complex_t IT_0326 = IT_0011*IT_0139;
    const complex_t IT_0327 = IT_0090*IT_0091*IT_0137;
    const complex_t IT_0328 = 0.101321183642338*IT_0327;
    const complex_t IT_0329 = IT_0011*IT_0328;
    const complex_t IT_0330 = IT_0033*IT_0328;
    const complex_t IT_0331 = U_sd_20*conjq(U_sd_24);
    const complex_t IT_0332 = U_sd_10*conjq(U_sd_14);
    const complex_t IT_0333 = U_sd_00*conjq(U_sd_04);
    const complex_t IT_0334 = IT_0331 + IT_0332 + IT_0333;
    const complex_t IT_0335 = (complex_t{0, 1})*e_em*IT_0010*((-1.5)*IT_0018
      *IT_0334 + IT_0006*IT_0007*((-0.5)*IT_0334 + U_sd_30*conjq(U_sd_34) +
       U_sd_40*conjq(U_sd_44) + U_sd_50*conjq(U_sd_54)));
    const complex_t IT_0336 = (-0.666666666666667)*IT_0335;
    const complex_t IT_0337 = mty::lt::C0iC(9, 0, 0, 0, IT_0021, IT_0022,
       IT_0101, mty::lt::reg_int);
    const complex_t IT_0338 = IT_0336*IT_0337;
    const complex_t IT_0339 = IT_0012*IT_0111*IT_0338;
    const complex_t IT_0340 = 0.101321183642338*IT_0339;
    const complex_t IT_0341 = IT_0011*IT_0340;
    const complex_t IT_0342 = IT_0033*IT_0340;
    const complex_t IT_0343 = IT_0035*IT_0100*IT_0338;
    const complex_t IT_0344 = 0.101321183642338*IT_0343;
    const complex_t IT_0345 = IT_0011*IT_0344;
    const complex_t IT_0346 = IT_0033*IT_0344;
    const complex_t IT_0347 = conjq(U_sd_20)*U_sd_25;
    const complex_t IT_0348 = conjq(U_sd_10)*U_sd_15;
    const complex_t IT_0349 = conjq(U_sd_00)*U_sd_05;
    const complex_t IT_0350 = IT_0347 + IT_0348 + IT_0349;
    const complex_t IT_0351 = (complex_t{0, 1})*e_em*IT_0010*((-1.5)*IT_0018
      *IT_0350 + IT_0006*IT_0007*((-0.5)*IT_0350 + conjq(U_sd_30)*U_sd_35 +
       conjq(U_sd_40)*U_sd_45 + conjq(U_sd_50)*U_sd_55));
    const complex_t IT_0352 = (-0.666666666666667)*IT_0351;
    const complex_t IT_0353 = mty::lt::C0iC(9, 0, 0, 0, IT_0021, IT_0117,
       IT_0022, mty::lt::reg_int);
    const complex_t IT_0354 = IT_0352*IT_0353;
    const complex_t IT_0355 = IT_0036*IT_0115*IT_0354;
    const complex_t IT_0356 = 0.101321183642338*IT_0355;
    const complex_t IT_0357 = IT_0011*IT_0356;
    const complex_t IT_0358 = IT_0033*IT_0356;
    const complex_t IT_0359 = IT_0013*IT_0170*IT_0354;
    const complex_t IT_0360 = 0.101321183642338*IT_0359;
    const complex_t IT_0361 = IT_0011*IT_0360;
    const complex_t IT_0362 = IT_0033*IT_0360;
    const complex_t IT_0363 = conjq(U_sd_20)*U_sd_22;
    const complex_t IT_0364 = conjq(U_sd_10)*U_sd_12;
    const complex_t IT_0365 = conjq(U_sd_00)*U_sd_02;
    const complex_t IT_0366 = IT_0363 + IT_0364 + IT_0365;
    const complex_t IT_0367 = (complex_t{0, 1})*e_em*IT_0010*((-1.5)*IT_0018
      *IT_0366 + IT_0006*IT_0007*((-0.5)*IT_0366 + conjq(U_sd_30)*U_sd_32 +
       conjq(U_sd_40)*U_sd_42 + conjq(U_sd_50)*U_sd_52));
    const complex_t IT_0368 = (-0.666666666666667)*IT_0367;
    const complex_t IT_0369 = mty::lt::C0iC(9, 0, 0, 0, IT_0021, IT_0083,
       IT_0022, mty::lt::reg_int);
    const complex_t IT_0370 = IT_0368*IT_0369;
    const complex_t IT_0371 = IT_0036*IT_0154*IT_0370;
    const complex_t IT_0372 = 0.101321183642338*IT_0371;
    const complex_t IT_0373 = IT_0011*IT_0372;
    const complex_t IT_0374 = IT_0033*IT_0372;
    const complex_t IT_0375 = IT_0013*IT_0081*IT_0370;
    const complex_t IT_0376 = 0.101321183642338*IT_0375;
    const complex_t IT_0377 = IT_0011*IT_0376;
    const complex_t IT_0378 = IT_0033*IT_0376;
    const complex_t IT_0379 = IT_0012*IT_0068*IT_0230;
    const complex_t IT_0380 = 0.101321183642338*IT_0379;
    const complex_t IT_0381 = IT_0011*IT_0380;
    const complex_t IT_0382 = IT_0033*IT_0380;
    const complex_t IT_0383 = IT_0011*IT_0232;
    const complex_t IT_0384 = IT_0033*IT_0168;
    const complex_t IT_0385 = IT_0060*IT_0075*IT_0166;
    const complex_t IT_0386 = 0.101321183642338*IT_0385;
    const complex_t IT_0387 = IT_0011*IT_0386;
    const complex_t IT_0388 = IT_0033*IT_0386;
    const complex_t IT_0389 = IT_0041*IT_0053;
    const complex_t IT_0390 = IT_0012*IT_0013*IT_0046*IT_0389;
    const complex_t IT_0391 = 0.101321183642338*IT_0043*IT_0390;
    const complex_t IT_0392 = IT_0011*IT_0391;
    const complex_t IT_0393 = IT_0033*IT_0391;
    const complex_t IT_0394 = mty::lt::B0iC(3, IT_0041, IT_0021, IT_0083,
       mty::lt::reg_int);
    const complex_t IT_0395 = IT_0041*IT_0394;
    const complex_t IT_0396 = IT_0046*IT_0081*IT_0141*IT_0395;
    const complex_t IT_0397 = 0.101321183642338*IT_0043*IT_0396;
    const complex_t IT_0398 = IT_0011*IT_0397;
    const complex_t IT_0399 = IT_0033*IT_0397;
    const complex_t IT_0400 = m_b*IT_0394;
    const complex_t IT_0401 = IT_0046*IT_0082*IT_0154*IT_0400;
    const complex_t IT_0402 = IT_0043*IT_0044*IT_0401;
    const complex_t IT_0403 = IT_0011*IT_0402;
    const complex_t IT_0404 = IT_0033*IT_0402;
    const complex_t IT_0405 = IT_0033*IT_0180;
    const complex_t IT_0406 = IT_0041*IT_0093;
    const complex_t IT_0407 = IT_0046*IT_0128*IT_0129*IT_0406;
    const complex_t IT_0408 = 0.101321183642338*IT_0043*IT_0407;
    const complex_t IT_0409 = IT_0011*IT_0408;
    const complex_t IT_0410 = IT_0033*IT_0408;
    const complex_t IT_0411 = IT_0041*IT_0234;
    const complex_t IT_0412 = IT_0046*IT_0099*IT_0111*IT_0411;
    const complex_t IT_0413 = 0.101321183642338*IT_0043*IT_0412;
    const complex_t IT_0414 = IT_0011*IT_0413;
    const complex_t IT_0415 = IT_0033*IT_0413;
    const complex_t IT_0416 = IT_0033*IT_0237;
    const complex_t IT_0417 = IT_0033*IT_0183;
    const complex_t IT_0418 = IT_0011*IT_0188;
    const complex_t IT_0419 = IT_0033*IT_0245;
    const complex_t IT_0420 = IT_0172*IT_0177;
    const complex_t IT_0421 = IT_0046*IT_0090*IT_0129*IT_0420;
    const complex_t IT_0422 = 0.101321183642338*IT_0043*IT_0421;
    const complex_t IT_0423 = IT_0011*IT_0422;
    const complex_t IT_0424 = IT_0033*IT_0422;
    const complex_t IT_0425 = IT_0046*IT_0110*IT_0111*IT_0255;
    const complex_t IT_0426 = 0.101321183642338*IT_0043*IT_0425;
    const complex_t IT_0427 = IT_0011*IT_0426;
    const complex_t IT_0428 = IT_0033*IT_0426;
    const complex_t IT_0429 = IT_0046*IT_0115*IT_0116*IT_0173;
    const complex_t IT_0430 = 0.101321183642338*IT_0043*IT_0429;
    const complex_t IT_0431 = IT_0011*IT_0430;
    const complex_t IT_0432 = IT_0033*IT_0430;
    const complex_t IT_0433 = IT_0033*IT_0214;
    const complex_t IT_0434 = mty::lt::B0iC(3, IT_0042, IT_0021, IT_0061,
       mty::lt::reg_int);
    const complex_t IT_0435 = IT_0042*IT_0434;
    const complex_t IT_0436 = IT_0046*IT_0059*IT_0068*IT_0435;
    const complex_t IT_0437 = 0.101321183642338*IT_0108*IT_0436;
    const complex_t IT_0438 = IT_0011*IT_0437;
    const complex_t IT_0439 = IT_0033*IT_0437;
    const complex_t IT_0440 = m_s*IT_0434;
    const complex_t IT_0441 = IT_0046*IT_0060*IT_0075*IT_0440;
    const complex_t IT_0442 = IT_0108*IT_0109*IT_0441;
    const complex_t IT_0443 = IT_0011*IT_0442;
    const complex_t IT_0444 = IT_0033*IT_0442;
    const complex_t IT_0445 = IT_0084*IT_0200;
    const complex_t IT_0446 = IT_0046*IT_0081*IT_0082*IT_0445;
    const complex_t IT_0447 = 0.101321183642338*IT_0108*IT_0446;
    const complex_t IT_0448 = IT_0011*IT_0447;
    const complex_t IT_0449 = IT_0033*IT_0447;
    const complex_t IT_0450 = mty::lt::B0iC(3, IT_0042, IT_0021, IT_0083,
       mty::lt::reg_int);
    const complex_t IT_0451 = IT_0042*IT_0450;
    const complex_t IT_0452 = IT_0046*IT_0081*IT_0141*IT_0451;
    const complex_t IT_0453 = 0.101321183642338*IT_0108*IT_0452;
    const complex_t IT_0454 = IT_0011*IT_0453;
    const complex_t IT_0455 = IT_0033*IT_0453;
    const complex_t IT_0456 = m_s*IT_0450;
    const complex_t IT_0457 = IT_0046*IT_0082*IT_0154*IT_0456;
    const complex_t IT_0458 = IT_0108*IT_0109*IT_0457;
    const complex_t IT_0459 = IT_0011*IT_0458;
    const complex_t IT_0460 = IT_0033*IT_0458;
    const complex_t IT_0461 = IT_0177*IT_0200;
    const complex_t IT_0462 = IT_0046*IT_0091*IT_0128*IT_0461;
    const complex_t IT_0463 = 0.101321183642338*IT_0108*IT_0462;
    const complex_t IT_0464 = IT_0011*IT_0463;
    const complex_t IT_0465 = IT_0033*IT_0463;
    const complex_t IT_0466 = mty::lt::B0iC(3, IT_0042, IT_0021, IT_0092,
       mty::lt::reg_int);
    const complex_t IT_0467 = IT_0042*IT_0466;
    const complex_t IT_0468 = IT_0046*IT_0128*IT_0129*IT_0467;
    const complex_t IT_0469 = 0.101321183642338*IT_0108*IT_0468;
    const complex_t IT_0470 = IT_0011*IT_0469;
    const complex_t IT_0471 = IT_0033*IT_0469;
    const complex_t IT_0472 = m_s*IT_0466;
    const complex_t IT_0473 = IT_0046*IT_0090*IT_0091*IT_0472;
    const complex_t IT_0474 = IT_0108*IT_0109*IT_0473;
    const complex_t IT_0475 = IT_0011*IT_0474;
    const complex_t IT_0476 = IT_0033*IT_0474;
    const complex_t IT_0477 = IT_0102*IT_0200;
    const complex_t IT_0478 = IT_0046*IT_0099*IT_0100*IT_0477;
    const complex_t IT_0479 = 0.101321183642338*IT_0108*IT_0478;
    const complex_t IT_0480 = IT_0011*IT_0479;
    const complex_t IT_0481 = IT_0033*IT_0479;
    const complex_t IT_0482 = mty::lt::B0iC(3, IT_0042, IT_0021, IT_0101,
       mty::lt::reg_int);
    const complex_t IT_0483 = IT_0042*IT_0482;
    const complex_t IT_0484 = IT_0046*IT_0099*IT_0111*IT_0483;
    const complex_t IT_0485 = 0.101321183642338*IT_0108*IT_0484;
    const complex_t IT_0486 = IT_0011*IT_0485;
    const complex_t IT_0487 = IT_0033*IT_0485;
    const complex_t IT_0488 = m_s*IT_0482;
    const complex_t IT_0489 = IT_0046*IT_0100*IT_0110*IT_0488;
    const complex_t IT_0490 = IT_0108*IT_0109*IT_0489;
    const complex_t IT_0491 = IT_0011*IT_0490;
    const complex_t IT_0492 = IT_0033*IT_0490;
    const complex_t IT_0493 = IT_0118*IT_0200;
    const complex_t IT_0494 = IT_0046*IT_0170*IT_0171*IT_0493;
    const complex_t IT_0495 = 0.101321183642338*IT_0108*IT_0494;
    const complex_t IT_0496 = IT_0011*IT_0495;
    const complex_t IT_0497 = IT_0033*IT_0495;
    const complex_t IT_0498 = mty::lt::B0iC(3, IT_0042, IT_0021, IT_0117,
       mty::lt::reg_int);
    const complex_t IT_0499 = IT_0042*IT_0498;
    const complex_t IT_0500 = IT_0046*IT_0116*IT_0170*IT_0499;
    const complex_t IT_0501 = 0.101321183642338*IT_0108*IT_0500;
    const complex_t IT_0502 = IT_0011*IT_0501;
    const complex_t IT_0503 = IT_0033*IT_0501;
    const complex_t IT_0504 = m_s*IT_0498;
    const complex_t IT_0505 = IT_0046*IT_0115*IT_0171*IT_0504;
    const complex_t IT_0506 = IT_0108*IT_0109*IT_0505;
    const complex_t IT_0507 = IT_0011*IT_0506;
    const complex_t IT_0508 = IT_0033*IT_0506;
    const complex_t IT_0509 = IT_0033*IT_0253;
    const complex_t IT_0510 = IT_0046*IT_0063*IT_0068*IT_0075;
    const complex_t IT_0511 = IT_0108*IT_0109*IT_0510;
    const complex_t IT_0512 = IT_0011*IT_0511;
    const complex_t IT_0513 = IT_0033*IT_0511;
    const complex_t IT_0514 = IT_0046*IT_0085*IT_0141*IT_0154;
    const complex_t IT_0515 = IT_0108*IT_0109*IT_0514;
    const complex_t IT_0516 = IT_0011*IT_0515;
    const complex_t IT_0517 = IT_0033*IT_0515;
    const complex_t IT_0518 = IT_0046*IT_0090*IT_0129*IT_0178;
    const complex_t IT_0519 = IT_0108*IT_0109*IT_0518;
    const complex_t IT_0520 = IT_0011*IT_0519;
    const complex_t IT_0521 = IT_0033*IT_0519;
    const complex_t IT_0522 = IT_0011*IT_0113;
    const complex_t IT_0523 = conjq(U_sd_20)*U_sd_21;
    const complex_t IT_0524 = conjq(U_sd_10)*U_sd_11;
    const complex_t IT_0525 = conjq(U_sd_00)*U_sd_01;
    const complex_t IT_0526 = IT_0523 + IT_0524 + IT_0525;
    const complex_t IT_0527 = (complex_t{0, 1})*e_em*IT_0010*((-1.5)*IT_0018
      *IT_0526 + IT_0006*IT_0007*((-0.5)*IT_0526 + conjq(U_sd_30)*U_sd_31 +
       conjq(U_sd_40)*U_sd_41 + conjq(U_sd_50)*U_sd_51));
    const complex_t IT_0528 = (-0.666666666666667)*IT_0527;
    const complex_t IT_0529 = IT_0229*IT_0528;
    const complex_t IT_0530 = IT_0036*IT_0075*IT_0529;
    const complex_t IT_0531 = 0.101321183642338*IT_0530;
    const complex_t IT_0532 = IT_0011*IT_0531;
    const complex_t IT_0533 = IT_0033*IT_0531;
    const complex_t IT_0534 = IT_0013*IT_0059*IT_0529;
    const complex_t IT_0535 = 0.101321183642338*IT_0534;
    const complex_t IT_0536 = IT_0011*IT_0535;
    const complex_t IT_0537 = IT_0033*IT_0535;
    const complex_t IT_0538 = U_sd_20*conjq(U_sd_22);
    const complex_t IT_0539 = U_sd_10*conjq(U_sd_12);
    const complex_t IT_0540 = U_sd_00*conjq(U_sd_02);
    const complex_t IT_0541 = IT_0538 + IT_0539 + IT_0540;
    const complex_t IT_0542 = (complex_t{0, 1})*e_em*IT_0010*((-1.5)*IT_0018
      *IT_0541 + IT_0006*IT_0007*((-0.5)*IT_0541 + U_sd_30*conjq(U_sd_32) +
       U_sd_40*conjq(U_sd_42) + U_sd_50*conjq(U_sd_52)));
    const complex_t IT_0543 = (-0.666666666666667)*IT_0542;
    const complex_t IT_0544 = IT_0369*IT_0543;
    const complex_t IT_0545 = IT_0012*IT_0141*IT_0544;
    const complex_t IT_0546 = 0.101321183642338*IT_0545;
    const complex_t IT_0547 = IT_0011*IT_0546;
    const complex_t IT_0548 = IT_0033*IT_0546;
    const complex_t IT_0549 = IT_0035*IT_0082*IT_0544;
    const complex_t IT_0550 = 0.101321183642338*IT_0549;
    const complex_t IT_0551 = IT_0011*IT_0550;
    const complex_t IT_0552 = IT_0033*IT_0550;
    const complex_t IT_0553 = U_sd_21*conjq(U_sd_22);
    const complex_t IT_0554 = U_sd_11*conjq(U_sd_12);
    const complex_t IT_0555 = U_sd_01*conjq(U_sd_02);
    const complex_t IT_0556 = IT_0553 + IT_0554 + IT_0555;
    const complex_t IT_0557 = (complex_t{0, 1})*e_em*IT_0010*((-1.5)*IT_0018
      *IT_0556 + IT_0006*IT_0007*((-0.5)*IT_0556 + U_sd_31*conjq(U_sd_32) +
       U_sd_41*conjq(U_sd_42) + U_sd_51*conjq(U_sd_52)));
    const complex_t IT_0558 = (-0.666666666666667)*IT_0557;
    const complex_t IT_0559 = IT_0300*IT_0558;
    const complex_t IT_0560 = IT_0059*IT_0141*IT_0559;
    const complex_t IT_0561 = 0.101321183642338*IT_0560;
    const complex_t IT_0562 = IT_0011*IT_0561;
    const complex_t IT_0563 = IT_0033*IT_0561;
    const complex_t IT_0564 = IT_0075*IT_0082*IT_0559;
    const complex_t IT_0565 = 0.101321183642338*IT_0564;
    const complex_t IT_0566 = IT_0011*IT_0565;
    const complex_t IT_0567 = IT_0033*IT_0565;
    const complex_t IT_0568 = U_sd_20*conjq(U_sd_23);
    const complex_t IT_0569 = U_sd_10*conjq(U_sd_13);
    const complex_t IT_0570 = U_sd_00*conjq(U_sd_03);
    const complex_t IT_0571 = IT_0568 + IT_0569 + IT_0570;
    const complex_t IT_0572 = (complex_t{0, 1})*e_em*IT_0010*((-1.5)*IT_0018
      *IT_0571 + IT_0006*IT_0007*((-0.5)*IT_0571 + U_sd_30*conjq(U_sd_33) +
       U_sd_40*conjq(U_sd_43) + U_sd_50*conjq(U_sd_53)));
    const complex_t IT_0573 = (-0.666666666666667)*IT_0572;
    const complex_t IT_0574 = mty::lt::C0iC(9, 0, 0, 0, IT_0021, IT_0022,
       IT_0092, mty::lt::reg_int);
    const complex_t IT_0575 = IT_0573*IT_0574;
    const complex_t IT_0576 = IT_0012*IT_0129*IT_0575;
    const complex_t IT_0577 = 0.101321183642338*IT_0576;
    const complex_t IT_0578 = IT_0011*IT_0577;
    const complex_t IT_0579 = IT_0033*IT_0577;
    const complex_t IT_0580 = IT_0035*IT_0091*IT_0575;
    const complex_t IT_0581 = 0.101321183642338*IT_0580;
    const complex_t IT_0582 = IT_0011*IT_0581;
    const complex_t IT_0583 = IT_0033*IT_0581;
    const complex_t IT_0584 = conjq(U_sd_20)*U_sd_23;
    const complex_t IT_0585 = conjq(U_sd_10)*U_sd_13;
    const complex_t IT_0586 = conjq(U_sd_00)*U_sd_03;
    const complex_t IT_0587 = IT_0584 + IT_0585 + IT_0586;
    const complex_t IT_0588 = (complex_t{0, 1})*e_em*IT_0010*((-1.5)*IT_0018
      *IT_0587 + IT_0006*IT_0007*((-0.5)*IT_0587 + conjq(U_sd_30)*U_sd_33 +
       conjq(U_sd_40)*U_sd_43 + conjq(U_sd_50)*U_sd_53));
    const complex_t IT_0589 = (-0.666666666666667)*IT_0588;
    const complex_t IT_0590 = IT_0574*IT_0589;
    const complex_t IT_0591 = IT_0036*IT_0090*IT_0590;
    const complex_t IT_0592 = 0.101321183642338*IT_0591;
    const complex_t IT_0593 = IT_0011*IT_0592;
    const complex_t IT_0594 = IT_0033*IT_0592;
    const complex_t IT_0595 = IT_0013*IT_0128*IT_0590;
    const complex_t IT_0596 = 0.101321183642338*IT_0595;
    const complex_t IT_0597 = IT_0011*IT_0596;
    const complex_t IT_0598 = IT_0033*IT_0596;
    const complex_t IT_0599 = U_sd_21*conjq(U_sd_23);
    const complex_t IT_0600 = U_sd_11*conjq(U_sd_13);
    const complex_t IT_0601 = U_sd_01*conjq(U_sd_03);
    const complex_t IT_0602 = IT_0599 + IT_0600 + IT_0601;
    const complex_t IT_0603 = (complex_t{0, 1})*e_em*IT_0010*((-1.5)*IT_0018
      *IT_0602 + IT_0006*IT_0007*((-0.5)*IT_0602 + U_sd_31*conjq(U_sd_33) +
       U_sd_41*conjq(U_sd_43) + U_sd_51*conjq(U_sd_53)));
    const complex_t IT_0604 = (-0.666666666666667)*IT_0603;
    const complex_t IT_0605 = mty::lt::C0iC(9, 0, 0, 0, IT_0021, IT_0092,
       IT_0061, mty::lt::reg_int);
    const complex_t IT_0606 = IT_0604*IT_0605;
    const complex_t IT_0607 = IT_0059*IT_0129*IT_0606;
    const complex_t IT_0608 = 0.101321183642338*IT_0607;
    const complex_t IT_0609 = IT_0011*IT_0608;
    const complex_t IT_0610 = IT_0033*IT_0608;
    const complex_t IT_0611 = IT_0075*IT_0091*IT_0606;
    const complex_t IT_0612 = 0.101321183642338*IT_0611;
    const complex_t IT_0613 = IT_0011*IT_0612;
    const complex_t IT_0614 = IT_0033*IT_0612;
    const complex_t IT_0615 = conjq(U_sd_21)*U_sd_23;
    const complex_t IT_0616 = conjq(U_sd_11)*U_sd_13;
    const complex_t IT_0617 = conjq(U_sd_01)*U_sd_03;
    const complex_t IT_0618 = IT_0615 + IT_0616 + IT_0617;
    const complex_t IT_0619 = (complex_t{0, 1})*e_em*IT_0010*((-1.5)*IT_0018
      *IT_0618 + IT_0006*IT_0007*((-0.5)*IT_0618 + conjq(U_sd_31)*U_sd_33 +
       conjq(U_sd_41)*U_sd_43 + conjq(U_sd_51)*U_sd_53));
    const complex_t IT_0620 = (-0.666666666666667)*IT_0619;
    const complex_t IT_0621 = IT_0605*IT_0620;
    const complex_t IT_0622 = IT_0060*IT_0090*IT_0621;
    const complex_t IT_0623 = 0.101321183642338*IT_0622;
    const complex_t IT_0624 = IT_0011*IT_0623;
    const complex_t IT_0625 = IT_0033*IT_0623;
    const complex_t IT_0626 = IT_0068*IT_0128*IT_0621;
    const complex_t IT_0627 = 0.101321183642338*IT_0626;
    const complex_t IT_0628 = IT_0011*IT_0627;
    const complex_t IT_0629 = IT_0033*IT_0627;
    const complex_t IT_0630 = U_sd_22*conjq(U_sd_23);
    const complex_t IT_0631 = U_sd_12*conjq(U_sd_13);
    const complex_t IT_0632 = U_sd_02*conjq(U_sd_03);
    const complex_t IT_0633 = IT_0630 + IT_0631 + IT_0632;
    const complex_t IT_0634 = (complex_t{0, 1})*e_em*IT_0010*((-1.5)*IT_0018
      *IT_0633 + IT_0006*IT_0007*((-0.5)*IT_0633 + U_sd_32*conjq(U_sd_33) +
       U_sd_42*conjq(U_sd_43) + U_sd_52*conjq(U_sd_53)));
    const complex_t IT_0635 = (-0.666666666666667)*IT_0634;
    const complex_t IT_0636 = mty::lt::C0iC(9, 0, 0, 0, IT_0021, IT_0083,
       IT_0092, mty::lt::reg_int);
    const complex_t IT_0637 = IT_0635*IT_0636;
    const complex_t IT_0638 = IT_0081*IT_0129*IT_0637;
    const complex_t IT_0639 = 0.101321183642338*IT_0638;
    const complex_t IT_0640 = IT_0011*IT_0639;
    const complex_t IT_0641 = IT_0033*IT_0639;
    const complex_t IT_0642 = IT_0091*IT_0154*IT_0637;
    const complex_t IT_0643 = 0.101321183642338*IT_0642;
    const complex_t IT_0644 = IT_0011*IT_0643;
    const complex_t IT_0645 = IT_0033*IT_0643;
    const complex_t IT_0646 = conjq(U_sd_22)*U_sd_23;
    const complex_t IT_0647 = conjq(U_sd_12)*U_sd_13;
    const complex_t IT_0648 = conjq(U_sd_02)*U_sd_03;
    const complex_t IT_0649 = IT_0646 + IT_0647 + IT_0648;
    const complex_t IT_0650 = (complex_t{0, 1})*e_em*IT_0010*((-1.5)*IT_0018
      *IT_0649 + IT_0006*IT_0007*((-0.5)*IT_0649 + conjq(U_sd_32)*U_sd_33 +
       conjq(U_sd_42)*U_sd_43 + conjq(U_sd_52)*U_sd_53));
    const complex_t IT_0651 = (-0.666666666666667)*IT_0650;
    const complex_t IT_0652 = IT_0636*IT_0651;
    const complex_t IT_0653 = IT_0082*IT_0090*IT_0652;
    const complex_t IT_0654 = 0.101321183642338*IT_0653;
    const complex_t IT_0655 = IT_0011*IT_0654;
    const complex_t IT_0656 = IT_0033*IT_0654;
    const complex_t IT_0657 = IT_0128*IT_0141*IT_0652;
    const complex_t IT_0658 = 0.101321183642338*IT_0657;
    const complex_t IT_0659 = IT_0011*IT_0658;
    const complex_t IT_0660 = IT_0033*IT_0658;
    const complex_t IT_0661 = U_sd_24*conjq(U_sd_24);
    const complex_t IT_0662 = U_sd_14*conjq(U_sd_14);
    const complex_t IT_0663 = U_sd_04*conjq(U_sd_04);
    const complex_t IT_0664 = IT_0661 + IT_0662 + IT_0663;
    const complex_t IT_0665 = (complex_t{0, 1})*e_em*IT_0010*((-1.5)*IT_0018
      *IT_0664 + IT_0006*IT_0007*((-0.5)*IT_0664 + U_sd_34*conjq(U_sd_34) +
       U_sd_44*conjq(U_sd_44) + U_sd_54*conjq(U_sd_54)));
    const complex_t IT_0666 = (-0.666666666666667)*IT_0665;
    const complex_t IT_0667 = mty::lt::C0iC(9, 0, 0, 0, IT_0021, IT_0101,
       IT_0101, mty::lt::reg_int);
    const complex_t IT_0668 = IT_0666*IT_0667;
    const complex_t IT_0669 = IT_0099*IT_0111*IT_0668;
    const complex_t IT_0670 = 0.101321183642338*IT_0669;
    const complex_t IT_0671 = IT_0011*IT_0670;
    const complex_t IT_0672 = IT_0033*IT_0670;
    const complex_t IT_0673 = IT_0100*IT_0110*IT_0668;
    const complex_t IT_0674 = 0.101321183642338*IT_0673;
    const complex_t IT_0675 = IT_0011*IT_0674;
    const complex_t IT_0676 = IT_0033*IT_0674;
    const complex_t IT_0677 = conjq(U_sd_20)*U_sd_24;
    const complex_t IT_0678 = conjq(U_sd_10)*U_sd_14;
    const complex_t IT_0679 = conjq(U_sd_00)*U_sd_04;
    const complex_t IT_0680 = IT_0677 + IT_0678 + IT_0679;
    const complex_t IT_0681 = (complex_t{0, 1})*e_em*IT_0010*((-1.5)*IT_0018
      *IT_0680 + IT_0006*IT_0007*((-0.5)*IT_0680 + conjq(U_sd_30)*U_sd_34 +
       conjq(U_sd_40)*U_sd_44 + conjq(U_sd_50)*U_sd_54));
    const complex_t IT_0682 = (-0.666666666666667)*IT_0681;
    const complex_t IT_0683 = IT_0337*IT_0682;
    const complex_t IT_0684 = IT_0036*IT_0110*IT_0683;
    const complex_t IT_0685 = 0.101321183642338*IT_0684;
    const complex_t IT_0686 = IT_0011*IT_0685;
    const complex_t IT_0687 = IT_0033*IT_0685;
    const complex_t IT_0688 = IT_0013*IT_0099*IT_0683;
    const complex_t IT_0689 = 0.101321183642338*IT_0688;
    const complex_t IT_0690 = IT_0011*IT_0689;
    const complex_t IT_0691 = IT_0033*IT_0689;
    const complex_t IT_0692 = conjq(U_sd_21)*U_sd_24;
    const complex_t IT_0693 = conjq(U_sd_11)*U_sd_14;
    const complex_t IT_0694 = conjq(U_sd_01)*U_sd_04;
    const complex_t IT_0695 = IT_0692 + IT_0693 + IT_0694;
    const complex_t IT_0696 = (complex_t{0, 1})*e_em*IT_0010*((-1.5)*IT_0018
      *IT_0695 + IT_0006*IT_0007*((-0.5)*IT_0695 + conjq(U_sd_31)*U_sd_34 +
       conjq(U_sd_41)*U_sd_44 + conjq(U_sd_51)*U_sd_54));
    const complex_t IT_0697 = (-0.666666666666667)*IT_0696;
    const complex_t IT_0698 = IT_0284*IT_0697;
    const complex_t IT_0699 = IT_0060*IT_0110*IT_0698;
    const complex_t IT_0700 = 0.101321183642338*IT_0699;
    const complex_t IT_0701 = IT_0011*IT_0700;
    const complex_t IT_0702 = IT_0033*IT_0700;
    const complex_t IT_0703 = IT_0068*IT_0099*IT_0698;
    const complex_t IT_0704 = 0.101321183642338*IT_0703;
    const complex_t IT_0705 = IT_0011*IT_0704;
    const complex_t IT_0706 = IT_0033*IT_0704;
    const complex_t IT_0707 = conjq(U_sd_22)*U_sd_24;
    const complex_t IT_0708 = conjq(U_sd_12)*U_sd_14;
    const complex_t IT_0709 = conjq(U_sd_02)*U_sd_04;
    const complex_t IT_0710 = IT_0707 + IT_0708 + IT_0709;
    const complex_t IT_0711 = (complex_t{0, 1})*e_em*IT_0010*((-1.5)*IT_0018
      *IT_0710 + IT_0006*IT_0007*((-0.5)*IT_0710 + conjq(U_sd_32)*U_sd_34 +
       conjq(U_sd_42)*U_sd_44 + conjq(U_sd_52)*U_sd_54));
    const complex_t IT_0712 = (-0.666666666666667)*IT_0711;
    const complex_t IT_0713 = IT_0316*IT_0712;
    const complex_t IT_0714 = IT_0082*IT_0110*IT_0713;
    const complex_t IT_0715 = 0.101321183642338*IT_0714;
    const complex_t IT_0716 = IT_0011*IT_0715;
    const complex_t IT_0717 = IT_0033*IT_0715;
    const complex_t IT_0718 = IT_0099*IT_0141*IT_0713;
    const complex_t IT_0719 = 0.101321183642338*IT_0718;
    const complex_t IT_0720 = IT_0011*IT_0719;
    const complex_t IT_0721 = IT_0033*IT_0719;
    const complex_t IT_0722 = conjq(U_sd_23)*U_sd_24;
    const complex_t IT_0723 = conjq(U_sd_13)*U_sd_14;
    const complex_t IT_0724 = conjq(U_sd_03)*U_sd_04;
    const complex_t IT_0725 = IT_0722 + IT_0723 + IT_0724;
    const complex_t IT_0726 = (complex_t{0, 1})*e_em*IT_0010*((-1.5)*IT_0018
      *IT_0725 + IT_0006*IT_0007*((-0.5)*IT_0725 + conjq(U_sd_33)*U_sd_34 +
       conjq(U_sd_43)*U_sd_44 + conjq(U_sd_53)*U_sd_54));
    const complex_t IT_0727 = (-0.666666666666667)*IT_0726;
    const complex_t IT_0728 = mty::lt::C0iC(9, 0, 0, 0, IT_0021, IT_0092,
       IT_0101, mty::lt::reg_int);
    const complex_t IT_0729 = IT_0727*IT_0728;
    const complex_t IT_0730 = IT_0091*IT_0110*IT_0729;
    const complex_t IT_0731 = 0.101321183642338*IT_0730;
    const complex_t IT_0732 = IT_0011*IT_0731;
    const complex_t IT_0733 = IT_0033*IT_0731;
    const complex_t IT_0734 = IT_0099*IT_0129*IT_0729;
    const complex_t IT_0735 = 0.101321183642338*IT_0734;
    const complex_t IT_0736 = IT_0011*IT_0735;
    const complex_t IT_0737 = IT_0033*IT_0735;
    const complex_t IT_0738 = U_sd_23*conjq(U_sd_24);
    const complex_t IT_0739 = U_sd_13*conjq(U_sd_14);
    const complex_t IT_0740 = U_sd_03*conjq(U_sd_04);
    const complex_t IT_0741 = IT_0738 + IT_0739 + IT_0740;
    const complex_t IT_0742 = (complex_t{0, 1})*e_em*IT_0010*((-1.5)*IT_0018
      *IT_0741 + IT_0006*IT_0007*((-0.5)*IT_0741 + U_sd_33*conjq(U_sd_34) +
       U_sd_43*conjq(U_sd_44) + U_sd_53*conjq(U_sd_54)));
    const complex_t IT_0743 = (-0.666666666666667)*IT_0742;
    const complex_t IT_0744 = IT_0728*IT_0743;
    const complex_t IT_0745 = IT_0111*IT_0128*IT_0744;
    const complex_t IT_0746 = 0.101321183642338*IT_0745;
    const complex_t IT_0747 = IT_0011*IT_0746;
    const complex_t IT_0748 = IT_0033*IT_0746;
    const complex_t IT_0749 = IT_0090*IT_0100*IT_0744;
    const complex_t IT_0750 = 0.101321183642338*IT_0749;
    const complex_t IT_0751 = IT_0011*IT_0750;
    const complex_t IT_0752 = IT_0033*IT_0750;
    const complex_t IT_0753 = U_sd_25*conjq(U_sd_25);
    const complex_t IT_0754 = U_sd_15*conjq(U_sd_15);
    const complex_t IT_0755 = U_sd_05*conjq(U_sd_05);
    const complex_t IT_0756 = IT_0753 + IT_0754 + IT_0755;
    const complex_t IT_0757 = (complex_t{0, 1})*e_em*IT_0010*((-1.5)*IT_0018
      *IT_0756 + IT_0006*IT_0007*((-0.5)*IT_0756 + U_sd_35*conjq(U_sd_35) +
       U_sd_45*conjq(U_sd_45) + U_sd_55*conjq(U_sd_55)));
    const complex_t IT_0758 = (-0.666666666666667)*IT_0757;
    const complex_t IT_0759 = mty::lt::C0iC(9, 0, 0, 0, IT_0021, IT_0117,
       IT_0117, mty::lt::reg_int);
    const complex_t IT_0760 = IT_0758*IT_0759;
    const complex_t IT_0761 = IT_0116*IT_0170*IT_0760;
    const complex_t IT_0762 = 0.101321183642338*IT_0761;
    const complex_t IT_0763 = IT_0011*IT_0762;
    const complex_t IT_0764 = IT_0033*IT_0762;
    const complex_t IT_0765 = IT_0115*IT_0171*IT_0760;
    const complex_t IT_0766 = 0.101321183642338*IT_0765;
    const complex_t IT_0767 = IT_0011*IT_0766;
    const complex_t IT_0768 = IT_0033*IT_0766;
    const complex_t IT_0769 = U_sd_20*conjq(U_sd_25);
    const complex_t IT_0770 = U_sd_10*conjq(U_sd_15);
    const complex_t IT_0771 = U_sd_00*conjq(U_sd_05);
    const complex_t IT_0772 = IT_0769 + IT_0770 + IT_0771;
    const complex_t IT_0773 = (complex_t{0, 1})*e_em*IT_0010*((-1.5)*IT_0018
      *IT_0772 + IT_0006*IT_0007*((-0.5)*IT_0772 + U_sd_30*conjq(U_sd_35) +
       U_sd_40*conjq(U_sd_45) + U_sd_50*conjq(U_sd_55)));
    const complex_t IT_0774 = (-0.666666666666667)*IT_0773;
    const complex_t IT_0775 = IT_0353*IT_0774;
    const complex_t IT_0776 = IT_0012*IT_0116*IT_0775;
    const complex_t IT_0777 = 0.101321183642338*IT_0776;
    const complex_t IT_0778 = IT_0011*IT_0777;
    const complex_t IT_0779 = IT_0033*IT_0777;
    const complex_t IT_0780 = IT_0035*IT_0171*IT_0775;
    const complex_t IT_0781 = 0.101321183642338*IT_0780;
    const complex_t IT_0782 = IT_0011*IT_0781;
    const complex_t IT_0783 = IT_0033*IT_0781;
    const complex_t IT_0784 = U_sd_21*conjq(U_sd_25);
    const complex_t IT_0785 = U_sd_11*conjq(U_sd_15);
    const complex_t IT_0786 = U_sd_01*conjq(U_sd_05);
    const complex_t IT_0787 = IT_0784 + IT_0785 + IT_0786;
    const complex_t IT_0788 = (complex_t{0, 1})*e_em*IT_0010*((-1.5)*IT_0018
      *IT_0787 + IT_0006*IT_0007*((-0.5)*IT_0787 + U_sd_31*conjq(U_sd_35) +
       U_sd_41*conjq(U_sd_45) + U_sd_51*conjq(U_sd_55)));
    const complex_t IT_0789 = (-0.666666666666667)*IT_0788;
    const complex_t IT_0790 = IT_0268*IT_0789;
    const complex_t IT_0791 = IT_0059*IT_0116*IT_0790;
    const complex_t IT_0792 = 0.101321183642338*IT_0791;
    const complex_t IT_0793 = IT_0011*IT_0792;
    const complex_t IT_0794 = IT_0033*IT_0792;
    const complex_t IT_0795 = IT_0075*IT_0171*IT_0790;
    const complex_t IT_0796 = 0.101321183642338*IT_0795;
    const complex_t IT_0797 = IT_0011*IT_0796;
    const complex_t IT_0798 = IT_0033*IT_0796;
    const complex_t IT_0799 = conjq(U_sd_22)*U_sd_25;
    const complex_t IT_0800 = conjq(U_sd_12)*U_sd_15;
    const complex_t IT_0801 = conjq(U_sd_02)*U_sd_05;
    const complex_t IT_0802 = IT_0799 + IT_0800 + IT_0801;
    const complex_t IT_0803 = (complex_t{0, 1})*e_em*IT_0010*((-1.5)*IT_0018
      *IT_0802 + IT_0006*IT_0007*((-0.5)*IT_0802 + conjq(U_sd_32)*U_sd_35 +
       conjq(U_sd_42)*U_sd_45 + conjq(U_sd_52)*U_sd_55));
    const complex_t IT_0804 = (-0.666666666666667)*IT_0803;
    const complex_t IT_0805 = mty::lt::C0iC(9, 0, 0, 0, IT_0021, IT_0083,
       IT_0117, mty::lt::reg_int);
    const complex_t IT_0806 = IT_0804*IT_0805;
    const complex_t IT_0807 = IT_0082*IT_0115*IT_0806;
    const complex_t IT_0808 = 0.101321183642338*IT_0807;
    const complex_t IT_0809 = IT_0011*IT_0808;
    const complex_t IT_0810 = IT_0033*IT_0808;
    const complex_t IT_0811 = IT_0141*IT_0170*IT_0806;
    const complex_t IT_0812 = 0.101321183642338*IT_0811;
    const complex_t IT_0813 = IT_0011*IT_0812;
    const complex_t IT_0814 = IT_0033*IT_0812;
    const complex_t IT_0815 = U_sd_22*conjq(U_sd_25);
    const complex_t IT_0816 = U_sd_12*conjq(U_sd_15);
    const complex_t IT_0817 = U_sd_02*conjq(U_sd_05);
    const complex_t IT_0818 = IT_0815 + IT_0816 + IT_0817;
    const complex_t IT_0819 = (complex_t{0, 1})*e_em*IT_0010*((-1.5)*IT_0018
      *IT_0818 + IT_0006*IT_0007*((-0.5)*IT_0818 + U_sd_32*conjq(U_sd_35) +
       U_sd_42*conjq(U_sd_45) + U_sd_52*conjq(U_sd_55)));
    const complex_t IT_0820 = (-0.666666666666667)*IT_0819;
    const complex_t IT_0821 = IT_0805*IT_0820;
    const complex_t IT_0822 = IT_0081*IT_0116*IT_0821;
    const complex_t IT_0823 = 0.101321183642338*IT_0822;
    const complex_t IT_0824 = IT_0011*IT_0823;
    const complex_t IT_0825 = IT_0033*IT_0823;
    const complex_t IT_0826 = IT_0154*IT_0171*IT_0821;
    const complex_t IT_0827 = 0.101321183642338*IT_0826;
    const complex_t IT_0828 = IT_0011*IT_0827;
    const complex_t IT_0829 = IT_0033*IT_0827;
    const complex_t IT_0830 = U_sd_23*conjq(U_sd_25);
    const complex_t IT_0831 = U_sd_13*conjq(U_sd_15);
    const complex_t IT_0832 = U_sd_03*conjq(U_sd_05);
    const complex_t IT_0833 = IT_0830 + IT_0831 + IT_0832;
    const complex_t IT_0834 = (complex_t{0, 1})*e_em*IT_0010*((-1.5)*IT_0018
      *IT_0833 + IT_0006*IT_0007*((-0.5)*IT_0833 + U_sd_33*conjq(U_sd_35) +
       U_sd_43*conjq(U_sd_45) + U_sd_53*conjq(U_sd_55)));
    const complex_t IT_0835 = (-0.666666666666667)*IT_0834;
    const complex_t IT_0836 = mty::lt::C0iC(9, 0, 0, 0, IT_0021, IT_0117,
       IT_0092, mty::lt::reg_int);
    const complex_t IT_0837 = IT_0835*IT_0836;
    const complex_t IT_0838 = IT_0116*IT_0128*IT_0837;
    const complex_t IT_0839 = 0.101321183642338*IT_0838;
    const complex_t IT_0840 = IT_0011*IT_0839;
    const complex_t IT_0841 = IT_0033*IT_0839;
    const complex_t IT_0842 = IT_0090*IT_0171*IT_0837;
    const complex_t IT_0843 = 0.101321183642338*IT_0842;
    const complex_t IT_0844 = IT_0011*IT_0843;
    const complex_t IT_0845 = IT_0033*IT_0843;
    const complex_t IT_0846 = conjq(U_sd_23)*U_sd_25;
    const complex_t IT_0847 = conjq(U_sd_13)*U_sd_15;
    const complex_t IT_0848 = conjq(U_sd_03)*U_sd_05;
    const complex_t IT_0849 = IT_0846 + IT_0847 + IT_0848;
    const complex_t IT_0850 = (complex_t{0, 1})*e_em*IT_0010*((-1.5)*IT_0018
      *IT_0849 + IT_0006*IT_0007*((-0.5)*IT_0849 + conjq(U_sd_33)*U_sd_35 +
       conjq(U_sd_43)*U_sd_45 + conjq(U_sd_53)*U_sd_55));
    const complex_t IT_0851 = (-0.666666666666667)*IT_0850;
    const complex_t IT_0852 = IT_0836*IT_0851;
    const complex_t IT_0853 = IT_0091*IT_0115*IT_0852;
    const complex_t IT_0854 = 0.101321183642338*IT_0853;
    const complex_t IT_0855 = IT_0011*IT_0854;
    const complex_t IT_0856 = IT_0033*IT_0854;
    const complex_t IT_0857 = IT_0129*IT_0170*IT_0852;
    const complex_t IT_0858 = 0.101321183642338*IT_0857;
    const complex_t IT_0859 = IT_0011*IT_0858;
    const complex_t IT_0860 = IT_0033*IT_0858;
    const complex_t IT_0861 = U_sd_24*conjq(U_sd_25);
    const complex_t IT_0862 = U_sd_14*conjq(U_sd_15);
    const complex_t IT_0863 = U_sd_04*conjq(U_sd_05);
    const complex_t IT_0864 = IT_0861 + IT_0862 + IT_0863;
    const complex_t IT_0865 = (complex_t{0, 1})*e_em*IT_0010*((-1.5)*IT_0018
      *IT_0864 + IT_0006*IT_0007*((-0.5)*IT_0864 + U_sd_34*conjq(U_sd_35) +
       U_sd_44*conjq(U_sd_45) + U_sd_54*conjq(U_sd_55)));
    const complex_t IT_0866 = (-0.666666666666667)*IT_0865;
    const complex_t IT_0867 = mty::lt::C0iC(9, 0, 0, 0, IT_0021, IT_0117,
       IT_0101, mty::lt::reg_int);
    const complex_t IT_0868 = IT_0866*IT_0867;
    const complex_t IT_0869 = IT_0099*IT_0116*IT_0868;
    const complex_t IT_0870 = 0.101321183642338*IT_0869;
    const complex_t IT_0871 = IT_0011*IT_0870;
    const complex_t IT_0872 = IT_0033*IT_0870;
    const complex_t IT_0873 = IT_0110*IT_0171*IT_0868;
    const complex_t IT_0874 = 0.101321183642338*IT_0873;
    const complex_t IT_0875 = IT_0011*IT_0874;
    const complex_t IT_0876 = IT_0033*IT_0874;
    const complex_t IT_0877 = conjq(U_sd_24)*U_sd_25;
    const complex_t IT_0878 = conjq(U_sd_14)*U_sd_15;
    const complex_t IT_0879 = conjq(U_sd_04)*U_sd_05;
    const complex_t IT_0880 = IT_0877 + IT_0878 + IT_0879;
    const complex_t IT_0881 = (complex_t{0, 1})*e_em*IT_0010*((-1.5)*IT_0018
      *IT_0880 + IT_0006*IT_0007*((-0.5)*IT_0880 + conjq(U_sd_34)*U_sd_35 +
       conjq(U_sd_44)*U_sd_45 + conjq(U_sd_54)*U_sd_55));
    const complex_t IT_0882 = (-0.666666666666667)*IT_0881;
    const complex_t IT_0883 = IT_0867*IT_0882;
    const complex_t IT_0884 = IT_0100*IT_0115*IT_0883;
    const complex_t IT_0885 = 0.101321183642338*IT_0884;
    const complex_t IT_0886 = IT_0011*IT_0885;
    const complex_t IT_0887 = IT_0033*IT_0885;
    const complex_t IT_0888 = IT_0111*IT_0170*IT_0883;
    const complex_t IT_0889 = 0.101321183642338*IT_0888;
    const complex_t IT_0890 = IT_0011*IT_0889;
    const complex_t IT_0891 = IT_0033*IT_0889;
    const complex_t IT_0892 = IT_0012*IT_0036*IT_0048*IT_0124;
    const complex_t IT_0893 = IT_0108*IT_0109*IT_0892;
    const complex_t IT_0894 = IT_0011*IT_0893;
    const complex_t IT_0895 = IT_0033*IT_0893;
    const complex_t IT_0896 = IT_0012*IT_0013*IT_0124*IT_0207;
    const complex_t IT_0897 = IT_0108*IT_0109*IT_0896;
    const complex_t IT_0898 = IT_0011*IT_0897;
    const complex_t IT_0899 = IT_0033*IT_0897;
    const complex_t IT_0900 = IT_0012*IT_0036*IT_0124*IT_0195;
    const complex_t IT_0901 = 0.101321183642338*IT_0043*IT_0900;
    const complex_t IT_0902 = IT_0011*IT_0901;
    const complex_t IT_0903 = IT_0033*IT_0901;
    const complex_t IT_0904 = IT_0012*IT_0013*IT_0054*IT_0124;
    const complex_t IT_0905 = IT_0043*IT_0044*IT_0904;
    const complex_t IT_0906 = IT_0011*IT_0905;
    const complex_t IT_0907 = IT_0033*IT_0905;
    const complex_t IT_0908 = IT_0035*IT_0036*IT_0124*IT_0247;
    const complex_t IT_0909 = 0.101321183642338*IT_0108*IT_0908;
    const complex_t IT_0910 = IT_0011*IT_0909;
    const complex_t IT_0911 = IT_0033*IT_0909;
    const complex_t IT_0912 = IT_0035*IT_0036*IT_0124*IT_0389;
    const complex_t IT_0913 = 0.101321183642338*IT_0043*IT_0912;
    const complex_t IT_0914 = IT_0011*IT_0913;
    const complex_t IT_0915 = IT_0033*IT_0913;
    const complex_t IT_0916 = IT_0059*IT_0060*IT_0063*IT_0124;
    const complex_t IT_0917 = IT_0108*IT_0109*IT_0916;
    const complex_t IT_0918 = IT_0011*IT_0917;
    const complex_t IT_0919 = IT_0033*IT_0917;
    const complex_t IT_0920 = IT_0059*IT_0068*IT_0124*IT_0440;
    const complex_t IT_0921 = IT_0108*IT_0109*IT_0920;
    const complex_t IT_0922 = IT_0011*IT_0921;
    const complex_t IT_0923 = IT_0033*IT_0921;
    const complex_t IT_0924 = IT_0011*IT_0218;
    const complex_t IT_0925 = IT_0059*IT_0068*IT_0076*IT_0124;
    const complex_t IT_0926 = IT_0043*IT_0044*IT_0925;
    const complex_t IT_0927 = IT_0011*IT_0926;
    const complex_t IT_0928 = IT_0033*IT_0926;
    const complex_t IT_0929 = IT_0060*IT_0075*IT_0124*IT_0435;
    const complex_t IT_0930 = 0.101321183642338*IT_0108*IT_0929;
    const complex_t IT_0931 = IT_0011*IT_0930;
    const complex_t IT_0932 = IT_0033*IT_0930;
    const complex_t IT_0933 = IT_0060*IT_0070*IT_0075*IT_0124;
    const complex_t IT_0934 = 0.101321183642338*IT_0043*IT_0933;
    const complex_t IT_0935 = IT_0011*IT_0934;
    const complex_t IT_0936 = IT_0033*IT_0934;
    const complex_t IT_0937 = IT_0081*IT_0082*IT_0085*IT_0124;
    const complex_t IT_0938 = IT_0108*IT_0109*IT_0937;
    const complex_t IT_0939 = IT_0011*IT_0938;
    const complex_t IT_0940 = IT_0033*IT_0938;
    const complex_t IT_0941 = IT_0081*IT_0124*IT_0141*IT_0456;
    const complex_t IT_0942 = IT_0108*IT_0109*IT_0941;
    const complex_t IT_0943 = IT_0011*IT_0942;
    const complex_t IT_0944 = IT_0033*IT_0942;
    const complex_t IT_0945 = IT_0081*IT_0082*IT_0124*IT_0243;
    const complex_t IT_0946 = 0.101321183642338*IT_0043*IT_0945;
    const complex_t IT_0947 = IT_0011*IT_0946;
    const complex_t IT_0948 = IT_0033*IT_0946;
    const complex_t IT_0949 = IT_0081*IT_0124*IT_0141*IT_0400;
    const complex_t IT_0950 = IT_0043*IT_0044*IT_0949;
    const complex_t IT_0951 = IT_0011*IT_0950;
    const complex_t IT_0952 = IT_0033*IT_0950;
    const complex_t IT_0953 = IT_0082*IT_0124*IT_0154*IT_0451;
    const complex_t IT_0954 = 0.101321183642338*IT_0108*IT_0953;
    const complex_t IT_0955 = IT_0011*IT_0954;
    const complex_t IT_0956 = IT_0033*IT_0954;
    const complex_t IT_0957 = IT_0082*IT_0124*IT_0154*IT_0395;
    const complex_t IT_0958 = 0.101321183642338*IT_0043*IT_0957;
    const complex_t IT_0959 = IT_0011*IT_0958;
    const complex_t IT_0960 = IT_0033*IT_0958;
    const complex_t IT_0961 = IT_0091*IT_0124*IT_0128*IT_0178;
    const complex_t IT_0962 = IT_0108*IT_0109*IT_0961;
    const complex_t IT_0963 = IT_0011*IT_0962;
    const complex_t IT_0964 = IT_0033*IT_0962;
    const complex_t IT_0965 = IT_0124*IT_0128*IT_0129*IT_0472;
    const complex_t IT_0966 = IT_0108*IT_0109*IT_0965;
    const complex_t IT_0967 = IT_0011*IT_0966;
    const complex_t IT_0968 = IT_0033*IT_0966;
    const complex_t IT_0969 = IT_0091*IT_0124*IT_0128*IT_0420;
    const complex_t IT_0970 = 0.101321183642338*IT_0043*IT_0969;
    const complex_t IT_0971 = IT_0011*IT_0970;
    const complex_t IT_0972 = IT_0033*IT_0970;
    const complex_t IT_0973 = IT_0094*IT_0124*IT_0128*IT_0129;
    const complex_t IT_0974 = IT_0043*IT_0044*IT_0973;
    const complex_t IT_0975 = IT_0011*IT_0974;
    const complex_t IT_0976 = IT_0033*IT_0974;
    const complex_t IT_0977 = IT_0090*IT_0091*IT_0124*IT_0467;
    const complex_t IT_0978 = 0.101321183642338*IT_0108*IT_0977;
    const complex_t IT_0979 = IT_0011*IT_0978;
    const complex_t IT_0980 = IT_0033*IT_0978;
    const complex_t IT_0981 = IT_0090*IT_0091*IT_0124*IT_0406;
    const complex_t IT_0982 = 0.101321183642338*IT_0043*IT_0981;
    const complex_t IT_0983 = IT_0011*IT_0982;
    const complex_t IT_0984 = IT_0033*IT_0982;
    const complex_t IT_0985 = IT_0099*IT_0100*IT_0103*IT_0124;
    const complex_t IT_0986 = IT_0108*IT_0109*IT_0985;
    const complex_t IT_0987 = IT_0011*IT_0986;
    const complex_t IT_0988 = IT_0033*IT_0986;
    const complex_t IT_0989 = IT_0099*IT_0111*IT_0124*IT_0488;
    const complex_t IT_0990 = IT_0108*IT_0109*IT_0989;
    const complex_t IT_0991 = IT_0011*IT_0990;
    const complex_t IT_0992 = IT_0033*IT_0990;
    const complex_t IT_0993 = IT_0011*IT_0257;
    const complex_t IT_0994 = IT_0099*IT_0111*IT_0124*IT_0235;
    const complex_t IT_0995 = IT_0043*IT_0044*IT_0994;
    const complex_t IT_0996 = IT_0011*IT_0995;
    const complex_t IT_0997 = IT_0033*IT_0995;
    const complex_t IT_0998 = IT_0100*IT_0110*IT_0124*IT_0483;
    const complex_t IT_0999 = 0.101321183642338*IT_0108*IT_0998;
    const complex_t IT_1000 = IT_0011*IT_0999;
    const complex_t IT_1001 = IT_0033*IT_0999;
    const complex_t IT_1002 = IT_0100*IT_0110*IT_0124*IT_0411;
    const complex_t IT_1003 = 0.101321183642338*IT_0043*IT_1002;
    const complex_t IT_1004 = IT_0011*IT_1003;
    const complex_t IT_1005 = IT_0033*IT_1003;
    const complex_t IT_1006 = IT_0119*IT_0124*IT_0170*IT_0171;
    const complex_t IT_1007 = IT_0108*IT_0109*IT_1006;
    const complex_t IT_1008 = IT_0011*IT_1007;
    const complex_t IT_1009 = IT_0033*IT_1007;
    const complex_t IT_1010 = IT_0116*IT_0124*IT_0170*IT_0504;
    const complex_t IT_1011 = IT_0108*IT_0109*IT_1010;
    const complex_t IT_1012 = IT_0011*IT_1011;
    const complex_t IT_1013 = IT_0033*IT_1011;
    const complex_t IT_1014 = IT_0033*IT_0175;
    const complex_t IT_1015 = IT_0033*IT_0221;
    const complex_t IT_1016 = IT_0115*IT_0124*IT_0171*IT_0499;
    const complex_t IT_1017 = 0.101321183642338*IT_0108*IT_1016;
    const complex_t IT_1018 = IT_0011*IT_1017;
    const complex_t IT_1019 = IT_0033*IT_1017;
    const complex_t IT_1020 = IT_0011*IT_0260;
    const complex_t IT_1021 = IT_0013*IT_0035*IT_0124*IT_0201;
    const complex_t IT_1022 = 0.101321183642338*IT_0108*IT_1021;
    const complex_t IT_1023 = IT_0011*IT_1022;
    const complex_t IT_1024 = IT_0033*IT_1022;
    const complex_t IT_1025 = IT_0068*IT_0075*IT_0124*IT_0212;
    const complex_t IT_1026 = 0.101321183642338*IT_0108*IT_1025;
    const complex_t IT_1027 = IT_0011*IT_1026;
    const complex_t IT_1028 = IT_0033*IT_1026;
    const complex_t IT_1029 = IT_0124*IT_0141*IT_0154*IT_0445;
    const complex_t IT_1030 = 0.101321183642338*IT_0108*IT_1029;
    const complex_t IT_1031 = IT_0011*IT_1030;
    const complex_t IT_1032 = IT_0033*IT_1030;
    const complex_t IT_1033 = IT_0090*IT_0124*IT_0129*IT_0461;
    const complex_t IT_1034 = 0.101321183642338*IT_0108*IT_1033;
    const complex_t IT_1035 = IT_0011*IT_1034;
    const complex_t IT_1036 = IT_0033*IT_1034;
    const complex_t IT_1037 = IT_0110*IT_0111*IT_0124*IT_0477;
    const complex_t IT_1038 = 0.101321183642338*IT_0108*IT_1037;
    const complex_t IT_1039 = IT_0011*IT_1038;
    const complex_t IT_1040 = IT_0033*IT_1038;
    const complex_t IT_1041 = IT_0115*IT_0116*IT_0124*IT_0493;
    const complex_t IT_1042 = 0.101321183642338*IT_0108*IT_1041;
    const complex_t IT_1043 = IT_0011*IT_1042;
    const complex_t IT_1044 = IT_0033*IT_1042;
    const complex_t IT_1045 = IT_0013*IT_0035*IT_0048*IT_0124;
    const complex_t IT_1046 = IT_0043*IT_0044*IT_1045;
    const complex_t IT_1047 = IT_0011*IT_1046;
    const complex_t IT_1048 = IT_0033*IT_1046;
    const complex_t IT_1049 = IT_0063*IT_0068*IT_0075*IT_0124;
    const complex_t IT_1050 = IT_0043*IT_0044*IT_1049;
    const complex_t IT_1051 = IT_0011*IT_1050;
    const complex_t IT_1052 = IT_0033*IT_1050;
    const complex_t IT_1053 = IT_0085*IT_0124*IT_0141*IT_0154;
    const complex_t IT_1054 = IT_0043*IT_0044*IT_1053;
    const complex_t IT_1055 = IT_0011*IT_1054;
    const complex_t IT_1056 = IT_0033*IT_1054;
    const complex_t IT_1057 = IT_0090*IT_0124*IT_0129*IT_0178;
    const complex_t IT_1058 = IT_0043*IT_0044*IT_1057;
    const complex_t IT_1059 = IT_0011*IT_1058;
    const complex_t IT_1060 = IT_0033*IT_1058;
    const complex_t IT_1061 = IT_0103*IT_0110*IT_0111*IT_0124;
    const complex_t IT_1062 = IT_0043*IT_0044*IT_1061;
    const complex_t IT_1063 = IT_0011*IT_1062;
    const complex_t IT_1064 = IT_0033*IT_1062;
    const complex_t IT_1065 = IT_0011*IT_0126;
    const complex_t IT_1066 = IT_0027 + IT_0034 + IT_0039 + IT_0040 + -IT_0051
       + -IT_0052 + -IT_0057 + -IT_0058 + -IT_0066 + -IT_0067 + -IT_0073 + 
      -IT_0074 + -IT_0079 + -IT_0080 + -IT_0088 + -IT_0089 + -IT_0097 + -IT_0098
       + -IT_0106 + -IT_0107 + IT_0114 + IT_0122 + IT_0123 + -IT_0127 + IT_0140 
      + IT_0152 + IT_0153 + IT_0157 + IT_0158 + IT_0169 + -IT_0176 + -IT_0181 + 
      -IT_0184 + -IT_0189 + -IT_0193 + -IT_0194 + -IT_0198 + -IT_0199 + IT_0204 
      + IT_0205 + IT_0210 + IT_0211 + IT_0215 + -IT_0219 + -IT_0222 + IT_0233 + 
      -IT_0238 + -IT_0241 + -IT_0242 + -IT_0246 + IT_0250 + IT_0251 + IT_0254 + 
      -IT_0258 + -IT_0261 + IT_0272 + IT_0273 + IT_0276 + IT_0277 + IT_0288 +
       IT_0289 + IT_0292 + IT_0293 + IT_0304 + IT_0305 + IT_0308 + IT_0309 +
       IT_0320 + IT_0321 + IT_0324 + IT_0325 + IT_0326 + IT_0329 + IT_0330 +
       IT_0341 + IT_0342 + IT_0345 + IT_0346 + IT_0357 + IT_0358 + IT_0361 +
       IT_0362 + IT_0373 + IT_0374 + IT_0377 + IT_0378 + IT_0381 + IT_0382 +
       IT_0383 + IT_0384 + IT_0387 + IT_0388 + -IT_0392 + -IT_0393 + -IT_0398 + 
      -IT_0399 + -IT_0403 + -IT_0404 + -IT_0405 + -IT_0409 + -IT_0410 + -IT_0414
       + -IT_0415 + -IT_0416 + -IT_0417 + -IT_0418 + -IT_0419 + -IT_0423 + 
      -IT_0424 + -IT_0427 + -IT_0428 + -IT_0431 + -IT_0432 + IT_0433 + IT_0438 +
       IT_0439 + IT_0443 + IT_0444 + IT_0448 + IT_0449 + IT_0454 + IT_0455 +
       IT_0459 + IT_0460 + IT_0464 + IT_0465 + IT_0470 + IT_0471 + IT_0475 +
       IT_0476 + IT_0480 + IT_0481 + IT_0486 + IT_0487 + IT_0491 + IT_0492 +
       IT_0496 + IT_0497 + IT_0502 + IT_0503 + IT_0507 + IT_0508 + IT_0509 +
       IT_0512 + IT_0513 + IT_0516 + IT_0517 + IT_0520 + IT_0521 + IT_0522 +
       IT_0532 + IT_0533 + IT_0536 + IT_0537 + IT_0547 + IT_0548 + IT_0551 +
       IT_0552 + IT_0562 + IT_0563 + IT_0566 + IT_0567 + IT_0578 + IT_0579 +
       IT_0582 + IT_0583 + IT_0593 + IT_0594 + IT_0597 + IT_0598 + IT_0609 +
       IT_0610 + IT_0613 + IT_0614 + IT_0624 + IT_0625 + IT_0628 + IT_0629 +
       IT_0640 + IT_0641 + IT_0644 + IT_0645 + IT_0655 + IT_0656 + IT_0659 +
       IT_0660 + IT_0671 + IT_0672 + IT_0675 + IT_0676 + IT_0686 + IT_0687 +
       IT_0690 + IT_0691 + IT_0701 + IT_0702 + IT_0705 + IT_0706 + IT_0716 +
       IT_0717 + IT_0720 + IT_0721 + IT_0732 + IT_0733 + IT_0736 + IT_0737 +
       IT_0747 + IT_0748 + IT_0751 + IT_0752 + IT_0763 + IT_0764 + IT_0767 +
       IT_0768 + IT_0778 + IT_0779 + IT_0782 + IT_0783 + IT_0793 + IT_0794 +
       IT_0797 + IT_0798 + IT_0809 + IT_0810 + IT_0813 + IT_0814 + IT_0824 +
       IT_0825 + IT_0828 + IT_0829 + IT_0840 + IT_0841 + IT_0844 + IT_0845 +
       IT_0855 + IT_0856 + IT_0859 + IT_0860 + IT_0871 + IT_0872 + IT_0875 +
       IT_0876 + IT_0886 + IT_0887 + IT_0890 + IT_0891 + IT_0894 + IT_0895 +
       IT_0898 + IT_0899 + -IT_0902 + -IT_0903 + -IT_0906 + -IT_0907 + IT_0910 +
       IT_0911 + -IT_0914 + -IT_0915 + IT_0918 + IT_0919 + IT_0922 + IT_0923 + 
      -IT_0924 + -IT_0927 + -IT_0928 + IT_0931 + IT_0932 + -IT_0935 + -IT_0936 +
       IT_0939 + IT_0940 + IT_0943 + IT_0944 + -IT_0947 + -IT_0948 + -IT_0951 + 
      -IT_0952 + IT_0955 + IT_0956 + -IT_0959 + -IT_0960 + IT_0963 + IT_0964 +
       IT_0967 + IT_0968 + -IT_0971 + -IT_0972 + -IT_0975 + -IT_0976 + IT_0979 +
       IT_0980 + -IT_0983 + -IT_0984 + IT_0987 + IT_0988 + IT_0991 + IT_0992 + 
      -IT_0993 + -IT_0996 + -IT_0997 + IT_1000 + IT_1001 + -IT_1004 + -IT_1005 +
       IT_1008 + IT_1009 + IT_1012 + IT_1013 + -IT_1014 + -IT_1015 + IT_1018 +
       IT_1019 + -IT_1020 + IT_1023 + IT_1024 + IT_1027 + IT_1028 + IT_1031 +
       IT_1032 + IT_1035 + IT_1036 + IT_1039 + IT_1040 + IT_1043 + IT_1044 + 
      -IT_1047 + -IT_1048 + -IT_1051 + -IT_1052 + -IT_1055 + -IT_1056 + -IT_1059
       + -IT_1060 + -IT_1063 + -IT_1064 + -IT_1065;
    const complex_t IT_1067 = powq(M_Z, 2);
    const complex_t IT_1068 = powq(m_mu, 2);
    const complex_t IT_1069 = cpowq((-2)*s_34 + IT_1067 + (-2)*IT_1068 + 
      -reg_prop, -1);
    const complex_t IT_1070 = cpowq(IT_0007, 2);
    const complex_t IT_1071 = IT_1066*IT_1069*IT_1070;
    const complex_t IT_1072 = IT_0004*IT_1071;
    const complex_t IT_1073 = (-0.666666666666667)*IT_1072;
    const complex_t IT_1074 = IT_0027 + IT_0034 + -IT_0039 + -IT_0040 + 
      -IT_0051 + -IT_0052 + -IT_0057 + -IT_0058 + -IT_0066 + -IT_0067 + -IT_0073
       + -IT_0074 + -IT_0079 + -IT_0080 + -IT_0088 + -IT_0089 + -IT_0097 + 
      -IT_0098 + -IT_0106 + -IT_0107 + IT_0114 + IT_0122 + IT_0123 + IT_0127 +
       IT_0140 + IT_0152 + IT_0153 + -IT_0157 + -IT_0158 + IT_0169 + IT_0176 + 
      -IT_0181 + -IT_0184 + -IT_0189 + -IT_0193 + -IT_0194 + -IT_0198 + -IT_0199
       + IT_0204 + IT_0205 + IT_0210 + IT_0211 + IT_0215 + IT_0219 + IT_0222 + 
      -IT_0233 + -IT_0238 + -IT_0241 + -IT_0242 + -IT_0246 + IT_0250 + IT_0251 +
       IT_0254 + IT_0258 + IT_0261 + -IT_0272 + -IT_0273 + IT_0276 + IT_0277 +
       IT_0288 + IT_0289 + -IT_0292 + -IT_0293 + -IT_0304 + -IT_0305 + IT_0308 +
       IT_0309 + IT_0320 + IT_0321 + -IT_0324 + -IT_0325 + IT_0326 + -IT_0329 + 
      -IT_0330 + IT_0341 + IT_0342 + -IT_0345 + -IT_0346 + -IT_0357 + -IT_0358 +
       IT_0361 + IT_0362 + -IT_0373 + -IT_0374 + IT_0377 + IT_0378 + IT_0381 +
       IT_0382 + -IT_0383 + IT_0384 + -IT_0387 + -IT_0388 + -IT_0392 + -IT_0393 
      + -IT_0398 + -IT_0399 + -IT_0403 + -IT_0404 + -IT_0405 + -IT_0409 + 
      -IT_0410 + -IT_0414 + -IT_0415 + -IT_0416 + -IT_0417 + -IT_0418 + -IT_0419
       + -IT_0423 + -IT_0424 + -IT_0427 + -IT_0428 + -IT_0431 + -IT_0432 +
       IT_0433 + IT_0438 + IT_0439 + IT_0443 + IT_0444 + IT_0448 + IT_0449 +
       IT_0454 + IT_0455 + IT_0459 + IT_0460 + IT_0464 + IT_0465 + IT_0470 +
       IT_0471 + IT_0475 + IT_0476 + IT_0480 + IT_0481 + IT_0486 + IT_0487 +
       IT_0491 + IT_0492 + IT_0496 + IT_0497 + IT_0502 + IT_0503 + IT_0507 +
       IT_0508 + IT_0509 + IT_0512 + IT_0513 + IT_0516 + IT_0517 + IT_0520 +
       IT_0521 + IT_0522 + -IT_0532 + -IT_0533 + IT_0536 + IT_0537 + IT_0547 +
       IT_0548 + -IT_0551 + -IT_0552 + IT_0562 + IT_0563 + -IT_0566 + -IT_0567 +
       IT_0578 + IT_0579 + -IT_0582 + -IT_0583 + -IT_0593 + -IT_0594 + IT_0597 +
       IT_0598 + IT_0609 + IT_0610 + -IT_0613 + -IT_0614 + -IT_0624 + -IT_0625 +
       IT_0628 + IT_0629 + IT_0640 + IT_0641 + -IT_0644 + -IT_0645 + -IT_0655 + 
      -IT_0656 + IT_0659 + IT_0660 + IT_0671 + IT_0672 + -IT_0675 + -IT_0676 + 
      -IT_0686 + -IT_0687 + IT_0690 + IT_0691 + -IT_0701 + -IT_0702 + IT_0705 +
       IT_0706 + -IT_0716 + -IT_0717 + IT_0720 + IT_0721 + -IT_0732 + -IT_0733 +
       IT_0736 + IT_0737 + IT_0747 + IT_0748 + -IT_0751 + -IT_0752 + IT_0763 +
       IT_0764 + -IT_0767 + -IT_0768 + IT_0778 + IT_0779 + -IT_0782 + -IT_0783 +
       IT_0793 + IT_0794 + -IT_0797 + -IT_0798 + -IT_0809 + -IT_0810 + IT_0813 +
       IT_0814 + IT_0824 + IT_0825 + -IT_0828 + -IT_0829 + IT_0840 + IT_0841 + 
      -IT_0844 + -IT_0845 + -IT_0855 + -IT_0856 + IT_0859 + IT_0860 + IT_0871 +
       IT_0872 + -IT_0875 + -IT_0876 + -IT_0886 + -IT_0887 + IT_0890 + IT_0891 +
       -IT_0894 + -IT_0895 + -IT_0898 + -IT_0899 + IT_0902 + IT_0903 + IT_0906 +
       IT_0907 + -IT_0910 + -IT_0911 + IT_0914 + IT_0915 + -IT_0918 + -IT_0919 +
       -IT_0922 + -IT_0923 + IT_0924 + IT_0927 + IT_0928 + -IT_0931 + -IT_0932 +
       IT_0935 + IT_0936 + -IT_0939 + -IT_0940 + -IT_0943 + -IT_0944 + IT_0947 +
       IT_0948 + IT_0951 + IT_0952 + -IT_0955 + -IT_0956 + IT_0959 + IT_0960 + 
      -IT_0963 + -IT_0964 + -IT_0967 + -IT_0968 + IT_0971 + IT_0972 + IT_0975 +
       IT_0976 + -IT_0979 + -IT_0980 + IT_0983 + IT_0984 + -IT_0987 + -IT_0988 +
       -IT_0991 + -IT_0992 + IT_0993 + IT_0996 + IT_0997 + -IT_1000 + -IT_1001 +
       IT_1004 + IT_1005 + -IT_1008 + -IT_1009 + -IT_1012 + -IT_1013 + IT_1014 +
       IT_1015 + -IT_1018 + -IT_1019 + IT_1020 + -IT_1023 + -IT_1024 + -IT_1027 
      + -IT_1028 + -IT_1031 + -IT_1032 + -IT_1035 + -IT_1036 + -IT_1039 + 
      -IT_1040 + -IT_1043 + -IT_1044 + IT_1047 + IT_1048 + IT_1051 + IT_1052 +
       IT_1055 + IT_1056 + IT_1059 + IT_1060 + IT_1063 + IT_1064 + IT_1065;
    const complex_t IT_1075 = IT_1069*IT_1070*IT_1074;
    const complex_t IT_1076 = IT_0004*IT_1075;
    const complex_t IT_1077 = 0.666666666666667*IT_1076;
    return -IT_1073 + IT_1077;
}
} // End of namespace c9_nmfv
