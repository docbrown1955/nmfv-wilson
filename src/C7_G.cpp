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
#include "C7_G.h"
#include "common.h"

#include "params.h"
#include "libcomplexop.h"

namespace c9_nmfv {

complex_t C7_G(
        param_t const &param
        )
{
    clearcache();
    const real_t M_W = param.M_W;
    const real_t g_s = param.g_s;
    const real_t m_b = param.m_b;
    const real_t m_s = param.m_s;
    const real_t V_tb = param.V_tb;
    const real_t e_em = param.e_em;
    const real_t m_sG = param.m_sG;
    const real_t s_12 = param.s_12;
    const real_t m_sb_L = param.m_sb_L;
    const real_t m_sb_R = param.m_sb_R;
    const real_t m_sd_L = param.m_sd_L;
    const real_t m_sd_R = param.m_sd_R;
    const real_t m_ss_L = param.m_ss_L;
    const real_t m_ss_R = param.m_ss_R;
    const real_t theta_W = param.theta_W;
    const complex_t V_ts = param.V_ts;
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
    const complex_t IT_0001 = powq(m_b, -1);
    const complex_t IT_0002 = powq(V_tb, -1);
    const complex_t IT_0003 = cpowq(conjq(V_ts), -1);
    const complex_t IT_0004 = powq(e_em, -3);
    const complex_t IT_0005 = 9.86960440108936*IT_0000*IT_0001*IT_0002*IT_0003
      *IT_0004;
    const complex_t IT_0006 = (complex_t{0, 1.4142135623731})*g_s*conjq
      (U_sd_20);
    const complex_t IT_0007 = (complex_t{0, 1.4142135623731})*g_s*U_sd_40;
    const complex_t IT_0008 = IT_0006*IT_0007;
    const complex_t IT_0009 = 0.101321183642338*m_sG;
    const complex_t IT_0010 = IT_0008*IT_0009;
    const complex_t IT_0011 = cosq(theta_W);
    const complex_t IT_0012 = cpowq(IT_0011, -1);
    const complex_t IT_0013 = tanq(theta_W);
    const complex_t IT_0014 = cpowq(IT_0013, 2);
    const complex_t IT_0015 = cpowq(1 + IT_0014, (-0.5));
    const complex_t IT_0016 = (complex_t{0, 1})*e_em*IT_0012*IT_0015;
    const complex_t IT_0017 = 0.333333333333333*IT_0016;
    const complex_t IT_0018 = powq(m_b, 2);
    const complex_t IT_0019 = powq(m_s, 2);
    const complex_t IT_0020 = powq(m_sG, 2);
    const complex_t IT_0021 = powq(m_sd_L, 2);
    const complex_t IT_0022 = mty::lt::C0iC(0, IT_0018, IT_0018 + IT_0019 + (
      -2)*s_12, IT_0019, IT_0020, IT_0021, IT_0021, mty::lt::reg_int);
    const complex_t IT_0023 = IT_0017*IT_0022;
    const complex_t IT_0024 = 0.666666666666667*IT_0016;
    const complex_t IT_0025 = mty::lt::C0iC(3, IT_0018, IT_0018 + IT_0019 + (
      -2)*s_12, IT_0019, IT_0020, IT_0021, IT_0021, mty::lt::reg_int);
    const complex_t IT_0026 = IT_0024*IT_0025;
    const complex_t IT_0027 = IT_0023 + IT_0026;
    const complex_t IT_0028 = IT_0010*IT_0027;
    const complex_t IT_0029 = (complex_t{0, 1.4142135623731})*g_s*conjq
      (U_sd_50);
    const complex_t IT_0030 = (complex_t{0, 1.4142135623731})*g_s*U_sd_10;
    const complex_t IT_0031 = IT_0029*IT_0030;
    const complex_t IT_0032 = IT_0009*IT_0031;
    const complex_t IT_0033 = IT_0027*IT_0032;
    const complex_t IT_0034 = (complex_t{0, 1.4142135623731})*g_s*conjq
      (U_sd_21);
    const complex_t IT_0035 = (complex_t{0, 1.4142135623731})*g_s*U_sd_41;
    const complex_t IT_0036 = IT_0034*IT_0035;
    const complex_t IT_0037 = IT_0009*IT_0036;
    const complex_t IT_0038 = powq(m_ss_L, 2);
    const complex_t IT_0039 = mty::lt::C0iC(0, IT_0018, IT_0018 + IT_0019 + (
      -2)*s_12, IT_0019, IT_0020, IT_0038, IT_0038, mty::lt::reg_int);
    const complex_t IT_0040 = IT_0017*IT_0039;
    const complex_t IT_0041 = mty::lt::C0iC(3, IT_0018, IT_0018 + IT_0019 + (
      -2)*s_12, IT_0019, IT_0020, IT_0038, IT_0038, mty::lt::reg_int);
    const complex_t IT_0042 = IT_0024*IT_0041;
    const complex_t IT_0043 = IT_0040 + IT_0042;
    const complex_t IT_0044 = IT_0037*IT_0043;
    const complex_t IT_0045 = (complex_t{0, 1.4142135623731})*g_s*conjq
      (U_sd_51);
    const complex_t IT_0046 = (complex_t{0, 1.4142135623731})*g_s*U_sd_11;
    const complex_t IT_0047 = IT_0045*IT_0046;
    const complex_t IT_0048 = IT_0009*IT_0047;
    const complex_t IT_0049 = IT_0043*IT_0048;
    const complex_t IT_0050 = (complex_t{0, 1.4142135623731})*g_s*conjq
      (U_sd_22);
    const complex_t IT_0051 = (complex_t{0, 1.4142135623731})*g_s*U_sd_42;
    const complex_t IT_0052 = IT_0050*IT_0051;
    const complex_t IT_0053 = IT_0009*IT_0052;
    const complex_t IT_0054 = powq(m_sb_L, 2);
    const complex_t IT_0055 = mty::lt::C0iC(0, IT_0018, IT_0018 + IT_0019 + (
      -2)*s_12, IT_0019, IT_0020, IT_0054, IT_0054, mty::lt::reg_int);
    const complex_t IT_0056 = IT_0017*IT_0055;
    const complex_t IT_0057 = mty::lt::C0iC(3, IT_0018, IT_0018 + IT_0019 + (
      -2)*s_12, IT_0019, IT_0020, IT_0054, IT_0054, mty::lt::reg_int);
    const complex_t IT_0058 = IT_0024*IT_0057;
    const complex_t IT_0059 = IT_0056 + IT_0058;
    const complex_t IT_0060 = IT_0053*IT_0059;
    const complex_t IT_0061 = (complex_t{0, 1.4142135623731})*g_s*conjq
      (U_sd_52);
    const complex_t IT_0062 = (complex_t{0, 1.4142135623731})*g_s*U_sd_12;
    const complex_t IT_0063 = IT_0061*IT_0062;
    const complex_t IT_0064 = IT_0009*IT_0063;
    const complex_t IT_0065 = IT_0059*IT_0064;
    const complex_t IT_0066 = (complex_t{0, 1.4142135623731})*g_s*conjq
      (U_sd_23);
    const complex_t IT_0067 = (complex_t{0, 1.4142135623731})*g_s*U_sd_43;
    const complex_t IT_0068 = IT_0066*IT_0067;
    const complex_t IT_0069 = IT_0009*IT_0068;
    const complex_t IT_0070 = powq(m_sd_R, 2);
    const complex_t IT_0071 = mty::lt::C0iC(0, IT_0018, IT_0018 + IT_0019 + (
      -2)*s_12, IT_0019, IT_0020, IT_0070, IT_0070, mty::lt::reg_int);
    const complex_t IT_0072 = IT_0017*IT_0071;
    const complex_t IT_0073 = mty::lt::C0iC(3, IT_0018, IT_0018 + IT_0019 + (
      -2)*s_12, IT_0019, IT_0020, IT_0070, IT_0070, mty::lt::reg_int);
    const complex_t IT_0074 = IT_0024*IT_0073;
    const complex_t IT_0075 = IT_0072 + IT_0074;
    const complex_t IT_0076 = IT_0069*IT_0075;
    const complex_t IT_0077 = (complex_t{0, 1.4142135623731})*g_s*conjq
      (U_sd_53);
    const complex_t IT_0078 = (complex_t{0, 1.4142135623731})*g_s*U_sd_13;
    const complex_t IT_0079 = IT_0077*IT_0078;
    const complex_t IT_0080 = IT_0009*IT_0079;
    const complex_t IT_0081 = IT_0075*IT_0080;
    const complex_t IT_0082 = (complex_t{0, 1.4142135623731})*g_s*conjq
      (U_sd_24);
    const complex_t IT_0083 = (complex_t{0, 1.4142135623731})*g_s*U_sd_44;
    const complex_t IT_0084 = IT_0082*IT_0083;
    const complex_t IT_0085 = IT_0009*IT_0084;
    const complex_t IT_0086 = powq(m_ss_R, 2);
    const complex_t IT_0087 = mty::lt::C0iC(0, IT_0018, IT_0018 + IT_0019 + (
      -2)*s_12, IT_0019, IT_0020, IT_0086, IT_0086, mty::lt::reg_int);
    const complex_t IT_0088 = IT_0017*IT_0087;
    const complex_t IT_0089 = mty::lt::C0iC(3, IT_0018, IT_0018 + IT_0019 + (
      -2)*s_12, IT_0019, IT_0020, IT_0086, IT_0086, mty::lt::reg_int);
    const complex_t IT_0090 = IT_0024*IT_0089;
    const complex_t IT_0091 = IT_0088 + IT_0090;
    const complex_t IT_0092 = IT_0085*IT_0091;
    const complex_t IT_0093 = IT_0034*IT_0046;
    const complex_t IT_0094 = 0.101321183642338*IT_0093;
    const complex_t IT_0095 = mty::lt::C0iC(6, IT_0018, IT_0018 + IT_0019 + (
      -2)*s_12, IT_0019, IT_0020, IT_0038, IT_0038, mty::lt::reg_int);
    const complex_t IT_0096 = IT_0017*IT_0095;
    const complex_t IT_0097 = m_s*IT_0096;
    const complex_t IT_0098 = mty::lt::C0iC(15, IT_0018, IT_0018 + IT_0019 + (
      -2)*s_12, IT_0019, IT_0020, IT_0038, IT_0038, mty::lt::reg_int);
    const complex_t IT_0099 = IT_0024*IT_0098;
    const complex_t IT_0100 = m_s*IT_0099;
    const complex_t IT_0101 = IT_0017*IT_0041;
    const complex_t IT_0102 = m_b*IT_0101;
    const complex_t IT_0103 = mty::lt::C0iC(12, IT_0018, IT_0018 + IT_0019 + (
      -2)*s_12, IT_0019, IT_0020, IT_0038, IT_0038, mty::lt::reg_int);
    const complex_t IT_0104 = IT_0024*IT_0103;
    const complex_t IT_0105 = m_b*IT_0104;
    const complex_t IT_0106 = IT_0097 + IT_0100 + IT_0102 + IT_0105;
    const complex_t IT_0107 = IT_0094*IT_0106;
    const complex_t IT_0108 = IT_0035*IT_0045;
    const complex_t IT_0109 = 0.101321183642338*IT_0108;
    const complex_t IT_0110 = IT_0106*IT_0109;
    const complex_t IT_0111 = IT_0050*IT_0062;
    const complex_t IT_0112 = 0.101321183642338*IT_0111;
    const complex_t IT_0113 = mty::lt::C0iC(6, IT_0018, IT_0018 + IT_0019 + (
      -2)*s_12, IT_0019, IT_0020, IT_0054, IT_0054, mty::lt::reg_int);
    const complex_t IT_0114 = IT_0017*IT_0113;
    const complex_t IT_0115 = m_s*IT_0114;
    const complex_t IT_0116 = mty::lt::C0iC(15, IT_0018, IT_0018 + IT_0019 + (
      -2)*s_12, IT_0019, IT_0020, IT_0054, IT_0054, mty::lt::reg_int);
    const complex_t IT_0117 = IT_0024*IT_0116;
    const complex_t IT_0118 = m_s*IT_0117;
    const complex_t IT_0119 = IT_0017*IT_0057;
    const complex_t IT_0120 = m_b*IT_0119;
    const complex_t IT_0121 = mty::lt::C0iC(12, IT_0018, IT_0018 + IT_0019 + (
      -2)*s_12, IT_0019, IT_0020, IT_0054, IT_0054, mty::lt::reg_int);
    const complex_t IT_0122 = IT_0024*IT_0121;
    const complex_t IT_0123 = m_b*IT_0122;
    const complex_t IT_0124 = IT_0115 + IT_0118 + IT_0120 + IT_0123;
    const complex_t IT_0125 = IT_0112*IT_0124;
    const complex_t IT_0126 = IT_0051*IT_0061;
    const complex_t IT_0127 = 0.101321183642338*IT_0126;
    const complex_t IT_0128 = IT_0124*IT_0127;
    const complex_t IT_0129 = IT_0067*IT_0077;
    const complex_t IT_0130 = 0.101321183642338*IT_0129;
    const complex_t IT_0131 = mty::lt::C0iC(6, IT_0018, IT_0018 + IT_0019 + (
      -2)*s_12, IT_0019, IT_0020, IT_0070, IT_0070, mty::lt::reg_int);
    const complex_t IT_0132 = IT_0017*IT_0131;
    const complex_t IT_0133 = m_s*IT_0132;
    const complex_t IT_0134 = mty::lt::C0iC(15, IT_0018, IT_0018 + IT_0019 + (
      -2)*s_12, IT_0019, IT_0020, IT_0070, IT_0070, mty::lt::reg_int);
    const complex_t IT_0135 = IT_0024*IT_0134;
    const complex_t IT_0136 = m_s*IT_0135;
    const complex_t IT_0137 = IT_0017*IT_0073;
    const complex_t IT_0138 = m_b*IT_0137;
    const complex_t IT_0139 = mty::lt::C0iC(12, IT_0018, IT_0018 + IT_0019 + (
      -2)*s_12, IT_0019, IT_0020, IT_0070, IT_0070, mty::lt::reg_int);
    const complex_t IT_0140 = IT_0024*IT_0139;
    const complex_t IT_0141 = m_b*IT_0140;
    const complex_t IT_0142 = IT_0133 + IT_0136 + IT_0138 + IT_0141;
    const complex_t IT_0143 = IT_0130*IT_0142;
    const complex_t IT_0144 = (complex_t{0, 1.4142135623731})*g_s*conjq
      (U_sd_55);
    const complex_t IT_0145 = (complex_t{0, 1.4142135623731})*g_s*U_sd_45;
    const complex_t IT_0146 = IT_0144*IT_0145;
    const complex_t IT_0147 = 0.101321183642338*IT_0146;
    const complex_t IT_0148 = powq(m_sb_R, 2);
    const complex_t IT_0149 = mty::lt::C0iC(6, IT_0018, IT_0018 + IT_0019 + (
      -2)*s_12, IT_0019, IT_0020, IT_0148, IT_0148, mty::lt::reg_int);
    const complex_t IT_0150 = IT_0017*IT_0149;
    const complex_t IT_0151 = m_s*IT_0150;
    const complex_t IT_0152 = mty::lt::C0iC(15, IT_0018, IT_0018 + IT_0019 + (
      -2)*s_12, IT_0019, IT_0020, IT_0148, IT_0148, mty::lt::reg_int);
    const complex_t IT_0153 = IT_0024*IT_0152;
    const complex_t IT_0154 = m_s*IT_0153;
    const complex_t IT_0155 = mty::lt::C0iC(3, IT_0018, IT_0018 + IT_0019 + (
      -2)*s_12, IT_0019, IT_0020, IT_0148, IT_0148, mty::lt::reg_int);
    const complex_t IT_0156 = IT_0017*IT_0155;
    const complex_t IT_0157 = m_b*IT_0156;
    const complex_t IT_0158 = mty::lt::C0iC(12, IT_0018, IT_0018 + IT_0019 + (
      -2)*s_12, IT_0019, IT_0020, IT_0148, IT_0148, mty::lt::reg_int);
    const complex_t IT_0159 = IT_0024*IT_0158;
    const complex_t IT_0160 = m_b*IT_0159;
    const complex_t IT_0161 = IT_0151 + IT_0154 + IT_0157 + IT_0160;
    const complex_t IT_0162 = IT_0147*IT_0161;
    const complex_t IT_0163 = IT_0006*IT_0030;
    const complex_t IT_0164 = 0.101321183642338*IT_0163;
    const complex_t IT_0165 = mty::lt::C0iC(15, IT_0018, IT_0018 + IT_0019 + (
      -2)*s_12, IT_0019, IT_0020, IT_0021, IT_0021, mty::lt::reg_int);
    const complex_t IT_0166 = IT_0024*IT_0165;
    const complex_t IT_0167 = m_s*IT_0166;
    const complex_t IT_0168 = IT_0017*IT_0025;
    const complex_t IT_0169 = m_b*IT_0168;
    const complex_t IT_0170 = mty::lt::C0iC(12, IT_0018, IT_0018 + IT_0019 + (
      -2)*s_12, IT_0019, IT_0020, IT_0021, IT_0021, mty::lt::reg_int);
    const complex_t IT_0171 = IT_0024*IT_0170;
    const complex_t IT_0172 = m_b*IT_0171;
    const complex_t IT_0173 = mty::lt::C0iC(6, IT_0018, IT_0018 + IT_0019 + (
      -2)*s_12, IT_0019, IT_0020, IT_0021, IT_0021, mty::lt::reg_int);
    const complex_t IT_0174 = IT_0017*IT_0173;
    const complex_t IT_0175 = m_s*IT_0174;
    const complex_t IT_0176 = IT_0167 + IT_0169 + IT_0172 + IT_0175;
    const complex_t IT_0177 = IT_0164*IT_0176;
    const complex_t IT_0178 = IT_0007*IT_0029;
    const complex_t IT_0179 = 0.101321183642338*IT_0178;
    const complex_t IT_0180 = IT_0176*IT_0179;
    const complex_t IT_0181 = IT_0066*IT_0078;
    const complex_t IT_0182 = 0.101321183642338*IT_0181;
    const complex_t IT_0183 = IT_0142*IT_0182;
    const complex_t IT_0184 = (complex_t{0, 1.4142135623731})*g_s*U_sd_14;
    const complex_t IT_0185 = IT_0082*IT_0184;
    const complex_t IT_0186 = 0.101321183642338*IT_0185;
    const complex_t IT_0187 = mty::lt::C0iC(6, IT_0018, IT_0018 + IT_0019 + (
      -2)*s_12, IT_0019, IT_0020, IT_0086, IT_0086, mty::lt::reg_int);
    const complex_t IT_0188 = IT_0017*IT_0187;
    const complex_t IT_0189 = m_s*IT_0188;
    const complex_t IT_0190 = mty::lt::C0iC(15, IT_0018, IT_0018 + IT_0019 + (
      -2)*s_12, IT_0019, IT_0020, IT_0086, IT_0086, mty::lt::reg_int);
    const complex_t IT_0191 = IT_0024*IT_0190;
    const complex_t IT_0192 = m_s*IT_0191;
    const complex_t IT_0193 = mty::lt::C0iC(12, IT_0018, IT_0018 + IT_0019 + (
      -2)*s_12, IT_0019, IT_0020, IT_0086, IT_0086, mty::lt::reg_int);
    const complex_t IT_0194 = IT_0024*IT_0193;
    const complex_t IT_0195 = m_b*IT_0194;
    const complex_t IT_0196 = IT_0017*IT_0089;
    const complex_t IT_0197 = m_b*IT_0196;
    const complex_t IT_0198 = IT_0189 + IT_0192 + IT_0195 + IT_0197;
    const complex_t IT_0199 = IT_0186*IT_0198;
    const complex_t IT_0200 = (complex_t{0, 1.4142135623731})*g_s*conjq
      (U_sd_54);
    const complex_t IT_0201 = IT_0083*IT_0200;
    const complex_t IT_0202 = 0.101321183642338*IT_0201;
    const complex_t IT_0203 = IT_0198*IT_0202;
    const complex_t IT_0204 = IT_0184*IT_0200;
    const complex_t IT_0205 = IT_0009*IT_0204;
    const complex_t IT_0206 = IT_0091*IT_0205;
    const complex_t IT_0207 = (complex_t{0, 1.4142135623731})*g_s*conjq
      (U_sd_25);
    const complex_t IT_0208 = IT_0145*IT_0207;
    const complex_t IT_0209 = IT_0009*IT_0208;
    const complex_t IT_0210 = mty::lt::C0iC(0, IT_0018, IT_0018 + IT_0019 + (
      -2)*s_12, IT_0019, IT_0020, IT_0148, IT_0148, mty::lt::reg_int);
    const complex_t IT_0211 = IT_0017*IT_0210;
    const complex_t IT_0212 = IT_0024*IT_0155;
    const complex_t IT_0213 = IT_0211 + IT_0212;
    const complex_t IT_0214 = IT_0209*IT_0213;
    const complex_t IT_0215 = (complex_t{0, 1.4142135623731})*g_s*U_sd_15;
    const complex_t IT_0216 = IT_0207*IT_0215;
    const complex_t IT_0217 = 0.101321183642338*IT_0216;
    const complex_t IT_0218 = IT_0161*IT_0217;
    const complex_t IT_0219 = IT_0144*IT_0215;
    const complex_t IT_0220 = IT_0009*IT_0219;
    const complex_t IT_0221 = IT_0213*IT_0220;
    const complex_t IT_0222 = IT_0028 + IT_0033 + IT_0044 + IT_0049 + IT_0060 
      + IT_0065 + IT_0076 + IT_0081 + IT_0092 + IT_0107 + IT_0110 + IT_0125 +
       IT_0128 + IT_0143 + IT_0162 + IT_0177 + IT_0180 + IT_0183 + IT_0199 +
       IT_0203 + IT_0206 + IT_0214 + IT_0218 + IT_0221;
    const complex_t IT_0223 = sinq(theta_W);
    const complex_t IT_0224 = cpowq(IT_0223, 2);
    const complex_t IT_0225 = IT_0222*IT_0224;
    const complex_t IT_0226 = IT_0005*IT_0225;
    const complex_t IT_0227 = IT_0022*IT_0024;
    const complex_t IT_0228 = IT_0024*IT_0173;
    const complex_t IT_0229 = -IT_0227 + -IT_0228;
    const complex_t IT_0230 = IT_0023 + IT_0229;
    const complex_t IT_0231 = IT_0032*IT_0230;
    const complex_t IT_0232 = IT_0024*IT_0039;
    const complex_t IT_0233 = IT_0024*IT_0095;
    const complex_t IT_0234 = -IT_0232 + -IT_0233;
    const complex_t IT_0235 = IT_0040 + IT_0234;
    const complex_t IT_0236 = IT_0037*IT_0235;
    const complex_t IT_0237 = IT_0097 + IT_0102;
    const complex_t IT_0238 = m_s*IT_0233;
    const complex_t IT_0239 = mty::lt::C0iC(18, IT_0018, IT_0018 + IT_0019 + (
      -2)*s_12, IT_0019, IT_0020, IT_0038, IT_0038, mty::lt::reg_int);
    const complex_t IT_0240 = IT_0024*IT_0239;
    const complex_t IT_0241 = m_s*IT_0240;
    const complex_t IT_0242 = m_b*IT_0042;
    const complex_t IT_0243 = m_b*IT_0099;
    const complex_t IT_0244 = -IT_0238 + -IT_0241 + -IT_0242 + -IT_0243;
    const complex_t IT_0245 = IT_0237 + IT_0244;
    const complex_t IT_0246 = IT_0094*IT_0245;
    const complex_t IT_0247 = IT_0109*IT_0245;
    const complex_t IT_0248 = IT_0048*IT_0235;
    const complex_t IT_0249 = IT_0024*IT_0055;
    const complex_t IT_0250 = IT_0024*IT_0113;
    const complex_t IT_0251 = -IT_0249 + -IT_0250;
    const complex_t IT_0252 = IT_0056 + IT_0251;
    const complex_t IT_0253 = IT_0053*IT_0252;
    const complex_t IT_0254 = IT_0115 + IT_0120;
    const complex_t IT_0255 = m_s*IT_0250;
    const complex_t IT_0256 = mty::lt::C0iC(18, IT_0018, IT_0018 + IT_0019 + (
      -2)*s_12, IT_0019, IT_0020, IT_0054, IT_0054, mty::lt::reg_int);
    const complex_t IT_0257 = IT_0024*IT_0256;
    const complex_t IT_0258 = m_s*IT_0257;
    const complex_t IT_0259 = m_b*IT_0058;
    const complex_t IT_0260 = m_b*IT_0117;
    const complex_t IT_0261 = -IT_0255 + -IT_0258 + -IT_0259 + -IT_0260;
    const complex_t IT_0262 = IT_0254 + IT_0261;
    const complex_t IT_0263 = IT_0112*IT_0262;
    const complex_t IT_0264 = IT_0127*IT_0262;
    const complex_t IT_0265 = IT_0064*IT_0252;
    const complex_t IT_0266 = IT_0024*IT_0071;
    const complex_t IT_0267 = IT_0024*IT_0131;
    const complex_t IT_0268 = -IT_0266 + -IT_0267;
    const complex_t IT_0269 = IT_0072 + IT_0268;
    const complex_t IT_0270 = IT_0069*IT_0269;
    const complex_t IT_0271 = IT_0133 + IT_0138;
    const complex_t IT_0272 = m_s*IT_0267;
    const complex_t IT_0273 = mty::lt::C0iC(18, IT_0018, IT_0018 + IT_0019 + (
      -2)*s_12, IT_0019, IT_0020, IT_0070, IT_0070, mty::lt::reg_int);
    const complex_t IT_0274 = IT_0024*IT_0273;
    const complex_t IT_0275 = m_s*IT_0274;
    const complex_t IT_0276 = m_b*IT_0074;
    const complex_t IT_0277 = m_b*IT_0135;
    const complex_t IT_0278 = -IT_0272 + -IT_0275 + -IT_0276 + -IT_0277;
    const complex_t IT_0279 = IT_0271 + IT_0278;
    const complex_t IT_0280 = IT_0130*IT_0279;
    const complex_t IT_0281 = IT_0080*IT_0269;
    const complex_t IT_0282 = IT_0024*IT_0087;
    const complex_t IT_0283 = IT_0024*IT_0187;
    const complex_t IT_0284 = -IT_0282 + -IT_0283;
    const complex_t IT_0285 = IT_0088 + IT_0284;
    const complex_t IT_0286 = IT_0085*IT_0285;
    const complex_t IT_0287 = IT_0151 + IT_0157;
    const complex_t IT_0288 = IT_0024*IT_0149;
    const complex_t IT_0289 = m_s*IT_0288;
    const complex_t IT_0290 = mty::lt::C0iC(18, IT_0018, IT_0018 + IT_0019 + (
      -2)*s_12, IT_0019, IT_0020, IT_0148, IT_0148, mty::lt::reg_int);
    const complex_t IT_0291 = IT_0024*IT_0290;
    const complex_t IT_0292 = m_s*IT_0291;
    const complex_t IT_0293 = m_b*IT_0212;
    const complex_t IT_0294 = m_b*IT_0153;
    const complex_t IT_0295 = -IT_0289 + -IT_0292 + -IT_0293 + -IT_0294;
    const complex_t IT_0296 = IT_0287 + IT_0295;
    const complex_t IT_0297 = IT_0217*IT_0296;
    const complex_t IT_0298 = IT_0147*IT_0296;
    const complex_t IT_0299 = IT_0024*IT_0210;
    const complex_t IT_0300 = -IT_0288 + -IT_0299;
    const complex_t IT_0301 = IT_0211 + IT_0300;
    const complex_t IT_0302 = IT_0220*IT_0301;
    const complex_t IT_0303 = IT_0010*IT_0230;
    const complex_t IT_0304 = IT_0169 + IT_0175;
    const complex_t IT_0305 = m_s*IT_0228;
    const complex_t IT_0306 = mty::lt::C0iC(18, IT_0018, IT_0018 + IT_0019 + (
      -2)*s_12, IT_0019, IT_0020, IT_0021, IT_0021, mty::lt::reg_int);
    const complex_t IT_0307 = IT_0024*IT_0306;
    const complex_t IT_0308 = m_s*IT_0307;
    const complex_t IT_0309 = m_b*IT_0026;
    const complex_t IT_0310 = m_b*IT_0166;
    const complex_t IT_0311 = -IT_0305 + -IT_0308 + -IT_0309 + -IT_0310;
    const complex_t IT_0312 = IT_0304 + IT_0311;
    const complex_t IT_0313 = IT_0164*IT_0312;
    const complex_t IT_0314 = IT_0179*IT_0312;
    const complex_t IT_0315 = IT_0182*IT_0279;
    const complex_t IT_0316 = m_s*IT_0283;
    const complex_t IT_0317 = mty::lt::C0iC(18, IT_0018, IT_0018 + IT_0019 + (
      -2)*s_12, IT_0019, IT_0020, IT_0086, IT_0086, mty::lt::reg_int);
    const complex_t IT_0318 = IT_0024*IT_0317;
    const complex_t IT_0319 = m_s*IT_0318;
    const complex_t IT_0320 = m_b*IT_0090;
    const complex_t IT_0321 = m_b*IT_0191;
    const complex_t IT_0322 = -IT_0316 + -IT_0319 + -IT_0320 + -IT_0321;
    const complex_t IT_0323 = IT_0189 + IT_0197;
    const complex_t IT_0324 = IT_0322 + IT_0323;
    const complex_t IT_0325 = IT_0186*IT_0324;
    const complex_t IT_0326 = IT_0202*IT_0324;
    const complex_t IT_0327 = IT_0205*IT_0285;
    const complex_t IT_0328 = IT_0209*IT_0301;
    const complex_t IT_0329 = IT_0231 + IT_0236 + IT_0246 + IT_0247 + IT_0248 
      + IT_0253 + IT_0263 + IT_0264 + IT_0265 + IT_0270 + IT_0280 + IT_0281 +
       IT_0286 + IT_0297 + IT_0298 + IT_0302 + IT_0303 + IT_0313 + IT_0314 +
       IT_0315 + IT_0325 + IT_0326 + IT_0327 + IT_0328;
    const complex_t IT_0330 = IT_0224*IT_0329;
    const complex_t IT_0331 = IT_0005*IT_0330;
    const complex_t IT_0332 = IT_0167 + IT_0175;
    const complex_t IT_0333 = -IT_0169 + -IT_0172;
    const complex_t IT_0334 = IT_0332 + IT_0333;
    const complex_t IT_0335 = IT_0164*IT_0334;
    const complex_t IT_0336 = IT_0097 + IT_0100;
    const complex_t IT_0337 = -IT_0102 + -IT_0105;
    const complex_t IT_0338 = IT_0336 + IT_0337;
    const complex_t IT_0339 = IT_0094*IT_0338;
    const complex_t IT_0340 = IT_0109*IT_0338;
    const complex_t IT_0341 = IT_0115 + IT_0118;
    const complex_t IT_0342 = -IT_0120 + -IT_0123;
    const complex_t IT_0343 = IT_0341 + IT_0342;
    const complex_t IT_0344 = IT_0112*IT_0343;
    const complex_t IT_0345 = IT_0127*IT_0343;
    const complex_t IT_0346 = IT_0179*IT_0334;
    const complex_t IT_0347 = IT_0133 + IT_0136;
    const complex_t IT_0348 = -IT_0138 + -IT_0141;
    const complex_t IT_0349 = IT_0347 + IT_0348;
    const complex_t IT_0350 = IT_0182*IT_0349;
    const complex_t IT_0351 = IT_0130*IT_0349;
    const complex_t IT_0352 = IT_0189 + IT_0192;
    const complex_t IT_0353 = -IT_0195 + -IT_0197;
    const complex_t IT_0354 = IT_0352 + IT_0353;
    const complex_t IT_0355 = IT_0186*IT_0354;
    const complex_t IT_0356 = IT_0202*IT_0354;
    const complex_t IT_0357 = -IT_0157 + -IT_0160;
    const complex_t IT_0358 = IT_0151 + IT_0154;
    const complex_t IT_0359 = IT_0357 + IT_0358;
    const complex_t IT_0360 = IT_0217*IT_0359;
    const complex_t IT_0361 = IT_0147*IT_0359;
    const complex_t IT_0362 = IT_0028 + -IT_0033 + IT_0044 + -IT_0049 +
       IT_0060 + -IT_0065 + IT_0076 + -IT_0081 + IT_0092 + -IT_0206 + IT_0214 + 
      -IT_0221 + IT_0335 + IT_0339 + -IT_0340 + IT_0344 + -IT_0345 + -IT_0346 +
       IT_0350 + -IT_0351 + IT_0355 + -IT_0356 + IT_0360 + -IT_0361;
    const complex_t IT_0363 = IT_0224*IT_0362;
    const complex_t IT_0364 = IT_0005*IT_0363;
    const complex_t IT_0365 = 1.33333333333333*IT_0364;
    const complex_t IT_0366 = -IT_0169 + -IT_0305 + -IT_0308;
    const complex_t IT_0367 = IT_0175 + IT_0309 + IT_0310;
    const complex_t IT_0368 = IT_0366 + IT_0367;
    const complex_t IT_0369 = IT_0179*IT_0368;
    const complex_t IT_0370 = IT_0097 + IT_0242 + IT_0243;
    const complex_t IT_0371 = -IT_0102 + -IT_0238 + -IT_0241;
    const complex_t IT_0372 = IT_0370 + IT_0371;
    const complex_t IT_0373 = IT_0094*IT_0372;
    const complex_t IT_0374 = IT_0109*IT_0372;
    const complex_t IT_0375 = IT_0115 + IT_0259 + IT_0260;
    const complex_t IT_0376 = -IT_0120 + -IT_0255 + -IT_0258;
    const complex_t IT_0377 = IT_0375 + IT_0376;
    const complex_t IT_0378 = IT_0112*IT_0377;
    const complex_t IT_0379 = IT_0127*IT_0377;
    const complex_t IT_0380 = IT_0133 + IT_0276 + IT_0277;
    const complex_t IT_0381 = -IT_0138 + -IT_0272 + -IT_0275;
    const complex_t IT_0382 = IT_0380 + IT_0381;
    const complex_t IT_0383 = IT_0130*IT_0382;
    const complex_t IT_0384 = IT_0182*IT_0382;
    const complex_t IT_0385 = IT_0164*IT_0368;
    const complex_t IT_0386 = IT_0189 + IT_0320 + IT_0321;
    const complex_t IT_0387 = -IT_0197 + -IT_0316 + -IT_0319;
    const complex_t IT_0388 = IT_0386 + IT_0387;
    const complex_t IT_0389 = IT_0186*IT_0388;
    const complex_t IT_0390 = IT_0202*IT_0388;
    const complex_t IT_0391 = IT_0151 + IT_0293 + IT_0294;
    const complex_t IT_0392 = -IT_0157 + -IT_0289 + -IT_0292;
    const complex_t IT_0393 = IT_0391 + IT_0392;
    const complex_t IT_0394 = IT_0217*IT_0393;
    const complex_t IT_0395 = IT_0147*IT_0393;
    const complex_t IT_0396 = -IT_0231 + IT_0236 + -IT_0248 + IT_0253 + 
      -IT_0265 + IT_0270 + -IT_0281 + IT_0286 + -IT_0302 + IT_0303 + -IT_0327 +
       IT_0328 + -IT_0369 + IT_0373 + -IT_0374 + IT_0378 + -IT_0379 + -IT_0383 +
       IT_0384 + IT_0385 + IT_0389 + -IT_0390 + IT_0394 + -IT_0395;
    const complex_t IT_0397 = IT_0224*IT_0396;
    const complex_t IT_0398 = -IT_0397;
    const complex_t IT_0399 = -IT_0398;
    const complex_t IT_0400 = IT_0005*IT_0399;
    const complex_t IT_0401 = (-1.33333333333333)*IT_0400;
    return (complex_t{0, (-0.333333333333333)})*IT_0226 + (complex_t{0,
       0.333333333333333})*IT_0331 + (complex_t{0, 0.25})*IT_0365 + 
      (complex_t{0, 0.25})*IT_0401;
}
} // End of namespace c9_nmfv
