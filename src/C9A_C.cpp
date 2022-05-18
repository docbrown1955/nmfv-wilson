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
#include "C9A_C.h"
#include "common.h"

#include "params.h"
#include "libcomplexop.h"

namespace c9_nmfv {

complex_t C9A_C(
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
    const real_t m_mu = param.m_mu;
    const real_t s_12 = param.s_12;
    const real_t s_34 = param.s_34;
    const real_t m_C1p = param.m_C1p;
    const real_t m_C2p = param.m_C2p;
    const real_t m_sc_L = param.m_sc_L;
    const real_t m_sc_R = param.m_sc_R;
    const real_t m_st_L = param.m_st_L;
    const real_t m_st_R = param.m_st_R;
    const real_t m_su_L = param.m_su_L;
    const real_t m_su_R = param.m_su_R;
    const real_t theta_W = param.theta_W;
    const real_t V_ub_mod = param.V_ub_mod;
    const real_t reg_prop = param.reg_prop;
    const real_t delta_wolf = param.delta_wolf;
    const complex_t U_d1 = param.U_d1;
    const complex_t U_d2 = param.U_d2;
    const complex_t V_cs = param.V_cs;
    const complex_t V_ts = param.V_ts;
    const complex_t V_u1 = param.V_u1;
    const complex_t V_u2 = param.V_u2;
    const complex_t U_Wm1 = param.U_Wm1;
    const complex_t U_Wm2 = param.U_Wm2;
    const complex_t V_Wp1 = param.V_Wp1;
    const complex_t V_Wp2 = param.V_Wp2;
    const complex_t U_su_00 = param.U_su_00;
    const complex_t U_su_01 = param.U_su_01;
    const complex_t U_su_02 = param.U_su_02;
    const complex_t U_su_03 = param.U_su_03;
    const complex_t U_su_04 = param.U_su_04;
    const complex_t U_su_05 = param.U_su_05;
    const complex_t U_su_10 = param.U_su_10;
    const complex_t U_su_11 = param.U_su_11;
    const complex_t U_su_12 = param.U_su_12;
    const complex_t U_su_13 = param.U_su_13;
    const complex_t U_su_14 = param.U_su_14;
    const complex_t U_su_15 = param.U_su_15;
    const complex_t U_su_20 = param.U_su_20;
    const complex_t U_su_21 = param.U_su_21;
    const complex_t U_su_22 = param.U_su_22;
    const complex_t U_su_23 = param.U_su_23;
    const complex_t U_su_24 = param.U_su_24;
    const complex_t U_su_25 = param.U_su_25;
    const complex_t U_su_30 = param.U_su_30;
    const complex_t U_su_31 = param.U_su_31;
    const complex_t U_su_32 = param.U_su_32;
    const complex_t U_su_33 = param.U_su_33;
    const complex_t U_su_34 = param.U_su_34;
    const complex_t U_su_35 = param.U_su_35;
    const complex_t U_su_40 = param.U_su_40;
    const complex_t U_su_41 = param.U_su_41;
    const complex_t U_su_42 = param.U_su_42;
    const complex_t U_su_43 = param.U_su_43;
    const complex_t U_su_44 = param.U_su_44;
    const complex_t U_su_45 = param.U_su_45;
    const complex_t U_su_50 = param.U_su_50;
    const complex_t U_su_51 = param.U_su_51;
    const complex_t U_su_52 = param.U_su_52;
    const complex_t U_su_53 = param.U_su_53;
    const complex_t U_su_54 = param.U_su_54;
    const complex_t U_su_55 = param.U_su_55;
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
    const complex_t IT_0017 = 0.101321183642338*m_C1p;
    const complex_t IT_0018 = sinq(theta_W);
    const complex_t IT_0019 = cpowq(IT_0018, -1);
    const complex_t IT_0020 = V_cb*e_em*V_Wp1*conjq(U_su_12);
    const complex_t IT_0021 = IT_0019*IT_0020;
    const complex_t IT_0022 = V_tb*e_em*V_Wp1*conjq(U_su_22);
    const complex_t IT_0023 = IT_0019*IT_0022;
    const complex_t IT_0024 = cexpq((complex_t{0, -1})*delta_wolf);
    const complex_t IT_0025 = IT_0019*IT_0024;
    const complex_t IT_0026 = e_em*V_Wp1*conjq(U_su_02)*V_ub_mod;
    const complex_t IT_0027 = IT_0025*IT_0026;
    const complex_t IT_0028 = sinq(beta);
    const complex_t IT_0029 = cpowq(IT_0028, -1);
    const complex_t IT_0030 = IT_0019*IT_0029;
    const complex_t IT_0031 = powq(M_W, -1);
    const complex_t IT_0032 = m_c*V_cb*V_u1*e_em*IT_0031*conjq(U_su_42);
    const complex_t IT_0033 = IT_0030*IT_0032;
    const complex_t IT_0034 = 1.4142135623731*IT_0033;
    const complex_t IT_0035 = m_t*V_tb*V_u1*e_em*IT_0031*conjq(U_su_52);
    const complex_t IT_0036 = IT_0030*IT_0035;
    const complex_t IT_0037 = 1.4142135623731*IT_0036;
    const complex_t IT_0038 = IT_0019*IT_0024*IT_0029;
    const complex_t IT_0039 = m_u*V_u1*e_em*IT_0031*conjq(U_su_32)*V_ub_mod;
    const complex_t IT_0040 = IT_0038*IT_0039;
    const complex_t IT_0041 = 1.4142135623731*IT_0040;
    const complex_t IT_0042 = (complex_t{0, 1})*(IT_0021 + IT_0023 + IT_0027 +
       (-0.5)*IT_0034 + (-0.5)*IT_0037 + (-0.5)*IT_0041);
    const complex_t IT_0043 = cosq(beta);
    const complex_t IT_0044 = cpowq(IT_0043, -1);
    const complex_t IT_0045 = IT_0019*IT_0044;
    const complex_t IT_0046 = m_s*U_d1*conjq(V_cs)*e_em*IT_0031*U_su_12;
    const complex_t IT_0047 = IT_0045*IT_0046;
    const complex_t IT_0048 = 1.4142135623731*IT_0047;
    const complex_t IT_0049 = m_s*U_d1*conjq(V_ts)*e_em*IT_0031*U_su_22;
    const complex_t IT_0050 = IT_0045*IT_0049;
    const complex_t IT_0051 = 1.4142135623731*IT_0050;
    const complex_t IT_0052 = m_s*U_d1*V_us*e_em*IT_0031*U_su_02;
    const complex_t IT_0053 = IT_0045*IT_0052;
    const complex_t IT_0054 = 1.4142135623731*IT_0053;
    const complex_t IT_0055 = (complex_t{0, 1})*(IT_0048 + IT_0051 + IT_0054);
    const complex_t IT_0056 = 0.5*IT_0055;
    const complex_t IT_0057 = IT_0042*IT_0056;
    const complex_t IT_0058 = IT_0017*IT_0057;
    const complex_t IT_0059 = (-0.666666666666667)*IT_0006;
    const complex_t IT_0060 = powq(m_C1p, 2);
    const complex_t IT_0061 = powq(m_st_L, 2);
    const complex_t IT_0062 = mty::lt::C0iC(0, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0060, IT_0061, IT_0061, mty::lt::reg_int);
    const complex_t IT_0063 = IT_0059*IT_0062;
    const complex_t IT_0064 = (-1.33333333333333)*IT_0006;
    const complex_t IT_0065 = mty::lt::C0iC(3, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0060, IT_0061, IT_0061, mty::lt::reg_int);
    const complex_t IT_0066 = IT_0064*IT_0065;
    const complex_t IT_0067 = IT_0063 + IT_0066;
    const complex_t IT_0068 = IT_0058*IT_0067;
    const complex_t IT_0069 = IT_0002*IT_0005;
    const complex_t IT_0070 = e_em*U_Wm1*conjq(U_Wm2);
    const complex_t IT_0071 = IT_0069*IT_0070;
    const complex_t IT_0072 = U_d1*conjq(U_d2)*e_em;
    const complex_t IT_0073 = IT_0069*IT_0072;
    const complex_t IT_0074 = (complex_t{0, 1})*(IT_0071 + IT_0073);
    const complex_t IT_0075 = -IT_0074;
    const complex_t IT_0076 = m_s*U_d2*conjq(V_cs)*e_em*IT_0031*U_su_12;
    const complex_t IT_0077 = IT_0045*IT_0076;
    const complex_t IT_0078 = 1.4142135623731*IT_0077;
    const complex_t IT_0079 = m_s*U_d2*conjq(V_ts)*e_em*IT_0031*U_su_22;
    const complex_t IT_0080 = IT_0045*IT_0079;
    const complex_t IT_0081 = 1.4142135623731*IT_0080;
    const complex_t IT_0082 = m_s*U_d2*V_us*e_em*IT_0031*U_su_02;
    const complex_t IT_0083 = IT_0045*IT_0082;
    const complex_t IT_0084 = 1.4142135623731*IT_0083;
    const complex_t IT_0085 = (complex_t{0, 1})*(IT_0078 + IT_0081 + IT_0084);
    const complex_t IT_0086 = 0.5*IT_0085;
    const complex_t IT_0087 = IT_0042*IT_0075*IT_0086;
    const complex_t IT_0088 = IT_0017*IT_0087;
    const complex_t IT_0089 = powq(m_C2p, 2);
    const complex_t IT_0090 = mty::lt::C0iC(3, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0089, IT_0061, mty::lt::reg_int);
    const complex_t IT_0091 = mty::lt::C0iC(6, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0089, IT_0061, mty::lt::reg_int);
    const complex_t IT_0092 = mty::lt::C0iC(0, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0089, IT_0061, mty::lt::reg_int);
    const complex_t IT_0093 = IT_0090 + IT_0091 + IT_0092;
    const complex_t IT_0094 = IT_0088*IT_0093;
    const complex_t IT_0095 = V_us*e_em*conjq(V_Wp1)*U_su_02;
    const complex_t IT_0096 = IT_0019*IT_0095;
    const complex_t IT_0097 = conjq(V_cs)*e_em*conjq(V_Wp1)*U_su_12;
    const complex_t IT_0098 = IT_0019*IT_0097;
    const complex_t IT_0099 = conjq(V_ts)*e_em*conjq(V_Wp1)*U_su_22;
    const complex_t IT_0100 = IT_0019*IT_0099;
    const complex_t IT_0101 = m_u*conjq(V_u1)*V_us*e_em*IT_0031*U_su_32;
    const complex_t IT_0102 = IT_0030*IT_0101;
    const complex_t IT_0103 = 1.4142135623731*IT_0102;
    const complex_t IT_0104 = m_c*conjq(V_cs)*conjq(V_u1)*e_em*IT_0031*U_su_42;
    const complex_t IT_0105 = IT_0030*IT_0104;
    const complex_t IT_0106 = 1.4142135623731*IT_0105;
    const complex_t IT_0107 = m_t*conjq(V_ts)*conjq(V_u1)*e_em*IT_0031*U_su_52;
    const complex_t IT_0108 = IT_0030*IT_0107;
    const complex_t IT_0109 = 1.4142135623731*IT_0108;
    const complex_t IT_0110 = (complex_t{0, 1})*(IT_0096 + IT_0098 + IT_0100 +
       (-0.5)*IT_0103 + (-0.5)*IT_0106 + (-0.5)*IT_0109);
    const complex_t IT_0111 = IT_0042*IT_0110;
    const complex_t IT_0112 = 0.101321183642338*IT_0111;
    const complex_t IT_0113 = mty::lt::C0iC(15, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0060, IT_0061, IT_0061, mty::lt::reg_int);
    const complex_t IT_0114 = IT_0064*IT_0113;
    const complex_t IT_0115 = m_s*IT_0114;
    const complex_t IT_0116 = IT_0059*IT_0065;
    const complex_t IT_0117 = m_b*IT_0116;
    const complex_t IT_0118 = mty::lt::C0iC(12, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0060, IT_0061, IT_0061, mty::lt::reg_int);
    const complex_t IT_0119 = IT_0064*IT_0118;
    const complex_t IT_0120 = m_b*IT_0119;
    const complex_t IT_0121 = mty::lt::C0iC(6, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0060, IT_0061, IT_0061, mty::lt::reg_int);
    const complex_t IT_0122 = IT_0059*IT_0121;
    const complex_t IT_0123 = m_s*IT_0122;
    const complex_t IT_0124 = IT_0115 + IT_0117 + IT_0120 + IT_0123;
    const complex_t IT_0125 = IT_0112*IT_0124;
    const complex_t IT_0126 = V_cb*e_em*V_Wp2*conjq(U_su_12);
    const complex_t IT_0127 = IT_0019*IT_0126;
    const complex_t IT_0128 = V_tb*e_em*V_Wp2*conjq(U_su_22);
    const complex_t IT_0129 = IT_0019*IT_0128;
    const complex_t IT_0130 = e_em*V_Wp2*conjq(U_su_02)*V_ub_mod;
    const complex_t IT_0131 = IT_0025*IT_0130;
    const complex_t IT_0132 = m_c*V_cb*V_u2*e_em*IT_0031*conjq(U_su_42);
    const complex_t IT_0133 = IT_0030*IT_0132;
    const complex_t IT_0134 = 1.4142135623731*IT_0133;
    const complex_t IT_0135 = m_t*V_tb*V_u2*e_em*IT_0031*conjq(U_su_52);
    const complex_t IT_0136 = IT_0030*IT_0135;
    const complex_t IT_0137 = 1.4142135623731*IT_0136;
    const complex_t IT_0138 = m_u*V_u2*e_em*IT_0031*conjq(U_su_32)*V_ub_mod;
    const complex_t IT_0139 = IT_0038*IT_0138;
    const complex_t IT_0140 = 1.4142135623731*IT_0139;
    const complex_t IT_0141 = (complex_t{0, 1})*(IT_0127 + IT_0129 + IT_0131 +
       (-0.5)*IT_0134 + (-0.5)*IT_0137 + (-0.5)*IT_0140);
    const complex_t IT_0142 = V_us*e_em*conjq(V_Wp2)*U_su_02;
    const complex_t IT_0143 = IT_0019*IT_0142;
    const complex_t IT_0144 = conjq(V_cs)*e_em*conjq(V_Wp2)*U_su_12;
    const complex_t IT_0145 = IT_0019*IT_0144;
    const complex_t IT_0146 = conjq(V_ts)*e_em*conjq(V_Wp2)*U_su_22;
    const complex_t IT_0147 = IT_0019*IT_0146;
    const complex_t IT_0148 = m_u*conjq(V_u2)*V_us*e_em*IT_0031*U_su_32;
    const complex_t IT_0149 = IT_0030*IT_0148;
    const complex_t IT_0150 = 1.4142135623731*IT_0149;
    const complex_t IT_0151 = m_c*conjq(V_cs)*conjq(V_u2)*e_em*IT_0031*U_su_42;
    const complex_t IT_0152 = IT_0030*IT_0151;
    const complex_t IT_0153 = 1.4142135623731*IT_0152;
    const complex_t IT_0154 = m_t*conjq(V_ts)*conjq(V_u2)*e_em*IT_0031*U_su_52;
    const complex_t IT_0155 = IT_0030*IT_0154;
    const complex_t IT_0156 = 1.4142135623731*IT_0155;
    const complex_t IT_0157 = (complex_t{0, 1})*(IT_0143 + IT_0145 + IT_0147 +
       (-0.5)*IT_0150 + (-0.5)*IT_0153 + (-0.5)*IT_0156);
    const complex_t IT_0158 = IT_0141*IT_0157;
    const complex_t IT_0159 = 0.101321183642338*IT_0158;
    const complex_t IT_0160 = mty::lt::C0iC(6, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0089, IT_0061, IT_0061, mty::lt::reg_int);
    const complex_t IT_0161 = IT_0059*IT_0160;
    const complex_t IT_0162 = m_s*IT_0161;
    const complex_t IT_0163 = mty::lt::C0iC(3, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0089, IT_0061, IT_0061, mty::lt::reg_int);
    const complex_t IT_0164 = IT_0059*IT_0163;
    const complex_t IT_0165 = m_b*IT_0164;
    const complex_t IT_0166 = mty::lt::C0iC(12, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0089, IT_0061, IT_0061, mty::lt::reg_int);
    const complex_t IT_0167 = IT_0064*IT_0166;
    const complex_t IT_0168 = m_b*IT_0167;
    const complex_t IT_0169 = mty::lt::C0iC(15, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0089, IT_0061, IT_0061, mty::lt::reg_int);
    const complex_t IT_0170 = IT_0064*IT_0169;
    const complex_t IT_0171 = m_s*IT_0170;
    const complex_t IT_0172 = IT_0162 + IT_0165 + IT_0168 + IT_0171;
    const complex_t IT_0173 = IT_0159*IT_0172;
    const complex_t IT_0174 = 0.101321183642338*m_C2p;
    const complex_t IT_0175 = IT_0086*IT_0141;
    const complex_t IT_0176 = IT_0174*IT_0175;
    const complex_t IT_0177 = IT_0064*IT_0163;
    const complex_t IT_0178 = mty::lt::C0iC(0, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0089, IT_0061, IT_0061, mty::lt::reg_int);
    const complex_t IT_0179 = IT_0059*IT_0178;
    const complex_t IT_0180 = IT_0177 + IT_0179;
    const complex_t IT_0181 = IT_0176*IT_0180;
    const complex_t IT_0182 = m_b*conjq(U_d1)*V_cb*e_em*IT_0031*conjq(U_su_10);
    const complex_t IT_0183 = IT_0045*IT_0182;
    const complex_t IT_0184 = 1.4142135623731*IT_0183;
    const complex_t IT_0185 = m_b*conjq(U_d1)*V_tb*e_em*IT_0031*conjq(U_su_20);
    const complex_t IT_0186 = IT_0045*IT_0185;
    const complex_t IT_0187 = 1.4142135623731*IT_0186;
    const complex_t IT_0188 = IT_0019*IT_0024*IT_0044;
    const complex_t IT_0189 = m_b*conjq(U_d1)*e_em*IT_0031*conjq(U_su_00)
      *V_ub_mod;
    const complex_t IT_0190 = IT_0188*IT_0189;
    const complex_t IT_0191 = 1.4142135623731*IT_0190;
    const complex_t IT_0192 = (complex_t{0, 1})*(IT_0184 + IT_0187 + IT_0191);
    const complex_t IT_0193 = 0.5*IT_0192;
    const complex_t IT_0194 = V_us*e_em*conjq(V_Wp1)*U_su_00;
    const complex_t IT_0195 = IT_0019*IT_0194;
    const complex_t IT_0196 = conjq(V_cs)*e_em*conjq(V_Wp1)*U_su_10;
    const complex_t IT_0197 = IT_0019*IT_0196;
    const complex_t IT_0198 = conjq(V_ts)*e_em*conjq(V_Wp1)*U_su_20;
    const complex_t IT_0199 = IT_0019*IT_0198;
    const complex_t IT_0200 = m_u*conjq(V_u1)*V_us*e_em*IT_0031*U_su_30;
    const complex_t IT_0201 = IT_0030*IT_0200;
    const complex_t IT_0202 = 1.4142135623731*IT_0201;
    const complex_t IT_0203 = m_c*conjq(V_cs)*conjq(V_u1)*e_em*IT_0031*U_su_40;
    const complex_t IT_0204 = IT_0030*IT_0203;
    const complex_t IT_0205 = 1.4142135623731*IT_0204;
    const complex_t IT_0206 = m_t*conjq(V_ts)*conjq(V_u1)*e_em*IT_0031*U_su_50;
    const complex_t IT_0207 = IT_0030*IT_0206;
    const complex_t IT_0208 = 1.4142135623731*IT_0207;
    const complex_t IT_0209 = (complex_t{0, 1})*(IT_0195 + IT_0197 + IT_0199 +
       (-0.5)*IT_0202 + (-0.5)*IT_0205 + (-0.5)*IT_0208);
    const complex_t IT_0210 = IT_0193*IT_0209;
    const complex_t IT_0211 = IT_0017*IT_0210;
    const complex_t IT_0212 = powq(m_su_L, 2);
    const complex_t IT_0213 = mty::lt::C0iC(3, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0060, IT_0212, IT_0212, mty::lt::reg_int);
    const complex_t IT_0214 = IT_0064*IT_0213;
    const complex_t IT_0215 = mty::lt::C0iC(0, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0060, IT_0212, IT_0212, mty::lt::reg_int);
    const complex_t IT_0216 = IT_0059*IT_0215;
    const complex_t IT_0217 = IT_0214 + IT_0216;
    const complex_t IT_0218 = IT_0211*IT_0217;
    const complex_t IT_0219 = m_b*conjq(U_d2)*V_cb*e_em*IT_0031*conjq(U_su_10);
    const complex_t IT_0220 = IT_0045*IT_0219;
    const complex_t IT_0221 = 1.4142135623731*IT_0220;
    const complex_t IT_0222 = m_b*conjq(U_d2)*V_tb*e_em*IT_0031*conjq(U_su_20);
    const complex_t IT_0223 = IT_0045*IT_0222;
    const complex_t IT_0224 = 1.4142135623731*IT_0223;
    const complex_t IT_0225 = m_b*conjq(U_d2)*e_em*IT_0031*conjq(U_su_00)
      *V_ub_mod;
    const complex_t IT_0226 = IT_0188*IT_0225;
    const complex_t IT_0227 = 1.4142135623731*IT_0226;
    const complex_t IT_0228 = (complex_t{0, 1})*(IT_0221 + IT_0224 + IT_0227);
    const complex_t IT_0229 = 0.5*IT_0228;
    const complex_t IT_0230 = V_us*e_em*conjq(V_Wp2)*U_su_00;
    const complex_t IT_0231 = IT_0019*IT_0230;
    const complex_t IT_0232 = conjq(V_cs)*e_em*conjq(V_Wp2)*U_su_10;
    const complex_t IT_0233 = IT_0019*IT_0232;
    const complex_t IT_0234 = conjq(V_ts)*e_em*conjq(V_Wp2)*U_su_20;
    const complex_t IT_0235 = IT_0019*IT_0234;
    const complex_t IT_0236 = m_u*conjq(V_u2)*V_us*e_em*IT_0031*U_su_30;
    const complex_t IT_0237 = IT_0030*IT_0236;
    const complex_t IT_0238 = 1.4142135623731*IT_0237;
    const complex_t IT_0239 = m_c*conjq(V_cs)*conjq(V_u2)*e_em*IT_0031*U_su_40;
    const complex_t IT_0240 = IT_0030*IT_0239;
    const complex_t IT_0241 = 1.4142135623731*IT_0240;
    const complex_t IT_0242 = m_t*conjq(V_ts)*conjq(V_u2)*e_em*IT_0031*U_su_50;
    const complex_t IT_0243 = IT_0030*IT_0242;
    const complex_t IT_0244 = 1.4142135623731*IT_0243;
    const complex_t IT_0245 = (complex_t{0, 1})*(IT_0231 + IT_0233 + IT_0235 +
       (-0.5)*IT_0238 + (-0.5)*IT_0241 + (-0.5)*IT_0244);
    const complex_t IT_0246 = IT_0229*IT_0245;
    const complex_t IT_0247 = IT_0174*IT_0246;
    const complex_t IT_0248 = mty::lt::C0iC(0, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0089, IT_0212, IT_0212, mty::lt::reg_int);
    const complex_t IT_0249 = IT_0059*IT_0248;
    const complex_t IT_0250 = mty::lt::C0iC(3, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0089, IT_0212, IT_0212, mty::lt::reg_int);
    const complex_t IT_0251 = IT_0064*IT_0250;
    const complex_t IT_0252 = IT_0249 + IT_0251;
    const complex_t IT_0253 = IT_0247*IT_0252;
    const complex_t IT_0254 = V_cb*e_em*V_Wp1*conjq(U_su_10);
    const complex_t IT_0255 = IT_0019*IT_0254;
    const complex_t IT_0256 = V_tb*e_em*V_Wp1*conjq(U_su_20);
    const complex_t IT_0257 = IT_0019*IT_0256;
    const complex_t IT_0258 = e_em*V_Wp1*conjq(U_su_00)*V_ub_mod;
    const complex_t IT_0259 = IT_0025*IT_0258;
    const complex_t IT_0260 = m_c*V_cb*V_u1*e_em*IT_0031*conjq(U_su_40);
    const complex_t IT_0261 = IT_0030*IT_0260;
    const complex_t IT_0262 = 1.4142135623731*IT_0261;
    const complex_t IT_0263 = m_t*V_tb*V_u1*e_em*IT_0031*conjq(U_su_50);
    const complex_t IT_0264 = IT_0030*IT_0263;
    const complex_t IT_0265 = 1.4142135623731*IT_0264;
    const complex_t IT_0266 = m_u*V_u1*e_em*IT_0031*conjq(U_su_30)*V_ub_mod;
    const complex_t IT_0267 = IT_0038*IT_0266;
    const complex_t IT_0268 = 1.4142135623731*IT_0267;
    const complex_t IT_0269 = (complex_t{0, 1})*(IT_0255 + IT_0257 + IT_0259 +
       (-0.5)*IT_0262 + (-0.5)*IT_0265 + (-0.5)*IT_0268);
    const complex_t IT_0270 = IT_0209*IT_0269;
    const complex_t IT_0271 = 0.101321183642338*IT_0270;
    const complex_t IT_0272 = IT_0059*IT_0213;
    const complex_t IT_0273 = m_b*IT_0272;
    const complex_t IT_0274 = mty::lt::C0iC(12, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0060, IT_0212, IT_0212, mty::lt::reg_int);
    const complex_t IT_0275 = IT_0064*IT_0274;
    const complex_t IT_0276 = m_b*IT_0275;
    const complex_t IT_0277 = mty::lt::C0iC(6, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0060, IT_0212, IT_0212, mty::lt::reg_int);
    const complex_t IT_0278 = IT_0059*IT_0277;
    const complex_t IT_0279 = m_s*IT_0278;
    const complex_t IT_0280 = mty::lt::C0iC(15, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0060, IT_0212, IT_0212, mty::lt::reg_int);
    const complex_t IT_0281 = IT_0064*IT_0280;
    const complex_t IT_0282 = m_s*IT_0281;
    const complex_t IT_0283 = IT_0273 + IT_0276 + IT_0279 + IT_0282;
    const complex_t IT_0284 = IT_0271*IT_0283;
    const complex_t IT_0285 = V_cb*e_em*V_Wp2*conjq(U_su_10);
    const complex_t IT_0286 = IT_0019*IT_0285;
    const complex_t IT_0287 = V_tb*e_em*V_Wp2*conjq(U_su_20);
    const complex_t IT_0288 = IT_0019*IT_0287;
    const complex_t IT_0289 = e_em*V_Wp2*conjq(U_su_00)*V_ub_mod;
    const complex_t IT_0290 = IT_0025*IT_0289;
    const complex_t IT_0291 = m_c*V_cb*V_u2*e_em*IT_0031*conjq(U_su_40);
    const complex_t IT_0292 = IT_0030*IT_0291;
    const complex_t IT_0293 = 1.4142135623731*IT_0292;
    const complex_t IT_0294 = m_t*V_tb*V_u2*e_em*IT_0031*conjq(U_su_50);
    const complex_t IT_0295 = IT_0030*IT_0294;
    const complex_t IT_0296 = 1.4142135623731*IT_0295;
    const complex_t IT_0297 = m_u*V_u2*e_em*IT_0031*conjq(U_su_30)*V_ub_mod;
    const complex_t IT_0298 = IT_0038*IT_0297;
    const complex_t IT_0299 = 1.4142135623731*IT_0298;
    const complex_t IT_0300 = (complex_t{0, 1})*(IT_0286 + IT_0288 + IT_0290 +
       (-0.5)*IT_0293 + (-0.5)*IT_0296 + (-0.5)*IT_0299);
    const complex_t IT_0301 = IT_0245*IT_0300;
    const complex_t IT_0302 = 0.101321183642338*IT_0301;
    const complex_t IT_0303 = mty::lt::C0iC(6, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0089, IT_0212, IT_0212, mty::lt::reg_int);
    const complex_t IT_0304 = IT_0059*IT_0303;
    const complex_t IT_0305 = m_s*IT_0304;
    const complex_t IT_0306 = IT_0059*IT_0250;
    const complex_t IT_0307 = m_b*IT_0306;
    const complex_t IT_0308 = mty::lt::C0iC(12, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0089, IT_0212, IT_0212, mty::lt::reg_int);
    const complex_t IT_0309 = IT_0064*IT_0308;
    const complex_t IT_0310 = m_b*IT_0309;
    const complex_t IT_0311 = mty::lt::C0iC(15, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0089, IT_0212, IT_0212, mty::lt::reg_int);
    const complex_t IT_0312 = IT_0064*IT_0311;
    const complex_t IT_0313 = m_s*IT_0312;
    const complex_t IT_0314 = IT_0305 + IT_0307 + IT_0310 + IT_0313;
    const complex_t IT_0315 = IT_0302*IT_0314;
    const complex_t IT_0316 = e_em*U_Wm1*conjq(U_Wm1);
    const complex_t IT_0317 = IT_0069*IT_0316;
    const complex_t IT_0318 = U_d1*conjq(U_d1)*e_em;
    const complex_t IT_0319 = IT_0069*IT_0318;
    const complex_t IT_0320 = (complex_t{0, 1})*(IT_0317 + IT_0319);
    const complex_t IT_0321 = -IT_0320;
    const complex_t IT_0322 = m_s*U_d1*conjq(V_cs)*e_em*IT_0031*U_su_14;
    const complex_t IT_0323 = IT_0045*IT_0322;
    const complex_t IT_0324 = 1.4142135623731*IT_0323;
    const complex_t IT_0325 = m_s*U_d1*conjq(V_ts)*e_em*IT_0031*U_su_24;
    const complex_t IT_0326 = IT_0045*IT_0325;
    const complex_t IT_0327 = 1.4142135623731*IT_0326;
    const complex_t IT_0328 = m_s*U_d1*V_us*e_em*IT_0031*U_su_04;
    const complex_t IT_0329 = IT_0045*IT_0328;
    const complex_t IT_0330 = 1.4142135623731*IT_0329;
    const complex_t IT_0331 = (complex_t{0, 1})*(IT_0324 + IT_0327 + IT_0330);
    const complex_t IT_0332 = 0.5*IT_0331;
    const complex_t IT_0333 = m_b*conjq(U_d1)*V_cb*e_em*IT_0031*conjq(U_su_14);
    const complex_t IT_0334 = IT_0045*IT_0333;
    const complex_t IT_0335 = 1.4142135623731*IT_0334;
    const complex_t IT_0336 = m_b*conjq(U_d1)*V_tb*e_em*IT_0031*conjq(U_su_24);
    const complex_t IT_0337 = IT_0045*IT_0336;
    const complex_t IT_0338 = 1.4142135623731*IT_0337;
    const complex_t IT_0339 = m_b*conjq(U_d1)*e_em*IT_0031*conjq(U_su_04)
      *V_ub_mod;
    const complex_t IT_0340 = IT_0188*IT_0339;
    const complex_t IT_0341 = 1.4142135623731*IT_0340;
    const complex_t IT_0342 = (complex_t{0, 1})*(IT_0335 + IT_0338 + IT_0341);
    const complex_t IT_0343 = 0.5*IT_0342;
    const complex_t IT_0344 = IT_0321*IT_0332*IT_0343;
    const complex_t IT_0345 = 0.101321183642338*IT_0344;
    const complex_t IT_0346 = powq(m_sc_R, 2);
    const complex_t IT_0347 = mty::lt::C0iC(18, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0060, IT_0346, mty::lt::reg_int);
    const complex_t IT_0348 = m_b*IT_0347;
    const complex_t IT_0349 = mty::lt::C0iC(6, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0060, IT_0346, mty::lt::reg_int);
    const complex_t IT_0350 = m_b*IT_0349;
    const complex_t IT_0351 = mty::lt::C0iC(3, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0060, IT_0346, mty::lt::reg_int);
    const complex_t IT_0352 = m_b*IT_0351;
    const complex_t IT_0353 = mty::lt::C0iC(12, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0060, IT_0346, mty::lt::reg_int);
    const complex_t IT_0354 = m_b*IT_0353;
    const complex_t IT_0355 = IT_0348 + IT_0350 + IT_0352 + IT_0354;
    const complex_t IT_0356 = m_s*IT_0351;
    const complex_t IT_0357 = m_s*IT_0353;
    const complex_t IT_0358 = mty::lt::C0iC(15, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0060, IT_0346, mty::lt::reg_int);
    const complex_t IT_0359 = m_b*IT_0358;
    const complex_t IT_0360 = m_s*IT_0358;
    const complex_t IT_0361 = -IT_0356 + -IT_0357 + 2*IT_0359 + -IT_0360;
    const complex_t IT_0362 = IT_0355 + IT_0361;
    const complex_t IT_0363 = IT_0345*IT_0362;
    const complex_t IT_0364 = m_s*U_d1*conjq(V_cs)*e_em*IT_0031*U_su_15;
    const complex_t IT_0365 = IT_0045*IT_0364;
    const complex_t IT_0366 = 1.4142135623731*IT_0365;
    const complex_t IT_0367 = m_s*U_d1*conjq(V_ts)*e_em*IT_0031*U_su_25;
    const complex_t IT_0368 = IT_0045*IT_0367;
    const complex_t IT_0369 = 1.4142135623731*IT_0368;
    const complex_t IT_0370 = m_s*U_d1*V_us*e_em*IT_0031*U_su_05;
    const complex_t IT_0371 = IT_0045*IT_0370;
    const complex_t IT_0372 = 1.4142135623731*IT_0371;
    const complex_t IT_0373 = (complex_t{0, 1})*(IT_0366 + IT_0369 + IT_0372);
    const complex_t IT_0374 = 0.5*IT_0373;
    const complex_t IT_0375 = m_b*conjq(U_d1)*V_cb*e_em*IT_0031*conjq(U_su_15);
    const complex_t IT_0376 = IT_0045*IT_0375;
    const complex_t IT_0377 = 1.4142135623731*IT_0376;
    const complex_t IT_0378 = m_b*conjq(U_d1)*V_tb*e_em*IT_0031*conjq(U_su_25);
    const complex_t IT_0379 = IT_0045*IT_0378;
    const complex_t IT_0380 = 1.4142135623731*IT_0379;
    const complex_t IT_0381 = m_b*conjq(U_d1)*e_em*IT_0031*conjq(U_su_05)
      *V_ub_mod;
    const complex_t IT_0382 = IT_0188*IT_0381;
    const complex_t IT_0383 = 1.4142135623731*IT_0382;
    const complex_t IT_0384 = (complex_t{0, 1})*(IT_0377 + IT_0380 + IT_0383);
    const complex_t IT_0385 = 0.5*IT_0384;
    const complex_t IT_0386 = IT_0321*IT_0374*IT_0385;
    const complex_t IT_0387 = 0.101321183642338*IT_0386;
    const complex_t IT_0388 = powq(m_st_R, 2);
    const complex_t IT_0389 = mty::lt::C0iC(6, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0060, IT_0388, mty::lt::reg_int);
    const complex_t IT_0390 = m_b*IT_0389;
    const complex_t IT_0391 = mty::lt::C0iC(3, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0060, IT_0388, mty::lt::reg_int);
    const complex_t IT_0392 = m_b*IT_0391;
    const complex_t IT_0393 = mty::lt::C0iC(12, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0060, IT_0388, mty::lt::reg_int);
    const complex_t IT_0394 = m_b*IT_0393;
    const complex_t IT_0395 = mty::lt::C0iC(18, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0060, IT_0388, mty::lt::reg_int);
    const complex_t IT_0396 = m_b*IT_0395;
    const complex_t IT_0397 = IT_0390 + IT_0392 + IT_0394 + IT_0396;
    const complex_t IT_0398 = mty::lt::C0iC(15, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0060, IT_0388, mty::lt::reg_int);
    const complex_t IT_0399 = m_b*IT_0398;
    const complex_t IT_0400 = m_s*IT_0393;
    const complex_t IT_0401 = m_s*IT_0398;
    const complex_t IT_0402 = m_s*IT_0391;
    const complex_t IT_0403 = 2*IT_0399 + -IT_0400 + -IT_0401 + -IT_0402;
    const complex_t IT_0404 = IT_0397 + IT_0403;
    const complex_t IT_0405 = IT_0387*IT_0404;
    const complex_t IT_0406 = m_b*conjq(U_d1)*V_cb*e_em*IT_0031*conjq(U_su_13);
    const complex_t IT_0407 = IT_0045*IT_0406;
    const complex_t IT_0408 = 1.4142135623731*IT_0407;
    const complex_t IT_0409 = m_b*conjq(U_d1)*V_tb*e_em*IT_0031*conjq(U_su_23);
    const complex_t IT_0410 = IT_0045*IT_0409;
    const complex_t IT_0411 = 1.4142135623731*IT_0410;
    const complex_t IT_0412 = m_b*conjq(U_d1)*e_em*IT_0031*conjq(U_su_03)
      *V_ub_mod;
    const complex_t IT_0413 = IT_0188*IT_0412;
    const complex_t IT_0414 = 1.4142135623731*IT_0413;
    const complex_t IT_0415 = (complex_t{0, 1})*(IT_0408 + IT_0411 + IT_0414);
    const complex_t IT_0416 = 0.5*IT_0415;
    const complex_t IT_0417 = m_s*U_d1*conjq(V_cs)*e_em*IT_0031*U_su_13;
    const complex_t IT_0418 = IT_0045*IT_0417;
    const complex_t IT_0419 = 1.4142135623731*IT_0418;
    const complex_t IT_0420 = m_s*U_d1*conjq(V_ts)*e_em*IT_0031*U_su_23;
    const complex_t IT_0421 = IT_0045*IT_0420;
    const complex_t IT_0422 = 1.4142135623731*IT_0421;
    const complex_t IT_0423 = m_s*U_d1*V_us*e_em*IT_0031*U_su_03;
    const complex_t IT_0424 = IT_0045*IT_0423;
    const complex_t IT_0425 = 1.4142135623731*IT_0424;
    const complex_t IT_0426 = (complex_t{0, 1})*(IT_0419 + IT_0422 + IT_0425);
    const complex_t IT_0427 = 0.5*IT_0426;
    const complex_t IT_0428 = IT_0416*IT_0427;
    const complex_t IT_0429 = 0.101321183642338*IT_0428;
    const complex_t IT_0430 = powq(m_su_R, 2);
    const complex_t IT_0431 = mty::lt::C0iC(3, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0060, IT_0430, IT_0430, mty::lt::reg_int);
    const complex_t IT_0432 = IT_0059*IT_0431;
    const complex_t IT_0433 = m_b*IT_0432;
    const complex_t IT_0434 = mty::lt::C0iC(6, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0060, IT_0430, IT_0430, mty::lt::reg_int);
    const complex_t IT_0435 = IT_0059*IT_0434;
    const complex_t IT_0436 = m_s*IT_0435;
    const complex_t IT_0437 = mty::lt::C0iC(12, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0060, IT_0430, IT_0430, mty::lt::reg_int);
    const complex_t IT_0438 = IT_0064*IT_0437;
    const complex_t IT_0439 = m_b*IT_0438;
    const complex_t IT_0440 = mty::lt::C0iC(15, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0060, IT_0430, IT_0430, mty::lt::reg_int);
    const complex_t IT_0441 = IT_0064*IT_0440;
    const complex_t IT_0442 = m_s*IT_0441;
    const complex_t IT_0443 = IT_0433 + IT_0436 + IT_0439 + IT_0442;
    const complex_t IT_0444 = IT_0429*IT_0443;
    const complex_t IT_0445 = m_b*conjq(U_d2)*V_cb*e_em*IT_0031*conjq(U_su_13);
    const complex_t IT_0446 = IT_0045*IT_0445;
    const complex_t IT_0447 = 1.4142135623731*IT_0446;
    const complex_t IT_0448 = m_b*conjq(U_d2)*V_tb*e_em*IT_0031*conjq(U_su_23);
    const complex_t IT_0449 = IT_0045*IT_0448;
    const complex_t IT_0450 = 1.4142135623731*IT_0449;
    const complex_t IT_0451 = m_b*conjq(U_d2)*e_em*IT_0031*conjq(U_su_03)
      *V_ub_mod;
    const complex_t IT_0452 = IT_0188*IT_0451;
    const complex_t IT_0453 = 1.4142135623731*IT_0452;
    const complex_t IT_0454 = (complex_t{0, 1})*(IT_0447 + IT_0450 + IT_0453);
    const complex_t IT_0455 = 0.5*IT_0454;
    const complex_t IT_0456 = m_s*U_d2*conjq(V_cs)*e_em*IT_0031*U_su_13;
    const complex_t IT_0457 = IT_0045*IT_0456;
    const complex_t IT_0458 = 1.4142135623731*IT_0457;
    const complex_t IT_0459 = m_s*U_d2*conjq(V_ts)*e_em*IT_0031*U_su_23;
    const complex_t IT_0460 = IT_0045*IT_0459;
    const complex_t IT_0461 = 1.4142135623731*IT_0460;
    const complex_t IT_0462 = m_s*U_d2*V_us*e_em*IT_0031*U_su_03;
    const complex_t IT_0463 = IT_0045*IT_0462;
    const complex_t IT_0464 = 1.4142135623731*IT_0463;
    const complex_t IT_0465 = (complex_t{0, 1})*(IT_0458 + IT_0461 + IT_0464);
    const complex_t IT_0466 = 0.5*IT_0465;
    const complex_t IT_0467 = IT_0455*IT_0466;
    const complex_t IT_0468 = 0.101321183642338*IT_0467;
    const complex_t IT_0469 = mty::lt::C0iC(3, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0089, IT_0430, IT_0430, mty::lt::reg_int);
    const complex_t IT_0470 = IT_0059*IT_0469;
    const complex_t IT_0471 = m_b*IT_0470;
    const complex_t IT_0472 = mty::lt::C0iC(12, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0089, IT_0430, IT_0430, mty::lt::reg_int);
    const complex_t IT_0473 = IT_0064*IT_0472;
    const complex_t IT_0474 = m_b*IT_0473;
    const complex_t IT_0475 = mty::lt::C0iC(6, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0089, IT_0430, IT_0430, mty::lt::reg_int);
    const complex_t IT_0476 = IT_0059*IT_0475;
    const complex_t IT_0477 = m_s*IT_0476;
    const complex_t IT_0478 = mty::lt::C0iC(15, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0089, IT_0430, IT_0430, mty::lt::reg_int);
    const complex_t IT_0479 = IT_0064*IT_0478;
    const complex_t IT_0480 = m_s*IT_0479;
    const complex_t IT_0481 = IT_0471 + IT_0474 + IT_0477 + IT_0480;
    const complex_t IT_0482 = IT_0468*IT_0481;
    const complex_t IT_0483 = m_b*conjq(U_d2)*V_cb*e_em*IT_0031*conjq(U_su_11);
    const complex_t IT_0484 = IT_0045*IT_0483;
    const complex_t IT_0485 = 1.4142135623731*IT_0484;
    const complex_t IT_0486 = m_b*conjq(U_d2)*V_tb*e_em*IT_0031*conjq(U_su_21);
    const complex_t IT_0487 = IT_0045*IT_0486;
    const complex_t IT_0488 = 1.4142135623731*IT_0487;
    const complex_t IT_0489 = m_b*conjq(U_d2)*e_em*IT_0031*conjq(U_su_01)
      *V_ub_mod;
    const complex_t IT_0490 = IT_0188*IT_0489;
    const complex_t IT_0491 = 1.4142135623731*IT_0490;
    const complex_t IT_0492 = (complex_t{0, 1})*(IT_0485 + IT_0488 + IT_0491);
    const complex_t IT_0493 = 0.5*IT_0492;
    const complex_t IT_0494 = m_s*U_d2*conjq(V_cs)*e_em*IT_0031*U_su_11;
    const complex_t IT_0495 = IT_0045*IT_0494;
    const complex_t IT_0496 = 1.4142135623731*IT_0495;
    const complex_t IT_0497 = m_s*U_d2*conjq(V_ts)*e_em*IT_0031*U_su_21;
    const complex_t IT_0498 = IT_0045*IT_0497;
    const complex_t IT_0499 = 1.4142135623731*IT_0498;
    const complex_t IT_0500 = m_s*U_d2*V_us*e_em*IT_0031*U_su_01;
    const complex_t IT_0501 = IT_0045*IT_0500;
    const complex_t IT_0502 = 1.4142135623731*IT_0501;
    const complex_t IT_0503 = (complex_t{0, 1})*(IT_0496 + IT_0499 + IT_0502);
    const complex_t IT_0504 = 0.5*IT_0503;
    const complex_t IT_0505 = IT_0493*IT_0504;
    const complex_t IT_0506 = 0.101321183642338*IT_0505;
    const complex_t IT_0507 = powq(m_sc_L, 2);
    const complex_t IT_0508 = mty::lt::C0iC(3, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0089, IT_0507, IT_0507, mty::lt::reg_int);
    const complex_t IT_0509 = IT_0059*IT_0508;
    const complex_t IT_0510 = m_b*IT_0509;
    const complex_t IT_0511 = mty::lt::C0iC(12, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0089, IT_0507, IT_0507, mty::lt::reg_int);
    const complex_t IT_0512 = IT_0064*IT_0511;
    const complex_t IT_0513 = m_b*IT_0512;
    const complex_t IT_0514 = mty::lt::C0iC(6, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0089, IT_0507, IT_0507, mty::lt::reg_int);
    const complex_t IT_0515 = IT_0059*IT_0514;
    const complex_t IT_0516 = m_s*IT_0515;
    const complex_t IT_0517 = mty::lt::C0iC(15, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0089, IT_0507, IT_0507, mty::lt::reg_int);
    const complex_t IT_0518 = IT_0064*IT_0517;
    const complex_t IT_0519 = m_s*IT_0518;
    const complex_t IT_0520 = IT_0510 + IT_0513 + IT_0516 + IT_0519;
    const complex_t IT_0521 = IT_0506*IT_0520;
    const complex_t IT_0522 = V_cb*e_em*V_Wp1*conjq(U_su_11);
    const complex_t IT_0523 = IT_0019*IT_0522;
    const complex_t IT_0524 = V_tb*e_em*V_Wp1*conjq(U_su_21);
    const complex_t IT_0525 = IT_0019*IT_0524;
    const complex_t IT_0526 = e_em*V_Wp1*conjq(U_su_01)*V_ub_mod;
    const complex_t IT_0527 = IT_0025*IT_0526;
    const complex_t IT_0528 = m_c*V_cb*V_u1*e_em*IT_0031*conjq(U_su_41);
    const complex_t IT_0529 = IT_0030*IT_0528;
    const complex_t IT_0530 = 1.4142135623731*IT_0529;
    const complex_t IT_0531 = m_t*V_tb*V_u1*e_em*IT_0031*conjq(U_su_51);
    const complex_t IT_0532 = IT_0030*IT_0531;
    const complex_t IT_0533 = 1.4142135623731*IT_0532;
    const complex_t IT_0534 = m_u*V_u1*e_em*IT_0031*conjq(U_su_31)*V_ub_mod;
    const complex_t IT_0535 = IT_0038*IT_0534;
    const complex_t IT_0536 = 1.4142135623731*IT_0535;
    const complex_t IT_0537 = (complex_t{0, 1})*(IT_0523 + IT_0525 + IT_0527 +
       (-0.5)*IT_0530 + (-0.5)*IT_0533 + (-0.5)*IT_0536);
    const complex_t IT_0538 = m_s*U_d1*conjq(V_cs)*e_em*IT_0031*U_su_11;
    const complex_t IT_0539 = IT_0045*IT_0538;
    const complex_t IT_0540 = 1.4142135623731*IT_0539;
    const complex_t IT_0541 = m_s*U_d1*conjq(V_ts)*e_em*IT_0031*U_su_21;
    const complex_t IT_0542 = IT_0045*IT_0541;
    const complex_t IT_0543 = 1.4142135623731*IT_0542;
    const complex_t IT_0544 = m_s*U_d1*V_us*e_em*IT_0031*U_su_01;
    const complex_t IT_0545 = IT_0045*IT_0544;
    const complex_t IT_0546 = 1.4142135623731*IT_0545;
    const complex_t IT_0547 = (complex_t{0, 1})*(IT_0540 + IT_0543 + IT_0546);
    const complex_t IT_0548 = 0.5*IT_0547;
    const complex_t IT_0549 = IT_0537*IT_0548;
    const complex_t IT_0550 = IT_0017*IT_0549;
    const complex_t IT_0551 = mty::lt::C0iC(3, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0060, IT_0507, IT_0507, mty::lt::reg_int);
    const complex_t IT_0552 = IT_0064*IT_0551;
    const complex_t IT_0553 = mty::lt::C0iC(0, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0060, IT_0507, IT_0507, mty::lt::reg_int);
    const complex_t IT_0554 = IT_0059*IT_0553;
    const complex_t IT_0555 = IT_0552 + IT_0554;
    const complex_t IT_0556 = IT_0550*IT_0555;
    const complex_t IT_0557 = e_em*V_Wp2*conjq(V_Wp2);
    const complex_t IT_0558 = IT_0069*IT_0557;
    const complex_t IT_0559 = V_u2*conjq(V_u2)*e_em;
    const complex_t IT_0560 = IT_0069*IT_0559;
    const complex_t IT_0561 = (complex_t{0, 1})*(IT_0558 + IT_0560);
    const complex_t IT_0562 = m_b*conjq(U_d2)*V_cb*e_em*IT_0031*conjq(U_su_12);
    const complex_t IT_0563 = IT_0045*IT_0562;
    const complex_t IT_0564 = 1.4142135623731*IT_0563;
    const complex_t IT_0565 = m_b*conjq(U_d2)*V_tb*e_em*IT_0031*conjq(U_su_22);
    const complex_t IT_0566 = IT_0045*IT_0565;
    const complex_t IT_0567 = 1.4142135623731*IT_0566;
    const complex_t IT_0568 = m_b*conjq(U_d2)*e_em*IT_0031*conjq(U_su_02)
      *V_ub_mod;
    const complex_t IT_0569 = IT_0188*IT_0568;
    const complex_t IT_0570 = 1.4142135623731*IT_0569;
    const complex_t IT_0571 = (complex_t{0, 1})*(IT_0564 + IT_0567 + IT_0570);
    const complex_t IT_0572 = 0.5*IT_0571;
    const complex_t IT_0573 = IT_0157*IT_0561*IT_0572;
    const complex_t IT_0574 = IT_0174*IT_0573;
    const complex_t IT_0575 = mty::lt::C0iC(0, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0089, IT_0061, mty::lt::reg_int);
    const complex_t IT_0576 = mty::lt::C0iC(3, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0089, IT_0061, mty::lt::reg_int);
    const complex_t IT_0577 = mty::lt::C0iC(6, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0089, IT_0061, mty::lt::reg_int);
    const complex_t IT_0578 = IT_0575 + IT_0576 + IT_0577;
    const complex_t IT_0579 = IT_0574*IT_0578;
    const complex_t IT_0580 = e_em*conjq(V_Wp1)*V_Wp2;
    const complex_t IT_0581 = IT_0069*IT_0580;
    const complex_t IT_0582 = conjq(V_u1)*V_u2*e_em;
    const complex_t IT_0583 = IT_0069*IT_0582;
    const complex_t IT_0584 = (complex_t{0, 1})*(IT_0581 + IT_0583);
    const complex_t IT_0585 = m_b*conjq(U_d1)*V_cb*e_em*IT_0031*conjq(U_su_12);
    const complex_t IT_0586 = IT_0045*IT_0585;
    const complex_t IT_0587 = 1.4142135623731*IT_0586;
    const complex_t IT_0588 = m_b*conjq(U_d1)*V_tb*e_em*IT_0031*conjq(U_su_22);
    const complex_t IT_0589 = IT_0045*IT_0588;
    const complex_t IT_0590 = 1.4142135623731*IT_0589;
    const complex_t IT_0591 = m_b*conjq(U_d1)*e_em*IT_0031*conjq(U_su_02)
      *V_ub_mod;
    const complex_t IT_0592 = IT_0188*IT_0591;
    const complex_t IT_0593 = 1.4142135623731*IT_0592;
    const complex_t IT_0594 = (complex_t{0, 1})*(IT_0587 + IT_0590 + IT_0593);
    const complex_t IT_0595 = 0.5*IT_0594;
    const complex_t IT_0596 = IT_0157*IT_0584*IT_0595;
    const complex_t IT_0597 = IT_0017*IT_0596;
    const complex_t IT_0598 = IT_0093*IT_0597;
    const complex_t IT_0599 = V_cb*e_em*V_Wp1*conjq(U_su_13);
    const complex_t IT_0600 = IT_0019*IT_0599;
    const complex_t IT_0601 = V_tb*e_em*V_Wp1*conjq(U_su_23);
    const complex_t IT_0602 = IT_0019*IT_0601;
    const complex_t IT_0603 = e_em*V_Wp1*conjq(U_su_03)*V_ub_mod;
    const complex_t IT_0604 = IT_0025*IT_0603;
    const complex_t IT_0605 = m_c*V_cb*V_u1*e_em*IT_0031*conjq(U_su_43);
    const complex_t IT_0606 = IT_0030*IT_0605;
    const complex_t IT_0607 = 1.4142135623731*IT_0606;
    const complex_t IT_0608 = m_t*V_tb*V_u1*e_em*IT_0031*conjq(U_su_53);
    const complex_t IT_0609 = IT_0030*IT_0608;
    const complex_t IT_0610 = 1.4142135623731*IT_0609;
    const complex_t IT_0611 = m_u*V_u1*e_em*IT_0031*conjq(U_su_33)*V_ub_mod;
    const complex_t IT_0612 = IT_0038*IT_0611;
    const complex_t IT_0613 = 1.4142135623731*IT_0612;
    const complex_t IT_0614 = (complex_t{0, 1})*(IT_0600 + IT_0602 + IT_0604 +
       (-0.5)*IT_0607 + (-0.5)*IT_0610 + (-0.5)*IT_0613);
    const complex_t IT_0615 = V_us*e_em*conjq(V_Wp2)*U_su_03;
    const complex_t IT_0616 = IT_0019*IT_0615;
    const complex_t IT_0617 = conjq(V_cs)*e_em*conjq(V_Wp2)*U_su_13;
    const complex_t IT_0618 = IT_0019*IT_0617;
    const complex_t IT_0619 = conjq(V_ts)*e_em*conjq(V_Wp2)*U_su_23;
    const complex_t IT_0620 = IT_0019*IT_0619;
    const complex_t IT_0621 = m_u*conjq(V_u2)*V_us*e_em*IT_0031*U_su_33;
    const complex_t IT_0622 = IT_0030*IT_0621;
    const complex_t IT_0623 = 1.4142135623731*IT_0622;
    const complex_t IT_0624 = m_c*conjq(V_cs)*conjq(V_u2)*e_em*IT_0031*U_su_43;
    const complex_t IT_0625 = IT_0030*IT_0624;
    const complex_t IT_0626 = 1.4142135623731*IT_0625;
    const complex_t IT_0627 = m_t*conjq(V_ts)*conjq(V_u2)*e_em*IT_0031*U_su_53;
    const complex_t IT_0628 = IT_0030*IT_0627;
    const complex_t IT_0629 = 1.4142135623731*IT_0628;
    const complex_t IT_0630 = (complex_t{0, 1})*(IT_0616 + IT_0618 + IT_0620 +
       (-0.5)*IT_0623 + (-0.5)*IT_0626 + (-0.5)*IT_0629);
    const complex_t IT_0631 = IT_0584*IT_0614*IT_0630;
    const complex_t IT_0632 = 0.101321183642338*IT_0631;
    const complex_t IT_0633 = mty::lt::C0iC(18, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0089, IT_0430, mty::lt::reg_int);
    const complex_t IT_0634 = m_b*IT_0633;
    const complex_t IT_0635 = mty::lt::C0iC(3, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0089, IT_0430, mty::lt::reg_int);
    const complex_t IT_0636 = m_b*IT_0635;
    const complex_t IT_0637 = mty::lt::C0iC(12, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0089, IT_0430, mty::lt::reg_int);
    const complex_t IT_0638 = m_b*IT_0637;
    const complex_t IT_0639 = mty::lt::C0iC(6, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0089, IT_0430, mty::lt::reg_int);
    const complex_t IT_0640 = m_b*IT_0639;
    const complex_t IT_0641 = IT_0634 + IT_0636 + IT_0638 + IT_0640;
    const complex_t IT_0642 = mty::lt::C0iC(15, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0089, IT_0430, mty::lt::reg_int);
    const complex_t IT_0643 = m_s*IT_0642;
    const complex_t IT_0644 = m_s*IT_0635;
    const complex_t IT_0645 = m_s*IT_0637;
    const complex_t IT_0646 = m_b*IT_0642;
    const complex_t IT_0647 = -IT_0643 + -IT_0644 + -IT_0645 + 2*IT_0646;
    const complex_t IT_0648 = IT_0641 + IT_0647;
    const complex_t IT_0649 = IT_0632*IT_0648;
    const complex_t IT_0650 = IT_0042*IT_0157*IT_0584;
    const complex_t IT_0651 = 0.101321183642338*IT_0650;
    const complex_t IT_0652 = m_b*IT_0091;
    const complex_t IT_0653 = m_b*IT_0090;
    const complex_t IT_0654 = mty::lt::C0iC(12, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0089, IT_0061, mty::lt::reg_int);
    const complex_t IT_0655 = m_b*IT_0654;
    const complex_t IT_0656 = mty::lt::C0iC(18, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0089, IT_0061, mty::lt::reg_int);
    const complex_t IT_0657 = m_b*IT_0656;
    const complex_t IT_0658 = IT_0652 + IT_0653 + IT_0655 + IT_0657;
    const complex_t IT_0659 = m_s*IT_0090;
    const complex_t IT_0660 = mty::lt::C0iC(15, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0089, IT_0061, mty::lt::reg_int);
    const complex_t IT_0661 = m_b*IT_0660;
    const complex_t IT_0662 = m_s*IT_0654;
    const complex_t IT_0663 = m_s*IT_0660;
    const complex_t IT_0664 = -IT_0659 + 2*IT_0661 + -IT_0662 + -IT_0663;
    const complex_t IT_0665 = IT_0658 + IT_0664;
    const complex_t IT_0666 = IT_0651*IT_0665;
    const complex_t IT_0667 = e_em*V_Wp1*conjq(V_Wp2);
    const complex_t IT_0668 = IT_0069*IT_0667;
    const complex_t IT_0669 = V_u1*conjq(V_u2)*e_em;
    const complex_t IT_0670 = IT_0069*IT_0669;
    const complex_t IT_0671 = (complex_t{0, 1})*(IT_0668 + IT_0670);
    const complex_t IT_0672 = IT_0110*IT_0572*IT_0671;
    const complex_t IT_0673 = IT_0174*IT_0672;
    const complex_t IT_0674 = mty::lt::C0iC(0, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0060, IT_0061, mty::lt::reg_int);
    const complex_t IT_0675 = mty::lt::C0iC(3, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0060, IT_0061, mty::lt::reg_int);
    const complex_t IT_0676 = mty::lt::C0iC(6, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0060, IT_0061, mty::lt::reg_int);
    const complex_t IT_0677 = IT_0674 + IT_0675 + IT_0676;
    const complex_t IT_0678 = IT_0673*IT_0677;
    const complex_t IT_0679 = V_us*e_em*conjq(V_Wp1)*U_su_03;
    const complex_t IT_0680 = IT_0019*IT_0679;
    const complex_t IT_0681 = conjq(V_cs)*e_em*conjq(V_Wp1)*U_su_13;
    const complex_t IT_0682 = IT_0019*IT_0681;
    const complex_t IT_0683 = conjq(V_ts)*e_em*conjq(V_Wp1)*U_su_23;
    const complex_t IT_0684 = IT_0019*IT_0683;
    const complex_t IT_0685 = m_u*conjq(V_u1)*V_us*e_em*IT_0031*U_su_33;
    const complex_t IT_0686 = IT_0030*IT_0685;
    const complex_t IT_0687 = 1.4142135623731*IT_0686;
    const complex_t IT_0688 = m_c*conjq(V_cs)*conjq(V_u1)*e_em*IT_0031*U_su_43;
    const complex_t IT_0689 = IT_0030*IT_0688;
    const complex_t IT_0690 = 1.4142135623731*IT_0689;
    const complex_t IT_0691 = m_t*conjq(V_ts)*conjq(V_u1)*e_em*IT_0031*U_su_53;
    const complex_t IT_0692 = IT_0030*IT_0691;
    const complex_t IT_0693 = 1.4142135623731*IT_0692;
    const complex_t IT_0694 = (complex_t{0, 1})*(IT_0680 + IT_0682 + IT_0684 +
       (-0.5)*IT_0687 + (-0.5)*IT_0690 + (-0.5)*IT_0693);
    const complex_t IT_0695 = IT_0455*IT_0671*IT_0694;
    const complex_t IT_0696 = IT_0174*IT_0695;
    const complex_t IT_0697 = mty::lt::C0iC(0, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0060, IT_0430, mty::lt::reg_int);
    const complex_t IT_0698 = mty::lt::C0iC(3, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0060, IT_0430, mty::lt::reg_int);
    const complex_t IT_0699 = mty::lt::C0iC(6, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0060, IT_0430, mty::lt::reg_int);
    const complex_t IT_0700 = IT_0697 + IT_0698 + IT_0699;
    const complex_t IT_0701 = IT_0696*IT_0700;
    const complex_t IT_0702 = m_b*conjq(U_d2)*V_cb*e_em*IT_0031*conjq(U_su_15);
    const complex_t IT_0703 = IT_0045*IT_0702;
    const complex_t IT_0704 = 1.4142135623731*IT_0703;
    const complex_t IT_0705 = m_b*conjq(U_d2)*V_tb*e_em*IT_0031*conjq(U_su_25);
    const complex_t IT_0706 = IT_0045*IT_0705;
    const complex_t IT_0707 = 1.4142135623731*IT_0706;
    const complex_t IT_0708 = m_b*conjq(U_d2)*e_em*IT_0031*conjq(U_su_05)
      *V_ub_mod;
    const complex_t IT_0709 = IT_0188*IT_0708;
    const complex_t IT_0710 = 1.4142135623731*IT_0709;
    const complex_t IT_0711 = (complex_t{0, 1})*(IT_0704 + IT_0707 + IT_0710);
    const complex_t IT_0712 = 0.5*IT_0711;
    const complex_t IT_0713 = V_us*e_em*conjq(V_Wp1)*U_su_05;
    const complex_t IT_0714 = IT_0019*IT_0713;
    const complex_t IT_0715 = conjq(V_cs)*e_em*conjq(V_Wp1)*U_su_15;
    const complex_t IT_0716 = IT_0019*IT_0715;
    const complex_t IT_0717 = conjq(V_ts)*e_em*conjq(V_Wp1)*U_su_25;
    const complex_t IT_0718 = IT_0019*IT_0717;
    const complex_t IT_0719 = m_u*conjq(V_u1)*V_us*e_em*IT_0031*U_su_35;
    const complex_t IT_0720 = IT_0030*IT_0719;
    const complex_t IT_0721 = 1.4142135623731*IT_0720;
    const complex_t IT_0722 = m_c*conjq(V_cs)*conjq(V_u1)*e_em*IT_0031*U_su_45;
    const complex_t IT_0723 = IT_0030*IT_0722;
    const complex_t IT_0724 = 1.4142135623731*IT_0723;
    const complex_t IT_0725 = m_t*conjq(V_ts)*conjq(V_u1)*e_em*IT_0031*U_su_55;
    const complex_t IT_0726 = IT_0030*IT_0725;
    const complex_t IT_0727 = 1.4142135623731*IT_0726;
    const complex_t IT_0728 = (complex_t{0, 1})*(IT_0714 + IT_0716 + IT_0718 +
       (-0.5)*IT_0721 + (-0.5)*IT_0724 + (-0.5)*IT_0727);
    const complex_t IT_0729 = IT_0671*IT_0712*IT_0728;
    const complex_t IT_0730 = IT_0174*IT_0729;
    const complex_t IT_0731 = mty::lt::C0iC(0, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0060, IT_0388, mty::lt::reg_int);
    const complex_t IT_0732 = mty::lt::C0iC(3, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0060, IT_0388, mty::lt::reg_int);
    const complex_t IT_0733 = mty::lt::C0iC(6, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0060, IT_0388, mty::lt::reg_int);
    const complex_t IT_0734 = IT_0731 + IT_0732 + IT_0733;
    const complex_t IT_0735 = IT_0730*IT_0734;
    const complex_t IT_0736 = e_em*U_Wm2*conjq(U_Wm2);
    const complex_t IT_0737 = IT_0069*IT_0736;
    const complex_t IT_0738 = U_d2*conjq(U_d2)*e_em;
    const complex_t IT_0739 = IT_0069*IT_0738;
    const complex_t IT_0740 = (complex_t{0, 1})*(IT_0737 + IT_0739);
    const complex_t IT_0741 = -IT_0740;
    const complex_t IT_0742 = IT_0493*IT_0504*IT_0741;
    const complex_t IT_0743 = 0.101321183642338*IT_0742;
    const complex_t IT_0744 = mty::lt::C0iC(6, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0089, IT_0507, mty::lt::reg_int);
    const complex_t IT_0745 = m_b*IT_0744;
    const complex_t IT_0746 = mty::lt::C0iC(3, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0089, IT_0507, mty::lt::reg_int);
    const complex_t IT_0747 = m_b*IT_0746;
    const complex_t IT_0748 = mty::lt::C0iC(12, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0089, IT_0507, mty::lt::reg_int);
    const complex_t IT_0749 = m_b*IT_0748;
    const complex_t IT_0750 = mty::lt::C0iC(18, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0089, IT_0507, mty::lt::reg_int);
    const complex_t IT_0751 = m_b*IT_0750;
    const complex_t IT_0752 = IT_0745 + IT_0747 + IT_0749 + IT_0751;
    const complex_t IT_0753 = m_s*IT_0746;
    const complex_t IT_0754 = m_s*IT_0748;
    const complex_t IT_0755 = mty::lt::C0iC(15, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0089, IT_0507, mty::lt::reg_int);
    const complex_t IT_0756 = m_s*IT_0755;
    const complex_t IT_0757 = m_b*IT_0755;
    const complex_t IT_0758 = -IT_0753 + -IT_0754 + -IT_0756 + 2*IT_0757;
    const complex_t IT_0759 = IT_0752 + IT_0758;
    const complex_t IT_0760 = IT_0743*IT_0759;
    const complex_t IT_0761 = IT_0086*IT_0572*IT_0741;
    const complex_t IT_0762 = 0.101321183642338*IT_0761;
    const complex_t IT_0763 = m_b*IT_0576;
    const complex_t IT_0764 = m_b*IT_0577;
    const complex_t IT_0765 = mty::lt::C0iC(12, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0089, IT_0061, mty::lt::reg_int);
    const complex_t IT_0766 = m_b*IT_0765;
    const complex_t IT_0767 = mty::lt::C0iC(18, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0089, IT_0061, mty::lt::reg_int);
    const complex_t IT_0768 = m_b*IT_0767;
    const complex_t IT_0769 = IT_0763 + IT_0764 + IT_0766 + IT_0768;
    const complex_t IT_0770 = m_s*IT_0576;
    const complex_t IT_0771 = m_s*IT_0765;
    const complex_t IT_0772 = mty::lt::C0iC(15, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0089, IT_0061, mty::lt::reg_int);
    const complex_t IT_0773 = m_s*IT_0772;
    const complex_t IT_0774 = m_b*IT_0772;
    const complex_t IT_0775 = -IT_0770 + -IT_0771 + -IT_0773 + 2*IT_0774;
    const complex_t IT_0776 = IT_0769 + IT_0775;
    const complex_t IT_0777 = IT_0762*IT_0776;
    const complex_t IT_0778 = IT_0455*IT_0466*IT_0741;
    const complex_t IT_0779 = 0.101321183642338*IT_0778;
    const complex_t IT_0780 = mty::lt::C0iC(6, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0089, IT_0430, mty::lt::reg_int);
    const complex_t IT_0781 = m_b*IT_0780;
    const complex_t IT_0782 = mty::lt::C0iC(3, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0089, IT_0430, mty::lt::reg_int);
    const complex_t IT_0783 = m_b*IT_0782;
    const complex_t IT_0784 = mty::lt::C0iC(12, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0089, IT_0430, mty::lt::reg_int);
    const complex_t IT_0785 = m_b*IT_0784;
    const complex_t IT_0786 = mty::lt::C0iC(18, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0089, IT_0430, mty::lt::reg_int);
    const complex_t IT_0787 = m_b*IT_0786;
    const complex_t IT_0788 = IT_0781 + IT_0783 + IT_0785 + IT_0787;
    const complex_t IT_0789 = m_s*IT_0782;
    const complex_t IT_0790 = m_s*IT_0784;
    const complex_t IT_0791 = mty::lt::C0iC(15, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0089, IT_0430, mty::lt::reg_int);
    const complex_t IT_0792 = m_s*IT_0791;
    const complex_t IT_0793 = m_b*IT_0791;
    const complex_t IT_0794 = -IT_0789 + -IT_0790 + -IT_0792 + 2*IT_0793;
    const complex_t IT_0795 = IT_0788 + IT_0794;
    const complex_t IT_0796 = IT_0779*IT_0795;
    const complex_t IT_0797 = m_s*U_d2*conjq(V_cs)*e_em*IT_0031*U_su_10;
    const complex_t IT_0798 = IT_0045*IT_0797;
    const complex_t IT_0799 = 1.4142135623731*IT_0798;
    const complex_t IT_0800 = m_s*U_d2*conjq(V_ts)*e_em*IT_0031*U_su_20;
    const complex_t IT_0801 = IT_0045*IT_0800;
    const complex_t IT_0802 = 1.4142135623731*IT_0801;
    const complex_t IT_0803 = m_s*U_d2*V_us*e_em*IT_0031*U_su_00;
    const complex_t IT_0804 = IT_0045*IT_0803;
    const complex_t IT_0805 = 1.4142135623731*IT_0804;
    const complex_t IT_0806 = (complex_t{0, 1})*(IT_0799 + IT_0802 + IT_0805);
    const complex_t IT_0807 = 0.5*IT_0806;
    const complex_t IT_0808 = IT_0229*IT_0741*IT_0807;
    const complex_t IT_0809 = 0.101321183642338*IT_0808;
    const complex_t IT_0810 = mty::lt::C0iC(6, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0089, IT_0212, mty::lt::reg_int);
    const complex_t IT_0811 = m_b*IT_0810;
    const complex_t IT_0812 = mty::lt::C0iC(3, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0089, IT_0212, mty::lt::reg_int);
    const complex_t IT_0813 = m_b*IT_0812;
    const complex_t IT_0814 = mty::lt::C0iC(12, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0089, IT_0212, mty::lt::reg_int);
    const complex_t IT_0815 = m_b*IT_0814;
    const complex_t IT_0816 = mty::lt::C0iC(18, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0089, IT_0212, mty::lt::reg_int);
    const complex_t IT_0817 = m_b*IT_0816;
    const complex_t IT_0818 = IT_0811 + IT_0813 + IT_0815 + IT_0817;
    const complex_t IT_0819 = m_s*IT_0812;
    const complex_t IT_0820 = m_s*IT_0814;
    const complex_t IT_0821 = mty::lt::C0iC(15, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0089, IT_0212, mty::lt::reg_int);
    const complex_t IT_0822 = m_s*IT_0821;
    const complex_t IT_0823 = m_b*IT_0821;
    const complex_t IT_0824 = -IT_0819 + -IT_0820 + -IT_0822 + 2*IT_0823;
    const complex_t IT_0825 = IT_0818 + IT_0824;
    const complex_t IT_0826 = IT_0809*IT_0825;
    const complex_t IT_0827 = m_b*conjq(U_d2)*V_cb*e_em*IT_0031*conjq(U_su_14);
    const complex_t IT_0828 = IT_0045*IT_0827;
    const complex_t IT_0829 = 1.4142135623731*IT_0828;
    const complex_t IT_0830 = m_b*conjq(U_d2)*V_tb*e_em*IT_0031*conjq(U_su_24);
    const complex_t IT_0831 = IT_0045*IT_0830;
    const complex_t IT_0832 = 1.4142135623731*IT_0831;
    const complex_t IT_0833 = m_b*conjq(U_d2)*e_em*IT_0031*conjq(U_su_04)
      *V_ub_mod;
    const complex_t IT_0834 = IT_0188*IT_0833;
    const complex_t IT_0835 = 1.4142135623731*IT_0834;
    const complex_t IT_0836 = (complex_t{0, 1})*(IT_0829 + IT_0832 + IT_0835);
    const complex_t IT_0837 = 0.5*IT_0836;
    const complex_t IT_0838 = m_s*U_d2*conjq(V_cs)*e_em*IT_0031*U_su_14;
    const complex_t IT_0839 = IT_0045*IT_0838;
    const complex_t IT_0840 = 1.4142135623731*IT_0839;
    const complex_t IT_0841 = m_s*U_d2*conjq(V_ts)*e_em*IT_0031*U_su_24;
    const complex_t IT_0842 = IT_0045*IT_0841;
    const complex_t IT_0843 = 1.4142135623731*IT_0842;
    const complex_t IT_0844 = m_s*U_d2*V_us*e_em*IT_0031*U_su_04;
    const complex_t IT_0845 = IT_0045*IT_0844;
    const complex_t IT_0846 = 1.4142135623731*IT_0845;
    const complex_t IT_0847 = (complex_t{0, 1})*(IT_0840 + IT_0843 + IT_0846);
    const complex_t IT_0848 = 0.5*IT_0847;
    const complex_t IT_0849 = IT_0741*IT_0837*IT_0848;
    const complex_t IT_0850 = 0.101321183642338*IT_0849;
    const complex_t IT_0851 = mty::lt::C0iC(6, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0089, IT_0346, mty::lt::reg_int);
    const complex_t IT_0852 = m_b*IT_0851;
    const complex_t IT_0853 = mty::lt::C0iC(3, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0089, IT_0346, mty::lt::reg_int);
    const complex_t IT_0854 = m_b*IT_0853;
    const complex_t IT_0855 = mty::lt::C0iC(12, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0089, IT_0346, mty::lt::reg_int);
    const complex_t IT_0856 = m_b*IT_0855;
    const complex_t IT_0857 = mty::lt::C0iC(18, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0089, IT_0346, mty::lt::reg_int);
    const complex_t IT_0858 = m_b*IT_0857;
    const complex_t IT_0859 = IT_0852 + IT_0854 + IT_0856 + IT_0858;
    const complex_t IT_0860 = m_s*IT_0853;
    const complex_t IT_0861 = m_s*IT_0855;
    const complex_t IT_0862 = mty::lt::C0iC(15, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0089, IT_0346, mty::lt::reg_int);
    const complex_t IT_0863 = m_s*IT_0862;
    const complex_t IT_0864 = m_b*IT_0862;
    const complex_t IT_0865 = -IT_0860 + -IT_0861 + -IT_0863 + 2*IT_0864;
    const complex_t IT_0866 = IT_0859 + IT_0865;
    const complex_t IT_0867 = IT_0850*IT_0866;
    const complex_t IT_0868 = m_s*U_d2*conjq(V_cs)*e_em*IT_0031*U_su_15;
    const complex_t IT_0869 = IT_0045*IT_0868;
    const complex_t IT_0870 = 1.4142135623731*IT_0869;
    const complex_t IT_0871 = m_s*U_d2*conjq(V_ts)*e_em*IT_0031*U_su_25;
    const complex_t IT_0872 = IT_0045*IT_0871;
    const complex_t IT_0873 = 1.4142135623731*IT_0872;
    const complex_t IT_0874 = m_s*U_d2*V_us*e_em*IT_0031*U_su_05;
    const complex_t IT_0875 = IT_0045*IT_0874;
    const complex_t IT_0876 = 1.4142135623731*IT_0875;
    const complex_t IT_0877 = (complex_t{0, 1})*(IT_0870 + IT_0873 + IT_0876);
    const complex_t IT_0878 = 0.5*IT_0877;
    const complex_t IT_0879 = IT_0712*IT_0741*IT_0878;
    const complex_t IT_0880 = 0.101321183642338*IT_0879;
    const complex_t IT_0881 = mty::lt::C0iC(18, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0089, IT_0388, mty::lt::reg_int);
    const complex_t IT_0882 = m_b*IT_0881;
    const complex_t IT_0883 = mty::lt::C0iC(6, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0089, IT_0388, mty::lt::reg_int);
    const complex_t IT_0884 = m_b*IT_0883;
    const complex_t IT_0885 = mty::lt::C0iC(3, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0089, IT_0388, mty::lt::reg_int);
    const complex_t IT_0886 = m_b*IT_0885;
    const complex_t IT_0887 = mty::lt::C0iC(12, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0089, IT_0388, mty::lt::reg_int);
    const complex_t IT_0888 = m_b*IT_0887;
    const complex_t IT_0889 = IT_0882 + IT_0884 + IT_0886 + IT_0888;
    const complex_t IT_0890 = mty::lt::C0iC(15, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0089, IT_0388, mty::lt::reg_int);
    const complex_t IT_0891 = m_b*IT_0890;
    const complex_t IT_0892 = m_s*IT_0885;
    const complex_t IT_0893 = m_s*IT_0887;
    const complex_t IT_0894 = m_s*IT_0890;
    const complex_t IT_0895 = 2*IT_0891 + -IT_0892 + -IT_0893 + -IT_0894;
    const complex_t IT_0896 = IT_0889 + IT_0895;
    const complex_t IT_0897 = IT_0880*IT_0896;
    const complex_t IT_0898 = e_em*V_Wp1*conjq(V_Wp1);
    const complex_t IT_0899 = IT_0069*IT_0898;
    const complex_t IT_0900 = V_u1*conjq(V_u1)*e_em;
    const complex_t IT_0901 = IT_0069*IT_0900;
    const complex_t IT_0902 = (complex_t{0, 1})*(IT_0899 + IT_0901);
    const complex_t IT_0903 = V_us*e_em*conjq(V_Wp1)*U_su_01;
    const complex_t IT_0904 = IT_0019*IT_0903;
    const complex_t IT_0905 = conjq(V_cs)*e_em*conjq(V_Wp1)*U_su_11;
    const complex_t IT_0906 = IT_0019*IT_0905;
    const complex_t IT_0907 = conjq(V_ts)*e_em*conjq(V_Wp1)*U_su_21;
    const complex_t IT_0908 = IT_0019*IT_0907;
    const complex_t IT_0909 = m_u*conjq(V_u1)*V_us*e_em*IT_0031*U_su_31;
    const complex_t IT_0910 = IT_0030*IT_0909;
    const complex_t IT_0911 = 1.4142135623731*IT_0910;
    const complex_t IT_0912 = m_c*conjq(V_cs)*conjq(V_u1)*e_em*IT_0031*U_su_41;
    const complex_t IT_0913 = IT_0030*IT_0912;
    const complex_t IT_0914 = 1.4142135623731*IT_0913;
    const complex_t IT_0915 = m_t*conjq(V_ts)*conjq(V_u1)*e_em*IT_0031*U_su_51;
    const complex_t IT_0916 = IT_0030*IT_0915;
    const complex_t IT_0917 = 1.4142135623731*IT_0916;
    const complex_t IT_0918 = (complex_t{0, 1})*(IT_0904 + IT_0906 + IT_0908 +
       (-0.5)*IT_0911 + (-0.5)*IT_0914 + (-0.5)*IT_0917);
    const complex_t IT_0919 = IT_0537*IT_0902*IT_0918;
    const complex_t IT_0920 = 0.101321183642338*IT_0919;
    const complex_t IT_0921 = mty::lt::C0iC(6, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0060, IT_0507, mty::lt::reg_int);
    const complex_t IT_0922 = m_b*IT_0921;
    const complex_t IT_0923 = mty::lt::C0iC(3, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0060, IT_0507, mty::lt::reg_int);
    const complex_t IT_0924 = m_b*IT_0923;
    const complex_t IT_0925 = mty::lt::C0iC(12, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0060, IT_0507, mty::lt::reg_int);
    const complex_t IT_0926 = m_b*IT_0925;
    const complex_t IT_0927 = mty::lt::C0iC(18, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0060, IT_0507, mty::lt::reg_int);
    const complex_t IT_0928 = m_b*IT_0927;
    const complex_t IT_0929 = IT_0922 + IT_0924 + IT_0926 + IT_0928;
    const complex_t IT_0930 = m_s*IT_0923;
    const complex_t IT_0931 = mty::lt::C0iC(15, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0060, IT_0507, mty::lt::reg_int);
    const complex_t IT_0932 = m_b*IT_0931;
    const complex_t IT_0933 = m_s*IT_0925;
    const complex_t IT_0934 = m_s*IT_0931;
    const complex_t IT_0935 = -IT_0930 + 2*IT_0932 + -IT_0933 + -IT_0934;
    const complex_t IT_0936 = IT_0929 + IT_0935;
    const complex_t IT_0937 = IT_0920*IT_0936;
    const complex_t IT_0938 = IT_0209*IT_0269*IT_0902;
    const complex_t IT_0939 = 0.101321183642338*IT_0938;
    const complex_t IT_0940 = mty::lt::C0iC(6, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0060, IT_0212, mty::lt::reg_int);
    const complex_t IT_0941 = m_b*IT_0940;
    const complex_t IT_0942 = mty::lt::C0iC(3, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0060, IT_0212, mty::lt::reg_int);
    const complex_t IT_0943 = m_b*IT_0942;
    const complex_t IT_0944 = mty::lt::C0iC(12, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0060, IT_0212, mty::lt::reg_int);
    const complex_t IT_0945 = m_b*IT_0944;
    const complex_t IT_0946 = mty::lt::C0iC(18, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0060, IT_0212, mty::lt::reg_int);
    const complex_t IT_0947 = m_b*IT_0946;
    const complex_t IT_0948 = IT_0941 + IT_0943 + IT_0945 + IT_0947;
    const complex_t IT_0949 = m_s*IT_0942;
    const complex_t IT_0950 = mty::lt::C0iC(15, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0060, IT_0212, mty::lt::reg_int);
    const complex_t IT_0951 = m_b*IT_0950;
    const complex_t IT_0952 = m_s*IT_0944;
    const complex_t IT_0953 = m_s*IT_0950;
    const complex_t IT_0954 = -IT_0949 + 2*IT_0951 + -IT_0952 + -IT_0953;
    const complex_t IT_0955 = IT_0948 + IT_0954;
    const complex_t IT_0956 = IT_0939*IT_0955;
    const complex_t IT_0957 = IT_0042*IT_0110*IT_0902;
    const complex_t IT_0958 = 0.101321183642338*IT_0957;
    const complex_t IT_0959 = mty::lt::C0iC(3, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0060, IT_0061, mty::lt::reg_int);
    const complex_t IT_0960 = m_s*IT_0959;
    const complex_t IT_0961 = mty::lt::C0iC(15, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0060, IT_0061, mty::lt::reg_int);
    const complex_t IT_0962 = m_b*IT_0961;
    const complex_t IT_0963 = mty::lt::C0iC(12, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0060, IT_0061, mty::lt::reg_int);
    const complex_t IT_0964 = m_s*IT_0963;
    const complex_t IT_0965 = m_s*IT_0961;
    const complex_t IT_0966 = -IT_0960 + 2*IT_0962 + -IT_0964 + -IT_0965;
    const complex_t IT_0967 = m_b*IT_0963;
    const complex_t IT_0968 = mty::lt::C0iC(18, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0060, IT_0061, mty::lt::reg_int);
    const complex_t IT_0969 = m_b*IT_0968;
    const complex_t IT_0970 = mty::lt::C0iC(6, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0060, IT_0061, mty::lt::reg_int);
    const complex_t IT_0971 = m_b*IT_0970;
    const complex_t IT_0972 = m_b*IT_0959;
    const complex_t IT_0973 = IT_0967 + IT_0969 + IT_0971 + IT_0972;
    const complex_t IT_0974 = IT_0966 + IT_0973;
    const complex_t IT_0975 = IT_0958*IT_0974;
    const complex_t IT_0976 = V_cb*e_em*V_Wp1*conjq(U_su_14);
    const complex_t IT_0977 = IT_0019*IT_0976;
    const complex_t IT_0978 = V_tb*e_em*V_Wp1*conjq(U_su_24);
    const complex_t IT_0979 = IT_0019*IT_0978;
    const complex_t IT_0980 = e_em*V_Wp1*conjq(U_su_04)*V_ub_mod;
    const complex_t IT_0981 = IT_0025*IT_0980;
    const complex_t IT_0982 = m_c*V_cb*V_u1*e_em*IT_0031*conjq(U_su_44);
    const complex_t IT_0983 = IT_0030*IT_0982;
    const complex_t IT_0984 = 1.4142135623731*IT_0983;
    const complex_t IT_0985 = m_t*V_tb*V_u1*e_em*IT_0031*conjq(U_su_54);
    const complex_t IT_0986 = IT_0030*IT_0985;
    const complex_t IT_0987 = 1.4142135623731*IT_0986;
    const complex_t IT_0988 = m_u*V_u1*e_em*IT_0031*conjq(U_su_34)*V_ub_mod;
    const complex_t IT_0989 = IT_0038*IT_0988;
    const complex_t IT_0990 = 1.4142135623731*IT_0989;
    const complex_t IT_0991 = (complex_t{0, 1})*(IT_0977 + IT_0979 + IT_0981 +
       (-0.5)*IT_0984 + (-0.5)*IT_0987 + (-0.5)*IT_0990);
    const complex_t IT_0992 = V_us*e_em*conjq(V_Wp1)*U_su_04;
    const complex_t IT_0993 = IT_0019*IT_0992;
    const complex_t IT_0994 = conjq(V_cs)*e_em*conjq(V_Wp1)*U_su_14;
    const complex_t IT_0995 = IT_0019*IT_0994;
    const complex_t IT_0996 = conjq(V_ts)*e_em*conjq(V_Wp1)*U_su_24;
    const complex_t IT_0997 = IT_0019*IT_0996;
    const complex_t IT_0998 = m_u*conjq(V_u1)*V_us*e_em*IT_0031*U_su_34;
    const complex_t IT_0999 = IT_0030*IT_0998;
    const complex_t IT_1000 = 1.4142135623731*IT_0999;
    const complex_t IT_1001 = m_c*conjq(V_cs)*conjq(V_u1)*e_em*IT_0031*U_su_44;
    const complex_t IT_1002 = IT_0030*IT_1001;
    const complex_t IT_1003 = 1.4142135623731*IT_1002;
    const complex_t IT_1004 = m_t*conjq(V_ts)*conjq(V_u1)*e_em*IT_0031*U_su_54;
    const complex_t IT_1005 = IT_0030*IT_1004;
    const complex_t IT_1006 = 1.4142135623731*IT_1005;
    const complex_t IT_1007 = (complex_t{0, 1})*(IT_0993 + IT_0995 + IT_0997 +
       (-0.5)*IT_1000 + (-0.5)*IT_1003 + (-0.5)*IT_1006);
    const complex_t IT_1008 = IT_0991*IT_1007;
    const complex_t IT_1009 = 0.101321183642338*IT_1008;
    const complex_t IT_1010 = mty::lt::C0iC(6, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0060, IT_0346, IT_0346, mty::lt::reg_int);
    const complex_t IT_1011 = IT_0059*IT_1010;
    const complex_t IT_1012 = m_s*IT_1011;
    const complex_t IT_1013 = mty::lt::C0iC(3, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0060, IT_0346, IT_0346, mty::lt::reg_int);
    const complex_t IT_1014 = IT_0059*IT_1013;
    const complex_t IT_1015 = m_b*IT_1014;
    const complex_t IT_1016 = mty::lt::C0iC(12, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0060, IT_0346, IT_0346, mty::lt::reg_int);
    const complex_t IT_1017 = IT_0064*IT_1016;
    const complex_t IT_1018 = m_b*IT_1017;
    const complex_t IT_1019 = mty::lt::C0iC(15, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0060, IT_0346, IT_0346, mty::lt::reg_int);
    const complex_t IT_1020 = IT_0064*IT_1019;
    const complex_t IT_1021 = m_s*IT_1020;
    const complex_t IT_1022 = IT_1012 + IT_1015 + IT_1018 + IT_1021;
    const complex_t IT_1023 = IT_1009*IT_1022;
    const complex_t IT_1024 = V_cb*e_em*V_Wp2*conjq(U_su_14);
    const complex_t IT_1025 = IT_0019*IT_1024;
    const complex_t IT_1026 = V_tb*e_em*V_Wp2*conjq(U_su_24);
    const complex_t IT_1027 = IT_0019*IT_1026;
    const complex_t IT_1028 = e_em*V_Wp2*conjq(U_su_04)*V_ub_mod;
    const complex_t IT_1029 = IT_0025*IT_1028;
    const complex_t IT_1030 = m_c*V_cb*V_u2*e_em*IT_0031*conjq(U_su_44);
    const complex_t IT_1031 = IT_0030*IT_1030;
    const complex_t IT_1032 = 1.4142135623731*IT_1031;
    const complex_t IT_1033 = m_t*V_tb*V_u2*e_em*IT_0031*conjq(U_su_54);
    const complex_t IT_1034 = IT_0030*IT_1033;
    const complex_t IT_1035 = 1.4142135623731*IT_1034;
    const complex_t IT_1036 = m_u*V_u2*e_em*IT_0031*conjq(U_su_34)*V_ub_mod;
    const complex_t IT_1037 = IT_0038*IT_1036;
    const complex_t IT_1038 = 1.4142135623731*IT_1037;
    const complex_t IT_1039 = (complex_t{0, 1})*(IT_1025 + IT_1027 + IT_1029 +
       (-0.5)*IT_1032 + (-0.5)*IT_1035 + (-0.5)*IT_1038);
    const complex_t IT_1040 = V_us*e_em*conjq(V_Wp2)*U_su_04;
    const complex_t IT_1041 = IT_0019*IT_1040;
    const complex_t IT_1042 = conjq(V_cs)*e_em*conjq(V_Wp2)*U_su_14;
    const complex_t IT_1043 = IT_0019*IT_1042;
    const complex_t IT_1044 = conjq(V_ts)*e_em*conjq(V_Wp2)*U_su_24;
    const complex_t IT_1045 = IT_0019*IT_1044;
    const complex_t IT_1046 = m_u*conjq(V_u2)*V_us*e_em*IT_0031*U_su_34;
    const complex_t IT_1047 = IT_0030*IT_1046;
    const complex_t IT_1048 = 1.4142135623731*IT_1047;
    const complex_t IT_1049 = m_c*conjq(V_cs)*conjq(V_u2)*e_em*IT_0031*U_su_44;
    const complex_t IT_1050 = IT_0030*IT_1049;
    const complex_t IT_1051 = 1.4142135623731*IT_1050;
    const complex_t IT_1052 = m_t*conjq(V_ts)*conjq(V_u2)*e_em*IT_0031*U_su_54;
    const complex_t IT_1053 = IT_0030*IT_1052;
    const complex_t IT_1054 = 1.4142135623731*IT_1053;
    const complex_t IT_1055 = (complex_t{0, 1})*(IT_1041 + IT_1043 + IT_1045 +
       (-0.5)*IT_1048 + (-0.5)*IT_1051 + (-0.5)*IT_1054);
    const complex_t IT_1056 = IT_1039*IT_1055;
    const complex_t IT_1057 = 0.101321183642338*IT_1056;
    const complex_t IT_1058 = mty::lt::C0iC(6, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0089, IT_0346, IT_0346, mty::lt::reg_int);
    const complex_t IT_1059 = IT_0059*IT_1058;
    const complex_t IT_1060 = m_s*IT_1059;
    const complex_t IT_1061 = mty::lt::C0iC(3, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0089, IT_0346, IT_0346, mty::lt::reg_int);
    const complex_t IT_1062 = IT_0059*IT_1061;
    const complex_t IT_1063 = m_b*IT_1062;
    const complex_t IT_1064 = mty::lt::C0iC(12, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0089, IT_0346, IT_0346, mty::lt::reg_int);
    const complex_t IT_1065 = IT_0064*IT_1064;
    const complex_t IT_1066 = m_b*IT_1065;
    const complex_t IT_1067 = mty::lt::C0iC(15, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0089, IT_0346, IT_0346, mty::lt::reg_int);
    const complex_t IT_1068 = IT_0064*IT_1067;
    const complex_t IT_1069 = m_s*IT_1068;
    const complex_t IT_1070 = IT_1060 + IT_1063 + IT_1066 + IT_1069;
    const complex_t IT_1071 = IT_1057*IT_1070;
    const complex_t IT_1072 = IT_0374*IT_0385;
    const complex_t IT_1073 = 0.101321183642338*IT_1072;
    const complex_t IT_1074 = mty::lt::C0iC(3, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0060, IT_0388, IT_0388, mty::lt::reg_int);
    const complex_t IT_1075 = IT_0059*IT_1074;
    const complex_t IT_1076 = m_b*IT_1075;
    const complex_t IT_1077 = mty::lt::C0iC(12, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0060, IT_0388, IT_0388, mty::lt::reg_int);
    const complex_t IT_1078 = IT_0064*IT_1077;
    const complex_t IT_1079 = m_b*IT_1078;
    const complex_t IT_1080 = mty::lt::C0iC(6, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0060, IT_0388, IT_0388, mty::lt::reg_int);
    const complex_t IT_1081 = IT_0059*IT_1080;
    const complex_t IT_1082 = m_s*IT_1081;
    const complex_t IT_1083 = mty::lt::C0iC(15, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0060, IT_0388, IT_0388, mty::lt::reg_int);
    const complex_t IT_1084 = IT_0064*IT_1083;
    const complex_t IT_1085 = m_s*IT_1084;
    const complex_t IT_1086 = IT_1076 + IT_1079 + IT_1082 + IT_1085;
    const complex_t IT_1087 = IT_1073*IT_1086;
    const complex_t IT_1088 = e_em*conjq(U_Wm1)*U_Wm2;
    const complex_t IT_1089 = IT_0069*IT_1088;
    const complex_t IT_1090 = conjq(U_d1)*U_d2*e_em;
    const complex_t IT_1091 = IT_0069*IT_1090;
    const complex_t IT_1092 = (complex_t{0, 1})*(IT_1089 + IT_1091);
    const complex_t IT_1093 = -IT_1092;
    const complex_t IT_1094 = m_s*U_d1*conjq(V_cs)*e_em*IT_0031*U_su_10;
    const complex_t IT_1095 = IT_0045*IT_1094;
    const complex_t IT_1096 = 1.4142135623731*IT_1095;
    const complex_t IT_1097 = m_s*U_d1*conjq(V_ts)*e_em*IT_0031*U_su_20;
    const complex_t IT_1098 = IT_0045*IT_1097;
    const complex_t IT_1099 = 1.4142135623731*IT_1098;
    const complex_t IT_1100 = m_s*U_d1*V_us*e_em*IT_0031*U_su_00;
    const complex_t IT_1101 = IT_0045*IT_1100;
    const complex_t IT_1102 = 1.4142135623731*IT_1101;
    const complex_t IT_1103 = (complex_t{0, 1})*(IT_1096 + IT_1099 + IT_1102);
    const complex_t IT_1104 = 0.5*IT_1103;
    const complex_t IT_1105 = IT_0300*IT_1093*IT_1104;
    const complex_t IT_1106 = IT_0174*IT_1105;
    const complex_t IT_1107 = mty::lt::C0iC(0, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0060, IT_0212, mty::lt::reg_int);
    const complex_t IT_1108 = mty::lt::C0iC(3, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0060, IT_0212, mty::lt::reg_int);
    const complex_t IT_1109 = mty::lt::C0iC(6, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0060, IT_0212, mty::lt::reg_int);
    const complex_t IT_1110 = IT_1107 + IT_1108 + IT_1109;
    const complex_t IT_1111 = IT_1106*IT_1110;
    const complex_t IT_1112 = V_cb*e_em*V_Wp2*conjq(U_su_13);
    const complex_t IT_1113 = IT_0019*IT_1112;
    const complex_t IT_1114 = V_tb*e_em*V_Wp2*conjq(U_su_23);
    const complex_t IT_1115 = IT_0019*IT_1114;
    const complex_t IT_1116 = e_em*V_Wp2*conjq(U_su_03)*V_ub_mod;
    const complex_t IT_1117 = IT_0025*IT_1116;
    const complex_t IT_1118 = m_c*V_cb*V_u2*e_em*IT_0031*conjq(U_su_43);
    const complex_t IT_1119 = IT_0030*IT_1118;
    const complex_t IT_1120 = 1.4142135623731*IT_1119;
    const complex_t IT_1121 = m_t*V_tb*V_u2*e_em*IT_0031*conjq(U_su_53);
    const complex_t IT_1122 = IT_0030*IT_1121;
    const complex_t IT_1123 = 1.4142135623731*IT_1122;
    const complex_t IT_1124 = m_u*V_u2*e_em*IT_0031*conjq(U_su_33)*V_ub_mod;
    const complex_t IT_1125 = IT_0038*IT_1124;
    const complex_t IT_1126 = 1.4142135623731*IT_1125;
    const complex_t IT_1127 = (complex_t{0, 1})*(IT_1113 + IT_1115 + IT_1117 +
       (-0.5)*IT_1120 + (-0.5)*IT_1123 + (-0.5)*IT_1126);
    const complex_t IT_1128 = IT_0427*IT_1093*IT_1127;
    const complex_t IT_1129 = IT_0174*IT_1128;
    const complex_t IT_1130 = IT_0700*IT_1129;
    const complex_t IT_1131 = IT_0332*IT_1039*IT_1093;
    const complex_t IT_1132 = IT_0174*IT_1131;
    const complex_t IT_1133 = mty::lt::C0iC(0, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0060, IT_0346, mty::lt::reg_int);
    const complex_t IT_1134 = mty::lt::C0iC(3, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0060, IT_0346, mty::lt::reg_int);
    const complex_t IT_1135 = mty::lt::C0iC(6, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0060, IT_0346, mty::lt::reg_int);
    const complex_t IT_1136 = IT_1133 + IT_1134 + IT_1135;
    const complex_t IT_1137 = IT_1132*IT_1136;
    const complex_t IT_1138 = IT_0042*IT_0056*IT_0321;
    const complex_t IT_1139 = IT_0017*IT_1138;
    const complex_t IT_1140 = mty::lt::C0iC(0, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0060, IT_0061, mty::lt::reg_int);
    const complex_t IT_1141 = IT_0959 + IT_0970 + IT_1140;
    const complex_t IT_1142 = IT_1139*IT_1141;
    const complex_t IT_1143 = IT_0321*IT_0537*IT_0548;
    const complex_t IT_1144 = IT_0017*IT_1143;
    const complex_t IT_1145 = mty::lt::C0iC(0, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0060, IT_0507, mty::lt::reg_int);
    const complex_t IT_1146 = IT_0921 + IT_0923 + IT_1145;
    const complex_t IT_1147 = IT_1144*IT_1146;
    const complex_t IT_1148 = IT_0321*IT_0427*IT_0614;
    const complex_t IT_1149 = IT_0017*IT_1148;
    const complex_t IT_1150 = mty::lt::C0iC(6, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0060, IT_0430, mty::lt::reg_int);
    const complex_t IT_1151 = mty::lt::C0iC(3, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0060, IT_0430, mty::lt::reg_int);
    const complex_t IT_1152 = mty::lt::C0iC(0, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0060, IT_0430, mty::lt::reg_int);
    const complex_t IT_1153 = IT_1150 + IT_1151 + IT_1152;
    const complex_t IT_1154 = IT_1149*IT_1153;
    const complex_t IT_1155 = IT_0269*IT_0321*IT_1104;
    const complex_t IT_1156 = IT_0017*IT_1155;
    const complex_t IT_1157 = mty::lt::C0iC(0, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0060, IT_0212, mty::lt::reg_int);
    const complex_t IT_1158 = IT_0940 + IT_0942 + IT_1157;
    const complex_t IT_1159 = IT_1156*IT_1158;
    const complex_t IT_1160 = IT_0321*IT_0332*IT_0991;
    const complex_t IT_1161 = IT_0017*IT_1160;
    const complex_t IT_1162 = mty::lt::C0iC(0, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0060, IT_0346, mty::lt::reg_int);
    const complex_t IT_1163 = IT_0349 + IT_0351 + IT_1162;
    const complex_t IT_1164 = IT_1161*IT_1163;
    const complex_t IT_1165 = V_cb*e_em*V_Wp1*conjq(U_su_15);
    const complex_t IT_1166 = IT_0019*IT_1165;
    const complex_t IT_1167 = V_tb*e_em*V_Wp1*conjq(U_su_25);
    const complex_t IT_1168 = IT_0019*IT_1167;
    const complex_t IT_1169 = e_em*V_Wp1*conjq(U_su_05)*V_ub_mod;
    const complex_t IT_1170 = IT_0025*IT_1169;
    const complex_t IT_1171 = m_c*V_cb*V_u1*e_em*IT_0031*conjq(U_su_45);
    const complex_t IT_1172 = IT_0030*IT_1171;
    const complex_t IT_1173 = 1.4142135623731*IT_1172;
    const complex_t IT_1174 = m_t*V_tb*V_u1*e_em*IT_0031*conjq(U_su_55);
    const complex_t IT_1175 = IT_0030*IT_1174;
    const complex_t IT_1176 = 1.4142135623731*IT_1175;
    const complex_t IT_1177 = m_u*V_u1*e_em*IT_0031*conjq(U_su_35)*V_ub_mod;
    const complex_t IT_1178 = IT_0038*IT_1177;
    const complex_t IT_1179 = 1.4142135623731*IT_1178;
    const complex_t IT_1180 = (complex_t{0, 1})*(IT_1166 + IT_1168 + IT_1170 +
       (-0.5)*IT_1173 + (-0.5)*IT_1176 + (-0.5)*IT_1179);
    const complex_t IT_1181 = IT_0321*IT_0374*IT_1180;
    const complex_t IT_1182 = IT_0017*IT_1181;
    const complex_t IT_1183 = mty::lt::C0iC(0, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0060, IT_0388, mty::lt::reg_int);
    const complex_t IT_1184 = IT_0389 + IT_0391 + IT_1183;
    const complex_t IT_1185 = IT_1182*IT_1184;
    const complex_t IT_1186 = IT_0455*IT_0630;
    const complex_t IT_1187 = IT_0174*IT_1186;
    const complex_t IT_1188 = mty::lt::C0iC(0, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0089, IT_0430, IT_0430, mty::lt::reg_int);
    const complex_t IT_1189 = IT_0059*IT_1188;
    const complex_t IT_1190 = IT_0064*IT_0469;
    const complex_t IT_1191 = IT_1189 + IT_1190;
    const complex_t IT_1192 = IT_1187*IT_1191;
    const complex_t IT_1193 = V_cb*e_em*V_Wp2*conjq(U_su_11);
    const complex_t IT_1194 = IT_0019*IT_1193;
    const complex_t IT_1195 = V_tb*e_em*V_Wp2*conjq(U_su_21);
    const complex_t IT_1196 = IT_0019*IT_1195;
    const complex_t IT_1197 = e_em*V_Wp2*conjq(U_su_01)*V_ub_mod;
    const complex_t IT_1198 = IT_0025*IT_1197;
    const complex_t IT_1199 = m_c*V_cb*V_u2*e_em*IT_0031*conjq(U_su_41);
    const complex_t IT_1200 = IT_0030*IT_1199;
    const complex_t IT_1201 = 1.4142135623731*IT_1200;
    const complex_t IT_1202 = m_t*V_tb*V_u2*e_em*IT_0031*conjq(U_su_51);
    const complex_t IT_1203 = IT_0030*IT_1202;
    const complex_t IT_1204 = 1.4142135623731*IT_1203;
    const complex_t IT_1205 = m_u*V_u2*e_em*IT_0031*conjq(U_su_31)*V_ub_mod;
    const complex_t IT_1206 = IT_0038*IT_1205;
    const complex_t IT_1207 = 1.4142135623731*IT_1206;
    const complex_t IT_1208 = (complex_t{0, 1})*(IT_1194 + IT_1196 + IT_1198 +
       (-0.5)*IT_1201 + (-0.5)*IT_1204 + (-0.5)*IT_1207);
    const complex_t IT_1209 = IT_0504*IT_1208;
    const complex_t IT_1210 = IT_0174*IT_1209;
    const complex_t IT_1211 = mty::lt::C0iC(0, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0089, IT_0507, IT_0507, mty::lt::reg_int);
    const complex_t IT_1212 = IT_0059*IT_1211;
    const complex_t IT_1213 = IT_0064*IT_0508;
    const complex_t IT_1214 = IT_1212 + IT_1213;
    const complex_t IT_1215 = IT_1210*IT_1214;
    const complex_t IT_1216 = IT_0193*IT_0245*IT_0584;
    const complex_t IT_1217 = IT_0017*IT_1216;
    const complex_t IT_1218 = mty::lt::C0iC(6, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0089, IT_0212, mty::lt::reg_int);
    const complex_t IT_1219 = mty::lt::C0iC(3, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0089, IT_0212, mty::lt::reg_int);
    const complex_t IT_1220 = mty::lt::C0iC(0, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0089, IT_0212, mty::lt::reg_int);
    const complex_t IT_1221 = IT_1218 + IT_1219 + IT_1220;
    const complex_t IT_1222 = IT_1217*IT_1221;
    const complex_t IT_1223 = IT_0343*IT_0584*IT_1055;
    const complex_t IT_1224 = IT_0017*IT_1223;
    const complex_t IT_1225 = mty::lt::C0iC(3, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0089, IT_0346, mty::lt::reg_int);
    const complex_t IT_1226 = mty::lt::C0iC(6, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0089, IT_0346, mty::lt::reg_int);
    const complex_t IT_1227 = mty::lt::C0iC(0, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0089, IT_0346, mty::lt::reg_int);
    const complex_t IT_1228 = IT_1225 + IT_1226 + IT_1227;
    const complex_t IT_1229 = IT_1224*IT_1228;
    const complex_t IT_1230 = V_us*e_em*conjq(V_Wp2)*U_su_05;
    const complex_t IT_1231 = IT_0019*IT_1230;
    const complex_t IT_1232 = conjq(V_cs)*e_em*conjq(V_Wp2)*U_su_15;
    const complex_t IT_1233 = IT_0019*IT_1232;
    const complex_t IT_1234 = conjq(V_ts)*e_em*conjq(V_Wp2)*U_su_25;
    const complex_t IT_1235 = IT_0019*IT_1234;
    const complex_t IT_1236 = m_u*conjq(V_u2)*V_us*e_em*IT_0031*U_su_35;
    const complex_t IT_1237 = IT_0030*IT_1236;
    const complex_t IT_1238 = 1.4142135623731*IT_1237;
    const complex_t IT_1239 = m_c*conjq(V_cs)*conjq(V_u2)*e_em*IT_0031*U_su_45;
    const complex_t IT_1240 = IT_0030*IT_1239;
    const complex_t IT_1241 = 1.4142135623731*IT_1240;
    const complex_t IT_1242 = m_t*conjq(V_ts)*conjq(V_u2)*e_em*IT_0031*U_su_55;
    const complex_t IT_1243 = IT_0030*IT_1242;
    const complex_t IT_1244 = 1.4142135623731*IT_1243;
    const complex_t IT_1245 = (complex_t{0, 1})*(IT_1231 + IT_1233 + IT_1235 +
       (-0.5)*IT_1238 + (-0.5)*IT_1241 + (-0.5)*IT_1244);
    const complex_t IT_1246 = IT_0385*IT_0584*IT_1245;
    const complex_t IT_1247 = IT_0017*IT_1246;
    const complex_t IT_1248 = mty::lt::C0iC(3, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0089, IT_0388, mty::lt::reg_int);
    const complex_t IT_1249 = mty::lt::C0iC(6, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0089, IT_0388, mty::lt::reg_int);
    const complex_t IT_1250 = mty::lt::C0iC(0, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0089, IT_0388, mty::lt::reg_int);
    const complex_t IT_1251 = IT_1248 + IT_1249 + IT_1250;
    const complex_t IT_1252 = IT_1247*IT_1251;
    const complex_t IT_1253 = IT_0504*IT_0741*IT_1208;
    const complex_t IT_1254 = IT_0174*IT_1253;
    const complex_t IT_1255 = mty::lt::C0iC(0, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0089, IT_0507, mty::lt::reg_int);
    const complex_t IT_1256 = IT_0744 + IT_0746 + IT_1255;
    const complex_t IT_1257 = IT_1254*IT_1256;
    const complex_t IT_1258 = IT_0466*IT_0741*IT_1127;
    const complex_t IT_1259 = IT_0174*IT_1258;
    const complex_t IT_1260 = mty::lt::C0iC(0, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0089, IT_0430, mty::lt::reg_int);
    const complex_t IT_1261 = IT_0780 + IT_0782 + IT_1260;
    const complex_t IT_1262 = IT_1259*IT_1261;
    const complex_t IT_1263 = IT_0300*IT_0741*IT_0807;
    const complex_t IT_1264 = IT_0174*IT_1263;
    const complex_t IT_1265 = mty::lt::C0iC(0, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0089, IT_0212, mty::lt::reg_int);
    const complex_t IT_1266 = IT_0810 + IT_0812 + IT_1265;
    const complex_t IT_1267 = IT_1264*IT_1266;
    const complex_t IT_1268 = IT_0741*IT_0848*IT_1039;
    const complex_t IT_1269 = IT_0174*IT_1268;
    const complex_t IT_1270 = mty::lt::C0iC(0, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0089, IT_0346, mty::lt::reg_int);
    const complex_t IT_1271 = IT_0851 + IT_0853 + IT_1270;
    const complex_t IT_1272 = IT_1269*IT_1271;
    const complex_t IT_1273 = V_cb*e_em*V_Wp2*conjq(U_su_15);
    const complex_t IT_1274 = IT_0019*IT_1273;
    const complex_t IT_1275 = V_tb*e_em*V_Wp2*conjq(U_su_25);
    const complex_t IT_1276 = IT_0019*IT_1275;
    const complex_t IT_1277 = e_em*V_Wp2*conjq(U_su_05)*V_ub_mod;
    const complex_t IT_1278 = IT_0025*IT_1277;
    const complex_t IT_1279 = m_c*V_cb*V_u2*e_em*IT_0031*conjq(U_su_45);
    const complex_t IT_1280 = IT_0030*IT_1279;
    const complex_t IT_1281 = 1.4142135623731*IT_1280;
    const complex_t IT_1282 = m_t*V_tb*V_u2*e_em*IT_0031*conjq(U_su_55);
    const complex_t IT_1283 = IT_0030*IT_1282;
    const complex_t IT_1284 = 1.4142135623731*IT_1283;
    const complex_t IT_1285 = m_u*V_u2*e_em*IT_0031*conjq(U_su_35)*V_ub_mod;
    const complex_t IT_1286 = IT_0038*IT_1285;
    const complex_t IT_1287 = 1.4142135623731*IT_1286;
    const complex_t IT_1288 = (complex_t{0, 1})*(IT_1274 + IT_1276 + IT_1278 +
       (-0.5)*IT_1281 + (-0.5)*IT_1284 + (-0.5)*IT_1287);
    const complex_t IT_1289 = IT_0741*IT_0878*IT_1288;
    const complex_t IT_1290 = IT_0174*IT_1289;
    const complex_t IT_1291 = mty::lt::C0iC(0, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0089, IT_0388, mty::lt::reg_int);
    const complex_t IT_1292 = IT_0883 + IT_0885 + IT_1291;
    const complex_t IT_1293 = IT_1290*IT_1292;
    const complex_t IT_1294 = IT_0548*IT_1093*IT_1208;
    const complex_t IT_1295 = IT_0174*IT_1294;
    const complex_t IT_1296 = mty::lt::C0iC(0, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0060, IT_0507, mty::lt::reg_int);
    const complex_t IT_1297 = mty::lt::C0iC(3, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0060, IT_0507, mty::lt::reg_int);
    const complex_t IT_1298 = mty::lt::C0iC(6, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0060, IT_0507, mty::lt::reg_int);
    const complex_t IT_1299 = IT_1296 + IT_1297 + IT_1298;
    const complex_t IT_1300 = IT_1295*IT_1299;
    const complex_t IT_1301 = IT_0374*IT_1093*IT_1288;
    const complex_t IT_1302 = IT_0174*IT_1301;
    const complex_t IT_1303 = IT_0734*IT_1302;
    const complex_t IT_1304 = IT_0343*IT_0902*IT_1007;
    const complex_t IT_1305 = IT_0017*IT_1304;
    const complex_t IT_1306 = IT_1163*IT_1305;
    const complex_t IT_1307 = IT_0343*IT_1007;
    const complex_t IT_1308 = IT_0017*IT_1307;
    const complex_t IT_1309 = IT_0064*IT_1013;
    const complex_t IT_1310 = mty::lt::C0iC(0, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0060, IT_0346, IT_0346, mty::lt::reg_int);
    const complex_t IT_1311 = IT_0059*IT_1310;
    const complex_t IT_1312 = IT_1309 + IT_1311;
    const complex_t IT_1313 = IT_1308*IT_1312;
    const complex_t IT_1314 = mty::lt::C0iC(0, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0089, IT_0346, IT_0346, mty::lt::reg_int);
    const complex_t IT_1315 = IT_0059*IT_1314;
    const complex_t IT_1316 = IT_0064*IT_1061;
    const complex_t IT_1317 = IT_1315 + IT_1316;
    const complex_t IT_1318 = IT_0837*IT_1055;
    const complex_t IT_1319 = IT_0174*IT_1318;
    const complex_t IT_1320 = IT_1317*IT_1319;
    const complex_t IT_1321 = IT_0385*IT_0728;
    const complex_t IT_1322 = IT_0017*IT_1321;
    const complex_t IT_1323 = mty::lt::C0iC(0, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0060, IT_0388, IT_0388, mty::lt::reg_int);
    const complex_t IT_1324 = IT_0059*IT_1323;
    const complex_t IT_1325 = IT_0064*IT_1074;
    const complex_t IT_1326 = IT_1324 + IT_1325;
    const complex_t IT_1327 = IT_1322*IT_1326;
    const complex_t IT_1328 = IT_0878*IT_1288;
    const complex_t IT_1329 = IT_0174*IT_1328;
    const complex_t IT_1330 = mty::lt::C0iC(0, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0089, IT_0388, IT_0388, mty::lt::reg_int);
    const complex_t IT_1331 = IT_0059*IT_1330;
    const complex_t IT_1332 = mty::lt::C0iC(3, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0089, IT_0388, IT_0388, mty::lt::reg_int);
    const complex_t IT_1333 = IT_0064*IT_1332;
    const complex_t IT_1334 = IT_1331 + IT_1333;
    const complex_t IT_1335 = IT_1329*IT_1334;
    const complex_t IT_1336 = m_b*conjq(U_d1)*V_cb*e_em*IT_0031*conjq(U_su_11);
    const complex_t IT_1337 = IT_0045*IT_1336;
    const complex_t IT_1338 = 1.4142135623731*IT_1337;
    const complex_t IT_1339 = m_b*conjq(U_d1)*V_tb*e_em*IT_0031*conjq(U_su_21);
    const complex_t IT_1340 = IT_0045*IT_1339;
    const complex_t IT_1341 = 1.4142135623731*IT_1340;
    const complex_t IT_1342 = m_b*conjq(U_d1)*e_em*IT_0031*conjq(U_su_01)
      *V_ub_mod;
    const complex_t IT_1343 = IT_0188*IT_1342;
    const complex_t IT_1344 = 1.4142135623731*IT_1343;
    const complex_t IT_1345 = (complex_t{0, 1})*(IT_1338 + IT_1341 + IT_1344);
    const complex_t IT_1346 = 0.5*IT_1345;
    const complex_t IT_1347 = IT_0075*IT_0504*IT_1346;
    const complex_t IT_1348 = 0.101321183642338*IT_1347;
    const complex_t IT_1349 = mty::lt::C0iC(3, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0089, IT_0507, mty::lt::reg_int);
    const complex_t IT_1350 = m_b*IT_1349;
    const complex_t IT_1351 = mty::lt::C0iC(6, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0089, IT_0507, mty::lt::reg_int);
    const complex_t IT_1352 = m_b*IT_1351;
    const complex_t IT_1353 = mty::lt::C0iC(12, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0089, IT_0507, mty::lt::reg_int);
    const complex_t IT_1354 = m_b*IT_1353;
    const complex_t IT_1355 = mty::lt::C0iC(18, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0089, IT_0507, mty::lt::reg_int);
    const complex_t IT_1356 = m_b*IT_1355;
    const complex_t IT_1357 = IT_1350 + IT_1352 + IT_1354 + IT_1356;
    const complex_t IT_1358 = m_s*IT_1349;
    const complex_t IT_1359 = mty::lt::C0iC(15, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0089, IT_0507, mty::lt::reg_int);
    const complex_t IT_1360 = m_b*IT_1359;
    const complex_t IT_1361 = m_s*IT_1353;
    const complex_t IT_1362 = m_s*IT_1359;
    const complex_t IT_1363 = -IT_1358 + 2*IT_1360 + -IT_1361 + -IT_1362;
    const complex_t IT_1364 = IT_1357 + IT_1363;
    const complex_t IT_1365 = IT_1348*IT_1364;
    const complex_t IT_1366 = IT_0548*IT_1346;
    const complex_t IT_1367 = 0.101321183642338*IT_1366;
    const complex_t IT_1368 = IT_0059*IT_0551;
    const complex_t IT_1369 = m_b*IT_1368;
    const complex_t IT_1370 = mty::lt::C0iC(12, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0060, IT_0507, IT_0507, mty::lt::reg_int);
    const complex_t IT_1371 = IT_0064*IT_1370;
    const complex_t IT_1372 = m_b*IT_1371;
    const complex_t IT_1373 = mty::lt::C0iC(6, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0060, IT_0507, IT_0507, mty::lt::reg_int);
    const complex_t IT_1374 = IT_0059*IT_1373;
    const complex_t IT_1375 = m_s*IT_1374;
    const complex_t IT_1376 = mty::lt::C0iC(15, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0060, IT_0507, IT_0507, mty::lt::reg_int);
    const complex_t IT_1377 = IT_0064*IT_1376;
    const complex_t IT_1378 = m_s*IT_1377;
    const complex_t IT_1379 = IT_1369 + IT_1372 + IT_1375 + IT_1378;
    const complex_t IT_1380 = IT_1367*IT_1379;
    const complex_t IT_1381 = V_us*e_em*conjq(V_Wp2)*U_su_01;
    const complex_t IT_1382 = IT_0019*IT_1381;
    const complex_t IT_1383 = conjq(V_cs)*e_em*conjq(V_Wp2)*U_su_11;
    const complex_t IT_1384 = IT_0019*IT_1383;
    const complex_t IT_1385 = conjq(V_ts)*e_em*conjq(V_Wp2)*U_su_21;
    const complex_t IT_1386 = IT_0019*IT_1385;
    const complex_t IT_1387 = m_u*conjq(V_u2)*V_us*e_em*IT_0031*U_su_31;
    const complex_t IT_1388 = IT_0030*IT_1387;
    const complex_t IT_1389 = 1.4142135623731*IT_1388;
    const complex_t IT_1390 = m_c*conjq(V_cs)*conjq(V_u2)*e_em*IT_0031*U_su_41;
    const complex_t IT_1391 = IT_0030*IT_1390;
    const complex_t IT_1392 = 1.4142135623731*IT_1391;
    const complex_t IT_1393 = m_t*conjq(V_ts)*conjq(V_u2)*e_em*IT_0031*U_su_51;
    const complex_t IT_1394 = IT_0030*IT_1393;
    const complex_t IT_1395 = 1.4142135623731*IT_1394;
    const complex_t IT_1396 = (complex_t{0, 1})*(IT_1382 + IT_1384 + IT_1386 +
       (-0.5)*IT_1389 + (-0.5)*IT_1392 + (-0.5)*IT_1395);
    const complex_t IT_1397 = IT_0537*IT_0584*IT_1396;
    const complex_t IT_1398 = 0.101321183642338*IT_1397;
    const complex_t IT_1399 = IT_1364*IT_1398;
    const complex_t IT_1400 = IT_0209*IT_0300*IT_0671;
    const complex_t IT_1401 = 0.101321183642338*IT_1400;
    const complex_t IT_1402 = m_b*IT_1108;
    const complex_t IT_1403 = m_b*IT_1109;
    const complex_t IT_1404 = mty::lt::C0iC(12, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0060, IT_0212, mty::lt::reg_int);
    const complex_t IT_1405 = m_b*IT_1404;
    const complex_t IT_1406 = mty::lt::C0iC(18, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0060, IT_0212, mty::lt::reg_int);
    const complex_t IT_1407 = m_b*IT_1406;
    const complex_t IT_1408 = IT_1402 + IT_1403 + IT_1405 + IT_1407;
    const complex_t IT_1409 = m_s*IT_1108;
    const complex_t IT_1410 = m_s*IT_1404;
    const complex_t IT_1411 = mty::lt::C0iC(15, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0060, IT_0212, mty::lt::reg_int);
    const complex_t IT_1412 = m_s*IT_1411;
    const complex_t IT_1413 = m_b*IT_1411;
    const complex_t IT_1414 = -IT_1409 + -IT_1410 + -IT_1412 + 2*IT_1413;
    const complex_t IT_1415 = IT_1408 + IT_1414;
    const complex_t IT_1416 = IT_1401*IT_1415;
    const complex_t IT_1417 = IT_0269*IT_1104;
    const complex_t IT_1418 = IT_0017*IT_1417;
    const complex_t IT_1419 = IT_0217*IT_1418;
    const complex_t IT_1420 = IT_0671*IT_0918*IT_1208;
    const complex_t IT_1421 = 0.101321183642338*IT_1420;
    const complex_t IT_1422 = m_b*IT_1298;
    const complex_t IT_1423 = m_b*IT_1297;
    const complex_t IT_1424 = mty::lt::C0iC(12, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0060, IT_0507, mty::lt::reg_int);
    const complex_t IT_1425 = m_b*IT_1424;
    const complex_t IT_1426 = mty::lt::C0iC(18, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0060, IT_0507, mty::lt::reg_int);
    const complex_t IT_1427 = m_b*IT_1426;
    const complex_t IT_1428 = IT_1422 + IT_1423 + IT_1425 + IT_1427;
    const complex_t IT_1429 = mty::lt::C0iC(15, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0060, IT_0507, mty::lt::reg_int);
    const complex_t IT_1430 = m_s*IT_1429;
    const complex_t IT_1431 = m_b*IT_1429;
    const complex_t IT_1432 = m_s*IT_1297;
    const complex_t IT_1433 = m_s*IT_1424;
    const complex_t IT_1434 = -IT_1430 + 2*IT_1431 + -IT_1432 + -IT_1433;
    const complex_t IT_1435 = IT_1428 + IT_1434;
    const complex_t IT_1436 = IT_1421*IT_1435;
    const complex_t IT_1437 = IT_0321*IT_0416*IT_0427;
    const complex_t IT_1438 = 0.101321183642338*IT_1437;
    const complex_t IT_1439 = m_b*IT_1150;
    const complex_t IT_1440 = m_b*IT_1151;
    const complex_t IT_1441 = mty::lt::C0iC(12, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0060, IT_0430, mty::lt::reg_int);
    const complex_t IT_1442 = m_b*IT_1441;
    const complex_t IT_1443 = mty::lt::C0iC(18, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0060, IT_0430, mty::lt::reg_int);
    const complex_t IT_1444 = m_b*IT_1443;
    const complex_t IT_1445 = IT_1439 + IT_1440 + IT_1442 + IT_1444;
    const complex_t IT_1446 = m_s*IT_1151;
    const complex_t IT_1447 = m_s*IT_1441;
    const complex_t IT_1448 = mty::lt::C0iC(15, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0060, IT_0430, mty::lt::reg_int);
    const complex_t IT_1449 = m_b*IT_1448;
    const complex_t IT_1450 = m_s*IT_1448;
    const complex_t IT_1451 = -IT_1446 + -IT_1447 + 2*IT_1449 + -IT_1450;
    const complex_t IT_1452 = IT_1445 + IT_1451;
    const complex_t IT_1453 = IT_1438*IT_1452;
    const complex_t IT_1454 = IT_0193*IT_0321*IT_1104;
    const complex_t IT_1455 = 0.101321183642338*IT_1454;
    const complex_t IT_1456 = IT_0955*IT_1455;
    const complex_t IT_1457 = IT_0427*IT_0614;
    const complex_t IT_1458 = IT_0017*IT_1457;
    const complex_t IT_1459 = mty::lt::C0iC(0, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0060, IT_0430, IT_0430, mty::lt::reg_int);
    const complex_t IT_1460 = IT_0059*IT_1459;
    const complex_t IT_1461 = IT_0064*IT_0431;
    const complex_t IT_1462 = IT_1460 + IT_1461;
    const complex_t IT_1463 = IT_1458*IT_1462;
    const complex_t IT_1464 = IT_0141*IT_0157*IT_0561;
    const complex_t IT_1465 = 0.101321183642338*IT_1464;
    const complex_t IT_1466 = IT_0776*IT_1465;
    const complex_t IT_1467 = IT_0561*IT_0712*IT_1245;
    const complex_t IT_1468 = IT_0174*IT_1467;
    const complex_t IT_1469 = IT_1292*IT_1468;
    const complex_t IT_1470 = IT_0332*IT_0991;
    const complex_t IT_1471 = IT_0017*IT_1470;
    const complex_t IT_1472 = IT_1312*IT_1471;
    const complex_t IT_1473 = IT_0848*IT_1039;
    const complex_t IT_1474 = IT_0174*IT_1473;
    const complex_t IT_1475 = IT_1317*IT_1474;
    const complex_t IT_1476 = mty::lt::C0iC(0, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0089, IT_0507, mty::lt::reg_int);
    const complex_t IT_1477 = IT_1349 + IT_1351 + IT_1476;
    const complex_t IT_1478 = IT_0584*IT_1346*IT_1396;
    const complex_t IT_1479 = IT_0017*IT_1478;
    const complex_t IT_1480 = IT_1477*IT_1479;
    const complex_t IT_1481 = IT_0075*IT_0086*IT_0595;
    const complex_t IT_1482 = 0.101321183642338*IT_1481;
    const complex_t IT_1483 = IT_0665*IT_1482;
    const complex_t IT_1484 = IT_0493*IT_0548*IT_1093;
    const complex_t IT_1485 = 0.101321183642338*IT_1484;
    const complex_t IT_1486 = IT_1435*IT_1485;
    const complex_t IT_1487 = IT_0056*IT_0572*IT_1093;
    const complex_t IT_1488 = 0.101321183642338*IT_1487;
    const complex_t IT_1489 = m_b*IT_0675;
    const complex_t IT_1490 = m_b*IT_0676;
    const complex_t IT_1491 = mty::lt::C0iC(12, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0060, IT_0061, mty::lt::reg_int);
    const complex_t IT_1492 = m_b*IT_1491;
    const complex_t IT_1493 = mty::lt::C0iC(18, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0060, IT_0061, mty::lt::reg_int);
    const complex_t IT_1494 = m_b*IT_1493;
    const complex_t IT_1495 = IT_1489 + IT_1490 + IT_1492 + IT_1494;
    const complex_t IT_1496 = m_s*IT_0675;
    const complex_t IT_1497 = m_s*IT_1491;
    const complex_t IT_1498 = mty::lt::C0iC(15, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0060, IT_0061, mty::lt::reg_int);
    const complex_t IT_1499 = m_s*IT_1498;
    const complex_t IT_1500 = m_b*IT_1498;
    const complex_t IT_1501 = -IT_1496 + -IT_1497 + -IT_1499 + 2*IT_1500;
    const complex_t IT_1502 = IT_1495 + IT_1501;
    const complex_t IT_1503 = IT_1488*IT_1502;
    const complex_t IT_1504 = IT_0427*IT_0455*IT_1093;
    const complex_t IT_1505 = 0.101321183642338*IT_1504;
    const complex_t IT_1506 = mty::lt::C0iC(12, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0060, IT_0430, mty::lt::reg_int);
    const complex_t IT_1507 = m_b*IT_1506;
    const complex_t IT_1508 = m_b*IT_0699;
    const complex_t IT_1509 = m_b*IT_0698;
    const complex_t IT_1510 = mty::lt::C0iC(18, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0060, IT_0430, mty::lt::reg_int);
    const complex_t IT_1511 = m_b*IT_1510;
    const complex_t IT_1512 = IT_1507 + IT_1508 + IT_1509 + IT_1511;
    const complex_t IT_1513 = m_s*IT_1506;
    const complex_t IT_1514 = m_s*IT_0698;
    const complex_t IT_1515 = mty::lt::C0iC(15, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0060, IT_0430, mty::lt::reg_int);
    const complex_t IT_1516 = m_s*IT_1515;
    const complex_t IT_1517 = m_b*IT_1515;
    const complex_t IT_1518 = -IT_1513 + -IT_1514 + -IT_1516 + 2*IT_1517;
    const complex_t IT_1519 = IT_1512 + IT_1518;
    const complex_t IT_1520 = IT_1505*IT_1519;
    const complex_t IT_1521 = IT_0229*IT_1093*IT_1104;
    const complex_t IT_1522 = 0.101321183642338*IT_1521;
    const complex_t IT_1523 = IT_1415*IT_1522;
    const complex_t IT_1524 = IT_0332*IT_0837*IT_1093;
    const complex_t IT_1525 = 0.101321183642338*IT_1524;
    const complex_t IT_1526 = m_b*IT_1134;
    const complex_t IT_1527 = m_b*IT_1135;
    const complex_t IT_1528 = mty::lt::C0iC(12, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0060, IT_0346, mty::lt::reg_int);
    const complex_t IT_1529 = m_b*IT_1528;
    const complex_t IT_1530 = mty::lt::C0iC(18, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0060, IT_0346, mty::lt::reg_int);
    const complex_t IT_1531 = m_b*IT_1530;
    const complex_t IT_1532 = IT_1526 + IT_1527 + IT_1529 + IT_1531;
    const complex_t IT_1533 = m_s*IT_1134;
    const complex_t IT_1534 = m_s*IT_1528;
    const complex_t IT_1535 = mty::lt::C0iC(15, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0060, IT_0346, mty::lt::reg_int);
    const complex_t IT_1536 = m_s*IT_1535;
    const complex_t IT_1537 = m_b*IT_1535;
    const complex_t IT_1538 = -IT_1533 + -IT_1534 + -IT_1536 + 2*IT_1537;
    const complex_t IT_1539 = IT_1532 + IT_1538;
    const complex_t IT_1540 = IT_1525*IT_1539;
    const complex_t IT_1541 = IT_0374*IT_0712*IT_1093;
    const complex_t IT_1542 = 0.101321183642338*IT_1541;
    const complex_t IT_1543 = mty::lt::C0iC(18, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0060, IT_0388, mty::lt::reg_int);
    const complex_t IT_1544 = m_b*IT_1543;
    const complex_t IT_1545 = m_b*IT_0733;
    const complex_t IT_1546 = m_b*IT_0732;
    const complex_t IT_1547 = mty::lt::C0iC(12, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0060, IT_0388, mty::lt::reg_int);
    const complex_t IT_1548 = m_b*IT_1547;
    const complex_t IT_1549 = IT_1544 + IT_1545 + IT_1546 + IT_1548;
    const complex_t IT_1550 = mty::lt::C0iC(15, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0089, IT_0060, IT_0388, mty::lt::reg_int);
    const complex_t IT_1551 = m_s*IT_1550;
    const complex_t IT_1552 = m_b*IT_1550;
    const complex_t IT_1553 = m_s*IT_1547;
    const complex_t IT_1554 = m_s*IT_0732;
    const complex_t IT_1555 = -IT_1551 + 2*IT_1552 + -IT_1553 + -IT_1554;
    const complex_t IT_1556 = IT_1549 + IT_1555;
    const complex_t IT_1557 = IT_1542*IT_1556;
    const complex_t IT_1558 = IT_0300*IT_0807;
    const complex_t IT_1559 = IT_0174*IT_1558;
    const complex_t IT_1560 = IT_0252*IT_1559;
    const complex_t IT_1561 = IT_0075*IT_0504*IT_0537;
    const complex_t IT_1562 = IT_0017*IT_1561;
    const complex_t IT_1563 = IT_1477*IT_1562;
    const complex_t IT_1564 = IT_0493*IT_0671*IT_0918;
    const complex_t IT_1565 = IT_0174*IT_1564;
    const complex_t IT_1566 = IT_1299*IT_1565;
    const complex_t IT_1567 = IT_0671*IT_0694*IT_1127;
    const complex_t IT_1568 = 0.101321183642338*IT_1567;
    const complex_t IT_1569 = IT_1519*IT_1568;
    const complex_t IT_1570 = IT_0110*IT_0141*IT_0671;
    const complex_t IT_1571 = 0.101321183642338*IT_1570;
    const complex_t IT_1572 = IT_1502*IT_1571;
    const complex_t IT_1573 = IT_0056*IT_0141*IT_1093;
    const complex_t IT_1574 = IT_0174*IT_1573;
    const complex_t IT_1575 = IT_0677*IT_1574;
    const complex_t IT_1576 = IT_0075*IT_0416*IT_0466;
    const complex_t IT_1577 = 0.101321183642338*IT_1576;
    const complex_t IT_1578 = IT_0648*IT_1577;
    const complex_t IT_1579 = IT_0075*IT_0193*IT_0807;
    const complex_t IT_1580 = 0.101321183642338*IT_1579;
    const complex_t IT_1581 = m_b*IT_1218;
    const complex_t IT_1582 = mty::lt::C0iC(12, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0089, IT_0212, mty::lt::reg_int);
    const complex_t IT_1583 = m_b*IT_1582;
    const complex_t IT_1584 = mty::lt::C0iC(18, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0089, IT_0212, mty::lt::reg_int);
    const complex_t IT_1585 = m_b*IT_1584;
    const complex_t IT_1586 = m_b*IT_1219;
    const complex_t IT_1587 = IT_1581 + IT_1583 + IT_1585 + IT_1586;
    const complex_t IT_1588 = m_s*IT_1219;
    const complex_t IT_1589 = mty::lt::C0iC(15, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0089, IT_0212, mty::lt::reg_int);
    const complex_t IT_1590 = m_b*IT_1589;
    const complex_t IT_1591 = m_s*IT_1582;
    const complex_t IT_1592 = m_s*IT_1589;
    const complex_t IT_1593 = -IT_1588 + 2*IT_1590 + -IT_1591 + -IT_1592;
    const complex_t IT_1594 = IT_1587 + IT_1593;
    const complex_t IT_1595 = IT_1580*IT_1594;
    const complex_t IT_1596 = IT_0075*IT_0385*IT_0878;
    const complex_t IT_1597 = 0.101321183642338*IT_1596;
    const complex_t IT_1598 = m_s*IT_1248;
    const complex_t IT_1599 = mty::lt::C0iC(15, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0089, IT_0388, mty::lt::reg_int);
    const complex_t IT_1600 = m_b*IT_1599;
    const complex_t IT_1601 = mty::lt::C0iC(12, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0089, IT_0388, mty::lt::reg_int);
    const complex_t IT_1602 = m_s*IT_1601;
    const complex_t IT_1603 = m_s*IT_1599;
    const complex_t IT_1604 = -IT_1598 + 2*IT_1600 + -IT_1602 + -IT_1603;
    const complex_t IT_1605 = m_b*IT_1248;
    const complex_t IT_1606 = m_b*IT_1601;
    const complex_t IT_1607 = m_b*IT_1249;
    const complex_t IT_1608 = mty::lt::C0iC(18, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0089, IT_0388, mty::lt::reg_int);
    const complex_t IT_1609 = m_b*IT_1608;
    const complex_t IT_1610 = IT_1605 + IT_1606 + IT_1607 + IT_1609;
    const complex_t IT_1611 = IT_1604 + IT_1610;
    const complex_t IT_1612 = IT_1597*IT_1611;
    const complex_t IT_1613 = IT_0416*IT_0694;
    const complex_t IT_1614 = IT_0017*IT_1613;
    const complex_t IT_1615 = IT_1462*IT_1614;
    const complex_t IT_1616 = IT_0614*IT_0694;
    const complex_t IT_1617 = 0.101321183642338*IT_1616;
    const complex_t IT_1618 = IT_0443*IT_1617;
    const complex_t IT_1619 = IT_0466*IT_1127;
    const complex_t IT_1620 = IT_0174*IT_1619;
    const complex_t IT_1621 = IT_1191*IT_1620;
    const complex_t IT_1622 = IT_0584*IT_0991*IT_1055;
    const complex_t IT_1623 = 0.101321183642338*IT_1622;
    const complex_t IT_1624 = m_b*IT_1225;
    const complex_t IT_1625 = m_b*IT_1226;
    const complex_t IT_1626 = mty::lt::C0iC(12, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0089, IT_0346, mty::lt::reg_int);
    const complex_t IT_1627 = m_b*IT_1626;
    const complex_t IT_1628 = mty::lt::C0iC(18, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0089, IT_0346, mty::lt::reg_int);
    const complex_t IT_1629 = m_b*IT_1628;
    const complex_t IT_1630 = IT_1624 + IT_1625 + IT_1627 + IT_1629;
    const complex_t IT_1631 = m_s*IT_1225;
    const complex_t IT_1632 = m_s*IT_1626;
    const complex_t IT_1633 = mty::lt::C0iC(15, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0089, IT_0346, mty::lt::reg_int);
    const complex_t IT_1634 = m_b*IT_1633;
    const complex_t IT_1635 = m_s*IT_1633;
    const complex_t IT_1636 = -IT_1631 + -IT_1632 + 2*IT_1634 + -IT_1635;
    const complex_t IT_1637 = IT_1630 + IT_1636;
    const complex_t IT_1638 = IT_1623*IT_1637;
    const complex_t IT_1639 = IT_0416*IT_0694*IT_0902;
    const complex_t IT_1640 = IT_0017*IT_1639;
    const complex_t IT_1641 = IT_1153*IT_1640;
    const complex_t IT_1642 = IT_0193*IT_0209*IT_0902;
    const complex_t IT_1643 = IT_0017*IT_1642;
    const complex_t IT_1644 = IT_1158*IT_1643;
    const complex_t IT_1645 = IT_0614*IT_0694*IT_0902;
    const complex_t IT_1646 = 0.101321183642338*IT_1645;
    const complex_t IT_1647 = IT_1452*IT_1646;
    const complex_t IT_1648 = IT_0837*IT_0848;
    const complex_t IT_1649 = 0.101321183642338*IT_1648;
    const complex_t IT_1650 = IT_1070*IT_1649;
    const complex_t IT_1651 = IT_0455*IT_0561*IT_0630;
    const complex_t IT_1652 = IT_0174*IT_1651;
    const complex_t IT_1653 = IT_1261*IT_1652;
    const complex_t IT_1654 = IT_0110*IT_0595;
    const complex_t IT_1655 = IT_0017*IT_1654;
    const complex_t IT_1656 = IT_0067*IT_1655;
    const complex_t IT_1657 = IT_0056*IT_0595;
    const complex_t IT_1658 = 0.101321183642338*IT_1657;
    const complex_t IT_1659 = IT_0124*IT_1658;
    const complex_t IT_1660 = IT_0086*IT_0572;
    const complex_t IT_1661 = 0.101321183642338*IT_1660;
    const complex_t IT_1662 = IT_0172*IT_1661;
    const complex_t IT_1663 = IT_0561*IT_0630*IT_1127;
    const complex_t IT_1664 = 0.101321183642338*IT_1663;
    const complex_t IT_1665 = IT_0795*IT_1664;
    const complex_t IT_1666 = IT_0245*IT_0269*IT_0584;
    const complex_t IT_1667 = 0.101321183642338*IT_1666;
    const complex_t IT_1668 = IT_1594*IT_1667;
    const complex_t IT_1669 = IT_0584*IT_1180*IT_1245;
    const complex_t IT_1670 = 0.101321183642338*IT_1669;
    const complex_t IT_1671 = IT_1611*IT_1670;
    const complex_t IT_1672 = IT_0110*IT_0595*IT_0902;
    const complex_t IT_1673 = IT_0017*IT_1672;
    const complex_t IT_1674 = IT_1141*IT_1673;
    const complex_t IT_1675 = IT_0728*IT_1180;
    const complex_t IT_1676 = 0.101321183642338*IT_1675;
    const complex_t IT_1677 = IT_1086*IT_1676;
    const complex_t IT_1678 = IT_0059*IT_1332;
    const complex_t IT_1679 = m_b*IT_1678;
    const complex_t IT_1680 = mty::lt::C0iC(12, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0089, IT_0388, IT_0388, mty::lt::reg_int);
    const complex_t IT_1681 = IT_0064*IT_1680;
    const complex_t IT_1682 = m_b*IT_1681;
    const complex_t IT_1683 = mty::lt::C0iC(6, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0089, IT_0388, IT_0388, mty::lt::reg_int);
    const complex_t IT_1684 = IT_0059*IT_1683;
    const complex_t IT_1685 = m_s*IT_1684;
    const complex_t IT_1686 = mty::lt::C0iC(15, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0089, IT_0388, IT_0388, mty::lt::reg_int);
    const complex_t IT_1687 = IT_0064*IT_1686;
    const complex_t IT_1688 = m_s*IT_1687;
    const complex_t IT_1689 = IT_1679 + IT_1682 + IT_1685 + IT_1688;
    const complex_t IT_1690 = IT_1245*IT_1288;
    const complex_t IT_1691 = 0.101321183642338*IT_1690;
    const complex_t IT_1692 = IT_1689*IT_1691;
    const complex_t IT_1693 = IT_0229*IT_0245*IT_0561;
    const complex_t IT_1694 = IT_0174*IT_1693;
    const complex_t IT_1695 = IT_1266*IT_1694;
    const complex_t IT_1696 = IT_0561*IT_0837*IT_1055;
    const complex_t IT_1697 = IT_0174*IT_1696;
    const complex_t IT_1698 = IT_1271*IT_1697;
    const complex_t IT_1699 = IT_0245*IT_0300*IT_0561;
    const complex_t IT_1700 = 0.101321183642338*IT_1699;
    const complex_t IT_1701 = IT_0825*IT_1700;
    const complex_t IT_1702 = IT_0671*IT_0728*IT_1288;
    const complex_t IT_1703 = 0.101321183642338*IT_1702;
    const complex_t IT_1704 = IT_1556*IT_1703;
    const complex_t IT_1705 = IT_0902*IT_0918*IT_1346;
    const complex_t IT_1706 = IT_0017*IT_1705;
    const complex_t IT_1707 = IT_1146*IT_1706;
    const complex_t IT_1708 = IT_0902*IT_0991*IT_1007;
    const complex_t IT_1709 = 0.101321183642338*IT_1708;
    const complex_t IT_1710 = IT_0362*IT_1709;
    const complex_t IT_1711 = IT_0728*IT_0902*IT_1180;
    const complex_t IT_1712 = 0.101321183642338*IT_1711;
    const complex_t IT_1713 = IT_0404*IT_1712;
    const complex_t IT_1714 = IT_0332*IT_0343;
    const complex_t IT_1715 = 0.101321183642338*IT_1714;
    const complex_t IT_1716 = IT_1022*IT_1715;
    const complex_t IT_1717 = IT_0712*IT_0878;
    const complex_t IT_1718 = 0.101321183642338*IT_1717;
    const complex_t IT_1719 = IT_1689*IT_1718;
    const complex_t IT_1720 = IT_0374*IT_1180;
    const complex_t IT_1721 = IT_0017*IT_1720;
    const complex_t IT_1722 = IT_1326*IT_1721;
    const complex_t IT_1723 = IT_0385*IT_0728*IT_0902;
    const complex_t IT_1724 = IT_0017*IT_1723;
    const complex_t IT_1725 = IT_1184*IT_1724;
    const complex_t IT_1726 = IT_0416*IT_0584*IT_0630;
    const complex_t IT_1727 = IT_0017*IT_1726;
    const complex_t IT_1728 = mty::lt::C0iC(0, IT_0008 + IT_0009 + (-2)*s_12,
       IT_0009, IT_0008, IT_0060, IT_0089, IT_0430, mty::lt::reg_int);
    const complex_t IT_1729 = IT_0635 + IT_0639 + IT_1728;
    const complex_t IT_1730 = IT_1727*IT_1729;
    const complex_t IT_1731 = IT_0075*IT_0848*IT_0991;
    const complex_t IT_1732 = IT_0017*IT_1731;
    const complex_t IT_1733 = IT_1228*IT_1732;
    const complex_t IT_1734 = IT_0712*IT_1245;
    const complex_t IT_1735 = IT_0174*IT_1734;
    const complex_t IT_1736 = IT_1334*IT_1735;
    const complex_t IT_1737 = IT_0157*IT_0572;
    const complex_t IT_1738 = IT_0174*IT_1737;
    const complex_t IT_1739 = IT_0180*IT_1738;
    const complex_t IT_1740 = IT_0193*IT_1104;
    const complex_t IT_1741 = 0.101321183642338*IT_1740;
    const complex_t IT_1742 = IT_0283*IT_1741;
    const complex_t IT_1743 = IT_0229*IT_0807;
    const complex_t IT_1744 = 0.101321183642338*IT_1743;
    const complex_t IT_1745 = IT_0314*IT_1744;
    const complex_t IT_1746 = IT_0075*IT_0343*IT_0848;
    const complex_t IT_1747 = 0.101321183642338*IT_1746;
    const complex_t IT_1748 = IT_1637*IT_1747;
    const complex_t IT_1749 = IT_0075*IT_0466*IT_0614;
    const complex_t IT_1750 = IT_0017*IT_1749;
    const complex_t IT_1751 = IT_1729*IT_1750;
    const complex_t IT_1752 = IT_0075*IT_0269*IT_0807;
    const complex_t IT_1753 = IT_0017*IT_1752;
    const complex_t IT_1754 = IT_1221*IT_1753;
    const complex_t IT_1755 = IT_0075*IT_0878*IT_1180;
    const complex_t IT_1756 = IT_0017*IT_1755;
    const complex_t IT_1757 = IT_1251*IT_1756;
    const complex_t IT_1758 = IT_0321*IT_0548*IT_1346;
    const complex_t IT_1759 = 0.101321183642338*IT_1758;
    const complex_t IT_1760 = IT_0936*IT_1759;
    const complex_t IT_1761 = IT_0056*IT_0321*IT_0595;
    const complex_t IT_1762 = 0.101321183642338*IT_1761;
    const complex_t IT_1763 = IT_0974*IT_1762;
    const complex_t IT_1764 = IT_0630*IT_1127;
    const complex_t IT_1765 = 0.101321183642338*IT_1764;
    const complex_t IT_1766 = IT_0481*IT_1765;
    const complex_t IT_1767 = IT_0918*IT_1346;
    const complex_t IT_1768 = IT_0017*IT_1767;
    const complex_t IT_1769 = IT_0555*IT_1768;
    const complex_t IT_1770 = IT_0493*IT_1396;
    const complex_t IT_1771 = IT_0174*IT_1770;
    const complex_t IT_1772 = IT_1214*IT_1771;
    const complex_t IT_1773 = IT_0537*IT_0918;
    const complex_t IT_1774 = 0.101321183642338*IT_1773;
    const complex_t IT_1775 = IT_1379*IT_1774;
    const complex_t IT_1776 = IT_1208*IT_1396;
    const complex_t IT_1777 = 0.101321183642338*IT_1776;
    const complex_t IT_1778 = IT_0520*IT_1777;
    const complex_t IT_1779 = IT_0493*IT_0561*IT_1396;
    const complex_t IT_1780 = IT_0174*IT_1779;
    const complex_t IT_1781 = IT_1256*IT_1780;
    const complex_t IT_1782 = IT_0561*IT_1208*IT_1396;
    const complex_t IT_1783 = 0.101321183642338*IT_1782;
    const complex_t IT_1784 = IT_0759*IT_1783;
    const complex_t IT_1785 = IT_0561*IT_1039*IT_1055;
    const complex_t IT_1786 = 0.101321183642338*IT_1785;
    const complex_t IT_1787 = IT_0866*IT_1786;
    const complex_t IT_1788 = IT_0561*IT_1245*IT_1288;
    const complex_t IT_1789 = 0.101321183642338*IT_1788;
    const complex_t IT_1790 = IT_0896*IT_1789;
    const complex_t IT_1791 = IT_0209*IT_0229*IT_0671;
    const complex_t IT_1792 = IT_0174*IT_1791;
    const complex_t IT_1793 = IT_1110*IT_1792;
    const complex_t IT_1794 = IT_0671*IT_0837*IT_1007;
    const complex_t IT_1795 = IT_0174*IT_1794;
    const complex_t IT_1796 = IT_1136*IT_1795;
    const complex_t IT_1797 = IT_0671*IT_1007*IT_1039;
    const complex_t IT_1798 = 0.101321183642338*IT_1797;
    const complex_t IT_1799 = IT_1539*IT_1798;
    const complex_t IT_1800 = IT_0086*IT_0141*IT_0741;
    const complex_t IT_1801 = IT_0174*IT_1800;
    const complex_t IT_1802 = IT_0578*IT_1801;
    const complex_t IT_1803 = IT_0068 + 2*IT_0094 + IT_0125 + IT_0173 +
       IT_0181 + IT_0218 + IT_0253 + IT_0284 + IT_0315 + 2*IT_0363 + 2*IT_0405 +
       IT_0444 + IT_0482 + IT_0521 + IT_0556 + (-2)*IT_0579 + (-2)*IT_0598 + (-2
      )*IT_0649 + (-2)*IT_0666 + (-2)*IT_0678 + (-2)*IT_0701 + (-2)*IT_0735 + 2
      *IT_0760 + 2*IT_0777 + 2*IT_0796 + 2*IT_0826 + 2*IT_0867 + 2*IT_0897 + (-2
      )*IT_0937 + (-2)*IT_0956 + (-2)*IT_0975 + IT_1023 + IT_1071 + IT_1087 + 2
      *IT_1111 + 2*IT_1130 + 2*IT_1137 + 2*IT_1142 + 2*IT_1147 + 2*IT_1154 + 2
      *IT_1159 + 2*IT_1164 + 2*IT_1185 + IT_1192 + IT_1215 + (-2)*IT_1222 + (-2)
      *IT_1229 + (-2)*IT_1252 + 2*IT_1257 + 2*IT_1262 + 2*IT_1267 + 2*IT_1272 +
       2*IT_1293 + 2*IT_1300 + 2*IT_1303 + (-2)*IT_1306 + IT_1313 + IT_1320 +
       IT_1327 + IT_1335 + 2*IT_1365 + IT_1380 + (-2)*IT_1399 + (-2)*IT_1416 +
       IT_1419 + (-2)*IT_1436 + 2*IT_1453 + 2*IT_1456 + IT_1463 + (-2)*IT_1466 +
       (-2)*IT_1469 + IT_1472 + IT_1475 + (-2)*IT_1480 + 2*IT_1483 + 2*IT_1486 +
       2*IT_1503 + 2*IT_1520 + 2*IT_1523 + 2*IT_1540 + 2*IT_1557 + IT_1560 + 2
      *IT_1563 + (-2)*IT_1566 + (-2)*IT_1569 + (-2)*IT_1572 + 2*IT_1575 + 2
      *IT_1578 + 2*IT_1595 + 2*IT_1612 + IT_1615 + IT_1618 + IT_1621 + (-2)
      *IT_1638 + (-2)*IT_1641 + (-2)*IT_1644 + (-2)*IT_1647 + IT_1650 + (-2)
      *IT_1653 + IT_1656 + IT_1659 + IT_1662 + (-2)*IT_1665 + (-2)*IT_1668 + (-2
      )*IT_1671 + (-2)*IT_1674 + IT_1677 + IT_1692 + (-2)*IT_1695 + (-2)*IT_1698
       + (-2)*IT_1701 + (-2)*IT_1704 + (-2)*IT_1707 + (-2)*IT_1710 + (-2)
      *IT_1713 + IT_1716 + IT_1719 + IT_1722 + (-2)*IT_1725 + (-2)*IT_1730 + 2
      *IT_1733 + IT_1736 + IT_1739 + IT_1742 + IT_1745 + 2*IT_1748 + 2*IT_1751 +
       2*IT_1754 + 2*IT_1757 + 2*IT_1760 + 2*IT_1763 + IT_1766 + IT_1769 +
       IT_1772 + IT_1775 + IT_1778 + (-2)*IT_1781 + (-2)*IT_1784 + (-2)*IT_1787 
      + (-2)*IT_1790 + (-2)*IT_1793 + (-2)*IT_1796 + (-2)*IT_1799 + 2*IT_1802;
    const complex_t IT_1804 = cpowq(IT_0018, 2);
    const complex_t IT_1805 = IT_1803*IT_1804;
    const complex_t IT_1806 = IT_0016*IT_1805;
    const complex_t IT_1807 = IT_0321*IT_0385*IT_0728;
    const complex_t IT_1808 = IT_0017*IT_1807;
    const complex_t IT_1809 = IT_0391*IT_1808;
    const complex_t IT_1810 = IT_0300*IT_0671*IT_1104;
    const complex_t IT_1811 = IT_0017*IT_1810;
    const complex_t IT_1812 = IT_1108*IT_1811;
    const complex_t IT_1813 = IT_0427*IT_0671*IT_1127;
    const complex_t IT_1814 = IT_0017*IT_1813;
    const complex_t IT_1815 = IT_0698*IT_1814;
    const complex_t IT_1816 = IT_0374*IT_0671*IT_1288;
    const complex_t IT_1817 = IT_0017*IT_1816;
    const complex_t IT_1818 = IT_0732*IT_1817;
    const complex_t IT_1819 = IT_0075*IT_1346*IT_1396;
    const complex_t IT_1820 = IT_0174*IT_1819;
    const complex_t IT_1821 = IT_1349*IT_1820;
    const complex_t IT_1822 = IT_0117 + IT_0123;
    const complex_t IT_1823 = m_b*IT_0066;
    const complex_t IT_1824 = m_b*IT_0114;
    const complex_t IT_1825 = IT_0064*IT_0121;
    const complex_t IT_1826 = m_s*IT_1825;
    const complex_t IT_1827 = mty::lt::C0iC(18, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0060, IT_0061, IT_0061, mty::lt::reg_int);
    const complex_t IT_1828 = IT_0064*IT_1827;
    const complex_t IT_1829 = m_s*IT_1828;
    const complex_t IT_1830 = -IT_1823 + -IT_1824 + -IT_1826 + -IT_1829;
    const complex_t IT_1831 = IT_1822 + IT_1830;
    const complex_t IT_1832 = IT_0112*IT_1831;
    const complex_t IT_1833 = IT_0162 + IT_0165;
    const complex_t IT_1834 = m_b*IT_0170;
    const complex_t IT_1835 = m_b*IT_0177;
    const complex_t IT_1836 = mty::lt::C0iC(18, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0089, IT_0061, IT_0061, mty::lt::reg_int);
    const complex_t IT_1837 = IT_0064*IT_1836;
    const complex_t IT_1838 = m_s*IT_1837;
    const complex_t IT_1839 = IT_0064*IT_0160;
    const complex_t IT_1840 = m_s*IT_1839;
    const complex_t IT_1841 = -IT_1834 + -IT_1835 + -IT_1838 + -IT_1840;
    const complex_t IT_1842 = IT_1833 + IT_1841;
    const complex_t IT_1843 = IT_0159*IT_1842;
    const complex_t IT_1844 = IT_0064*IT_0178;
    const complex_t IT_1845 = -IT_1839 + -IT_1844;
    const complex_t IT_1846 = IT_0179 + IT_1845;
    const complex_t IT_1847 = IT_0176*IT_1846;
    const complex_t IT_1848 = IT_0064*IT_0215;
    const complex_t IT_1849 = IT_0064*IT_0277;
    const complex_t IT_1850 = -IT_1848 + -IT_1849;
    const complex_t IT_1851 = IT_0216 + IT_1850;
    const complex_t IT_1852 = IT_0211*IT_1851;
    const complex_t IT_1853 = IT_0943 + IT_0945 + IT_0951;
    const complex_t IT_1854 = -IT_0949 + -IT_0952;
    const complex_t IT_1855 = IT_1853 + IT_1854;
    const complex_t IT_1856 = IT_1455*IT_1855;
    const complex_t IT_1857 = IT_0352 + IT_0354 + IT_0359;
    const complex_t IT_1858 = -IT_0356 + -IT_0357;
    const complex_t IT_1859 = IT_1857 + IT_1858;
    const complex_t IT_1860 = IT_0345*IT_1859;
    const complex_t IT_1861 = IT_0392 + IT_0394 + IT_0399;
    const complex_t IT_1862 = -IT_0400 + -IT_0402;
    const complex_t IT_1863 = IT_1861 + IT_1862;
    const complex_t IT_1864 = IT_0387*IT_1863;
    const complex_t IT_1865 = IT_1369 + IT_1375;
    const complex_t IT_1866 = m_b*IT_0552;
    const complex_t IT_1867 = m_b*IT_1377;
    const complex_t IT_1868 = IT_0064*IT_1373;
    const complex_t IT_1869 = m_s*IT_1868;
    const complex_t IT_1870 = mty::lt::C0iC(18, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0060, IT_0507, IT_0507, mty::lt::reg_int);
    const complex_t IT_1871 = IT_0064*IT_1870;
    const complex_t IT_1872 = m_s*IT_1871;
    const complex_t IT_1873 = -IT_1866 + -IT_1867 + -IT_1869 + -IT_1872;
    const complex_t IT_1874 = IT_1865 + IT_1873;
    const complex_t IT_1875 = IT_1367*IT_1874;
    const complex_t IT_1876 = IT_0510 + IT_0516;
    const complex_t IT_1877 = m_b*IT_1213;
    const complex_t IT_1878 = m_b*IT_0518;
    const complex_t IT_1879 = IT_0064*IT_0514;
    const complex_t IT_1880 = m_s*IT_1879;
    const complex_t IT_1881 = mty::lt::C0iC(18, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0089, IT_0507, IT_0507, mty::lt::reg_int);
    const complex_t IT_1882 = IT_0064*IT_1881;
    const complex_t IT_1883 = m_s*IT_1882;
    const complex_t IT_1884 = -IT_1877 + -IT_1878 + -IT_1880 + -IT_1883;
    const complex_t IT_1885 = IT_1876 + IT_1884;
    const complex_t IT_1886 = IT_0506*IT_1885;
    const complex_t IT_1887 = IT_0064*IT_0553;
    const complex_t IT_1888 = -IT_1868 + -IT_1887;
    const complex_t IT_1889 = IT_0554 + IT_1888;
    const complex_t IT_1890 = IT_0550*IT_1889;
    const complex_t IT_1891 = IT_0504*IT_0561*IT_1208;
    const complex_t IT_1892 = IT_0174*IT_1891;
    const complex_t IT_1893 = IT_0746*IT_1892;
    const complex_t IT_1894 = IT_0466*IT_0561*IT_1127;
    const complex_t IT_1895 = IT_0174*IT_1894;
    const complex_t IT_1896 = IT_0782*IT_1895;
    const complex_t IT_1897 = IT_0300*IT_0561*IT_0807;
    const complex_t IT_1898 = IT_0174*IT_1897;
    const complex_t IT_1899 = IT_0812*IT_1898;
    const complex_t IT_1900 = IT_0504*IT_0537*IT_0584;
    const complex_t IT_1901 = IT_0174*IT_1900;
    const complex_t IT_1902 = IT_1349*IT_1901;
    const complex_t IT_1903 = IT_0584*IT_0848*IT_0991;
    const complex_t IT_1904 = IT_0174*IT_1903;
    const complex_t IT_1905 = IT_1225*IT_1904;
    const complex_t IT_1906 = IT_1402 + IT_1405 + IT_1413;
    const complex_t IT_1907 = -IT_1409 + -IT_1410;
    const complex_t IT_1908 = IT_1906 + IT_1907;
    const complex_t IT_1909 = IT_1401*IT_1908;
    const complex_t IT_1910 = IT_1489 + IT_1492 + IT_1500;
    const complex_t IT_1911 = -IT_1496 + -IT_1497;
    const complex_t IT_1912 = IT_1910 + IT_1911;
    const complex_t IT_1913 = IT_1571*IT_1912;
    const complex_t IT_1914 = IT_1488*IT_1912;
    const complex_t IT_1915 = IT_1507 + IT_1509 + IT_1517;
    const complex_t IT_1916 = -IT_1513 + -IT_1514;
    const complex_t IT_1917 = IT_1915 + IT_1916;
    const complex_t IT_1918 = IT_1505*IT_1917;
    const complex_t IT_1919 = IT_1709*IT_1859;
    const complex_t IT_1920 = IT_1712*IT_1863;
    const complex_t IT_1921 = IT_0064*IT_1310;
    const complex_t IT_1922 = IT_0064*IT_1010;
    const complex_t IT_1923 = -IT_1921 + -IT_1922;
    const complex_t IT_1924 = IT_1311 + IT_1923;
    const complex_t IT_1925 = IT_1308*IT_1924;
    const complex_t IT_1926 = IT_0075*IT_0416*IT_0630;
    const complex_t IT_1927 = IT_0174*IT_1926;
    const complex_t IT_1928 = IT_0635*IT_1927;
    const complex_t IT_1929 = IT_0075*IT_0343*IT_1055;
    const complex_t IT_1930 = IT_0174*IT_1929;
    const complex_t IT_1931 = IT_1225*IT_1930;
    const complex_t IT_1932 = IT_0075*IT_0385*IT_1245;
    const complex_t IT_1933 = IT_0174*IT_1932;
    const complex_t IT_1934 = IT_1248*IT_1933;
    const complex_t IT_1935 = IT_0110*IT_0321*IT_0595;
    const complex_t IT_1936 = IT_0017*IT_1935;
    const complex_t IT_1937 = IT_0959*IT_1936;
    const complex_t IT_1938 = IT_0321*IT_0343*IT_1007;
    const complex_t IT_1939 = IT_0017*IT_1938;
    const complex_t IT_1940 = IT_0351*IT_1939;
    const complex_t IT_1941 = IT_0064*IT_1211;
    const complex_t IT_1942 = -IT_1879 + -IT_1941;
    const complex_t IT_1943 = IT_1212 + IT_1942;
    const complex_t IT_1944 = IT_1771*IT_1943;
    const complex_t IT_1945 = IT_1210*IT_1943;
    const complex_t IT_1946 = IT_0561*IT_0848*IT_1039;
    const complex_t IT_1947 = IT_0174*IT_1946;
    const complex_t IT_1948 = IT_0853*IT_1947;
    const complex_t IT_1949 = IT_0561*IT_0878*IT_1288;
    const complex_t IT_1950 = IT_0174*IT_1949;
    const complex_t IT_1951 = IT_0885*IT_1950;
    const complex_t IT_1952 = IT_1568*IT_1917;
    const complex_t IT_1953 = IT_0455*IT_0694*IT_1093;
    const complex_t IT_1954 = IT_0017*IT_1953;
    const complex_t IT_1955 = IT_0698*IT_1954;
    const complex_t IT_1956 = IT_0056*IT_0141*IT_0671;
    const complex_t IT_1957 = IT_0017*IT_1956;
    const complex_t IT_1958 = IT_0675*IT_1957;
    const complex_t IT_1959 = IT_0763 + IT_0766 + IT_0774;
    const complex_t IT_1960 = -IT_0770 + -IT_0771;
    const complex_t IT_1961 = IT_1959 + IT_1960;
    const complex_t IT_1962 = IT_1465*IT_1961;
    const complex_t IT_1963 = IT_1583 + IT_1586 + IT_1590;
    const complex_t IT_1964 = -IT_1588 + -IT_1591;
    const complex_t IT_1965 = IT_1963 + IT_1964;
    const complex_t IT_1966 = IT_1667*IT_1965;
    const complex_t IT_1967 = IT_0636 + IT_0638 + IT_0646;
    const complex_t IT_1968 = -IT_0644 + -IT_0645;
    const complex_t IT_1969 = IT_1967 + IT_1968;
    const complex_t IT_1970 = IT_0632*IT_1969;
    const complex_t IT_1971 = IT_1423 + IT_1425 + IT_1431;
    const complex_t IT_1972 = -IT_1432 + -IT_1433;
    const complex_t IT_1973 = IT_1971 + IT_1972;
    const complex_t IT_1974 = IT_1421*IT_1973;
    const complex_t IT_1975 = IT_0886 + IT_0888 + IT_0891;
    const complex_t IT_1976 = -IT_0892 + -IT_0893;
    const complex_t IT_1977 = IT_1975 + IT_1976;
    const complex_t IT_1978 = IT_0880*IT_1977;
    const complex_t IT_1979 = IT_1522*IT_1908;
    const complex_t IT_1980 = IT_1526 + IT_1529 + IT_1537;
    const complex_t IT_1981 = -IT_1533 + -IT_1534;
    const complex_t IT_1982 = IT_1980 + IT_1981;
    const complex_t IT_1983 = IT_1525*IT_1982;
    const complex_t IT_1984 = IT_1546 + IT_1548 + IT_1552;
    const complex_t IT_1985 = -IT_1553 + -IT_1554;
    const complex_t IT_1986 = IT_1984 + IT_1985;
    const complex_t IT_1987 = IT_1542*IT_1986;
    const complex_t IT_1988 = IT_0064*IT_1058;
    const complex_t IT_1989 = IT_0064*IT_1314;
    const complex_t IT_1990 = -IT_1988 + -IT_1989;
    const complex_t IT_1991 = IT_1315 + IT_1990;
    const complex_t IT_1992 = IT_1319*IT_1991;
    const complex_t IT_1993 = IT_1012 + IT_1015;
    const complex_t IT_1994 = mty::lt::C0iC(18, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0060, IT_0346, IT_0346, mty::lt::reg_int);
    const complex_t IT_1995 = IT_0064*IT_1994;
    const complex_t IT_1996 = m_s*IT_1995;
    const complex_t IT_1997 = m_b*IT_1309;
    const complex_t IT_1998 = m_b*IT_1020;
    const complex_t IT_1999 = m_s*IT_1922;
    const complex_t IT_2000 = -IT_1996 + -IT_1997 + -IT_1998 + -IT_1999;
    const complex_t IT_2001 = IT_1993 + IT_2000;
    const complex_t IT_2002 = IT_1009*IT_2001;
    const complex_t IT_2003 = IT_1060 + IT_1063;
    const complex_t IT_2004 = m_b*IT_1068;
    const complex_t IT_2005 = m_s*IT_1988;
    const complex_t IT_2006 = mty::lt::C0iC(18, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0089, IT_0346, IT_0346, mty::lt::reg_int);
    const complex_t IT_2007 = IT_0064*IT_2006;
    const complex_t IT_2008 = m_s*IT_2007;
    const complex_t IT_2009 = m_b*IT_1316;
    const complex_t IT_2010 = -IT_2004 + -IT_2005 + -IT_2008 + -IT_2009;
    const complex_t IT_2011 = IT_2003 + IT_2010;
    const complex_t IT_2012 = IT_1057*IT_2011;
    const complex_t IT_2013 = IT_0332*IT_0902*IT_0991;
    const complex_t IT_2014 = IT_0017*IT_2013;
    const complex_t IT_2015 = IT_0351*IT_2014;
    const complex_t IT_2016 = IT_0584*IT_0878*IT_1180;
    const complex_t IT_2017 = IT_0174*IT_2016;
    const complex_t IT_2018 = IT_1248*IT_2017;
    const complex_t IT_2019 = IT_1471*IT_1924;
    const complex_t IT_2020 = IT_1474*IT_1991;
    const complex_t IT_2021 = IT_0157*IT_0572*IT_0741;
    const complex_t IT_2022 = IT_0174*IT_2021;
    const complex_t IT_2023 = IT_0576*IT_2022;
    const complex_t IT_2024 = IT_0427*IT_0614*IT_0902;
    const complex_t IT_2025 = IT_0017*IT_2024;
    const complex_t IT_2026 = IT_1151*IT_2025;
    const complex_t IT_2027 = IT_0075*IT_0193*IT_0245;
    const complex_t IT_2028 = IT_0174*IT_2027;
    const complex_t IT_2029 = IT_1219*IT_2028;
    const complex_t IT_2030 = IT_1774*IT_1874;
    const complex_t IT_2031 = IT_0064*IT_1323;
    const complex_t IT_2032 = IT_0064*IT_1080;
    const complex_t IT_2033 = -IT_2031 + -IT_2032;
    const complex_t IT_2034 = IT_1324 + IT_2033;
    const complex_t IT_2035 = IT_1721*IT_2034;
    const complex_t IT_2036 = IT_0062*IT_0064;
    const complex_t IT_2037 = -IT_1825 + -IT_2036;
    const complex_t IT_2038 = IT_0063 + IT_2037;
    const complex_t IT_2039 = IT_0058*IT_2038;
    const complex_t IT_2040 = IT_0273 + IT_0279;
    const complex_t IT_2041 = m_b*IT_0214;
    const complex_t IT_2042 = m_b*IT_0281;
    const complex_t IT_2043 = mty::lt::C0iC(18, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0060, IT_0212, IT_0212, mty::lt::reg_int);
    const complex_t IT_2044 = IT_0064*IT_2043;
    const complex_t IT_2045 = m_s*IT_2044;
    const complex_t IT_2046 = m_s*IT_1849;
    const complex_t IT_2047 = -IT_2041 + -IT_2042 + -IT_2045 + -IT_2046;
    const complex_t IT_2048 = IT_2040 + IT_2047;
    const complex_t IT_2049 = IT_0271*IT_2048;
    const complex_t IT_2050 = IT_0305 + IT_0307;
    const complex_t IT_2051 = m_b*IT_0251;
    const complex_t IT_2052 = m_b*IT_0312;
    const complex_t IT_2053 = mty::lt::C0iC(18, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0089, IT_0212, IT_0212, mty::lt::reg_int);
    const complex_t IT_2054 = IT_0064*IT_2053;
    const complex_t IT_2055 = m_s*IT_2054;
    const complex_t IT_2056 = IT_0064*IT_0303;
    const complex_t IT_2057 = m_s*IT_2056;
    const complex_t IT_2058 = -IT_2051 + -IT_2052 + -IT_2055 + -IT_2057;
    const complex_t IT_2059 = IT_2050 + IT_2058;
    const complex_t IT_2060 = IT_0302*IT_2059;
    const complex_t IT_2061 = IT_0075*IT_0157*IT_0595;
    const complex_t IT_2062 = IT_0174*IT_2061;
    const complex_t IT_2063 = IT_0090*IT_2062;
    const complex_t IT_2064 = IT_1580*IT_1965;
    const complex_t IT_2065 = IT_1624 + IT_1627 + IT_1634;
    const complex_t IT_2066 = -IT_1631 + -IT_1632;
    const complex_t IT_2067 = IT_2065 + IT_2066;
    const complex_t IT_2068 = IT_1747*IT_2067;
    const complex_t IT_2069 = IT_1600 + IT_1605 + IT_1606;
    const complex_t IT_2070 = -IT_1598 + -IT_1602;
    const complex_t IT_2071 = IT_2069 + IT_2070;
    const complex_t IT_2072 = IT_1597*IT_2071;
    const complex_t IT_2073 = IT_0924 + IT_0926 + IT_0932;
    const complex_t IT_2074 = -IT_0930 + -IT_0933;
    const complex_t IT_2075 = IT_2073 + IT_2074;
    const complex_t IT_2076 = IT_1759*IT_2075;
    const complex_t IT_2077 = IT_0962 + IT_0967 + IT_0972;
    const complex_t IT_2078 = -IT_0960 + -IT_0964;
    const complex_t IT_2079 = IT_2077 + IT_2078;
    const complex_t IT_2080 = IT_1762*IT_2079;
    const complex_t IT_2081 = IT_1440 + IT_1442 + IT_1449;
    const complex_t IT_2082 = -IT_1446 + -IT_1447;
    const complex_t IT_2083 = IT_2081 + IT_2082;
    const complex_t IT_2084 = IT_1438*IT_2083;
    const complex_t IT_2085 = IT_0433 + IT_0436;
    const complex_t IT_2086 = m_b*IT_1461;
    const complex_t IT_2087 = m_b*IT_0441;
    const complex_t IT_2088 = IT_0064*IT_0434;
    const complex_t IT_2089 = m_s*IT_2088;
    const complex_t IT_2090 = mty::lt::C0iC(18, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0060, IT_0430, IT_0430, mty::lt::reg_int);
    const complex_t IT_2091 = IT_0064*IT_2090;
    const complex_t IT_2092 = m_s*IT_2091;
    const complex_t IT_2093 = -IT_2086 + -IT_2087 + -IT_2089 + -IT_2092;
    const complex_t IT_2094 = IT_2085 + IT_2093;
    const complex_t IT_2095 = IT_0429*IT_2094;
    const complex_t IT_2096 = IT_0471 + IT_0477;
    const complex_t IT_2097 = m_b*IT_1190;
    const complex_t IT_2098 = m_b*IT_0479;
    const complex_t IT_2099 = IT_0064*IT_0475;
    const complex_t IT_2100 = m_s*IT_2099;
    const complex_t IT_2101 = mty::lt::C0iC(18, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0089, IT_0430, IT_0430, mty::lt::reg_int);
    const complex_t IT_2102 = IT_0064*IT_2101;
    const complex_t IT_2103 = m_s*IT_2102;
    const complex_t IT_2104 = -IT_2097 + -IT_2098 + -IT_2100 + -IT_2103;
    const complex_t IT_2105 = IT_2096 + IT_2104;
    const complex_t IT_2106 = IT_0468*IT_2105;
    const complex_t IT_2107 = -IT_0753 + -IT_0754;
    const complex_t IT_2108 = IT_0747 + IT_0749 + IT_0757;
    const complex_t IT_2109 = IT_2107 + IT_2108;
    const complex_t IT_2110 = IT_0743*IT_2109;
    const complex_t IT_2111 = IT_0762*IT_1961;
    const complex_t IT_2112 = IT_0783 + IT_0785 + IT_0793;
    const complex_t IT_2113 = -IT_0789 + -IT_0790;
    const complex_t IT_2114 = IT_2112 + IT_2113;
    const complex_t IT_2115 = IT_0779*IT_2114;
    const complex_t IT_2116 = IT_0813 + IT_0815 + IT_0823;
    const complex_t IT_2117 = -IT_0819 + -IT_0820;
    const complex_t IT_2118 = IT_2116 + IT_2117;
    const complex_t IT_2119 = IT_0809*IT_2118;
    const complex_t IT_2120 = IT_0854 + IT_0856 + IT_0864;
    const complex_t IT_2121 = -IT_0860 + -IT_0861;
    const complex_t IT_2122 = IT_2120 + IT_2121;
    const complex_t IT_2123 = IT_0850*IT_2122;
    const complex_t IT_2124 = IT_0958*IT_2079;
    const complex_t IT_2125 = IT_0064*IT_0248;
    const complex_t IT_2126 = -IT_2056 + -IT_2125;
    const complex_t IT_2127 = IT_0249 + IT_2126;
    const complex_t IT_2128 = IT_0247*IT_2127;
    const complex_t IT_2129 = IT_0064*IT_1459;
    const complex_t IT_2130 = -IT_2088 + -IT_2129;
    const complex_t IT_2131 = IT_1460 + IT_2130;
    const complex_t IT_2132 = IT_1458*IT_2131;
    const complex_t IT_2133 = IT_0466*IT_0584*IT_0614;
    const complex_t IT_2134 = IT_0174*IT_2133;
    const complex_t IT_2135 = IT_0635*IT_2134;
    const complex_t IT_2136 = IT_0269*IT_0584*IT_0807;
    const complex_t IT_2137 = IT_0174*IT_2136;
    const complex_t IT_2138 = IT_1219*IT_2137;
    const complex_t IT_2139 = IT_0548*IT_0671*IT_1208;
    const complex_t IT_2140 = IT_0017*IT_2139;
    const complex_t IT_2141 = IT_1297*IT_2140;
    const complex_t IT_2142 = IT_0374*IT_0902*IT_1180;
    const complex_t IT_2143 = IT_0017*IT_2142;
    const complex_t IT_2144 = IT_0391*IT_2143;
    const complex_t IT_2145 = IT_0064*IT_1188;
    const complex_t IT_2146 = -IT_2099 + -IT_2145;
    const complex_t IT_2147 = IT_1189 + IT_2146;
    const complex_t IT_2148 = IT_1620*IT_2147;
    const complex_t IT_2149 = IT_1798*IT_1982;
    const complex_t IT_2150 = IT_0920*IT_2075;
    const complex_t IT_2151 = IT_0939*IT_1855;
    const complex_t IT_2152 = IT_1646*IT_2083;
    const complex_t IT_2153 = IT_1076 + IT_1082;
    const complex_t IT_2154 = m_b*IT_1325;
    const complex_t IT_2155 = m_b*IT_1084;
    const complex_t IT_2156 = m_s*IT_2032;
    const complex_t IT_2157 = mty::lt::C0iC(18, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0060, IT_0388, IT_0388, mty::lt::reg_int);
    const complex_t IT_2158 = IT_0064*IT_2157;
    const complex_t IT_2159 = m_s*IT_2158;
    const complex_t IT_2160 = -IT_2154 + -IT_2155 + -IT_2156 + -IT_2159;
    const complex_t IT_2161 = IT_2153 + IT_2160;
    const complex_t IT_2162 = IT_1073*IT_2161;
    const complex_t IT_2163 = IT_0332*IT_0671*IT_1039;
    const complex_t IT_2164 = IT_0017*IT_2163;
    const complex_t IT_2165 = IT_1134*IT_2164;
    const complex_t IT_2166 = IT_0193*IT_0209*IT_0321;
    const complex_t IT_2167 = IT_0017*IT_2166;
    const complex_t IT_2168 = IT_0942*IT_2167;
    const complex_t IT_2169 = IT_0321*IT_0918*IT_1346;
    const complex_t IT_2170 = IT_0017*IT_2169;
    const complex_t IT_2171 = IT_0923*IT_2170;
    const complex_t IT_2172 = IT_0209*IT_0229*IT_1093;
    const complex_t IT_2173 = IT_0017*IT_2172;
    const complex_t IT_2174 = IT_1108*IT_2173;
    const complex_t IT_2175 = IT_0537*IT_0548*IT_0902;
    const complex_t IT_2176 = IT_0017*IT_2175;
    const complex_t IT_2177 = IT_0923*IT_2176;
    const complex_t IT_2178 = IT_0321*IT_0416*IT_0694;
    const complex_t IT_2179 = IT_0017*IT_2178;
    const complex_t IT_2180 = IT_1151*IT_2179;
    const complex_t IT_2181 = IT_1700*IT_2118;
    const complex_t IT_2182 = IT_1786*IT_2122;
    const complex_t IT_2183 = IT_1789*IT_1977;
    const complex_t IT_2184 = -IT_0659 + -IT_0662;
    const complex_t IT_2185 = IT_0653 + IT_0655 + IT_0661;
    const complex_t IT_2186 = IT_2184 + IT_2185;
    const complex_t IT_2187 = IT_0651*IT_2186;
    const complex_t IT_2188 = IT_0042*IT_0086*IT_0584;
    const complex_t IT_2189 = IT_0174*IT_2188;
    const complex_t IT_2190 = IT_0090*IT_2189;
    const complex_t IT_2191 = IT_1485*IT_1973;
    const complex_t IT_2192 = IT_1676*IT_2161;
    const complex_t IT_2193 = IT_1679 + IT_1685;
    const complex_t IT_2194 = m_b*IT_1333;
    const complex_t IT_2195 = IT_0064*IT_1683;
    const complex_t IT_2196 = m_s*IT_2195;
    const complex_t IT_2197 = m_b*IT_1687;
    const complex_t IT_2198 = mty::lt::C0iC(18, IT_0008, IT_0008 + IT_0009 + (
      -2)*s_12, IT_0009, IT_0089, IT_0388, IT_0388, mty::lt::reg_int);
    const complex_t IT_2199 = IT_0064*IT_2198;
    const complex_t IT_2200 = m_s*IT_2199;
    const complex_t IT_2201 = -IT_2194 + -IT_2196 + -IT_2197 + -IT_2200;
    const complex_t IT_2202 = IT_2193 + IT_2201;
    const complex_t IT_2203 = IT_1691*IT_2202;
    const complex_t IT_2204 = IT_0086*IT_0141*IT_0561;
    const complex_t IT_2205 = IT_0174*IT_2204;
    const complex_t IT_2206 = IT_0576*IT_2205;
    const complex_t IT_2207 = IT_1187*IT_2147;
    const complex_t IT_2208 = IT_1664*IT_2114;
    const complex_t IT_2209 = IT_1715*IT_2001;
    const complex_t IT_2210 = IT_1649*IT_2011;
    const complex_t IT_2211 = IT_1768*IT_1889;
    const complex_t IT_2212 = IT_0493*IT_0741*IT_1396;
    const complex_t IT_2213 = IT_0174*IT_2212;
    const complex_t IT_2214 = IT_0746*IT_2213;
    const complex_t IT_2215 = IT_0455*IT_0630*IT_0741;
    const complex_t IT_2216 = IT_0174*IT_2215;
    const complex_t IT_2217 = IT_0782*IT_2216;
    const complex_t IT_2218 = IT_0229*IT_0245*IT_0741;
    const complex_t IT_2219 = IT_0174*IT_2218;
    const complex_t IT_2220 = IT_0812*IT_2219;
    const complex_t IT_2221 = IT_0741*IT_0837*IT_1055;
    const complex_t IT_2222 = IT_0174*IT_2221;
    const complex_t IT_2223 = IT_0853*IT_2222;
    const complex_t IT_2224 = IT_0712*IT_0741*IT_1245;
    const complex_t IT_2225 = IT_0174*IT_2224;
    const complex_t IT_2226 = IT_0885*IT_2225;
    const complex_t IT_2227 = IT_0493*IT_0918*IT_1093;
    const complex_t IT_2228 = IT_0017*IT_2227;
    const complex_t IT_2229 = IT_1297*IT_2228;
    const complex_t IT_2230 = IT_0837*IT_1007*IT_1093;
    const complex_t IT_2231 = IT_0017*IT_2230;
    const complex_t IT_2232 = IT_1134*IT_2231;
    const complex_t IT_2233 = IT_1322*IT_2034;
    const complex_t IT_2234 = IT_0269*IT_0902*IT_1104;
    const complex_t IT_2235 = IT_0017*IT_2234;
    const complex_t IT_2236 = IT_0942*IT_2235;
    const complex_t IT_2237 = IT_1655*IT_2038;
    const complex_t IT_2238 = IT_1658*IT_1831;
    const complex_t IT_2239 = IT_1741*IT_2048;
    const complex_t IT_2240 = IT_1744*IT_2059;
    const complex_t IT_2241 = IT_1418*IT_1851;
    const complex_t IT_2242 = IT_1559*IT_2127;
    const complex_t IT_2243 = IT_1703*IT_1986;
    const complex_t IT_2244 = IT_0712*IT_0728*IT_1093;
    const complex_t IT_2245 = IT_0017*IT_2244;
    const complex_t IT_2246 = IT_0732*IT_2245;
    const complex_t IT_2247 = IT_1350 + IT_1354 + IT_1360;
    const complex_t IT_2248 = -IT_1358 + -IT_1361;
    const complex_t IT_2249 = IT_2247 + IT_2248;
    const complex_t IT_2250 = IT_1348*IT_2249;
    const complex_t IT_2251 = IT_1482*IT_2186;
    const complex_t IT_2252 = IT_1577*IT_1969;
    const complex_t IT_2253 = IT_1765*IT_2105;
    const complex_t IT_2254 = IT_1783*IT_2109;
    const complex_t IT_2255 = IT_0042*IT_0056*IT_0902;
    const complex_t IT_2256 = IT_0017*IT_2255;
    const complex_t IT_2257 = IT_0959*IT_2256;
    const complex_t IT_2258 = IT_1738*IT_1846;
    const complex_t IT_2259 = IT_1661*IT_1842;
    const complex_t IT_2260 = IT_1614*IT_2131;
    const complex_t IT_2261 = IT_1617*IT_2094;
    const complex_t IT_2262 = IT_1777*IT_1885;
    const complex_t IT_2263 = IT_1398*IT_2249;
    const complex_t IT_2264 = IT_1623*IT_2067;
    const complex_t IT_2265 = IT_1670*IT_2071;
    const complex_t IT_2266 = IT_0110*IT_0572*IT_1093;
    const complex_t IT_2267 = IT_0017*IT_2266;
    const complex_t IT_2268 = IT_0675*IT_2267;
    const complex_t IT_2269 = IT_0064*IT_1330;
    const complex_t IT_2270 = -IT_2195 + -IT_2269;
    const complex_t IT_2271 = IT_1331 + IT_2270;
    const complex_t IT_2272 = IT_1735*IT_2271;
    const complex_t IT_2273 = IT_1718*IT_2202;
    const complex_t IT_2274 = IT_1329*IT_2271;
    const complex_t IT_2275 = IT_1809 + -IT_1812 + -IT_1815 + -IT_1818 +
       IT_1821 + 0.5*IT_1832 + 0.5*IT_1843 + 0.5*IT_1847 + 0.5*IT_1852 + IT_1856
       + IT_1860 + IT_1864 + 0.5*IT_1875 + 0.5*IT_1886 + 0.5*IT_1890 + -IT_1893 
      + -IT_1896 + -IT_1899 + -IT_1902 + -IT_1905 + -IT_1909 + -IT_1913 +
       IT_1914 + IT_1918 + -IT_1919 + -IT_1920 + 0.5*IT_1925 + IT_1928 + IT_1931
       + IT_1934 + IT_1937 + IT_1940 + 0.5*IT_1944 + 0.5*IT_1945 + -IT_1948 + 
      -IT_1951 + -IT_1952 + IT_1955 + -IT_1958 + -IT_1962 + -IT_1966 + -IT_1970 
      + -IT_1974 + IT_1978 + IT_1979 + IT_1983 + IT_1987 + 0.5*IT_1992 + 0.5
      *IT_2002 + 0.5*IT_2012 + -IT_2015 + -IT_2018 + 0.5*IT_2019 + 0.5*IT_2020 +
       IT_2023 + -IT_2026 + IT_2029 + 0.5*IT_2030 + 0.5*IT_2035 + 0.5*IT_2039 +
       0.5*IT_2049 + 0.5*IT_2060 + IT_2063 + IT_2064 + IT_2068 + IT_2072 +
       IT_2076 + IT_2080 + IT_2084 + 0.5*IT_2095 + 0.5*IT_2106 + IT_2110 +
       IT_2111 + IT_2115 + IT_2119 + IT_2123 + -IT_2124 + 0.5*IT_2128 + 0.5
      *IT_2132 + -IT_2135 + -IT_2138 + -IT_2141 + -IT_2144 + 0.5*IT_2148 + 
      -IT_2149 + -IT_2150 + -IT_2151 + -IT_2152 + 0.5*IT_2162 + -IT_2165 +
       IT_2168 + IT_2171 + IT_2174 + -IT_2177 + IT_2180 + -IT_2181 + -IT_2182 + 
      -IT_2183 + -IT_2187 + -IT_2190 + IT_2191 + 0.5*IT_2192 + 0.5*IT_2203 + 
      -IT_2206 + 0.5*IT_2207 + -IT_2208 + 0.5*IT_2209 + 0.5*IT_2210 + 0.5
      *IT_2211 + IT_2214 + IT_2217 + IT_2220 + IT_2223 + IT_2226 + IT_2229 +
       IT_2232 + 0.5*IT_2233 + -IT_2236 + 0.5*IT_2237 + 0.5*IT_2238 + 0.5
      *IT_2239 + 0.5*IT_2240 + 0.5*IT_2241 + 0.5*IT_2242 + -IT_2243 + IT_2246 +
       IT_2250 + IT_2251 + IT_2252 + 0.5*IT_2253 + -IT_2254 + -IT_2257 + 0.5
      *IT_2258 + 0.5*IT_2259 + 0.5*IT_2260 + 0.5*IT_2261 + 0.5*IT_2262 + 
      -IT_2263 + -IT_2264 + -IT_2265 + IT_2268 + 0.5*IT_2272 + 0.5*IT_2273 + 0.5
      *IT_2274;
    const complex_t IT_2276 = IT_1804*IT_2275;
    const complex_t IT_2277 = -IT_2276;
    const complex_t IT_2278 = (-2)*IT_2277;
    const complex_t IT_2279 = IT_0016*IT_2278;
    const complex_t IT_2280 = -IT_2279;
    const complex_t IT_2281 = powq(m_b + m_s, -1);
    const complex_t IT_2282 = IT_1624 + IT_1625 + IT_1627 + IT_1629 + IT_1631 
      + IT_1632 + IT_1635;
    const complex_t IT_2283 = 2*IT_1634;
    const complex_t IT_2284 = IT_2282 + IT_2283;
    const complex_t IT_2285 = IT_1747*IT_2284;
    const complex_t IT_2286 = IT_1402 + IT_1403 + IT_1405 + IT_1407 + IT_1409 
      + IT_1410 + IT_1412;
    const complex_t IT_2287 = 2*IT_1413;
    const complex_t IT_2288 = IT_2286 + IT_2287;
    const complex_t IT_2289 = IT_1522*IT_2288;
    const complex_t IT_2290 = IT_0922 + IT_0924 + IT_0926 + IT_0928 + IT_0930 
      + IT_0933 + IT_0934;
    const complex_t IT_2291 = 2*IT_0932;
    const complex_t IT_2292 = IT_2290 + IT_2291;
    const complex_t IT_2293 = IT_0920*IT_2292;
    const complex_t IT_2294 = IT_0165 + IT_0168;
    const complex_t IT_2295 = -IT_0162 + -IT_0171;
    const complex_t IT_2296 = IT_2294 + IT_2295;
    const complex_t IT_2297 = IT_1661*IT_2296;
    const complex_t IT_2298 = IT_1350 + IT_1352 + IT_1354 + IT_1356 + IT_1358 
      + IT_1361 + IT_1362;
    const complex_t IT_2299 = 2*IT_1360;
    const complex_t IT_2300 = IT_2298 + IT_2299;
    const complex_t IT_2301 = IT_1348*IT_2300;
    const complex_t IT_2302 = IT_0652 + IT_0653 + IT_0655 + IT_0657 + IT_0659 
      + IT_0662 + IT_0663;
    const complex_t IT_2303 = 2*IT_0661;
    const complex_t IT_2304 = IT_2302 + IT_2303;
    const complex_t IT_2305 = IT_1482*IT_2304;
    const complex_t IT_2306 = IT_0634 + IT_0636 + IT_0638 + IT_0640 + IT_0643 
      + IT_0644 + IT_0645;
    const complex_t IT_2307 = 2*IT_0646;
    const complex_t IT_2308 = IT_2306 + IT_2307;
    const complex_t IT_2309 = IT_1577*IT_2308;
    const complex_t IT_2310 = IT_0307 + IT_0310;
    const complex_t IT_2311 = -IT_0305 + -IT_0313;
    const complex_t IT_2312 = IT_2310 + IT_2311;
    const complex_t IT_2313 = IT_0302*IT_2312;
    const complex_t IT_2314 = IT_1581 + IT_1583 + IT_1585 + IT_1586 + IT_1588 
      + IT_1591 + IT_1592;
    const complex_t IT_2315 = 2*IT_1590;
    const complex_t IT_2316 = IT_2314 + IT_2315;
    const complex_t IT_2317 = IT_1580*IT_2316;
    const complex_t IT_2318 = 2*IT_0962;
    const complex_t IT_2319 = IT_0960 + IT_0964 + IT_0965 + IT_0967 + IT_0969 
      + IT_0971 + IT_0972;
    const complex_t IT_2320 = IT_2318 + IT_2319;
    const complex_t IT_2321 = IT_1762*IT_2320;
    const complex_t IT_2322 = IT_1439 + IT_1440 + IT_1442 + IT_1444 + IT_1446 
      + IT_1447 + IT_1450;
    const complex_t IT_2323 = 2*IT_1449;
    const complex_t IT_2324 = IT_2322 + IT_2323;
    const complex_t IT_2325 = IT_1438*IT_2324;
    const complex_t IT_2326 = IT_0941 + IT_0943 + IT_0945 + IT_0947 + IT_0949 
      + IT_0952 + IT_0953;
    const complex_t IT_2327 = 2*IT_0951;
    const complex_t IT_2328 = IT_2326 + IT_2327;
    const complex_t IT_2329 = IT_1455*IT_2328;
    const complex_t IT_2330 = IT_0390 + IT_0392 + IT_0394 + IT_0396 + IT_0400 
      + IT_0401 + IT_0402;
    const complex_t IT_2331 = 2*IT_0399;
    const complex_t IT_2332 = IT_2330 + IT_2331;
    const complex_t IT_2333 = IT_0387*IT_2332;
    const complex_t IT_2334 = IT_0433 + IT_0439;
    const complex_t IT_2335 = -IT_0436 + -IT_0442;
    const complex_t IT_2336 = IT_2334 + IT_2335;
    const complex_t IT_2337 = IT_0429*IT_2336;
    const complex_t IT_2338 = IT_1369 + IT_1372;
    const complex_t IT_2339 = -IT_1375 + -IT_1378;
    const complex_t IT_2340 = IT_2338 + IT_2339;
    const complex_t IT_2341 = IT_1367*IT_2340;
    const complex_t IT_2342 = IT_0510 + IT_0513;
    const complex_t IT_2343 = -IT_0516 + -IT_0519;
    const complex_t IT_2344 = IT_2342 + IT_2343;
    const complex_t IT_2345 = IT_0506*IT_2344;
    const complex_t IT_2346 = IT_0348 + IT_0350 + IT_0352 + IT_0354 + IT_0356 
      + IT_0357 + IT_0360;
    const complex_t IT_2347 = 2*IT_0359;
    const complex_t IT_2348 = IT_2346 + IT_2347;
    const complex_t IT_2349 = IT_0345*IT_2348;
    const complex_t IT_2350 = IT_0117 + IT_0120;
    const complex_t IT_2351 = -IT_0115 + -IT_0123;
    const complex_t IT_2352 = IT_2350 + IT_2351;
    const complex_t IT_2353 = IT_0112*IT_2352;
    const complex_t IT_2354 = IT_0159*IT_2296;
    const complex_t IT_2355 = IT_0273 + IT_0276;
    const complex_t IT_2356 = -IT_0279 + -IT_0282;
    const complex_t IT_2357 = IT_2355 + IT_2356;
    const complex_t IT_2358 = IT_0271*IT_2357;
    const complex_t IT_2359 = IT_1598 + IT_1602 + IT_1603 + IT_1605 + IT_1606 
      + IT_1607 + IT_1609;
    const complex_t IT_2360 = 2*IT_1600;
    const complex_t IT_2361 = IT_2359 + IT_2360;
    const complex_t IT_2362 = IT_1597*IT_2361;
    const complex_t IT_2363 = IT_1759*IT_2292;
    const complex_t IT_2364 = IT_0471 + IT_0474;
    const complex_t IT_2365 = -IT_0477 + -IT_0480;
    const complex_t IT_2366 = IT_2364 + IT_2365;
    const complex_t IT_2367 = IT_0468*IT_2366;
    const complex_t IT_2368 = 2*IT_0757;
    const complex_t IT_2369 = IT_0745 + IT_0747 + IT_0749 + IT_0751 + IT_0753 
      + IT_0754 + IT_0756;
    const complex_t IT_2370 = IT_2368 + IT_2369;
    const complex_t IT_2371 = IT_1783*IT_2370;
    const complex_t IT_2372 = 2*IT_1431;
    const complex_t IT_2373 = IT_1422 + IT_1423 + IT_1425 + IT_1427 + IT_1430 
      + IT_1432 + IT_1433;
    const complex_t IT_2374 = IT_2372 + IT_2373;
    const complex_t IT_2375 = IT_1421*IT_2374;
    const complex_t IT_2376 = IT_1401*IT_2288;
    const complex_t IT_2377 = IT_1489 + IT_1490 + IT_1492 + IT_1494 + IT_1496 
      + IT_1497 + IT_1499;
    const complex_t IT_2378 = 2*IT_1500;
    const complex_t IT_2379 = IT_2377 + IT_2378;
    const complex_t IT_2380 = IT_1571*IT_2379;
    const complex_t IT_2381 = IT_0743*IT_2370;
    const complex_t IT_2382 = IT_0763 + IT_0764 + IT_0766 + IT_0768 + IT_0770 
      + IT_0771 + IT_0773;
    const complex_t IT_2383 = 2*IT_0774;
    const complex_t IT_2384 = IT_2382 + IT_2383;
    const complex_t IT_2385 = IT_0762*IT_2384;
    const complex_t IT_2386 = 2*IT_0793;
    const complex_t IT_2387 = IT_0781 + IT_0783 + IT_0785 + IT_0787 + IT_0789 
      + IT_0790 + IT_0792;
    const complex_t IT_2388 = IT_2386 + IT_2387;
    const complex_t IT_2389 = IT_0779*IT_2388;
    const complex_t IT_2390 = 2*IT_0823;
    const complex_t IT_2391 = IT_0811 + IT_0813 + IT_0815 + IT_0817 + IT_0819 
      + IT_0820 + IT_0822;
    const complex_t IT_2392 = IT_2390 + IT_2391;
    const complex_t IT_2393 = IT_0809*IT_2392;
    const complex_t IT_2394 = 2*IT_0864;
    const complex_t IT_2395 = IT_0852 + IT_0854 + IT_0856 + IT_0858 + IT_0860 
      + IT_0861 + IT_0863;
    const complex_t IT_2396 = IT_2394 + IT_2395;
    const complex_t IT_2397 = IT_0850*IT_2396;
    const complex_t IT_2398 = 2*IT_0891;
    const complex_t IT_2399 = IT_0882 + IT_0884 + IT_0886 + IT_0888 + IT_0892 
      + IT_0893 + IT_0894;
    const complex_t IT_2400 = IT_2398 + IT_2399;
    const complex_t IT_2401 = IT_0880*IT_2400;
    const complex_t IT_2402 = IT_1485*IT_2374;
    const complex_t IT_2403 = IT_1488*IT_2379;
    const complex_t IT_2404 = IT_1507 + IT_1508 + IT_1509 + IT_1511 + IT_1513 
      + IT_1514 + IT_1516;
    const complex_t IT_2405 = 2*IT_1517;
    const complex_t IT_2406 = IT_2404 + IT_2405;
    const complex_t IT_2407 = IT_1505*IT_2406;
    const complex_t IT_2408 = IT_1526 + IT_1527 + IT_1529 + IT_1531 + IT_1533 
      + IT_1534 + IT_1536;
    const complex_t IT_2409 = 2*IT_1537;
    const complex_t IT_2410 = IT_2408 + IT_2409;
    const complex_t IT_2411 = IT_1525*IT_2410;
    const complex_t IT_2412 = IT_1544 + IT_1545 + IT_1546 + IT_1548 + IT_1551 
      + IT_1553 + IT_1554;
    const complex_t IT_2413 = 2*IT_1552;
    const complex_t IT_2414 = IT_2412 + IT_2413;
    const complex_t IT_2415 = IT_1542*IT_2414;
    const complex_t IT_2416 = IT_0958*IT_2320;
    const complex_t IT_2417 = IT_1712*IT_2332;
    const complex_t IT_2418 = IT_1015 + IT_1018;
    const complex_t IT_2419 = -IT_1012 + -IT_1021;
    const complex_t IT_2420 = IT_2418 + IT_2419;
    const complex_t IT_2421 = IT_1009*IT_2420;
    const complex_t IT_2422 = IT_1063 + IT_1066;
    const complex_t IT_2423 = -IT_1060 + -IT_1069;
    const complex_t IT_2424 = IT_2422 + IT_2423;
    const complex_t IT_2425 = IT_1057*IT_2424;
    const complex_t IT_2426 = IT_1658*IT_2352;
    const complex_t IT_2427 = IT_1744*IT_2312;
    const complex_t IT_2428 = IT_1700*IT_2392;
    const complex_t IT_2429 = IT_1664*IT_2388;
    const complex_t IT_2430 = IT_1623*IT_2284;
    const complex_t IT_2431 = IT_1670*IT_2361;
    const complex_t IT_2432 = IT_1703*IT_2414;
    const complex_t IT_2433 = IT_0939*IT_2328;
    const complex_t IT_2434 = IT_1646*IT_2324;
    const complex_t IT_2435 = IT_1709*IT_2348;
    const complex_t IT_2436 = IT_1715*IT_2420;
    const complex_t IT_2437 = IT_1076 + IT_1079;
    const complex_t IT_2438 = -IT_1082 + -IT_1085;
    const complex_t IT_2439 = IT_2437 + IT_2438;
    const complex_t IT_2440 = IT_1676*IT_2439;
    const complex_t IT_2441 = IT_1465*IT_2384;
    const complex_t IT_2442 = IT_1786*IT_2396;
    const complex_t IT_2443 = IT_1649*IT_2424;
    const complex_t IT_2444 = IT_1073*IT_2439;
    const complex_t IT_2445 = IT_1741*IT_2357;
    const complex_t IT_2446 = IT_1617*IT_2336;
    const complex_t IT_2447 = IT_1765*IT_2366;
    const complex_t IT_2448 = IT_1774*IT_2340;
    const complex_t IT_2449 = IT_1777*IT_2344;
    const complex_t IT_2450 = IT_1789*IT_2400;
    const complex_t IT_2451 = IT_1398*IT_2300;
    const complex_t IT_2452 = IT_1667*IT_2316;
    const complex_t IT_2453 = IT_0632*IT_2308;
    const complex_t IT_2454 = IT_0651*IT_2304;
    const complex_t IT_2455 = IT_1568*IT_2406;
    const complex_t IT_2456 = IT_1798*IT_2410;
    const complex_t IT_2457 = IT_1679 + IT_1682;
    const complex_t IT_2458 = -IT_1685 + -IT_1688;
    const complex_t IT_2459 = IT_2457 + IT_2458;
    const complex_t IT_2460 = IT_1718*IT_2459;
    const complex_t IT_2461 = IT_1691*IT_2459;
    const complex_t IT_2462 = IT_0068 + 2*IT_0094 + IT_0181 + -IT_0218 + 
      -IT_0253 + IT_0556 + 2*IT_0579 + 2*IT_0598 + 2*IT_0678 + 2*IT_0701 + 2
      *IT_0735 + 2*IT_1111 + 2*IT_1130 + 2*IT_1137 + 2*IT_1142 + 2*IT_1147 + 2
      *IT_1154 + 2*IT_1159 + 2*IT_1164 + 2*IT_1185 + -IT_1192 + IT_1215 + 2
      *IT_1222 + 2*IT_1229 + 2*IT_1252 + 2*IT_1257 + 2*IT_1262 + 2*IT_1267 + 2
      *IT_1272 + 2*IT_1293 + 2*IT_1300 + 2*IT_1303 + 2*IT_1306 + -IT_1313 + 
      -IT_1320 + -IT_1327 + IT_1335 + IT_1419 + IT_1463 + 2*IT_1469 + IT_1472 +
       IT_1475 + 2*IT_1480 + IT_1560 + 2*IT_1563 + 2*IT_1566 + 2*IT_1575 + 
      -IT_1615 + IT_1621 + 2*IT_1641 + 2*IT_1644 + 2*IT_1653 + -IT_1656 + 2
      *IT_1674 + 2*IT_1695 + 2*IT_1698 + 2*IT_1707 + IT_1722 + 2*IT_1725 + 2
      *IT_1730 + 2*IT_1733 + -IT_1736 + -IT_1739 + 2*IT_1751 + 2*IT_1754 + 2
      *IT_1757 + -IT_1769 + -IT_1772 + 2*IT_1781 + 2*IT_1793 + 2*IT_1796 + 2
      *IT_1802 + 2*IT_2285 + 2*IT_2289 + 2*IT_2293 + IT_2297 + 2*IT_2301 + 2
      *IT_2305 + 2*IT_2309 + -IT_2313 + 2*IT_2317 + 2*IT_2321 + 2*IT_2325 + 2
      *IT_2329 + 2*IT_2333 + IT_2337 + IT_2341 + IT_2345 + 2*IT_2349 + -IT_2353 
      + -IT_2354 + -IT_2358 + 2*IT_2362 + 2*IT_2363 + IT_2367 + 2*IT_2371 + 2
      *IT_2375 + 2*IT_2376 + 2*IT_2380 + 2*IT_2381 + 2*IT_2385 + 2*IT_2389 + 2
      *IT_2393 + 2*IT_2397 + 2*IT_2401 + 2*IT_2402 + 2*IT_2403 + 2*IT_2407 + 2
      *IT_2411 + 2*IT_2415 + 2*IT_2416 + 2*IT_2417 + -IT_2421 + -IT_2425 +
       IT_2426 + IT_2427 + 2*IT_2428 + 2*IT_2429 + 2*IT_2430 + 2*IT_2431 + 2
      *IT_2432 + 2*IT_2433 + 2*IT_2434 + 2*IT_2435 + IT_2436 + -IT_2440 + 2
      *IT_2441 + 2*IT_2442 + IT_2443 + IT_2444 + IT_2445 + -IT_2446 + -IT_2447 +
       -IT_2448 + -IT_2449 + 2*IT_2450 + 2*IT_2451 + 2*IT_2452 + 2*IT_2453 + 2
      *IT_2454 + 2*IT_2455 + 2*IT_2456 + IT_2460 + -IT_2461;
    const complex_t IT_2463 = IT_1804*IT_2462;
    const complex_t IT_2464 = 0.5*IT_2463;
    const complex_t IT_2465 = 2*IT_2464;
    const complex_t IT_2466 = IT_0016*IT_2465;
    const complex_t IT_2467 = -IT_2466;
    const complex_t IT_2468 = IT_0117 + IT_1826 + IT_1829;
    const complex_t IT_2469 = -IT_0123 + -IT_1823 + -IT_1824;
    const complex_t IT_2470 = IT_2468 + IT_2469;
    const complex_t IT_2471 = IT_0112*IT_2470;
    const complex_t IT_2472 = IT_0165 + IT_1838 + IT_1840;
    const complex_t IT_2473 = -IT_0162 + -IT_1834 + -IT_1835;
    const complex_t IT_2474 = IT_2472 + IT_2473;
    const complex_t IT_2475 = IT_0159*IT_2474;
    const complex_t IT_2476 = IT_0273 + IT_2045 + IT_2046;
    const complex_t IT_2477 = -IT_0279 + -IT_2041 + -IT_2042;
    const complex_t IT_2478 = IT_2476 + IT_2477;
    const complex_t IT_2479 = IT_0271*IT_2478;
    const complex_t IT_2480 = IT_0307 + IT_2055 + IT_2057;
    const complex_t IT_2481 = -IT_0305 + -IT_2051 + -IT_2052;
    const complex_t IT_2482 = IT_2480 + IT_2481;
    const complex_t IT_2483 = IT_0302*IT_2482;
    const complex_t IT_2484 = IT_1583 + IT_1586 + IT_1588 + IT_1590 + IT_1591;
    const complex_t IT_2485 = IT_1580*IT_2484;
    const complex_t IT_2486 = IT_1598 + IT_1600 + IT_1602 + IT_1605 + IT_1606;
    const complex_t IT_2487 = IT_1597*IT_2486;
    const complex_t IT_2488 = IT_0924 + IT_0926 + IT_0930 + IT_0932 + IT_0933;
    const complex_t IT_2489 = IT_1759*IT_2488;
    const complex_t IT_2490 = IT_0960 + IT_0962 + IT_0964 + IT_0967 + IT_0972;
    const complex_t IT_2491 = IT_1762*IT_2490;
    const complex_t IT_2492 = IT_1440 + IT_1442 + IT_1446 + IT_1447 + IT_1449;
    const complex_t IT_2493 = IT_1438*IT_2492;
    const complex_t IT_2494 = IT_0943 + IT_0945 + IT_0949 + IT_0951 + IT_0952;
    const complex_t IT_2495 = IT_1455*IT_2494;
    const complex_t IT_2496 = IT_0352 + IT_0354 + IT_0356 + IT_0357 + IT_0359;
    const complex_t IT_2497 = IT_0345*IT_2496;
    const complex_t IT_2498 = IT_0392 + IT_0394 + IT_0399 + IT_0400 + IT_0402;
    const complex_t IT_2499 = IT_0387*IT_2498;
    const complex_t IT_2500 = IT_0433 + IT_2089 + IT_2092;
    const complex_t IT_2501 = -IT_0436 + -IT_2086 + -IT_2087;
    const complex_t IT_2502 = IT_2500 + IT_2501;
    const complex_t IT_2503 = IT_0429*IT_2502;
    const complex_t IT_2504 = IT_0471 + IT_2100 + IT_2103;
    const complex_t IT_2505 = -IT_0477 + -IT_2097 + -IT_2098;
    const complex_t IT_2506 = IT_2504 + IT_2505;
    const complex_t IT_2507 = IT_0468*IT_2506;
    const complex_t IT_2508 = IT_1423 + IT_1425 + IT_1431 + IT_1432 + IT_1433;
    const complex_t IT_2509 = IT_1421*IT_2508;
    const complex_t IT_2510 = IT_1402 + IT_1405 + IT_1409 + IT_1410 + IT_1413;
    const complex_t IT_2511 = IT_1401*IT_2510;
    const complex_t IT_2512 = IT_1489 + IT_1492 + IT_1496 + IT_1497 + IT_1500;
    const complex_t IT_2513 = IT_1571*IT_2512;
    const complex_t IT_2514 = IT_1546 + IT_1548 + IT_1552 + IT_1553 + IT_1554;
    const complex_t IT_2515 = IT_1703*IT_2514;
    const complex_t IT_2516 = IT_0747 + IT_0749 + IT_0753 + IT_0754 + IT_0757;
    const complex_t IT_2517 = IT_0743*IT_2516;
    const complex_t IT_2518 = IT_0763 + IT_0766 + IT_0770 + IT_0771 + IT_0774;
    const complex_t IT_2519 = IT_0762*IT_2518;
    const complex_t IT_2520 = IT_0783 + IT_0785 + IT_0789 + IT_0790 + IT_0793;
    const complex_t IT_2521 = IT_0779*IT_2520;
    const complex_t IT_2522 = IT_0813 + IT_0815 + IT_0819 + IT_0820 + IT_0823;
    const complex_t IT_2523 = IT_0809*IT_2522;
    const complex_t IT_2524 = IT_0854 + IT_0856 + IT_0860 + IT_0861 + IT_0864;
    const complex_t IT_2525 = IT_0850*IT_2524;
    const complex_t IT_2526 = IT_0886 + IT_0888 + IT_0891 + IT_0892 + IT_0893;
    const complex_t IT_2527 = IT_0880*IT_2526;
    const complex_t IT_2528 = IT_1485*IT_2508;
    const complex_t IT_2529 = IT_1488*IT_2512;
    const complex_t IT_2530 = IT_1507 + IT_1509 + IT_1513 + IT_1514 + IT_1517;
    const complex_t IT_2531 = IT_1505*IT_2530;
    const complex_t IT_2532 = IT_1522*IT_2510;
    const complex_t IT_2533 = IT_1526 + IT_1529 + IT_1533 + IT_1534 + IT_1537;
    const complex_t IT_2534 = IT_1525*IT_2533;
    const complex_t IT_2535 = IT_1542*IT_2514;
    const complex_t IT_2536 = IT_1712*IT_2498;
    const complex_t IT_2537 = IT_1015 + IT_1996 + IT_1999;
    const complex_t IT_2538 = -IT_1012 + -IT_1997 + -IT_1998;
    const complex_t IT_2539 = IT_2537 + IT_2538;
    const complex_t IT_2540 = IT_1009*IT_2539;
    const complex_t IT_2541 = IT_1063 + IT_2005 + IT_2008;
    const complex_t IT_2542 = -IT_1060 + -IT_2004 + -IT_2009;
    const complex_t IT_2543 = IT_2541 + IT_2542;
    const complex_t IT_2544 = IT_1057*IT_2543;
    const complex_t IT_2545 = IT_1568*IT_2530;
    const complex_t IT_2546 = IT_1715*IT_2539;
    const complex_t IT_2547 = IT_1649*IT_2543;
    const complex_t IT_2548 = IT_1076 + IT_2156 + IT_2159;
    const complex_t IT_2549 = -IT_1082 + -IT_2154 + -IT_2155;
    const complex_t IT_2550 = IT_2548 + IT_2549;
    const complex_t IT_2551 = IT_1073*IT_2550;
    const complex_t IT_2552 = IT_1679 + IT_2196 + IT_2200;
    const complex_t IT_2553 = -IT_1685 + -IT_2194 + -IT_2197;
    const complex_t IT_2554 = IT_2552 + IT_2553;
    const complex_t IT_2555 = IT_1718*IT_2554;
    const complex_t IT_2556 = IT_1350 + IT_1354 + IT_1358 + IT_1360 + IT_1361;
    const complex_t IT_2557 = IT_1348*IT_2556;
    const complex_t IT_2558 = IT_0653 + IT_0655 + IT_0659 + IT_0661 + IT_0662;
    const complex_t IT_2559 = IT_1482*IT_2558;
    const complex_t IT_2560 = IT_0636 + IT_0638 + IT_0644 + IT_0645 + IT_0646;
    const complex_t IT_2561 = IT_1577*IT_2560;
    const complex_t IT_2562 = IT_1624 + IT_1627 + IT_1631 + IT_1632 + IT_1634;
    const complex_t IT_2563 = IT_1747*IT_2562;
    const complex_t IT_2564 = IT_1765*IT_2506;
    const complex_t IT_2565 = IT_1369 + IT_1869 + IT_1872;
    const complex_t IT_2566 = -IT_1375 + -IT_1866 + -IT_1867;
    const complex_t IT_2567 = IT_2565 + IT_2566;
    const complex_t IT_2568 = IT_1367*IT_2567;
    const complex_t IT_2569 = IT_0510 + IT_1880 + IT_1883;
    const complex_t IT_2570 = -IT_0516 + -IT_1877 + -IT_1878;
    const complex_t IT_2571 = IT_2569 + IT_2570;
    const complex_t IT_2572 = IT_0506*IT_2571;
    const complex_t IT_2573 = IT_1774*IT_2567;
    const complex_t IT_2574 = IT_1783*IT_2516;
    const complex_t IT_2575 = IT_1664*IT_2520;
    const complex_t IT_2576 = IT_1465*IT_2518;
    const complex_t IT_2577 = IT_1670*IT_2486;
    const complex_t IT_2578 = IT_1798*IT_2533;
    const complex_t IT_2579 = IT_1646*IT_2492;
    const complex_t IT_2580 = IT_0958*IT_2490;
    const complex_t IT_2581 = IT_1709*IT_2496;
    const complex_t IT_2582 = IT_1658*IT_2470;
    const complex_t IT_2583 = IT_1661*IT_2474;
    const complex_t IT_2584 = IT_1741*IT_2478;
    const complex_t IT_2585 = IT_1744*IT_2482;
    const complex_t IT_2586 = IT_1617*IT_2502;
    const complex_t IT_2587 = IT_1777*IT_2571;
    const complex_t IT_2588 = IT_1700*IT_2522;
    const complex_t IT_2589 = IT_1786*IT_2524;
    const complex_t IT_2590 = IT_1789*IT_2526;
    const complex_t IT_2591 = IT_1398*IT_2556;
    const complex_t IT_2592 = IT_1667*IT_2484;
    const complex_t IT_2593 = IT_0632*IT_2560;
    const complex_t IT_2594 = IT_0651*IT_2558;
    const complex_t IT_2595 = IT_1623*IT_2562;
    const complex_t IT_2596 = IT_0920*IT_2488;
    const complex_t IT_2597 = IT_0939*IT_2494;
    const complex_t IT_2598 = IT_1676*IT_2550;
    const complex_t IT_2599 = IT_1691*IT_2554;
    const complex_t IT_2600 = IT_1809 + IT_1812 + IT_1815 + IT_1818 + IT_1821 
      + (-0.5)*IT_1847 + 0.5*IT_1852 + (-0.5)*IT_1890 + IT_1893 + IT_1896 +
       IT_1899 + IT_1902 + IT_1905 + 0.5*IT_1925 + IT_1928 + IT_1931 + IT_1934 +
       IT_1937 + IT_1940 + 0.5*IT_1944 + (-0.5)*IT_1945 + IT_1948 + IT_1951 +
       IT_1955 + IT_1958 + 0.5*IT_1992 + IT_2015 + IT_2018 + (-0.5)*IT_2019 + (
      -0.5)*IT_2020 + IT_2023 + IT_2026 + IT_2029 + (-0.5)*IT_2035 + (-0.5)
      *IT_2039 + IT_2063 + 0.5*IT_2128 + (-0.5)*IT_2132 + IT_2135 + IT_2138 +
       IT_2141 + IT_2144 + (-0.5)*IT_2148 + IT_2165 + IT_2168 + IT_2171 +
       IT_2174 + IT_2177 + IT_2180 + IT_2190 + IT_2206 + 0.5*IT_2207 + 0.5
      *IT_2211 + IT_2214 + IT_2217 + IT_2220 + IT_2223 + IT_2226 + IT_2229 +
       IT_2232 + 0.5*IT_2233 + IT_2236 + 0.5*IT_2237 + (-0.5)*IT_2241 + (-0.5)
      *IT_2242 + IT_2246 + IT_2257 + 0.5*IT_2258 + 0.5*IT_2260 + IT_2268 + 0.5
      *IT_2272 + (-0.5)*IT_2274 + 0.5*IT_2471 + 0.5*IT_2475 + 0.5*IT_2479 + 0.5
      *IT_2483 + -IT_2485 + -IT_2487 + -IT_2489 + -IT_2491 + -IT_2493 + -IT_2495
       + -IT_2497 + -IT_2499 + (-0.5)*IT_2503 + (-0.5)*IT_2507 + -IT_2509 + 
      -IT_2511 + -IT_2513 + -IT_2515 + -IT_2517 + -IT_2519 + -IT_2521 + -IT_2523
       + -IT_2525 + -IT_2527 + -IT_2528 + -IT_2529 + -IT_2531 + -IT_2532 + 
      -IT_2534 + -IT_2535 + -IT_2536 + 0.5*IT_2540 + 0.5*IT_2544 + -IT_2545 + (
      -0.5)*IT_2546 + (-0.5)*IT_2547 + (-0.5)*IT_2551 + (-0.5)*IT_2555 + 
      -IT_2557 + -IT_2559 + -IT_2561 + -IT_2563 + 0.5*IT_2564 + (-0.5)*IT_2568 +
       (-0.5)*IT_2572 + 0.5*IT_2573 + -IT_2574 + -IT_2575 + -IT_2576 + -IT_2577 
      + -IT_2578 + -IT_2579 + -IT_2580 + -IT_2581 + (-0.5)*IT_2582 + (-0.5)
      *IT_2583 + (-0.5)*IT_2584 + (-0.5)*IT_2585 + 0.5*IT_2586 + 0.5*IT_2587 + 
      -IT_2588 + -IT_2589 + -IT_2590 + -IT_2591 + -IT_2592 + -IT_2593 + -IT_2594
       + -IT_2595 + -IT_2596 + -IT_2597 + 0.5*IT_2598 + 0.5*IT_2599;
    const complex_t IT_2601 = IT_1804*IT_2600;
    const complex_t IT_2602 = IT_0016*IT_2601;
    const complex_t IT_2603 = (-2)*IT_2602;
    const complex_t IT_2604 = (-0.5)*IT_0000*IT_0007*IT_0011*(IT_1806 + 
      -IT_2280)*(s_12 + -1./2*IT_0008 + -1./2*IT_0009 + -1./2*reg_prop) + (-0.5)
      *IT_0007*IT_0011*IT_2281*(IT_2467 + -IT_2603)*(s_12 + -1./2*IT_0008 + -1.
      /2*IT_0009 + -1./2*reg_prop);
    return IT_2604;
}
} // End of namespace c9_nmfv
