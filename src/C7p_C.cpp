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
#include "C7p_C.h"
#include "common.h"

#include "params.h"
#include "libcomplexop.h"

namespace c9_nmfv {

complex_t C7p_C(
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
    const real_t s_12 = param.s_12;
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
    const complex_t IT_0000 = powq(M_W, 2);
    const complex_t IT_0001 = powq(m_b, -1);
    const complex_t IT_0002 = powq(V_tb, -1);
    const complex_t IT_0003 = cpowq(conjq(V_ts), -1);
    const complex_t IT_0004 = powq(e_em, -3);
    const complex_t IT_0005 = 9.86960440108936*IT_0000*IT_0001*IT_0002*IT_0003
      *IT_0004;
    const complex_t IT_0006 = cosq(theta_W);
    const complex_t IT_0007 = cpowq(IT_0006, -1);
    const complex_t IT_0008 = tanq(theta_W);
    const complex_t IT_0009 = cpowq(IT_0008, 2);
    const complex_t IT_0010 = cpowq(1 + IT_0009, (-0.5));
    const complex_t IT_0011 = (complex_t{0, 1})*e_em*IT_0007*IT_0010;
    const complex_t IT_0012 = (-0.666666666666667)*IT_0011;
    const complex_t IT_0013 = powq(m_b, 2);
    const complex_t IT_0014 = powq(m_s, 2);
    const complex_t IT_0015 = powq(m_C2p, 2);
    const complex_t IT_0016 = powq(m_st_L, 2);
    const complex_t IT_0017 = mty::lt::C0iC(0, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0015, IT_0016, IT_0016, mty::lt::reg_int);
    const complex_t IT_0018 = IT_0012*IT_0017;
    const complex_t IT_0019 = (-1.33333333333333)*IT_0011;
    const complex_t IT_0020 = mty::lt::C0iC(3, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0015, IT_0016, IT_0016, mty::lt::reg_int);
    const complex_t IT_0021 = IT_0019*IT_0020;
    const complex_t IT_0022 = IT_0018 + IT_0021;
    const complex_t IT_0023 = 0.101321183642338*m_C2p;
    const complex_t IT_0024 = sinq(theta_W);
    const complex_t IT_0025 = cpowq(IT_0024, -1);
    const complex_t IT_0026 = V_us*e_em*conjq(V_Wp2)*U_su_02;
    const complex_t IT_0027 = IT_0025*IT_0026;
    const complex_t IT_0028 = conjq(V_cs)*e_em*conjq(V_Wp2)*U_su_12;
    const complex_t IT_0029 = IT_0025*IT_0028;
    const complex_t IT_0030 = conjq(V_ts)*e_em*conjq(V_Wp2)*U_su_22;
    const complex_t IT_0031 = IT_0025*IT_0030;
    const complex_t IT_0032 = sinq(beta);
    const complex_t IT_0033 = cpowq(IT_0032, -1);
    const complex_t IT_0034 = IT_0025*IT_0033;
    const complex_t IT_0035 = powq(M_W, -1);
    const complex_t IT_0036 = m_u*conjq(V_u2)*V_us*e_em*IT_0035*U_su_32;
    const complex_t IT_0037 = IT_0034*IT_0036;
    const complex_t IT_0038 = 1.4142135623731*IT_0037;
    const complex_t IT_0039 = m_c*conjq(V_cs)*conjq(V_u2)*e_em*IT_0035*U_su_42;
    const complex_t IT_0040 = IT_0034*IT_0039;
    const complex_t IT_0041 = 1.4142135623731*IT_0040;
    const complex_t IT_0042 = m_t*conjq(V_ts)*conjq(V_u2)*e_em*IT_0035*U_su_52;
    const complex_t IT_0043 = IT_0034*IT_0042;
    const complex_t IT_0044 = 1.4142135623731*IT_0043;
    const complex_t IT_0045 = (complex_t{0, 1})*(IT_0027 + IT_0029 + IT_0031 +
       (-0.5)*IT_0038 + (-0.5)*IT_0041 + (-0.5)*IT_0044);
    const complex_t IT_0046 = cosq(beta);
    const complex_t IT_0047 = cpowq(IT_0046, -1);
    const complex_t IT_0048 = IT_0025*IT_0047;
    const complex_t IT_0049 = m_b*conjq(U_d2)*V_cb*e_em*IT_0035*conjq(U_su_12);
    const complex_t IT_0050 = IT_0048*IT_0049;
    const complex_t IT_0051 = 1.4142135623731*IT_0050;
    const complex_t IT_0052 = m_b*conjq(U_d2)*V_tb*e_em*IT_0035*conjq(U_su_22);
    const complex_t IT_0053 = IT_0048*IT_0052;
    const complex_t IT_0054 = 1.4142135623731*IT_0053;
    const complex_t IT_0055 = cexpq((complex_t{0, -1})*delta_wolf);
    const complex_t IT_0056 = IT_0025*IT_0047*IT_0055;
    const complex_t IT_0057 = m_b*conjq(U_d2)*e_em*IT_0035*conjq(U_su_02)
      *V_ub_mod;
    const complex_t IT_0058 = IT_0056*IT_0057;
    const complex_t IT_0059 = 1.4142135623731*IT_0058;
    const complex_t IT_0060 = (complex_t{0, 1})*(IT_0051 + IT_0054 + IT_0059);
    const complex_t IT_0061 = 0.5*IT_0060;
    const complex_t IT_0062 = IT_0045*IT_0061;
    const complex_t IT_0063 = IT_0023*IT_0062;
    const complex_t IT_0064 = IT_0022*IT_0063;
    const complex_t IT_0065 = powq(m_C1p, 2);
    const complex_t IT_0066 = mty::lt::C0iC(0, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0065, IT_0016, IT_0016, mty::lt::reg_int);
    const complex_t IT_0067 = IT_0012*IT_0066;
    const complex_t IT_0068 = mty::lt::C0iC(3, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0065, IT_0016, IT_0016, mty::lt::reg_int);
    const complex_t IT_0069 = IT_0019*IT_0068;
    const complex_t IT_0070 = IT_0067 + IT_0069;
    const complex_t IT_0071 = 0.101321183642338*m_C1p;
    const complex_t IT_0072 = V_cb*e_em*V_Wp1*conjq(U_su_12);
    const complex_t IT_0073 = IT_0025*IT_0072;
    const complex_t IT_0074 = V_tb*e_em*V_Wp1*conjq(U_su_22);
    const complex_t IT_0075 = IT_0025*IT_0074;
    const complex_t IT_0076 = IT_0025*IT_0055;
    const complex_t IT_0077 = e_em*V_Wp1*conjq(U_su_02)*V_ub_mod;
    const complex_t IT_0078 = IT_0076*IT_0077;
    const complex_t IT_0079 = m_c*V_cb*V_u1*e_em*IT_0035*conjq(U_su_42);
    const complex_t IT_0080 = IT_0034*IT_0079;
    const complex_t IT_0081 = 1.4142135623731*IT_0080;
    const complex_t IT_0082 = m_t*V_tb*V_u1*e_em*IT_0035*conjq(U_su_52);
    const complex_t IT_0083 = IT_0034*IT_0082;
    const complex_t IT_0084 = 1.4142135623731*IT_0083;
    const complex_t IT_0085 = IT_0025*IT_0033*IT_0055;
    const complex_t IT_0086 = m_u*V_u1*e_em*IT_0035*conjq(U_su_32)*V_ub_mod;
    const complex_t IT_0087 = IT_0085*IT_0086;
    const complex_t IT_0088 = 1.4142135623731*IT_0087;
    const complex_t IT_0089 = (complex_t{0, 1})*(IT_0073 + IT_0075 + IT_0078 +
       (-0.5)*IT_0081 + (-0.5)*IT_0084 + (-0.5)*IT_0088);
    const complex_t IT_0090 = m_s*U_d1*conjq(V_cs)*e_em*IT_0035*U_su_12;
    const complex_t IT_0091 = IT_0048*IT_0090;
    const complex_t IT_0092 = 1.4142135623731*IT_0091;
    const complex_t IT_0093 = m_s*U_d1*conjq(V_ts)*e_em*IT_0035*U_su_22;
    const complex_t IT_0094 = IT_0048*IT_0093;
    const complex_t IT_0095 = 1.4142135623731*IT_0094;
    const complex_t IT_0096 = m_s*U_d1*V_us*e_em*IT_0035*U_su_02;
    const complex_t IT_0097 = IT_0048*IT_0096;
    const complex_t IT_0098 = 1.4142135623731*IT_0097;
    const complex_t IT_0099 = (complex_t{0, 1})*(IT_0092 + IT_0095 + IT_0098);
    const complex_t IT_0100 = 0.5*IT_0099;
    const complex_t IT_0101 = IT_0089*IT_0100;
    const complex_t IT_0102 = IT_0071*IT_0101;
    const complex_t IT_0103 = IT_0070*IT_0102;
    const complex_t IT_0104 = V_cb*e_em*V_Wp2*conjq(U_su_10);
    const complex_t IT_0105 = IT_0025*IT_0104;
    const complex_t IT_0106 = V_tb*e_em*V_Wp2*conjq(U_su_20);
    const complex_t IT_0107 = IT_0025*IT_0106;
    const complex_t IT_0108 = e_em*V_Wp2*conjq(U_su_00)*V_ub_mod;
    const complex_t IT_0109 = IT_0076*IT_0108;
    const complex_t IT_0110 = m_c*V_cb*V_u2*e_em*IT_0035*conjq(U_su_40);
    const complex_t IT_0111 = IT_0034*IT_0110;
    const complex_t IT_0112 = 1.4142135623731*IT_0111;
    const complex_t IT_0113 = m_t*V_tb*V_u2*e_em*IT_0035*conjq(U_su_50);
    const complex_t IT_0114 = IT_0034*IT_0113;
    const complex_t IT_0115 = 1.4142135623731*IT_0114;
    const complex_t IT_0116 = m_u*V_u2*e_em*IT_0035*conjq(U_su_30)*V_ub_mod;
    const complex_t IT_0117 = IT_0085*IT_0116;
    const complex_t IT_0118 = 1.4142135623731*IT_0117;
    const complex_t IT_0119 = (complex_t{0, 1})*(IT_0105 + IT_0107 + IT_0109 +
       (-0.5)*IT_0112 + (-0.5)*IT_0115 + (-0.5)*IT_0118);
    const complex_t IT_0120 = m_s*U_d2*conjq(V_cs)*e_em*IT_0035*U_su_10;
    const complex_t IT_0121 = IT_0048*IT_0120;
    const complex_t IT_0122 = 1.4142135623731*IT_0121;
    const complex_t IT_0123 = m_s*U_d2*conjq(V_ts)*e_em*IT_0035*U_su_20;
    const complex_t IT_0124 = IT_0048*IT_0123;
    const complex_t IT_0125 = 1.4142135623731*IT_0124;
    const complex_t IT_0126 = m_s*U_d2*V_us*e_em*IT_0035*U_su_00;
    const complex_t IT_0127 = IT_0048*IT_0126;
    const complex_t IT_0128 = 1.4142135623731*IT_0127;
    const complex_t IT_0129 = (complex_t{0, 1})*(IT_0122 + IT_0125 + IT_0128);
    const complex_t IT_0130 = 0.5*IT_0129;
    const complex_t IT_0131 = IT_0119*IT_0130;
    const complex_t IT_0132 = IT_0023*IT_0131;
    const complex_t IT_0133 = powq(m_su_L, 2);
    const complex_t IT_0134 = mty::lt::C0iC(0, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0015, IT_0133, IT_0133, mty::lt::reg_int);
    const complex_t IT_0135 = IT_0012*IT_0134;
    const complex_t IT_0136 = mty::lt::C0iC(3, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0015, IT_0133, IT_0133, mty::lt::reg_int);
    const complex_t IT_0137 = IT_0019*IT_0136;
    const complex_t IT_0138 = IT_0135 + IT_0137;
    const complex_t IT_0139 = IT_0132*IT_0138;
    const complex_t IT_0140 = powq(m_su_R, 2);
    const complex_t IT_0141 = mty::lt::C0iC(3, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0065, IT_0140, mty::lt::reg_int);
    const complex_t IT_0142 = mty::lt::C0iC(6, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0065, IT_0140, mty::lt::reg_int);
    const complex_t IT_0143 = mty::lt::C0iC(0, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0065, IT_0140, mty::lt::reg_int);
    const complex_t IT_0144 = IT_0141 + IT_0142 + IT_0143;
    const complex_t IT_0145 = IT_0007*IT_0010;
    const complex_t IT_0146 = e_em*U_Wm1*conjq(U_Wm1);
    const complex_t IT_0147 = IT_0145*IT_0146;
    const complex_t IT_0148 = U_d1*conjq(U_d1)*e_em;
    const complex_t IT_0149 = IT_0145*IT_0148;
    const complex_t IT_0150 = (complex_t{0, 1})*(IT_0147 + IT_0149);
    const complex_t IT_0151 = -IT_0150;
    const complex_t IT_0152 = V_cb*e_em*V_Wp1*conjq(U_su_13);
    const complex_t IT_0153 = IT_0025*IT_0152;
    const complex_t IT_0154 = V_tb*e_em*V_Wp1*conjq(U_su_23);
    const complex_t IT_0155 = IT_0025*IT_0154;
    const complex_t IT_0156 = e_em*V_Wp1*conjq(U_su_03)*V_ub_mod;
    const complex_t IT_0157 = IT_0076*IT_0156;
    const complex_t IT_0158 = m_c*V_cb*V_u1*e_em*IT_0035*conjq(U_su_43);
    const complex_t IT_0159 = IT_0034*IT_0158;
    const complex_t IT_0160 = 1.4142135623731*IT_0159;
    const complex_t IT_0161 = m_t*V_tb*V_u1*e_em*IT_0035*conjq(U_su_53);
    const complex_t IT_0162 = IT_0034*IT_0161;
    const complex_t IT_0163 = 1.4142135623731*IT_0162;
    const complex_t IT_0164 = m_u*V_u1*e_em*IT_0035*conjq(U_su_33)*V_ub_mod;
    const complex_t IT_0165 = IT_0085*IT_0164;
    const complex_t IT_0166 = 1.4142135623731*IT_0165;
    const complex_t IT_0167 = (complex_t{0, 1})*(IT_0153 + IT_0155 + IT_0157 +
       (-0.5)*IT_0160 + (-0.5)*IT_0163 + (-0.5)*IT_0166);
    const complex_t IT_0168 = m_s*U_d1*conjq(V_cs)*e_em*IT_0035*U_su_13;
    const complex_t IT_0169 = IT_0048*IT_0168;
    const complex_t IT_0170 = 1.4142135623731*IT_0169;
    const complex_t IT_0171 = m_s*U_d1*conjq(V_ts)*e_em*IT_0035*U_su_23;
    const complex_t IT_0172 = IT_0048*IT_0171;
    const complex_t IT_0173 = 1.4142135623731*IT_0172;
    const complex_t IT_0174 = m_s*U_d1*V_us*e_em*IT_0035*U_su_03;
    const complex_t IT_0175 = IT_0048*IT_0174;
    const complex_t IT_0176 = 1.4142135623731*IT_0175;
    const complex_t IT_0177 = (complex_t{0, 1})*(IT_0170 + IT_0173 + IT_0176);
    const complex_t IT_0178 = 0.5*IT_0177;
    const complex_t IT_0179 = IT_0151*IT_0167*IT_0178;
    const complex_t IT_0180 = IT_0071*IT_0179;
    const complex_t IT_0181 = IT_0144*IT_0180;
    const complex_t IT_0182 = mty::lt::C0iC(3, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0065, IT_0133, mty::lt::reg_int);
    const complex_t IT_0183 = mty::lt::C0iC(6, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0065, IT_0133, mty::lt::reg_int);
    const complex_t IT_0184 = mty::lt::C0iC(0, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0065, IT_0133, mty::lt::reg_int);
    const complex_t IT_0185 = IT_0182 + IT_0183 + IT_0184;
    const complex_t IT_0186 = V_cb*e_em*V_Wp1*conjq(U_su_10);
    const complex_t IT_0187 = IT_0025*IT_0186;
    const complex_t IT_0188 = V_tb*e_em*V_Wp1*conjq(U_su_20);
    const complex_t IT_0189 = IT_0025*IT_0188;
    const complex_t IT_0190 = e_em*V_Wp1*conjq(U_su_00)*V_ub_mod;
    const complex_t IT_0191 = IT_0076*IT_0190;
    const complex_t IT_0192 = m_c*V_cb*V_u1*e_em*IT_0035*conjq(U_su_40);
    const complex_t IT_0193 = IT_0034*IT_0192;
    const complex_t IT_0194 = 1.4142135623731*IT_0193;
    const complex_t IT_0195 = m_t*V_tb*V_u1*e_em*IT_0035*conjq(U_su_50);
    const complex_t IT_0196 = IT_0034*IT_0195;
    const complex_t IT_0197 = 1.4142135623731*IT_0196;
    const complex_t IT_0198 = m_u*V_u1*e_em*IT_0035*conjq(U_su_30)*V_ub_mod;
    const complex_t IT_0199 = IT_0085*IT_0198;
    const complex_t IT_0200 = 1.4142135623731*IT_0199;
    const complex_t IT_0201 = (complex_t{0, 1})*(IT_0187 + IT_0189 + IT_0191 +
       (-0.5)*IT_0194 + (-0.5)*IT_0197 + (-0.5)*IT_0200);
    const complex_t IT_0202 = m_s*U_d1*conjq(V_cs)*e_em*IT_0035*U_su_10;
    const complex_t IT_0203 = IT_0048*IT_0202;
    const complex_t IT_0204 = 1.4142135623731*IT_0203;
    const complex_t IT_0205 = m_s*U_d1*conjq(V_ts)*e_em*IT_0035*U_su_20;
    const complex_t IT_0206 = IT_0048*IT_0205;
    const complex_t IT_0207 = 1.4142135623731*IT_0206;
    const complex_t IT_0208 = m_s*U_d1*V_us*e_em*IT_0035*U_su_00;
    const complex_t IT_0209 = IT_0048*IT_0208;
    const complex_t IT_0210 = 1.4142135623731*IT_0209;
    const complex_t IT_0211 = (complex_t{0, 1})*(IT_0204 + IT_0207 + IT_0210);
    const complex_t IT_0212 = 0.5*IT_0211;
    const complex_t IT_0213 = IT_0151*IT_0201*IT_0212;
    const complex_t IT_0214 = IT_0071*IT_0213;
    const complex_t IT_0215 = IT_0185*IT_0214;
    const complex_t IT_0216 = e_em*V_Wp2*conjq(V_Wp2);
    const complex_t IT_0217 = IT_0145*IT_0216;
    const complex_t IT_0218 = V_u2*conjq(V_u2)*e_em;
    const complex_t IT_0219 = IT_0145*IT_0218;
    const complex_t IT_0220 = (complex_t{0, 1})*(IT_0217 + IT_0219);
    const complex_t IT_0221 = m_b*conjq(U_d2)*V_cb*e_em*IT_0035*conjq(U_su_11);
    const complex_t IT_0222 = IT_0048*IT_0221;
    const complex_t IT_0223 = 1.4142135623731*IT_0222;
    const complex_t IT_0224 = m_b*conjq(U_d2)*V_tb*e_em*IT_0035*conjq(U_su_21);
    const complex_t IT_0225 = IT_0048*IT_0224;
    const complex_t IT_0226 = 1.4142135623731*IT_0225;
    const complex_t IT_0227 = m_b*conjq(U_d2)*e_em*IT_0035*conjq(U_su_01)
      *V_ub_mod;
    const complex_t IT_0228 = IT_0056*IT_0227;
    const complex_t IT_0229 = 1.4142135623731*IT_0228;
    const complex_t IT_0230 = (complex_t{0, 1})*(IT_0223 + IT_0226 + IT_0229);
    const complex_t IT_0231 = 0.5*IT_0230;
    const complex_t IT_0232 = V_us*e_em*conjq(V_Wp2)*U_su_01;
    const complex_t IT_0233 = IT_0025*IT_0232;
    const complex_t IT_0234 = conjq(V_cs)*e_em*conjq(V_Wp2)*U_su_11;
    const complex_t IT_0235 = IT_0025*IT_0234;
    const complex_t IT_0236 = conjq(V_ts)*e_em*conjq(V_Wp2)*U_su_21;
    const complex_t IT_0237 = IT_0025*IT_0236;
    const complex_t IT_0238 = m_u*conjq(V_u2)*V_us*e_em*IT_0035*U_su_31;
    const complex_t IT_0239 = IT_0034*IT_0238;
    const complex_t IT_0240 = 1.4142135623731*IT_0239;
    const complex_t IT_0241 = m_c*conjq(V_cs)*conjq(V_u2)*e_em*IT_0035*U_su_41;
    const complex_t IT_0242 = IT_0034*IT_0241;
    const complex_t IT_0243 = 1.4142135623731*IT_0242;
    const complex_t IT_0244 = m_t*conjq(V_ts)*conjq(V_u2)*e_em*IT_0035*U_su_51;
    const complex_t IT_0245 = IT_0034*IT_0244;
    const complex_t IT_0246 = 1.4142135623731*IT_0245;
    const complex_t IT_0247 = (complex_t{0, 1})*(IT_0233 + IT_0235 + IT_0237 +
       (-0.5)*IT_0240 + (-0.5)*IT_0243 + (-0.5)*IT_0246);
    const complex_t IT_0248 = IT_0220*IT_0231*IT_0247;
    const complex_t IT_0249 = IT_0023*IT_0248;
    const complex_t IT_0250 = powq(m_sc_L, 2);
    const complex_t IT_0251 = mty::lt::C0iC(0, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0015, IT_0250, mty::lt::reg_int);
    const complex_t IT_0252 = mty::lt::C0iC(3, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0015, IT_0250, mty::lt::reg_int);
    const complex_t IT_0253 = mty::lt::C0iC(6, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0015, IT_0250, mty::lt::reg_int);
    const complex_t IT_0254 = IT_0251 + IT_0252 + IT_0253;
    const complex_t IT_0255 = IT_0249*IT_0254;
    const complex_t IT_0256 = V_us*e_em*conjq(V_Wp1)*U_su_02;
    const complex_t IT_0257 = IT_0025*IT_0256;
    const complex_t IT_0258 = conjq(V_cs)*e_em*conjq(V_Wp1)*U_su_12;
    const complex_t IT_0259 = IT_0025*IT_0258;
    const complex_t IT_0260 = conjq(V_ts)*e_em*conjq(V_Wp1)*U_su_22;
    const complex_t IT_0261 = IT_0025*IT_0260;
    const complex_t IT_0262 = m_u*conjq(V_u1)*V_us*e_em*IT_0035*U_su_32;
    const complex_t IT_0263 = IT_0034*IT_0262;
    const complex_t IT_0264 = 1.4142135623731*IT_0263;
    const complex_t IT_0265 = m_c*conjq(V_cs)*conjq(V_u1)*e_em*IT_0035*U_su_42;
    const complex_t IT_0266 = IT_0034*IT_0265;
    const complex_t IT_0267 = 1.4142135623731*IT_0266;
    const complex_t IT_0268 = m_t*conjq(V_ts)*conjq(V_u1)*e_em*IT_0035*U_su_52;
    const complex_t IT_0269 = IT_0034*IT_0268;
    const complex_t IT_0270 = 1.4142135623731*IT_0269;
    const complex_t IT_0271 = (complex_t{0, 1})*(IT_0257 + IT_0259 + IT_0261 +
       (-0.5)*IT_0264 + (-0.5)*IT_0267 + (-0.5)*IT_0270);
    const complex_t IT_0272 = m_b*conjq(U_d1)*V_cb*e_em*IT_0035*conjq(U_su_12);
    const complex_t IT_0273 = IT_0048*IT_0272;
    const complex_t IT_0274 = 1.4142135623731*IT_0273;
    const complex_t IT_0275 = m_b*conjq(U_d1)*V_tb*e_em*IT_0035*conjq(U_su_22);
    const complex_t IT_0276 = IT_0048*IT_0275;
    const complex_t IT_0277 = 1.4142135623731*IT_0276;
    const complex_t IT_0278 = m_b*conjq(U_d1)*e_em*IT_0035*conjq(U_su_02)
      *V_ub_mod;
    const complex_t IT_0279 = IT_0056*IT_0278;
    const complex_t IT_0280 = 1.4142135623731*IT_0279;
    const complex_t IT_0281 = (complex_t{0, 1})*(IT_0274 + IT_0277 + IT_0280);
    const complex_t IT_0282 = 0.5*IT_0281;
    const complex_t IT_0283 = IT_0271*IT_0282;
    const complex_t IT_0284 = IT_0071*IT_0283;
    const complex_t IT_0285 = IT_0070*IT_0284;
    const complex_t IT_0286 = IT_0012*IT_0068;
    const complex_t IT_0287 = m_b*IT_0286;
    const complex_t IT_0288 = mty::lt::C0iC(12, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0065, IT_0016, IT_0016, mty::lt::reg_int);
    const complex_t IT_0289 = IT_0019*IT_0288;
    const complex_t IT_0290 = m_b*IT_0289;
    const complex_t IT_0291 = mty::lt::C0iC(6, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0065, IT_0016, IT_0016, mty::lt::reg_int);
    const complex_t IT_0292 = IT_0012*IT_0291;
    const complex_t IT_0293 = m_s*IT_0292;
    const complex_t IT_0294 = mty::lt::C0iC(15, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0065, IT_0016, IT_0016, mty::lt::reg_int);
    const complex_t IT_0295 = IT_0019*IT_0294;
    const complex_t IT_0296 = m_s*IT_0295;
    const complex_t IT_0297 = IT_0287 + IT_0290 + IT_0293 + IT_0296;
    const complex_t IT_0298 = IT_0100*IT_0282;
    const complex_t IT_0299 = 0.101321183642338*IT_0298;
    const complex_t IT_0300 = IT_0297*IT_0299;
    const complex_t IT_0301 = IT_0012*IT_0020;
    const complex_t IT_0302 = m_b*IT_0301;
    const complex_t IT_0303 = mty::lt::C0iC(12, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0015, IT_0016, IT_0016, mty::lt::reg_int);
    const complex_t IT_0304 = IT_0019*IT_0303;
    const complex_t IT_0305 = m_b*IT_0304;
    const complex_t IT_0306 = mty::lt::C0iC(6, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0015, IT_0016, IT_0016, mty::lt::reg_int);
    const complex_t IT_0307 = IT_0012*IT_0306;
    const complex_t IT_0308 = m_s*IT_0307;
    const complex_t IT_0309 = mty::lt::C0iC(15, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0015, IT_0016, IT_0016, mty::lt::reg_int);
    const complex_t IT_0310 = IT_0019*IT_0309;
    const complex_t IT_0311 = m_s*IT_0310;
    const complex_t IT_0312 = IT_0302 + IT_0305 + IT_0308 + IT_0311;
    const complex_t IT_0313 = m_s*U_d2*conjq(V_cs)*e_em*IT_0035*U_su_12;
    const complex_t IT_0314 = IT_0048*IT_0313;
    const complex_t IT_0315 = 1.4142135623731*IT_0314;
    const complex_t IT_0316 = m_s*U_d2*conjq(V_ts)*e_em*IT_0035*U_su_22;
    const complex_t IT_0317 = IT_0048*IT_0316;
    const complex_t IT_0318 = 1.4142135623731*IT_0317;
    const complex_t IT_0319 = m_s*U_d2*V_us*e_em*IT_0035*U_su_02;
    const complex_t IT_0320 = IT_0048*IT_0319;
    const complex_t IT_0321 = 1.4142135623731*IT_0320;
    const complex_t IT_0322 = (complex_t{0, 1})*(IT_0315 + IT_0318 + IT_0321);
    const complex_t IT_0323 = 0.5*IT_0322;
    const complex_t IT_0324 = IT_0061*IT_0323;
    const complex_t IT_0325 = 0.101321183642338*IT_0324;
    const complex_t IT_0326 = IT_0312*IT_0325;
    const complex_t IT_0327 = IT_0089*IT_0271;
    const complex_t IT_0328 = 0.101321183642338*IT_0327;
    const complex_t IT_0329 = IT_0297*IT_0328;
    const complex_t IT_0330 = V_cb*e_em*V_Wp2*conjq(U_su_12);
    const complex_t IT_0331 = IT_0025*IT_0330;
    const complex_t IT_0332 = V_tb*e_em*V_Wp2*conjq(U_su_22);
    const complex_t IT_0333 = IT_0025*IT_0332;
    const complex_t IT_0334 = e_em*V_Wp2*conjq(U_su_02)*V_ub_mod;
    const complex_t IT_0335 = IT_0076*IT_0334;
    const complex_t IT_0336 = m_c*V_cb*V_u2*e_em*IT_0035*conjq(U_su_42);
    const complex_t IT_0337 = IT_0034*IT_0336;
    const complex_t IT_0338 = 1.4142135623731*IT_0337;
    const complex_t IT_0339 = m_t*V_tb*V_u2*e_em*IT_0035*conjq(U_su_52);
    const complex_t IT_0340 = IT_0034*IT_0339;
    const complex_t IT_0341 = 1.4142135623731*IT_0340;
    const complex_t IT_0342 = m_u*V_u2*e_em*IT_0035*conjq(U_su_32)*V_ub_mod;
    const complex_t IT_0343 = IT_0085*IT_0342;
    const complex_t IT_0344 = 1.4142135623731*IT_0343;
    const complex_t IT_0345 = (complex_t{0, 1})*(IT_0331 + IT_0333 + IT_0335 +
       (-0.5)*IT_0338 + (-0.5)*IT_0341 + (-0.5)*IT_0344);
    const complex_t IT_0346 = IT_0045*IT_0345;
    const complex_t IT_0347 = 0.101321183642338*IT_0346;
    const complex_t IT_0348 = IT_0312*IT_0347;
    const complex_t IT_0349 = IT_0323*IT_0345;
    const complex_t IT_0350 = IT_0023*IT_0349;
    const complex_t IT_0351 = IT_0022*IT_0350;
    const complex_t IT_0352 = mty::lt::C0iC(3, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0065, IT_0133, IT_0133, mty::lt::reg_int);
    const complex_t IT_0353 = IT_0019*IT_0352;
    const complex_t IT_0354 = mty::lt::C0iC(0, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0065, IT_0133, IT_0133, mty::lt::reg_int);
    const complex_t IT_0355 = IT_0012*IT_0354;
    const complex_t IT_0356 = IT_0353 + IT_0355;
    const complex_t IT_0357 = m_b*conjq(U_d1)*V_cb*e_em*IT_0035*conjq(U_su_10);
    const complex_t IT_0358 = IT_0048*IT_0357;
    const complex_t IT_0359 = 1.4142135623731*IT_0358;
    const complex_t IT_0360 = m_b*conjq(U_d1)*V_tb*e_em*IT_0035*conjq(U_su_20);
    const complex_t IT_0361 = IT_0048*IT_0360;
    const complex_t IT_0362 = 1.4142135623731*IT_0361;
    const complex_t IT_0363 = m_b*conjq(U_d1)*e_em*IT_0035*conjq(U_su_00)
      *V_ub_mod;
    const complex_t IT_0364 = IT_0056*IT_0363;
    const complex_t IT_0365 = 1.4142135623731*IT_0364;
    const complex_t IT_0366 = (complex_t{0, 1})*(IT_0359 + IT_0362 + IT_0365);
    const complex_t IT_0367 = 0.5*IT_0366;
    const complex_t IT_0368 = V_us*e_em*conjq(V_Wp1)*U_su_00;
    const complex_t IT_0369 = IT_0025*IT_0368;
    const complex_t IT_0370 = conjq(V_cs)*e_em*conjq(V_Wp1)*U_su_10;
    const complex_t IT_0371 = IT_0025*IT_0370;
    const complex_t IT_0372 = conjq(V_ts)*e_em*conjq(V_Wp1)*U_su_20;
    const complex_t IT_0373 = IT_0025*IT_0372;
    const complex_t IT_0374 = m_u*conjq(V_u1)*V_us*e_em*IT_0035*U_su_30;
    const complex_t IT_0375 = IT_0034*IT_0374;
    const complex_t IT_0376 = 1.4142135623731*IT_0375;
    const complex_t IT_0377 = m_c*conjq(V_cs)*conjq(V_u1)*e_em*IT_0035*U_su_40;
    const complex_t IT_0378 = IT_0034*IT_0377;
    const complex_t IT_0379 = 1.4142135623731*IT_0378;
    const complex_t IT_0380 = m_t*conjq(V_ts)*conjq(V_u1)*e_em*IT_0035*U_su_50;
    const complex_t IT_0381 = IT_0034*IT_0380;
    const complex_t IT_0382 = 1.4142135623731*IT_0381;
    const complex_t IT_0383 = (complex_t{0, 1})*(IT_0369 + IT_0371 + IT_0373 +
       (-0.5)*IT_0376 + (-0.5)*IT_0379 + (-0.5)*IT_0382);
    const complex_t IT_0384 = IT_0367*IT_0383;
    const complex_t IT_0385 = IT_0071*IT_0384;
    const complex_t IT_0386 = IT_0356*IT_0385;
    const complex_t IT_0387 = IT_0212*IT_0367;
    const complex_t IT_0388 = 0.101321183642338*IT_0387;
    const complex_t IT_0389 = IT_0012*IT_0352;
    const complex_t IT_0390 = m_b*IT_0389;
    const complex_t IT_0391 = mty::lt::C0iC(12, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0065, IT_0133, IT_0133, mty::lt::reg_int);
    const complex_t IT_0392 = IT_0019*IT_0391;
    const complex_t IT_0393 = m_b*IT_0392;
    const complex_t IT_0394 = mty::lt::C0iC(6, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0065, IT_0133, IT_0133, mty::lt::reg_int);
    const complex_t IT_0395 = IT_0012*IT_0394;
    const complex_t IT_0396 = m_s*IT_0395;
    const complex_t IT_0397 = mty::lt::C0iC(15, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0065, IT_0133, IT_0133, mty::lt::reg_int);
    const complex_t IT_0398 = IT_0019*IT_0397;
    const complex_t IT_0399 = m_s*IT_0398;
    const complex_t IT_0400 = IT_0390 + IT_0393 + IT_0396 + IT_0399;
    const complex_t IT_0401 = IT_0388*IT_0400;
    const complex_t IT_0402 = m_b*conjq(U_d2)*V_cb*e_em*IT_0035*conjq(U_su_10);
    const complex_t IT_0403 = IT_0048*IT_0402;
    const complex_t IT_0404 = 1.4142135623731*IT_0403;
    const complex_t IT_0405 = m_b*conjq(U_d2)*V_tb*e_em*IT_0035*conjq(U_su_20);
    const complex_t IT_0406 = IT_0048*IT_0405;
    const complex_t IT_0407 = 1.4142135623731*IT_0406;
    const complex_t IT_0408 = m_b*conjq(U_d2)*e_em*IT_0035*conjq(U_su_00)
      *V_ub_mod;
    const complex_t IT_0409 = IT_0056*IT_0408;
    const complex_t IT_0410 = 1.4142135623731*IT_0409;
    const complex_t IT_0411 = (complex_t{0, 1})*(IT_0404 + IT_0407 + IT_0410);
    const complex_t IT_0412 = 0.5*IT_0411;
    const complex_t IT_0413 = V_us*e_em*conjq(V_Wp2)*U_su_00;
    const complex_t IT_0414 = IT_0025*IT_0413;
    const complex_t IT_0415 = conjq(V_cs)*e_em*conjq(V_Wp2)*U_su_10;
    const complex_t IT_0416 = IT_0025*IT_0415;
    const complex_t IT_0417 = conjq(V_ts)*e_em*conjq(V_Wp2)*U_su_20;
    const complex_t IT_0418 = IT_0025*IT_0417;
    const complex_t IT_0419 = m_u*conjq(V_u2)*V_us*e_em*IT_0035*U_su_30;
    const complex_t IT_0420 = IT_0034*IT_0419;
    const complex_t IT_0421 = 1.4142135623731*IT_0420;
    const complex_t IT_0422 = m_c*conjq(V_cs)*conjq(V_u2)*e_em*IT_0035*U_su_40;
    const complex_t IT_0423 = IT_0034*IT_0422;
    const complex_t IT_0424 = 1.4142135623731*IT_0423;
    const complex_t IT_0425 = m_t*conjq(V_ts)*conjq(V_u2)*e_em*IT_0035*U_su_50;
    const complex_t IT_0426 = IT_0034*IT_0425;
    const complex_t IT_0427 = 1.4142135623731*IT_0426;
    const complex_t IT_0428 = (complex_t{0, 1})*(IT_0414 + IT_0416 + IT_0418 +
       (-0.5)*IT_0421 + (-0.5)*IT_0424 + (-0.5)*IT_0427);
    const complex_t IT_0429 = IT_0412*IT_0428;
    const complex_t IT_0430 = IT_0023*IT_0429;
    const complex_t IT_0431 = IT_0138*IT_0430;
    const complex_t IT_0432 = IT_0130*IT_0412;
    const complex_t IT_0433 = 0.101321183642338*IT_0432;
    const complex_t IT_0434 = mty::lt::C0iC(12, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0015, IT_0133, IT_0133, mty::lt::reg_int);
    const complex_t IT_0435 = IT_0019*IT_0434;
    const complex_t IT_0436 = m_b*IT_0435;
    const complex_t IT_0437 = mty::lt::C0iC(15, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0015, IT_0133, IT_0133, mty::lt::reg_int);
    const complex_t IT_0438 = IT_0019*IT_0437;
    const complex_t IT_0439 = m_s*IT_0438;
    const complex_t IT_0440 = IT_0012*IT_0136;
    const complex_t IT_0441 = m_b*IT_0440;
    const complex_t IT_0442 = mty::lt::C0iC(6, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0015, IT_0133, IT_0133, mty::lt::reg_int);
    const complex_t IT_0443 = IT_0012*IT_0442;
    const complex_t IT_0444 = m_s*IT_0443;
    const complex_t IT_0445 = IT_0436 + IT_0439 + IT_0441 + IT_0444;
    const complex_t IT_0446 = IT_0433*IT_0445;
    const complex_t IT_0447 = IT_0201*IT_0383;
    const complex_t IT_0448 = 0.101321183642338*IT_0447;
    const complex_t IT_0449 = IT_0400*IT_0448;
    const complex_t IT_0450 = IT_0119*IT_0428;
    const complex_t IT_0451 = 0.101321183642338*IT_0450;
    const complex_t IT_0452 = IT_0445*IT_0451;
    const complex_t IT_0453 = IT_0201*IT_0212;
    const complex_t IT_0454 = IT_0071*IT_0453;
    const complex_t IT_0455 = IT_0356*IT_0454;
    const complex_t IT_0456 = e_em*U_Wm1*conjq(U_Wm2);
    const complex_t IT_0457 = IT_0145*IT_0456;
    const complex_t IT_0458 = U_d1*conjq(U_d2)*e_em;
    const complex_t IT_0459 = IT_0145*IT_0458;
    const complex_t IT_0460 = (complex_t{0, 1})*(IT_0457 + IT_0459);
    const complex_t IT_0461 = -IT_0460;
    const complex_t IT_0462 = m_b*conjq(U_d1)*V_cb*e_em*IT_0035*conjq(U_su_11);
    const complex_t IT_0463 = IT_0048*IT_0462;
    const complex_t IT_0464 = 1.4142135623731*IT_0463;
    const complex_t IT_0465 = m_b*conjq(U_d1)*V_tb*e_em*IT_0035*conjq(U_su_21);
    const complex_t IT_0466 = IT_0048*IT_0465;
    const complex_t IT_0467 = 1.4142135623731*IT_0466;
    const complex_t IT_0468 = m_b*conjq(U_d1)*e_em*IT_0035*conjq(U_su_01)
      *V_ub_mod;
    const complex_t IT_0469 = IT_0056*IT_0468;
    const complex_t IT_0470 = 1.4142135623731*IT_0469;
    const complex_t IT_0471 = (complex_t{0, 1})*(IT_0464 + IT_0467 + IT_0470);
    const complex_t IT_0472 = 0.5*IT_0471;
    const complex_t IT_0473 = m_s*U_d2*conjq(V_cs)*e_em*IT_0035*U_su_11;
    const complex_t IT_0474 = IT_0048*IT_0473;
    const complex_t IT_0475 = 1.4142135623731*IT_0474;
    const complex_t IT_0476 = m_s*U_d2*conjq(V_ts)*e_em*IT_0035*U_su_21;
    const complex_t IT_0477 = IT_0048*IT_0476;
    const complex_t IT_0478 = 1.4142135623731*IT_0477;
    const complex_t IT_0479 = m_s*U_d2*V_us*e_em*IT_0035*U_su_01;
    const complex_t IT_0480 = IT_0048*IT_0479;
    const complex_t IT_0481 = 1.4142135623731*IT_0480;
    const complex_t IT_0482 = (complex_t{0, 1})*(IT_0475 + IT_0478 + IT_0481);
    const complex_t IT_0483 = 0.5*IT_0482;
    const complex_t IT_0484 = IT_0461*IT_0472*IT_0483;
    const complex_t IT_0485 = 0.101321183642338*IT_0484;
    const complex_t IT_0486 = mty::lt::C0iC(6, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0015, IT_0250, mty::lt::reg_int);
    const complex_t IT_0487 = m_b*IT_0486;
    const complex_t IT_0488 = mty::lt::C0iC(3, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0015, IT_0250, mty::lt::reg_int);
    const complex_t IT_0489 = m_b*IT_0488;
    const complex_t IT_0490 = mty::lt::C0iC(12, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0015, IT_0250, mty::lt::reg_int);
    const complex_t IT_0491 = m_b*IT_0490;
    const complex_t IT_0492 = mty::lt::C0iC(18, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0015, IT_0250, mty::lt::reg_int);
    const complex_t IT_0493 = m_b*IT_0492;
    const complex_t IT_0494 = IT_0487 + IT_0489 + IT_0491 + IT_0493;
    const complex_t IT_0495 = m_s*IT_0488;
    const complex_t IT_0496 = mty::lt::C0iC(15, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0015, IT_0250, mty::lt::reg_int);
    const complex_t IT_0497 = m_b*IT_0496;
    const complex_t IT_0498 = m_s*IT_0490;
    const complex_t IT_0499 = m_s*IT_0496;
    const complex_t IT_0500 = -IT_0495 + 2*IT_0497 + -IT_0498 + -IT_0499;
    const complex_t IT_0501 = IT_0494 + IT_0500;
    const complex_t IT_0502 = IT_0485*IT_0501;
    const complex_t IT_0503 = IT_0282*IT_0323*IT_0461;
    const complex_t IT_0504 = 0.101321183642338*IT_0503;
    const complex_t IT_0505 = mty::lt::C0iC(6, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0015, IT_0016, mty::lt::reg_int);
    const complex_t IT_0506 = m_b*IT_0505;
    const complex_t IT_0507 = mty::lt::C0iC(3, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0015, IT_0016, mty::lt::reg_int);
    const complex_t IT_0508 = m_b*IT_0507;
    const complex_t IT_0509 = mty::lt::C0iC(12, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0015, IT_0016, mty::lt::reg_int);
    const complex_t IT_0510 = m_b*IT_0509;
    const complex_t IT_0511 = mty::lt::C0iC(18, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0015, IT_0016, mty::lt::reg_int);
    const complex_t IT_0512 = m_b*IT_0511;
    const complex_t IT_0513 = IT_0506 + IT_0508 + IT_0510 + IT_0512;
    const complex_t IT_0514 = m_s*IT_0507;
    const complex_t IT_0515 = mty::lt::C0iC(15, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0015, IT_0016, mty::lt::reg_int);
    const complex_t IT_0516 = m_b*IT_0515;
    const complex_t IT_0517 = m_s*IT_0509;
    const complex_t IT_0518 = m_s*IT_0515;
    const complex_t IT_0519 = -IT_0514 + 2*IT_0516 + -IT_0517 + -IT_0518;
    const complex_t IT_0520 = IT_0513 + IT_0519;
    const complex_t IT_0521 = IT_0504*IT_0520;
    const complex_t IT_0522 = m_b*conjq(U_d1)*V_cb*e_em*IT_0035*conjq(U_su_13);
    const complex_t IT_0523 = IT_0048*IT_0522;
    const complex_t IT_0524 = 1.4142135623731*IT_0523;
    const complex_t IT_0525 = m_b*conjq(U_d1)*V_tb*e_em*IT_0035*conjq(U_su_23);
    const complex_t IT_0526 = IT_0048*IT_0525;
    const complex_t IT_0527 = 1.4142135623731*IT_0526;
    const complex_t IT_0528 = m_b*conjq(U_d1)*e_em*IT_0035*conjq(U_su_03)
      *V_ub_mod;
    const complex_t IT_0529 = IT_0056*IT_0528;
    const complex_t IT_0530 = 1.4142135623731*IT_0529;
    const complex_t IT_0531 = (complex_t{0, 1})*(IT_0524 + IT_0527 + IT_0530);
    const complex_t IT_0532 = 0.5*IT_0531;
    const complex_t IT_0533 = m_s*U_d2*conjq(V_cs)*e_em*IT_0035*U_su_13;
    const complex_t IT_0534 = IT_0048*IT_0533;
    const complex_t IT_0535 = 1.4142135623731*IT_0534;
    const complex_t IT_0536 = m_s*U_d2*conjq(V_ts)*e_em*IT_0035*U_su_23;
    const complex_t IT_0537 = IT_0048*IT_0536;
    const complex_t IT_0538 = 1.4142135623731*IT_0537;
    const complex_t IT_0539 = m_s*U_d2*V_us*e_em*IT_0035*U_su_03;
    const complex_t IT_0540 = IT_0048*IT_0539;
    const complex_t IT_0541 = 1.4142135623731*IT_0540;
    const complex_t IT_0542 = (complex_t{0, 1})*(IT_0535 + IT_0538 + IT_0541);
    const complex_t IT_0543 = 0.5*IT_0542;
    const complex_t IT_0544 = IT_0461*IT_0532*IT_0543;
    const complex_t IT_0545 = 0.101321183642338*IT_0544;
    const complex_t IT_0546 = mty::lt::C0iC(3, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0015, IT_0140, mty::lt::reg_int);
    const complex_t IT_0547 = m_b*IT_0546;
    const complex_t IT_0548 = mty::lt::C0iC(12, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0015, IT_0140, mty::lt::reg_int);
    const complex_t IT_0549 = m_b*IT_0548;
    const complex_t IT_0550 = mty::lt::C0iC(6, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0015, IT_0140, mty::lt::reg_int);
    const complex_t IT_0551 = m_b*IT_0550;
    const complex_t IT_0552 = mty::lt::C0iC(18, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0015, IT_0140, mty::lt::reg_int);
    const complex_t IT_0553 = m_b*IT_0552;
    const complex_t IT_0554 = IT_0547 + IT_0549 + IT_0551 + IT_0553;
    const complex_t IT_0555 = m_s*IT_0546;
    const complex_t IT_0556 = mty::lt::C0iC(15, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0015, IT_0140, mty::lt::reg_int);
    const complex_t IT_0557 = m_b*IT_0556;
    const complex_t IT_0558 = m_s*IT_0548;
    const complex_t IT_0559 = m_s*IT_0556;
    const complex_t IT_0560 = -IT_0555 + 2*IT_0557 + -IT_0558 + -IT_0559;
    const complex_t IT_0561 = IT_0554 + IT_0560;
    const complex_t IT_0562 = IT_0545*IT_0561;
    const complex_t IT_0563 = IT_0130*IT_0367*IT_0461;
    const complex_t IT_0564 = 0.101321183642338*IT_0563;
    const complex_t IT_0565 = mty::lt::C0iC(3, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0015, IT_0133, mty::lt::reg_int);
    const complex_t IT_0566 = m_b*IT_0565;
    const complex_t IT_0567 = mty::lt::C0iC(12, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0015, IT_0133, mty::lt::reg_int);
    const complex_t IT_0568 = m_b*IT_0567;
    const complex_t IT_0569 = mty::lt::C0iC(6, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0015, IT_0133, mty::lt::reg_int);
    const complex_t IT_0570 = m_b*IT_0569;
    const complex_t IT_0571 = mty::lt::C0iC(18, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0015, IT_0133, mty::lt::reg_int);
    const complex_t IT_0572 = m_b*IT_0571;
    const complex_t IT_0573 = IT_0566 + IT_0568 + IT_0570 + IT_0572;
    const complex_t IT_0574 = m_s*IT_0565;
    const complex_t IT_0575 = mty::lt::C0iC(15, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0015, IT_0133, mty::lt::reg_int);
    const complex_t IT_0576 = m_b*IT_0575;
    const complex_t IT_0577 = m_s*IT_0567;
    const complex_t IT_0578 = m_s*IT_0575;
    const complex_t IT_0579 = -IT_0574 + 2*IT_0576 + -IT_0577 + -IT_0578;
    const complex_t IT_0580 = IT_0573 + IT_0579;
    const complex_t IT_0581 = IT_0564*IT_0580;
    const complex_t IT_0582 = m_s*U_d2*conjq(V_cs)*e_em*IT_0035*U_su_14;
    const complex_t IT_0583 = IT_0048*IT_0582;
    const complex_t IT_0584 = 1.4142135623731*IT_0583;
    const complex_t IT_0585 = m_s*U_d2*conjq(V_ts)*e_em*IT_0035*U_su_24;
    const complex_t IT_0586 = IT_0048*IT_0585;
    const complex_t IT_0587 = 1.4142135623731*IT_0586;
    const complex_t IT_0588 = m_s*U_d2*V_us*e_em*IT_0035*U_su_04;
    const complex_t IT_0589 = IT_0048*IT_0588;
    const complex_t IT_0590 = 1.4142135623731*IT_0589;
    const complex_t IT_0591 = (complex_t{0, 1})*(IT_0584 + IT_0587 + IT_0590);
    const complex_t IT_0592 = 0.5*IT_0591;
    const complex_t IT_0593 = m_b*conjq(U_d1)*V_cb*e_em*IT_0035*conjq(U_su_14);
    const complex_t IT_0594 = IT_0048*IT_0593;
    const complex_t IT_0595 = 1.4142135623731*IT_0594;
    const complex_t IT_0596 = m_b*conjq(U_d1)*V_tb*e_em*IT_0035*conjq(U_su_24);
    const complex_t IT_0597 = IT_0048*IT_0596;
    const complex_t IT_0598 = 1.4142135623731*IT_0597;
    const complex_t IT_0599 = m_b*conjq(U_d1)*e_em*IT_0035*conjq(U_su_04)
      *V_ub_mod;
    const complex_t IT_0600 = IT_0056*IT_0599;
    const complex_t IT_0601 = 1.4142135623731*IT_0600;
    const complex_t IT_0602 = (complex_t{0, 1})*(IT_0595 + IT_0598 + IT_0601);
    const complex_t IT_0603 = 0.5*IT_0602;
    const complex_t IT_0604 = IT_0461*IT_0592*IT_0603;
    const complex_t IT_0605 = 0.101321183642338*IT_0604;
    const complex_t IT_0606 = powq(m_sc_R, 2);
    const complex_t IT_0607 = mty::lt::C0iC(3, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0015, IT_0606, mty::lt::reg_int);
    const complex_t IT_0608 = m_b*IT_0607;
    const complex_t IT_0609 = mty::lt::C0iC(12, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0015, IT_0606, mty::lt::reg_int);
    const complex_t IT_0610 = m_b*IT_0609;
    const complex_t IT_0611 = mty::lt::C0iC(6, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0015, IT_0606, mty::lt::reg_int);
    const complex_t IT_0612 = m_b*IT_0611;
    const complex_t IT_0613 = mty::lt::C0iC(18, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0015, IT_0606, mty::lt::reg_int);
    const complex_t IT_0614 = m_b*IT_0613;
    const complex_t IT_0615 = IT_0608 + IT_0610 + IT_0612 + IT_0614;
    const complex_t IT_0616 = m_s*IT_0607;
    const complex_t IT_0617 = mty::lt::C0iC(15, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0015, IT_0606, mty::lt::reg_int);
    const complex_t IT_0618 = m_b*IT_0617;
    const complex_t IT_0619 = m_s*IT_0609;
    const complex_t IT_0620 = m_s*IT_0617;
    const complex_t IT_0621 = -IT_0616 + 2*IT_0618 + -IT_0619 + -IT_0620;
    const complex_t IT_0622 = IT_0615 + IT_0621;
    const complex_t IT_0623 = IT_0605*IT_0622;
    const complex_t IT_0624 = m_b*conjq(U_d1)*V_cb*e_em*IT_0035*conjq(U_su_15);
    const complex_t IT_0625 = IT_0048*IT_0624;
    const complex_t IT_0626 = 1.4142135623731*IT_0625;
    const complex_t IT_0627 = m_b*conjq(U_d1)*V_tb*e_em*IT_0035*conjq(U_su_25);
    const complex_t IT_0628 = IT_0048*IT_0627;
    const complex_t IT_0629 = 1.4142135623731*IT_0628;
    const complex_t IT_0630 = m_b*conjq(U_d1)*e_em*IT_0035*conjq(U_su_05)
      *V_ub_mod;
    const complex_t IT_0631 = IT_0056*IT_0630;
    const complex_t IT_0632 = 1.4142135623731*IT_0631;
    const complex_t IT_0633 = (complex_t{0, 1})*(IT_0626 + IT_0629 + IT_0632);
    const complex_t IT_0634 = 0.5*IT_0633;
    const complex_t IT_0635 = m_s*U_d2*conjq(V_cs)*e_em*IT_0035*U_su_15;
    const complex_t IT_0636 = IT_0048*IT_0635;
    const complex_t IT_0637 = 1.4142135623731*IT_0636;
    const complex_t IT_0638 = m_s*U_d2*conjq(V_ts)*e_em*IT_0035*U_su_25;
    const complex_t IT_0639 = IT_0048*IT_0638;
    const complex_t IT_0640 = 1.4142135623731*IT_0639;
    const complex_t IT_0641 = m_s*U_d2*V_us*e_em*IT_0035*U_su_05;
    const complex_t IT_0642 = IT_0048*IT_0641;
    const complex_t IT_0643 = 1.4142135623731*IT_0642;
    const complex_t IT_0644 = (complex_t{0, 1})*(IT_0637 + IT_0640 + IT_0643);
    const complex_t IT_0645 = 0.5*IT_0644;
    const complex_t IT_0646 = IT_0461*IT_0634*IT_0645;
    const complex_t IT_0647 = 0.101321183642338*IT_0646;
    const complex_t IT_0648 = powq(m_st_R, 2);
    const complex_t IT_0649 = mty::lt::C0iC(3, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0015, IT_0648, mty::lt::reg_int);
    const complex_t IT_0650 = m_b*IT_0649;
    const complex_t IT_0651 = mty::lt::C0iC(12, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0015, IT_0648, mty::lt::reg_int);
    const complex_t IT_0652 = m_b*IT_0651;
    const complex_t IT_0653 = mty::lt::C0iC(6, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0015, IT_0648, mty::lt::reg_int);
    const complex_t IT_0654 = m_b*IT_0653;
    const complex_t IT_0655 = mty::lt::C0iC(18, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0015, IT_0648, mty::lt::reg_int);
    const complex_t IT_0656 = m_b*IT_0655;
    const complex_t IT_0657 = IT_0650 + IT_0652 + IT_0654 + IT_0656;
    const complex_t IT_0658 = m_s*IT_0649;
    const complex_t IT_0659 = mty::lt::C0iC(15, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0015, IT_0648, mty::lt::reg_int);
    const complex_t IT_0660 = m_b*IT_0659;
    const complex_t IT_0661 = m_s*IT_0651;
    const complex_t IT_0662 = m_s*IT_0659;
    const complex_t IT_0663 = -IT_0658 + 2*IT_0660 + -IT_0661 + -IT_0662;
    const complex_t IT_0664 = IT_0657 + IT_0663;
    const complex_t IT_0665 = IT_0647*IT_0664;
    const complex_t IT_0666 = IT_0089*IT_0323*IT_0461;
    const complex_t IT_0667 = IT_0071*IT_0666;
    const complex_t IT_0668 = mty::lt::C0iC(0, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0015, IT_0016, mty::lt::reg_int);
    const complex_t IT_0669 = IT_0505 + IT_0507 + IT_0668;
    const complex_t IT_0670 = IT_0667*IT_0669;
    const complex_t IT_0671 = mty::lt::C0iC(0, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0015, IT_0250, mty::lt::reg_int);
    const complex_t IT_0672 = IT_0486 + IT_0488 + IT_0671;
    const complex_t IT_0673 = V_cb*e_em*V_Wp1*conjq(U_su_11);
    const complex_t IT_0674 = IT_0025*IT_0673;
    const complex_t IT_0675 = V_tb*e_em*V_Wp1*conjq(U_su_21);
    const complex_t IT_0676 = IT_0025*IT_0675;
    const complex_t IT_0677 = e_em*V_Wp1*conjq(U_su_01)*V_ub_mod;
    const complex_t IT_0678 = IT_0076*IT_0677;
    const complex_t IT_0679 = m_c*V_cb*V_u1*e_em*IT_0035*conjq(U_su_41);
    const complex_t IT_0680 = IT_0034*IT_0679;
    const complex_t IT_0681 = 1.4142135623731*IT_0680;
    const complex_t IT_0682 = m_t*V_tb*V_u1*e_em*IT_0035*conjq(U_su_51);
    const complex_t IT_0683 = IT_0034*IT_0682;
    const complex_t IT_0684 = 1.4142135623731*IT_0683;
    const complex_t IT_0685 = m_u*V_u1*e_em*IT_0035*conjq(U_su_31)*V_ub_mod;
    const complex_t IT_0686 = IT_0085*IT_0685;
    const complex_t IT_0687 = 1.4142135623731*IT_0686;
    const complex_t IT_0688 = (complex_t{0, 1})*(IT_0674 + IT_0676 + IT_0678 +
       (-0.5)*IT_0681 + (-0.5)*IT_0684 + (-0.5)*IT_0687);
    const complex_t IT_0689 = IT_0461*IT_0483*IT_0688;
    const complex_t IT_0690 = IT_0071*IT_0689;
    const complex_t IT_0691 = IT_0672*IT_0690;
    const complex_t IT_0692 = IT_0167*IT_0461*IT_0543;
    const complex_t IT_0693 = IT_0071*IT_0692;
    const complex_t IT_0694 = mty::lt::C0iC(0, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0015, IT_0140, mty::lt::reg_int);
    const complex_t IT_0695 = IT_0546 + IT_0550 + IT_0694;
    const complex_t IT_0696 = IT_0693*IT_0695;
    const complex_t IT_0697 = IT_0130*IT_0201*IT_0461;
    const complex_t IT_0698 = IT_0071*IT_0697;
    const complex_t IT_0699 = mty::lt::C0iC(0, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0015, IT_0133, mty::lt::reg_int);
    const complex_t IT_0700 = IT_0565 + IT_0569 + IT_0699;
    const complex_t IT_0701 = IT_0698*IT_0700;
    const complex_t IT_0702 = V_cb*e_em*V_Wp1*conjq(U_su_14);
    const complex_t IT_0703 = IT_0025*IT_0702;
    const complex_t IT_0704 = V_tb*e_em*V_Wp1*conjq(U_su_24);
    const complex_t IT_0705 = IT_0025*IT_0704;
    const complex_t IT_0706 = e_em*V_Wp1*conjq(U_su_04)*V_ub_mod;
    const complex_t IT_0707 = IT_0076*IT_0706;
    const complex_t IT_0708 = m_c*V_cb*V_u1*e_em*IT_0035*conjq(U_su_44);
    const complex_t IT_0709 = IT_0034*IT_0708;
    const complex_t IT_0710 = 1.4142135623731*IT_0709;
    const complex_t IT_0711 = m_t*V_tb*V_u1*e_em*IT_0035*conjq(U_su_54);
    const complex_t IT_0712 = IT_0034*IT_0711;
    const complex_t IT_0713 = 1.4142135623731*IT_0712;
    const complex_t IT_0714 = m_u*V_u1*e_em*IT_0035*conjq(U_su_34)*V_ub_mod;
    const complex_t IT_0715 = IT_0085*IT_0714;
    const complex_t IT_0716 = 1.4142135623731*IT_0715;
    const complex_t IT_0717 = (complex_t{0, 1})*(IT_0703 + IT_0705 + IT_0707 +
       (-0.5)*IT_0710 + (-0.5)*IT_0713 + (-0.5)*IT_0716);
    const complex_t IT_0718 = IT_0461*IT_0592*IT_0717;
    const complex_t IT_0719 = IT_0071*IT_0718;
    const complex_t IT_0720 = mty::lt::C0iC(0, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0015, IT_0606, mty::lt::reg_int);
    const complex_t IT_0721 = IT_0607 + IT_0611 + IT_0720;
    const complex_t IT_0722 = IT_0719*IT_0721;
    const complex_t IT_0723 = V_cb*e_em*V_Wp1*conjq(U_su_15);
    const complex_t IT_0724 = IT_0025*IT_0723;
    const complex_t IT_0725 = V_tb*e_em*V_Wp1*conjq(U_su_25);
    const complex_t IT_0726 = IT_0025*IT_0725;
    const complex_t IT_0727 = e_em*V_Wp1*conjq(U_su_05)*V_ub_mod;
    const complex_t IT_0728 = IT_0076*IT_0727;
    const complex_t IT_0729 = m_c*V_cb*V_u1*e_em*IT_0035*conjq(U_su_45);
    const complex_t IT_0730 = IT_0034*IT_0729;
    const complex_t IT_0731 = 1.4142135623731*IT_0730;
    const complex_t IT_0732 = m_t*V_tb*V_u1*e_em*IT_0035*conjq(U_su_55);
    const complex_t IT_0733 = IT_0034*IT_0732;
    const complex_t IT_0734 = 1.4142135623731*IT_0733;
    const complex_t IT_0735 = m_u*V_u1*e_em*IT_0035*conjq(U_su_35)*V_ub_mod;
    const complex_t IT_0736 = IT_0085*IT_0735;
    const complex_t IT_0737 = 1.4142135623731*IT_0736;
    const complex_t IT_0738 = (complex_t{0, 1})*(IT_0724 + IT_0726 + IT_0728 +
       (-0.5)*IT_0731 + (-0.5)*IT_0734 + (-0.5)*IT_0737);
    const complex_t IT_0739 = IT_0461*IT_0645*IT_0738;
    const complex_t IT_0740 = IT_0071*IT_0739;
    const complex_t IT_0741 = mty::lt::C0iC(0, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0015, IT_0648, mty::lt::reg_int);
    const complex_t IT_0742 = IT_0649 + IT_0653 + IT_0741;
    const complex_t IT_0743 = IT_0740*IT_0742;
    const complex_t IT_0744 = m_s*U_d1*conjq(V_cs)*e_em*IT_0035*U_su_11;
    const complex_t IT_0745 = IT_0048*IT_0744;
    const complex_t IT_0746 = 1.4142135623731*IT_0745;
    const complex_t IT_0747 = m_s*U_d1*conjq(V_ts)*e_em*IT_0035*U_su_21;
    const complex_t IT_0748 = IT_0048*IT_0747;
    const complex_t IT_0749 = 1.4142135623731*IT_0748;
    const complex_t IT_0750 = m_s*U_d1*V_us*e_em*IT_0035*U_su_01;
    const complex_t IT_0751 = IT_0048*IT_0750;
    const complex_t IT_0752 = 1.4142135623731*IT_0751;
    const complex_t IT_0753 = (complex_t{0, 1})*(IT_0746 + IT_0749 + IT_0752);
    const complex_t IT_0754 = 0.5*IT_0753;
    const complex_t IT_0755 = IT_0151*IT_0472*IT_0754;
    const complex_t IT_0756 = 0.101321183642338*IT_0755;
    const complex_t IT_0757 = mty::lt::C0iC(6, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0065, IT_0250, mty::lt::reg_int);
    const complex_t IT_0758 = m_b*IT_0757;
    const complex_t IT_0759 = mty::lt::C0iC(3, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0065, IT_0250, mty::lt::reg_int);
    const complex_t IT_0760 = m_b*IT_0759;
    const complex_t IT_0761 = mty::lt::C0iC(12, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0065, IT_0250, mty::lt::reg_int);
    const complex_t IT_0762 = m_b*IT_0761;
    const complex_t IT_0763 = mty::lt::C0iC(18, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0065, IT_0250, mty::lt::reg_int);
    const complex_t IT_0764 = m_b*IT_0763;
    const complex_t IT_0765 = IT_0758 + IT_0760 + IT_0762 + IT_0764;
    const complex_t IT_0766 = m_s*IT_0759;
    const complex_t IT_0767 = mty::lt::C0iC(15, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0065, IT_0250, mty::lt::reg_int);
    const complex_t IT_0768 = m_b*IT_0767;
    const complex_t IT_0769 = m_s*IT_0767;
    const complex_t IT_0770 = m_s*IT_0761;
    const complex_t IT_0771 = -IT_0766 + 2*IT_0768 + -IT_0769 + -IT_0770;
    const complex_t IT_0772 = IT_0765 + IT_0771;
    const complex_t IT_0773 = IT_0756*IT_0772;
    const complex_t IT_0774 = IT_0100*IT_0151*IT_0282;
    const complex_t IT_0775 = 0.101321183642338*IT_0774;
    const complex_t IT_0776 = mty::lt::C0iC(6, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0065, IT_0016, mty::lt::reg_int);
    const complex_t IT_0777 = m_b*IT_0776;
    const complex_t IT_0778 = mty::lt::C0iC(3, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0065, IT_0016, mty::lt::reg_int);
    const complex_t IT_0779 = m_b*IT_0778;
    const complex_t IT_0780 = mty::lt::C0iC(12, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0065, IT_0016, mty::lt::reg_int);
    const complex_t IT_0781 = m_b*IT_0780;
    const complex_t IT_0782 = mty::lt::C0iC(18, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0065, IT_0016, mty::lt::reg_int);
    const complex_t IT_0783 = m_b*IT_0782;
    const complex_t IT_0784 = IT_0777 + IT_0779 + IT_0781 + IT_0783;
    const complex_t IT_0785 = m_s*IT_0778;
    const complex_t IT_0786 = mty::lt::C0iC(15, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0065, IT_0016, mty::lt::reg_int);
    const complex_t IT_0787 = m_b*IT_0786;
    const complex_t IT_0788 = m_s*IT_0780;
    const complex_t IT_0789 = m_s*IT_0786;
    const complex_t IT_0790 = -IT_0785 + 2*IT_0787 + -IT_0788 + -IT_0789;
    const complex_t IT_0791 = IT_0784 + IT_0790;
    const complex_t IT_0792 = IT_0775*IT_0791;
    const complex_t IT_0793 = IT_0151*IT_0178*IT_0532;
    const complex_t IT_0794 = 0.101321183642338*IT_0793;
    const complex_t IT_0795 = m_b*IT_0142;
    const complex_t IT_0796 = mty::lt::C0iC(12, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0065, IT_0140, mty::lt::reg_int);
    const complex_t IT_0797 = m_b*IT_0796;
    const complex_t IT_0798 = mty::lt::C0iC(18, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0065, IT_0140, mty::lt::reg_int);
    const complex_t IT_0799 = m_b*IT_0798;
    const complex_t IT_0800 = m_b*IT_0141;
    const complex_t IT_0801 = IT_0795 + IT_0797 + IT_0799 + IT_0800;
    const complex_t IT_0802 = m_s*IT_0141;
    const complex_t IT_0803 = mty::lt::C0iC(15, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0065, IT_0140, mty::lt::reg_int);
    const complex_t IT_0804 = m_b*IT_0803;
    const complex_t IT_0805 = m_s*IT_0796;
    const complex_t IT_0806 = m_s*IT_0803;
    const complex_t IT_0807 = -IT_0802 + 2*IT_0804 + -IT_0805 + -IT_0806;
    const complex_t IT_0808 = IT_0801 + IT_0807;
    const complex_t IT_0809 = IT_0794*IT_0808;
    const complex_t IT_0810 = IT_0151*IT_0212*IT_0367;
    const complex_t IT_0811 = 0.101321183642338*IT_0810;
    const complex_t IT_0812 = m_b*IT_0183;
    const complex_t IT_0813 = m_b*IT_0182;
    const complex_t IT_0814 = mty::lt::C0iC(12, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0065, IT_0133, mty::lt::reg_int);
    const complex_t IT_0815 = m_b*IT_0814;
    const complex_t IT_0816 = mty::lt::C0iC(18, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0065, IT_0133, mty::lt::reg_int);
    const complex_t IT_0817 = m_b*IT_0816;
    const complex_t IT_0818 = IT_0812 + IT_0813 + IT_0815 + IT_0817;
    const complex_t IT_0819 = m_s*IT_0182;
    const complex_t IT_0820 = mty::lt::C0iC(15, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0065, IT_0133, mty::lt::reg_int);
    const complex_t IT_0821 = m_b*IT_0820;
    const complex_t IT_0822 = m_s*IT_0814;
    const complex_t IT_0823 = m_s*IT_0820;
    const complex_t IT_0824 = -IT_0819 + 2*IT_0821 + -IT_0822 + -IT_0823;
    const complex_t IT_0825 = IT_0818 + IT_0824;
    const complex_t IT_0826 = IT_0811*IT_0825;
    const complex_t IT_0827 = m_s*U_d1*conjq(V_cs)*e_em*IT_0035*U_su_14;
    const complex_t IT_0828 = IT_0048*IT_0827;
    const complex_t IT_0829 = 1.4142135623731*IT_0828;
    const complex_t IT_0830 = m_s*U_d1*conjq(V_ts)*e_em*IT_0035*U_su_24;
    const complex_t IT_0831 = IT_0048*IT_0830;
    const complex_t IT_0832 = 1.4142135623731*IT_0831;
    const complex_t IT_0833 = m_s*U_d1*V_us*e_em*IT_0035*U_su_04;
    const complex_t IT_0834 = IT_0048*IT_0833;
    const complex_t IT_0835 = 1.4142135623731*IT_0834;
    const complex_t IT_0836 = (complex_t{0, 1})*(IT_0829 + IT_0832 + IT_0835);
    const complex_t IT_0837 = 0.5*IT_0836;
    const complex_t IT_0838 = IT_0151*IT_0603*IT_0837;
    const complex_t IT_0839 = 0.101321183642338*IT_0838;
    const complex_t IT_0840 = mty::lt::C0iC(6, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0065, IT_0606, mty::lt::reg_int);
    const complex_t IT_0841 = m_b*IT_0840;
    const complex_t IT_0842 = mty::lt::C0iC(3, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0065, IT_0606, mty::lt::reg_int);
    const complex_t IT_0843 = m_b*IT_0842;
    const complex_t IT_0844 = mty::lt::C0iC(12, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0065, IT_0606, mty::lt::reg_int);
    const complex_t IT_0845 = m_b*IT_0844;
    const complex_t IT_0846 = mty::lt::C0iC(18, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0065, IT_0606, mty::lt::reg_int);
    const complex_t IT_0847 = m_b*IT_0846;
    const complex_t IT_0848 = IT_0841 + IT_0843 + IT_0845 + IT_0847;
    const complex_t IT_0849 = m_s*IT_0842;
    const complex_t IT_0850 = mty::lt::C0iC(15, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0065, IT_0606, mty::lt::reg_int);
    const complex_t IT_0851 = m_b*IT_0850;
    const complex_t IT_0852 = m_s*IT_0844;
    const complex_t IT_0853 = m_s*IT_0850;
    const complex_t IT_0854 = -IT_0849 + 2*IT_0851 + -IT_0852 + -IT_0853;
    const complex_t IT_0855 = IT_0848 + IT_0854;
    const complex_t IT_0856 = IT_0839*IT_0855;
    const complex_t IT_0857 = m_s*U_d1*conjq(V_cs)*e_em*IT_0035*U_su_15;
    const complex_t IT_0858 = IT_0048*IT_0857;
    const complex_t IT_0859 = 1.4142135623731*IT_0858;
    const complex_t IT_0860 = m_s*U_d1*conjq(V_ts)*e_em*IT_0035*U_su_25;
    const complex_t IT_0861 = IT_0048*IT_0860;
    const complex_t IT_0862 = 1.4142135623731*IT_0861;
    const complex_t IT_0863 = m_s*U_d1*V_us*e_em*IT_0035*U_su_05;
    const complex_t IT_0864 = IT_0048*IT_0863;
    const complex_t IT_0865 = 1.4142135623731*IT_0864;
    const complex_t IT_0866 = (complex_t{0, 1})*(IT_0859 + IT_0862 + IT_0865);
    const complex_t IT_0867 = 0.5*IT_0866;
    const complex_t IT_0868 = IT_0151*IT_0634*IT_0867;
    const complex_t IT_0869 = 0.101321183642338*IT_0868;
    const complex_t IT_0870 = mty::lt::C0iC(6, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0065, IT_0648, mty::lt::reg_int);
    const complex_t IT_0871 = m_b*IT_0870;
    const complex_t IT_0872 = mty::lt::C0iC(3, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0065, IT_0648, mty::lt::reg_int);
    const complex_t IT_0873 = m_b*IT_0872;
    const complex_t IT_0874 = mty::lt::C0iC(12, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0065, IT_0648, mty::lt::reg_int);
    const complex_t IT_0875 = m_b*IT_0874;
    const complex_t IT_0876 = mty::lt::C0iC(18, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0065, IT_0648, mty::lt::reg_int);
    const complex_t IT_0877 = m_b*IT_0876;
    const complex_t IT_0878 = IT_0871 + IT_0873 + IT_0875 + IT_0877;
    const complex_t IT_0879 = m_s*IT_0872;
    const complex_t IT_0880 = mty::lt::C0iC(15, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0065, IT_0648, mty::lt::reg_int);
    const complex_t IT_0881 = m_b*IT_0880;
    const complex_t IT_0882 = m_s*IT_0874;
    const complex_t IT_0883 = m_s*IT_0880;
    const complex_t IT_0884 = -IT_0879 + 2*IT_0881 + -IT_0882 + -IT_0883;
    const complex_t IT_0885 = IT_0878 + IT_0884;
    const complex_t IT_0886 = IT_0869*IT_0885;
    const complex_t IT_0887 = IT_0089*IT_0100*IT_0151;
    const complex_t IT_0888 = IT_0071*IT_0887;
    const complex_t IT_0889 = mty::lt::C0iC(0, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0065, IT_0016, mty::lt::reg_int);
    const complex_t IT_0890 = IT_0776 + IT_0778 + IT_0889;
    const complex_t IT_0891 = IT_0888*IT_0890;
    const complex_t IT_0892 = mty::lt::C0iC(0, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0065, IT_0250, mty::lt::reg_int);
    const complex_t IT_0893 = IT_0757 + IT_0759 + IT_0892;
    const complex_t IT_0894 = IT_0151*IT_0688*IT_0754;
    const complex_t IT_0895 = IT_0071*IT_0894;
    const complex_t IT_0896 = IT_0893*IT_0895;
    const complex_t IT_0897 = mty::lt::C0iC(0, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0065, IT_0606, mty::lt::reg_int);
    const complex_t IT_0898 = IT_0840 + IT_0842 + IT_0897;
    const complex_t IT_0899 = IT_0151*IT_0717*IT_0837;
    const complex_t IT_0900 = IT_0071*IT_0899;
    const complex_t IT_0901 = IT_0898*IT_0900;
    const complex_t IT_0902 = mty::lt::C0iC(0, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0065, IT_0065, IT_0648, mty::lt::reg_int);
    const complex_t IT_0903 = IT_0870 + IT_0872 + IT_0902;
    const complex_t IT_0904 = IT_0151*IT_0738*IT_0867;
    const complex_t IT_0905 = IT_0071*IT_0904;
    const complex_t IT_0906 = IT_0903*IT_0905;
    const complex_t IT_0907 = V_us*e_em*conjq(V_Wp1)*U_su_03;
    const complex_t IT_0908 = IT_0025*IT_0907;
    const complex_t IT_0909 = conjq(V_cs)*e_em*conjq(V_Wp1)*U_su_13;
    const complex_t IT_0910 = IT_0025*IT_0909;
    const complex_t IT_0911 = conjq(V_ts)*e_em*conjq(V_Wp1)*U_su_23;
    const complex_t IT_0912 = IT_0025*IT_0911;
    const complex_t IT_0913 = m_u*conjq(V_u1)*V_us*e_em*IT_0035*U_su_33;
    const complex_t IT_0914 = IT_0034*IT_0913;
    const complex_t IT_0915 = 1.4142135623731*IT_0914;
    const complex_t IT_0916 = m_c*conjq(V_cs)*conjq(V_u1)*e_em*IT_0035*U_su_43;
    const complex_t IT_0917 = IT_0034*IT_0916;
    const complex_t IT_0918 = 1.4142135623731*IT_0917;
    const complex_t IT_0919 = m_t*conjq(V_ts)*conjq(V_u1)*e_em*IT_0035*U_su_53;
    const complex_t IT_0920 = IT_0034*IT_0919;
    const complex_t IT_0921 = 1.4142135623731*IT_0920;
    const complex_t IT_0922 = (complex_t{0, 1})*(IT_0908 + IT_0910 + IT_0912 +
       (-0.5)*IT_0915 + (-0.5)*IT_0918 + (-0.5)*IT_0921);
    const complex_t IT_0923 = IT_0532*IT_0922;
    const complex_t IT_0924 = IT_0071*IT_0923;
    const complex_t IT_0925 = mty::lt::C0iC(0, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0065, IT_0140, IT_0140, mty::lt::reg_int);
    const complex_t IT_0926 = IT_0012*IT_0925;
    const complex_t IT_0927 = mty::lt::C0iC(3, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0065, IT_0140, IT_0140, mty::lt::reg_int);
    const complex_t IT_0928 = IT_0019*IT_0927;
    const complex_t IT_0929 = IT_0926 + IT_0928;
    const complex_t IT_0930 = IT_0924*IT_0929;
    const complex_t IT_0931 = IT_0178*IT_0532;
    const complex_t IT_0932 = 0.101321183642338*IT_0931;
    const complex_t IT_0933 = IT_0012*IT_0927;
    const complex_t IT_0934 = m_b*IT_0933;
    const complex_t IT_0935 = mty::lt::C0iC(12, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0065, IT_0140, IT_0140, mty::lt::reg_int);
    const complex_t IT_0936 = IT_0019*IT_0935;
    const complex_t IT_0937 = m_b*IT_0936;
    const complex_t IT_0938 = mty::lt::C0iC(6, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0065, IT_0140, IT_0140, mty::lt::reg_int);
    const complex_t IT_0939 = IT_0012*IT_0938;
    const complex_t IT_0940 = m_s*IT_0939;
    const complex_t IT_0941 = mty::lt::C0iC(15, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0065, IT_0140, IT_0140, mty::lt::reg_int);
    const complex_t IT_0942 = IT_0019*IT_0941;
    const complex_t IT_0943 = m_s*IT_0942;
    const complex_t IT_0944 = IT_0934 + IT_0937 + IT_0940 + IT_0943;
    const complex_t IT_0945 = IT_0932*IT_0944;
    const complex_t IT_0946 = m_b*conjq(U_d2)*V_cb*e_em*IT_0035*conjq(U_su_13);
    const complex_t IT_0947 = IT_0048*IT_0946;
    const complex_t IT_0948 = 1.4142135623731*IT_0947;
    const complex_t IT_0949 = m_b*conjq(U_d2)*V_tb*e_em*IT_0035*conjq(U_su_23);
    const complex_t IT_0950 = IT_0048*IT_0949;
    const complex_t IT_0951 = 1.4142135623731*IT_0950;
    const complex_t IT_0952 = m_b*conjq(U_d2)*e_em*IT_0035*conjq(U_su_03)
      *V_ub_mod;
    const complex_t IT_0953 = IT_0056*IT_0952;
    const complex_t IT_0954 = 1.4142135623731*IT_0953;
    const complex_t IT_0955 = (complex_t{0, 1})*(IT_0948 + IT_0951 + IT_0954);
    const complex_t IT_0956 = 0.5*IT_0955;
    const complex_t IT_0957 = V_us*e_em*conjq(V_Wp2)*U_su_03;
    const complex_t IT_0958 = IT_0025*IT_0957;
    const complex_t IT_0959 = conjq(V_cs)*e_em*conjq(V_Wp2)*U_su_13;
    const complex_t IT_0960 = IT_0025*IT_0959;
    const complex_t IT_0961 = conjq(V_ts)*e_em*conjq(V_Wp2)*U_su_23;
    const complex_t IT_0962 = IT_0025*IT_0961;
    const complex_t IT_0963 = m_u*conjq(V_u2)*V_us*e_em*IT_0035*U_su_33;
    const complex_t IT_0964 = IT_0034*IT_0963;
    const complex_t IT_0965 = 1.4142135623731*IT_0964;
    const complex_t IT_0966 = m_c*conjq(V_cs)*conjq(V_u2)*e_em*IT_0035*U_su_43;
    const complex_t IT_0967 = IT_0034*IT_0966;
    const complex_t IT_0968 = 1.4142135623731*IT_0967;
    const complex_t IT_0969 = m_t*conjq(V_ts)*conjq(V_u2)*e_em*IT_0035*U_su_53;
    const complex_t IT_0970 = IT_0034*IT_0969;
    const complex_t IT_0971 = 1.4142135623731*IT_0970;
    const complex_t IT_0972 = (complex_t{0, 1})*(IT_0958 + IT_0960 + IT_0962 +
       (-0.5)*IT_0965 + (-0.5)*IT_0968 + (-0.5)*IT_0971);
    const complex_t IT_0973 = IT_0956*IT_0972;
    const complex_t IT_0974 = IT_0023*IT_0973;
    const complex_t IT_0975 = mty::lt::C0iC(0, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0015, IT_0140, IT_0140, mty::lt::reg_int);
    const complex_t IT_0976 = IT_0012*IT_0975;
    const complex_t IT_0977 = mty::lt::C0iC(3, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0015, IT_0140, IT_0140, mty::lt::reg_int);
    const complex_t IT_0978 = IT_0019*IT_0977;
    const complex_t IT_0979 = IT_0976 + IT_0978;
    const complex_t IT_0980 = IT_0974*IT_0979;
    const complex_t IT_0981 = IT_0543*IT_0956;
    const complex_t IT_0982 = 0.101321183642338*IT_0981;
    const complex_t IT_0983 = IT_0012*IT_0977;
    const complex_t IT_0984 = m_b*IT_0983;
    const complex_t IT_0985 = mty::lt::C0iC(12, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0015, IT_0140, IT_0140, mty::lt::reg_int);
    const complex_t IT_0986 = IT_0019*IT_0985;
    const complex_t IT_0987 = m_b*IT_0986;
    const complex_t IT_0988 = mty::lt::C0iC(6, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0015, IT_0140, IT_0140, mty::lt::reg_int);
    const complex_t IT_0989 = IT_0012*IT_0988;
    const complex_t IT_0990 = m_s*IT_0989;
    const complex_t IT_0991 = mty::lt::C0iC(15, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0015, IT_0140, IT_0140, mty::lt::reg_int);
    const complex_t IT_0992 = IT_0019*IT_0991;
    const complex_t IT_0993 = m_s*IT_0992;
    const complex_t IT_0994 = IT_0984 + IT_0987 + IT_0990 + IT_0993;
    const complex_t IT_0995 = IT_0982*IT_0994;
    const complex_t IT_0996 = IT_0167*IT_0922;
    const complex_t IT_0997 = 0.101321183642338*IT_0996;
    const complex_t IT_0998 = IT_0944*IT_0997;
    const complex_t IT_0999 = V_cb*e_em*V_Wp2*conjq(U_su_13);
    const complex_t IT_1000 = IT_0025*IT_0999;
    const complex_t IT_1001 = V_tb*e_em*V_Wp2*conjq(U_su_23);
    const complex_t IT_1002 = IT_0025*IT_1001;
    const complex_t IT_1003 = e_em*V_Wp2*conjq(U_su_03)*V_ub_mod;
    const complex_t IT_1004 = IT_0076*IT_1003;
    const complex_t IT_1005 = m_c*V_cb*V_u2*e_em*IT_0035*conjq(U_su_43);
    const complex_t IT_1006 = IT_0034*IT_1005;
    const complex_t IT_1007 = 1.4142135623731*IT_1006;
    const complex_t IT_1008 = m_t*V_tb*V_u2*e_em*IT_0035*conjq(U_su_53);
    const complex_t IT_1009 = IT_0034*IT_1008;
    const complex_t IT_1010 = 1.4142135623731*IT_1009;
    const complex_t IT_1011 = m_u*V_u2*e_em*IT_0035*conjq(U_su_33)*V_ub_mod;
    const complex_t IT_1012 = IT_0085*IT_1011;
    const complex_t IT_1013 = 1.4142135623731*IT_1012;
    const complex_t IT_1014 = (complex_t{0, 1})*(IT_1000 + IT_1002 + IT_1004 +
       (-0.5)*IT_1007 + (-0.5)*IT_1010 + (-0.5)*IT_1013);
    const complex_t IT_1015 = IT_0972*IT_1014;
    const complex_t IT_1016 = 0.101321183642338*IT_1015;
    const complex_t IT_1017 = IT_0994*IT_1016;
    const complex_t IT_1018 = IT_0167*IT_0178;
    const complex_t IT_1019 = IT_0071*IT_1018;
    const complex_t IT_1020 = IT_0929*IT_1019;
    const complex_t IT_1021 = IT_0543*IT_1014;
    const complex_t IT_1022 = IT_0023*IT_1021;
    const complex_t IT_1023 = IT_0979*IT_1022;
    const complex_t IT_1024 = V_us*e_em*conjq(V_Wp1)*U_su_01;
    const complex_t IT_1025 = IT_0025*IT_1024;
    const complex_t IT_1026 = conjq(V_cs)*e_em*conjq(V_Wp1)*U_su_11;
    const complex_t IT_1027 = IT_0025*IT_1026;
    const complex_t IT_1028 = conjq(V_ts)*e_em*conjq(V_Wp1)*U_su_21;
    const complex_t IT_1029 = IT_0025*IT_1028;
    const complex_t IT_1030 = m_u*conjq(V_u1)*V_us*e_em*IT_0035*U_su_31;
    const complex_t IT_1031 = IT_0034*IT_1030;
    const complex_t IT_1032 = 1.4142135623731*IT_1031;
    const complex_t IT_1033 = m_c*conjq(V_cs)*conjq(V_u1)*e_em*IT_0035*U_su_41;
    const complex_t IT_1034 = IT_0034*IT_1033;
    const complex_t IT_1035 = 1.4142135623731*IT_1034;
    const complex_t IT_1036 = m_t*conjq(V_ts)*conjq(V_u1)*e_em*IT_0035*U_su_51;
    const complex_t IT_1037 = IT_0034*IT_1036;
    const complex_t IT_1038 = 1.4142135623731*IT_1037;
    const complex_t IT_1039 = (complex_t{0, 1})*(IT_1025 + IT_1027 + IT_1029 +
       (-0.5)*IT_1032 + (-0.5)*IT_1035 + (-0.5)*IT_1038);
    const complex_t IT_1040 = IT_0472*IT_1039;
    const complex_t IT_1041 = IT_0071*IT_1040;
    const complex_t IT_1042 = mty::lt::C0iC(3, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0065, IT_0250, IT_0250, mty::lt::reg_int);
    const complex_t IT_1043 = IT_0019*IT_1042;
    const complex_t IT_1044 = mty::lt::C0iC(0, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0065, IT_0250, IT_0250, mty::lt::reg_int);
    const complex_t IT_1045 = IT_0012*IT_1044;
    const complex_t IT_1046 = IT_1043 + IT_1045;
    const complex_t IT_1047 = IT_1041*IT_1046;
    const complex_t IT_1048 = IT_0472*IT_0754;
    const complex_t IT_1049 = 0.101321183642338*IT_1048;
    const complex_t IT_1050 = IT_0012*IT_1042;
    const complex_t IT_1051 = m_b*IT_1050;
    const complex_t IT_1052 = mty::lt::C0iC(12, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0065, IT_0250, IT_0250, mty::lt::reg_int);
    const complex_t IT_1053 = IT_0019*IT_1052;
    const complex_t IT_1054 = m_b*IT_1053;
    const complex_t IT_1055 = mty::lt::C0iC(6, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0065, IT_0250, IT_0250, mty::lt::reg_int);
    const complex_t IT_1056 = IT_0012*IT_1055;
    const complex_t IT_1057 = m_s*IT_1056;
    const complex_t IT_1058 = mty::lt::C0iC(15, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0065, IT_0250, IT_0250, mty::lt::reg_int);
    const complex_t IT_1059 = IT_0019*IT_1058;
    const complex_t IT_1060 = m_s*IT_1059;
    const complex_t IT_1061 = IT_1051 + IT_1054 + IT_1057 + IT_1060;
    const complex_t IT_1062 = IT_1049*IT_1061;
    const complex_t IT_1063 = IT_0231*IT_0247;
    const complex_t IT_1064 = IT_0023*IT_1063;
    const complex_t IT_1065 = mty::lt::C0iC(0, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0015, IT_0250, IT_0250, mty::lt::reg_int);
    const complex_t IT_1066 = IT_0012*IT_1065;
    const complex_t IT_1067 = mty::lt::C0iC(3, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0015, IT_0250, IT_0250, mty::lt::reg_int);
    const complex_t IT_1068 = IT_0019*IT_1067;
    const complex_t IT_1069 = IT_1066 + IT_1068;
    const complex_t IT_1070 = IT_1064*IT_1069;
    const complex_t IT_1071 = IT_0231*IT_0483;
    const complex_t IT_1072 = 0.101321183642338*IT_1071;
    const complex_t IT_1073 = IT_0012*IT_1067;
    const complex_t IT_1074 = m_b*IT_1073;
    const complex_t IT_1075 = mty::lt::C0iC(6, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0015, IT_0250, IT_0250, mty::lt::reg_int);
    const complex_t IT_1076 = IT_0012*IT_1075;
    const complex_t IT_1077 = m_s*IT_1076;
    const complex_t IT_1078 = mty::lt::C0iC(12, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0015, IT_0250, IT_0250, mty::lt::reg_int);
    const complex_t IT_1079 = IT_0019*IT_1078;
    const complex_t IT_1080 = m_b*IT_1079;
    const complex_t IT_1081 = mty::lt::C0iC(15, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0015, IT_0250, IT_0250, mty::lt::reg_int);
    const complex_t IT_1082 = IT_0019*IT_1081;
    const complex_t IT_1083 = m_s*IT_1082;
    const complex_t IT_1084 = IT_1074 + IT_1077 + IT_1080 + IT_1083;
    const complex_t IT_1085 = IT_1072*IT_1084;
    const complex_t IT_1086 = IT_0688*IT_1039;
    const complex_t IT_1087 = 0.101321183642338*IT_1086;
    const complex_t IT_1088 = IT_1061*IT_1087;
    const complex_t IT_1089 = V_cb*e_em*V_Wp2*conjq(U_su_11);
    const complex_t IT_1090 = IT_0025*IT_1089;
    const complex_t IT_1091 = V_tb*e_em*V_Wp2*conjq(U_su_21);
    const complex_t IT_1092 = IT_0025*IT_1091;
    const complex_t IT_1093 = e_em*V_Wp2*conjq(U_su_01)*V_ub_mod;
    const complex_t IT_1094 = IT_0076*IT_1093;
    const complex_t IT_1095 = m_c*V_cb*V_u2*e_em*IT_0035*conjq(U_su_41);
    const complex_t IT_1096 = IT_0034*IT_1095;
    const complex_t IT_1097 = 1.4142135623731*IT_1096;
    const complex_t IT_1098 = m_t*V_tb*V_u2*e_em*IT_0035*conjq(U_su_51);
    const complex_t IT_1099 = IT_0034*IT_1098;
    const complex_t IT_1100 = 1.4142135623731*IT_1099;
    const complex_t IT_1101 = m_u*V_u2*e_em*IT_0035*conjq(U_su_31)*V_ub_mod;
    const complex_t IT_1102 = IT_0085*IT_1101;
    const complex_t IT_1103 = 1.4142135623731*IT_1102;
    const complex_t IT_1104 = (complex_t{0, 1})*(IT_1090 + IT_1092 + IT_1094 +
       (-0.5)*IT_1097 + (-0.5)*IT_1100 + (-0.5)*IT_1103);
    const complex_t IT_1105 = IT_0247*IT_1104;
    const complex_t IT_1106 = 0.101321183642338*IT_1105;
    const complex_t IT_1107 = IT_1084*IT_1106;
    const complex_t IT_1108 = IT_0688*IT_0754;
    const complex_t IT_1109 = IT_0071*IT_1108;
    const complex_t IT_1110 = IT_1046*IT_1109;
    const complex_t IT_1111 = IT_0483*IT_1104;
    const complex_t IT_1112 = IT_0023*IT_1111;
    const complex_t IT_1113 = IT_1069*IT_1112;
    const complex_t IT_1114 = mty::lt::C0iC(0, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0015, IT_0016, mty::lt::reg_int);
    const complex_t IT_1115 = mty::lt::C0iC(3, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0015, IT_0016, mty::lt::reg_int);
    const complex_t IT_1116 = mty::lt::C0iC(6, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0015, IT_0016, mty::lt::reg_int);
    const complex_t IT_1117 = IT_1114 + IT_1115 + IT_1116;
    const complex_t IT_1118 = IT_0045*IT_0061*IT_0220;
    const complex_t IT_1119 = IT_0023*IT_1118;
    const complex_t IT_1120 = IT_1117*IT_1119;
    const complex_t IT_1121 = IT_0220*IT_0956*IT_0972;
    const complex_t IT_1122 = IT_0023*IT_1121;
    const complex_t IT_1123 = mty::lt::C0iC(0, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0015, IT_0140, mty::lt::reg_int);
    const complex_t IT_1124 = mty::lt::C0iC(3, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0015, IT_0140, mty::lt::reg_int);
    const complex_t IT_1125 = mty::lt::C0iC(6, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0015, IT_0140, mty::lt::reg_int);
    const complex_t IT_1126 = IT_1123 + IT_1124 + IT_1125;
    const complex_t IT_1127 = IT_1122*IT_1126;
    const complex_t IT_1128 = IT_0220*IT_0412*IT_0428;
    const complex_t IT_1129 = IT_0023*IT_1128;
    const complex_t IT_1130 = mty::lt::C0iC(0, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0015, IT_0133, mty::lt::reg_int);
    const complex_t IT_1131 = mty::lt::C0iC(3, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0015, IT_0133, mty::lt::reg_int);
    const complex_t IT_1132 = mty::lt::C0iC(6, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0015, IT_0133, mty::lt::reg_int);
    const complex_t IT_1133 = IT_1130 + IT_1131 + IT_1132;
    const complex_t IT_1134 = IT_1129*IT_1133;
    const complex_t IT_1135 = mty::lt::C0iC(0, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0015, IT_0606, mty::lt::reg_int);
    const complex_t IT_1136 = mty::lt::C0iC(3, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0015, IT_0606, mty::lt::reg_int);
    const complex_t IT_1137 = mty::lt::C0iC(6, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0015, IT_0606, mty::lt::reg_int);
    const complex_t IT_1138 = IT_1135 + IT_1136 + IT_1137;
    const complex_t IT_1139 = V_us*e_em*conjq(V_Wp2)*U_su_04;
    const complex_t IT_1140 = IT_0025*IT_1139;
    const complex_t IT_1141 = conjq(V_cs)*e_em*conjq(V_Wp2)*U_su_14;
    const complex_t IT_1142 = IT_0025*IT_1141;
    const complex_t IT_1143 = conjq(V_ts)*e_em*conjq(V_Wp2)*U_su_24;
    const complex_t IT_1144 = IT_0025*IT_1143;
    const complex_t IT_1145 = m_u*conjq(V_u2)*V_us*e_em*IT_0035*U_su_34;
    const complex_t IT_1146 = IT_0034*IT_1145;
    const complex_t IT_1147 = 1.4142135623731*IT_1146;
    const complex_t IT_1148 = m_c*conjq(V_cs)*conjq(V_u2)*e_em*IT_0035*U_su_44;
    const complex_t IT_1149 = IT_0034*IT_1148;
    const complex_t IT_1150 = 1.4142135623731*IT_1149;
    const complex_t IT_1151 = m_t*conjq(V_ts)*conjq(V_u2)*e_em*IT_0035*U_su_54;
    const complex_t IT_1152 = IT_0034*IT_1151;
    const complex_t IT_1153 = 1.4142135623731*IT_1152;
    const complex_t IT_1154 = (complex_t{0, 1})*(IT_1140 + IT_1142 + IT_1144 +
       (-0.5)*IT_1147 + (-0.5)*IT_1150 + (-0.5)*IT_1153);
    const complex_t IT_1155 = m_b*conjq(U_d2)*V_cb*e_em*IT_0035*conjq(U_su_14);
    const complex_t IT_1156 = IT_0048*IT_1155;
    const complex_t IT_1157 = 1.4142135623731*IT_1156;
    const complex_t IT_1158 = m_b*conjq(U_d2)*V_tb*e_em*IT_0035*conjq(U_su_24);
    const complex_t IT_1159 = IT_0048*IT_1158;
    const complex_t IT_1160 = 1.4142135623731*IT_1159;
    const complex_t IT_1161 = m_b*conjq(U_d2)*e_em*IT_0035*conjq(U_su_04)
      *V_ub_mod;
    const complex_t IT_1162 = IT_0056*IT_1161;
    const complex_t IT_1163 = 1.4142135623731*IT_1162;
    const complex_t IT_1164 = (complex_t{0, 1})*(IT_1157 + IT_1160 + IT_1163);
    const complex_t IT_1165 = 0.5*IT_1164;
    const complex_t IT_1166 = IT_0220*IT_1154*IT_1165;
    const complex_t IT_1167 = IT_0023*IT_1166;
    const complex_t IT_1168 = IT_1138*IT_1167;
    const complex_t IT_1169 = m_b*conjq(U_d2)*V_cb*e_em*IT_0035*conjq(U_su_15);
    const complex_t IT_1170 = IT_0048*IT_1169;
    const complex_t IT_1171 = 1.4142135623731*IT_1170;
    const complex_t IT_1172 = m_b*conjq(U_d2)*V_tb*e_em*IT_0035*conjq(U_su_25);
    const complex_t IT_1173 = IT_0048*IT_1172;
    const complex_t IT_1174 = 1.4142135623731*IT_1173;
    const complex_t IT_1175 = m_b*conjq(U_d2)*e_em*IT_0035*conjq(U_su_05)
      *V_ub_mod;
    const complex_t IT_1176 = IT_0056*IT_1175;
    const complex_t IT_1177 = 1.4142135623731*IT_1176;
    const complex_t IT_1178 = (complex_t{0, 1})*(IT_1171 + IT_1174 + IT_1177);
    const complex_t IT_1179 = 0.5*IT_1178;
    const complex_t IT_1180 = V_us*e_em*conjq(V_Wp2)*U_su_05;
    const complex_t IT_1181 = IT_0025*IT_1180;
    const complex_t IT_1182 = conjq(V_cs)*e_em*conjq(V_Wp2)*U_su_15;
    const complex_t IT_1183 = IT_0025*IT_1182;
    const complex_t IT_1184 = conjq(V_ts)*e_em*conjq(V_Wp2)*U_su_25;
    const complex_t IT_1185 = IT_0025*IT_1184;
    const complex_t IT_1186 = m_u*conjq(V_u2)*V_us*e_em*IT_0035*U_su_35;
    const complex_t IT_1187 = IT_0034*IT_1186;
    const complex_t IT_1188 = 1.4142135623731*IT_1187;
    const complex_t IT_1189 = m_c*conjq(V_cs)*conjq(V_u2)*e_em*IT_0035*U_su_45;
    const complex_t IT_1190 = IT_0034*IT_1189;
    const complex_t IT_1191 = 1.4142135623731*IT_1190;
    const complex_t IT_1192 = m_t*conjq(V_ts)*conjq(V_u2)*e_em*IT_0035*U_su_55;
    const complex_t IT_1193 = IT_0034*IT_1192;
    const complex_t IT_1194 = 1.4142135623731*IT_1193;
    const complex_t IT_1195 = (complex_t{0, 1})*(IT_1181 + IT_1183 + IT_1185 +
       (-0.5)*IT_1188 + (-0.5)*IT_1191 + (-0.5)*IT_1194);
    const complex_t IT_1196 = IT_0220*IT_1179*IT_1195;
    const complex_t IT_1197 = IT_0023*IT_1196;
    const complex_t IT_1198 = mty::lt::C0iC(0, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0015, IT_0648, mty::lt::reg_int);
    const complex_t IT_1199 = mty::lt::C0iC(3, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0015, IT_0648, mty::lt::reg_int);
    const complex_t IT_1200 = mty::lt::C0iC(6, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0015, IT_0648, mty::lt::reg_int);
    const complex_t IT_1201 = IT_1198 + IT_1199 + IT_1200;
    const complex_t IT_1202 = IT_1197*IT_1201;
    const complex_t IT_1203 = IT_0220*IT_0247*IT_1104;
    const complex_t IT_1204 = 0.101321183642338*IT_1203;
    const complex_t IT_1205 = mty::lt::C0iC(18, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0015, IT_0250, mty::lt::reg_int);
    const complex_t IT_1206 = m_b*IT_1205;
    const complex_t IT_1207 = m_b*IT_0253;
    const complex_t IT_1208 = m_b*IT_0252;
    const complex_t IT_1209 = mty::lt::C0iC(12, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0015, IT_0250, mty::lt::reg_int);
    const complex_t IT_1210 = m_b*IT_1209;
    const complex_t IT_1211 = IT_1206 + IT_1207 + IT_1208 + IT_1210;
    const complex_t IT_1212 = mty::lt::C0iC(15, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0015, IT_0250, mty::lt::reg_int);
    const complex_t IT_1213 = m_s*IT_1212;
    const complex_t IT_1214 = m_b*IT_1212;
    const complex_t IT_1215 = m_s*IT_1209;
    const complex_t IT_1216 = m_s*IT_0252;
    const complex_t IT_1217 = -IT_1213 + 2*IT_1214 + -IT_1215 + -IT_1216;
    const complex_t IT_1218 = IT_1211 + IT_1217;
    const complex_t IT_1219 = IT_1204*IT_1218;
    const complex_t IT_1220 = IT_0119*IT_0220*IT_0428;
    const complex_t IT_1221 = 0.101321183642338*IT_1220;
    const complex_t IT_1222 = m_b*IT_1132;
    const complex_t IT_1223 = m_b*IT_1131;
    const complex_t IT_1224 = mty::lt::C0iC(18, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0015, IT_0133, mty::lt::reg_int);
    const complex_t IT_1225 = m_b*IT_1224;
    const complex_t IT_1226 = mty::lt::C0iC(12, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0015, IT_0133, mty::lt::reg_int);
    const complex_t IT_1227 = m_b*IT_1226;
    const complex_t IT_1228 = IT_1222 + IT_1223 + IT_1225 + IT_1227;
    const complex_t IT_1229 = m_s*IT_1131;
    const complex_t IT_1230 = mty::lt::C0iC(15, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0015, IT_0133, mty::lt::reg_int);
    const complex_t IT_1231 = m_s*IT_1230;
    const complex_t IT_1232 = m_b*IT_1230;
    const complex_t IT_1233 = m_s*IT_1226;
    const complex_t IT_1234 = -IT_1229 + -IT_1231 + 2*IT_1232 + -IT_1233;
    const complex_t IT_1235 = IT_1228 + IT_1234;
    const complex_t IT_1236 = IT_1221*IT_1235;
    const complex_t IT_1237 = IT_0220*IT_0972*IT_1014;
    const complex_t IT_1238 = 0.101321183642338*IT_1237;
    const complex_t IT_1239 = mty::lt::C0iC(12, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0015, IT_0140, mty::lt::reg_int);
    const complex_t IT_1240 = m_b*IT_1239;
    const complex_t IT_1241 = mty::lt::C0iC(18, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0015, IT_0140, mty::lt::reg_int);
    const complex_t IT_1242 = m_b*IT_1241;
    const complex_t IT_1243 = m_b*IT_1125;
    const complex_t IT_1244 = m_b*IT_1124;
    const complex_t IT_1245 = IT_1240 + IT_1242 + IT_1243 + IT_1244;
    const complex_t IT_1246 = m_s*IT_1124;
    const complex_t IT_1247 = m_s*IT_1239;
    const complex_t IT_1248 = mty::lt::C0iC(15, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0015, IT_0140, mty::lt::reg_int);
    const complex_t IT_1249 = m_s*IT_1248;
    const complex_t IT_1250 = m_b*IT_1248;
    const complex_t IT_1251 = -IT_1246 + -IT_1247 + -IT_1249 + 2*IT_1250;
    const complex_t IT_1252 = IT_1245 + IT_1251;
    const complex_t IT_1253 = IT_1238*IT_1252;
    const complex_t IT_1254 = IT_0045*IT_0220*IT_0345;
    const complex_t IT_1255 = 0.101321183642338*IT_1254;
    const complex_t IT_1256 = mty::lt::C0iC(12, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0015, IT_0016, mty::lt::reg_int);
    const complex_t IT_1257 = m_b*IT_1256;
    const complex_t IT_1258 = mty::lt::C0iC(18, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0015, IT_0016, mty::lt::reg_int);
    const complex_t IT_1259 = m_b*IT_1258;
    const complex_t IT_1260 = m_b*IT_1116;
    const complex_t IT_1261 = m_b*IT_1115;
    const complex_t IT_1262 = IT_1257 + IT_1259 + IT_1260 + IT_1261;
    const complex_t IT_1263 = m_s*IT_1256;
    const complex_t IT_1264 = mty::lt::C0iC(15, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0015, IT_0016, mty::lt::reg_int);
    const complex_t IT_1265 = m_s*IT_1264;
    const complex_t IT_1266 = m_b*IT_1264;
    const complex_t IT_1267 = m_s*IT_1115;
    const complex_t IT_1268 = -IT_1263 + -IT_1265 + 2*IT_1266 + -IT_1267;
    const complex_t IT_1269 = IT_1262 + IT_1268;
    const complex_t IT_1270 = IT_1255*IT_1269;
    const complex_t IT_1271 = V_cb*e_em*V_Wp2*conjq(U_su_14);
    const complex_t IT_1272 = IT_0025*IT_1271;
    const complex_t IT_1273 = V_tb*e_em*V_Wp2*conjq(U_su_24);
    const complex_t IT_1274 = IT_0025*IT_1273;
    const complex_t IT_1275 = e_em*V_Wp2*conjq(U_su_04)*V_ub_mod;
    const complex_t IT_1276 = IT_0076*IT_1275;
    const complex_t IT_1277 = m_c*V_cb*V_u2*e_em*IT_0035*conjq(U_su_44);
    const complex_t IT_1278 = IT_0034*IT_1277;
    const complex_t IT_1279 = 1.4142135623731*IT_1278;
    const complex_t IT_1280 = m_t*V_tb*V_u2*e_em*IT_0035*conjq(U_su_54);
    const complex_t IT_1281 = IT_0034*IT_1280;
    const complex_t IT_1282 = 1.4142135623731*IT_1281;
    const complex_t IT_1283 = m_u*V_u2*e_em*IT_0035*conjq(U_su_34)*V_ub_mod;
    const complex_t IT_1284 = IT_0085*IT_1283;
    const complex_t IT_1285 = 1.4142135623731*IT_1284;
    const complex_t IT_1286 = (complex_t{0, 1})*(IT_1272 + IT_1274 + IT_1276 +
       (-0.5)*IT_1279 + (-0.5)*IT_1282 + (-0.5)*IT_1285);
    const complex_t IT_1287 = IT_0220*IT_1154*IT_1286;
    const complex_t IT_1288 = 0.101321183642338*IT_1287;
    const complex_t IT_1289 = m_b*IT_1137;
    const complex_t IT_1290 = m_b*IT_1136;
    const complex_t IT_1291 = mty::lt::C0iC(12, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0015, IT_0606, mty::lt::reg_int);
    const complex_t IT_1292 = m_b*IT_1291;
    const complex_t IT_1293 = mty::lt::C0iC(18, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0015, IT_0606, mty::lt::reg_int);
    const complex_t IT_1294 = m_b*IT_1293;
    const complex_t IT_1295 = IT_1289 + IT_1290 + IT_1292 + IT_1294;
    const complex_t IT_1296 = m_s*IT_1136;
    const complex_t IT_1297 = m_s*IT_1291;
    const complex_t IT_1298 = mty::lt::C0iC(15, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0015, IT_0606, mty::lt::reg_int);
    const complex_t IT_1299 = m_s*IT_1298;
    const complex_t IT_1300 = m_b*IT_1298;
    const complex_t IT_1301 = -IT_1296 + -IT_1297 + -IT_1299 + 2*IT_1300;
    const complex_t IT_1302 = IT_1295 + IT_1301;
    const complex_t IT_1303 = IT_1288*IT_1302;
    const complex_t IT_1304 = V_cb*e_em*V_Wp2*conjq(U_su_15);
    const complex_t IT_1305 = IT_0025*IT_1304;
    const complex_t IT_1306 = V_tb*e_em*V_Wp2*conjq(U_su_25);
    const complex_t IT_1307 = IT_0025*IT_1306;
    const complex_t IT_1308 = e_em*V_Wp2*conjq(U_su_05)*V_ub_mod;
    const complex_t IT_1309 = IT_0076*IT_1308;
    const complex_t IT_1310 = m_c*V_cb*V_u2*e_em*IT_0035*conjq(U_su_45);
    const complex_t IT_1311 = IT_0034*IT_1310;
    const complex_t IT_1312 = 1.4142135623731*IT_1311;
    const complex_t IT_1313 = m_t*V_tb*V_u2*e_em*IT_0035*conjq(U_su_55);
    const complex_t IT_1314 = IT_0034*IT_1313;
    const complex_t IT_1315 = 1.4142135623731*IT_1314;
    const complex_t IT_1316 = m_u*V_u2*e_em*IT_0035*conjq(U_su_35)*V_ub_mod;
    const complex_t IT_1317 = IT_0085*IT_1316;
    const complex_t IT_1318 = 1.4142135623731*IT_1317;
    const complex_t IT_1319 = (complex_t{0, 1})*(IT_1305 + IT_1307 + IT_1309 +
       (-0.5)*IT_1312 + (-0.5)*IT_1315 + (-0.5)*IT_1318);
    const complex_t IT_1320 = IT_0220*IT_1195*IT_1319;
    const complex_t IT_1321 = 0.101321183642338*IT_1320;
    const complex_t IT_1322 = m_b*IT_1199;
    const complex_t IT_1323 = mty::lt::C0iC(18, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0015, IT_0648, mty::lt::reg_int);
    const complex_t IT_1324 = m_b*IT_1323;
    const complex_t IT_1325 = m_b*IT_1200;
    const complex_t IT_1326 = mty::lt::C0iC(12, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0015, IT_0648, mty::lt::reg_int);
    const complex_t IT_1327 = m_b*IT_1326;
    const complex_t IT_1328 = IT_1322 + IT_1324 + IT_1325 + IT_1327;
    const complex_t IT_1329 = m_s*IT_1199;
    const complex_t IT_1330 = m_s*IT_1326;
    const complex_t IT_1331 = mty::lt::C0iC(15, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0015, IT_0648, mty::lt::reg_int);
    const complex_t IT_1332 = m_s*IT_1331;
    const complex_t IT_1333 = m_b*IT_1331;
    const complex_t IT_1334 = -IT_1329 + -IT_1330 + -IT_1332 + 2*IT_1333;
    const complex_t IT_1335 = IT_1328 + IT_1334;
    const complex_t IT_1336 = IT_1321*IT_1335;
    const complex_t IT_1337 = e_em*conjq(V_Wp1)*V_Wp2;
    const complex_t IT_1338 = IT_0145*IT_1337;
    const complex_t IT_1339 = conjq(V_u1)*V_u2*e_em;
    const complex_t IT_1340 = IT_0145*IT_1339;
    const complex_t IT_1341 = (complex_t{0, 1})*(IT_1338 + IT_1340);
    const complex_t IT_1342 = IT_0247*IT_0472*IT_1341;
    const complex_t IT_1343 = IT_0071*IT_1342;
    const complex_t IT_1344 = IT_0672*IT_1343;
    const complex_t IT_1345 = IT_0045*IT_0282*IT_1341;
    const complex_t IT_1346 = IT_0071*IT_1345;
    const complex_t IT_1347 = IT_0669*IT_1346;
    const complex_t IT_1348 = IT_0532*IT_0972*IT_1341;
    const complex_t IT_1349 = IT_0071*IT_1348;
    const complex_t IT_1350 = IT_0695*IT_1349;
    const complex_t IT_1351 = IT_0367*IT_0428*IT_1341;
    const complex_t IT_1352 = IT_0071*IT_1351;
    const complex_t IT_1353 = IT_0700*IT_1352;
    const complex_t IT_1354 = IT_0603*IT_1154*IT_1341;
    const complex_t IT_1355 = IT_0071*IT_1354;
    const complex_t IT_1356 = IT_0721*IT_1355;
    const complex_t IT_1357 = IT_0634*IT_1195*IT_1341;
    const complex_t IT_1358 = IT_0071*IT_1357;
    const complex_t IT_1359 = IT_0742*IT_1358;
    const complex_t IT_1360 = IT_0247*IT_0688*IT_1341;
    const complex_t IT_1361 = 0.101321183642338*IT_1360;
    const complex_t IT_1362 = IT_0501*IT_1361;
    const complex_t IT_1363 = IT_0201*IT_0428*IT_1341;
    const complex_t IT_1364 = 0.101321183642338*IT_1363;
    const complex_t IT_1365 = IT_0580*IT_1364;
    const complex_t IT_1366 = IT_0167*IT_0972*IT_1341;
    const complex_t IT_1367 = 0.101321183642338*IT_1366;
    const complex_t IT_1368 = IT_0561*IT_1367;
    const complex_t IT_1369 = IT_0045*IT_0089*IT_1341;
    const complex_t IT_1370 = 0.101321183642338*IT_1369;
    const complex_t IT_1371 = IT_0520*IT_1370;
    const complex_t IT_1372 = IT_0717*IT_1154*IT_1341;
    const complex_t IT_1373 = 0.101321183642338*IT_1372;
    const complex_t IT_1374 = IT_0622*IT_1373;
    const complex_t IT_1375 = IT_0738*IT_1195*IT_1341;
    const complex_t IT_1376 = 0.101321183642338*IT_1375;
    const complex_t IT_1377 = IT_0664*IT_1376;
    const complex_t IT_1378 = e_em*V_Wp1*conjq(V_Wp2);
    const complex_t IT_1379 = IT_0145*IT_1378;
    const complex_t IT_1380 = V_u1*conjq(V_u2)*e_em;
    const complex_t IT_1381 = IT_0145*IT_1380;
    const complex_t IT_1382 = (complex_t{0, 1})*(IT_1379 + IT_1381);
    const complex_t IT_1383 = IT_0231*IT_1039*IT_1382;
    const complex_t IT_1384 = IT_0023*IT_1383;
    const complex_t IT_1385 = mty::lt::C0iC(0, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0065, IT_0250, mty::lt::reg_int);
    const complex_t IT_1386 = mty::lt::C0iC(3, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0065, IT_0250, mty::lt::reg_int);
    const complex_t IT_1387 = mty::lt::C0iC(6, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0065, IT_0250, mty::lt::reg_int);
    const complex_t IT_1388 = IT_1385 + IT_1386 + IT_1387;
    const complex_t IT_1389 = IT_1384*IT_1388;
    const complex_t IT_1390 = mty::lt::C0iC(0, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0065, IT_0016, mty::lt::reg_int);
    const complex_t IT_1391 = mty::lt::C0iC(3, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0065, IT_0016, mty::lt::reg_int);
    const complex_t IT_1392 = mty::lt::C0iC(6, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0065, IT_0016, mty::lt::reg_int);
    const complex_t IT_1393 = IT_1390 + IT_1391 + IT_1392;
    const complex_t IT_1394 = IT_0061*IT_0271*IT_1382;
    const complex_t IT_1395 = IT_0023*IT_1394;
    const complex_t IT_1396 = IT_1393*IT_1395;
    const complex_t IT_1397 = IT_0922*IT_0956*IT_1382;
    const complex_t IT_1398 = IT_0023*IT_1397;
    const complex_t IT_1399 = mty::lt::C0iC(0, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0065, IT_0140, mty::lt::reg_int);
    const complex_t IT_1400 = mty::lt::C0iC(3, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0065, IT_0140, mty::lt::reg_int);
    const complex_t IT_1401 = mty::lt::C0iC(6, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0065, IT_0140, mty::lt::reg_int);
    const complex_t IT_1402 = IT_1399 + IT_1400 + IT_1401;
    const complex_t IT_1403 = IT_1398*IT_1402;
    const complex_t IT_1404 = IT_0383*IT_0412*IT_1382;
    const complex_t IT_1405 = IT_0023*IT_1404;
    const complex_t IT_1406 = mty::lt::C0iC(0, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0065, IT_0133, mty::lt::reg_int);
    const complex_t IT_1407 = mty::lt::C0iC(3, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0065, IT_0133, mty::lt::reg_int);
    const complex_t IT_1408 = mty::lt::C0iC(6, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0065, IT_0133, mty::lt::reg_int);
    const complex_t IT_1409 = IT_1406 + IT_1407 + IT_1408;
    const complex_t IT_1410 = IT_1405*IT_1409;
    const complex_t IT_1411 = mty::lt::C0iC(0, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0065, IT_0606, mty::lt::reg_int);
    const complex_t IT_1412 = mty::lt::C0iC(3, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0065, IT_0606, mty::lt::reg_int);
    const complex_t IT_1413 = mty::lt::C0iC(6, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0065, IT_0606, mty::lt::reg_int);
    const complex_t IT_1414 = IT_1411 + IT_1412 + IT_1413;
    const complex_t IT_1415 = V_us*e_em*conjq(V_Wp1)*U_su_04;
    const complex_t IT_1416 = IT_0025*IT_1415;
    const complex_t IT_1417 = conjq(V_cs)*e_em*conjq(V_Wp1)*U_su_14;
    const complex_t IT_1418 = IT_0025*IT_1417;
    const complex_t IT_1419 = conjq(V_ts)*e_em*conjq(V_Wp1)*U_su_24;
    const complex_t IT_1420 = IT_0025*IT_1419;
    const complex_t IT_1421 = m_u*conjq(V_u1)*V_us*e_em*IT_0035*U_su_34;
    const complex_t IT_1422 = IT_0034*IT_1421;
    const complex_t IT_1423 = 1.4142135623731*IT_1422;
    const complex_t IT_1424 = m_c*conjq(V_cs)*conjq(V_u1)*e_em*IT_0035*U_su_44;
    const complex_t IT_1425 = IT_0034*IT_1424;
    const complex_t IT_1426 = 1.4142135623731*IT_1425;
    const complex_t IT_1427 = m_t*conjq(V_ts)*conjq(V_u1)*e_em*IT_0035*U_su_54;
    const complex_t IT_1428 = IT_0034*IT_1427;
    const complex_t IT_1429 = 1.4142135623731*IT_1428;
    const complex_t IT_1430 = (complex_t{0, 1})*(IT_1416 + IT_1418 + IT_1420 +
       (-0.5)*IT_1423 + (-0.5)*IT_1426 + (-0.5)*IT_1429);
    const complex_t IT_1431 = IT_1165*IT_1382*IT_1430;
    const complex_t IT_1432 = IT_0023*IT_1431;
    const complex_t IT_1433 = IT_1414*IT_1432;
    const complex_t IT_1434 = V_us*e_em*conjq(V_Wp1)*U_su_05;
    const complex_t IT_1435 = IT_0025*IT_1434;
    const complex_t IT_1436 = conjq(V_cs)*e_em*conjq(V_Wp1)*U_su_15;
    const complex_t IT_1437 = IT_0025*IT_1436;
    const complex_t IT_1438 = conjq(V_ts)*e_em*conjq(V_Wp1)*U_su_25;
    const complex_t IT_1439 = IT_0025*IT_1438;
    const complex_t IT_1440 = m_u*conjq(V_u1)*V_us*e_em*IT_0035*U_su_35;
    const complex_t IT_1441 = IT_0034*IT_1440;
    const complex_t IT_1442 = 1.4142135623731*IT_1441;
    const complex_t IT_1443 = m_c*conjq(V_cs)*conjq(V_u1)*e_em*IT_0035*U_su_45;
    const complex_t IT_1444 = IT_0034*IT_1443;
    const complex_t IT_1445 = 1.4142135623731*IT_1444;
    const complex_t IT_1446 = m_t*conjq(V_ts)*conjq(V_u1)*e_em*IT_0035*U_su_55;
    const complex_t IT_1447 = IT_0034*IT_1446;
    const complex_t IT_1448 = 1.4142135623731*IT_1447;
    const complex_t IT_1449 = (complex_t{0, 1})*(IT_1435 + IT_1437 + IT_1439 +
       (-0.5)*IT_1442 + (-0.5)*IT_1445 + (-0.5)*IT_1448);
    const complex_t IT_1450 = IT_1179*IT_1382*IT_1449;
    const complex_t IT_1451 = IT_0023*IT_1450;
    const complex_t IT_1452 = mty::lt::C0iC(0, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0065, IT_0648, mty::lt::reg_int);
    const complex_t IT_1453 = mty::lt::C0iC(3, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0065, IT_0648, mty::lt::reg_int);
    const complex_t IT_1454 = mty::lt::C0iC(6, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0065, IT_0648, mty::lt::reg_int);
    const complex_t IT_1455 = IT_1452 + IT_1453 + IT_1454;
    const complex_t IT_1456 = IT_1451*IT_1455;
    const complex_t IT_1457 = IT_1039*IT_1104*IT_1382;
    const complex_t IT_1458 = 0.101321183642338*IT_1457;
    const complex_t IT_1459 = m_b*IT_1387;
    const complex_t IT_1460 = m_b*IT_1386;
    const complex_t IT_1461 = mty::lt::C0iC(12, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0065, IT_0250, mty::lt::reg_int);
    const complex_t IT_1462 = m_b*IT_1461;
    const complex_t IT_1463 = mty::lt::C0iC(18, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0065, IT_0250, mty::lt::reg_int);
    const complex_t IT_1464 = m_b*IT_1463;
    const complex_t IT_1465 = IT_1459 + IT_1460 + IT_1462 + IT_1464;
    const complex_t IT_1466 = m_s*IT_1386;
    const complex_t IT_1467 = m_s*IT_1461;
    const complex_t IT_1468 = mty::lt::C0iC(15, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0065, IT_0250, mty::lt::reg_int);
    const complex_t IT_1469 = m_s*IT_1468;
    const complex_t IT_1470 = m_b*IT_1468;
    const complex_t IT_1471 = -IT_1466 + -IT_1467 + -IT_1469 + 2*IT_1470;
    const complex_t IT_1472 = IT_1465 + IT_1471;
    const complex_t IT_1473 = IT_1458*IT_1472;
    const complex_t IT_1474 = IT_0119*IT_0383*IT_1382;
    const complex_t IT_1475 = 0.101321183642338*IT_1474;
    const complex_t IT_1476 = mty::lt::C0iC(18, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0065, IT_0133, mty::lt::reg_int);
    const complex_t IT_1477 = m_b*IT_1476;
    const complex_t IT_1478 = m_b*IT_1408;
    const complex_t IT_1479 = m_b*IT_1407;
    const complex_t IT_1480 = mty::lt::C0iC(12, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0065, IT_0133, mty::lt::reg_int);
    const complex_t IT_1481 = m_b*IT_1480;
    const complex_t IT_1482 = IT_1477 + IT_1478 + IT_1479 + IT_1481;
    const complex_t IT_1483 = m_s*IT_1407;
    const complex_t IT_1484 = m_s*IT_1480;
    const complex_t IT_1485 = mty::lt::C0iC(15, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0065, IT_0133, mty::lt::reg_int);
    const complex_t IT_1486 = m_s*IT_1485;
    const complex_t IT_1487 = m_b*IT_1485;
    const complex_t IT_1488 = -IT_1483 + -IT_1484 + -IT_1486 + 2*IT_1487;
    const complex_t IT_1489 = IT_1482 + IT_1488;
    const complex_t IT_1490 = IT_1475*IT_1489;
    const complex_t IT_1491 = IT_0922*IT_1014*IT_1382;
    const complex_t IT_1492 = 0.101321183642338*IT_1491;
    const complex_t IT_1493 = m_b*IT_1401;
    const complex_t IT_1494 = mty::lt::C0iC(12, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0065, IT_0140, mty::lt::reg_int);
    const complex_t IT_1495 = m_b*IT_1494;
    const complex_t IT_1496 = mty::lt::C0iC(18, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0065, IT_0140, mty::lt::reg_int);
    const complex_t IT_1497 = m_b*IT_1496;
    const complex_t IT_1498 = m_b*IT_1400;
    const complex_t IT_1499 = IT_1493 + IT_1495 + IT_1497 + IT_1498;
    const complex_t IT_1500 = mty::lt::C0iC(15, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0065, IT_0140, mty::lt::reg_int);
    const complex_t IT_1501 = m_b*IT_1500;
    const complex_t IT_1502 = m_s*IT_1494;
    const complex_t IT_1503 = m_s*IT_1500;
    const complex_t IT_1504 = m_s*IT_1400;
    const complex_t IT_1505 = 2*IT_1501 + -IT_1502 + -IT_1503 + -IT_1504;
    const complex_t IT_1506 = IT_1499 + IT_1505;
    const complex_t IT_1507 = IT_1492*IT_1506;
    const complex_t IT_1508 = IT_0271*IT_0345*IT_1382;
    const complex_t IT_1509 = 0.101321183642338*IT_1508;
    const complex_t IT_1510 = m_b*IT_1392;
    const complex_t IT_1511 = mty::lt::C0iC(18, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0065, IT_0016, mty::lt::reg_int);
    const complex_t IT_1512 = m_b*IT_1511;
    const complex_t IT_1513 = m_b*IT_1391;
    const complex_t IT_1514 = mty::lt::C0iC(12, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0065, IT_0016, mty::lt::reg_int);
    const complex_t IT_1515 = m_b*IT_1514;
    const complex_t IT_1516 = IT_1510 + IT_1512 + IT_1513 + IT_1515;
    const complex_t IT_1517 = mty::lt::C0iC(15, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0065, IT_0016, mty::lt::reg_int);
    const complex_t IT_1518 = m_s*IT_1517;
    const complex_t IT_1519 = m_b*IT_1517;
    const complex_t IT_1520 = m_s*IT_1391;
    const complex_t IT_1521 = m_s*IT_1514;
    const complex_t IT_1522 = -IT_1518 + 2*IT_1519 + -IT_1520 + -IT_1521;
    const complex_t IT_1523 = IT_1516 + IT_1522;
    const complex_t IT_1524 = IT_1509*IT_1523;
    const complex_t IT_1525 = IT_1286*IT_1382*IT_1430;
    const complex_t IT_1526 = 0.101321183642338*IT_1525;
    const complex_t IT_1527 = m_b*IT_1413;
    const complex_t IT_1528 = m_b*IT_1412;
    const complex_t IT_1529 = mty::lt::C0iC(12, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0065, IT_0606, mty::lt::reg_int);
    const complex_t IT_1530 = m_b*IT_1529;
    const complex_t IT_1531 = mty::lt::C0iC(18, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0065, IT_0606, mty::lt::reg_int);
    const complex_t IT_1532 = m_b*IT_1531;
    const complex_t IT_1533 = IT_1527 + IT_1528 + IT_1530 + IT_1532;
    const complex_t IT_1534 = m_s*IT_1412;
    const complex_t IT_1535 = m_s*IT_1529;
    const complex_t IT_1536 = mty::lt::C0iC(15, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0065, IT_0606, mty::lt::reg_int);
    const complex_t IT_1537 = m_s*IT_1536;
    const complex_t IT_1538 = m_b*IT_1536;
    const complex_t IT_1539 = -IT_1534 + -IT_1535 + -IT_1537 + 2*IT_1538;
    const complex_t IT_1540 = IT_1533 + IT_1539;
    const complex_t IT_1541 = IT_1526*IT_1540;
    const complex_t IT_1542 = IT_1319*IT_1382*IT_1449;
    const complex_t IT_1543 = 0.101321183642338*IT_1542;
    const complex_t IT_1544 = m_b*IT_1454;
    const complex_t IT_1545 = m_b*IT_1453;
    const complex_t IT_1546 = mty::lt::C0iC(12, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0065, IT_0648, mty::lt::reg_int);
    const complex_t IT_1547 = m_b*IT_1546;
    const complex_t IT_1548 = mty::lt::C0iC(18, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0065, IT_0648, mty::lt::reg_int);
    const complex_t IT_1549 = m_b*IT_1548;
    const complex_t IT_1550 = IT_1544 + IT_1545 + IT_1547 + IT_1549;
    const complex_t IT_1551 = m_s*IT_1453;
    const complex_t IT_1552 = m_s*IT_1546;
    const complex_t IT_1553 = mty::lt::C0iC(15, IT_0013 + IT_0014 + (-2)*s_12,
       IT_0014, IT_0013, IT_0015, IT_0065, IT_0648, mty::lt::reg_int);
    const complex_t IT_1554 = m_s*IT_1553;
    const complex_t IT_1555 = m_b*IT_1553;
    const complex_t IT_1556 = -IT_1551 + -IT_1552 + -IT_1554 + 2*IT_1555;
    const complex_t IT_1557 = IT_1550 + IT_1556;
    const complex_t IT_1558 = IT_1543*IT_1557;
    const complex_t IT_1559 = e_em*U_Wm2*conjq(U_Wm2);
    const complex_t IT_1560 = IT_0145*IT_1559;
    const complex_t IT_1561 = U_d2*conjq(U_d2)*e_em;
    const complex_t IT_1562 = IT_0145*IT_1561;
    const complex_t IT_1563 = (complex_t{0, 1})*(IT_1560 + IT_1562);
    const complex_t IT_1564 = -IT_1563;
    const complex_t IT_1565 = IT_0231*IT_0483*IT_1564;
    const complex_t IT_1566 = 0.101321183642338*IT_1565;
    const complex_t IT_1567 = IT_1218*IT_1566;
    const complex_t IT_1568 = IT_0061*IT_0323*IT_1564;
    const complex_t IT_1569 = 0.101321183642338*IT_1568;
    const complex_t IT_1570 = IT_1269*IT_1569;
    const complex_t IT_1571 = IT_0543*IT_0956*IT_1564;
    const complex_t IT_1572 = 0.101321183642338*IT_1571;
    const complex_t IT_1573 = IT_1252*IT_1572;
    const complex_t IT_1574 = IT_0130*IT_0412*IT_1564;
    const complex_t IT_1575 = 0.101321183642338*IT_1574;
    const complex_t IT_1576 = IT_1235*IT_1575;
    const complex_t IT_1577 = IT_0592*IT_1165*IT_1564;
    const complex_t IT_1578 = 0.101321183642338*IT_1577;
    const complex_t IT_1579 = IT_1302*IT_1578;
    const complex_t IT_1580 = IT_0645*IT_1179*IT_1564;
    const complex_t IT_1581 = 0.101321183642338*IT_1580;
    const complex_t IT_1582 = IT_1335*IT_1581;
    const complex_t IT_1583 = IT_0323*IT_0345*IT_1564;
    const complex_t IT_1584 = IT_0023*IT_1583;
    const complex_t IT_1585 = IT_1117*IT_1584;
    const complex_t IT_1586 = IT_0483*IT_1104*IT_1564;
    const complex_t IT_1587 = IT_0023*IT_1586;
    const complex_t IT_1588 = IT_0254*IT_1587;
    const complex_t IT_1589 = IT_0543*IT_1014*IT_1564;
    const complex_t IT_1590 = IT_0023*IT_1589;
    const complex_t IT_1591 = IT_1126*IT_1590;
    const complex_t IT_1592 = IT_0119*IT_0130*IT_1564;
    const complex_t IT_1593 = IT_0023*IT_1592;
    const complex_t IT_1594 = IT_1133*IT_1593;
    const complex_t IT_1595 = IT_0592*IT_1286*IT_1564;
    const complex_t IT_1596 = IT_0023*IT_1595;
    const complex_t IT_1597 = IT_1138*IT_1596;
    const complex_t IT_1598 = IT_0645*IT_1319*IT_1564;
    const complex_t IT_1599 = IT_0023*IT_1598;
    const complex_t IT_1600 = IT_1201*IT_1599;
    const complex_t IT_1601 = e_em*conjq(U_Wm1)*U_Wm2;
    const complex_t IT_1602 = IT_0145*IT_1601;
    const complex_t IT_1603 = conjq(U_d1)*U_d2*e_em;
    const complex_t IT_1604 = IT_0145*IT_1603;
    const complex_t IT_1605 = (complex_t{0, 1})*(IT_1602 + IT_1604);
    const complex_t IT_1606 = -IT_1605;
    const complex_t IT_1607 = IT_0231*IT_0754*IT_1606;
    const complex_t IT_1608 = 0.101321183642338*IT_1607;
    const complex_t IT_1609 = IT_1472*IT_1608;
    const complex_t IT_1610 = IT_0061*IT_0100*IT_1606;
    const complex_t IT_1611 = 0.101321183642338*IT_1610;
    const complex_t IT_1612 = IT_1523*IT_1611;
    const complex_t IT_1613 = IT_0178*IT_0956*IT_1606;
    const complex_t IT_1614 = 0.101321183642338*IT_1613;
    const complex_t IT_1615 = IT_1506*IT_1614;
    const complex_t IT_1616 = IT_0212*IT_0412*IT_1606;
    const complex_t IT_1617 = 0.101321183642338*IT_1616;
    const complex_t IT_1618 = IT_1489*IT_1617;
    const complex_t IT_1619 = IT_0837*IT_1165*IT_1606;
    const complex_t IT_1620 = 0.101321183642338*IT_1619;
    const complex_t IT_1621 = IT_1540*IT_1620;
    const complex_t IT_1622 = IT_0867*IT_1179*IT_1606;
    const complex_t IT_1623 = 0.101321183642338*IT_1622;
    const complex_t IT_1624 = IT_1557*IT_1623;
    const complex_t IT_1625 = IT_0100*IT_0345*IT_1606;
    const complex_t IT_1626 = IT_0023*IT_1625;
    const complex_t IT_1627 = IT_1393*IT_1626;
    const complex_t IT_1628 = IT_0754*IT_1104*IT_1606;
    const complex_t IT_1629 = IT_0023*IT_1628;
    const complex_t IT_1630 = IT_1388*IT_1629;
    const complex_t IT_1631 = IT_0119*IT_0212*IT_1606;
    const complex_t IT_1632 = IT_0023*IT_1631;
    const complex_t IT_1633 = IT_1409*IT_1632;
    const complex_t IT_1634 = IT_0178*IT_1014*IT_1606;
    const complex_t IT_1635 = IT_0023*IT_1634;
    const complex_t IT_1636 = IT_1402*IT_1635;
    const complex_t IT_1637 = IT_0837*IT_1286*IT_1606;
    const complex_t IT_1638 = IT_0023*IT_1637;
    const complex_t IT_1639 = IT_1414*IT_1638;
    const complex_t IT_1640 = IT_0867*IT_1319*IT_1606;
    const complex_t IT_1641 = IT_0023*IT_1640;
    const complex_t IT_1642 = IT_1455*IT_1641;
    const complex_t IT_1643 = e_em*V_Wp1*conjq(V_Wp1);
    const complex_t IT_1644 = IT_0145*IT_1643;
    const complex_t IT_1645 = V_u1*conjq(V_u1)*e_em;
    const complex_t IT_1646 = IT_0145*IT_1645;
    const complex_t IT_1647 = (complex_t{0, 1})*(IT_1644 + IT_1646);
    const complex_t IT_1648 = IT_0472*IT_1039*IT_1647;
    const complex_t IT_1649 = IT_0071*IT_1648;
    const complex_t IT_1650 = IT_0893*IT_1649;
    const complex_t IT_1651 = IT_0271*IT_0282*IT_1647;
    const complex_t IT_1652 = IT_0071*IT_1651;
    const complex_t IT_1653 = IT_0890*IT_1652;
    const complex_t IT_1654 = IT_0532*IT_0922*IT_1647;
    const complex_t IT_1655 = IT_0071*IT_1654;
    const complex_t IT_1656 = IT_0144*IT_1655;
    const complex_t IT_1657 = IT_0367*IT_0383*IT_1647;
    const complex_t IT_1658 = IT_0071*IT_1657;
    const complex_t IT_1659 = IT_0185*IT_1658;
    const complex_t IT_1660 = IT_0603*IT_1430*IT_1647;
    const complex_t IT_1661 = IT_0071*IT_1660;
    const complex_t IT_1662 = IT_0898*IT_1661;
    const complex_t IT_1663 = IT_0634*IT_1449*IT_1647;
    const complex_t IT_1664 = IT_0071*IT_1663;
    const complex_t IT_1665 = IT_0903*IT_1664;
    const complex_t IT_1666 = IT_0688*IT_1039*IT_1647;
    const complex_t IT_1667 = 0.101321183642338*IT_1666;
    const complex_t IT_1668 = IT_0772*IT_1667;
    const complex_t IT_1669 = IT_0201*IT_0383*IT_1647;
    const complex_t IT_1670 = 0.101321183642338*IT_1669;
    const complex_t IT_1671 = IT_0825*IT_1670;
    const complex_t IT_1672 = IT_0167*IT_0922*IT_1647;
    const complex_t IT_1673 = 0.101321183642338*IT_1672;
    const complex_t IT_1674 = IT_0808*IT_1673;
    const complex_t IT_1675 = IT_0089*IT_0271*IT_1647;
    const complex_t IT_1676 = 0.101321183642338*IT_1675;
    const complex_t IT_1677 = IT_0791*IT_1676;
    const complex_t IT_1678 = IT_0717*IT_1430*IT_1647;
    const complex_t IT_1679 = 0.101321183642338*IT_1678;
    const complex_t IT_1680 = IT_0855*IT_1679;
    const complex_t IT_1681 = IT_0738*IT_1449*IT_1647;
    const complex_t IT_1682 = 0.101321183642338*IT_1681;
    const complex_t IT_1683 = IT_0885*IT_1682;
    const complex_t IT_1684 = IT_0603*IT_1430;
    const complex_t IT_1685 = IT_0071*IT_1684;
    const complex_t IT_1686 = mty::lt::C0iC(0, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0065, IT_0606, IT_0606, mty::lt::reg_int);
    const complex_t IT_1687 = IT_0012*IT_1686;
    const complex_t IT_1688 = mty::lt::C0iC(3, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0065, IT_0606, IT_0606, mty::lt::reg_int);
    const complex_t IT_1689 = IT_0019*IT_1688;
    const complex_t IT_1690 = IT_1687 + IT_1689;
    const complex_t IT_1691 = IT_1685*IT_1690;
    const complex_t IT_1692 = IT_0603*IT_0837;
    const complex_t IT_1693 = 0.101321183642338*IT_1692;
    const complex_t IT_1694 = IT_0012*IT_1688;
    const complex_t IT_1695 = m_b*IT_1694;
    const complex_t IT_1696 = mty::lt::C0iC(12, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0065, IT_0606, IT_0606, mty::lt::reg_int);
    const complex_t IT_1697 = IT_0019*IT_1696;
    const complex_t IT_1698 = m_b*IT_1697;
    const complex_t IT_1699 = mty::lt::C0iC(6, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0065, IT_0606, IT_0606, mty::lt::reg_int);
    const complex_t IT_1700 = IT_0012*IT_1699;
    const complex_t IT_1701 = m_s*IT_1700;
    const complex_t IT_1702 = mty::lt::C0iC(15, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0065, IT_0606, IT_0606, mty::lt::reg_int);
    const complex_t IT_1703 = IT_0019*IT_1702;
    const complex_t IT_1704 = m_s*IT_1703;
    const complex_t IT_1705 = IT_1695 + IT_1698 + IT_1701 + IT_1704;
    const complex_t IT_1706 = IT_1693*IT_1705;
    const complex_t IT_1707 = IT_1154*IT_1165;
    const complex_t IT_1708 = IT_0023*IT_1707;
    const complex_t IT_1709 = mty::lt::C0iC(0, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0015, IT_0606, IT_0606, mty::lt::reg_int);
    const complex_t IT_1710 = IT_0012*IT_1709;
    const complex_t IT_1711 = mty::lt::C0iC(3, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0015, IT_0606, IT_0606, mty::lt::reg_int);
    const complex_t IT_1712 = IT_0019*IT_1711;
    const complex_t IT_1713 = IT_1710 + IT_1712;
    const complex_t IT_1714 = IT_1708*IT_1713;
    const complex_t IT_1715 = IT_0592*IT_1165;
    const complex_t IT_1716 = 0.101321183642338*IT_1715;
    const complex_t IT_1717 = mty::lt::C0iC(15, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0015, IT_0606, IT_0606, mty::lt::reg_int);
    const complex_t IT_1718 = IT_0019*IT_1717;
    const complex_t IT_1719 = m_s*IT_1718;
    const complex_t IT_1720 = IT_0012*IT_1711;
    const complex_t IT_1721 = m_b*IT_1720;
    const complex_t IT_1722 = mty::lt::C0iC(12, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0015, IT_0606, IT_0606, mty::lt::reg_int);
    const complex_t IT_1723 = IT_0019*IT_1722;
    const complex_t IT_1724 = m_b*IT_1723;
    const complex_t IT_1725 = mty::lt::C0iC(6, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0015, IT_0606, IT_0606, mty::lt::reg_int);
    const complex_t IT_1726 = IT_0012*IT_1725;
    const complex_t IT_1727 = m_s*IT_1726;
    const complex_t IT_1728 = IT_1719 + IT_1721 + IT_1724 + IT_1727;
    const complex_t IT_1729 = IT_1716*IT_1728;
    const complex_t IT_1730 = IT_0717*IT_1430;
    const complex_t IT_1731 = 0.101321183642338*IT_1730;
    const complex_t IT_1732 = IT_1705*IT_1731;
    const complex_t IT_1733 = IT_1154*IT_1286;
    const complex_t IT_1734 = 0.101321183642338*IT_1733;
    const complex_t IT_1735 = IT_1728*IT_1734;
    const complex_t IT_1736 = IT_0717*IT_0837;
    const complex_t IT_1737 = IT_0071*IT_1736;
    const complex_t IT_1738 = IT_1690*IT_1737;
    const complex_t IT_1739 = IT_0592*IT_1286;
    const complex_t IT_1740 = IT_0023*IT_1739;
    const complex_t IT_1741 = IT_1713*IT_1740;
    const complex_t IT_1742 = IT_0634*IT_1449;
    const complex_t IT_1743 = IT_0071*IT_1742;
    const complex_t IT_1744 = mty::lt::C0iC(0, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0065, IT_0648, IT_0648, mty::lt::reg_int);
    const complex_t IT_1745 = IT_0012*IT_1744;
    const complex_t IT_1746 = mty::lt::C0iC(3, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0065, IT_0648, IT_0648, mty::lt::reg_int);
    const complex_t IT_1747 = IT_0019*IT_1746;
    const complex_t IT_1748 = IT_1745 + IT_1747;
    const complex_t IT_1749 = IT_1743*IT_1748;
    const complex_t IT_1750 = IT_0634*IT_0867;
    const complex_t IT_1751 = 0.101321183642338*IT_1750;
    const complex_t IT_1752 = IT_0012*IT_1746;
    const complex_t IT_1753 = m_b*IT_1752;
    const complex_t IT_1754 = mty::lt::C0iC(12, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0065, IT_0648, IT_0648, mty::lt::reg_int);
    const complex_t IT_1755 = IT_0019*IT_1754;
    const complex_t IT_1756 = m_b*IT_1755;
    const complex_t IT_1757 = mty::lt::C0iC(6, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0065, IT_0648, IT_0648, mty::lt::reg_int);
    const complex_t IT_1758 = IT_0012*IT_1757;
    const complex_t IT_1759 = m_s*IT_1758;
    const complex_t IT_1760 = mty::lt::C0iC(15, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0065, IT_0648, IT_0648, mty::lt::reg_int);
    const complex_t IT_1761 = IT_0019*IT_1760;
    const complex_t IT_1762 = m_s*IT_1761;
    const complex_t IT_1763 = IT_1753 + IT_1756 + IT_1759 + IT_1762;
    const complex_t IT_1764 = IT_1751*IT_1763;
    const complex_t IT_1765 = IT_1179*IT_1195;
    const complex_t IT_1766 = IT_0023*IT_1765;
    const complex_t IT_1767 = mty::lt::C0iC(0, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0015, IT_0648, IT_0648, mty::lt::reg_int);
    const complex_t IT_1768 = IT_0012*IT_1767;
    const complex_t IT_1769 = mty::lt::C0iC(3, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0015, IT_0648, IT_0648, mty::lt::reg_int);
    const complex_t IT_1770 = IT_0019*IT_1769;
    const complex_t IT_1771 = IT_1768 + IT_1770;
    const complex_t IT_1772 = IT_1766*IT_1771;
    const complex_t IT_1773 = IT_0645*IT_1179;
    const complex_t IT_1774 = 0.101321183642338*IT_1773;
    const complex_t IT_1775 = IT_0012*IT_1769;
    const complex_t IT_1776 = m_b*IT_1775;
    const complex_t IT_1777 = mty::lt::C0iC(12, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0015, IT_0648, IT_0648, mty::lt::reg_int);
    const complex_t IT_1778 = IT_0019*IT_1777;
    const complex_t IT_1779 = m_b*IT_1778;
    const complex_t IT_1780 = mty::lt::C0iC(6, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0015, IT_0648, IT_0648, mty::lt::reg_int);
    const complex_t IT_1781 = IT_0012*IT_1780;
    const complex_t IT_1782 = m_s*IT_1781;
    const complex_t IT_1783 = mty::lt::C0iC(15, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0015, IT_0648, IT_0648, mty::lt::reg_int);
    const complex_t IT_1784 = IT_0019*IT_1783;
    const complex_t IT_1785 = m_s*IT_1784;
    const complex_t IT_1786 = IT_1776 + IT_1779 + IT_1782 + IT_1785;
    const complex_t IT_1787 = IT_1774*IT_1786;
    const complex_t IT_1788 = IT_0738*IT_1449;
    const complex_t IT_1789 = 0.101321183642338*IT_1788;
    const complex_t IT_1790 = IT_1763*IT_1789;
    const complex_t IT_1791 = IT_1195*IT_1319;
    const complex_t IT_1792 = 0.101321183642338*IT_1791;
    const complex_t IT_1793 = IT_1786*IT_1792;
    const complex_t IT_1794 = IT_0738*IT_0867;
    const complex_t IT_1795 = IT_0071*IT_1794;
    const complex_t IT_1796 = IT_1748*IT_1795;
    const complex_t IT_1797 = IT_0645*IT_1319;
    const complex_t IT_1798 = IT_0023*IT_1797;
    const complex_t IT_1799 = IT_1771*IT_1798;
    const complex_t IT_1800 = IT_0064 + IT_0103 + IT_0139 + 2*IT_0181 + 2
      *IT_0215 + (-2)*IT_0255 + IT_0285 + IT_0300 + IT_0326 + IT_0329 + IT_0348 
      + IT_0351 + IT_0386 + IT_0401 + IT_0431 + IT_0446 + IT_0449 + IT_0452 +
       IT_0455 + 2*IT_0502 + 2*IT_0521 + 2*IT_0562 + 2*IT_0581 + 2*IT_0623 + 2
      *IT_0665 + 2*IT_0670 + 2*IT_0691 + 2*IT_0696 + 2*IT_0701 + 2*IT_0722 + 2
      *IT_0743 + 2*IT_0773 + 2*IT_0792 + 2*IT_0809 + 2*IT_0826 + 2*IT_0856 + 2
      *IT_0886 + 2*IT_0891 + 2*IT_0896 + 2*IT_0901 + 2*IT_0906 + IT_0930 +
       IT_0945 + IT_0980 + IT_0995 + IT_0998 + IT_1017 + IT_1020 + IT_1023 +
       IT_1047 + IT_1062 + IT_1070 + IT_1085 + IT_1088 + IT_1107 + IT_1110 +
       IT_1113 + (-2)*IT_1120 + (-2)*IT_1127 + (-2)*IT_1134 + (-2)*IT_1168 + (-2
      )*IT_1202 + (-2)*IT_1219 + (-2)*IT_1236 + (-2)*IT_1253 + (-2)*IT_1270 + (
      -2)*IT_1303 + (-2)*IT_1336 + (-2)*IT_1344 + (-2)*IT_1347 + (-2)*IT_1350 + 
      (-2)*IT_1353 + (-2)*IT_1356 + (-2)*IT_1359 + (-2)*IT_1362 + (-2)*IT_1365 +
       (-2)*IT_1368 + (-2)*IT_1371 + (-2)*IT_1374 + (-2)*IT_1377 + (-2)*IT_1389 
      + (-2)*IT_1396 + (-2)*IT_1403 + (-2)*IT_1410 + (-2)*IT_1433 + (-2)*IT_1456
       + (-2)*IT_1473 + (-2)*IT_1490 + (-2)*IT_1507 + (-2)*IT_1524 + (-2)
      *IT_1541 + (-2)*IT_1558 + 2*IT_1567 + 2*IT_1570 + 2*IT_1573 + 2*IT_1576 +
       2*IT_1579 + 2*IT_1582 + 2*IT_1585 + 2*IT_1588 + 2*IT_1591 + 2*IT_1594 + 2
      *IT_1597 + 2*IT_1600 + 2*IT_1609 + 2*IT_1612 + 2*IT_1615 + 2*IT_1618 + 2
      *IT_1621 + 2*IT_1624 + 2*IT_1627 + 2*IT_1630 + 2*IT_1633 + 2*IT_1636 + 2
      *IT_1639 + 2*IT_1642 + (-2)*IT_1650 + (-2)*IT_1653 + (-2)*IT_1656 + (-2)
      *IT_1659 + (-2)*IT_1662 + (-2)*IT_1665 + (-2)*IT_1668 + (-2)*IT_1671 + (-2
      )*IT_1674 + (-2)*IT_1677 + (-2)*IT_1680 + (-2)*IT_1683 + IT_1691 + IT_1706
       + IT_1714 + IT_1729 + IT_1732 + IT_1735 + IT_1738 + IT_1741 + IT_1749 +
       IT_1764 + IT_1772 + IT_1787 + IT_1790 + IT_1793 + IT_1796 + IT_1799;
    const complex_t IT_1801 = cpowq(IT_0024, 2);
    const complex_t IT_1802 = IT_1800*IT_1801;
    const complex_t IT_1803 = IT_0005*IT_1802;
    const complex_t IT_1804 = -IT_1803;
    const complex_t IT_1805 = IT_0201*IT_0212*IT_1647;
    const complex_t IT_1806 = IT_0071*IT_1805;
    const complex_t IT_1807 = IT_0182*IT_1806;
    const complex_t IT_1808 = IT_0019*IT_0354;
    const complex_t IT_1809 = IT_0019*IT_0394;
    const complex_t IT_1810 = -IT_1808 + -IT_1809;
    const complex_t IT_1811 = IT_0355 + IT_1810;
    const complex_t IT_1812 = IT_0454*IT_1811;
    const complex_t IT_1813 = IT_0045*IT_0061*IT_1564;
    const complex_t IT_1814 = IT_0023*IT_1813;
    const complex_t IT_1815 = IT_1115*IT_1814;
    const complex_t IT_1816 = IT_0061*IT_0271*IT_1606;
    const complex_t IT_1817 = IT_0071*IT_1816;
    const complex_t IT_1818 = IT_1391*IT_1817;
    const complex_t IT_1819 = IT_0922*IT_0956*IT_1606;
    const complex_t IT_1820 = IT_0071*IT_1819;
    const complex_t IT_1821 = IT_1400*IT_1820;
    const complex_t IT_1822 = IT_0383*IT_0412*IT_1606;
    const complex_t IT_1823 = IT_0071*IT_1822;
    const complex_t IT_1824 = IT_1407*IT_1823;
    const complex_t IT_1825 = IT_0247*IT_0461*IT_0472;
    const complex_t IT_1826 = IT_0023*IT_1825;
    const complex_t IT_1827 = IT_0488*IT_1826;
    const complex_t IT_1828 = IT_0367*IT_0428*IT_0461;
    const complex_t IT_1829 = IT_0023*IT_1828;
    const complex_t IT_1830 = IT_0565*IT_1829;
    const complex_t IT_1831 = IT_0461*IT_0634*IT_1195;
    const complex_t IT_1832 = IT_0023*IT_1831;
    const complex_t IT_1833 = IT_0649*IT_1832;
    const complex_t IT_1834 = IT_0151*IT_0472*IT_1039;
    const complex_t IT_1835 = IT_0071*IT_1834;
    const complex_t IT_1836 = IT_0759*IT_1835;
    const complex_t IT_1837 = IT_0151*IT_0367*IT_0383;
    const complex_t IT_1838 = IT_0071*IT_1837;
    const complex_t IT_1839 = IT_0182*IT_1838;
    const complex_t IT_1840 = IT_0151*IT_0603*IT_1430;
    const complex_t IT_1841 = IT_0071*IT_1840;
    const complex_t IT_1842 = IT_0842*IT_1841;
    const complex_t IT_1843 = IT_0220*IT_0323*IT_0345;
    const complex_t IT_1844 = IT_0023*IT_1843;
    const complex_t IT_1845 = IT_1115*IT_1844;
    const complex_t IT_1846 = IT_0119*IT_0130*IT_0220;
    const complex_t IT_1847 = IT_0023*IT_1846;
    const complex_t IT_1848 = IT_1131*IT_1847;
    const complex_t IT_1849 = IT_0220*IT_0592*IT_1286;
    const complex_t IT_1850 = IT_0023*IT_1849;
    const complex_t IT_1851 = IT_1136*IT_1850;
    const complex_t IT_1852 = IT_0100*IT_0345*IT_1382;
    const complex_t IT_1853 = IT_0071*IT_1852;
    const complex_t IT_1854 = IT_1391*IT_1853;
    const complex_t IT_1855 = IT_0231*IT_0247*IT_1564;
    const complex_t IT_1856 = IT_0023*IT_1855;
    const complex_t IT_1857 = IT_0252*IT_1856;
    const complex_t IT_1858 = IT_0019*IT_0066;
    const complex_t IT_1859 = IT_0019*IT_0291;
    const complex_t IT_1860 = -IT_1858 + -IT_1859;
    const complex_t IT_1861 = IT_0067 + IT_1860;
    const complex_t IT_1862 = IT_0284*IT_1861;
    const complex_t IT_1863 = IT_0287 + IT_0293;
    const complex_t IT_1864 = m_b*IT_0069;
    const complex_t IT_1865 = m_s*IT_1859;
    const complex_t IT_1866 = m_b*IT_0295;
    const complex_t IT_1867 = mty::lt::C0iC(18, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0065, IT_0016, IT_0016, mty::lt::reg_int);
    const complex_t IT_1868 = IT_0019*IT_1867;
    const complex_t IT_1869 = m_s*IT_1868;
    const complex_t IT_1870 = -IT_1864 + -IT_1865 + -IT_1866 + -IT_1869;
    const complex_t IT_1871 = IT_1863 + IT_1870;
    const complex_t IT_1872 = IT_0299*IT_1871;
    const complex_t IT_1873 = IT_0017*IT_0019;
    const complex_t IT_1874 = IT_0019*IT_0306;
    const complex_t IT_1875 = -IT_1873 + -IT_1874;
    const complex_t IT_1876 = IT_0018 + IT_1875;
    const complex_t IT_1877 = IT_0063*IT_1876;
    const complex_t IT_1878 = IT_0302 + IT_0308;
    const complex_t IT_1879 = m_b*IT_0021;
    const complex_t IT_1880 = m_s*IT_1874;
    const complex_t IT_1881 = m_b*IT_0310;
    const complex_t IT_1882 = mty::lt::C0iC(18, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0015, IT_0016, IT_0016, mty::lt::reg_int);
    const complex_t IT_1883 = IT_0019*IT_1882;
    const complex_t IT_1884 = m_s*IT_1883;
    const complex_t IT_1885 = -IT_1879 + -IT_1880 + -IT_1881 + -IT_1884;
    const complex_t IT_1886 = IT_1878 + IT_1885;
    const complex_t IT_1887 = IT_0325*IT_1886;
    const complex_t IT_1888 = IT_0328*IT_1871;
    const complex_t IT_1889 = IT_0347*IT_1886;
    const complex_t IT_1890 = IT_0102*IT_1861;
    const complex_t IT_1891 = IT_0350*IT_1876;
    const complex_t IT_1892 = IT_0385*IT_1811;
    const complex_t IT_1893 = IT_0390 + IT_0396;
    const complex_t IT_1894 = m_b*IT_0353;
    const complex_t IT_1895 = m_b*IT_0398;
    const complex_t IT_1896 = m_s*IT_1809;
    const complex_t IT_1897 = mty::lt::C0iC(18, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0065, IT_0133, IT_0133, mty::lt::reg_int);
    const complex_t IT_1898 = IT_0019*IT_1897;
    const complex_t IT_1899 = m_s*IT_1898;
    const complex_t IT_1900 = -IT_1894 + -IT_1895 + -IT_1896 + -IT_1899;
    const complex_t IT_1901 = IT_1893 + IT_1900;
    const complex_t IT_1902 = IT_0388*IT_1901;
    const complex_t IT_1903 = IT_0019*IT_0442;
    const complex_t IT_1904 = IT_0019*IT_0134;
    const complex_t IT_1905 = -IT_1903 + -IT_1904;
    const complex_t IT_1906 = IT_0135 + IT_1905;
    const complex_t IT_1907 = IT_0430*IT_1906;
    const complex_t IT_1908 = IT_0441 + IT_0444;
    const complex_t IT_1909 = m_b*IT_0137;
    const complex_t IT_1910 = m_b*IT_0438;
    const complex_t IT_1911 = m_s*IT_1903;
    const complex_t IT_1912 = mty::lt::C0iC(18, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0015, IT_0133, IT_0133, mty::lt::reg_int);
    const complex_t IT_1913 = IT_0019*IT_1912;
    const complex_t IT_1914 = m_s*IT_1913;
    const complex_t IT_1915 = -IT_1909 + -IT_1910 + -IT_1911 + -IT_1914;
    const complex_t IT_1916 = IT_1908 + IT_1915;
    const complex_t IT_1917 = IT_0433*IT_1916;
    const complex_t IT_1918 = IT_0448*IT_1901;
    const complex_t IT_1919 = IT_0451*IT_1916;
    const complex_t IT_1920 = IT_0132*IT_1906;
    const complex_t IT_1921 = IT_0489 + IT_0491 + IT_0497;
    const complex_t IT_1922 = -IT_0495 + -IT_0498;
    const complex_t IT_1923 = IT_1921 + IT_1922;
    const complex_t IT_1924 = IT_0485*IT_1923;
    const complex_t IT_1925 = IT_0045*IT_0282*IT_0461;
    const complex_t IT_1926 = IT_0023*IT_1925;
    const complex_t IT_1927 = IT_0507*IT_1926;
    const complex_t IT_1928 = IT_0508 + IT_0510 + IT_0516;
    const complex_t IT_1929 = -IT_0514 + -IT_0517;
    const complex_t IT_1930 = IT_1928 + IT_1929;
    const complex_t IT_1931 = IT_0504*IT_1930;
    const complex_t IT_1932 = IT_0461*IT_0532*IT_0972;
    const complex_t IT_1933 = IT_0023*IT_1932;
    const complex_t IT_1934 = IT_0546*IT_1933;
    const complex_t IT_1935 = IT_0547 + IT_0549 + IT_0557;
    const complex_t IT_1936 = -IT_0555 + -IT_0558;
    const complex_t IT_1937 = IT_1935 + IT_1936;
    const complex_t IT_1938 = IT_0545*IT_1937;
    const complex_t IT_1939 = IT_0566 + IT_0568 + IT_0576;
    const complex_t IT_1940 = -IT_0574 + -IT_0577;
    const complex_t IT_1941 = IT_1939 + IT_1940;
    const complex_t IT_1942 = IT_0564*IT_1941;
    const complex_t IT_1943 = IT_0461*IT_0603*IT_1154;
    const complex_t IT_1944 = IT_0023*IT_1943;
    const complex_t IT_1945 = IT_0607*IT_1944;
    const complex_t IT_1946 = IT_0608 + IT_0610 + IT_0618;
    const complex_t IT_1947 = -IT_0616 + -IT_0619;
    const complex_t IT_1948 = IT_1946 + IT_1947;
    const complex_t IT_1949 = IT_0605*IT_1948;
    const complex_t IT_1950 = IT_0650 + IT_0652 + IT_0660;
    const complex_t IT_1951 = -IT_0658 + -IT_0661;
    const complex_t IT_1952 = IT_1950 + IT_1951;
    const complex_t IT_1953 = IT_0647*IT_1952;
    const complex_t IT_1954 = IT_0760 + IT_0762 + IT_0768;
    const complex_t IT_1955 = -IT_0766 + -IT_0770;
    const complex_t IT_1956 = IT_1954 + IT_1955;
    const complex_t IT_1957 = IT_0756*IT_1956;
    const complex_t IT_1958 = IT_0151*IT_0271*IT_0282;
    const complex_t IT_1959 = IT_0071*IT_1958;
    const complex_t IT_1960 = IT_0778*IT_1959;
    const complex_t IT_1961 = IT_0779 + IT_0781 + IT_0787;
    const complex_t IT_1962 = -IT_0785 + -IT_0788;
    const complex_t IT_1963 = IT_1961 + IT_1962;
    const complex_t IT_1964 = IT_0775*IT_1963;
    const complex_t IT_1965 = IT_0151*IT_0532*IT_0922;
    const complex_t IT_1966 = IT_0071*IT_1965;
    const complex_t IT_1967 = IT_0141*IT_1966;
    const complex_t IT_1968 = IT_0797 + IT_0800 + IT_0804;
    const complex_t IT_1969 = -IT_0802 + -IT_0805;
    const complex_t IT_1970 = IT_1968 + IT_1969;
    const complex_t IT_1971 = IT_0794*IT_1970;
    const complex_t IT_1972 = IT_0813 + IT_0815 + IT_0821;
    const complex_t IT_1973 = -IT_0819 + -IT_0822;
    const complex_t IT_1974 = IT_1972 + IT_1973;
    const complex_t IT_1975 = IT_0811*IT_1974;
    const complex_t IT_1976 = IT_0843 + IT_0845 + IT_0851;
    const complex_t IT_1977 = -IT_0849 + -IT_0852;
    const complex_t IT_1978 = IT_1976 + IT_1977;
    const complex_t IT_1979 = IT_0839*IT_1978;
    const complex_t IT_1980 = IT_0151*IT_0634*IT_1449;
    const complex_t IT_1981 = IT_0071*IT_1980;
    const complex_t IT_1982 = IT_0872*IT_1981;
    const complex_t IT_1983 = IT_0873 + IT_0875 + IT_0881;
    const complex_t IT_1984 = -IT_0879 + -IT_0882;
    const complex_t IT_1985 = IT_1983 + IT_1984;
    const complex_t IT_1986 = IT_0869*IT_1985;
    const complex_t IT_1987 = IT_0019*IT_0925;
    const complex_t IT_1988 = IT_0019*IT_0938;
    const complex_t IT_1989 = -IT_1987 + -IT_1988;
    const complex_t IT_1990 = IT_0926 + IT_1989;
    const complex_t IT_1991 = IT_0924*IT_1990;
    const complex_t IT_1992 = IT_0934 + IT_0940;
    const complex_t IT_1993 = m_b*IT_0928;
    const complex_t IT_1994 = m_b*IT_0942;
    const complex_t IT_1995 = m_s*IT_1988;
    const complex_t IT_1996 = mty::lt::C0iC(18, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0065, IT_0140, IT_0140, mty::lt::reg_int);
    const complex_t IT_1997 = IT_0019*IT_1996;
    const complex_t IT_1998 = m_s*IT_1997;
    const complex_t IT_1999 = -IT_1993 + -IT_1994 + -IT_1995 + -IT_1998;
    const complex_t IT_2000 = IT_1992 + IT_1999;
    const complex_t IT_2001 = IT_0932*IT_2000;
    const complex_t IT_2002 = IT_0019*IT_0975;
    const complex_t IT_2003 = IT_0019*IT_0988;
    const complex_t IT_2004 = -IT_2002 + -IT_2003;
    const complex_t IT_2005 = IT_0976 + IT_2004;
    const complex_t IT_2006 = IT_0974*IT_2005;
    const complex_t IT_2007 = IT_0984 + IT_0990;
    const complex_t IT_2008 = m_b*IT_0978;
    const complex_t IT_2009 = m_b*IT_0992;
    const complex_t IT_2010 = m_s*IT_2003;
    const complex_t IT_2011 = mty::lt::C0iC(18, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0015, IT_0140, IT_0140, mty::lt::reg_int);
    const complex_t IT_2012 = IT_0019*IT_2011;
    const complex_t IT_2013 = m_s*IT_2012;
    const complex_t IT_2014 = -IT_2008 + -IT_2009 + -IT_2010 + -IT_2013;
    const complex_t IT_2015 = IT_2007 + IT_2014;
    const complex_t IT_2016 = IT_0982*IT_2015;
    const complex_t IT_2017 = IT_0997*IT_2000;
    const complex_t IT_2018 = IT_1016*IT_2015;
    const complex_t IT_2019 = IT_1019*IT_1990;
    const complex_t IT_2020 = IT_1022*IT_2005;
    const complex_t IT_2021 = IT_0019*IT_1055;
    const complex_t IT_2022 = IT_0019*IT_1044;
    const complex_t IT_2023 = -IT_2021 + -IT_2022;
    const complex_t IT_2024 = IT_1045 + IT_2023;
    const complex_t IT_2025 = IT_1041*IT_2024;
    const complex_t IT_2026 = IT_1051 + IT_1057;
    const complex_t IT_2027 = m_b*IT_1043;
    const complex_t IT_2028 = m_b*IT_1059;
    const complex_t IT_2029 = m_s*IT_2021;
    const complex_t IT_2030 = mty::lt::C0iC(18, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0065, IT_0250, IT_0250, mty::lt::reg_int);
    const complex_t IT_2031 = IT_0019*IT_2030;
    const complex_t IT_2032 = m_s*IT_2031;
    const complex_t IT_2033 = -IT_2027 + -IT_2028 + -IT_2029 + -IT_2032;
    const complex_t IT_2034 = IT_2026 + IT_2033;
    const complex_t IT_2035 = IT_1049*IT_2034;
    const complex_t IT_2036 = IT_0019*IT_1065;
    const complex_t IT_2037 = IT_0019*IT_1075;
    const complex_t IT_2038 = -IT_2036 + -IT_2037;
    const complex_t IT_2039 = IT_1066 + IT_2038;
    const complex_t IT_2040 = IT_1064*IT_2039;
    const complex_t IT_2041 = IT_1074 + IT_1077;
    const complex_t IT_2042 = m_b*IT_1068;
    const complex_t IT_2043 = m_b*IT_1082;
    const complex_t IT_2044 = m_s*IT_2037;
    const complex_t IT_2045 = mty::lt::C0iC(18, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0015, IT_0250, IT_0250, mty::lt::reg_int);
    const complex_t IT_2046 = IT_0019*IT_2045;
    const complex_t IT_2047 = m_s*IT_2046;
    const complex_t IT_2048 = -IT_2042 + -IT_2043 + -IT_2044 + -IT_2047;
    const complex_t IT_2049 = IT_2041 + IT_2048;
    const complex_t IT_2050 = IT_1072*IT_2049;
    const complex_t IT_2051 = IT_1087*IT_2034;
    const complex_t IT_2052 = IT_1106*IT_2049;
    const complex_t IT_2053 = IT_1109*IT_2024;
    const complex_t IT_2054 = IT_1112*IT_2039;
    const complex_t IT_2055 = IT_1208 + IT_1210 + IT_1214;
    const complex_t IT_2056 = -IT_1215 + -IT_1216;
    const complex_t IT_2057 = IT_2055 + IT_2056;
    const complex_t IT_2058 = IT_1204*IT_2057;
    const complex_t IT_2059 = IT_1223 + IT_1227 + IT_1232;
    const complex_t IT_2060 = -IT_1229 + -IT_1233;
    const complex_t IT_2061 = IT_2059 + IT_2060;
    const complex_t IT_2062 = IT_1221*IT_2061;
    const complex_t IT_2063 = IT_1240 + IT_1244 + IT_1250;
    const complex_t IT_2064 = -IT_1246 + -IT_1247;
    const complex_t IT_2065 = IT_2063 + IT_2064;
    const complex_t IT_2066 = IT_1238*IT_2065;
    const complex_t IT_2067 = IT_1257 + IT_1261 + IT_1266;
    const complex_t IT_2068 = -IT_1263 + -IT_1267;
    const complex_t IT_2069 = IT_2067 + IT_2068;
    const complex_t IT_2070 = IT_1255*IT_2069;
    const complex_t IT_2071 = IT_1290 + IT_1292 + IT_1300;
    const complex_t IT_2072 = -IT_1296 + -IT_1297;
    const complex_t IT_2073 = IT_2071 + IT_2072;
    const complex_t IT_2074 = IT_1288*IT_2073;
    const complex_t IT_2075 = IT_1322 + IT_1327 + IT_1333;
    const complex_t IT_2076 = -IT_1329 + -IT_1330;
    const complex_t IT_2077 = IT_2075 + IT_2076;
    const complex_t IT_2078 = IT_1321*IT_2077;
    const complex_t IT_2079 = IT_0220*IT_0483*IT_1104;
    const complex_t IT_2080 = IT_0023*IT_2079;
    const complex_t IT_2081 = IT_0252*IT_2080;
    const complex_t IT_2082 = IT_0220*IT_0543*IT_1014;
    const complex_t IT_2083 = IT_0023*IT_2082;
    const complex_t IT_2084 = IT_1124*IT_2083;
    const complex_t IT_2085 = IT_0220*IT_0645*IT_1319;
    const complex_t IT_2086 = IT_0023*IT_2085;
    const complex_t IT_2087 = IT_1199*IT_2086;
    const complex_t IT_2088 = IT_1361*IT_1923;
    const complex_t IT_2089 = IT_1364*IT_1941;
    const complex_t IT_2090 = IT_1367*IT_1937;
    const complex_t IT_2091 = IT_1370*IT_1930;
    const complex_t IT_2092 = IT_1373*IT_1948;
    const complex_t IT_2093 = IT_1376*IT_1952;
    const complex_t IT_2094 = IT_0089*IT_0323*IT_1341;
    const complex_t IT_2095 = IT_0023*IT_2094;
    const complex_t IT_2096 = IT_0507*IT_2095;
    const complex_t IT_2097 = IT_0483*IT_0688*IT_1341;
    const complex_t IT_2098 = IT_0023*IT_2097;
    const complex_t IT_2099 = IT_0488*IT_2098;
    const complex_t IT_2100 = IT_0167*IT_0543*IT_1341;
    const complex_t IT_2101 = IT_0023*IT_2100;
    const complex_t IT_2102 = IT_0546*IT_2101;
    const complex_t IT_2103 = IT_0130*IT_0201*IT_1341;
    const complex_t IT_2104 = IT_0023*IT_2103;
    const complex_t IT_2105 = IT_0565*IT_2104;
    const complex_t IT_2106 = IT_0592*IT_0717*IT_1341;
    const complex_t IT_2107 = IT_0023*IT_2106;
    const complex_t IT_2108 = IT_0607*IT_2107;
    const complex_t IT_2109 = IT_0645*IT_0738*IT_1341;
    const complex_t IT_2110 = IT_0023*IT_2109;
    const complex_t IT_2111 = IT_0649*IT_2110;
    const complex_t IT_2112 = IT_1460 + IT_1462 + IT_1470;
    const complex_t IT_2113 = -IT_1466 + -IT_1467;
    const complex_t IT_2114 = IT_2112 + IT_2113;
    const complex_t IT_2115 = IT_1458*IT_2114;
    const complex_t IT_2116 = IT_1479 + IT_1481 + IT_1487;
    const complex_t IT_2117 = -IT_1483 + -IT_1484;
    const complex_t IT_2118 = IT_2116 + IT_2117;
    const complex_t IT_2119 = IT_1475*IT_2118;
    const complex_t IT_2120 = IT_1495 + IT_1498 + IT_1501;
    const complex_t IT_2121 = -IT_1502 + -IT_1504;
    const complex_t IT_2122 = IT_2120 + IT_2121;
    const complex_t IT_2123 = IT_1492*IT_2122;
    const complex_t IT_2124 = IT_1513 + IT_1515 + IT_1519;
    const complex_t IT_2125 = -IT_1520 + -IT_1521;
    const complex_t IT_2126 = IT_2124 + IT_2125;
    const complex_t IT_2127 = IT_1509*IT_2126;
    const complex_t IT_2128 = IT_1528 + IT_1530 + IT_1538;
    const complex_t IT_2129 = -IT_1534 + -IT_1535;
    const complex_t IT_2130 = IT_2128 + IT_2129;
    const complex_t IT_2131 = IT_1526*IT_2130;
    const complex_t IT_2132 = IT_1545 + IT_1547 + IT_1555;
    const complex_t IT_2133 = -IT_1551 + -IT_1552;
    const complex_t IT_2134 = IT_2132 + IT_2133;
    const complex_t IT_2135 = IT_1543*IT_2134;
    const complex_t IT_2136 = IT_0754*IT_1104*IT_1382;
    const complex_t IT_2137 = IT_0071*IT_2136;
    const complex_t IT_2138 = IT_1386*IT_2137;
    const complex_t IT_2139 = IT_0119*IT_0212*IT_1382;
    const complex_t IT_2140 = IT_0071*IT_2139;
    const complex_t IT_2141 = IT_1407*IT_2140;
    const complex_t IT_2142 = IT_0178*IT_1014*IT_1382;
    const complex_t IT_2143 = IT_0071*IT_2142;
    const complex_t IT_2144 = IT_1400*IT_2143;
    const complex_t IT_2145 = IT_0837*IT_1286*IT_1382;
    const complex_t IT_2146 = IT_0071*IT_2145;
    const complex_t IT_2147 = IT_1412*IT_2146;
    const complex_t IT_2148 = IT_0867*IT_1319*IT_1382;
    const complex_t IT_2149 = IT_0071*IT_2148;
    const complex_t IT_2150 = IT_1453*IT_2149;
    const complex_t IT_2151 = IT_1566*IT_2057;
    const complex_t IT_2152 = IT_1569*IT_2069;
    const complex_t IT_2153 = IT_0956*IT_0972*IT_1564;
    const complex_t IT_2154 = IT_0023*IT_2153;
    const complex_t IT_2155 = IT_1124*IT_2154;
    const complex_t IT_2156 = IT_1572*IT_2065;
    const complex_t IT_2157 = IT_0412*IT_0428*IT_1564;
    const complex_t IT_2158 = IT_0023*IT_2157;
    const complex_t IT_2159 = IT_1131*IT_2158;
    const complex_t IT_2160 = IT_1575*IT_2061;
    const complex_t IT_2161 = IT_1154*IT_1165*IT_1564;
    const complex_t IT_2162 = IT_0023*IT_2161;
    const complex_t IT_2163 = IT_1136*IT_2162;
    const complex_t IT_2164 = IT_1578*IT_2073;
    const complex_t IT_2165 = IT_1179*IT_1195*IT_1564;
    const complex_t IT_2166 = IT_0023*IT_2165;
    const complex_t IT_2167 = IT_1199*IT_2166;
    const complex_t IT_2168 = IT_1581*IT_2077;
    const complex_t IT_2169 = IT_0231*IT_1039*IT_1606;
    const complex_t IT_2170 = IT_0071*IT_2169;
    const complex_t IT_2171 = IT_1386*IT_2170;
    const complex_t IT_2172 = IT_1608*IT_2114;
    const complex_t IT_2173 = IT_1611*IT_2126;
    const complex_t IT_2174 = IT_1614*IT_2122;
    const complex_t IT_2175 = IT_1617*IT_2118;
    const complex_t IT_2176 = IT_1165*IT_1430*IT_1606;
    const complex_t IT_2177 = IT_0071*IT_2176;
    const complex_t IT_2178 = IT_1412*IT_2177;
    const complex_t IT_2179 = IT_1620*IT_2130;
    const complex_t IT_2180 = IT_1179*IT_1449*IT_1606;
    const complex_t IT_2181 = IT_0071*IT_2180;
    const complex_t IT_2182 = IT_1453*IT_2181;
    const complex_t IT_2183 = IT_1623*IT_2134;
    const complex_t IT_2184 = IT_1667*IT_1956;
    const complex_t IT_2185 = IT_1670*IT_1974;
    const complex_t IT_2186 = IT_1673*IT_1970;
    const complex_t IT_2187 = IT_1676*IT_1963;
    const complex_t IT_2188 = IT_1679*IT_1978;
    const complex_t IT_2189 = IT_1682*IT_1985;
    const complex_t IT_2190 = IT_0089*IT_0100*IT_1647;
    const complex_t IT_2191 = IT_0071*IT_2190;
    const complex_t IT_2192 = IT_0778*IT_2191;
    const complex_t IT_2193 = IT_0688*IT_0754*IT_1647;
    const complex_t IT_2194 = IT_0071*IT_2193;
    const complex_t IT_2195 = IT_0759*IT_2194;
    const complex_t IT_2196 = IT_0167*IT_0178*IT_1647;
    const complex_t IT_2197 = IT_0071*IT_2196;
    const complex_t IT_2198 = IT_0141*IT_2197;
    const complex_t IT_2199 = IT_0717*IT_0837*IT_1647;
    const complex_t IT_2200 = IT_0071*IT_2199;
    const complex_t IT_2201 = IT_0842*IT_2200;
    const complex_t IT_2202 = IT_0738*IT_0867*IT_1647;
    const complex_t IT_2203 = IT_0071*IT_2202;
    const complex_t IT_2204 = IT_0872*IT_2203;
    const complex_t IT_2205 = IT_0019*IT_1699;
    const complex_t IT_2206 = IT_0019*IT_1686;
    const complex_t IT_2207 = -IT_2205 + -IT_2206;
    const complex_t IT_2208 = IT_1687 + IT_2207;
    const complex_t IT_2209 = IT_1685*IT_2208;
    const complex_t IT_2210 = IT_1695 + IT_1701;
    const complex_t IT_2211 = mty::lt::C0iC(18, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0065, IT_0606, IT_0606, mty::lt::reg_int);
    const complex_t IT_2212 = IT_0019*IT_2211;
    const complex_t IT_2213 = m_s*IT_2212;
    const complex_t IT_2214 = m_b*IT_1689;
    const complex_t IT_2215 = m_b*IT_1703;
    const complex_t IT_2216 = m_s*IT_2205;
    const complex_t IT_2217 = -IT_2213 + -IT_2214 + -IT_2215 + -IT_2216;
    const complex_t IT_2218 = IT_2210 + IT_2217;
    const complex_t IT_2219 = IT_1693*IT_2218;
    const complex_t IT_2220 = IT_0019*IT_1725;
    const complex_t IT_2221 = IT_0019*IT_1709;
    const complex_t IT_2222 = -IT_2220 + -IT_2221;
    const complex_t IT_2223 = IT_1710 + IT_2222;
    const complex_t IT_2224 = IT_1708*IT_2223;
    const complex_t IT_2225 = IT_1721 + IT_1727;
    const complex_t IT_2226 = mty::lt::C0iC(18, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0015, IT_0606, IT_0606, mty::lt::reg_int);
    const complex_t IT_2227 = IT_0019*IT_2226;
    const complex_t IT_2228 = m_s*IT_2227;
    const complex_t IT_2229 = m_b*IT_1712;
    const complex_t IT_2230 = m_b*IT_1718;
    const complex_t IT_2231 = m_s*IT_2220;
    const complex_t IT_2232 = -IT_2228 + -IT_2229 + -IT_2230 + -IT_2231;
    const complex_t IT_2233 = IT_2225 + IT_2232;
    const complex_t IT_2234 = IT_1716*IT_2233;
    const complex_t IT_2235 = IT_1731*IT_2218;
    const complex_t IT_2236 = IT_1734*IT_2233;
    const complex_t IT_2237 = IT_1737*IT_2208;
    const complex_t IT_2238 = IT_1740*IT_2223;
    const complex_t IT_2239 = IT_0019*IT_1757;
    const complex_t IT_2240 = IT_0019*IT_1744;
    const complex_t IT_2241 = -IT_2239 + -IT_2240;
    const complex_t IT_2242 = IT_1745 + IT_2241;
    const complex_t IT_2243 = IT_1743*IT_2242;
    const complex_t IT_2244 = IT_1753 + IT_1759;
    const complex_t IT_2245 = m_b*IT_1747;
    const complex_t IT_2246 = m_b*IT_1761;
    const complex_t IT_2247 = m_s*IT_2239;
    const complex_t IT_2248 = mty::lt::C0iC(18, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0065, IT_0648, IT_0648, mty::lt::reg_int);
    const complex_t IT_2249 = IT_0019*IT_2248;
    const complex_t IT_2250 = m_s*IT_2249;
    const complex_t IT_2251 = -IT_2245 + -IT_2246 + -IT_2247 + -IT_2250;
    const complex_t IT_2252 = IT_2244 + IT_2251;
    const complex_t IT_2253 = IT_1751*IT_2252;
    const complex_t IT_2254 = IT_0019*IT_1780;
    const complex_t IT_2255 = IT_0019*IT_1767;
    const complex_t IT_2256 = -IT_2254 + -IT_2255;
    const complex_t IT_2257 = IT_1768 + IT_2256;
    const complex_t IT_2258 = IT_1766*IT_2257;
    const complex_t IT_2259 = IT_1776 + IT_1782;
    const complex_t IT_2260 = m_b*IT_1770;
    const complex_t IT_2261 = m_b*IT_1784;
    const complex_t IT_2262 = m_s*IT_2254;
    const complex_t IT_2263 = mty::lt::C0iC(18, IT_0013, IT_0013 + IT_0014 + (
      -2)*s_12, IT_0014, IT_0015, IT_0648, IT_0648, mty::lt::reg_int);
    const complex_t IT_2264 = IT_0019*IT_2263;
    const complex_t IT_2265 = m_s*IT_2264;
    const complex_t IT_2266 = -IT_2260 + -IT_2261 + -IT_2262 + -IT_2265;
    const complex_t IT_2267 = IT_2259 + IT_2266;
    const complex_t IT_2268 = IT_1774*IT_2267;
    const complex_t IT_2269 = IT_1789*IT_2252;
    const complex_t IT_2270 = IT_1792*IT_2267;
    const complex_t IT_2271 = IT_1795*IT_2242;
    const complex_t IT_2272 = IT_1798*IT_2257;
    const complex_t IT_2273 = IT_1807 + (-0.5)*IT_1812 + -IT_1815 + -IT_1818 +
       -IT_1821 + -IT_1824 + -IT_1827 + -IT_1830 + -IT_1833 + -IT_1836 + 
      -IT_1839 + -IT_1842 + IT_1845 + IT_1848 + IT_1851 + IT_1854 + -IT_1857 + (
      -0.5)*IT_1862 + (-0.5)*IT_1872 + (-0.5)*IT_1877 + (-0.5)*IT_1887 + (-0.5)
      *IT_1888 + (-0.5)*IT_1889 + (-0.5)*IT_1890 + (-0.5)*IT_1891 + (-0.5)
      *IT_1892 + (-0.5)*IT_1902 + (-0.5)*IT_1907 + (-0.5)*IT_1917 + (-0.5)
      *IT_1918 + (-0.5)*IT_1919 + (-0.5)*IT_1920 + -IT_1924 + -IT_1927 + 
      -IT_1931 + -IT_1934 + -IT_1938 + -IT_1942 + -IT_1945 + -IT_1949 + -IT_1953
       + -IT_1957 + -IT_1960 + -IT_1964 + -IT_1967 + -IT_1971 + -IT_1975 + 
      -IT_1979 + -IT_1982 + -IT_1986 + (-0.5)*IT_1991 + (-0.5)*IT_2001 + (-0.5)
      *IT_2006 + (-0.5)*IT_2016 + (-0.5)*IT_2017 + (-0.5)*IT_2018 + (-0.5)
      *IT_2019 + (-0.5)*IT_2020 + (-0.5)*IT_2025 + (-0.5)*IT_2035 + (-0.5)
      *IT_2040 + (-0.5)*IT_2050 + (-0.5)*IT_2051 + (-0.5)*IT_2052 + (-0.5)
      *IT_2053 + (-0.5)*IT_2054 + IT_2058 + IT_2062 + IT_2066 + IT_2070 +
       IT_2074 + IT_2078 + IT_2081 + IT_2084 + IT_2087 + IT_2088 + IT_2089 +
       IT_2090 + IT_2091 + IT_2092 + IT_2093 + IT_2096 + IT_2099 + IT_2102 +
       IT_2105 + IT_2108 + IT_2111 + IT_2115 + IT_2119 + IT_2123 + IT_2127 +
       IT_2131 + IT_2135 + IT_2138 + IT_2141 + IT_2144 + IT_2147 + IT_2150 + 
      -IT_2151 + -IT_2152 + -IT_2155 + -IT_2156 + -IT_2159 + -IT_2160 + -IT_2163
       + -IT_2164 + -IT_2167 + -IT_2168 + -IT_2171 + -IT_2172 + -IT_2173 + 
      -IT_2174 + -IT_2175 + -IT_2178 + -IT_2179 + -IT_2182 + -IT_2183 + IT_2184 
      + IT_2185 + IT_2186 + IT_2187 + IT_2188 + IT_2189 + IT_2192 + IT_2195 +
       IT_2198 + IT_2201 + IT_2204 + (-0.5)*IT_2209 + (-0.5)*IT_2219 + (-0.5)
      *IT_2224 + (-0.5)*IT_2234 + (-0.5)*IT_2235 + (-0.5)*IT_2236 + (-0.5)
      *IT_2237 + (-0.5)*IT_2238 + (-0.5)*IT_2243 + (-0.5)*IT_2253 + (-0.5)
      *IT_2258 + (-0.5)*IT_2268 + (-0.5)*IT_2269 + (-0.5)*IT_2270 + (-0.5)
      *IT_2271 + (-0.5)*IT_2272;
    const complex_t IT_2274 = IT_1801*IT_2273;
    const complex_t IT_2275 = -IT_2274;
    const complex_t IT_2276 = 2*IT_2275;
    const complex_t IT_2277 = IT_0005*IT_2276;
    const complex_t IT_2278 = -IT_0293 + -IT_0296;
    const complex_t IT_2279 = IT_0287 + IT_0290;
    const complex_t IT_2280 = IT_2278 + IT_2279;
    const complex_t IT_2281 = IT_0299*IT_2280;
    const complex_t IT_2282 = IT_0302 + IT_0305;
    const complex_t IT_2283 = -IT_0308 + -IT_0311;
    const complex_t IT_2284 = IT_2282 + IT_2283;
    const complex_t IT_2285 = IT_0325*IT_2284;
    const complex_t IT_2286 = IT_0328*IT_2280;
    const complex_t IT_2287 = IT_0347*IT_2284;
    const complex_t IT_2288 = IT_0390 + IT_0393;
    const complex_t IT_2289 = -IT_0396 + -IT_0399;
    const complex_t IT_2290 = IT_2288 + IT_2289;
    const complex_t IT_2291 = IT_0388*IT_2290;
    const complex_t IT_2292 = IT_0436 + IT_0441;
    const complex_t IT_2293 = -IT_0439 + -IT_0444;
    const complex_t IT_2294 = IT_2292 + IT_2293;
    const complex_t IT_2295 = IT_0433*IT_2294;
    const complex_t IT_2296 = IT_0448*IT_2290;
    const complex_t IT_2297 = IT_0451*IT_2294;
    const complex_t IT_2298 = IT_0487 + IT_0489 + IT_0491 + IT_0493 + IT_0495 
      + IT_0498 + IT_0499;
    const complex_t IT_2299 = 2*IT_0497;
    const complex_t IT_2300 = IT_2298 + IT_2299;
    const complex_t IT_2301 = IT_0485*IT_2300;
    const complex_t IT_2302 = IT_0506 + IT_0508 + IT_0510 + IT_0512 + IT_0514 
      + IT_0517 + IT_0518;
    const complex_t IT_2303 = 2*IT_0516;
    const complex_t IT_2304 = IT_2302 + IT_2303;
    const complex_t IT_2305 = IT_0504*IT_2304;
    const complex_t IT_2306 = IT_0547 + IT_0549 + IT_0551 + IT_0553 + IT_0555 
      + IT_0558 + IT_0559;
    const complex_t IT_2307 = 2*IT_0557;
    const complex_t IT_2308 = IT_2306 + IT_2307;
    const complex_t IT_2309 = IT_0545*IT_2308;
    const complex_t IT_2310 = IT_0566 + IT_0568 + IT_0570 + IT_0572 + IT_0574 
      + IT_0577 + IT_0578;
    const complex_t IT_2311 = 2*IT_0576;
    const complex_t IT_2312 = IT_2310 + IT_2311;
    const complex_t IT_2313 = IT_0564*IT_2312;
    const complex_t IT_2314 = IT_0608 + IT_0610 + IT_0612 + IT_0614 + IT_0616 
      + IT_0619 + IT_0620;
    const complex_t IT_2315 = 2*IT_0618;
    const complex_t IT_2316 = IT_2314 + IT_2315;
    const complex_t IT_2317 = IT_0605*IT_2316;
    const complex_t IT_2318 = IT_0650 + IT_0652 + IT_0654 + IT_0656 + IT_0658 
      + IT_0661 + IT_0662;
    const complex_t IT_2319 = 2*IT_0660;
    const complex_t IT_2320 = IT_2318 + IT_2319;
    const complex_t IT_2321 = IT_0647*IT_2320;
    const complex_t IT_2322 = 2*IT_0768;
    const complex_t IT_2323 = IT_0758 + IT_0760 + IT_0762 + IT_0764 + IT_0766 
      + IT_0769 + IT_0770;
    const complex_t IT_2324 = IT_2322 + IT_2323;
    const complex_t IT_2325 = IT_0756*IT_2324;
    const complex_t IT_2326 = 2*IT_0787;
    const complex_t IT_2327 = IT_0777 + IT_0779 + IT_0781 + IT_0783 + IT_0785 
      + IT_0788 + IT_0789;
    const complex_t IT_2328 = IT_2326 + IT_2327;
    const complex_t IT_2329 = IT_0775*IT_2328;
    const complex_t IT_2330 = IT_0795 + IT_0797 + IT_0799 + IT_0800 + IT_0802 
      + IT_0805 + IT_0806;
    const complex_t IT_2331 = 2*IT_0804;
    const complex_t IT_2332 = IT_2330 + IT_2331;
    const complex_t IT_2333 = IT_0794*IT_2332;
    const complex_t IT_2334 = 2*IT_0821;
    const complex_t IT_2335 = IT_0812 + IT_0813 + IT_0815 + IT_0817 + IT_0819 
      + IT_0822 + IT_0823;
    const complex_t IT_2336 = IT_2334 + IT_2335;
    const complex_t IT_2337 = IT_0811*IT_2336;
    const complex_t IT_2338 = IT_0841 + IT_0843 + IT_0845 + IT_0847 + IT_0849 
      + IT_0852 + IT_0853;
    const complex_t IT_2339 = 2*IT_0851;
    const complex_t IT_2340 = IT_2338 + IT_2339;
    const complex_t IT_2341 = IT_0839*IT_2340;
    const complex_t IT_2342 = IT_0871 + IT_0873 + IT_0875 + IT_0877 + IT_0879 
      + IT_0882 + IT_0883;
    const complex_t IT_2343 = 2*IT_0881;
    const complex_t IT_2344 = IT_2342 + IT_2343;
    const complex_t IT_2345 = IT_0869*IT_2344;
    const complex_t IT_2346 = IT_0934 + IT_0937;
    const complex_t IT_2347 = -IT_0940 + -IT_0943;
    const complex_t IT_2348 = IT_2346 + IT_2347;
    const complex_t IT_2349 = IT_0932*IT_2348;
    const complex_t IT_2350 = IT_0984 + IT_0987;
    const complex_t IT_2351 = -IT_0990 + -IT_0993;
    const complex_t IT_2352 = IT_2350 + IT_2351;
    const complex_t IT_2353 = IT_0982*IT_2352;
    const complex_t IT_2354 = IT_0997*IT_2348;
    const complex_t IT_2355 = IT_1016*IT_2352;
    const complex_t IT_2356 = IT_1051 + IT_1054;
    const complex_t IT_2357 = -IT_1057 + -IT_1060;
    const complex_t IT_2358 = IT_2356 + IT_2357;
    const complex_t IT_2359 = IT_1049*IT_2358;
    const complex_t IT_2360 = IT_1074 + IT_1080;
    const complex_t IT_2361 = -IT_1077 + -IT_1083;
    const complex_t IT_2362 = IT_2360 + IT_2361;
    const complex_t IT_2363 = IT_1072*IT_2362;
    const complex_t IT_2364 = IT_1087*IT_2358;
    const complex_t IT_2365 = IT_1106*IT_2362;
    const complex_t IT_2366 = IT_1206 + IT_1207 + IT_1208 + IT_1210 + IT_1213 
      + IT_1215 + IT_1216;
    const complex_t IT_2367 = 2*IT_1214;
    const complex_t IT_2368 = IT_2366 + IT_2367;
    const complex_t IT_2369 = IT_1204*IT_2368;
    const complex_t IT_2370 = IT_1222 + IT_1223 + IT_1225 + IT_1227 + IT_1229 
      + IT_1231 + IT_1233;
    const complex_t IT_2371 = 2*IT_1232;
    const complex_t IT_2372 = IT_2370 + IT_2371;
    const complex_t IT_2373 = IT_1221*IT_2372;
    const complex_t IT_2374 = IT_1240 + IT_1242 + IT_1243 + IT_1244 + IT_1246 
      + IT_1247 + IT_1249;
    const complex_t IT_2375 = 2*IT_1250;
    const complex_t IT_2376 = IT_2374 + IT_2375;
    const complex_t IT_2377 = IT_1238*IT_2376;
    const complex_t IT_2378 = 2*IT_1266;
    const complex_t IT_2379 = IT_1257 + IT_1259 + IT_1260 + IT_1261 + IT_1263 
      + IT_1265 + IT_1267;
    const complex_t IT_2380 = IT_2378 + IT_2379;
    const complex_t IT_2381 = IT_1255*IT_2380;
    const complex_t IT_2382 = IT_1289 + IT_1290 + IT_1292 + IT_1294 + IT_1296 
      + IT_1297 + IT_1299;
    const complex_t IT_2383 = 2*IT_1300;
    const complex_t IT_2384 = IT_2382 + IT_2383;
    const complex_t IT_2385 = IT_1288*IT_2384;
    const complex_t IT_2386 = IT_1322 + IT_1324 + IT_1325 + IT_1327 + IT_1329 
      + IT_1330 + IT_1332;
    const complex_t IT_2387 = 2*IT_1333;
    const complex_t IT_2388 = IT_2386 + IT_2387;
    const complex_t IT_2389 = IT_1321*IT_2388;
    const complex_t IT_2390 = IT_1361*IT_2300;
    const complex_t IT_2391 = IT_1364*IT_2312;
    const complex_t IT_2392 = IT_1367*IT_2308;
    const complex_t IT_2393 = IT_1370*IT_2304;
    const complex_t IT_2394 = IT_1373*IT_2316;
    const complex_t IT_2395 = IT_1376*IT_2320;
    const complex_t IT_2396 = IT_1459 + IT_1460 + IT_1462 + IT_1464 + IT_1466 
      + IT_1467 + IT_1469;
    const complex_t IT_2397 = 2*IT_1470;
    const complex_t IT_2398 = IT_2396 + IT_2397;
    const complex_t IT_2399 = IT_1458*IT_2398;
    const complex_t IT_2400 = IT_1477 + IT_1478 + IT_1479 + IT_1481 + IT_1483 
      + IT_1484 + IT_1486;
    const complex_t IT_2401 = 2*IT_1487;
    const complex_t IT_2402 = IT_2400 + IT_2401;
    const complex_t IT_2403 = IT_1475*IT_2402;
    const complex_t IT_2404 = IT_1493 + IT_1495 + IT_1497 + IT_1498 + IT_1502 
      + IT_1503 + IT_1504;
    const complex_t IT_2405 = 2*IT_1501;
    const complex_t IT_2406 = IT_2404 + IT_2405;
    const complex_t IT_2407 = IT_1492*IT_2406;
    const complex_t IT_2408 = IT_1510 + IT_1512 + IT_1513 + IT_1515 + IT_1518 
      + IT_1520 + IT_1521;
    const complex_t IT_2409 = 2*IT_1519;
    const complex_t IT_2410 = IT_2408 + IT_2409;
    const complex_t IT_2411 = IT_1509*IT_2410;
    const complex_t IT_2412 = 2*IT_1538;
    const complex_t IT_2413 = IT_1527 + IT_1528 + IT_1530 + IT_1532 + IT_1534 
      + IT_1535 + IT_1537;
    const complex_t IT_2414 = IT_2412 + IT_2413;
    const complex_t IT_2415 = IT_1526*IT_2414;
    const complex_t IT_2416 = IT_1544 + IT_1545 + IT_1547 + IT_1549 + IT_1551 
      + IT_1552 + IT_1554;
    const complex_t IT_2417 = 2*IT_1555;
    const complex_t IT_2418 = IT_2416 + IT_2417;
    const complex_t IT_2419 = IT_1543*IT_2418;
    const complex_t IT_2420 = IT_1566*IT_2368;
    const complex_t IT_2421 = IT_1569*IT_2380;
    const complex_t IT_2422 = IT_1572*IT_2376;
    const complex_t IT_2423 = IT_1575*IT_2372;
    const complex_t IT_2424 = IT_1578*IT_2384;
    const complex_t IT_2425 = IT_1581*IT_2388;
    const complex_t IT_2426 = IT_1608*IT_2398;
    const complex_t IT_2427 = IT_1611*IT_2410;
    const complex_t IT_2428 = IT_1614*IT_2406;
    const complex_t IT_2429 = IT_1617*IT_2402;
    const complex_t IT_2430 = IT_1620*IT_2414;
    const complex_t IT_2431 = IT_1623*IT_2418;
    const complex_t IT_2432 = IT_1667*IT_2324;
    const complex_t IT_2433 = IT_1670*IT_2336;
    const complex_t IT_2434 = IT_1673*IT_2332;
    const complex_t IT_2435 = IT_1676*IT_2328;
    const complex_t IT_2436 = IT_1679*IT_2340;
    const complex_t IT_2437 = IT_1682*IT_2344;
    const complex_t IT_2438 = IT_1695 + IT_1698;
    const complex_t IT_2439 = -IT_1701 + -IT_1704;
    const complex_t IT_2440 = IT_2438 + IT_2439;
    const complex_t IT_2441 = IT_1693*IT_2440;
    const complex_t IT_2442 = IT_1721 + IT_1724;
    const complex_t IT_2443 = -IT_1719 + -IT_1727;
    const complex_t IT_2444 = IT_2442 + IT_2443;
    const complex_t IT_2445 = IT_1716*IT_2444;
    const complex_t IT_2446 = IT_1731*IT_2440;
    const complex_t IT_2447 = IT_1734*IT_2444;
    const complex_t IT_2448 = IT_1753 + IT_1756;
    const complex_t IT_2449 = -IT_1759 + -IT_1762;
    const complex_t IT_2450 = IT_2448 + IT_2449;
    const complex_t IT_2451 = IT_1751*IT_2450;
    const complex_t IT_2452 = IT_1776 + IT_1779;
    const complex_t IT_2453 = -IT_1782 + -IT_1785;
    const complex_t IT_2454 = IT_2452 + IT_2453;
    const complex_t IT_2455 = IT_1774*IT_2454;
    const complex_t IT_2456 = IT_1789*IT_2450;
    const complex_t IT_2457 = IT_1792*IT_2454;
    const complex_t IT_2458 = IT_0064 + -IT_0103 + -IT_0139 + (-2)*IT_0181 + (
      -2)*IT_0215 + (-2)*IT_0255 + IT_0285 + -IT_0351 + IT_0386 + IT_0431 + 
      -IT_0455 + (-2)*IT_0670 + (-2)*IT_0691 + (-2)*IT_0696 + (-2)*IT_0701 + (-2
      )*IT_0722 + (-2)*IT_0743 + (-2)*IT_0891 + (-2)*IT_0896 + (-2)*IT_0901 + (
      -2)*IT_0906 + IT_0930 + IT_0980 + -IT_1020 + -IT_1023 + IT_1047 + IT_1070 
      + -IT_1110 + -IT_1113 + (-2)*IT_1120 + (-2)*IT_1127 + (-2)*IT_1134 + (-2)
      *IT_1168 + (-2)*IT_1202 + (-2)*IT_1344 + (-2)*IT_1347 + (-2)*IT_1350 + (-2
      )*IT_1353 + (-2)*IT_1356 + (-2)*IT_1359 + (-2)*IT_1389 + (-2)*IT_1396 + (
      -2)*IT_1403 + (-2)*IT_1410 + (-2)*IT_1433 + (-2)*IT_1456 + (-2)*IT_1585 + 
      (-2)*IT_1588 + (-2)*IT_1591 + (-2)*IT_1594 + (-2)*IT_1597 + (-2)*IT_1600 +
       (-2)*IT_1627 + (-2)*IT_1630 + (-2)*IT_1633 + (-2)*IT_1636 + (-2)*IT_1639 
      + (-2)*IT_1642 + (-2)*IT_1650 + (-2)*IT_1653 + (-2)*IT_1656 + (-2)*IT_1659
       + (-2)*IT_1662 + (-2)*IT_1665 + IT_1691 + IT_1714 + -IT_1738 + -IT_1741 +
       IT_1749 + IT_1772 + -IT_1796 + -IT_1799 + -IT_2281 + -IT_2285 + IT_2286 +
       IT_2287 + -IT_2291 + -IT_2295 + IT_2296 + IT_2297 + (-2)*IT_2301 + (-2)
      *IT_2305 + (-2)*IT_2309 + (-2)*IT_2313 + (-2)*IT_2317 + (-2)*IT_2321 + (-2
      )*IT_2325 + (-2)*IT_2329 + (-2)*IT_2333 + (-2)*IT_2337 + (-2)*IT_2341 + (
      -2)*IT_2345 + -IT_2349 + -IT_2353 + IT_2354 + IT_2355 + -IT_2359 + 
      -IT_2363 + IT_2364 + IT_2365 + (-2)*IT_2369 + (-2)*IT_2373 + (-2)*IT_2377 
      + (-2)*IT_2381 + (-2)*IT_2385 + (-2)*IT_2389 + (-2)*IT_2390 + (-2)*IT_2391
       + (-2)*IT_2392 + (-2)*IT_2393 + (-2)*IT_2394 + (-2)*IT_2395 + (-2)
      *IT_2399 + (-2)*IT_2403 + (-2)*IT_2407 + (-2)*IT_2411 + (-2)*IT_2415 + (-2
      )*IT_2419 + (-2)*IT_2420 + (-2)*IT_2421 + (-2)*IT_2422 + (-2)*IT_2423 + (
      -2)*IT_2424 + (-2)*IT_2425 + (-2)*IT_2426 + (-2)*IT_2427 + (-2)*IT_2428 + 
      (-2)*IT_2429 + (-2)*IT_2430 + (-2)*IT_2431 + (-2)*IT_2432 + (-2)*IT_2433 +
       (-2)*IT_2434 + (-2)*IT_2435 + (-2)*IT_2436 + (-2)*IT_2437 + -IT_2441 + 
      -IT_2445 + IT_2446 + IT_2447 + -IT_2451 + -IT_2455 + IT_2456 + IT_2457;
    const complex_t IT_2459 = IT_1801*IT_2458;
    const complex_t IT_2460 = (-0.5)*IT_2459;
    const complex_t IT_2461 = IT_0005*IT_2460;
    const complex_t IT_2462 = 2*IT_2461;
    const complex_t IT_2463 = IT_0287 + IT_1865 + IT_1869;
    const complex_t IT_2464 = -IT_0293 + -IT_1864 + -IT_1866;
    const complex_t IT_2465 = IT_2463 + IT_2464;
    const complex_t IT_2466 = IT_0299*IT_2465;
    const complex_t IT_2467 = IT_0302 + IT_1880 + IT_1884;
    const complex_t IT_2468 = -IT_0308 + -IT_1879 + -IT_1881;
    const complex_t IT_2469 = IT_2467 + IT_2468;
    const complex_t IT_2470 = IT_0325*IT_2469;
    const complex_t IT_2471 = IT_0328*IT_2465;
    const complex_t IT_2472 = IT_0347*IT_2469;
    const complex_t IT_2473 = IT_0390 + IT_1896 + IT_1899;
    const complex_t IT_2474 = -IT_0396 + -IT_1894 + -IT_1895;
    const complex_t IT_2475 = IT_2473 + IT_2474;
    const complex_t IT_2476 = IT_0388*IT_2475;
    const complex_t IT_2477 = IT_0441 + IT_1911 + IT_1914;
    const complex_t IT_2478 = -IT_0444 + -IT_1909 + -IT_1910;
    const complex_t IT_2479 = IT_2477 + IT_2478;
    const complex_t IT_2480 = IT_0433*IT_2479;
    const complex_t IT_2481 = IT_0448*IT_2475;
    const complex_t IT_2482 = IT_0451*IT_2479;
    const complex_t IT_2483 = IT_0489 + IT_0491 + IT_0495 + IT_0497 + IT_0498;
    const complex_t IT_2484 = IT_0485*IT_2483;
    const complex_t IT_2485 = IT_0508 + IT_0510 + IT_0514 + IT_0516 + IT_0517;
    const complex_t IT_2486 = IT_0504*IT_2485;
    const complex_t IT_2487 = IT_0547 + IT_0549 + IT_0555 + IT_0557 + IT_0558;
    const complex_t IT_2488 = IT_0545*IT_2487;
    const complex_t IT_2489 = IT_0566 + IT_0568 + IT_0574 + IT_0576 + IT_0577;
    const complex_t IT_2490 = IT_0564*IT_2489;
    const complex_t IT_2491 = IT_0608 + IT_0610 + IT_0616 + IT_0618 + IT_0619;
    const complex_t IT_2492 = IT_0605*IT_2491;
    const complex_t IT_2493 = IT_0650 + IT_0652 + IT_0658 + IT_0660 + IT_0661;
    const complex_t IT_2494 = IT_0647*IT_2493;
    const complex_t IT_2495 = IT_0760 + IT_0762 + IT_0766 + IT_0768 + IT_0770;
    const complex_t IT_2496 = IT_0756*IT_2495;
    const complex_t IT_2497 = IT_0779 + IT_0781 + IT_0785 + IT_0787 + IT_0788;
    const complex_t IT_2498 = IT_0775*IT_2497;
    const complex_t IT_2499 = IT_0797 + IT_0800 + IT_0802 + IT_0804 + IT_0805;
    const complex_t IT_2500 = IT_0794*IT_2499;
    const complex_t IT_2501 = IT_0813 + IT_0815 + IT_0819 + IT_0821 + IT_0822;
    const complex_t IT_2502 = IT_0811*IT_2501;
    const complex_t IT_2503 = IT_0843 + IT_0845 + IT_0849 + IT_0851 + IT_0852;
    const complex_t IT_2504 = IT_0839*IT_2503;
    const complex_t IT_2505 = IT_0873 + IT_0875 + IT_0879 + IT_0881 + IT_0882;
    const complex_t IT_2506 = IT_0869*IT_2505;
    const complex_t IT_2507 = IT_0934 + IT_1995 + IT_1998;
    const complex_t IT_2508 = -IT_0940 + -IT_1993 + -IT_1994;
    const complex_t IT_2509 = IT_2507 + IT_2508;
    const complex_t IT_2510 = IT_0932*IT_2509;
    const complex_t IT_2511 = IT_0984 + IT_2010 + IT_2013;
    const complex_t IT_2512 = -IT_0990 + -IT_2008 + -IT_2009;
    const complex_t IT_2513 = IT_2511 + IT_2512;
    const complex_t IT_2514 = IT_0982*IT_2513;
    const complex_t IT_2515 = IT_0997*IT_2509;
    const complex_t IT_2516 = IT_1016*IT_2513;
    const complex_t IT_2517 = IT_1051 + IT_2029 + IT_2032;
    const complex_t IT_2518 = -IT_1057 + -IT_2027 + -IT_2028;
    const complex_t IT_2519 = IT_2517 + IT_2518;
    const complex_t IT_2520 = IT_1049*IT_2519;
    const complex_t IT_2521 = IT_1074 + IT_2044 + IT_2047;
    const complex_t IT_2522 = -IT_1077 + -IT_2042 + -IT_2043;
    const complex_t IT_2523 = IT_2521 + IT_2522;
    const complex_t IT_2524 = IT_1072*IT_2523;
    const complex_t IT_2525 = IT_1087*IT_2519;
    const complex_t IT_2526 = IT_1106*IT_2523;
    const complex_t IT_2527 = IT_1208 + IT_1210 + IT_1214 + IT_1215 + IT_1216;
    const complex_t IT_2528 = IT_1204*IT_2527;
    const complex_t IT_2529 = IT_1223 + IT_1227 + IT_1229 + IT_1232 + IT_1233;
    const complex_t IT_2530 = IT_1221*IT_2529;
    const complex_t IT_2531 = IT_1240 + IT_1244 + IT_1246 + IT_1247 + IT_1250;
    const complex_t IT_2532 = IT_1238*IT_2531;
    const complex_t IT_2533 = IT_1257 + IT_1261 + IT_1263 + IT_1266 + IT_1267;
    const complex_t IT_2534 = IT_1255*IT_2533;
    const complex_t IT_2535 = IT_1290 + IT_1292 + IT_1296 + IT_1297 + IT_1300;
    const complex_t IT_2536 = IT_1288*IT_2535;
    const complex_t IT_2537 = IT_1322 + IT_1327 + IT_1329 + IT_1330 + IT_1333;
    const complex_t IT_2538 = IT_1321*IT_2537;
    const complex_t IT_2539 = IT_1361*IT_2483;
    const complex_t IT_2540 = IT_1364*IT_2489;
    const complex_t IT_2541 = IT_1367*IT_2487;
    const complex_t IT_2542 = IT_1370*IT_2485;
    const complex_t IT_2543 = IT_1373*IT_2491;
    const complex_t IT_2544 = IT_1376*IT_2493;
    const complex_t IT_2545 = IT_1460 + IT_1462 + IT_1466 + IT_1467 + IT_1470;
    const complex_t IT_2546 = IT_1458*IT_2545;
    const complex_t IT_2547 = IT_1479 + IT_1481 + IT_1483 + IT_1484 + IT_1487;
    const complex_t IT_2548 = IT_1475*IT_2547;
    const complex_t IT_2549 = IT_1495 + IT_1498 + IT_1501 + IT_1502 + IT_1504;
    const complex_t IT_2550 = IT_1492*IT_2549;
    const complex_t IT_2551 = IT_1513 + IT_1515 + IT_1519 + IT_1520 + IT_1521;
    const complex_t IT_2552 = IT_1509*IT_2551;
    const complex_t IT_2553 = IT_1528 + IT_1530 + IT_1534 + IT_1535 + IT_1538;
    const complex_t IT_2554 = IT_1526*IT_2553;
    const complex_t IT_2555 = IT_1545 + IT_1547 + IT_1551 + IT_1552 + IT_1555;
    const complex_t IT_2556 = IT_1543*IT_2555;
    const complex_t IT_2557 = IT_1566*IT_2527;
    const complex_t IT_2558 = IT_1569*IT_2533;
    const complex_t IT_2559 = IT_1572*IT_2531;
    const complex_t IT_2560 = IT_1575*IT_2529;
    const complex_t IT_2561 = IT_1578*IT_2535;
    const complex_t IT_2562 = IT_1581*IT_2537;
    const complex_t IT_2563 = IT_1608*IT_2545;
    const complex_t IT_2564 = IT_1611*IT_2551;
    const complex_t IT_2565 = IT_1614*IT_2549;
    const complex_t IT_2566 = IT_1617*IT_2547;
    const complex_t IT_2567 = IT_1620*IT_2553;
    const complex_t IT_2568 = IT_1623*IT_2555;
    const complex_t IT_2569 = IT_1667*IT_2495;
    const complex_t IT_2570 = IT_1670*IT_2501;
    const complex_t IT_2571 = IT_1673*IT_2499;
    const complex_t IT_2572 = IT_1676*IT_2497;
    const complex_t IT_2573 = IT_1679*IT_2503;
    const complex_t IT_2574 = IT_1682*IT_2505;
    const complex_t IT_2575 = IT_1695 + IT_2213 + IT_2216;
    const complex_t IT_2576 = -IT_1701 + -IT_2214 + -IT_2215;
    const complex_t IT_2577 = IT_2575 + IT_2576;
    const complex_t IT_2578 = IT_1693*IT_2577;
    const complex_t IT_2579 = IT_1721 + IT_2228 + IT_2231;
    const complex_t IT_2580 = -IT_1727 + -IT_2229 + -IT_2230;
    const complex_t IT_2581 = IT_2579 + IT_2580;
    const complex_t IT_2582 = IT_1716*IT_2581;
    const complex_t IT_2583 = IT_1731*IT_2577;
    const complex_t IT_2584 = IT_1734*IT_2581;
    const complex_t IT_2585 = IT_1753 + IT_2247 + IT_2250;
    const complex_t IT_2586 = -IT_1759 + -IT_2245 + -IT_2246;
    const complex_t IT_2587 = IT_2585 + IT_2586;
    const complex_t IT_2588 = IT_1751*IT_2587;
    const complex_t IT_2589 = IT_1776 + IT_2262 + IT_2265;
    const complex_t IT_2590 = -IT_1782 + -IT_2260 + -IT_2261;
    const complex_t IT_2591 = IT_2589 + IT_2590;
    const complex_t IT_2592 = IT_1774*IT_2591;
    const complex_t IT_2593 = IT_1789*IT_2587;
    const complex_t IT_2594 = IT_1792*IT_2591;
    const complex_t IT_2595 = (-2)*IT_1807 + IT_1812 + (-2)*IT_1815 + (-2)
      *IT_1818 + (-2)*IT_1821 + (-2)*IT_1824 + (-2)*IT_1827 + (-2)*IT_1830 + (-2
      )*IT_1833 + (-2)*IT_1836 + (-2)*IT_1839 + (-2)*IT_1842 + (-2)*IT_1845 + (
      -2)*IT_1848 + (-2)*IT_1851 + (-2)*IT_1854 + (-2)*IT_1857 + -IT_1862 + 
      -IT_1877 + IT_1890 + IT_1891 + -IT_1892 + -IT_1907 + IT_1920 + (-2)
      *IT_1927 + (-2)*IT_1934 + (-2)*IT_1945 + (-2)*IT_1960 + (-2)*IT_1967 + (-2
      )*IT_1982 + -IT_1991 + -IT_2006 + IT_2019 + IT_2020 + -IT_2025 + -IT_2040 
      + IT_2053 + IT_2054 + (-2)*IT_2081 + (-2)*IT_2084 + (-2)*IT_2087 + (-2)
      *IT_2096 + (-2)*IT_2099 + (-2)*IT_2102 + (-2)*IT_2105 + (-2)*IT_2108 + (-2
      )*IT_2111 + (-2)*IT_2138 + (-2)*IT_2141 + (-2)*IT_2144 + (-2)*IT_2147 + (
      -2)*IT_2150 + (-2)*IT_2155 + (-2)*IT_2159 + (-2)*IT_2163 + (-2)*IT_2167 + 
      (-2)*IT_2171 + (-2)*IT_2178 + (-2)*IT_2182 + (-2)*IT_2192 + (-2)*IT_2195 +
       (-2)*IT_2198 + (-2)*IT_2201 + (-2)*IT_2204 + -IT_2209 + -IT_2224 +
       IT_2237 + IT_2238 + -IT_2243 + -IT_2258 + IT_2271 + IT_2272 + IT_2466 +
       IT_2470 + -IT_2471 + -IT_2472 + IT_2476 + IT_2480 + -IT_2481 + -IT_2482 +
       2*IT_2484 + 2*IT_2486 + 2*IT_2488 + 2*IT_2490 + 2*IT_2492 + 2*IT_2494 + 2
      *IT_2496 + 2*IT_2498 + 2*IT_2500 + 2*IT_2502 + 2*IT_2504 + 2*IT_2506 +
       IT_2510 + IT_2514 + -IT_2515 + -IT_2516 + IT_2520 + IT_2524 + -IT_2525 + 
      -IT_2526 + 2*IT_2528 + 2*IT_2530 + 2*IT_2532 + 2*IT_2534 + 2*IT_2536 + 2
      *IT_2538 + 2*IT_2539 + 2*IT_2540 + 2*IT_2541 + 2*IT_2542 + 2*IT_2543 + 2
      *IT_2544 + 2*IT_2546 + 2*IT_2548 + 2*IT_2550 + 2*IT_2552 + 2*IT_2554 + 2
      *IT_2556 + 2*IT_2557 + 2*IT_2558 + 2*IT_2559 + 2*IT_2560 + 2*IT_2561 + 2
      *IT_2562 + 2*IT_2563 + 2*IT_2564 + 2*IT_2565 + 2*IT_2566 + 2*IT_2567 + 2
      *IT_2568 + 2*IT_2569 + 2*IT_2570 + 2*IT_2571 + 2*IT_2572 + 2*IT_2573 + 2
      *IT_2574 + IT_2578 + IT_2582 + -IT_2583 + -IT_2584 + IT_2588 + IT_2592 + 
      -IT_2593 + -IT_2594;
    const complex_t IT_2596 = IT_1801*IT_2595;
    const complex_t IT_2597 = (-0.5)*IT_2596;
    const complex_t IT_2598 = 2*IT_2597;
    const complex_t IT_2599 = IT_0005*IT_2598;
    return (complex_t{0, 0.25})*IT_1804 + (complex_t{0, 0.25})*IT_2277 + 
      (complex_t{0, (-0.25)})*IT_2462 + (complex_t{0, (-0.25)})*IT_2599;
}
} // End of namespace c9_nmfv
