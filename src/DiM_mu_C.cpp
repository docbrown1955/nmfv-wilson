#include <quadmath.h>
#include "clooptools.h"
#include "marty/looptools_init.h"
#include "marty/looptools_quad_extension.h"
#include <cmath>
#include "marty/looptools_interface.h"
#include "DiM_mu_C.h"
#include "common.h"

#include "params.h"
#include "libcomplexop.h"

namespace c9_nmfv {

complex_t DiM_mu_C(
        param_t const &param
        )
{
    clearcache();
    const real_t M_W = param.M_W;
    const real_t beta = param.beta;
    const real_t e_em = param.e_em;
    const real_t m_mu = param.m_mu;
    const real_t s_12 = param.s_12;
    const real_t m_C1p = param.m_C1p;
    const real_t m_C2p = param.m_C2p;
    const real_t m_snu_e = param.m_snu_e;
    const real_t theta_W = param.theta_W;
    const real_t m_snu_mu = param.m_snu_mu;
    const real_t m_snu_tau = param.m_snu_tau;
    const complex_t U_d1 = param.U_d1;
    const complex_t U_d2 = param.U_d2;
    const complex_t V_u1 = param.V_u1;
    const complex_t V_u2 = param.V_u2;
    const complex_t U_Wm1 = param.U_Wm1;
    const complex_t U_Wm2 = param.U_Wm2;
    const complex_t V_Wp1 = param.V_Wp1;
    const complex_t V_Wp2 = param.V_Wp2;
    const complex_t U_snu_10 = param.U_snu_10;
    const complex_t U_snu_11 = param.U_snu_11;
    const complex_t U_snu_12 = param.U_snu_12;
    const complex_t IT_0000 = powq(e_em, -1);
    const complex_t IT_0001 = m_mu*IT_0000;
    const complex_t IT_0002 = cosq(theta_W);
    const complex_t IT_0003 = cpowq(IT_0002, -1);
    const complex_t IT_0004 = tanq(theta_W);
    const complex_t IT_0005 = cpowq(IT_0004, 2);
    const complex_t IT_0006 = cpowq(1 + IT_0005, (-0.5));
    const complex_t IT_0007 = IT_0003*IT_0006;
    const complex_t IT_0008 = e_em*U_Wm1*conjq(U_Wm2);
    const complex_t IT_0009 = IT_0007*IT_0008;
    const complex_t IT_0010 = U_d1*conjq(U_d2)*e_em;
    const complex_t IT_0011 = IT_0007*IT_0010;
    const complex_t IT_0012 = (complex_t{0, 1})*(IT_0009 + IT_0011);
    const complex_t IT_0013 = -IT_0012;
    const complex_t IT_0014 = sinq(theta_W);
    const complex_t IT_0015 = cpowq(IT_0014, -1);
    const complex_t IT_0016 = (complex_t{0, 1})*e_em*V_Wp1*IT_0015*conjq
      (U_snu_10);
    const complex_t IT_0017 = powq(M_W, -1);
    const complex_t IT_0018 = cosq(beta);
    const complex_t IT_0019 = cpowq(IT_0018, -1);
    const complex_t IT_0020 = (complex_t{0, 1.4142135623731})*U_d2*e_em*m_mu
      *IT_0015*IT_0017*IT_0019*U_snu_10;
    const complex_t IT_0021 = 0.5*IT_0020;
    const complex_t IT_0022 = IT_0013*IT_0016*IT_0021;
    const complex_t IT_0023 = 0.101321183642338*m_C1p;
    const complex_t IT_0024 = IT_0022*IT_0023;
    const complex_t IT_0025 = powq(m_mu, 2);
    const complex_t IT_0026 = powq(m_C1p, 2);
    const complex_t IT_0027 = powq(m_C2p, 2);
    const complex_t IT_0028 = powq(m_snu_e, 2);
    const complex_t IT_0029 = mty::lt::C0iC(3, 2*IT_0025 + (-2)*s_12, IT_0025,
       IT_0025, IT_0026, IT_0027, IT_0028, mty::lt::reg_int);
    const complex_t IT_0030 = mty::lt::C0iC(6, 2*IT_0025 + (-2)*s_12, IT_0025,
       IT_0025, IT_0026, IT_0027, IT_0028, mty::lt::reg_int);
    const complex_t IT_0031 = mty::lt::C0iC(0, 2*IT_0025 + (-2)*s_12, IT_0025,
       IT_0025, IT_0026, IT_0027, IT_0028, mty::lt::reg_int);
    const complex_t IT_0032 = IT_0029 + IT_0030 + IT_0031;
    const complex_t IT_0033 = IT_0024*IT_0032;
    const complex_t IT_0034 = (complex_t{0, 1})*e_em*V_Wp1*IT_0015*conjq
      (U_snu_11);
    const complex_t IT_0035 = (complex_t{0, 1.4142135623731})*U_d2*e_em*m_mu
      *IT_0015*IT_0017*IT_0019*U_snu_11;
    const complex_t IT_0036 = 0.5*IT_0035;
    const complex_t IT_0037 = IT_0013*IT_0034*IT_0036;
    const complex_t IT_0038 = IT_0023*IT_0037;
    const complex_t IT_0039 = powq(m_snu_mu, 2);
    const complex_t IT_0040 = mty::lt::C0iC(3, 2*IT_0025 + (-2)*s_12, IT_0025,
       IT_0025, IT_0026, IT_0027, IT_0039, mty::lt::reg_int);
    const complex_t IT_0041 = mty::lt::C0iC(6, 2*IT_0025 + (-2)*s_12, IT_0025,
       IT_0025, IT_0026, IT_0027, IT_0039, mty::lt::reg_int);
    const complex_t IT_0042 = mty::lt::C0iC(0, 2*IT_0025 + (-2)*s_12, IT_0025,
       IT_0025, IT_0026, IT_0027, IT_0039, mty::lt::reg_int);
    const complex_t IT_0043 = IT_0040 + IT_0041 + IT_0042;
    const complex_t IT_0044 = IT_0038*IT_0043;
    const complex_t IT_0045 = (complex_t{0, 1})*e_em*V_Wp1*IT_0015*conjq
      (U_snu_12);
    const complex_t IT_0046 = (complex_t{0, 1.4142135623731})*U_d2*e_em*m_mu
      *IT_0015*IT_0017*IT_0019*U_snu_12;
    const complex_t IT_0047 = 0.5*IT_0046;
    const complex_t IT_0048 = IT_0013*IT_0045*IT_0047;
    const complex_t IT_0049 = IT_0023*IT_0048;
    const complex_t IT_0050 = powq(m_snu_tau, 2);
    const complex_t IT_0051 = mty::lt::C0iC(3, 2*IT_0025 + (-2)*s_12, IT_0025,
       IT_0025, IT_0026, IT_0027, IT_0050, mty::lt::reg_int);
    const complex_t IT_0052 = mty::lt::C0iC(6, 2*IT_0025 + (-2)*s_12, IT_0025,
       IT_0025, IT_0026, IT_0027, IT_0050, mty::lt::reg_int);
    const complex_t IT_0053 = mty::lt::C0iC(0, 2*IT_0025 + (-2)*s_12, IT_0025,
       IT_0025, IT_0026, IT_0027, IT_0050, mty::lt::reg_int);
    const complex_t IT_0054 = IT_0051 + IT_0052 + IT_0053;
    const complex_t IT_0055 = IT_0049*IT_0054;
    const complex_t IT_0056 = e_em*U_Wm1*conjq(U_Wm1);
    const complex_t IT_0057 = IT_0007*IT_0056;
    const complex_t IT_0058 = U_d1*conjq(U_d1)*e_em;
    const complex_t IT_0059 = IT_0007*IT_0058;
    const complex_t IT_0060 = (complex_t{0, 1})*(IT_0057 + IT_0059);
    const complex_t IT_0061 = -IT_0060;
    const complex_t IT_0062 = (complex_t{0, 1.4142135623731})*U_d1*e_em*m_mu
      *IT_0015*IT_0017*IT_0019*U_snu_10;
    const complex_t IT_0063 = 0.5*IT_0062;
    const complex_t IT_0064 = IT_0016*IT_0061*IT_0063;
    const complex_t IT_0065 = IT_0023*IT_0064;
    const complex_t IT_0066 = mty::lt::C0iC(3, 2*IT_0025 + (-2)*s_12, IT_0025,
       IT_0025, IT_0026, IT_0026, IT_0028, mty::lt::reg_int);
    const complex_t IT_0067 = mty::lt::C0iC(6, 2*IT_0025 + (-2)*s_12, IT_0025,
       IT_0025, IT_0026, IT_0026, IT_0028, mty::lt::reg_int);
    const complex_t IT_0068 = mty::lt::C0iC(0, 2*IT_0025 + (-2)*s_12, IT_0025,
       IT_0025, IT_0026, IT_0026, IT_0028, mty::lt::reg_int);
    const complex_t IT_0069 = IT_0066 + IT_0067 + IT_0068;
    const complex_t IT_0070 = IT_0065*IT_0069;
    const complex_t IT_0071 = (complex_t{0, 1.4142135623731})*U_d1*e_em*m_mu
      *IT_0015*IT_0017*IT_0019*U_snu_11;
    const complex_t IT_0072 = 0.5*IT_0071;
    const complex_t IT_0073 = IT_0034*IT_0061*IT_0072;
    const complex_t IT_0074 = IT_0023*IT_0073;
    const complex_t IT_0075 = mty::lt::C0iC(3, 2*IT_0025 + (-2)*s_12, IT_0025,
       IT_0025, IT_0026, IT_0026, IT_0039, mty::lt::reg_int);
    const complex_t IT_0076 = mty::lt::C0iC(6, 2*IT_0025 + (-2)*s_12, IT_0025,
       IT_0025, IT_0026, IT_0026, IT_0039, mty::lt::reg_int);
    const complex_t IT_0077 = mty::lt::C0iC(0, 2*IT_0025 + (-2)*s_12, IT_0025,
       IT_0025, IT_0026, IT_0026, IT_0039, mty::lt::reg_int);
    const complex_t IT_0078 = IT_0075 + IT_0076 + IT_0077;
    const complex_t IT_0079 = IT_0074*IT_0078;
    const complex_t IT_0080 = (complex_t{0, 1.4142135623731})*U_d1*e_em*m_mu
      *IT_0015*IT_0017*IT_0019*U_snu_12;
    const complex_t IT_0081 = 0.5*IT_0080;
    const complex_t IT_0082 = IT_0045*IT_0061*IT_0081;
    const complex_t IT_0083 = IT_0023*IT_0082;
    const complex_t IT_0084 = mty::lt::C0iC(3, 2*IT_0025 + (-2)*s_12, IT_0025,
       IT_0025, IT_0026, IT_0026, IT_0050, mty::lt::reg_int);
    const complex_t IT_0085 = mty::lt::C0iC(6, 2*IT_0025 + (-2)*s_12, IT_0025,
       IT_0025, IT_0026, IT_0026, IT_0050, mty::lt::reg_int);
    const complex_t IT_0086 = mty::lt::C0iC(0, 2*IT_0025 + (-2)*s_12, IT_0025,
       IT_0025, IT_0026, IT_0026, IT_0050, mty::lt::reg_int);
    const complex_t IT_0087 = IT_0084 + IT_0085 + IT_0086;
    const complex_t IT_0088 = IT_0083*IT_0087;
    const complex_t IT_0089 = 0.101321183642338*m_C2p;
    const complex_t IT_0090 = e_em*V_Wp2*conjq(V_Wp2);
    const complex_t IT_0091 = IT_0007*IT_0090;
    const complex_t IT_0092 = V_u2*conjq(V_u2)*e_em;
    const complex_t IT_0093 = IT_0007*IT_0092;
    const complex_t IT_0094 = (complex_t{0, 1})*(IT_0091 + IT_0093);
    const complex_t IT_0095 = (complex_t{0, 1.4142135623731})*conjq(U_d2)*e_em
      *m_mu*IT_0015*IT_0017*IT_0019*conjq(U_snu_10);
    const complex_t IT_0096 = 0.5*IT_0095;
    const complex_t IT_0097 = (complex_t{0, 1})*e_em*conjq(V_Wp2)*IT_0015
      *U_snu_10;
    const complex_t IT_0098 = IT_0094*IT_0096*IT_0097;
    const complex_t IT_0099 = IT_0089*IT_0098;
    const complex_t IT_0100 = mty::lt::C0iC(0, 2*IT_0025 + (-2)*s_12, IT_0025,
       IT_0025, IT_0027, IT_0027, IT_0028, mty::lt::reg_int);
    const complex_t IT_0101 = mty::lt::C0iC(3, 2*IT_0025 + (-2)*s_12, IT_0025,
       IT_0025, IT_0027, IT_0027, IT_0028, mty::lt::reg_int);
    const complex_t IT_0102 = mty::lt::C0iC(6, 2*IT_0025 + (-2)*s_12, IT_0025,
       IT_0025, IT_0027, IT_0027, IT_0028, mty::lt::reg_int);
    const complex_t IT_0103 = IT_0100 + IT_0101 + IT_0102;
    const complex_t IT_0104 = IT_0099*IT_0103;
    const complex_t IT_0105 = (complex_t{0, 1.4142135623731})*conjq(U_d2)*e_em
      *m_mu*IT_0015*IT_0017*IT_0019*conjq(U_snu_11);
    const complex_t IT_0106 = 0.5*IT_0105;
    const complex_t IT_0107 = (complex_t{0, 1})*e_em*conjq(V_Wp2)*IT_0015
      *U_snu_11;
    const complex_t IT_0108 = IT_0094*IT_0106*IT_0107;
    const complex_t IT_0109 = IT_0089*IT_0108;
    const complex_t IT_0110 = mty::lt::C0iC(0, 2*IT_0025 + (-2)*s_12, IT_0025,
       IT_0025, IT_0027, IT_0027, IT_0039, mty::lt::reg_int);
    const complex_t IT_0111 = mty::lt::C0iC(3, 2*IT_0025 + (-2)*s_12, IT_0025,
       IT_0025, IT_0027, IT_0027, IT_0039, mty::lt::reg_int);
    const complex_t IT_0112 = mty::lt::C0iC(6, 2*IT_0025 + (-2)*s_12, IT_0025,
       IT_0025, IT_0027, IT_0027, IT_0039, mty::lt::reg_int);
    const complex_t IT_0113 = IT_0110 + IT_0111 + IT_0112;
    const complex_t IT_0114 = IT_0109*IT_0113;
    const complex_t IT_0115 = (complex_t{0, 1.4142135623731})*conjq(U_d2)*e_em
      *m_mu*IT_0015*IT_0017*IT_0019*conjq(U_snu_12);
    const complex_t IT_0116 = 0.5*IT_0115;
    const complex_t IT_0117 = (complex_t{0, 1})*e_em*conjq(V_Wp2)*IT_0015
      *U_snu_12;
    const complex_t IT_0118 = IT_0094*IT_0116*IT_0117;
    const complex_t IT_0119 = IT_0089*IT_0118;
    const complex_t IT_0120 = mty::lt::C0iC(0, 2*IT_0025 + (-2)*s_12, IT_0025,
       IT_0025, IT_0027, IT_0027, IT_0050, mty::lt::reg_int);
    const complex_t IT_0121 = mty::lt::C0iC(3, 2*IT_0025 + (-2)*s_12, IT_0025,
       IT_0025, IT_0027, IT_0027, IT_0050, mty::lt::reg_int);
    const complex_t IT_0122 = mty::lt::C0iC(6, 2*IT_0025 + (-2)*s_12, IT_0025,
       IT_0025, IT_0027, IT_0027, IT_0050, mty::lt::reg_int);
    const complex_t IT_0123 = IT_0120 + IT_0121 + IT_0122;
    const complex_t IT_0124 = IT_0119*IT_0123;
    const complex_t IT_0125 = e_em*conjq(V_Wp1)*V_Wp2;
    const complex_t IT_0126 = IT_0007*IT_0125;
    const complex_t IT_0127 = conjq(V_u1)*V_u2*e_em;
    const complex_t IT_0128 = IT_0007*IT_0127;
    const complex_t IT_0129 = (complex_t{0, 1})*(IT_0126 + IT_0128);
    const complex_t IT_0130 = (complex_t{0, 1.4142135623731})*conjq(U_d1)*e_em
      *m_mu*IT_0015*IT_0017*IT_0019*conjq(U_snu_10);
    const complex_t IT_0131 = 0.5*IT_0130;
    const complex_t IT_0132 = IT_0097*IT_0129*IT_0131;
    const complex_t IT_0133 = IT_0023*IT_0132;
    const complex_t IT_0134 = IT_0032*IT_0133;
    const complex_t IT_0135 = (complex_t{0, 1.4142135623731})*conjq(U_d1)*e_em
      *m_mu*IT_0015*IT_0017*IT_0019*conjq(U_snu_11);
    const complex_t IT_0136 = 0.5*IT_0135;
    const complex_t IT_0137 = IT_0107*IT_0129*IT_0136;
    const complex_t IT_0138 = IT_0023*IT_0137;
    const complex_t IT_0139 = IT_0043*IT_0138;
    const complex_t IT_0140 = (complex_t{0, 1.4142135623731})*conjq(U_d1)*e_em
      *m_mu*IT_0015*IT_0017*IT_0019*conjq(U_snu_12);
    const complex_t IT_0141 = 0.5*IT_0140;
    const complex_t IT_0142 = IT_0117*IT_0129*IT_0141;
    const complex_t IT_0143 = IT_0023*IT_0142;
    const complex_t IT_0144 = IT_0054*IT_0143;
    const complex_t IT_0145 = e_em*V_Wp1*conjq(V_Wp2);
    const complex_t IT_0146 = IT_0007*IT_0145;
    const complex_t IT_0147 = V_u1*conjq(V_u2)*e_em;
    const complex_t IT_0148 = IT_0007*IT_0147;
    const complex_t IT_0149 = (complex_t{0, 1})*(IT_0146 + IT_0148);
    const complex_t IT_0150 = (complex_t{0, 1})*e_em*conjq(V_Wp1)*IT_0015
      *U_snu_10;
    const complex_t IT_0151 = IT_0096*IT_0149*IT_0150;
    const complex_t IT_0152 = IT_0089*IT_0151;
    const complex_t IT_0153 = mty::lt::C0iC(0, 2*IT_0025 + (-2)*s_12, IT_0025,
       IT_0025, IT_0027, IT_0026, IT_0028, mty::lt::reg_int);
    const complex_t IT_0154 = mty::lt::C0iC(3, 2*IT_0025 + (-2)*s_12, IT_0025,
       IT_0025, IT_0027, IT_0026, IT_0028, mty::lt::reg_int);
    const complex_t IT_0155 = mty::lt::C0iC(6, 2*IT_0025 + (-2)*s_12, IT_0025,
       IT_0025, IT_0027, IT_0026, IT_0028, mty::lt::reg_int);
    const complex_t IT_0156 = IT_0153 + IT_0154 + IT_0155;
    const complex_t IT_0157 = IT_0152*IT_0156;
    const complex_t IT_0158 = (complex_t{0, 1})*e_em*conjq(V_Wp1)*IT_0015
      *U_snu_11;
    const complex_t IT_0159 = IT_0106*IT_0149*IT_0158;
    const complex_t IT_0160 = IT_0089*IT_0159;
    const complex_t IT_0161 = mty::lt::C0iC(0, 2*IT_0025 + (-2)*s_12, IT_0025,
       IT_0025, IT_0027, IT_0026, IT_0039, mty::lt::reg_int);
    const complex_t IT_0162 = mty::lt::C0iC(3, 2*IT_0025 + (-2)*s_12, IT_0025,
       IT_0025, IT_0027, IT_0026, IT_0039, mty::lt::reg_int);
    const complex_t IT_0163 = mty::lt::C0iC(6, 2*IT_0025 + (-2)*s_12, IT_0025,
       IT_0025, IT_0027, IT_0026, IT_0039, mty::lt::reg_int);
    const complex_t IT_0164 = IT_0161 + IT_0162 + IT_0163;
    const complex_t IT_0165 = IT_0160*IT_0164;
    const complex_t IT_0166 = (complex_t{0, 1})*e_em*conjq(V_Wp1)*IT_0015
      *U_snu_12;
    const complex_t IT_0167 = IT_0116*IT_0149*IT_0166;
    const complex_t IT_0168 = IT_0089*IT_0167;
    const complex_t IT_0169 = mty::lt::C0iC(0, 2*IT_0025 + (-2)*s_12, IT_0025,
       IT_0025, IT_0027, IT_0026, IT_0050, mty::lt::reg_int);
    const complex_t IT_0170 = mty::lt::C0iC(3, 2*IT_0025 + (-2)*s_12, IT_0025,
       IT_0025, IT_0027, IT_0026, IT_0050, mty::lt::reg_int);
    const complex_t IT_0171 = mty::lt::C0iC(6, 2*IT_0025 + (-2)*s_12, IT_0025,
       IT_0025, IT_0027, IT_0026, IT_0050, mty::lt::reg_int);
    const complex_t IT_0172 = IT_0169 + IT_0170 + IT_0171;
    const complex_t IT_0173 = IT_0168*IT_0172;
    const complex_t IT_0174 = e_em*U_Wm2*conjq(U_Wm2);
    const complex_t IT_0175 = IT_0007*IT_0174;
    const complex_t IT_0176 = U_d2*conjq(U_d2)*e_em;
    const complex_t IT_0177 = IT_0007*IT_0176;
    const complex_t IT_0178 = (complex_t{0, 1})*(IT_0175 + IT_0177);
    const complex_t IT_0179 = -IT_0178;
    const complex_t IT_0180 = (complex_t{0, 1})*e_em*V_Wp2*IT_0015*conjq
      (U_snu_10);
    const complex_t IT_0181 = IT_0021*IT_0179*IT_0180;
    const complex_t IT_0182 = IT_0089*IT_0181;
    const complex_t IT_0183 = IT_0103*IT_0182;
    const complex_t IT_0184 = (complex_t{0, 1})*e_em*V_Wp2*IT_0015*conjq
      (U_snu_11);
    const complex_t IT_0185 = IT_0036*IT_0179*IT_0184;
    const complex_t IT_0186 = IT_0089*IT_0185;
    const complex_t IT_0187 = IT_0113*IT_0186;
    const complex_t IT_0188 = (complex_t{0, 1})*e_em*V_Wp2*IT_0015*conjq
      (U_snu_12);
    const complex_t IT_0189 = IT_0047*IT_0179*IT_0188;
    const complex_t IT_0190 = IT_0089*IT_0189;
    const complex_t IT_0191 = IT_0123*IT_0190;
    const complex_t IT_0192 = e_em*conjq(U_Wm1)*U_Wm2;
    const complex_t IT_0193 = IT_0007*IT_0192;
    const complex_t IT_0194 = conjq(U_d1)*U_d2*e_em;
    const complex_t IT_0195 = IT_0007*IT_0194;
    const complex_t IT_0196 = (complex_t{0, 1})*(IT_0193 + IT_0195);
    const complex_t IT_0197 = -IT_0196;
    const complex_t IT_0198 = IT_0063*IT_0180*IT_0197;
    const complex_t IT_0199 = IT_0089*IT_0198;
    const complex_t IT_0200 = IT_0156*IT_0199;
    const complex_t IT_0201 = IT_0072*IT_0184*IT_0197;
    const complex_t IT_0202 = IT_0089*IT_0201;
    const complex_t IT_0203 = IT_0164*IT_0202;
    const complex_t IT_0204 = IT_0081*IT_0188*IT_0197;
    const complex_t IT_0205 = IT_0089*IT_0204;
    const complex_t IT_0206 = IT_0172*IT_0205;
    const complex_t IT_0207 = e_em*V_Wp1*conjq(V_Wp1);
    const complex_t IT_0208 = IT_0007*IT_0207;
    const complex_t IT_0209 = V_u1*conjq(V_u1)*e_em;
    const complex_t IT_0210 = IT_0007*IT_0209;
    const complex_t IT_0211 = (complex_t{0, 1})*(IT_0208 + IT_0210);
    const complex_t IT_0212 = IT_0131*IT_0150*IT_0211;
    const complex_t IT_0213 = IT_0023*IT_0212;
    const complex_t IT_0214 = IT_0069*IT_0213;
    const complex_t IT_0215 = IT_0136*IT_0158*IT_0211;
    const complex_t IT_0216 = IT_0023*IT_0215;
    const complex_t IT_0217 = IT_0078*IT_0216;
    const complex_t IT_0218 = IT_0141*IT_0166*IT_0211;
    const complex_t IT_0219 = IT_0023*IT_0218;
    const complex_t IT_0220 = IT_0087*IT_0219;
    const complex_t IT_0221 = IT_0033 + IT_0044 + IT_0055 + IT_0070 + IT_0079 
      + IT_0088 + IT_0104 + IT_0114 + IT_0124 + IT_0134 + IT_0139 + IT_0144 +
       IT_0157 + IT_0165 + IT_0173 + IT_0183 + IT_0187 + IT_0191 + IT_0200 +
       IT_0203 + IT_0206 + IT_0214 + IT_0217 + IT_0220;
    const complex_t IT_0222 = IT_0013*IT_0021*IT_0131;
    const complex_t IT_0223 = 0.101321183642338*IT_0222;
    const complex_t IT_0224 = m_mu*IT_0029;
    const complex_t IT_0225 = mty::lt::C0iC(12, 2*IT_0025 + (-2)*s_12,
       IT_0025, IT_0025, IT_0026, IT_0027, IT_0028, mty::lt::reg_int);
    const complex_t IT_0226 = m_mu*IT_0225;
    const complex_t IT_0227 = IT_0224 + IT_0226;
    const complex_t IT_0228 = m_mu*IT_0030;
    const complex_t IT_0229 = mty::lt::C0iC(15, 2*IT_0025 + (-2)*s_12,
       IT_0025, IT_0025, IT_0026, IT_0027, IT_0028, mty::lt::reg_int);
    const complex_t IT_0230 = m_mu*IT_0229;
    const complex_t IT_0231 = mty::lt::C0iC(18, 2*IT_0025 + (-2)*s_12,
       IT_0025, IT_0025, IT_0026, IT_0027, IT_0028, mty::lt::reg_int);
    const complex_t IT_0232 = m_mu*IT_0231;
    const complex_t IT_0233 = 0.5*IT_0228 + 1.5*IT_0230 + 0.5*IT_0232;
    const complex_t IT_0234 = IT_0227 + IT_0233;
    const complex_t IT_0235 = IT_0223*IT_0234;
    const complex_t IT_0236 = IT_0013*IT_0036*IT_0136;
    const complex_t IT_0237 = 0.101321183642338*IT_0236;
    const complex_t IT_0238 = m_mu*IT_0040;
    const complex_t IT_0239 = mty::lt::C0iC(12, 2*IT_0025 + (-2)*s_12,
       IT_0025, IT_0025, IT_0026, IT_0027, IT_0039, mty::lt::reg_int);
    const complex_t IT_0240 = m_mu*IT_0239;
    const complex_t IT_0241 = IT_0238 + IT_0240;
    const complex_t IT_0242 = m_mu*IT_0041;
    const complex_t IT_0243 = mty::lt::C0iC(15, 2*IT_0025 + (-2)*s_12,
       IT_0025, IT_0025, IT_0026, IT_0027, IT_0039, mty::lt::reg_int);
    const complex_t IT_0244 = m_mu*IT_0243;
    const complex_t IT_0245 = mty::lt::C0iC(18, 2*IT_0025 + (-2)*s_12,
       IT_0025, IT_0025, IT_0026, IT_0027, IT_0039, mty::lt::reg_int);
    const complex_t IT_0246 = m_mu*IT_0245;
    const complex_t IT_0247 = 0.5*IT_0242 + 1.5*IT_0244 + 0.5*IT_0246;
    const complex_t IT_0248 = IT_0241 + IT_0247;
    const complex_t IT_0249 = IT_0237*IT_0248;
    const complex_t IT_0250 = IT_0013*IT_0047*IT_0141;
    const complex_t IT_0251 = 0.101321183642338*IT_0250;
    const complex_t IT_0252 = m_mu*IT_0051;
    const complex_t IT_0253 = mty::lt::C0iC(12, 2*IT_0025 + (-2)*s_12,
       IT_0025, IT_0025, IT_0026, IT_0027, IT_0050, mty::lt::reg_int);
    const complex_t IT_0254 = m_mu*IT_0253;
    const complex_t IT_0255 = IT_0252 + IT_0254;
    const complex_t IT_0256 = m_mu*IT_0052;
    const complex_t IT_0257 = mty::lt::C0iC(15, 2*IT_0025 + (-2)*s_12,
       IT_0025, IT_0025, IT_0026, IT_0027, IT_0050, mty::lt::reg_int);
    const complex_t IT_0258 = m_mu*IT_0257;
    const complex_t IT_0259 = mty::lt::C0iC(18, 2*IT_0025 + (-2)*s_12,
       IT_0025, IT_0025, IT_0026, IT_0027, IT_0050, mty::lt::reg_int);
    const complex_t IT_0260 = m_mu*IT_0259;
    const complex_t IT_0261 = 0.5*IT_0256 + 1.5*IT_0258 + 0.5*IT_0260;
    const complex_t IT_0262 = IT_0255 + IT_0261;
    const complex_t IT_0263 = IT_0251*IT_0262;
    const complex_t IT_0264 = IT_0061*IT_0063*IT_0131;
    const complex_t IT_0265 = 0.101321183642338*IT_0264;
    const complex_t IT_0266 = m_mu*IT_0066;
    const complex_t IT_0267 = mty::lt::C0iC(12, 2*IT_0025 + (-2)*s_12,
       IT_0025, IT_0025, IT_0026, IT_0026, IT_0028, mty::lt::reg_int);
    const complex_t IT_0268 = m_mu*IT_0267;
    const complex_t IT_0269 = IT_0266 + IT_0268;
    const complex_t IT_0270 = m_mu*IT_0067;
    const complex_t IT_0271 = mty::lt::C0iC(15, 2*IT_0025 + (-2)*s_12,
       IT_0025, IT_0025, IT_0026, IT_0026, IT_0028, mty::lt::reg_int);
    const complex_t IT_0272 = m_mu*IT_0271;
    const complex_t IT_0273 = mty::lt::C0iC(18, 2*IT_0025 + (-2)*s_12,
       IT_0025, IT_0025, IT_0026, IT_0026, IT_0028, mty::lt::reg_int);
    const complex_t IT_0274 = m_mu*IT_0273;
    const complex_t IT_0275 = 0.5*IT_0270 + 1.5*IT_0272 + 0.5*IT_0274;
    const complex_t IT_0276 = IT_0269 + IT_0275;
    const complex_t IT_0277 = IT_0265*IT_0276;
    const complex_t IT_0278 = IT_0061*IT_0072*IT_0136;
    const complex_t IT_0279 = 0.101321183642338*IT_0278;
    const complex_t IT_0280 = m_mu*IT_0075;
    const complex_t IT_0281 = mty::lt::C0iC(12, 2*IT_0025 + (-2)*s_12,
       IT_0025, IT_0025, IT_0026, IT_0026, IT_0039, mty::lt::reg_int);
    const complex_t IT_0282 = m_mu*IT_0281;
    const complex_t IT_0283 = IT_0280 + IT_0282;
    const complex_t IT_0284 = m_mu*IT_0076;
    const complex_t IT_0285 = mty::lt::C0iC(15, 2*IT_0025 + (-2)*s_12,
       IT_0025, IT_0025, IT_0026, IT_0026, IT_0039, mty::lt::reg_int);
    const complex_t IT_0286 = m_mu*IT_0285;
    const complex_t IT_0287 = mty::lt::C0iC(18, 2*IT_0025 + (-2)*s_12,
       IT_0025, IT_0025, IT_0026, IT_0026, IT_0039, mty::lt::reg_int);
    const complex_t IT_0288 = m_mu*IT_0287;
    const complex_t IT_0289 = 0.5*IT_0284 + 1.5*IT_0286 + 0.5*IT_0288;
    const complex_t IT_0290 = IT_0283 + IT_0289;
    const complex_t IT_0291 = IT_0279*IT_0290;
    const complex_t IT_0292 = IT_0061*IT_0081*IT_0141;
    const complex_t IT_0293 = 0.101321183642338*IT_0292;
    const complex_t IT_0294 = m_mu*IT_0084;
    const complex_t IT_0295 = mty::lt::C0iC(12, 2*IT_0025 + (-2)*s_12,
       IT_0025, IT_0025, IT_0026, IT_0026, IT_0050, mty::lt::reg_int);
    const complex_t IT_0296 = m_mu*IT_0295;
    const complex_t IT_0297 = IT_0294 + IT_0296;
    const complex_t IT_0298 = m_mu*IT_0085;
    const complex_t IT_0299 = mty::lt::C0iC(15, 2*IT_0025 + (-2)*s_12,
       IT_0025, IT_0025, IT_0026, IT_0026, IT_0050, mty::lt::reg_int);
    const complex_t IT_0300 = m_mu*IT_0299;
    const complex_t IT_0301 = mty::lt::C0iC(18, 2*IT_0025 + (-2)*s_12,
       IT_0025, IT_0025, IT_0026, IT_0026, IT_0050, mty::lt::reg_int);
    const complex_t IT_0302 = m_mu*IT_0301;
    const complex_t IT_0303 = 0.5*IT_0298 + 1.5*IT_0300 + 0.5*IT_0302;
    const complex_t IT_0304 = IT_0297 + IT_0303;
    const complex_t IT_0305 = IT_0293*IT_0304;
    const complex_t IT_0306 = IT_0094*IT_0097*IT_0180;
    const complex_t IT_0307 = 0.101321183642338*IT_0306;
    const complex_t IT_0308 = m_mu*IT_0101;
    const complex_t IT_0309 = mty::lt::C0iC(12, 2*IT_0025 + (-2)*s_12,
       IT_0025, IT_0025, IT_0027, IT_0027, IT_0028, mty::lt::reg_int);
    const complex_t IT_0310 = m_mu*IT_0309;
    const complex_t IT_0311 = IT_0308 + IT_0310;
    const complex_t IT_0312 = m_mu*IT_0102;
    const complex_t IT_0313 = mty::lt::C0iC(15, 2*IT_0025 + (-2)*s_12,
       IT_0025, IT_0025, IT_0027, IT_0027, IT_0028, mty::lt::reg_int);
    const complex_t IT_0314 = m_mu*IT_0313;
    const complex_t IT_0315 = mty::lt::C0iC(18, 2*IT_0025 + (-2)*s_12,
       IT_0025, IT_0025, IT_0027, IT_0027, IT_0028, mty::lt::reg_int);
    const complex_t IT_0316 = m_mu*IT_0315;
    const complex_t IT_0317 = 0.5*IT_0312 + 1.5*IT_0314 + 0.5*IT_0316;
    const complex_t IT_0318 = IT_0311 + IT_0317;
    const complex_t IT_0319 = IT_0307*IT_0318;
    const complex_t IT_0320 = IT_0094*IT_0107*IT_0184;
    const complex_t IT_0321 = 0.101321183642338*IT_0320;
    const complex_t IT_0322 = m_mu*IT_0111;
    const complex_t IT_0323 = mty::lt::C0iC(12, 2*IT_0025 + (-2)*s_12,
       IT_0025, IT_0025, IT_0027, IT_0027, IT_0039, mty::lt::reg_int);
    const complex_t IT_0324 = m_mu*IT_0323;
    const complex_t IT_0325 = IT_0322 + IT_0324;
    const complex_t IT_0326 = m_mu*IT_0112;
    const complex_t IT_0327 = mty::lt::C0iC(15, 2*IT_0025 + (-2)*s_12,
       IT_0025, IT_0025, IT_0027, IT_0027, IT_0039, mty::lt::reg_int);
    const complex_t IT_0328 = m_mu*IT_0327;
    const complex_t IT_0329 = mty::lt::C0iC(18, 2*IT_0025 + (-2)*s_12,
       IT_0025, IT_0025, IT_0027, IT_0027, IT_0039, mty::lt::reg_int);
    const complex_t IT_0330 = m_mu*IT_0329;
    const complex_t IT_0331 = 0.5*IT_0326 + 1.5*IT_0328 + 0.5*IT_0330;
    const complex_t IT_0332 = IT_0325 + IT_0331;
    const complex_t IT_0333 = IT_0321*IT_0332;
    const complex_t IT_0334 = IT_0094*IT_0117*IT_0188;
    const complex_t IT_0335 = 0.101321183642338*IT_0334;
    const complex_t IT_0336 = m_mu*IT_0121;
    const complex_t IT_0337 = mty::lt::C0iC(12, 2*IT_0025 + (-2)*s_12,
       IT_0025, IT_0025, IT_0027, IT_0027, IT_0050, mty::lt::reg_int);
    const complex_t IT_0338 = m_mu*IT_0337;
    const complex_t IT_0339 = IT_0336 + IT_0338;
    const complex_t IT_0340 = m_mu*IT_0122;
    const complex_t IT_0341 = mty::lt::C0iC(15, 2*IT_0025 + (-2)*s_12,
       IT_0025, IT_0025, IT_0027, IT_0027, IT_0050, mty::lt::reg_int);
    const complex_t IT_0342 = m_mu*IT_0341;
    const complex_t IT_0343 = mty::lt::C0iC(18, 2*IT_0025 + (-2)*s_12,
       IT_0025, IT_0025, IT_0027, IT_0027, IT_0050, mty::lt::reg_int);
    const complex_t IT_0344 = m_mu*IT_0343;
    const complex_t IT_0345 = 0.5*IT_0340 + 1.5*IT_0342 + 0.5*IT_0344;
    const complex_t IT_0346 = IT_0339 + IT_0345;
    const complex_t IT_0347 = IT_0335*IT_0346;
    const complex_t IT_0348 = IT_0016*IT_0097*IT_0129;
    const complex_t IT_0349 = 0.101321183642338*IT_0348;
    const complex_t IT_0350 = IT_0234*IT_0349;
    const complex_t IT_0351 = IT_0034*IT_0107*IT_0129;
    const complex_t IT_0352 = 0.101321183642338*IT_0351;
    const complex_t IT_0353 = IT_0248*IT_0352;
    const complex_t IT_0354 = IT_0045*IT_0117*IT_0129;
    const complex_t IT_0355 = 0.101321183642338*IT_0354;
    const complex_t IT_0356 = IT_0262*IT_0355;
    const complex_t IT_0357 = IT_0149*IT_0150*IT_0180;
    const complex_t IT_0358 = 0.101321183642338*IT_0357;
    const complex_t IT_0359 = m_mu*IT_0154;
    const complex_t IT_0360 = mty::lt::C0iC(12, 2*IT_0025 + (-2)*s_12,
       IT_0025, IT_0025, IT_0027, IT_0026, IT_0028, mty::lt::reg_int);
    const complex_t IT_0361 = m_mu*IT_0360;
    const complex_t IT_0362 = IT_0359 + IT_0361;
    const complex_t IT_0363 = m_mu*IT_0155;
    const complex_t IT_0364 = mty::lt::C0iC(15, 2*IT_0025 + (-2)*s_12,
       IT_0025, IT_0025, IT_0027, IT_0026, IT_0028, mty::lt::reg_int);
    const complex_t IT_0365 = m_mu*IT_0364;
    const complex_t IT_0366 = mty::lt::C0iC(18, 2*IT_0025 + (-2)*s_12,
       IT_0025, IT_0025, IT_0027, IT_0026, IT_0028, mty::lt::reg_int);
    const complex_t IT_0367 = m_mu*IT_0366;
    const complex_t IT_0368 = 0.5*IT_0363 + 1.5*IT_0365 + 0.5*IT_0367;
    const complex_t IT_0369 = IT_0362 + IT_0368;
    const complex_t IT_0370 = IT_0358*IT_0369;
    const complex_t IT_0371 = IT_0149*IT_0158*IT_0184;
    const complex_t IT_0372 = 0.101321183642338*IT_0371;
    const complex_t IT_0373 = m_mu*IT_0162;
    const complex_t IT_0374 = mty::lt::C0iC(12, 2*IT_0025 + (-2)*s_12,
       IT_0025, IT_0025, IT_0027, IT_0026, IT_0039, mty::lt::reg_int);
    const complex_t IT_0375 = m_mu*IT_0374;
    const complex_t IT_0376 = IT_0373 + IT_0375;
    const complex_t IT_0377 = m_mu*IT_0163;
    const complex_t IT_0378 = mty::lt::C0iC(15, 2*IT_0025 + (-2)*s_12,
       IT_0025, IT_0025, IT_0027, IT_0026, IT_0039, mty::lt::reg_int);
    const complex_t IT_0379 = m_mu*IT_0378;
    const complex_t IT_0380 = mty::lt::C0iC(18, 2*IT_0025 + (-2)*s_12,
       IT_0025, IT_0025, IT_0027, IT_0026, IT_0039, mty::lt::reg_int);
    const complex_t IT_0381 = m_mu*IT_0380;
    const complex_t IT_0382 = 0.5*IT_0377 + 1.5*IT_0379 + 0.5*IT_0381;
    const complex_t IT_0383 = IT_0376 + IT_0382;
    const complex_t IT_0384 = IT_0372*IT_0383;
    const complex_t IT_0385 = IT_0149*IT_0166*IT_0188;
    const complex_t IT_0386 = 0.101321183642338*IT_0385;
    const complex_t IT_0387 = m_mu*IT_0170;
    const complex_t IT_0388 = mty::lt::C0iC(12, 2*IT_0025 + (-2)*s_12,
       IT_0025, IT_0025, IT_0027, IT_0026, IT_0050, mty::lt::reg_int);
    const complex_t IT_0389 = m_mu*IT_0388;
    const complex_t IT_0390 = IT_0387 + IT_0389;
    const complex_t IT_0391 = m_mu*IT_0171;
    const complex_t IT_0392 = mty::lt::C0iC(15, 2*IT_0025 + (-2)*s_12,
       IT_0025, IT_0025, IT_0027, IT_0026, IT_0050, mty::lt::reg_int);
    const complex_t IT_0393 = m_mu*IT_0392;
    const complex_t IT_0394 = mty::lt::C0iC(18, 2*IT_0025 + (-2)*s_12,
       IT_0025, IT_0025, IT_0027, IT_0026, IT_0050, mty::lt::reg_int);
    const complex_t IT_0395 = m_mu*IT_0394;
    const complex_t IT_0396 = 0.5*IT_0391 + 1.5*IT_0393 + 0.5*IT_0395;
    const complex_t IT_0397 = IT_0390 + IT_0396;
    const complex_t IT_0398 = IT_0386*IT_0397;
    const complex_t IT_0399 = IT_0021*IT_0096*IT_0179;
    const complex_t IT_0400 = 0.101321183642338*IT_0399;
    const complex_t IT_0401 = IT_0318*IT_0400;
    const complex_t IT_0402 = IT_0036*IT_0106*IT_0179;
    const complex_t IT_0403 = 0.101321183642338*IT_0402;
    const complex_t IT_0404 = IT_0332*IT_0403;
    const complex_t IT_0405 = IT_0047*IT_0116*IT_0179;
    const complex_t IT_0406 = 0.101321183642338*IT_0405;
    const complex_t IT_0407 = IT_0346*IT_0406;
    const complex_t IT_0408 = IT_0063*IT_0096*IT_0197;
    const complex_t IT_0409 = 0.101321183642338*IT_0408;
    const complex_t IT_0410 = IT_0369*IT_0409;
    const complex_t IT_0411 = IT_0072*IT_0106*IT_0197;
    const complex_t IT_0412 = 0.101321183642338*IT_0411;
    const complex_t IT_0413 = IT_0383*IT_0412;
    const complex_t IT_0414 = IT_0081*IT_0116*IT_0197;
    const complex_t IT_0415 = 0.101321183642338*IT_0414;
    const complex_t IT_0416 = IT_0397*IT_0415;
    const complex_t IT_0417 = IT_0016*IT_0150*IT_0211;
    const complex_t IT_0418 = 0.101321183642338*IT_0417;
    const complex_t IT_0419 = IT_0276*IT_0418;
    const complex_t IT_0420 = IT_0034*IT_0158*IT_0211;
    const complex_t IT_0421 = 0.101321183642338*IT_0420;
    const complex_t IT_0422 = IT_0290*IT_0421;
    const complex_t IT_0423 = IT_0045*IT_0166*IT_0211;
    const complex_t IT_0424 = 0.101321183642338*IT_0423;
    const complex_t IT_0425 = IT_0304*IT_0424;
    const complex_t IT_0426 = 2*IT_0235 + 2*IT_0249 + 2*IT_0263 + 2*IT_0277 +
       2*IT_0291 + 2*IT_0305 + 2*IT_0319 + 2*IT_0333 + 2*IT_0347 + 2*IT_0350 + 2
      *IT_0353 + 2*IT_0356 + 2*IT_0370 + 2*IT_0384 + 2*IT_0398 + 2*IT_0401 + 2
      *IT_0404 + 2*IT_0407 + 2*IT_0410 + 2*IT_0413 + 2*IT_0416 + 2*IT_0419 + 2
      *IT_0422 + 2*IT_0425;
    const complex_t IT_0427 = IT_0221 + IT_0426;
    const complex_t IT_0428 = IT_0001*IT_0427;
    const complex_t IT_0429 = 0.25*IT_0428;
    const complex_t IT_0430 = IT_0013*IT_0097*IT_0131;
    const complex_t IT_0431 = IT_0089*IT_0430;
    const complex_t IT_0432 = IT_0029*IT_0431;
    const complex_t IT_0433 = IT_0013*IT_0107*IT_0136;
    const complex_t IT_0434 = IT_0089*IT_0433;
    const complex_t IT_0435 = IT_0040*IT_0434;
    const complex_t IT_0436 = IT_0013*IT_0117*IT_0141;
    const complex_t IT_0437 = IT_0089*IT_0436;
    const complex_t IT_0438 = IT_0051*IT_0437;
    const complex_t IT_0439 = IT_0061*IT_0131*IT_0150;
    const complex_t IT_0440 = IT_0023*IT_0439;
    const complex_t IT_0441 = IT_0066*IT_0440;
    const complex_t IT_0442 = IT_0061*IT_0136*IT_0158;
    const complex_t IT_0443 = IT_0023*IT_0442;
    const complex_t IT_0444 = IT_0075*IT_0443;
    const complex_t IT_0445 = IT_0061*IT_0141*IT_0166;
    const complex_t IT_0446 = IT_0023*IT_0445;
    const complex_t IT_0447 = IT_0084*IT_0446;
    const complex_t IT_0448 = IT_0021*IT_0094*IT_0180;
    const complex_t IT_0449 = IT_0089*IT_0448;
    const complex_t IT_0450 = IT_0101*IT_0449;
    const complex_t IT_0451 = IT_0036*IT_0094*IT_0184;
    const complex_t IT_0452 = IT_0089*IT_0451;
    const complex_t IT_0453 = IT_0111*IT_0452;
    const complex_t IT_0454 = IT_0047*IT_0094*IT_0188;
    const complex_t IT_0455 = IT_0089*IT_0454;
    const complex_t IT_0456 = IT_0121*IT_0455;
    const complex_t IT_0457 = IT_0016*IT_0021*IT_0129;
    const complex_t IT_0458 = IT_0089*IT_0457;
    const complex_t IT_0459 = IT_0029*IT_0458;
    const complex_t IT_0460 = IT_0034*IT_0036*IT_0129;
    const complex_t IT_0461 = IT_0089*IT_0460;
    const complex_t IT_0462 = IT_0040*IT_0461;
    const complex_t IT_0463 = IT_0045*IT_0047*IT_0129;
    const complex_t IT_0464 = IT_0089*IT_0463;
    const complex_t IT_0465 = IT_0051*IT_0464;
    const complex_t IT_0466 = IT_0063*IT_0149*IT_0180;
    const complex_t IT_0467 = IT_0023*IT_0466;
    const complex_t IT_0468 = IT_0154*IT_0467;
    const complex_t IT_0469 = IT_0072*IT_0149*IT_0184;
    const complex_t IT_0470 = IT_0023*IT_0469;
    const complex_t IT_0471 = IT_0162*IT_0470;
    const complex_t IT_0472 = IT_0081*IT_0149*IT_0188;
    const complex_t IT_0473 = IT_0023*IT_0472;
    const complex_t IT_0474 = IT_0170*IT_0473;
    const complex_t IT_0475 = IT_0096*IT_0097*IT_0179;
    const complex_t IT_0476 = IT_0089*IT_0475;
    const complex_t IT_0477 = IT_0101*IT_0476;
    const complex_t IT_0478 = IT_0106*IT_0107*IT_0179;
    const complex_t IT_0479 = IT_0089*IT_0478;
    const complex_t IT_0480 = IT_0111*IT_0479;
    const complex_t IT_0481 = IT_0116*IT_0117*IT_0179;
    const complex_t IT_0482 = IT_0089*IT_0481;
    const complex_t IT_0483 = IT_0121*IT_0482;
    const complex_t IT_0484 = IT_0096*IT_0150*IT_0197;
    const complex_t IT_0485 = IT_0023*IT_0484;
    const complex_t IT_0486 = IT_0154*IT_0485;
    const complex_t IT_0487 = IT_0106*IT_0158*IT_0197;
    const complex_t IT_0488 = IT_0023*IT_0487;
    const complex_t IT_0489 = IT_0162*IT_0488;
    const complex_t IT_0490 = IT_0116*IT_0166*IT_0197;
    const complex_t IT_0491 = IT_0023*IT_0490;
    const complex_t IT_0492 = IT_0170*IT_0491;
    const complex_t IT_0493 = IT_0016*IT_0063*IT_0211;
    const complex_t IT_0494 = IT_0023*IT_0493;
    const complex_t IT_0495 = IT_0066*IT_0494;
    const complex_t IT_0496 = IT_0034*IT_0072*IT_0211;
    const complex_t IT_0497 = IT_0023*IT_0496;
    const complex_t IT_0498 = IT_0075*IT_0497;
    const complex_t IT_0499 = IT_0045*IT_0081*IT_0211;
    const complex_t IT_0500 = IT_0023*IT_0499;
    const complex_t IT_0501 = IT_0084*IT_0500;
    const complex_t IT_0502 = IT_0432 + IT_0435 + IT_0438 + IT_0441 + IT_0444 
      + IT_0447 + IT_0450 + IT_0453 + IT_0456 + IT_0459 + IT_0462 + IT_0465 +
       IT_0468 + IT_0471 + IT_0474 + IT_0477 + IT_0480 + IT_0483 + IT_0486 +
       IT_0489 + IT_0492 + IT_0495 + IT_0498 + IT_0501;
    const complex_t IT_0503 = 0.5*IT_0230;
    const complex_t IT_0504 = IT_0227 + IT_0503;
    const complex_t IT_0505 = IT_0223*IT_0504;
    const complex_t IT_0506 = 0.5*IT_0244;
    const complex_t IT_0507 = IT_0241 + IT_0506;
    const complex_t IT_0508 = IT_0237*IT_0507;
    const complex_t IT_0509 = 0.5*IT_0258;
    const complex_t IT_0510 = IT_0255 + IT_0509;
    const complex_t IT_0511 = IT_0251*IT_0510;
    const complex_t IT_0512 = 0.5*IT_0272;
    const complex_t IT_0513 = IT_0269 + IT_0512;
    const complex_t IT_0514 = IT_0265*IT_0513;
    const complex_t IT_0515 = 0.5*IT_0286;
    const complex_t IT_0516 = IT_0283 + IT_0515;
    const complex_t IT_0517 = IT_0279*IT_0516;
    const complex_t IT_0518 = 0.5*IT_0300;
    const complex_t IT_0519 = IT_0297 + IT_0518;
    const complex_t IT_0520 = IT_0293*IT_0519;
    const complex_t IT_0521 = 0.5*IT_0314;
    const complex_t IT_0522 = IT_0311 + IT_0521;
    const complex_t IT_0523 = IT_0307*IT_0522;
    const complex_t IT_0524 = 0.5*IT_0328;
    const complex_t IT_0525 = IT_0325 + IT_0524;
    const complex_t IT_0526 = IT_0321*IT_0525;
    const complex_t IT_0527 = 0.5*IT_0342;
    const complex_t IT_0528 = IT_0339 + IT_0527;
    const complex_t IT_0529 = IT_0335*IT_0528;
    const complex_t IT_0530 = IT_0349*IT_0504;
    const complex_t IT_0531 = IT_0352*IT_0507;
    const complex_t IT_0532 = IT_0355*IT_0510;
    const complex_t IT_0533 = 0.5*IT_0365;
    const complex_t IT_0534 = IT_0362 + IT_0533;
    const complex_t IT_0535 = IT_0358*IT_0534;
    const complex_t IT_0536 = 0.5*IT_0379;
    const complex_t IT_0537 = IT_0376 + IT_0536;
    const complex_t IT_0538 = IT_0372*IT_0537;
    const complex_t IT_0539 = 0.5*IT_0393;
    const complex_t IT_0540 = IT_0390 + IT_0539;
    const complex_t IT_0541 = IT_0386*IT_0540;
    const complex_t IT_0542 = IT_0400*IT_0522;
    const complex_t IT_0543 = IT_0403*IT_0525;
    const complex_t IT_0544 = IT_0406*IT_0528;
    const complex_t IT_0545 = IT_0409*IT_0534;
    const complex_t IT_0546 = IT_0412*IT_0537;
    const complex_t IT_0547 = IT_0415*IT_0540;
    const complex_t IT_0548 = IT_0418*IT_0513;
    const complex_t IT_0549 = IT_0421*IT_0516;
    const complex_t IT_0550 = IT_0424*IT_0519;
    const complex_t IT_0551 = (-2)*IT_0505 + (-2)*IT_0508 + (-2)*IT_0511 + (-2
      )*IT_0514 + (-2)*IT_0517 + (-2)*IT_0520 + (-2)*IT_0523 + (-2)*IT_0526 + (
      -2)*IT_0529 + (-2)*IT_0530 + (-2)*IT_0531 + (-2)*IT_0532 + (-2)*IT_0535 + 
      (-2)*IT_0538 + (-2)*IT_0541 + (-2)*IT_0542 + (-2)*IT_0543 + (-2)*IT_0544 +
       (-2)*IT_0545 + (-2)*IT_0546 + (-2)*IT_0547 + (-2)*IT_0548 + (-2)*IT_0549 
      + (-2)*IT_0550;
    const complex_t IT_0552 = IT_0502 + IT_0551;
    const complex_t IT_0553 = IT_0001*IT_0552;
    const complex_t IT_0554 = 0.25*IT_0553;
    const complex_t IT_0555 = IT_0228 + IT_0230 + IT_0232;
    const complex_t IT_0556 = IT_0223*IT_0555;
    const complex_t IT_0557 = IT_0242 + IT_0244 + IT_0246;
    const complex_t IT_0558 = IT_0237*IT_0557;
    const complex_t IT_0559 = IT_0256 + IT_0258 + IT_0260;
    const complex_t IT_0560 = IT_0251*IT_0559;
    const complex_t IT_0561 = IT_0270 + IT_0272 + IT_0274;
    const complex_t IT_0562 = IT_0265*IT_0561;
    const complex_t IT_0563 = IT_0284 + IT_0286 + IT_0288;
    const complex_t IT_0564 = IT_0279*IT_0563;
    const complex_t IT_0565 = IT_0298 + IT_0300 + IT_0302;
    const complex_t IT_0566 = IT_0293*IT_0565;
    const complex_t IT_0567 = IT_0312 + IT_0314 + IT_0316;
    const complex_t IT_0568 = IT_0400*IT_0567;
    const complex_t IT_0569 = IT_0326 + IT_0328 + IT_0330;
    const complex_t IT_0570 = IT_0403*IT_0569;
    const complex_t IT_0571 = IT_0340 + IT_0342 + IT_0344;
    const complex_t IT_0572 = IT_0406*IT_0571;
    const complex_t IT_0573 = IT_0363 + IT_0365 + IT_0367;
    const complex_t IT_0574 = IT_0409*IT_0573;
    const complex_t IT_0575 = IT_0377 + IT_0379 + IT_0381;
    const complex_t IT_0576 = IT_0412*IT_0575;
    const complex_t IT_0577 = IT_0391 + IT_0393 + IT_0395;
    const complex_t IT_0578 = IT_0415*IT_0577;
    const complex_t IT_0579 = IT_0033 + IT_0044 + IT_0055 + IT_0070 + IT_0079 
      + IT_0088 + IT_0183 + IT_0187 + IT_0191 + IT_0200 + IT_0203 + IT_0206 +
       IT_0556 + IT_0558 + IT_0560 + IT_0562 + IT_0564 + IT_0566 + IT_0568 +
       IT_0570 + IT_0572 + IT_0574 + IT_0576 + IT_0578;
    const complex_t IT_0580 = IT_0307*IT_0567;
    const complex_t IT_0581 = IT_0321*IT_0569;
    const complex_t IT_0582 = IT_0335*IT_0571;
    const complex_t IT_0583 = IT_0349*IT_0555;
    const complex_t IT_0584 = IT_0352*IT_0557;
    const complex_t IT_0585 = IT_0355*IT_0559;
    const complex_t IT_0586 = IT_0358*IT_0573;
    const complex_t IT_0587 = IT_0372*IT_0575;
    const complex_t IT_0588 = IT_0386*IT_0577;
    const complex_t IT_0589 = IT_0418*IT_0561;
    const complex_t IT_0590 = IT_0421*IT_0563;
    const complex_t IT_0591 = IT_0424*IT_0565;
    const complex_t IT_0592 = -IT_0104 + -IT_0114 + -IT_0124 + -IT_0134 + 
      -IT_0139 + -IT_0144 + -IT_0157 + -IT_0165 + -IT_0173 + -IT_0214 + -IT_0217
       + -IT_0220 + -IT_0580 + -IT_0581 + -IT_0582 + -IT_0583 + -IT_0584 + 
      -IT_0585 + -IT_0586 + -IT_0587 + -IT_0588 + -IT_0589 + -IT_0590 + -IT_0591;
    const complex_t IT_0593 = IT_0579 + IT_0592;
    const complex_t IT_0594 = IT_0001*IT_0593;
    const complex_t IT_0595 = (-0.25)*IT_0594;
    const complex_t IT_0596 = IT_0223*IT_0230;
    const complex_t IT_0597 = IT_0237*IT_0244;
    const complex_t IT_0598 = IT_0251*IT_0258;
    const complex_t IT_0599 = IT_0265*IT_0272;
    const complex_t IT_0600 = IT_0279*IT_0286;
    const complex_t IT_0601 = IT_0293*IT_0300;
    const complex_t IT_0602 = IT_0314*IT_0400;
    const complex_t IT_0603 = IT_0328*IT_0403;
    const complex_t IT_0604 = IT_0342*IT_0406;
    const complex_t IT_0605 = IT_0365*IT_0409;
    const complex_t IT_0606 = IT_0379*IT_0412;
    const complex_t IT_0607 = IT_0393*IT_0415;
    const complex_t IT_0608 = IT_0432 + IT_0435 + IT_0438 + IT_0441 + IT_0444 
      + IT_0447 + IT_0477 + IT_0480 + IT_0483 + IT_0486 + IT_0489 + IT_0492 +
       IT_0596 + IT_0597 + IT_0598 + IT_0599 + IT_0600 + IT_0601 + IT_0602 +
       IT_0603 + IT_0604 + IT_0605 + IT_0606 + IT_0607;
    const complex_t IT_0609 = IT_0307*IT_0314;
    const complex_t IT_0610 = IT_0321*IT_0328;
    const complex_t IT_0611 = IT_0335*IT_0342;
    const complex_t IT_0612 = IT_0230*IT_0349;
    const complex_t IT_0613 = IT_0244*IT_0352;
    const complex_t IT_0614 = IT_0258*IT_0355;
    const complex_t IT_0615 = IT_0358*IT_0365;
    const complex_t IT_0616 = IT_0372*IT_0379;
    const complex_t IT_0617 = IT_0386*IT_0393;
    const complex_t IT_0618 = IT_0272*IT_0418;
    const complex_t IT_0619 = IT_0286*IT_0421;
    const complex_t IT_0620 = IT_0300*IT_0424;
    const complex_t IT_0621 = -IT_0450 + -IT_0453 + -IT_0456 + -IT_0459 + 
      -IT_0462 + -IT_0465 + -IT_0468 + -IT_0471 + -IT_0474 + -IT_0495 + -IT_0498
       + -IT_0501 + -IT_0609 + -IT_0610 + -IT_0611 + -IT_0612 + -IT_0613 + 
      -IT_0614 + -IT_0615 + -IT_0616 + -IT_0617 + -IT_0618 + -IT_0619 + -IT_0620;
    const complex_t IT_0622 = IT_0608 + IT_0621;
    const complex_t IT_0623 = IT_0001*IT_0622;
    const complex_t IT_0624 = 0.25*IT_0623;
    return (complex_t{0, 0.25})*IT_0429 + (complex_t{0, 0.25})*IT_0554 + 
      (complex_t{0, 0.25})*IT_0595 + (complex_t{0, 0.25})*IT_0624;
}
} // End of namespace c9_nmfv
