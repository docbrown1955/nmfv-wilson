#include "global.h"
#include "libdiagonalization.h"
#include "c9_nmfv.h"
#include "libcomplexop.h"

namespace c9_nmfv {


void updateSpectrum(param_t &params)
{
    updateDiagonalization(params);
    updateMassExpressions(params);
}

void updateDiagonalization(param_t &params)
{
    SpectrumInput  inputs;
    readDiagonalizationInputs (inputs,  params);
    SpectrumOutput outputs = updateDiagonalization(inputs);
    readDiagonalizationOutputs(outputs, params);
}

SpectrumOutput updateDiagonalization(SpectrumInput const &inputs)
{
    auto const &A_b = inputs.A_b;
    auto const &A_t = inputs.A_t;
    auto const &A_tau = inputs.A_tau;
    auto const &M_1 = inputs.M_1;
    auto const &M_2 = inputs.M_2;
    auto const &M_W = inputs.M_W;
    auto const &M_eL = inputs.M_eL;
    auto const &M_eR = inputs.M_eR;
    auto const &M_q1L = inputs.M_q1L;
    auto const &M_q3L = inputs.M_q3L;
    auto const &M_qbR = inputs.M_qbR;
    auto const &M_qdR = inputs.M_qdR;
    auto const &M_qtR = inputs.M_qtR;
    auto const &M_quR = inputs.M_quR;
    auto const &M_tauL = inputs.M_tauL;
    auto const &M_tauR = inputs.M_tauR;
    auto const &V_cb = inputs.V_cb;
    auto const &V_cd = inputs.V_cd;
    auto const &V_cs = inputs.V_cs;
    auto const &V_tb = inputs.V_tb;
    auto const &V_td = inputs.V_td;
    auto const &V_ts = inputs.V_ts;
    auto const &V_ub_mod = inputs.V_ub_mod;
    auto const &V_ud = inputs.V_ud;
    auto const &V_us = inputs.V_us;
    auto const &beta = inputs.beta;
    auto const &del_DR_1 = inputs.del_DR_1;
    auto const &del_DR_12 = inputs.del_DR_12;
    auto const &del_DR_2 = inputs.del_DR_2;
    auto const &del_ER_1 = inputs.del_ER_1;
    auto const &del_ER_12 = inputs.del_ER_12;
    auto const &del_ER_2 = inputs.del_ER_2;
    auto const &del_LL_1 = inputs.del_LL_1;
    auto const &del_LL_12 = inputs.del_LL_12;
    auto const &del_LL_2 = inputs.del_LL_2;
    auto const &del_QL_1 = inputs.del_QL_1;
    auto const &del_QL_12 = inputs.del_QL_12;
    auto const &del_QL_2 = inputs.del_QL_2;
    auto const &del_UR_1 = inputs.del_UR_1;
    auto const &del_UR_12 = inputs.del_UR_12;
    auto const &del_UR_2 = inputs.del_UR_2;
    auto const &delta_wolf = inputs.delta_wolf;
    auto const &e_em = inputs.e_em;
    auto const &m_b = inputs.m_b;
    auto const &m_c = inputs.m_c;
    auto const &m_d = inputs.m_d;
    auto const &m_e = inputs.m_e;
    auto const &m_mu = inputs.m_mu;
    auto const &m_s = inputs.m_s;
    auto const &m_t = inputs.m_t;
    auto const &m_tau = inputs.m_tau;
    auto const &m_u = inputs.m_u;
    auto const &mu_h = inputs.mu_h;
    auto const &theta_W = inputs.theta_W;

    SpectrumOutput outputs;

    Diagonalizer::applyDiagonalization(
        {
            M_1,
            0,
            (-2)*M_W*ccosq(beta)*csinq(theta_W)/ccosq(theta_W),
            2*M_W*csinq(beta)*csinq(theta_W)/ccosq(theta_W),
            0,
            M_2,
            2*M_W*ccosq(beta),
            (-2)*M_W*csinq(beta),
            0,
            0,
            0,
            (-2)*mu_h,
            0,
            0,
            0,
            0,
        },
        {&outputs.N_B1, &outputs.N_B2, &outputs.N_B3, &outputs.N_B4, &outputs.N_W1, &outputs.N_W2, &outputs.N_W3, &outputs.N_W4, &outputs.N_d1, &outputs.N_d2, &outputs.N_d3, &outputs.N_d4, &outputs.N_u1, &outputs.N_u2, &outputs.N_u3, &outputs.N_u4, },
        {&outputs.m_N_1, &outputs.m_N_2, &outputs.m_N_3, &outputs.m_N_4, }
        );

    Diagonalizer::applyBiDiagonalization(
        {
            M_2,
            1.41421*M_W*csinq(beta),
            1.41421*M_W*ccosq(beta),
            mu_h,
        },
        {&outputs.V_Wp1, &outputs.V_Wp2, &outputs.V_u1, &outputs.V_u2, },
        {&outputs.U_Wm1, &outputs.U_Wm2, &outputs.U_d1, &outputs.U_d2, },
        {&outputs.m_C1p, &outputs.m_C2p, }
        );

    Diagonalizer::applyDiagonalization(
        {
            cpowq(m_u, 2) + cpowq(M_q1L, 2) + 0.5*cpowq(M_W, 2)*cpowq(ccosq(beta), 2) + (-0.5)*cpowq(M_W, 2)*cpowq(csinq(beta), 2) + (-0.166667)*cpowq(M_W, 2)*cpowq(ccosq(beta), 2)*cpowq(ccosq(theta_W), -2)*cpowq(csinq(theta_W), 2) + 0.166667*cpowq(M_W, 2)*cpowq(ccosq(theta_W), -2)*cpowq(csinq(beta), 2)*cpowq(csinq(theta_W), 2),
            del_QL_1,
            del_QL_2,
            (-2)*m_u*mu_h*ccosq(beta)/csinq(beta),
            0,
            0,
            0,
            cpowq(m_c, 2) + cpowq(M_q1L, 2) + 0.5*cpowq(M_W, 2)*cpowq(ccosq(beta), 2) + (-0.5)*cpowq(M_W, 2)*cpowq(csinq(beta), 2) + (-0.166667)*cpowq(M_W, 2)*cpowq(ccosq(beta), 2)*cpowq(ccosq(theta_W), -2)*cpowq(csinq(theta_W), 2) + 0.166667*cpowq(M_W, 2)*cpowq(ccosq(theta_W), -2)*cpowq(csinq(beta), 2)*cpowq(csinq(theta_W), 2),
            del_QL_12,
            0,
            (-2)*m_c*mu_h*ccosq(beta)/csinq(beta),
            0,
            0,
            0,
            cpowq(m_t, 2) + cpowq(M_q3L, 2) + 0.5*cpowq(M_W, 2)*cpowq(ccosq(beta), 2) + (-0.5)*cpowq(M_W, 2)*cpowq(csinq(beta), 2) + (-0.166667)*cpowq(M_W, 2)*cpowq(ccosq(beta), 2)*cpowq(ccosq(theta_W), -2)*cpowq(csinq(theta_W), 2) + 0.166667*cpowq(M_W, 2)*cpowq(ccosq(theta_W), -2)*cpowq(csinq(beta), 2)*cpowq(csinq(theta_W), 2),
            0,
            0,
            (-2)*m_t*mu_h*ccosq(beta)/csinq(beta) + 2.82843*A_t*M_W*csinq(beta)*csinq(theta_W)/e_em,
            0,
            0,
            0,
            cpowq(m_u, 2) + cpowq(M_quR, 2) + 0.666667*cpowq(M_W, 2)*cpowq(ccosq(beta), 2)*cpowq(ccosq(theta_W), -2)*cpowq(csinq(theta_W), 2) + (-0.666667)*cpowq(M_W, 2)*cpowq(ccosq(theta_W), -2)*cpowq(csinq(beta), 2)*cpowq(csinq(theta_W), 2),
            del_UR_1,
            del_UR_2,
            0,
            0,
            0,
            0,
            cpowq(m_c, 2) + cpowq(M_quR, 2) + 0.666667*cpowq(M_W, 2)*cpowq(ccosq(beta), 2)*cpowq(ccosq(theta_W), -2)*cpowq(csinq(theta_W), 2) + (-0.666667)*cpowq(M_W, 2)*cpowq(ccosq(theta_W), -2)*cpowq(csinq(beta), 2)*cpowq(csinq(theta_W), 2),
            del_UR_12,
            0,
            0,
            0,
            0,
            0,
            cpowq(m_t, 2) + cpowq(M_qtR, 2) + 0.666667*cpowq(M_W, 2)*cpowq(ccosq(beta), 2)*cpowq(ccosq(theta_W), -2)*cpowq(csinq(theta_W), 2) + (-0.666667)*cpowq(M_W, 2)*cpowq(ccosq(theta_W), -2)*cpowq(csinq(beta), 2)*cpowq(csinq(theta_W), 2),
        },
        {&outputs.U_su_00, &outputs.U_su_01, &outputs.U_su_02, &outputs.U_su_03, &outputs.U_su_04, &outputs.U_su_05, &outputs.U_su_10, &outputs.U_su_11, &outputs.U_su_12, &outputs.U_su_13, &outputs.U_su_14, &outputs.U_su_15, &outputs.U_su_20, &outputs.U_su_21, &outputs.U_su_22, &outputs.U_su_23, &outputs.U_su_24, &outputs.U_su_25, &outputs.U_su_30, &outputs.U_su_31, &outputs.U_su_32, &outputs.U_su_33, &outputs.U_su_34, &outputs.U_su_35, &outputs.U_su_40, &outputs.U_su_41, &outputs.U_su_42, &outputs.U_su_43, &outputs.U_su_44, &outputs.U_su_45, &outputs.U_su_50, &outputs.U_su_51, &outputs.U_su_52, &outputs.U_su_53, &outputs.U_su_54, &outputs.U_su_55, },
        {&outputs.m_su_L, &outputs.m_sc_L, &outputs.m_st_L, &outputs.m_su_R, &outputs.m_sc_R, &outputs.m_st_R, }
        );
    if (0 > outputs.m_su_L) {
        std::cerr << "Warning: negative squared mass for " << "m_su_L" << ".\n";
    }
    outputs.m_su_L = sqrtq(outputs.m_su_L);
    if (0 > outputs.m_sc_L) {
        std::cerr << "Warning: negative squared mass for " << "m_sc_L" << ".\n";
    }
    outputs.m_sc_L = sqrtq(outputs.m_sc_L);
    if (0 > outputs.m_st_L) {
        std::cerr << "Warning: negative squared mass for " << "m_st_L" << ".\n";
    }
    outputs.m_st_L = sqrtq(outputs.m_st_L);
    if (0 > outputs.m_su_R) {
        std::cerr << "Warning: negative squared mass for " << "m_su_R" << ".\n";
    }
    outputs.m_su_R = sqrtq(outputs.m_su_R);
    if (0 > outputs.m_sc_R) {
        std::cerr << "Warning: negative squared mass for " << "m_sc_R" << ".\n";
    }
    outputs.m_sc_R = sqrtq(outputs.m_sc_R);
    if (0 > outputs.m_st_R) {
        std::cerr << "Warning: negative squared mass for " << "m_st_R" << ".\n";
    }
    outputs.m_st_R = sqrtq(outputs.m_st_R);

    Diagonalizer::applyDiagonalization(
        {
            cpowq(m_d, 2) + V_cd*conjq(V_cd)*cpowq(M_q1L, 2) + cpowq(V_ud, 2)*cpowq(M_q1L, 2) + V_td*conjq(V_td)*cpowq(M_q3L, 2) + 0.5*V_cd*V_ud*del_QL_1 + 0.5*conjq(V_cd)*V_ud*del_QL_1 + 0.5*V_td*V_ud*del_QL_2 + 0.5*conjq(V_td)*V_ud*del_QL_2 + 0.5*conjq(V_cd)*V_td*del_QL_12 + 0.5*V_cd*conjq(V_td)*del_QL_12 + (-0.5)*cpowq(M_W, 2)*cpowq(ccosq(beta), 2) + 0.5*cpowq(M_W, 2)*cpowq(csinq(beta), 2) + (-0.166667)*cpowq(M_W, 2)*cpowq(ccosq(beta), 2)*cpowq(ccosq(theta_W), -2)*cpowq(csinq(theta_W), 2) + 0.166667*cpowq(M_W, 2)*cpowq(ccosq(theta_W), -2)*cpowq(csinq(beta), 2)*cpowq(csinq(theta_W), 2),
            conjq(V_cd)*V_cs*cpowq(M_q1L, 2) + V_cd*conjq(V_cs)*cpowq(M_q1L, 2) + 2*V_ud*V_us*cpowq(M_q1L, 2) + conjq(V_td)*V_ts*cpowq(M_q3L, 2) + V_td*conjq(V_ts)*cpowq(M_q3L, 2) + 0.5*V_cs*V_ud*del_QL_1 + 0.5*conjq(V_cs)*V_ud*del_QL_1 + 0.5*V_cd*V_us*del_QL_1 + 0.5*conjq(V_cd)*V_us*del_QL_1 + 0.5*V_ts*V_ud*del_QL_2 + 0.5*conjq(V_ts)*V_ud*del_QL_2 + 0.5*V_td*V_us*del_QL_2 + 0.5*conjq(V_td)*V_us*del_QL_2 + 0.5*conjq(V_cs)*V_td*del_QL_12 + 0.5*V_cs*conjq(V_td)*del_QL_12 + 0.5*conjq(V_cd)*V_ts*del_QL_12 + 0.5*V_cd*conjq(V_ts)*del_QL_12,
            V_cb*V_cd*cpowq(M_q1L, 2) + V_cb*conjq(V_cd)*cpowq(M_q1L, 2) + V_tb*V_td*cpowq(M_q3L, 2) + V_tb*conjq(V_td)*cpowq(M_q3L, 2) + V_cb*V_ud*del_QL_1 + V_tb*V_ud*del_QL_2 + 0.5*V_cd*V_tb*del_QL_12 + 0.5*conjq(V_cd)*V_tb*del_QL_12 + 0.5*V_cb*V_td*del_QL_12 + 0.5*V_cb*conjq(V_td)*del_QL_12 + V_ud*cpowq(M_q1L, 2)*V_ub_mod*cexpq((complex_t{0, 1})*delta_wolf) + 0.5*V_cd*V_ub_mod*del_QL_1*cexpq((complex_t{0, 1})*delta_wolf) + 0.5*V_td*V_ub_mod*del_QL_2*cexpq((complex_t{0, 1})*delta_wolf) + V_ud*cpowq(M_q1L, 2)*V_ub_mod*cexpq((complex_t{0, -1})*delta_wolf) + 0.5*conjq(V_cd)*V_ub_mod*del_QL_1*cexpq((complex_t{0, -1})*delta_wolf) + 0.5*conjq(V_td)*V_ub_mod*del_QL_2*cexpq((complex_t{0, -1})*delta_wolf),
            (-2)*m_d*mu_h*csinq(beta)/ccosq(beta) + 2.82843*A_b*M_W*V_td*conjq(V_td)*ccosq(beta)*csinq(theta_W)/e_em,
            1.41421*A_b*M_W*conjq(V_td)*V_ts*ccosq(beta)*csinq(theta_W)/e_em + 1.41421*A_b*M_W*V_td*conjq(V_ts)*ccosq(beta)*csinq(theta_W)/e_em,
            1.41421*A_b*M_W*V_tb*V_td*ccosq(beta)*csinq(theta_W)/e_em + 1.41421*A_b*M_W*V_tb*conjq(V_td)*ccosq(beta)*csinq(theta_W)/e_em,
            0,
            cpowq(m_s, 2) + V_cs*conjq(V_cs)*cpowq(M_q1L, 2) + cpowq(V_us, 2)*cpowq(M_q1L, 2) + V_ts*conjq(V_ts)*cpowq(M_q3L, 2) + 0.5*V_cs*V_us*del_QL_1 + 0.5*conjq(V_cs)*V_us*del_QL_1 + 0.5*V_ts*V_us*del_QL_2 + 0.5*conjq(V_ts)*V_us*del_QL_2 + 0.5*conjq(V_cs)*V_ts*del_QL_12 + 0.5*V_cs*conjq(V_ts)*del_QL_12 + (-0.5)*cpowq(M_W, 2)*cpowq(ccosq(beta), 2) + 0.5*cpowq(M_W, 2)*cpowq(csinq(beta), 2) + (-0.166667)*cpowq(M_W, 2)*cpowq(ccosq(beta), 2)*cpowq(ccosq(theta_W), -2)*cpowq(csinq(theta_W), 2) + 0.166667*cpowq(M_W, 2)*cpowq(ccosq(theta_W), -2)*cpowq(csinq(beta), 2)*cpowq(csinq(theta_W), 2),
            V_cb*V_cs*cpowq(M_q1L, 2) + V_cb*conjq(V_cs)*cpowq(M_q1L, 2) + V_tb*V_ts*cpowq(M_q3L, 2) + V_tb*conjq(V_ts)*cpowq(M_q3L, 2) + V_cb*V_us*del_QL_1 + V_tb*V_us*del_QL_2 + 0.5*V_cs*V_tb*del_QL_12 + 0.5*conjq(V_cs)*V_tb*del_QL_12 + 0.5*V_cb*V_ts*del_QL_12 + 0.5*V_cb*conjq(V_ts)*del_QL_12 + V_us*cpowq(M_q1L, 2)*V_ub_mod*cexpq((complex_t{0, 1})*delta_wolf) + 0.5*V_cs*V_ub_mod*del_QL_1*cexpq((complex_t{0, 1})*delta_wolf) + 0.5*V_ts*V_ub_mod*del_QL_2*cexpq((complex_t{0, 1})*delta_wolf) + V_us*cpowq(M_q1L, 2)*V_ub_mod*cexpq((complex_t{0, -1})*delta_wolf) + 0.5*conjq(V_cs)*V_ub_mod*del_QL_1*cexpq((complex_t{0, -1})*delta_wolf) + 0.5*conjq(V_ts)*V_ub_mod*del_QL_2*cexpq((complex_t{0, -1})*delta_wolf),
            1.41421*A_b*M_W*conjq(V_td)*V_ts*ccosq(beta)*csinq(theta_W)/e_em + 1.41421*A_b*M_W*V_td*conjq(V_ts)*ccosq(beta)*csinq(theta_W)/e_em,
            (-2)*m_s*mu_h*csinq(beta)/ccosq(beta) + 2.82843*A_b*M_W*V_ts*conjq(V_ts)*ccosq(beta)*csinq(theta_W)/e_em,
            1.41421*A_b*M_W*V_tb*V_ts*ccosq(beta)*csinq(theta_W)/e_em + 1.41421*A_b*M_W*V_tb*conjq(V_ts)*ccosq(beta)*csinq(theta_W)/e_em,
            0,
            0,
            cpowq(m_b, 2) + cpowq(V_cb, 2)*cpowq(M_q1L, 2) + cpowq(V_tb, 2)*cpowq(M_q3L, 2) + V_cb*V_tb*del_QL_12 + 0.5*V_cb*V_ub_mod*del_QL_1*cexpq((complex_t{0, 1})*delta_wolf) + 0.5*V_tb*V_ub_mod*del_QL_2*cexpq((complex_t{0, 1})*delta_wolf) + 0.5*V_cb*V_ub_mod*del_QL_1*cexpq((complex_t{0, -1})*delta_wolf) + 0.5*V_tb*V_ub_mod*del_QL_2*cexpq((complex_t{0, -1})*delta_wolf) + cpowq(M_q1L, 2)*cpowq(V_ub_mod, 2)*cexpq((complex_t{0, 1})*delta_wolf)*cexpq((complex_t{0, -1})*delta_wolf) + (-0.5)*cpowq(M_W, 2)*cpowq(ccosq(beta), 2) + 0.5*cpowq(M_W, 2)*cpowq(csinq(beta), 2) + (-0.166667)*cpowq(M_W, 2)*cpowq(ccosq(beta), 2)*cpowq(ccosq(theta_W), -2)*cpowq(csinq(theta_W), 2) + 0.166667*cpowq(M_W, 2)*cpowq(ccosq(theta_W), -2)*cpowq(csinq(beta), 2)*cpowq(csinq(theta_W), 2),
            1.41421*A_b*M_W*V_tb*V_td*ccosq(beta)*csinq(theta_W)/e_em + 1.41421*A_b*M_W*V_tb*conjq(V_td)*ccosq(beta)*csinq(theta_W)/e_em,
            1.41421*A_b*M_W*V_tb*V_ts*ccosq(beta)*csinq(theta_W)/e_em + 1.41421*A_b*M_W*V_tb*conjq(V_ts)*ccosq(beta)*csinq(theta_W)/e_em,
            (-2)*m_b*mu_h*csinq(beta)/ccosq(beta) + 2.82843*A_b*M_W*cpowq(V_tb, 2)*ccosq(beta)*csinq(theta_W)/e_em,
            0,
            0,
            0,
            cpowq(m_d, 2) + V_td*conjq(V_td)*cpowq(M_qbR, 2) + V_cd*conjq(V_cd)*cpowq(M_qdR, 2) + cpowq(V_ud, 2)*cpowq(M_qdR, 2) + 0.5*V_cd*V_ud*del_DR_1 + 0.5*conjq(V_cd)*V_ud*del_DR_1 + 0.5*V_td*V_ud*del_DR_2 + 0.5*conjq(V_td)*V_ud*del_DR_2 + 0.5*conjq(V_cd)*V_td*del_DR_12 + 0.5*V_cd*conjq(V_td)*del_DR_12 + (-0.333333)*cpowq(M_W, 2)*cpowq(ccosq(beta), 2)*cpowq(ccosq(theta_W), -2)*cpowq(csinq(theta_W), 2) + 0.333333*cpowq(M_W, 2)*cpowq(ccosq(theta_W), -2)*cpowq(csinq(beta), 2)*cpowq(csinq(theta_W), 2),
            conjq(V_td)*V_ts*cpowq(M_qbR, 2) + V_td*conjq(V_ts)*cpowq(M_qbR, 2) + conjq(V_cd)*V_cs*cpowq(M_qdR, 2) + V_cd*conjq(V_cs)*cpowq(M_qdR, 2) + 2*V_ud*V_us*cpowq(M_qdR, 2) + 0.5*V_cs*V_ud*del_DR_1 + 0.5*conjq(V_cs)*V_ud*del_DR_1 + 0.5*V_cd*V_us*del_DR_1 + 0.5*conjq(V_cd)*V_us*del_DR_1 + 0.5*V_ts*V_ud*del_DR_2 + 0.5*conjq(V_ts)*V_ud*del_DR_2 + 0.5*V_td*V_us*del_DR_2 + 0.5*conjq(V_td)*V_us*del_DR_2 + 0.5*conjq(V_cs)*V_td*del_DR_12 + 0.5*V_cs*conjq(V_td)*del_DR_12 + 0.5*conjq(V_cd)*V_ts*del_DR_12 + 0.5*V_cd*conjq(V_ts)*del_DR_12,
            V_tb*V_td*cpowq(M_qbR, 2) + V_tb*conjq(V_td)*cpowq(M_qbR, 2) + V_cb*V_cd*cpowq(M_qdR, 2) + V_cb*conjq(V_cd)*cpowq(M_qdR, 2) + V_cb*V_ud*del_DR_1 + V_tb*V_ud*del_DR_2 + 0.5*V_cd*V_tb*del_DR_12 + 0.5*conjq(V_cd)*V_tb*del_DR_12 + 0.5*V_cb*V_td*del_DR_12 + 0.5*V_cb*conjq(V_td)*del_DR_12 + V_ud*cpowq(M_qdR, 2)*V_ub_mod*cexpq((complex_t{0, 1})*delta_wolf) + 0.5*V_cd*V_ub_mod*del_DR_1*cexpq((complex_t{0, 1})*delta_wolf) + 0.5*V_td*V_ub_mod*del_DR_2*cexpq((complex_t{0, 1})*delta_wolf) + V_ud*cpowq(M_qdR, 2)*V_ub_mod*cexpq((complex_t{0, -1})*delta_wolf) + 0.5*conjq(V_cd)*V_ub_mod*del_DR_1*cexpq((complex_t{0, -1})*delta_wolf) + 0.5*conjq(V_td)*V_ub_mod*del_DR_2*cexpq((complex_t{0, -1})*delta_wolf),
            0,
            0,
            0,
            0,
            cpowq(m_s, 2) + V_ts*conjq(V_ts)*cpowq(M_qbR, 2) + V_cs*conjq(V_cs)*cpowq(M_qdR, 2) + cpowq(V_us, 2)*cpowq(M_qdR, 2) + 0.5*V_cs*V_us*del_DR_1 + 0.5*conjq(V_cs)*V_us*del_DR_1 + 0.5*V_ts*V_us*del_DR_2 + 0.5*conjq(V_ts)*V_us*del_DR_2 + 0.5*conjq(V_cs)*V_ts*del_DR_12 + 0.5*V_cs*conjq(V_ts)*del_DR_12 + (-0.333333)*cpowq(M_W, 2)*cpowq(ccosq(beta), 2)*cpowq(ccosq(theta_W), -2)*cpowq(csinq(theta_W), 2) + 0.333333*cpowq(M_W, 2)*cpowq(ccosq(theta_W), -2)*cpowq(csinq(beta), 2)*cpowq(csinq(theta_W), 2),
            V_tb*V_ts*cpowq(M_qbR, 2) + V_tb*conjq(V_ts)*cpowq(M_qbR, 2) + V_cb*V_cs*cpowq(M_qdR, 2) + V_cb*conjq(V_cs)*cpowq(M_qdR, 2) + V_cb*V_us*del_DR_1 + V_tb*V_us*del_DR_2 + 0.5*V_cs*V_tb*del_DR_12 + 0.5*conjq(V_cs)*V_tb*del_DR_12 + 0.5*V_cb*V_ts*del_DR_12 + 0.5*V_cb*conjq(V_ts)*del_DR_12 + V_us*cpowq(M_qdR, 2)*V_ub_mod*cexpq((complex_t{0, 1})*delta_wolf) + 0.5*V_cs*V_ub_mod*del_DR_1*cexpq((complex_t{0, 1})*delta_wolf) + 0.5*V_ts*V_ub_mod*del_DR_2*cexpq((complex_t{0, 1})*delta_wolf) + V_us*cpowq(M_qdR, 2)*V_ub_mod*cexpq((complex_t{0, -1})*delta_wolf) + 0.5*conjq(V_cs)*V_ub_mod*del_DR_1*cexpq((complex_t{0, -1})*delta_wolf) + 0.5*conjq(V_ts)*V_ub_mod*del_DR_2*cexpq((complex_t{0, -1})*delta_wolf),
            0,
            0,
            0,
            0,
            0,
            cpowq(m_b, 2) + cpowq(V_tb, 2)*cpowq(M_qbR, 2) + cpowq(V_cb, 2)*cpowq(M_qdR, 2) + V_cb*V_tb*del_DR_12 + 0.5*V_cb*V_ub_mod*del_DR_1*cexpq((complex_t{0, 1})*delta_wolf) + 0.5*V_tb*V_ub_mod*del_DR_2*cexpq((complex_t{0, 1})*delta_wolf) + 0.5*V_cb*V_ub_mod*del_DR_1*cexpq((complex_t{0, -1})*delta_wolf) + 0.5*V_tb*V_ub_mod*del_DR_2*cexpq((complex_t{0, -1})*delta_wolf) + cpowq(M_qdR, 2)*cpowq(V_ub_mod, 2)*cexpq((complex_t{0, 1})*delta_wolf)*cexpq((complex_t{0, -1})*delta_wolf) + (-0.333333)*cpowq(M_W, 2)*cpowq(ccosq(beta), 2)*cpowq(ccosq(theta_W), -2)*cpowq(csinq(theta_W), 2) + 0.333333*cpowq(M_W, 2)*cpowq(ccosq(theta_W), -2)*cpowq(csinq(beta), 2)*cpowq(csinq(theta_W), 2),
        },
        {&outputs.U_sd_00, &outputs.U_sd_01, &outputs.U_sd_02, &outputs.U_sd_03, &outputs.U_sd_04, &outputs.U_sd_05, &outputs.U_sd_10, &outputs.U_sd_11, &outputs.U_sd_12, &outputs.U_sd_13, &outputs.U_sd_14, &outputs.U_sd_15, &outputs.U_sd_20, &outputs.U_sd_21, &outputs.U_sd_22, &outputs.U_sd_23, &outputs.U_sd_24, &outputs.U_sd_25, &outputs.U_sd_30, &outputs.U_sd_31, &outputs.U_sd_32, &outputs.U_sd_33, &outputs.U_sd_34, &outputs.U_sd_35, &outputs.U_sd_40, &outputs.U_sd_41, &outputs.U_sd_42, &outputs.U_sd_43, &outputs.U_sd_44, &outputs.U_sd_45, &outputs.U_sd_50, &outputs.U_sd_51, &outputs.U_sd_52, &outputs.U_sd_53, &outputs.U_sd_54, &outputs.U_sd_55, },
        {&outputs.m_sd_L, &outputs.m_ss_L, &outputs.m_sb_L, &outputs.m_sd_R, &outputs.m_ss_R, &outputs.m_sb_R, }
        );
    if (0 > outputs.m_sd_L) {
        std::cerr << "Warning: negative squared mass for " << "m_sd_L" << ".\n";
    }
    outputs.m_sd_L = sqrtq(outputs.m_sd_L);
    if (0 > outputs.m_ss_L) {
        std::cerr << "Warning: negative squared mass for " << "m_ss_L" << ".\n";
    }
    outputs.m_ss_L = sqrtq(outputs.m_ss_L);
    if (0 > outputs.m_sb_L) {
        std::cerr << "Warning: negative squared mass for " << "m_sb_L" << ".\n";
    }
    outputs.m_sb_L = sqrtq(outputs.m_sb_L);
    if (0 > outputs.m_sd_R) {
        std::cerr << "Warning: negative squared mass for " << "m_sd_R" << ".\n";
    }
    outputs.m_sd_R = sqrtq(outputs.m_sd_R);
    if (0 > outputs.m_ss_R) {
        std::cerr << "Warning: negative squared mass for " << "m_ss_R" << ".\n";
    }
    outputs.m_ss_R = sqrtq(outputs.m_ss_R);
    if (0 > outputs.m_sb_R) {
        std::cerr << "Warning: negative squared mass for " << "m_sb_R" << ".\n";
    }
    outputs.m_sb_R = sqrtq(outputs.m_sb_R);

    Diagonalizer::applyDiagonalization(
        {
            cpowq(m_e, 2) + cpowq(M_eL, 2) + (-0.5)*cpowq(M_W, 2)*cpowq(ccosq(beta), 2) + 0.5*cpowq(M_W, 2)*cpowq(csinq(beta), 2) + 0.5*cpowq(M_W, 2)*cpowq(ccosq(beta), 2)*cpowq(ccosq(theta_W), -2)*cpowq(csinq(theta_W), 2) + (-0.5)*cpowq(M_W, 2)*cpowq(ccosq(theta_W), -2)*cpowq(csinq(beta), 2)*cpowq(csinq(theta_W), 2),
            del_LL_1,
            del_LL_2,
            (-2)*m_e*mu_h*csinq(beta)/ccosq(beta),
            0,
            0,
            0,
            cpowq(M_eL, 2) + cpowq(m_mu, 2) + (-0.5)*cpowq(M_W, 2)*cpowq(ccosq(beta), 2) + 0.5*cpowq(M_W, 2)*cpowq(csinq(beta), 2) + 0.5*cpowq(M_W, 2)*cpowq(ccosq(beta), 2)*cpowq(ccosq(theta_W), -2)*cpowq(csinq(theta_W), 2) + (-0.5)*cpowq(M_W, 2)*cpowq(ccosq(theta_W), -2)*cpowq(csinq(beta), 2)*cpowq(csinq(theta_W), 2),
            del_LL_12,
            0,
            (-2)*m_mu*mu_h*csinq(beta)/ccosq(beta),
            0,
            0,
            0,
            cpowq(m_tau, 2) + cpowq(M_tauL, 2) + (-0.5)*cpowq(M_W, 2)*cpowq(ccosq(beta), 2) + 0.5*cpowq(M_W, 2)*cpowq(csinq(beta), 2) + 0.5*cpowq(M_W, 2)*cpowq(ccosq(beta), 2)*cpowq(ccosq(theta_W), -2)*cpowq(csinq(theta_W), 2) + (-0.5)*cpowq(M_W, 2)*cpowq(ccosq(theta_W), -2)*cpowq(csinq(beta), 2)*cpowq(csinq(theta_W), 2),
            0,
            0,
            (-2)*mu_h*m_tau*csinq(beta)/ccosq(beta) + 2.82843*M_W*A_tau*ccosq(beta)*csinq(theta_W)/e_em,
            0,
            0,
            0,
            cpowq(m_e, 2) + cpowq(M_eR, 2) + -cpowq(M_W, 2)*cpowq(ccosq(beta), 2)*cpowq(ccosq(theta_W), -2)*cpowq(csinq(theta_W), 2) + cpowq(M_W, 2)*cpowq(ccosq(theta_W), -2)*cpowq(csinq(beta), 2)*cpowq(csinq(theta_W), 2),
            del_ER_1,
            del_ER_2,
            0,
            0,
            0,
            0,
            cpowq(M_eR, 2) + cpowq(m_mu, 2) + -cpowq(M_W, 2)*cpowq(ccosq(beta), 2)*cpowq(ccosq(theta_W), -2)*cpowq(csinq(theta_W), 2) + cpowq(M_W, 2)*cpowq(ccosq(theta_W), -2)*cpowq(csinq(beta), 2)*cpowq(csinq(theta_W), 2),
            del_ER_12,
            0,
            0,
            0,
            0,
            0,
            cpowq(m_tau, 2) + cpowq(M_tauR, 2) + -cpowq(M_W, 2)*cpowq(ccosq(beta), 2)*cpowq(ccosq(theta_W), -2)*cpowq(csinq(theta_W), 2) + cpowq(M_W, 2)*cpowq(ccosq(theta_W), -2)*cpowq(csinq(beta), 2)*cpowq(csinq(theta_W), 2),
        },
        {&outputs.U_se_00, &outputs.U_se_01, &outputs.U_se_02, &outputs.U_se_03, &outputs.U_se_04, &outputs.U_se_05, &outputs.U_se_10, &outputs.U_se_11, &outputs.U_se_12, &outputs.U_se_13, &outputs.U_se_14, &outputs.U_se_15, &outputs.U_se_20, &outputs.U_se_21, &outputs.U_se_22, &outputs.U_se_23, &outputs.U_se_24, &outputs.U_se_25, &outputs.U_se_30, &outputs.U_se_31, &outputs.U_se_32, &outputs.U_se_33, &outputs.U_se_34, &outputs.U_se_35, &outputs.U_se_40, &outputs.U_se_41, &outputs.U_se_42, &outputs.U_se_43, &outputs.U_se_44, &outputs.U_se_45, &outputs.U_se_50, &outputs.U_se_51, &outputs.U_se_52, &outputs.U_se_53, &outputs.U_se_54, &outputs.U_se_55, },
        {&outputs.m_se_L, &outputs.m_smu_L, &outputs.m_stau_L, &outputs.m_se_R, &outputs.m_smu_R, &outputs.m_stau_R, }
        );
    if (0 > outputs.m_se_L) {
        std::cerr << "Warning: negative squared mass for " << "m_se_L" << ".\n";
    }
    outputs.m_se_L = sqrtq(outputs.m_se_L);
    if (0 > outputs.m_smu_L) {
        std::cerr << "Warning: negative squared mass for " << "m_smu_L" << ".\n";
    }
    outputs.m_smu_L = sqrtq(outputs.m_smu_L);
    if (0 > outputs.m_stau_L) {
        std::cerr << "Warning: negative squared mass for " << "m_stau_L" << ".\n";
    }
    outputs.m_stau_L = sqrtq(outputs.m_stau_L);
    if (0 > outputs.m_se_R) {
        std::cerr << "Warning: negative squared mass for " << "m_se_R" << ".\n";
    }
    outputs.m_se_R = sqrtq(outputs.m_se_R);
    if (0 > outputs.m_smu_R) {
        std::cerr << "Warning: negative squared mass for " << "m_smu_R" << ".\n";
    }
    outputs.m_smu_R = sqrtq(outputs.m_smu_R);
    if (0 > outputs.m_stau_R) {
        std::cerr << "Warning: negative squared mass for " << "m_stau_R" << ".\n";
    }
    outputs.m_stau_R = sqrtq(outputs.m_stau_R);

    Diagonalizer::applyDiagonalization(
        {
            cpowq(M_eL, 2) + 0.5*cpowq(M_W, 2)*cpowq(ccosq(beta), 2) + (-0.5)*cpowq(M_W, 2)*cpowq(csinq(beta), 2) + 0.5*cpowq(M_W, 2)*cpowq(ccosq(beta), 2)*cpowq(ccosq(theta_W), -2)*cpowq(csinq(theta_W), 2) + (-0.5)*cpowq(M_W, 2)*cpowq(ccosq(theta_W), -2)*cpowq(csinq(beta), 2)*cpowq(csinq(theta_W), 2),
            del_LL_1,
            del_LL_2,
            0,
            cpowq(M_eL, 2) + 0.5*cpowq(M_W, 2)*cpowq(ccosq(beta), 2) + (-0.5)*cpowq(M_W, 2)*cpowq(csinq(beta), 2) + 0.5*cpowq(M_W, 2)*cpowq(ccosq(beta), 2)*cpowq(ccosq(theta_W), -2)*cpowq(csinq(theta_W), 2) + (-0.5)*cpowq(M_W, 2)*cpowq(ccosq(theta_W), -2)*cpowq(csinq(beta), 2)*cpowq(csinq(theta_W), 2),
            del_LL_12,
            0,
            0,
            cpowq(M_tauL, 2) + 0.5*cpowq(M_W, 2)*cpowq(ccosq(beta), 2) + (-0.5)*cpowq(M_W, 2)*cpowq(csinq(beta), 2) + 0.5*cpowq(M_W, 2)*cpowq(ccosq(beta), 2)*cpowq(ccosq(theta_W), -2)*cpowq(csinq(theta_W), 2) + (-0.5)*cpowq(M_W, 2)*cpowq(ccosq(theta_W), -2)*cpowq(csinq(beta), 2)*cpowq(csinq(theta_W), 2),
        },
        {&outputs.U_snu_00, &outputs.U_snu_01, &outputs.U_snu_02, &outputs.U_snu_10, &outputs.U_snu_11, &outputs.U_snu_12, &outputs.U_snu_20, &outputs.U_snu_21, &outputs.U_snu_22, },
        {&outputs.m_snu_e, &outputs.m_snu_mu, &outputs.m_snu_tau, }
        );
    if (0 > outputs.m_snu_e) {
        std::cerr << "Warning: negative squared mass for " << "m_snu_e" << ".\n";
    }
    outputs.m_snu_e = sqrtq(outputs.m_snu_e);
    if (0 > outputs.m_snu_mu) {
        std::cerr << "Warning: negative squared mass for " << "m_snu_mu" << ".\n";
    }
    outputs.m_snu_mu = sqrtq(outputs.m_snu_mu);
    if (0 > outputs.m_snu_tau) {
        std::cerr << "Warning: negative squared mass for " << "m_snu_tau" << ".\n";
    }
    outputs.m_snu_tau = sqrtq(outputs.m_snu_tau);

    return outputs;
}

void updateMassExpressions(param_t &params)
{
    params.M_Z = crealq(m_Z(params));
    params.m_sG = crealq(m_sG(params));
    params.m_N_1 = crealq(m_N_1(params));
    params.m_N_2 = crealq(m_N_2(params));
    params.m_N_3 = crealq(m_N_3(params));
    params.m_N_4 = crealq(m_N_4(params));
    params.m_Gp = crealq(m_Gp(params));
    params.m_G0 = crealq(m_G0(params));
}

} // End of namespace c9_nmfv

