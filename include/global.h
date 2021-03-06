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

#ifndef CSL_LIB_GLOBAL
#define CSL_LIB_GLOBAL
#include "params.h"
#include "common.h"

namespace c9_nmfv {

void updateSpectrum(param_t &params);

struct SpectrumInput;
struct SpectrumOutput;

SpectrumOutput updateDiagonalization(SpectrumInput const&);

void updateDiagonalization(param_t &params);

////////////////////////////////////////////////////
// Here are the parameters to set before calling    
// updateDiagonalization()                          
////////////////////////////////////////////////////
struct SpectrumInput {
    complex_t A_b;
    complex_t A_t;
    complex_t A_tau;
    complex_t M_1;
    complex_t M_2;
    complex_t M_W;
    complex_t M_eL;
    complex_t M_eR;
    complex_t M_q1L;
    complex_t M_q3L;
    complex_t M_qbR;
    complex_t M_qdR;
    complex_t M_qtR;
    complex_t M_quR;
    complex_t M_tauL;
    complex_t M_tauR;
    complex_t V_cb;
    complex_t V_cd;
    complex_t V_cs;
    complex_t V_tb;
    complex_t V_td;
    complex_t V_ts;
    complex_t V_ub_mod;
    complex_t V_ud;
    complex_t V_us;
    complex_t beta;
    complex_t del_DR_1;
    complex_t del_DR_12;
    complex_t del_DR_2;
    complex_t del_ER_1;
    complex_t del_ER_12;
    complex_t del_ER_2;
    complex_t del_LL_1;
    complex_t del_LL_12;
    complex_t del_LL_2;
    complex_t del_QL_1;
    complex_t del_QL_12;
    complex_t del_QL_2;
    complex_t del_UR_1;
    complex_t del_UR_12;
    complex_t del_UR_2;
    complex_t delta_wolf;
    complex_t e_em;
    complex_t m_b;
    complex_t m_c;
    complex_t m_d;
    complex_t m_e;
    complex_t m_mu;
    complex_t m_s;
    complex_t m_t;
    complex_t m_tau;
    complex_t m_u;
    complex_t mu_h;
    complex_t theta_W;
};

////////////////////////////////////////////////////
// Here are the masses and mixings                 
// result of the diagonalization                   
////////////////////////////////////////////////////
struct SpectrumOutput {
    real_t m_N_1;
    real_t m_N_2;
    real_t m_N_3;
    real_t m_N_4;
    real_t m_C1p;
    real_t m_C2p;
    real_t m_su_L;
    real_t m_sc_L;
    real_t m_st_L;
    real_t m_su_R;
    real_t m_sc_R;
    real_t m_st_R;
    real_t m_sd_L;
    real_t m_ss_L;
    real_t m_sb_L;
    real_t m_sd_R;
    real_t m_ss_R;
    real_t m_sb_R;
    real_t m_se_L;
    real_t m_smu_L;
    real_t m_stau_L;
    real_t m_se_R;
    real_t m_smu_R;
    real_t m_stau_R;
    real_t m_snu_e;
    real_t m_snu_mu;
    real_t m_snu_tau;

    complex_t N_B1;
    complex_t N_B2;
    complex_t N_B3;
    complex_t N_B4;
    complex_t N_W1;
    complex_t N_W2;
    complex_t N_W3;
    complex_t N_W4;
    complex_t N_d1;
    complex_t N_d2;
    complex_t N_d3;
    complex_t N_d4;
    complex_t N_u1;
    complex_t N_u2;
    complex_t N_u3;
    complex_t N_u4;
    complex_t U_Wm1;
    complex_t U_Wm2;
    complex_t U_d1;
    complex_t U_d2;
    complex_t U_sd_00;
    complex_t U_sd_01;
    complex_t U_sd_02;
    complex_t U_sd_03;
    complex_t U_sd_04;
    complex_t U_sd_05;
    complex_t U_sd_10;
    complex_t U_sd_11;
    complex_t U_sd_12;
    complex_t U_sd_13;
    complex_t U_sd_14;
    complex_t U_sd_15;
    complex_t U_sd_20;
    complex_t U_sd_21;
    complex_t U_sd_22;
    complex_t U_sd_23;
    complex_t U_sd_24;
    complex_t U_sd_25;
    complex_t U_sd_30;
    complex_t U_sd_31;
    complex_t U_sd_32;
    complex_t U_sd_33;
    complex_t U_sd_34;
    complex_t U_sd_35;
    complex_t U_sd_40;
    complex_t U_sd_41;
    complex_t U_sd_42;
    complex_t U_sd_43;
    complex_t U_sd_44;
    complex_t U_sd_45;
    complex_t U_sd_50;
    complex_t U_sd_51;
    complex_t U_sd_52;
    complex_t U_sd_53;
    complex_t U_sd_54;
    complex_t U_sd_55;
    complex_t U_se_00;
    complex_t U_se_01;
    complex_t U_se_02;
    complex_t U_se_03;
    complex_t U_se_04;
    complex_t U_se_05;
    complex_t U_se_10;
    complex_t U_se_11;
    complex_t U_se_12;
    complex_t U_se_13;
    complex_t U_se_14;
    complex_t U_se_15;
    complex_t U_se_20;
    complex_t U_se_21;
    complex_t U_se_22;
    complex_t U_se_23;
    complex_t U_se_24;
    complex_t U_se_25;
    complex_t U_se_30;
    complex_t U_se_31;
    complex_t U_se_32;
    complex_t U_se_33;
    complex_t U_se_34;
    complex_t U_se_35;
    complex_t U_se_40;
    complex_t U_se_41;
    complex_t U_se_42;
    complex_t U_se_43;
    complex_t U_se_44;
    complex_t U_se_45;
    complex_t U_se_50;
    complex_t U_se_51;
    complex_t U_se_52;
    complex_t U_se_53;
    complex_t U_se_54;
    complex_t U_se_55;
    complex_t U_snu_00;
    complex_t U_snu_01;
    complex_t U_snu_02;
    complex_t U_snu_10;
    complex_t U_snu_11;
    complex_t U_snu_12;
    complex_t U_snu_20;
    complex_t U_snu_21;
    complex_t U_snu_22;
    complex_t U_su_00;
    complex_t U_su_01;
    complex_t U_su_02;
    complex_t U_su_03;
    complex_t U_su_04;
    complex_t U_su_05;
    complex_t U_su_10;
    complex_t U_su_11;
    complex_t U_su_12;
    complex_t U_su_13;
    complex_t U_su_14;
    complex_t U_su_15;
    complex_t U_su_20;
    complex_t U_su_21;
    complex_t U_su_22;
    complex_t U_su_23;
    complex_t U_su_24;
    complex_t U_su_25;
    complex_t U_su_30;
    complex_t U_su_31;
    complex_t U_su_32;
    complex_t U_su_33;
    complex_t U_su_34;
    complex_t U_su_35;
    complex_t U_su_40;
    complex_t U_su_41;
    complex_t U_su_42;
    complex_t U_su_43;
    complex_t U_su_44;
    complex_t U_su_45;
    complex_t U_su_50;
    complex_t U_su_51;
    complex_t U_su_52;
    complex_t U_su_53;
    complex_t U_su_54;
    complex_t U_su_55;
    complex_t V_Wp1;
    complex_t V_Wp2;
    complex_t V_u1;
    complex_t V_u2;
};

////////////////////////////////////////////////////
// Here is a generic function to read results      
// of the diagonalization in a corresponding struct
////////////////////////////////////////////////////

template<class Type>
void readDiagonalizationInputs(
        SpectrumInput &diagData,
        Type    const &input
        )
{
    diagData.A_b = input.A_b;
    diagData.A_t = input.A_t;
    diagData.A_tau = input.A_tau;
    diagData.M_1 = input.M_1;
    diagData.M_2 = input.M_2;
    diagData.M_W = input.M_W;
    diagData.M_eL = input.M_eL;
    diagData.M_eR = input.M_eR;
    diagData.M_q1L = input.M_q1L;
    diagData.M_q3L = input.M_q3L;
    diagData.M_qbR = input.M_qbR;
    diagData.M_qdR = input.M_qdR;
    diagData.M_qtR = input.M_qtR;
    diagData.M_quR = input.M_quR;
    diagData.M_tauL = input.M_tauL;
    diagData.M_tauR = input.M_tauR;
    diagData.V_cb = input.V_cb;
    diagData.V_cd = input.V_cd;
    diagData.V_cs = input.V_cs;
    diagData.V_tb = input.V_tb;
    diagData.V_td = input.V_td;
    diagData.V_ts = input.V_ts;
    diagData.V_ub_mod = input.V_ub_mod;
    diagData.V_ud = input.V_ud;
    diagData.V_us = input.V_us;
    diagData.beta = input.beta;
    diagData.del_DR_1 = input.del_DR_1;
    diagData.del_DR_12 = input.del_DR_12;
    diagData.del_DR_2 = input.del_DR_2;
    diagData.del_ER_1 = input.del_ER_1;
    diagData.del_ER_12 = input.del_ER_12;
    diagData.del_ER_2 = input.del_ER_2;
    diagData.del_LL_1 = input.del_LL_1;
    diagData.del_LL_12 = input.del_LL_12;
    diagData.del_LL_2 = input.del_LL_2;
    diagData.del_QL_1 = input.del_QL_1;
    diagData.del_QL_12 = input.del_QL_12;
    diagData.del_QL_2 = input.del_QL_2;
    diagData.del_UR_1 = input.del_UR_1;
    diagData.del_UR_12 = input.del_UR_12;
    diagData.del_UR_2 = input.del_UR_2;
    diagData.delta_wolf = input.delta_wolf;
    diagData.e_em = input.e_em;
    diagData.m_b = input.m_b;
    diagData.m_c = input.m_c;
    diagData.m_d = input.m_d;
    diagData.m_e = input.m_e;
    diagData.m_mu = input.m_mu;
    diagData.m_s = input.m_s;
    diagData.m_t = input.m_t;
    diagData.m_tau = input.m_tau;
    diagData.m_u = input.m_u;
    diagData.mu_h = input.mu_h;
    diagData.theta_W = input.theta_W;
}

template<class Type>
void readDiagonalizationOutputs(
        SpectrumOutput const &diagData,
        Type                 &output
        )
{
    output.m_N_1 = diagData.m_N_1;
    output.m_N_2 = diagData.m_N_2;
    output.m_N_3 = diagData.m_N_3;
    output.m_N_4 = diagData.m_N_4;
    output.m_C1p = diagData.m_C1p;
    output.m_C2p = diagData.m_C2p;
    output.m_su_L = diagData.m_su_L;
    output.m_sc_L = diagData.m_sc_L;
    output.m_st_L = diagData.m_st_L;
    output.m_su_R = diagData.m_su_R;
    output.m_sc_R = diagData.m_sc_R;
    output.m_st_R = diagData.m_st_R;
    output.m_sd_L = diagData.m_sd_L;
    output.m_ss_L = diagData.m_ss_L;
    output.m_sb_L = diagData.m_sb_L;
    output.m_sd_R = diagData.m_sd_R;
    output.m_ss_R = diagData.m_ss_R;
    output.m_sb_R = diagData.m_sb_R;
    output.m_se_L = diagData.m_se_L;
    output.m_smu_L = diagData.m_smu_L;
    output.m_stau_L = diagData.m_stau_L;
    output.m_se_R = diagData.m_se_R;
    output.m_smu_R = diagData.m_smu_R;
    output.m_stau_R = diagData.m_stau_R;
    output.m_snu_e = diagData.m_snu_e;
    output.m_snu_mu = diagData.m_snu_mu;
    output.m_snu_tau = diagData.m_snu_tau;
    output.N_B1 = diagData.N_B1;
    output.N_B2 = diagData.N_B2;
    output.N_B3 = diagData.N_B3;
    output.N_B4 = diagData.N_B4;
    output.N_W1 = diagData.N_W1;
    output.N_W2 = diagData.N_W2;
    output.N_W3 = diagData.N_W3;
    output.N_W4 = diagData.N_W4;
    output.N_d1 = diagData.N_d1;
    output.N_d2 = diagData.N_d2;
    output.N_d3 = diagData.N_d3;
    output.N_d4 = diagData.N_d4;
    output.N_u1 = diagData.N_u1;
    output.N_u2 = diagData.N_u2;
    output.N_u3 = diagData.N_u3;
    output.N_u4 = diagData.N_u4;
    output.U_Wm1 = diagData.U_Wm1;
    output.U_Wm2 = diagData.U_Wm2;
    output.U_d1 = diagData.U_d1;
    output.U_d2 = diagData.U_d2;
    output.U_sd_00 = diagData.U_sd_00;
    output.U_sd_01 = diagData.U_sd_01;
    output.U_sd_02 = diagData.U_sd_02;
    output.U_sd_03 = diagData.U_sd_03;
    output.U_sd_04 = diagData.U_sd_04;
    output.U_sd_05 = diagData.U_sd_05;
    output.U_sd_10 = diagData.U_sd_10;
    output.U_sd_11 = diagData.U_sd_11;
    output.U_sd_12 = diagData.U_sd_12;
    output.U_sd_13 = diagData.U_sd_13;
    output.U_sd_14 = diagData.U_sd_14;
    output.U_sd_15 = diagData.U_sd_15;
    output.U_sd_20 = diagData.U_sd_20;
    output.U_sd_21 = diagData.U_sd_21;
    output.U_sd_22 = diagData.U_sd_22;
    output.U_sd_23 = diagData.U_sd_23;
    output.U_sd_24 = diagData.U_sd_24;
    output.U_sd_25 = diagData.U_sd_25;
    output.U_sd_30 = diagData.U_sd_30;
    output.U_sd_31 = diagData.U_sd_31;
    output.U_sd_32 = diagData.U_sd_32;
    output.U_sd_33 = diagData.U_sd_33;
    output.U_sd_34 = diagData.U_sd_34;
    output.U_sd_35 = diagData.U_sd_35;
    output.U_sd_40 = diagData.U_sd_40;
    output.U_sd_41 = diagData.U_sd_41;
    output.U_sd_42 = diagData.U_sd_42;
    output.U_sd_43 = diagData.U_sd_43;
    output.U_sd_44 = diagData.U_sd_44;
    output.U_sd_45 = diagData.U_sd_45;
    output.U_sd_50 = diagData.U_sd_50;
    output.U_sd_51 = diagData.U_sd_51;
    output.U_sd_52 = diagData.U_sd_52;
    output.U_sd_53 = diagData.U_sd_53;
    output.U_sd_54 = diagData.U_sd_54;
    output.U_sd_55 = diagData.U_sd_55;
    output.U_se_00 = diagData.U_se_00;
    output.U_se_01 = diagData.U_se_01;
    output.U_se_02 = diagData.U_se_02;
    output.U_se_03 = diagData.U_se_03;
    output.U_se_04 = diagData.U_se_04;
    output.U_se_05 = diagData.U_se_05;
    output.U_se_10 = diagData.U_se_10;
    output.U_se_11 = diagData.U_se_11;
    output.U_se_12 = diagData.U_se_12;
    output.U_se_13 = diagData.U_se_13;
    output.U_se_14 = diagData.U_se_14;
    output.U_se_15 = diagData.U_se_15;
    output.U_se_20 = diagData.U_se_20;
    output.U_se_21 = diagData.U_se_21;
    output.U_se_22 = diagData.U_se_22;
    output.U_se_23 = diagData.U_se_23;
    output.U_se_24 = diagData.U_se_24;
    output.U_se_25 = diagData.U_se_25;
    output.U_se_30 = diagData.U_se_30;
    output.U_se_31 = diagData.U_se_31;
    output.U_se_32 = diagData.U_se_32;
    output.U_se_33 = diagData.U_se_33;
    output.U_se_34 = diagData.U_se_34;
    output.U_se_35 = diagData.U_se_35;
    output.U_se_40 = diagData.U_se_40;
    output.U_se_41 = diagData.U_se_41;
    output.U_se_42 = diagData.U_se_42;
    output.U_se_43 = diagData.U_se_43;
    output.U_se_44 = diagData.U_se_44;
    output.U_se_45 = diagData.U_se_45;
    output.U_se_50 = diagData.U_se_50;
    output.U_se_51 = diagData.U_se_51;
    output.U_se_52 = diagData.U_se_52;
    output.U_se_53 = diagData.U_se_53;
    output.U_se_54 = diagData.U_se_54;
    output.U_se_55 = diagData.U_se_55;
    output.U_snu_00 = diagData.U_snu_00;
    output.U_snu_01 = diagData.U_snu_01;
    output.U_snu_02 = diagData.U_snu_02;
    output.U_snu_10 = diagData.U_snu_10;
    output.U_snu_11 = diagData.U_snu_11;
    output.U_snu_12 = diagData.U_snu_12;
    output.U_snu_20 = diagData.U_snu_20;
    output.U_snu_21 = diagData.U_snu_21;
    output.U_snu_22 = diagData.U_snu_22;
    output.U_su_00 = diagData.U_su_00;
    output.U_su_01 = diagData.U_su_01;
    output.U_su_02 = diagData.U_su_02;
    output.U_su_03 = diagData.U_su_03;
    output.U_su_04 = diagData.U_su_04;
    output.U_su_05 = diagData.U_su_05;
    output.U_su_10 = diagData.U_su_10;
    output.U_su_11 = diagData.U_su_11;
    output.U_su_12 = diagData.U_su_12;
    output.U_su_13 = diagData.U_su_13;
    output.U_su_14 = diagData.U_su_14;
    output.U_su_15 = diagData.U_su_15;
    output.U_su_20 = diagData.U_su_20;
    output.U_su_21 = diagData.U_su_21;
    output.U_su_22 = diagData.U_su_22;
    output.U_su_23 = diagData.U_su_23;
    output.U_su_24 = diagData.U_su_24;
    output.U_su_25 = diagData.U_su_25;
    output.U_su_30 = diagData.U_su_30;
    output.U_su_31 = diagData.U_su_31;
    output.U_su_32 = diagData.U_su_32;
    output.U_su_33 = diagData.U_su_33;
    output.U_su_34 = diagData.U_su_34;
    output.U_su_35 = diagData.U_su_35;
    output.U_su_40 = diagData.U_su_40;
    output.U_su_41 = diagData.U_su_41;
    output.U_su_42 = diagData.U_su_42;
    output.U_su_43 = diagData.U_su_43;
    output.U_su_44 = diagData.U_su_44;
    output.U_su_45 = diagData.U_su_45;
    output.U_su_50 = diagData.U_su_50;
    output.U_su_51 = diagData.U_su_51;
    output.U_su_52 = diagData.U_su_52;
    output.U_su_53 = diagData.U_su_53;
    output.U_su_54 = diagData.U_su_54;
    output.U_su_55 = diagData.U_su_55;
    output.V_Wp1 = diagData.V_Wp1;
    output.V_Wp2 = diagData.V_Wp2;
    output.V_u1 = diagData.V_u1;
    output.V_u2 = diagData.V_u2;
}
void updateMassExpressions(param_t &params);


} // End of namespace c9_nmfv

#endif
