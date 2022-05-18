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

#ifndef CSL_LIB_PARAM_H_INCLUDED
#define CSL_LIB_PARAM_H_INCLUDED

#include <map>
#include <array>
#include "common.h"
#include "libcomplexop.h"
#include "csl/initSanitizer.h"

namespace c9_nmfv {

struct param_t {

    ///////////////////////////////////////
    // Elementary parameters to be defined 
    ///////////////////////////////////////

    csl::InitSanitizer<real_t> M_3 { "M_3" };
    csl::InitSanitizer<real_t> g_s { "g_s" };
    csl::InitSanitizer<real_t> V_cb { "V_cb" };
    csl::InitSanitizer<real_t> V_tb { "V_tb" };
    csl::InitSanitizer<real_t> V_us { "V_us" };
    csl::InitSanitizer<real_t> beta { "beta" };
    csl::InitSanitizer<real_t> e_em { "e_em" };
    csl::InitSanitizer<real_t> s_12 { "s_12" };
    csl::InitSanitizer<real_t> s_34 { "s_34" };
    csl::InitSanitizer<real_t> alpha { "alpha" };
    csl::InitSanitizer<real_t> Finite { "Finite" };
    csl::InitSanitizer<real_t> theta_W { "theta_W" };
    csl::InitSanitizer<real_t> V_ub_mod { "V_ub_mod" };
    csl::InitSanitizer<real_t> reg_prop { "reg_prop" };
    csl::InitSanitizer<real_t> delta_wolf { "delta_wolf" };
    csl::InitSanitizer<complex_t> A_b { "A_b" };
    csl::InitSanitizer<complex_t> A_t { "A_t" };
    csl::InitSanitizer<complex_t> M_1 { "M_1" };
    csl::InitSanitizer<complex_t> M_2 { "M_2" };
    csl::InitSanitizer<complex_t> M_eL { "M_eL" };
    csl::InitSanitizer<complex_t> M_eR { "M_eR" };
    csl::InitSanitizer<complex_t> V_cd { "V_cd" };
    csl::InitSanitizer<complex_t> V_cs { "V_cs" };
    csl::InitSanitizer<complex_t> V_td { "V_td" };
    csl::InitSanitizer<complex_t> V_ts { "V_ts" };
    csl::InitSanitizer<complex_t> V_ud { "V_ud" };
    csl::InitSanitizer<complex_t> mu_h { "mu_h" };
    csl::InitSanitizer<complex_t> A_tau { "A_tau" };
    csl::InitSanitizer<complex_t> M_q1L { "M_q1L" };
    csl::InitSanitizer<complex_t> M_q3L { "M_q3L" };
    csl::InitSanitizer<complex_t> M_qbR { "M_qbR" };
    csl::InitSanitizer<complex_t> M_qdR { "M_qdR" };
    csl::InitSanitizer<complex_t> M_qtR { "M_qtR" };
    csl::InitSanitizer<complex_t> M_quR { "M_quR" };
    csl::InitSanitizer<complex_t> M_tauL { "M_tauL" };
    csl::InitSanitizer<complex_t> M_tauR { "M_tauR" };
    csl::InitSanitizer<complex_t> del_DR_1 { "del_DR_1" };
    csl::InitSanitizer<complex_t> del_DR_2 { "del_DR_2" };
    csl::InitSanitizer<complex_t> del_ER_1 { "del_ER_1" };
    csl::InitSanitizer<complex_t> del_ER_2 { "del_ER_2" };
    csl::InitSanitizer<complex_t> del_LL_1 { "del_LL_1" };
    csl::InitSanitizer<complex_t> del_LL_2 { "del_LL_2" };
    csl::InitSanitizer<complex_t> del_QL_1 { "del_QL_1" };
    csl::InitSanitizer<complex_t> del_QL_2 { "del_QL_2" };
    csl::InitSanitizer<complex_t> del_UR_1 { "del_UR_1" };
    csl::InitSanitizer<complex_t> del_UR_2 { "del_UR_2" };
    csl::InitSanitizer<complex_t> del_DR_12 { "del_DR_12" };
    csl::InitSanitizer<complex_t> del_ER_12 { "del_ER_12" };
    csl::InitSanitizer<complex_t> del_LL_12 { "del_LL_12" };
    csl::InitSanitizer<complex_t> del_QL_12 { "del_QL_12" };
    csl::InitSanitizer<complex_t> del_UR_12 { "del_UR_12" };


    ///////////////////////////////////////
    // Parameters functions of others  
    // through diagonalization or mass 
    // expressions, see updateSpectrum()  
    // in global.h or set them by hand  
    // 
    // And other default parameters  
    ///////////////////////////////////////

    csl::InitSanitizer<real_t> M_W { "M_W" };
    csl::InitSanitizer<real_t> M_Z { "M_Z" };
    csl::InitSanitizer<real_t> m_b { "m_b" };
    csl::InitSanitizer<real_t> m_c { "m_c" };
    csl::InitSanitizer<real_t> m_d { "m_d" };
    csl::InitSanitizer<real_t> m_e { "m_e" };
    csl::InitSanitizer<real_t> m_s { "m_s" };
    csl::InitSanitizer<real_t> m_t { "m_t" };
    csl::InitSanitizer<real_t> m_u { "m_u" };
    csl::InitSanitizer<real_t> m_A0 { "m_A0" };
    csl::InitSanitizer<real_t> m_G0 { "m_G0" };
    csl::InitSanitizer<real_t> m_Gp { "m_Gp" };
    csl::InitSanitizer<real_t> m_H0 { "m_H0" };
    csl::InitSanitizer<real_t> m_Hp { "m_Hp" };
    csl::InitSanitizer<real_t> m_h0 { "m_h0" };
    csl::InitSanitizer<real_t> m_mu { "m_mu" };
    csl::InitSanitizer<real_t> m_sG { "m_sG" };
    csl::InitSanitizer<real_t> m_C1p { "m_C1p" };
    csl::InitSanitizer<real_t> m_C2p { "m_C2p" };
    csl::InitSanitizer<real_t> m_N_1 { "m_N_1" };
    csl::InitSanitizer<real_t> m_N_2 { "m_N_2" };
    csl::InitSanitizer<real_t> m_N_3 { "m_N_3" };
    csl::InitSanitizer<real_t> m_N_4 { "m_N_4" };
    csl::InitSanitizer<real_t> m_tau { "m_tau" };
    csl::InitSanitizer<real_t> m_sb_L { "m_sb_L" };
    csl::InitSanitizer<real_t> m_sb_R { "m_sb_R" };
    csl::InitSanitizer<real_t> m_sc_L { "m_sc_L" };
    csl::InitSanitizer<real_t> m_sc_R { "m_sc_R" };
    csl::InitSanitizer<real_t> m_sd_L { "m_sd_L" };
    csl::InitSanitizer<real_t> m_sd_R { "m_sd_R" };
    csl::InitSanitizer<real_t> m_se_L { "m_se_L" };
    csl::InitSanitizer<real_t> m_se_R { "m_se_R" };
    csl::InitSanitizer<real_t> m_ss_L { "m_ss_L" };
    csl::InitSanitizer<real_t> m_ss_R { "m_ss_R" };
    csl::InitSanitizer<real_t> m_st_L { "m_st_L" };
    csl::InitSanitizer<real_t> m_st_R { "m_st_R" };
    csl::InitSanitizer<real_t> m_su_L { "m_su_L" };
    csl::InitSanitizer<real_t> m_su_R { "m_su_R" };
    csl::InitSanitizer<real_t> m_smu_L { "m_smu_L" };
    csl::InitSanitizer<real_t> m_smu_R { "m_smu_R" };
    csl::InitSanitizer<real_t> m_snu_e { "m_snu_e" };
    csl::InitSanitizer<real_t> m_snu_mu { "m_snu_mu" };
    csl::InitSanitizer<real_t> m_stau_L { "m_stau_L" };
    csl::InitSanitizer<real_t> m_stau_R { "m_stau_R" };
    csl::InitSanitizer<real_t> m_snu_tau { "m_snu_tau" };
    csl::InitSanitizer<complex_t> N_B1 { "N_B1" };
    csl::InitSanitizer<complex_t> N_B2 { "N_B2" };
    csl::InitSanitizer<complex_t> N_B3 { "N_B3" };
    csl::InitSanitizer<complex_t> N_B4 { "N_B4" };
    csl::InitSanitizer<complex_t> N_W1 { "N_W1" };
    csl::InitSanitizer<complex_t> N_W2 { "N_W2" };
    csl::InitSanitizer<complex_t> N_W3 { "N_W3" };
    csl::InitSanitizer<complex_t> N_W4 { "N_W4" };
    csl::InitSanitizer<complex_t> N_d1 { "N_d1" };
    csl::InitSanitizer<complex_t> N_d2 { "N_d2" };
    csl::InitSanitizer<complex_t> N_d3 { "N_d3" };
    csl::InitSanitizer<complex_t> N_d4 { "N_d4" };
    csl::InitSanitizer<complex_t> N_u1 { "N_u1" };
    csl::InitSanitizer<complex_t> N_u2 { "N_u2" };
    csl::InitSanitizer<complex_t> N_u3 { "N_u3" };
    csl::InitSanitizer<complex_t> N_u4 { "N_u4" };
    csl::InitSanitizer<complex_t> U_d1 { "U_d1" };
    csl::InitSanitizer<complex_t> U_d2 { "U_d2" };
    csl::InitSanitizer<complex_t> V_u1 { "V_u1" };
    csl::InitSanitizer<complex_t> V_u2 { "V_u2" };
    csl::InitSanitizer<complex_t> U_Wm1 { "U_Wm1" };
    csl::InitSanitizer<complex_t> U_Wm2 { "U_Wm2" };
    csl::InitSanitizer<complex_t> V_Wp1 { "V_Wp1" };
    csl::InitSanitizer<complex_t> V_Wp2 { "V_Wp2" };
    csl::InitSanitizer<complex_t> U_sd_00 { "U_sd_00" };
    csl::InitSanitizer<complex_t> U_sd_01 { "U_sd_01" };
    csl::InitSanitizer<complex_t> U_sd_02 { "U_sd_02" };
    csl::InitSanitizer<complex_t> U_sd_03 { "U_sd_03" };
    csl::InitSanitizer<complex_t> U_sd_04 { "U_sd_04" };
    csl::InitSanitizer<complex_t> U_sd_05 { "U_sd_05" };
    csl::InitSanitizer<complex_t> U_sd_10 { "U_sd_10" };
    csl::InitSanitizer<complex_t> U_sd_11 { "U_sd_11" };
    csl::InitSanitizer<complex_t> U_sd_12 { "U_sd_12" };
    csl::InitSanitizer<complex_t> U_sd_13 { "U_sd_13" };
    csl::InitSanitizer<complex_t> U_sd_14 { "U_sd_14" };
    csl::InitSanitizer<complex_t> U_sd_15 { "U_sd_15" };
    csl::InitSanitizer<complex_t> U_sd_20 { "U_sd_20" };
    csl::InitSanitizer<complex_t> U_sd_21 { "U_sd_21" };
    csl::InitSanitizer<complex_t> U_sd_22 { "U_sd_22" };
    csl::InitSanitizer<complex_t> U_sd_23 { "U_sd_23" };
    csl::InitSanitizer<complex_t> U_sd_24 { "U_sd_24" };
    csl::InitSanitizer<complex_t> U_sd_25 { "U_sd_25" };
    csl::InitSanitizer<complex_t> U_sd_30 { "U_sd_30" };
    csl::InitSanitizer<complex_t> U_sd_31 { "U_sd_31" };
    csl::InitSanitizer<complex_t> U_sd_32 { "U_sd_32" };
    csl::InitSanitizer<complex_t> U_sd_33 { "U_sd_33" };
    csl::InitSanitizer<complex_t> U_sd_34 { "U_sd_34" };
    csl::InitSanitizer<complex_t> U_sd_35 { "U_sd_35" };
    csl::InitSanitizer<complex_t> U_sd_40 { "U_sd_40" };
    csl::InitSanitizer<complex_t> U_sd_41 { "U_sd_41" };
    csl::InitSanitizer<complex_t> U_sd_42 { "U_sd_42" };
    csl::InitSanitizer<complex_t> U_sd_43 { "U_sd_43" };
    csl::InitSanitizer<complex_t> U_sd_44 { "U_sd_44" };
    csl::InitSanitizer<complex_t> U_sd_45 { "U_sd_45" };
    csl::InitSanitizer<complex_t> U_sd_50 { "U_sd_50" };
    csl::InitSanitizer<complex_t> U_sd_51 { "U_sd_51" };
    csl::InitSanitizer<complex_t> U_sd_52 { "U_sd_52" };
    csl::InitSanitizer<complex_t> U_sd_53 { "U_sd_53" };
    csl::InitSanitizer<complex_t> U_sd_54 { "U_sd_54" };
    csl::InitSanitizer<complex_t> U_sd_55 { "U_sd_55" };
    csl::InitSanitizer<complex_t> U_se_00 { "U_se_00" };
    csl::InitSanitizer<complex_t> U_se_01 { "U_se_01" };
    csl::InitSanitizer<complex_t> U_se_02 { "U_se_02" };
    csl::InitSanitizer<complex_t> U_se_03 { "U_se_03" };
    csl::InitSanitizer<complex_t> U_se_04 { "U_se_04" };
    csl::InitSanitizer<complex_t> U_se_05 { "U_se_05" };
    csl::InitSanitizer<complex_t> U_se_10 { "U_se_10" };
    csl::InitSanitizer<complex_t> U_se_11 { "U_se_11" };
    csl::InitSanitizer<complex_t> U_se_12 { "U_se_12" };
    csl::InitSanitizer<complex_t> U_se_13 { "U_se_13" };
    csl::InitSanitizer<complex_t> U_se_14 { "U_se_14" };
    csl::InitSanitizer<complex_t> U_se_15 { "U_se_15" };
    csl::InitSanitizer<complex_t> U_se_20 { "U_se_20" };
    csl::InitSanitizer<complex_t> U_se_21 { "U_se_21" };
    csl::InitSanitizer<complex_t> U_se_22 { "U_se_22" };
    csl::InitSanitizer<complex_t> U_se_23 { "U_se_23" };
    csl::InitSanitizer<complex_t> U_se_24 { "U_se_24" };
    csl::InitSanitizer<complex_t> U_se_25 { "U_se_25" };
    csl::InitSanitizer<complex_t> U_se_30 { "U_se_30" };
    csl::InitSanitizer<complex_t> U_se_31 { "U_se_31" };
    csl::InitSanitizer<complex_t> U_se_32 { "U_se_32" };
    csl::InitSanitizer<complex_t> U_se_33 { "U_se_33" };
    csl::InitSanitizer<complex_t> U_se_34 { "U_se_34" };
    csl::InitSanitizer<complex_t> U_se_35 { "U_se_35" };
    csl::InitSanitizer<complex_t> U_se_40 { "U_se_40" };
    csl::InitSanitizer<complex_t> U_se_41 { "U_se_41" };
    csl::InitSanitizer<complex_t> U_se_42 { "U_se_42" };
    csl::InitSanitizer<complex_t> U_se_43 { "U_se_43" };
    csl::InitSanitizer<complex_t> U_se_44 { "U_se_44" };
    csl::InitSanitizer<complex_t> U_se_45 { "U_se_45" };
    csl::InitSanitizer<complex_t> U_se_50 { "U_se_50" };
    csl::InitSanitizer<complex_t> U_se_51 { "U_se_51" };
    csl::InitSanitizer<complex_t> U_se_52 { "U_se_52" };
    csl::InitSanitizer<complex_t> U_se_53 { "U_se_53" };
    csl::InitSanitizer<complex_t> U_se_54 { "U_se_54" };
    csl::InitSanitizer<complex_t> U_se_55 { "U_se_55" };
    csl::InitSanitizer<complex_t> U_su_00 { "U_su_00" };
    csl::InitSanitizer<complex_t> U_su_01 { "U_su_01" };
    csl::InitSanitizer<complex_t> U_su_02 { "U_su_02" };
    csl::InitSanitizer<complex_t> U_su_03 { "U_su_03" };
    csl::InitSanitizer<complex_t> U_su_04 { "U_su_04" };
    csl::InitSanitizer<complex_t> U_su_05 { "U_su_05" };
    csl::InitSanitizer<complex_t> U_su_10 { "U_su_10" };
    csl::InitSanitizer<complex_t> U_su_11 { "U_su_11" };
    csl::InitSanitizer<complex_t> U_su_12 { "U_su_12" };
    csl::InitSanitizer<complex_t> U_su_13 { "U_su_13" };
    csl::InitSanitizer<complex_t> U_su_14 { "U_su_14" };
    csl::InitSanitizer<complex_t> U_su_15 { "U_su_15" };
    csl::InitSanitizer<complex_t> U_su_20 { "U_su_20" };
    csl::InitSanitizer<complex_t> U_su_21 { "U_su_21" };
    csl::InitSanitizer<complex_t> U_su_22 { "U_su_22" };
    csl::InitSanitizer<complex_t> U_su_23 { "U_su_23" };
    csl::InitSanitizer<complex_t> U_su_24 { "U_su_24" };
    csl::InitSanitizer<complex_t> U_su_25 { "U_su_25" };
    csl::InitSanitizer<complex_t> U_su_30 { "U_su_30" };
    csl::InitSanitizer<complex_t> U_su_31 { "U_su_31" };
    csl::InitSanitizer<complex_t> U_su_32 { "U_su_32" };
    csl::InitSanitizer<complex_t> U_su_33 { "U_su_33" };
    csl::InitSanitizer<complex_t> U_su_34 { "U_su_34" };
    csl::InitSanitizer<complex_t> U_su_35 { "U_su_35" };
    csl::InitSanitizer<complex_t> U_su_40 { "U_su_40" };
    csl::InitSanitizer<complex_t> U_su_41 { "U_su_41" };
    csl::InitSanitizer<complex_t> U_su_42 { "U_su_42" };
    csl::InitSanitizer<complex_t> U_su_43 { "U_su_43" };
    csl::InitSanitizer<complex_t> U_su_44 { "U_su_44" };
    csl::InitSanitizer<complex_t> U_su_45 { "U_su_45" };
    csl::InitSanitizer<complex_t> U_su_50 { "U_su_50" };
    csl::InitSanitizer<complex_t> U_su_51 { "U_su_51" };
    csl::InitSanitizer<complex_t> U_su_52 { "U_su_52" };
    csl::InitSanitizer<complex_t> U_su_53 { "U_su_53" };
    csl::InitSanitizer<complex_t> U_su_54 { "U_su_54" };
    csl::InitSanitizer<complex_t> U_su_55 { "U_su_55" };
    csl::InitSanitizer<complex_t> U_snu_00 { "U_snu_00" };
    csl::InitSanitizer<complex_t> U_snu_01 { "U_snu_01" };
    csl::InitSanitizer<complex_t> U_snu_02 { "U_snu_02" };
    csl::InitSanitizer<complex_t> U_snu_10 { "U_snu_10" };
    csl::InitSanitizer<complex_t> U_snu_11 { "U_snu_11" };
    csl::InitSanitizer<complex_t> U_snu_12 { "U_snu_12" };
    csl::InitSanitizer<complex_t> U_snu_20 { "U_snu_20" };
    csl::InitSanitizer<complex_t> U_snu_21 { "U_snu_21" };
    csl::InitSanitizer<complex_t> U_snu_22 { "U_snu_22" };

    void reset()
    {
        using real_params = std::array<csl::InitSanitizer<real_t>*, 60>;
        using complex_params = std::array<csl::InitSanitizer<complex_t>*, 177>;

        for (auto &par : real_params{
                &M_3, &g_s, &V_cb, &V_tb, &V_us, 
                &beta, &e_em, &s_12, &s_34, &alpha, &Finite, 
                &theta_W, &V_ub_mod, &reg_prop, &delta_wolf, &M_W, &M_Z, 
                &m_b, &m_c, &m_d, &m_e, &m_s, &m_t, 
                &m_u, &m_A0, &m_G0, &m_Gp, &m_H0, &m_Hp, 
                &m_h0, &m_mu, &m_sG, &m_C1p, &m_C2p, &m_N_1, 
                &m_N_2, &m_N_3, &m_N_4, &m_tau, &m_sb_L, &m_sb_R, 
                &m_sc_L, &m_sc_R, &m_sd_L, &m_sd_R, &m_se_L, &m_se_R, 
                &m_ss_L, &m_ss_R, &m_st_L, &m_st_R, &m_su_L, &m_su_R, 
                &m_smu_L, &m_smu_R, &m_snu_e, &m_snu_mu, &m_stau_L, &m_stau_R, 
                &m_snu_tau, })
        {
            par->reset();
        }

        for (auto &par : complex_params{
                &A_b, &A_t, &M_1, &M_2, &M_eL, 
                &M_eR, &V_cd, &V_cs, &V_td, &V_ts, &V_ud, 
                &mu_h, &A_tau, &M_q1L, &M_q3L, &M_qbR, &M_qdR, 
                &M_qtR, &M_quR, &M_tauL, &M_tauR, &del_DR_1, &del_DR_2, 
                &del_ER_1, &del_ER_2, &del_LL_1, &del_LL_2, &del_QL_1, &del_QL_2, 
                &del_UR_1, &del_UR_2, &del_DR_12, &del_ER_12, &del_LL_12, &del_QL_12, 
                &del_UR_12, &N_B1, &N_B2, &N_B3, &N_B4, &N_W1, 
                &N_W2, &N_W3, &N_W4, &N_d1, &N_d2, &N_d3, 
                &N_d4, &N_u1, &N_u2, &N_u3, &N_u4, &U_d1, 
                &U_d2, &V_u1, &V_u2, &U_Wm1, &U_Wm2, &V_Wp1, 
                &V_Wp2, &U_sd_00, &U_sd_01, &U_sd_02, &U_sd_03, &U_sd_04, 
                &U_sd_05, &U_sd_10, &U_sd_11, &U_sd_12, &U_sd_13, &U_sd_14, 
                &U_sd_15, &U_sd_20, &U_sd_21, &U_sd_22, &U_sd_23, &U_sd_24, 
                &U_sd_25, &U_sd_30, &U_sd_31, &U_sd_32, &U_sd_33, &U_sd_34, 
                &U_sd_35, &U_sd_40, &U_sd_41, &U_sd_42, &U_sd_43, &U_sd_44, 
                &U_sd_45, &U_sd_50, &U_sd_51, &U_sd_52, &U_sd_53, &U_sd_54, 
                &U_sd_55, &U_se_00, &U_se_01, &U_se_02, &U_se_03, &U_se_04, 
                &U_se_05, &U_se_10, &U_se_11, &U_se_12, &U_se_13, &U_se_14, 
                &U_se_15, &U_se_20, &U_se_21, &U_se_22, &U_se_23, &U_se_24, 
                &U_se_25, &U_se_30, &U_se_31, &U_se_32, &U_se_33, &U_se_34, 
                &U_se_35, &U_se_40, &U_se_41, &U_se_42, &U_se_43, &U_se_44, 
                &U_se_45, &U_se_50, &U_se_51, &U_se_52, &U_se_53, &U_se_54, 
                &U_se_55, &U_su_00, &U_su_01, &U_su_02, &U_su_03, &U_su_04, 
                &U_su_05, &U_su_10, &U_su_11, &U_su_12, &U_su_13, &U_su_14, 
                &U_su_15, &U_su_20, &U_su_21, &U_su_22, &U_su_23, &U_su_24, 
                &U_su_25, &U_su_30, &U_su_31, &U_su_32, &U_su_33, &U_su_34, 
                &U_su_35, &U_su_40, &U_su_41, &U_su_42, &U_su_43, &U_su_44, 
                &U_su_45, &U_su_50, &U_su_51, &U_su_52, &U_su_53, &U_su_54, 
                &U_su_55, &U_snu_00, &U_snu_01, &U_snu_02, &U_snu_10, &U_snu_11, 
                &U_snu_12, &U_snu_20, &U_snu_21, &U_snu_22, })
        {
            par->reset();
        }
    }

    void print(std::ostream &out = std::cout)
    {
        using real_params = std::array<csl::InitSanitizer<real_t> const*, 60>;
        using complex_params = std::array<csl::InitSanitizer<complex_t> const*, 177>;

        out << "param_t struct:\n";
        out << "Real parameters\n";
        for (auto const &par : real_params{
                &M_3, &g_s, &V_cb, &V_tb, &V_us, 
                &beta, &e_em, &s_12, &s_34, &alpha, &Finite, 
                &theta_W, &V_ub_mod, &reg_prop, &delta_wolf, &M_W, &M_Z, 
                &m_b, &m_c, &m_d, &m_e, &m_s, &m_t, 
                &m_u, &m_A0, &m_G0, &m_Gp, &m_H0, &m_Hp, 
                &m_h0, &m_mu, &m_sG, &m_C1p, &m_C2p, &m_N_1, 
                &m_N_2, &m_N_3, &m_N_4, &m_tau, &m_sb_L, &m_sb_R, 
                &m_sc_L, &m_sc_R, &m_sd_L, &m_sd_R, &m_se_L, &m_se_R, 
                &m_ss_L, &m_ss_R, &m_st_L, &m_st_R, &m_su_L, &m_su_R, 
                &m_smu_L, &m_smu_R, &m_snu_e, &m_snu_mu, &m_stau_L, &m_stau_R, 
                &m_snu_tau, })
        {
            out << "  -> ";
            par->print(out);
        }

        out << "Complex parameters\n";
        for (auto const &par : complex_params{
                &A_b, &A_t, &M_1, &M_2, &M_eL, 
                &M_eR, &V_cd, &V_cs, &V_td, &V_ts, &V_ud, 
                &mu_h, &A_tau, &M_q1L, &M_q3L, &M_qbR, &M_qdR, 
                &M_qtR, &M_quR, &M_tauL, &M_tauR, &del_DR_1, &del_DR_2, 
                &del_ER_1, &del_ER_2, &del_LL_1, &del_LL_2, &del_QL_1, &del_QL_2, 
                &del_UR_1, &del_UR_2, &del_DR_12, &del_ER_12, &del_LL_12, &del_QL_12, 
                &del_UR_12, &N_B1, &N_B2, &N_B3, &N_B4, &N_W1, 
                &N_W2, &N_W3, &N_W4, &N_d1, &N_d2, &N_d3, 
                &N_d4, &N_u1, &N_u2, &N_u3, &N_u4, &U_d1, 
                &U_d2, &V_u1, &V_u2, &U_Wm1, &U_Wm2, &V_Wp1, 
                &V_Wp2, &U_sd_00, &U_sd_01, &U_sd_02, &U_sd_03, &U_sd_04, 
                &U_sd_05, &U_sd_10, &U_sd_11, &U_sd_12, &U_sd_13, &U_sd_14, 
                &U_sd_15, &U_sd_20, &U_sd_21, &U_sd_22, &U_sd_23, &U_sd_24, 
                &U_sd_25, &U_sd_30, &U_sd_31, &U_sd_32, &U_sd_33, &U_sd_34, 
                &U_sd_35, &U_sd_40, &U_sd_41, &U_sd_42, &U_sd_43, &U_sd_44, 
                &U_sd_45, &U_sd_50, &U_sd_51, &U_sd_52, &U_sd_53, &U_sd_54, 
                &U_sd_55, &U_se_00, &U_se_01, &U_se_02, &U_se_03, &U_se_04, 
                &U_se_05, &U_se_10, &U_se_11, &U_se_12, &U_se_13, &U_se_14, 
                &U_se_15, &U_se_20, &U_se_21, &U_se_22, &U_se_23, &U_se_24, 
                &U_se_25, &U_se_30, &U_se_31, &U_se_32, &U_se_33, &U_se_34, 
                &U_se_35, &U_se_40, &U_se_41, &U_se_42, &U_se_43, &U_se_44, 
                &U_se_45, &U_se_50, &U_se_51, &U_se_52, &U_se_53, &U_se_54, 
                &U_se_55, &U_su_00, &U_su_01, &U_su_02, &U_su_03, &U_su_04, 
                &U_su_05, &U_su_10, &U_su_11, &U_su_12, &U_su_13, &U_su_14, 
                &U_su_15, &U_su_20, &U_su_21, &U_su_22, &U_su_23, &U_su_24, 
                &U_su_25, &U_su_30, &U_su_31, &U_su_32, &U_su_33, &U_su_34, 
                &U_su_35, &U_su_40, &U_su_41, &U_su_42, &U_su_43, &U_su_44, 
                &U_su_45, &U_su_50, &U_su_51, &U_su_52, &U_su_53, &U_su_54, 
                &U_su_55, &U_snu_00, &U_snu_01, &U_snu_02, &U_snu_10, &U_snu_11, 
                &U_snu_12, &U_snu_20, &U_snu_21, &U_snu_22, })
        {
            out << "  -> ";
            par->print(out);
        }
        out << "\n";
    }

    std::map<std::string, csl::InitSanitizer<real_t>*> realParams {
        {"M_3", &M_3},
        {"g_s", &g_s},
        {"V_cb", &V_cb},
        {"V_tb", &V_tb},
        {"V_us", &V_us},
        {"beta", &beta},
        {"e_em", &e_em},
        {"s_12", &s_12},
        {"s_34", &s_34},
        {"alpha", &alpha},
        {"Finite", &Finite},
        {"theta_W", &theta_W},
        {"V_ub_mod", &V_ub_mod},
        {"reg_prop", &reg_prop},
        {"delta_wolf", &delta_wolf},
        {"M_W", &M_W},
        {"M_Z", &M_Z},
        {"m_b", &m_b},
        {"m_c", &m_c},
        {"m_d", &m_d},
        {"m_e", &m_e},
        {"m_s", &m_s},
        {"m_t", &m_t},
        {"m_u", &m_u},
        {"m_A0", &m_A0},
        {"m_G0", &m_G0},
        {"m_Gp", &m_Gp},
        {"m_H0", &m_H0},
        {"m_Hp", &m_Hp},
        {"m_h0", &m_h0},
        {"m_mu", &m_mu},
        {"m_sG", &m_sG},
        {"m_C1p", &m_C1p},
        {"m_C2p", &m_C2p},
        {"m_N_1", &m_N_1},
        {"m_N_2", &m_N_2},
        {"m_N_3", &m_N_3},
        {"m_N_4", &m_N_4},
        {"m_tau", &m_tau},
        {"m_sb_L", &m_sb_L},
        {"m_sb_R", &m_sb_R},
        {"m_sc_L", &m_sc_L},
        {"m_sc_R", &m_sc_R},
        {"m_sd_L", &m_sd_L},
        {"m_sd_R", &m_sd_R},
        {"m_se_L", &m_se_L},
        {"m_se_R", &m_se_R},
        {"m_ss_L", &m_ss_L},
        {"m_ss_R", &m_ss_R},
        {"m_st_L", &m_st_L},
        {"m_st_R", &m_st_R},
        {"m_su_L", &m_su_L},
        {"m_su_R", &m_su_R},
        {"m_smu_L", &m_smu_L},
        {"m_smu_R", &m_smu_R},
        {"m_snu_e", &m_snu_e},
        {"m_snu_mu", &m_snu_mu},
        {"m_stau_L", &m_stau_L},
        {"m_stau_R", &m_stau_R},
        {"m_snu_tau", &m_snu_tau},
    };

    std::map<std::string, csl::InitSanitizer<complex_t>*> complexParams {
        {"A_b", &A_b},
        {"A_t", &A_t},
        {"M_1", &M_1},
        {"M_2", &M_2},
        {"M_eL", &M_eL},
        {"M_eR", &M_eR},
        {"V_cd", &V_cd},
        {"V_cs", &V_cs},
        {"V_td", &V_td},
        {"V_ts", &V_ts},
        {"V_ud", &V_ud},
        {"mu_h", &mu_h},
        {"A_tau", &A_tau},
        {"M_q1L", &M_q1L},
        {"M_q3L", &M_q3L},
        {"M_qbR", &M_qbR},
        {"M_qdR", &M_qdR},
        {"M_qtR", &M_qtR},
        {"M_quR", &M_quR},
        {"M_tauL", &M_tauL},
        {"M_tauR", &M_tauR},
        {"del_DR_1", &del_DR_1},
        {"del_DR_2", &del_DR_2},
        {"del_ER_1", &del_ER_1},
        {"del_ER_2", &del_ER_2},
        {"del_LL_1", &del_LL_1},
        {"del_LL_2", &del_LL_2},
        {"del_QL_1", &del_QL_1},
        {"del_QL_2", &del_QL_2},
        {"del_UR_1", &del_UR_1},
        {"del_UR_2", &del_UR_2},
        {"del_DR_12", &del_DR_12},
        {"del_ER_12", &del_ER_12},
        {"del_LL_12", &del_LL_12},
        {"del_QL_12", &del_QL_12},
        {"del_UR_12", &del_UR_12},
        {"N_B1", &N_B1},
        {"N_B2", &N_B2},
        {"N_B3", &N_B3},
        {"N_B4", &N_B4},
        {"N_W1", &N_W1},
        {"N_W2", &N_W2},
        {"N_W3", &N_W3},
        {"N_W4", &N_W4},
        {"N_d1", &N_d1},
        {"N_d2", &N_d2},
        {"N_d3", &N_d3},
        {"N_d4", &N_d4},
        {"N_u1", &N_u1},
        {"N_u2", &N_u2},
        {"N_u3", &N_u3},
        {"N_u4", &N_u4},
        {"U_d1", &U_d1},
        {"U_d2", &U_d2},
        {"V_u1", &V_u1},
        {"V_u2", &V_u2},
        {"U_Wm1", &U_Wm1},
        {"U_Wm2", &U_Wm2},
        {"V_Wp1", &V_Wp1},
        {"V_Wp2", &V_Wp2},
        {"U_sd_00", &U_sd_00},
        {"U_sd_01", &U_sd_01},
        {"U_sd_02", &U_sd_02},
        {"U_sd_03", &U_sd_03},
        {"U_sd_04", &U_sd_04},
        {"U_sd_05", &U_sd_05},
        {"U_sd_10", &U_sd_10},
        {"U_sd_11", &U_sd_11},
        {"U_sd_12", &U_sd_12},
        {"U_sd_13", &U_sd_13},
        {"U_sd_14", &U_sd_14},
        {"U_sd_15", &U_sd_15},
        {"U_sd_20", &U_sd_20},
        {"U_sd_21", &U_sd_21},
        {"U_sd_22", &U_sd_22},
        {"U_sd_23", &U_sd_23},
        {"U_sd_24", &U_sd_24},
        {"U_sd_25", &U_sd_25},
        {"U_sd_30", &U_sd_30},
        {"U_sd_31", &U_sd_31},
        {"U_sd_32", &U_sd_32},
        {"U_sd_33", &U_sd_33},
        {"U_sd_34", &U_sd_34},
        {"U_sd_35", &U_sd_35},
        {"U_sd_40", &U_sd_40},
        {"U_sd_41", &U_sd_41},
        {"U_sd_42", &U_sd_42},
        {"U_sd_43", &U_sd_43},
        {"U_sd_44", &U_sd_44},
        {"U_sd_45", &U_sd_45},
        {"U_sd_50", &U_sd_50},
        {"U_sd_51", &U_sd_51},
        {"U_sd_52", &U_sd_52},
        {"U_sd_53", &U_sd_53},
        {"U_sd_54", &U_sd_54},
        {"U_sd_55", &U_sd_55},
        {"U_se_00", &U_se_00},
        {"U_se_01", &U_se_01},
        {"U_se_02", &U_se_02},
        {"U_se_03", &U_se_03},
        {"U_se_04", &U_se_04},
        {"U_se_05", &U_se_05},
        {"U_se_10", &U_se_10},
        {"U_se_11", &U_se_11},
        {"U_se_12", &U_se_12},
        {"U_se_13", &U_se_13},
        {"U_se_14", &U_se_14},
        {"U_se_15", &U_se_15},
        {"U_se_20", &U_se_20},
        {"U_se_21", &U_se_21},
        {"U_se_22", &U_se_22},
        {"U_se_23", &U_se_23},
        {"U_se_24", &U_se_24},
        {"U_se_25", &U_se_25},
        {"U_se_30", &U_se_30},
        {"U_se_31", &U_se_31},
        {"U_se_32", &U_se_32},
        {"U_se_33", &U_se_33},
        {"U_se_34", &U_se_34},
        {"U_se_35", &U_se_35},
        {"U_se_40", &U_se_40},
        {"U_se_41", &U_se_41},
        {"U_se_42", &U_se_42},
        {"U_se_43", &U_se_43},
        {"U_se_44", &U_se_44},
        {"U_se_45", &U_se_45},
        {"U_se_50", &U_se_50},
        {"U_se_51", &U_se_51},
        {"U_se_52", &U_se_52},
        {"U_se_53", &U_se_53},
        {"U_se_54", &U_se_54},
        {"U_se_55", &U_se_55},
        {"U_su_00", &U_su_00},
        {"U_su_01", &U_su_01},
        {"U_su_02", &U_su_02},
        {"U_su_03", &U_su_03},
        {"U_su_04", &U_su_04},
        {"U_su_05", &U_su_05},
        {"U_su_10", &U_su_10},
        {"U_su_11", &U_su_11},
        {"U_su_12", &U_su_12},
        {"U_su_13", &U_su_13},
        {"U_su_14", &U_su_14},
        {"U_su_15", &U_su_15},
        {"U_su_20", &U_su_20},
        {"U_su_21", &U_su_21},
        {"U_su_22", &U_su_22},
        {"U_su_23", &U_su_23},
        {"U_su_24", &U_su_24},
        {"U_su_25", &U_su_25},
        {"U_su_30", &U_su_30},
        {"U_su_31", &U_su_31},
        {"U_su_32", &U_su_32},
        {"U_su_33", &U_su_33},
        {"U_su_34", &U_su_34},
        {"U_su_35", &U_su_35},
        {"U_su_40", &U_su_40},
        {"U_su_41", &U_su_41},
        {"U_su_42", &U_su_42},
        {"U_su_43", &U_su_43},
        {"U_su_44", &U_su_44},
        {"U_su_45", &U_su_45},
        {"U_su_50", &U_su_50},
        {"U_su_51", &U_su_51},
        {"U_su_52", &U_su_52},
        {"U_su_53", &U_su_53},
        {"U_su_54", &U_su_54},
        {"U_su_55", &U_su_55},
        {"U_snu_00", &U_snu_00},
        {"U_snu_01", &U_snu_01},
        {"U_snu_02", &U_snu_02},
        {"U_snu_10", &U_snu_10},
        {"U_snu_11", &U_snu_11},
        {"U_snu_12", &U_snu_12},
        {"U_snu_20", &U_snu_20},
        {"U_snu_21", &U_snu_21},
        {"U_snu_22", &U_snu_22},
    };

};


}

#endif
