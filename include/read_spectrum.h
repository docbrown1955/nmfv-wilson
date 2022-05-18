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

#pragma once

#include "params.h"
#include "global.h"
#include "lhaData.h"

inline void readSpectrum(
        c9_nmfv::param_t         &param,
        mty::lha::LHAFileData const &data
        )
{
//    param.e_em = std::sqrt(4*M_PI * data.getValue("SMINPUTS", 1).value());
    param.e_em = 1.27931417E+02 ; 
    param.M_W = 80.385;
//    param.M_Z = data.getValue("SMINPUTS", 4).value(); 
    param.M_Z = 9.11876000E+01;
    param.m_b = data.getValue("SMINPUTS", 5).value();
    param.m_s = 0.150;
    param.m_t = data.getValue("SMINPUTS", 6).value();
    param.m_c = data.getValue("SMINPUTS",24).value();
    param.m_u = data.getValue("SMINPUTS",22).value();
    param.m_tau = data.getValue("SMINPUTS",7).value();
    param.m_u = data.getValue("SMINPUTS",22).value();

//    param.m_mu = 0; // 0.1;
    param.theta_W = crealq(cacosq(param.M_W / param.M_Z));
    param.alpha = data.getValue("alpha",0).value();
    param.beta = std::atan(data.getValue("MINPAR", 3).value());
    param.mu_h = data.getValue("EXTPAR", 23).value();
    param.M_1 = data.getValue("EXTPAR", 1).value();
    param.M_2 = data.getValue("EXTPAR", 2).value();
    param.M_3 = data.getValue("EXTPAR", 3).value();
    param.A_b = data.getValue("Td", 3,3).value()*data.getValue("Yd",3,3).value();
    param.A_t = data.getValue("Tu", 3,3).value()*data.getValue("Yu",3,3).value();
    param.A_tau = data.getValue("Te", 3,3).value()*data.getValue("Ye",3,3).value();
    param.m_A0 = data.getValue("EXTPAR", 26).value();
    // Block Gauge !
    param.g_s = data.getValue("gauge",3).value();
    // Block MaSS
    param.m_sG = data.getValue("MASS",1000021).value();
    param.m_Hp = data.getValue("MASS",37).value();
    param.m_H0 = data.getValue("MASS",35).value();
    param.m_su_L = data.getValue("MASS", 1000001).value();
    param.m_su_R = data.getValue("MASS", 2000001).value();
    param.m_sc_L = data.getValue("MASS", 1000003).value();
    param.m_sc_R = data.getValue("MASS", 2000003).value();
    param.m_st_L = data.getValue("MASS", 1000005).value();
    param.m_st_R = data.getValue("MASS", 2000005).value();
   
    param.m_sd_L = data.getValue("MASS", 1000002).value();
    param.m_sd_R = data.getValue("MASS", 2000002).value();
    param.m_ss_L = data.getValue("MASS", 1000004).value();
    param.m_ss_R = data.getValue("MASS", 2000004).value();
    param.m_sb_L = data.getValue("MASS", 1000006).value();
    param.m_sb_R = data.getValue("MASS", 2000006).value();
    
    param.m_snu_e   = data.getValue("MASS", 1000012).value();
    param.m_snu_mu  = data.getValue("MASS", 1000014).value();
    param.m_snu_tau = data.getValue("MASS", 1000016).value();

    param.m_se_L   = data.getValue("MASS", 1000011).value();
    param.m_smu_L  = data.getValue("MASS", 1000013).value();
    param.m_stau_L = data.getValue("MASS", 1000015).value();
    param.m_se_R   = data.getValue("MASS", 2000011).value();
    param.m_smu_R  = data.getValue("MASS", 2000013).value();
    param.m_stau_R = data.getValue("MASS", 2000015).value();

    param.m_C1p = data.getValue("MASS", 1000024).value();
    param.m_C2p = data.getValue("MASS", 1000037).value();
    
    param.m_N_1 = data.getValue("MASS", 1000022).value();
    param.m_N_2 = data.getValue("MASS", 1000023).value();
    param.m_N_3 = data.getValue("MASS", 1000025).value();
    param.m_N_4 = data.getValue("MASS", 1000035).value();


    param.Finite = 1;
    param.reg_prop =1e-3;
    // param.reg_int = 1e-2;


    param.U_Wm1 = data.getValue("Umix", 1, 1).value();
    param.U_d1  = data.getValue("Umix", 1, 2).value();
    param.U_Wm2 = data.getValue("Umix", 2, 1).value();
    param.U_d2  = data.getValue("Umix", 2, 2).value();

    param.V_Wp1 = data.getValue("Vmix", 1, 1).value();
    param.V_u1  = data.getValue("Vmix", 1, 2).value();
    param.V_Wp2 = data.getValue("Vmix", 2, 1).value();
    param.V_u2  = data.getValue("Vmix", 2, 2).value();

 //   param.U_st_00 = data.getValue("stopmix", 1, 1).value();
//    param.U_st_01 = data.getValue("stopmix", 2, 1).value();
//    param.U_st_10 = data.getValue("stopmix", 1, 2).value();
//    param.U_st_11 = data.getValue("stopmix", 2, 2).value();
//  USQmix : 
    param.U_su_00 = data.getValue("USQmix",1,1).value();
    param.U_su_01 = data.getValue("USQmix",2,1).value();
    param.U_su_02 = data.getValue("USQmix",3,1).value();
    param.U_su_03 = data.getValue("USQmix",4,1).value();
    param.U_su_04 = data.getValue("USQmix",5,1).value();
    param.U_su_05 = data.getValue("USQmix",6,1).value();
    param.U_su_10 = data.getValue("USQmix",1,2).value();
    param.U_su_11 = data.getValue("USQmix",2,2).value();
    param.U_su_12 = data.getValue("USQmix",3,2).value();
    param.U_su_13 = data.getValue("USQmix",4,2).value();
    param.U_su_14 = data.getValue("USQmix",5,2).value();
    param.U_su_15 = data.getValue("USQmix",6,2).value();
    param.U_su_20 = data.getValue("USQmix",1,3).value();
    param.U_su_21 = data.getValue("USQmix",2,3).value();
    param.U_su_22 = data.getValue("USQmix",3,3).value();
    param.U_su_23 = data.getValue("USQmix",4,3).value();
    param.U_su_24 = data.getValue("USQmix",5,3).value();
    param.U_su_25 = data.getValue("USQmix",6,3).value();
    param.U_su_30 = data.getValue("USQmix",1,4).value();
    param.U_su_31 = data.getValue("USQmix",2,4).value();
    param.U_su_32 = data.getValue("USQmix",3,4).value();
    param.U_su_33 = data.getValue("USQmix",4,4).value();
    param.U_su_34 = data.getValue("USQmix",5,4).value();
    param.U_su_35 = data.getValue("USQmix",6,4).value();
    param.U_su_40 = data.getValue("USQmix",1,5).value();
    param.U_su_41 = data.getValue("USQmix",2,5).value();
    param.U_su_42 = data.getValue("USQmix",3,5).value();
    param.U_su_43 = data.getValue("USQmix",4,5).value();
    param.U_su_44 = data.getValue("USQmix",5,5).value();
    param.U_su_45 = data.getValue("USQmix",6,5).value();
    param.U_su_50 = data.getValue("USQmix",1,6).value();
    param.U_su_51 = data.getValue("USQmix",2,6).value();
    param.U_su_52 = data.getValue("USQmix",3,6).value();
    param.U_su_53 = data.getValue("USQmix",4,6).value();
    param.U_su_54 = data.getValue("USQmix",5,6).value();
    param.U_su_55 = data.getValue("USQmix",6,6).value();
    // DSQmix: 
   param.U_sd_00 = data.getValue("DSQmix",1,1).value();
   param.U_sd_01 = data.getValue("DSQmix",2,1).value();
   param.U_sd_02 = data.getValue("DSQmix",3,1).value();
   param.U_sd_03 = data.getValue("DSQmix",4,1).value();
   param.U_sd_04 = data.getValue("DSQmix",5,1).value();
   param.U_sd_05 = data.getValue("DSQmix",6,1).value();
   param.U_sd_10 = data.getValue("DSQmix",1,2).value();
   param.U_sd_11 = data.getValue("DSQmix",2,2).value();
   param.U_sd_12 = data.getValue("DSQmix",3,2).value();
   param.U_sd_13 = data.getValue("DSQmix",4,2).value();
   param.U_sd_14 = data.getValue("DSQmix",5,2).value();
   param.U_sd_15 = data.getValue("DSQmix",6,2).value();
   param.U_sd_20 = data.getValue("DSQmix",1,3).value();
   param.U_sd_21 = data.getValue("DSQmix",2,3).value();
   param.U_sd_22 = data.getValue("DSQmix",3,3).value();
   param.U_sd_23 = data.getValue("DSQmix",4,3).value();
   param.U_sd_24 = data.getValue("DSQmix",5,3).value();
   param.U_sd_25 = data.getValue("DSQmix",6,3).value();
   param.U_sd_30 = data.getValue("DSQmix",1,4).value();
   param.U_sd_31 = data.getValue("DSQmix",2,4).value();
   param.U_sd_32 = data.getValue("DSQmix",3,4).value();
   param.U_sd_33 = data.getValue("DSQmix",4,4).value();
   param.U_sd_34 = data.getValue("DSQmix",5,4).value();
   param.U_sd_35 = data.getValue("DSQmix",6,4).value();
   param.U_sd_40 = data.getValue("DSQmix",1,5).value();
   param.U_sd_41 = data.getValue("DSQmix",2,5).value();
   param.U_sd_42 = data.getValue("DSQmix",3,5).value();
   param.U_sd_43 = data.getValue("DSQmix",4,5).value();
   param.U_sd_44 = data.getValue("DSQmix",5,5).value();
   param.U_sd_45 = data.getValue("DSQmix",6,5).value();
   param.U_sd_50 = data.getValue("DSQmix",1,6).value();
   param.U_sd_51 = data.getValue("DSQmix",2,6).value();
   param.U_sd_52 = data.getValue("DSQmix",3,6).value();
   param.U_sd_53 = data.getValue("DSQmix",4,6).value();
   param.U_sd_54 = data.getValue("DSQmix",5,6).value();
   param.U_sd_55 = data.getValue("DSQmix",6,6).value();
   //SELmix
   //
   param.U_se_00 = data.getValue("SELmix",1,1).value();
   param.U_se_01 = data.getValue("SELmix",2,1).value();
   param.U_se_02 = data.getValue("SELmix",3,1).value();
   param.U_se_03 = data.getValue("SELmix",4,1).value();
   param.U_se_04 = data.getValue("SELmix",5,1).value();
   param.U_se_05 = data.getValue("SELmix",6,1).value();
   param.U_se_10 = data.getValue("SELmix",1,2).value();
   param.U_se_11 = data.getValue("SELmix",2,2).value();
   param.U_se_12 = data.getValue("SELmix",3,2).value();
   param.U_se_13 = data.getValue("SELmix",4,2).value();
   param.U_se_14 = data.getValue("SELmix",5,2).value();
   param.U_se_15 = data.getValue("SELmix",6,2).value();
   param.U_se_20 = data.getValue("SELmix",1,3).value();
   param.U_se_21 = data.getValue("SELmix",2,3).value();
   param.U_se_22 = data.getValue("SELmix",3,3).value();
   param.U_se_23 = data.getValue("SELmix",4,3).value();
   param.U_se_24 = data.getValue("SELmix",5,3).value();
   param.U_se_25 = data.getValue("SELmix",6,3).value();
   param.U_se_30 = data.getValue("SELmix",1,4).value();
   param.U_se_31 = data.getValue("SELmix",2,4).value();
   param.U_se_32 = data.getValue("SELmix",3,4).value();
   param.U_se_33 = data.getValue("SELmix",4,4).value();
   param.U_se_34 = data.getValue("SELmix",5,4).value();
   param.U_se_35 = data.getValue("SELmix",6,4).value();
   param.U_se_40 = data.getValue("SELmix",1,5).value();
   param.U_se_41 = data.getValue("SELmix",2,5).value();
   param.U_se_42 = data.getValue("SELmix",3,5).value();
   param.U_se_43 = data.getValue("SELmix",4,5).value();
   param.U_se_44 = data.getValue("SELmix",5,5).value();
   param.U_se_45 = data.getValue("SELmix",6,5).value();
   param.U_se_50 = data.getValue("SELmix",1,6).value();
   param.U_se_51 = data.getValue("SELmix",2,6).value();
   param.U_se_52 = data.getValue("SELmix",3,6).value();
   param.U_se_53 = data.getValue("SELmix",4,6).value();
   param.U_se_54 = data.getValue("SELmix",5,6).value();
   param.U_se_55 = data.getValue("SELmix",6,6).value();
  // SNUmix :
   param.U_snu_00 = data.getValue("SNUmix",1,1).value();
   param.U_snu_01 = data.getValue("SNUmix",2,1).value();
   param.U_snu_02 = data.getValue("SNUmix",3,1).value();
   param.U_snu_10 = data.getValue("SNUmix",1,2).value();
   param.U_snu_11 = data.getValue("SNUmix",2,2).value();
   param.U_snu_12 = data.getValue("SNUmix",3,2).value();
   param.U_snu_20 = data.getValue("SNUmix",1,3).value();
   param.U_snu_21 = data.getValue("SNUmix",2,3).value();
   param.U_snu_22 = data.getValue("SNUmix",3,3).value();
   param.N_B1 = data.getValue("Nmix",1,1).value(); 
   param.N_B2 = data.getValue("Nmix",2,1).value(); 
   param.N_B3 = data.getValue("Nmix",3,1).value(); 
   param.N_B4 = data.getValue("Nmix",4,1).value(); 
   param.N_W1 = data.getValue("Nmix",1,2).value(); 
   param.N_W2 = data.getValue("Nmix",2,2).value(); 
   param.N_W3 = data.getValue("Nmix",3,2).value(); 
   param.N_W4 = data.getValue("Nmix",4,2).value(); 
   param.N_d1 = data.getValue("Nmix",1,3).value(); 
   param.N_d2 = data.getValue("Nmix",2,3).value(); 
   param.N_d3 = data.getValue("Nmix",3,3).value(); 
   param.N_d4 = data.getValue("Nmix",4,3).value(); 
   param.N_u1 = data.getValue("Nmix",1,4).value(); 
   param.N_u2 = data.getValue("Nmix",2,4).value(); 
   param.N_u3 = data.getValue("Nmix",3,4).value(); 
   param.N_u4 = data.getValue("Nmix",4,4).value(); 

//    updateSpectrum(param);
}

bool comp(
        std::string const &name,
        double a,
        double b
        )
{
    if (b == 0) {
        if (a != 0) {
            std::cout << name << " not good:" << '\n';
            std::cout << a << " != " << b << '\n';
            return false;
        }
    }
    else if (std::abs((b - a) / b) > 0.001) {
        std::cout << name << " has a discrepancy:" << '\n';
        std::cout << a << " != " << b << '\n';
        return false;
    }
    return true;
}

//inline void readSpectrum2(
//        c9_nmfv::param_t         &param,
//        mty::lha::LHAFileData const &data
//        )
//{
//    double MZ = data.getValue("ADDIT", 0).value();
//    // double MW = data.getValue("ADDIT", 1).value();
//    // double gL = data.getValue("ADDIT", 2).value();
//    // double gY = data.getValue("ADDIT", 3).value();
//    double sw2 = data.getValue("ADDIT", 4).value();
//    double m_st_L = data.getValue("ADDIT", 5).value();
//    double m_st_R = data.getValue("ADDIT", 6).value();
//
//    double m_C1p = data.getValue("ADDIT", 7).value();
//    double m_C2p = data.getValue("ADDIT", 8).value();
//    double U_Wm1  = data.getValue("ADDIT", 9).value();
//    double U_d1  = data.getValue("ADDIT", 10).value();
//    double U_Wm2  = data.getValue("ADDIT", 11).value();
//    double U_d2  = data.getValue("ADDIT", 12).value();
//
//    double V_Wp1  = data.getValue("ADDIT", 13).value();
//    double V_u1  = data.getValue("ADDIT", 14).value();
//    double V_Wp2  = data.getValue("ADDIT", 15).value();
//    double V_u2  = data.getValue("ADDIT", 16).value();
//    double U_st_00  = data.getValue("ADDIT", 17).value();
//    double U_st_10  = data.getValue("ADDIT", 18).value();
//    double U_st_01  = data.getValue("ADDIT", 19).value();
//
//    double U_st_11  = data.getValue("ADDIT", 20).value();
//    double m_u1 = data.getValue("ADDIT", 21).value();
//    double m_u2 = data.getValue("ADDIT", 22).value();
//    double m_sc_L = data.getValue("ADDIT", 23).value();
//    double m_sc_R = data.getValue("ADDIT", 24).value();
//    double tanb = data.getValue("ADDIT", 25).value();
//    // double GF = data.getValue("SMINPUTS", 2).value();
//    double alpha_em = data.getValue("SMINPUTS", 1).value();
//    double e_em = std::sqrt(4*M_PI/alpha_em);
//    double thetaW = std::asin(std::sqrt(sw2));
//
//
//    param.M_A = data.getValue("EXTPAR", 26).value();
//    param.e_em = e_em;
//    //param.M_W = MW;
//    param.M_W = MZ * std::cos(thetaW);
//    param.M_Z = MZ;
//    param.m_b = data.getValue("SMINPUTS", 5).value();
//    param.m_s = 0.150;
//    param.m_t = data.getValue("SMINPUTS", 6).value();
//    param.m_mu = 0; // 0.1;
//    param.theta_W = thetaW;
//
//    param.beta = std::atan(tanb);
//
//    param.Finite = 1;
//    param.reg_prop = 0;
//    // param.reg_int = 1e-2;
//
//    param.m_snu_e = data.getValue("MASS", 1000012).value();
//    // param.m_snu_mu = param.m_snu_e;
//
//    param.m_su_L = m_u1;
//    param.m_su_R = m_u2;
//
//    param.m_sc_L = m_sc_L;
//    param.m_sc_R = m_sc_R;
//
//    param.m_C1p = m_C1p;
//    param.m_C2p = m_C2p;
//
//    param.m_st_L = m_st_L;
//    param.m_st_R = m_st_R;
//
//    param.U_Wm1 = U_Wm1;
//    param.U_d1  = U_d1;
//    param.U_Wm2 = U_Wm2;
//    param.U_d2  = U_d2;
//
//    param.V_Wp1 = V_Wp1;
//    param.V_u1  = V_u1;
//    param.V_Wp2 = V_Wp2;
//    param.V_u2  = V_u2;
//
//    param.U_st_00 = U_st_00;
//    param.U_st_01 = U_st_01;
//    param.U_st_10 = U_st_10;
//    param.U_st_11 = U_st_11;
//    // updateSpectrum(param);
//}
