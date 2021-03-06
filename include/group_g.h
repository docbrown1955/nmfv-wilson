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

#ifndef CSL_LIB_C9_NMFV_G_H_INCLUDED
#define CSL_LIB_C9_NMFV_G_H_INCLUDED

#include <array>
#include "common.h"
#include "librarytensor.h"
#include "callable.h"
#include "csl/initSanitizer.h"
#include "params.h"
#include "m_Z.h"
#include "m_sG.h"
#include "m_N_1.h"
#include "m_N_2.h"
#include "m_N_3.h"
#include "m_N_4.h"
#include "m_sG.h"
#include "m_N_1.h"
#include "m_N_2.h"
#include "m_N_3.h"
#include "m_N_4.h"
#include "m_Gp.h"
#include "m_G0.h"
#include "C9B_N.h"
#include "C10B_N.h"
#include "C9Z_N.h"
#include "C10Z_N.h"
#include "C9A_N.h"
#include "C10A_N.h"
#include "C9B_G.h"
#include "C10B_G.h"
#include "C9Z_G.h"
#include "C10Z_G.h"
#include "C9A_G.h"
#include "C10A_G.h"
#include "C9B_C.h"
#include "C10B_C.h"
#include "C9Z_C.h"
#include "C10Z_C.h"
#include "C9A_C.h"
#include "C10A_C.h"
#include "C9B_H.h"
#include "C10B_H.h"
#include "C9Z_H.h"
#include "C10Z_H.h"
#include "C9A_H.h"
#include "C10A_H.h"
#include "C7_G.h"
#include "C7p_G.h"
#include "C7_C.h"
#include "C7p_C.h"
#include "C7_H.h"
#include "C7p_H.h"
#include "C7_N.h"
#include "C7p_N.h"
#include "DiM_mu_N.h"
#include "DiM_mup_N.h"
#include "DiM_mu_G.h"
#include "DiM_mup_G.h"
#include "DiM_mu_C.h"
#include "DiM_mup_C.h"
#include "DiM_mu_H.h"
#include "DiM_mup_H.h"
#include "DiM_mu_NH.h"
#include "DiM_mup_NH.h"

namespace c9_nmfv {


inline std::array<Callable<complex_t, param_t>, 55> f_G = {
    Callable{"m_Z", m_Z},
    Callable{"m_sG", m_sG},
    Callable{"m_N_1", m_N_1},
    Callable{"m_N_2", m_N_2},
    Callable{"m_N_3", m_N_3},
    Callable{"m_N_4", m_N_4},
    Callable{"m_sG", m_sG},
    Callable{"m_N_1", m_N_1},
    Callable{"m_N_2", m_N_2},
    Callable{"m_N_3", m_N_3},
    Callable{"m_N_4", m_N_4},
    Callable{"m_Gp", m_Gp},
    Callable{"m_G0", m_G0},
    Callable{"C9B_N", C9B_N},
    Callable{"C10B_N", C10B_N},
    Callable{"C9Z_N", C9Z_N},
    Callable{"C10Z_N", C10Z_N},
    Callable{"C9A_N", C9A_N},
    Callable{"C10A_N", C10A_N},
    Callable{"C9B_G", C9B_G},
    Callable{"C10B_G", C10B_G},
    Callable{"C9Z_G", C9Z_G},
    Callable{"C10Z_G", C10Z_G},
    Callable{"C9A_G", C9A_G},
    Callable{"C10A_G", C10A_G},
    Callable{"C9B_C", C9B_C},
    Callable{"C10B_C", C10B_C},
    Callable{"C9Z_C", C9Z_C},
    Callable{"C10Z_C", C10Z_C},
    Callable{"C9A_C", C9A_C},
    Callable{"C10A_C", C10A_C},
    Callable{"C9B_H", C9B_H},
    Callable{"C10B_H", C10B_H},
    Callable{"C9Z_H", C9Z_H},
    Callable{"C10Z_H", C10Z_H},
    Callable{"C9A_H", C9A_H},
    Callable{"C10A_H", C10A_H},
    Callable{"C7_G", C7_G},
    Callable{"C7p_G", C7p_G},
    Callable{"C7_C", C7_C},
    Callable{"C7p_C", C7p_C},
    Callable{"C7_H", C7_H},
    Callable{"C7p_H", C7p_H},
    Callable{"C7_N", C7_N},
    Callable{"C7p_N", C7p_N},
    Callable{"DiM_mu_N", DiM_mu_N},
    Callable{"DiM_mup_N", DiM_mup_N},
    Callable{"DiM_mu_G", DiM_mu_G},
    Callable{"DiM_mup_G", DiM_mup_G},
    Callable{"DiM_mu_C", DiM_mu_C},
    Callable{"DiM_mup_C", DiM_mup_C},
    Callable{"DiM_mu_H", DiM_mu_H},
    Callable{"DiM_mup_H", DiM_mup_H},
    Callable{"DiM_mu_NH", DiM_mu_NH},
    Callable{"DiM_mup_NH", DiM_mup_NH},
};

inline std::map<std::string, Callable<complex_t, param_t>> fmap_G {
    {"m_Z", f_G[0]},
    {"m_sG", f_G[1]},
    {"m_N_1", f_G[2]},
    {"m_N_2", f_G[3]},
    {"m_N_3", f_G[4]},
    {"m_N_4", f_G[5]},
    {"m_sG", f_G[6]},
    {"m_N_1", f_G[7]},
    {"m_N_2", f_G[8]},
    {"m_N_3", f_G[9]},
    {"m_N_4", f_G[10]},
    {"m_Gp", f_G[11]},
    {"m_G0", f_G[12]},
    {"C9B_N", f_G[13]},
    {"C10B_N", f_G[14]},
    {"C9Z_N", f_G[15]},
    {"C10Z_N", f_G[16]},
    {"C9A_N", f_G[17]},
    {"C10A_N", f_G[18]},
    {"C9B_G", f_G[19]},
    {"C10B_G", f_G[20]},
    {"C9Z_G", f_G[21]},
    {"C10Z_G", f_G[22]},
    {"C9A_G", f_G[23]},
    {"C10A_G", f_G[24]},
    {"C9B_C", f_G[25]},
    {"C10B_C", f_G[26]},
    {"C9Z_C", f_G[27]},
    {"C10Z_C", f_G[28]},
    {"C9A_C", f_G[29]},
    {"C10A_C", f_G[30]},
    {"C9B_H", f_G[31]},
    {"C10B_H", f_G[32]},
    {"C9Z_H", f_G[33]},
    {"C10Z_H", f_G[34]},
    {"C9A_H", f_G[35]},
    {"C10A_H", f_G[36]},
    {"C7_G", f_G[37]},
    {"C7p_G", f_G[38]},
    {"C7_C", f_G[39]},
    {"C7p_C", f_G[40]},
    {"C7_H", f_G[41]},
    {"C7p_H", f_G[42]},
    {"C7_N", f_G[43]},
    {"C7p_N", f_G[44]},
    {"DiM_mu_N", f_G[45]},
    {"DiM_mup_N", f_G[46]},
    {"DiM_mu_G", f_G[47]},
    {"DiM_mup_G", f_G[48]},
    {"DiM_mu_C", f_G[49]},
    {"DiM_mup_C", f_G[50]},
    {"DiM_mu_H", f_G[51]},
    {"DiM_mup_H", f_G[52]},
    {"DiM_mu_NH", f_G[53]},
    {"DiM_mup_NH", f_G[54]},
};


}
 // End of namespace c9_nmfv

#endif
