#include <quadmath.h>
#include "clooptools.h"
#include "marty/looptools_init.h"
#include "marty/looptools_quad_extension.h"
#include <cmath>
#include "marty/looptools_interface.h"
#include "m_N_1.h"
#include "common.h"

#include "params.h"
#include "libcomplexop.h"

namespace c9_nmfv {

complex_t m_N_1(
        param_t const &param
        )
{
    clearcache();
    const real_t m_N_1 = param.m_N_1;
    const complex_t IT_0000 = 0.5*m_N_1;
    return IT_0000;
}
} // End of namespace c9_nmfv
