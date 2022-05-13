#include <quadmath.h>
#include "clooptools.h"
#include "marty/looptools_init.h"
#include "marty/looptools_quad_extension.h"
#include <cmath>
#include "marty/looptools_interface.h"
#include "m_sG.h"
#include "common.h"

#include "params.h"
#include "libcomplexop.h"

namespace c9_nmfv {

complex_t m_sG(
        param_t const &param
        )
{
    clearcache();
    const real_t M_3 = param.M_3;
    const complex_t IT_0000 = 0.5*M_3;
    return IT_0000;
}
} // End of namespace c9_nmfv
