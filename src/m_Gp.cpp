#include <quadmath.h>
#include "clooptools.h"
#include "marty/looptools_init.h"
#include "marty/looptools_quad_extension.h"
#include <cmath>
#include "marty/looptools_interface.h"
#include "m_Gp.h"
#include "common.h"

#include "params.h"
#include "libcomplexop.h"

namespace c9_nmfv {

complex_t m_Gp(
        param_t const &param
        )
{
    clearcache();
    const real_t M_W = param.M_W;
    return M_W;
}
} // End of namespace c9_nmfv
