#include <quadmath.h>
#include "clooptools.h"
#include "marty/looptools_init.h"
#include "marty/looptools_quad_extension.h"
#include <cmath>
#include "marty/looptools_interface.h"
#include "m_G0.h"
#include "common.h"

#include "params.h"
#include "libcomplexop.h"

namespace c9_nmfv {

complex_t m_G0(
        param_t const &param
        )
{
    clearcache();
    const real_t M_Z = param.M_Z;
    return M_Z;
}
} // End of namespace c9_nmfv
