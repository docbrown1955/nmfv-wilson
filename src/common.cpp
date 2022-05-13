#include <quadmath.h>
#include "clooptools.h"
#include "marty/looptools_init.h"
#include "marty/looptools_quad_extension.h"
#include <cmath>
#include "marty/looptools_interface.h"
#include "common.h"

namespace c9_nmfv {

void setMu(const double mu)
{
    setmudim(mu * mu);
}

void setLambda2(const double lambda2)
{
    setlambda(lambda2);
}

void setUVDiv(const double x)
{
    setuvdiv(x);
}

} // End of namespace c9_nmfv

