#pragma once
#include "Force/ForceSpringMultiPointGeneric.h"

namespace redmax {

class ForceCable : public ForceSpringMultiPointGeneric {
public:
    dtype _l_rest; // resting length
    dtype _stiffness;
    dtype _damping;

    ForceCable(Simulation* sim, dtype stiffness = 1., dtype damping = 1.);

    void set_stiffness(dtype stiffness);
    void set_damping(dtype damping);

    void init();

    void computeSpringForce(dtype l, dtype ldot, dtype& f);
    void computeSpringForceWithDerivative(dtype l, dtype ldot, dtype& f, dtype& df_dl, dtype& df_dldot, bool verbose = false);
};

}