#pragma once
#include "Force/Force.h"

namespace redmax {

class Body;

class ForceSpringMultiPointGeneric : public Force {
public:
    std::vector<Body*> _bodies; // bodies to apply the forces.
    std::vector<Vector3> _xls; // application locations in the local coordinate.

    ForceSpringMultiPointGeneric(Simulation* sim);

    void addBodyPoint(Body* body, Vector3 loc);

    void computeForce(VectorX& fm, VectorX& fr, bool verbose = false);
    void computeForceWithDerivative(VectorX& fm, VectorX& fr, MatrixX& Km, MatrixX& Dm, MatrixX& Kr, MatrixX& Dr, bool verbose = false);

    virtual void computeSpringForce(dtype l, dtype l_dot, dtype& f) = 0;
    virtual void computeSpringForceWithDerivative(dtype l, dtype l_dot, dtype& f, dtype& df_dl, dtype& df_dldot, bool verbose = false) = 0;
};

}