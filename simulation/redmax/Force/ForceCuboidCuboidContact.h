#pragma once
#include "Common.h"
#include "Utils.h"
#include "Force/Force.h"

namespace redmax {

class BodyCuboid;

class ForceCuboidCuboidContact : public Force {
public:
    const BodyCuboid* _cuboid1; // cuboid bodies to consider contact
    const BodyCuboid* _cuboid2;
    dtype _stiffness;
    dtype _damping;

    ForceCuboidCuboidContact(
        Simulation* sim,
        const BodyCuboid* cuboid1, const BodyCuboid* cuboid2,
        dtype stiffness = 1., dtype damping = 0.);
    
    void set_stiffness(dtype stiffness);
    void set_damping(dtype damping);

    bool on_cuboid(const BodyCuboid* cuboid, const Vector3& xw);

    void computeForce(VectorX& fm, VectorX& fr, bool verbose = false);
    void computeForceWithDerivative(VectorX& fm, VectorX& fr, MatrixX& Km, MatrixX& Dm, MatrixX& Kr, MatrixX& Dr, bool verbose = false);
};

}