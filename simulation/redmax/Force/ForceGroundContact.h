#pragma once
#include "Force/Force.h"
#include "CollisionDetection/Contact.h"

namespace redmax {

class Body;

class ForceGroundContact : public Force {
public:
    Body* _contact_body;    // the body to contact with ground.
    Matrix4 _E_g;           // the ground frame, with Z-up.
    dtype _kn;              // normal stiffness
    dtype _kt;              // tangential stiffness
    dtype _mu;              // coefficient of friction
    dtype _damping;         // damping of the contact force

    ForceGroundContact(Simulation* sim, Body* contact_body, Matrix4 E_g, dtype kn = 1., dtype kt = 0., dtype mu = 0., dtype damping = 0.);

    void set_transform(Matrix4 E_g);
    void set_stiffness(dtype kn, dtype kt);
    void set_friction(dtype mu);
    void set_damping(dtype damping);

    void computeForce(VectorX& fm, VectorX& fr, bool verbose = false);
    void computeForceWithDerivative(
        VectorX& fm, VectorX& fr, 
        MatrixX& Km, MatrixX& Dm, 
        MatrixX& Kr, MatrixX& Dr, 
        bool verbose = false);
    void computeForceWithDerivative(
        VectorX& fm, VectorX& fr, 
        MatrixX& Km, MatrixX& Dm, 
        MatrixX& Kr, MatrixX& Dr, 
        MatrixX& dfm_dp, MatrixX& dfr_dp,
        bool verbose = false);

    void test_derivatives_runtime();
    void test_design_derivatives_runtime();

private:
    void computeForce(std::vector<Contact> &contacts, VectorX& fm, bool verbose = false);
    void computeForceWithDerivative(std::vector<Contact> &contacts, VectorX& fm, MatrixX& Km, MatrixX& Dm, bool verbose = false);
    void computeForceWithDerivative(std::vector<Contact> &contacts, VectorX& fm, MatrixX& Km, MatrixX& Dm, MatrixX& dfm_dp, bool verbose = false);
};

}