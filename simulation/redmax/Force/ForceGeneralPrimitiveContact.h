#pragma once
#include "Common.h"
#include "Utils.h"
#include "Force/Force.h"
#include "CollisionDetection/Contact.h"

namespace redmax {

class Body;
class BodyPrimitiveShape;

class ForceGeneralPrimitiveContact : public Force {
public:
    Body* _contact_body; // the general contact body should be able to give a list of contact points on the surface.
    BodyPrimitiveShape* _primitive_body; // the primitive body is supposed to have an anlytical distance field
    dtype _kn;              // normal stiffness
    dtype _kt;              // tangential stiffness
    dtype _mu;              // coefficient of friction
    dtype _damping;         // damping of the contact force
    dtype _scale;           // a scale parameter to control the stiffness (for continuation method)

    ForceGeneralPrimitiveContact(
        Simulation* sim,
        Body* contact_body, Body* primitive_body,
        dtype kn = 1., dtype kt = 0.,
        dtype mu = 0., dtype damping = 0.);
    
    void set_stiffness(dtype kn, dtype kt);
    void set_friction(dtype mu);
    void set_damping(dtype damping);
    void set_scale(dtype scale);

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

    void test_derivatives();
    void test_derivatives_runtime();
    void test_design_derivatives_runtime();

private:
    void computeForce(std::vector<Contact> &contacts, VectorX& fm, bool verbose = false);
    void computeForceWithDerivative(std::vector<Contact> &contacts, VectorX& fm, MatrixX& Km, MatrixX& Dm, bool verbose = false);
    void computeForceWithDerivative(std::vector<Contact> &contacts, VectorX& fm, MatrixX& Km, MatrixX& Dm, MatrixX& dfm_dp, bool verbose = false);
};

}