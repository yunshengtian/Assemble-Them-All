#pragma once
#include "Common.h"
#include "Utils.h"
#include "Force/ForceGeneralSDFContact.h"
#include "CollisionDetection/Contact.h"

namespace redmax {

class Body;
class BodySDFObj;

class ForceGeneralSDFContactWithSuction : public ForceGeneralSDFContact {
public:
    dtype _skn;
    bool _suction_active;
    std::vector<int> _suction_contact_ids;
    std::vector<Vector3> _surface_contact_pts;

    ForceGeneralSDFContactWithSuction(
        Simulation* sim,
        Body* contact_body, Body* SDF_body,
        dtype kn = 1., dtype kt = 0.,
        dtype mu = 0., dtype damping = 0., dtype skn = 1.);

    void enable_suction();
    void disable_suction();

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

private:
    void collision_detection(std::vector<Contact>& contacts, std::vector<Contact>& suction_contacts);

    void computeForce(std::vector<Contact> &contacts, VectorX& fm, bool verbose = false);
    void computeForceWithDerivative(std::vector<Contact> &contacts, VectorX& fm, MatrixX& Km, MatrixX& Dm, bool verbose = false);
    void computeForceWithDerivative(std::vector<Contact> &contacts, VectorX& fm, MatrixX& Km, MatrixX& Dm, MatrixX& dfm_dp, bool verbose = false);
};

}