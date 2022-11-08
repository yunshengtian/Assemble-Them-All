#include "Force/ForceCuboidCuboidContact.h"
#include "Body/BodyCuboid.h"
#include "CollisionDetection/CollisionDetection.h"
#include "CollisionDetection/Contact.h"
#include "Body/Body.h"
#include "Simulation.h"
#include "Joint/Joint.h"

namespace redmax {

ForceCuboidCuboidContact::ForceCuboidCuboidContact(
    Simulation* sim,
    const BodyCuboid* cuboid1, const BodyCuboid* cuboid2,
    dtype stiffness, dtype damping) : Force(sim), _cuboid1(cuboid1), _cuboid2(cuboid2) {
    _stiffness = stiffness;
    _damping = damping;
}

void ForceCuboidCuboidContact::set_stiffness(dtype stiffness) {
    _stiffness = stiffness;
}

void ForceCuboidCuboidContact::set_damping(dtype damping) {
    _damping = damping;
}

bool ForceCuboidCuboidContact::on_cuboid(const BodyCuboid* cuboid, const Vector3& xw) {
    Vector3 xl = cuboid->_E_i0.topLeftCorner(3, 3) * xw + cuboid->_E_i0.topRightCorner(3, 1);
    for (int i = 0;i < 3;i++) {
        dtype dist = fabs(fabs(xl(i)) - cuboid->_length(i) / 2.);
        if (dist < 1e-5) {
            return true;
        }
    }
    return false;
}

void ForceCuboidCuboidContact::computeForce(VectorX& fm, VectorX& fr, bool verbose) {
    Matrix4 E1 = _cuboid1->_E_0i;
    Matrix3 R1 = E1.topLeftCorner(3, 3);
    Vector3 p1 = E1.topRightCorner(3, 1);
    Matrix4 E2 = _cuboid2->_E_0i;
    Matrix3 R2 = E2.topLeftCorner(3, 3);
    Vector3 p2 = E2.topRightCorner(3, 1);

    Vector3 s1 = _cuboid1->_length;
    Vector3 s2 = _cuboid2->_length;

    std::vector<Contact> contacts;
    contacts.clear();
    collision_detection_cuboid_cuboid(_cuboid1, _cuboid2, contacts);

    for (int i = 0;i < contacts.size();i++) {
        Vector3 xw_c = contacts[i]._xw;
        Vector3 nw = -contacts[i]._normal * contacts[i]._d;
        Vector3 xw1, xw2;
        if (on_cuboid(_cuboid1, xw_c)) {
            xw1 = xw_c;
            xw2 = xw_c + nw;
        } else {
            xw2 = xw_c;
            xw1 = xw_c - nw;
        }
        Vector3 xl1 = _cuboid1->_E_i0.topLeftCorner(3, 3) * xw1 + _cuboid1->_E_i0.topRightCorner(3, 1);
        Vector3 xl2 = _cuboid2->_E_i0.topLeftCorner(3, 3) * xw2 + _cuboid2->_E_i0.topRightCorner(3, 1);

        assert(on_cuboid(_cuboid1, xw1));
        assert(on_cuboid(_cuboid2, xw2));
        
        MatrixX G1 = math::gamma(xl1);
        Vector3 vl1 = G1 * _cuboid1->_phi;
        Vector3 vw1 = R1 * vl1;
        MatrixX G2 = math::gamma(xl2);
        Vector3 vl2 = G2 * _cuboid2->_phi;
        Vector3 vw2 = R2 * vl2;

        Vector3 dx = xw2 - xw1;
        Vector3 dv = vw2 - vw1;
        Matrix3 I = Matrix3::Identity();
        Matrix3 Z = Matrix3::Zero();
        dtype ks = _stiffness;
        dtype kd = _damping;
        Vector3 f = ks * dx + kd * dv;

        fm.segment(_cuboid1->_index[0], 6) += G1.transpose() * (R1.transpose() * f);
        fm.segment(_cuboid2->_index[0], 6) -= G2.transpose() * (R2.transpose() * f);
    }
}

void ForceCuboidCuboidContact::computeForceWithDerivative(
    VectorX& fm, VectorX& fr, MatrixX& Km, MatrixX& Dm, MatrixX& Kr, MatrixX& Dr, bool verbose) {

    Matrix4 E1 = _cuboid1->_E_0i;
    Matrix3 R1 = E1.topLeftCorner(3, 3);
    Vector3 p1 = E1.topRightCorner(3, 1);
    Matrix4 E2 = _cuboid2->_E_0i;
    Matrix3 R2 = E2.topLeftCorner(3, 3);
    Vector3 p2 = E2.topRightCorner(3, 1);

    Vector3 s1 = _cuboid1->_length;
    Vector3 s2 = _cuboid2->_length;

    std::vector<Contact> contacts;
    contacts.clear();
    collision_detection_cuboid_cuboid(_cuboid1, _cuboid2, contacts);
    
    for (int i = 0;i < contacts.size();i++) {
        Vector3 xw_c = contacts[i]._xw;
        Vector3 nw = -contacts[i]._normal * contacts[i]._d;
        Vector3 xw1, xw2;
        if (on_cuboid(_cuboid1, xw_c)) {
            xw1 = xw_c;
            xw2 = xw_c + nw;
        } else {
            xw2 = xw_c;
            xw1 = xw_c - nw;
        }
        Vector3 xl1 = _cuboid1->_E_i0.topLeftCorner(3, 3) * xw1 + _cuboid1->_E_i0.topRightCorner(3, 1);
        Vector3 xl2 = _cuboid2->_E_i0.topLeftCorner(3, 3) * xw2 + _cuboid2->_E_i0.topRightCorner(3, 1);

        assert(on_cuboid(_cuboid1, xw1));
        assert(on_cuboid(_cuboid2, xw2));
        
        MatrixX G1 = math::gamma(xl1);
        Vector3 vl1 = G1 * _cuboid1->_phi;
        Vector3 vw1 = R1 * vl1;
        MatrixX G2 = math::gamma(xl2);
        Vector3 vl2 = G2 * _cuboid2->_phi;
        Vector3 vw2 = R2 * vl2;

        Vector3 dx = xw2 - xw1;
        Vector3 dv = vw2 - vw1;
        Matrix3 I = Matrix3::Identity();
        Matrix3 Z = Matrix3::Zero();
        dtype ks = _stiffness;
        dtype kd = _damping;
        Vector3 f = ks * dx + kd * dv;

        fm.segment(_cuboid1->_index[0], 6) += G1.transpose() * (R1.transpose() * f);
        fm.segment(_cuboid2->_index[0], 6) -= G2.transpose() * (R2.transpose() * f);

        Km.block(_cuboid1->_index[0], _cuboid1->_index[0], 6, 6) += 
            ks * G1.transpose() * (Matrix<dtype, 3, 6>() << math::skew(R1.transpose() * (xw2 - p1)), -I).finished()
            + kd * G1.transpose() * (Matrix<dtype, 3, 6>() << math::skew(R1.transpose() * xw2), Z).finished();
        Dm.block(_cuboid1->_index[0], _cuboid1->_index[0], 6, 6) -=
            kd * G1.transpose() * G1;
        Km.block(_cuboid2->_index[0], _cuboid2->_index[0], 6, 6) +=
            ks * G2.transpose() * (Matrix<dtype, 3, 6>() << math::skew(R2.transpose() * (xw1 - p2)), -I).finished()
            + kd * G2.transpose() * (Matrix<dtype, 3, 6>() << math::skew(R2.transpose() * xw1), Z).finished();
        Dm.block(_cuboid2->_index[0], _cuboid2->_index[0], 6, 6) -=
            kd * G2.transpose() * G2;
        Km.block(_cuboid1->_index[0], _cuboid2->_index[0], 6, 6) +=
            ks * G1.transpose() * (R1.transpose() * (R2 * (Matrix<dtype, 3, 6>() << -math::skew(xl2), I).finished()))
            - kd * G1.transpose() * (R1.transpose() * (R2 * (Matrix<dtype, 3, 6>() << math::skew(vl2), Z).finished()));
        Km.block(_cuboid2->_index[0], _cuboid1->_index[0], 6, 6) += 
            ks * G2.transpose() * (R2.transpose() * (R1 * (Matrix<dtype, 3, 6>() << -math::skew(xl1), I).finished()))
            - kd * G2.transpose() * (R2.transpose() * (R1 * (Matrix<dtype, 3, 6>() << -math::skew(vl1), Z).finished()));
        Dm.block(_cuboid1->_index[0], _cuboid2->_index[0], 6, 6) += 
            kd * G1.transpose() * R1.transpose() * R2 * G2;
        Dm.block(_cuboid2->_index[0], _cuboid1->_index[0], 6, 6) +=
            kd * G2.transpose() * R2.transpose() * R1 * G1;
    }
}

}