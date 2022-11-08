#include "Force/ForceGeneralSDFContactWithSuction.h"
#include "CollisionDetection/Contact.h"
#include "CollisionDetection/CollisionDetection.h"
#include "Body/Body.h"
#include "Body/BodySDFObj.h"
#include "Body/BodyCuboid.h"
#include "Simulation.h"
#include "Joint/Joint.h"

namespace redmax {

ForceGeneralSDFContactWithSuction::ForceGeneralSDFContactWithSuction(
    Simulation* sim,
    Body* contact_body, Body* SDF_body,
    dtype kn, dtype kt, dtype mu, dtype damping, dtype skn) : ForceGeneralSDFContact(sim, contact_body, SDF_body, kn, kt, mu, damping) {

    _suction_active = false;
    _skn = skn;
}

void ForceGeneralSDFContactWithSuction::enable_suction() {
    _suction_active = true;
}

void ForceGeneralSDFContactWithSuction::disable_suction() {
    _suction_active = false;
    _suction_contact_ids.clear();
    _surface_contact_pts.clear();
}

void ForceGeneralSDFContactWithSuction::collision_detection(std::vector<Contact>& contacts, std::vector<Contact>& suction_contacts) {

    collision_detection_general_SDF(_contact_body, _SDF_body, contacts);
    if (_suction_active) {
        collision_detection_general_SDF_with_suction(_contact_body, _SDF_body, suction_contacts, _suction_contact_ids, _surface_contact_pts);
    }
}

void ForceGeneralSDFContactWithSuction::computeForce(VectorX& fm, VectorX& fr, bool verbose) {
    auto t0 = clock();

    // detect the contact points between two bodies
    std::vector<Contact> contacts, suction_contacts;
    contacts.clear();
    suction_contacts.clear();
    collision_detection(contacts, suction_contacts);

    auto t1 = clock();

    VectorX fm_sub;
    ForceGeneralSDFContact::computeForce(contacts, fm_sub, verbose);
    if (_suction_active) {
        computeForce(suction_contacts, fm_sub, verbose);
    }
    
    fm.segment(_contact_body->_index[0], 6) += fm_sub.head(6);
    fm.segment(_SDF_body->_index[0], 6) += fm_sub.tail(6);

    auto t2 = clock();

    _time.add("collision detection", t1 - t0);
    _time.add("compute force", t2 - t1);
}

void ForceGeneralSDFContactWithSuction::computeForceWithDerivative(
    VectorX& fm, VectorX& fr, 
    MatrixX& Km, MatrixX& Dm, 
    MatrixX& Kr, MatrixX& Dr,
    bool verbose) {
    
    auto t0 = clock();

    // detect the contact points between two bodies
    std::vector<Contact> contacts, suction_contacts;
    contacts.clear();
    suction_contacts.clear();
    collision_detection(contacts, suction_contacts);

    auto t1 = clock();

    VectorX fm_sub;
    MatrixX Km_sub, Dm_sub;
    ForceGeneralSDFContact::computeForceWithDerivative(contacts, fm_sub, Km_sub, Dm_sub, verbose);
    if (_suction_active) {
        computeForceWithDerivative(suction_contacts, fm_sub, Km_sub, Dm_sub, verbose);
    }

    fm.segment(_contact_body->_index[0], 6) += fm_sub.head(6);
    fm.segment(_SDF_body->_index[0], 6) += fm_sub.tail(6);

    Km.block(_contact_body->_index[0], _contact_body->_index[0], 6, 6) += Km_sub.topLeftCorner(6, 6);
    Km.block(_contact_body->_index[0], _SDF_body->_index[0], 6, 6) += Km_sub.topRightCorner(6, 6);
    Km.block(_SDF_body->_index[0], _contact_body->_index[0], 6, 6) += Km_sub.bottomLeftCorner(6, 6);
    Km.block(_SDF_body->_index[0], _SDF_body->_index[0], 6, 6) += Km_sub.bottomRightCorner(6, 6);
    
    Dm.block(_contact_body->_index[0], _contact_body->_index[0], 6, 6) += Dm_sub.topLeftCorner(6, 6);
    Dm.block(_contact_body->_index[0], _SDF_body->_index[0], 6, 6) += Dm_sub.topRightCorner(6, 6);
    Dm.block(_SDF_body->_index[0], _contact_body->_index[0], 6, 6) += Dm_sub.bottomLeftCorner(6, 6);
    Dm.block(_SDF_body->_index[0], _SDF_body->_index[0], 6, 6) += Dm_sub.bottomRightCorner(6, 6);

    auto t2 = clock();

    _time.add("collision detection", t1 - t0);
    _time.add("compute force", t2 - t1);
}

void ForceGeneralSDFContactWithSuction::computeForceWithDerivative(
    VectorX& fm, VectorX& fr, 
    MatrixX& Km, MatrixX& Dm, 
    MatrixX& Kr, MatrixX& Dr,
    MatrixX& dfm_dp, MatrixX& dfr_dp,
    bool verbose) {

    auto t0 = clock();

    // detect the contact points between two bodies
    std::vector<Contact> contacts, suction_contacts;
    contacts.clear();
    suction_contacts.clear();
    collision_detection(contacts, suction_contacts);

    auto t1 = clock();

    VectorX fm_sub;
    MatrixX Km_sub, Dm_sub;
    ForceGeneralSDFContact::computeForceWithDerivative(contacts, fm_sub, Km_sub, Dm_sub, dfm_dp, verbose);
    if (_suction_active) {
        computeForceWithDerivative(suction_contacts, fm_sub, Km_sub, Dm_sub, dfm_dp, verbose);
    }

    fm.segment(_contact_body->_index[0], 6) += fm_sub.head(6);
    fm.segment(_SDF_body->_index[0], 6) += fm_sub.tail(6);

    Km.block(_contact_body->_index[0], _contact_body->_index[0], 6, 6) += Km_sub.topLeftCorner(6, 6);
    Km.block(_contact_body->_index[0], _SDF_body->_index[0], 6, 6) += Km_sub.topRightCorner(6, 6);
    Km.block(_SDF_body->_index[0], _contact_body->_index[0], 6, 6) += Km_sub.bottomLeftCorner(6, 6);
    Km.block(_SDF_body->_index[0], _SDF_body->_index[0], 6, 6) += Km_sub.bottomRightCorner(6, 6);
    
    Dm.block(_contact_body->_index[0], _contact_body->_index[0], 6, 6) += Dm_sub.topLeftCorner(6, 6);
    Dm.block(_contact_body->_index[0], _SDF_body->_index[0], 6, 6) += Dm_sub.topRightCorner(6, 6);
    Dm.block(_SDF_body->_index[0], _contact_body->_index[0], 6, 6) += Dm_sub.bottomLeftCorner(6, 6);
    Dm.block(_SDF_body->_index[0], _SDF_body->_index[0], 6, 6) += Dm_sub.bottomRightCorner(6, 6);

    auto t2 = clock();

    _time.add("collision detection", t1 - t0);
    _time.add("compute force", t2 - t1);
}

void ForceGeneralSDFContactWithSuction::computeForce(std::vector<Contact> &contacts, VectorX& fm, bool verbose) {
    Matrix3 R1 = _contact_body->_E_0i.topLeftCorner(3, 3);
    Vector3 p1 = _contact_body->_E_0i.topRightCorner(3, 1);
    Vector6 phi1 = _contact_body->_phi;
    Matrix3 R2 = _SDF_body->_E_0i.topLeftCorner(3, 3);
    Vector3 p2 = _SDF_body->_E_0i.topRightCorner(3, 1);
    Vector6 phi2 = _SDF_body->_phi;
    
    // fm = VectorX::Zero(12);

    dtype kn = _skn * _scale * _contact_body->_contact_scale;
    dtype kt = _kt * _scale * _contact_body->_contact_scale;
    dtype ks = kn;
    dtype kd = _damping;

    for (int i = 0;i < contacts.size();i++) {
        Vector3 xi1 = contacts[i]._xi;
        Vector3 xw1 = R1 * xi1 + p1;
        Vector3 xi2 = _surface_contact_pts[i];
        Vector3 xw2 = R2 * xi2 + p2;
        Vector3 dx = xw2 - xw1;
        
        MatrixX G1 = math::gamma(xi1);
        MatrixX G2 = math::gamma(xi2);
        MatrixX GTRT1 = G1.transpose() * R1.transpose();
        MatrixX GTRT2 = G2.transpose() * R2.transpose();

        Vector3 xi1_dot = G1 * phi1;
        Vector3 xw1_dot = R1 * xi1_dot;
        Vector3 xi2_dot = G2 * phi2;
        Vector3 xw2_dot = R2 * xi2_dot;
        Vector3 dv = xw2_dot - xw1_dot;

        // contact force
        Vector6 fc1 = GTRT1 * (ks * dx + kd * dv);
        fm.head(6) += fc1;
        Vector6 fc2 = -GTRT2 * (ks * dx + kd * dv);
        fm.tail(6) += fc2;
    }
}

void ForceGeneralSDFContactWithSuction::computeForceWithDerivative(
    std::vector<Contact> &contacts,
    VectorX& fm, MatrixX& Km, MatrixX& Dm,
    bool verbose) {
    
    // fm = VectorX::Zero(12);
    // Km = MatrixX::Zero(12, 12);
    // Dm = MatrixX::Zero(12, 12);

    Matrix3 R1 = _contact_body->_E_0i.topLeftCorner(3, 3);
    Vector3 p1 = _contact_body->_E_0i.topRightCorner(3, 1);
    Vector6 phi1 = _contact_body->_phi;
    Matrix3 R2 = _SDF_body->_E_0i.topLeftCorner(3, 3);
    Vector3 p2 = _SDF_body->_E_0i.topRightCorner(3, 1);
    Vector6 phi2 = _SDF_body->_phi;

    dtype kn = _skn * _scale * _contact_body->_contact_scale;
    dtype kt = _kt * _scale * _contact_body->_contact_scale;
    dtype ks = kn;
    dtype kd = _damping;

    for (int i = 0;i < contacts.size();i++) {
        Vector3 xi1 = contacts[i]._xi;
        Vector3 xw1 = R1 * xi1 + p1;
        Vector3 xi2 = _surface_contact_pts[i];
        Vector3 xw2 = R2 * xi2 + p2;
        Vector3 dx = xw2 - xw1;
        
        MatrixX G1 = math::gamma(xi1);
        MatrixX G2 = math::gamma(xi2);
        MatrixX GTRT1 = G1.transpose() * R1.transpose();
        MatrixX GTRT2 = G2.transpose() * R2.transpose();

        Vector3 xi1_dot = G1 * phi1;
        Vector3 xw1_dot = R1 * xi1_dot;
        Vector3 xi2_dot = G2 * phi2;
        Vector3 xw2_dot = R2 * xi2_dot;
        Vector3 dv = xw2_dot - xw1_dot;

        // contact force
        Vector6 fc1 = GTRT1 * (ks * dx + kd * dv);
        fm.head(6) += fc1;
        Vector6 fc2 = -GTRT2 * (ks * dx + kd * dv);
        fm.tail(6) += fc2;

        Matrix3 I = Matrix3::Identity();
        Matrix3 Z = Matrix3::Zero();

        Vector3 vi1 = G1 * phi1;
        Vector3 vi2 = G2 * phi2;

        // derivatives
        Km.topLeftCorner(6, 6) += 
            ks * G1.transpose() * (Matrix<dtype, 3, 6>() << math::skew(R1.transpose() * (xw2 - p1)), -I).finished()
            + kd * G1.transpose() * (Matrix<dtype, 3, 6>() << math::skew(R1.transpose() * xw2), Z).finished();
        Dm.topLeftCorner(6, 6) -=
            kd * G1.transpose() * G1;
        Km.bottomRightCorner(6, 6) +=
            ks * G2.transpose() * (Matrix<dtype, 3, 6>() << math::skew(R2.transpose() * (xw1 - p2)), -I).finished()
            + kd * G2.transpose() * (Matrix<dtype, 3, 6>() << math::skew(R2.transpose() * xw1), Z).finished();
        Dm.bottomRightCorner(6, 6) -=
            kd * G2.transpose() * G2;
        Km.topRightCorner(6, 6) +=
            ks * G1.transpose() * (R1.transpose() * (R2 * (Matrix<dtype, 3, 6>() << -math::skew(xi2), I).finished()))
            - kd * G1.transpose() * (R1.transpose() * (R2 * (Matrix<dtype, 3, 6>() << math::skew(vi2), Z).finished()));
        Km.bottomLeftCorner(6, 6) += 
            ks * G2.transpose() * (R2.transpose() * (R1 * (Matrix<dtype, 3, 6>() << -math::skew(xi1), I).finished()))
            - kd * G2.transpose() * (R2.transpose() * (R1 * (Matrix<dtype, 3, 6>() << -math::skew(vi1), Z).finished()));
        Dm.topRightCorner(6, 6) += 
            kd * G1.transpose() * R1.transpose() * R2 * G2;
        Dm.bottomLeftCorner(6, 6) +=
            kd * G2.transpose() * R2.transpose() * R1 * G1;

    }
}

void ForceGeneralSDFContactWithSuction::computeForceWithDerivative(
    std::vector<Contact> &contacts,
    VectorX& fm, MatrixX& Km, MatrixX& Dm, MatrixX& dfm_dp, 
    bool verbose) {
    
    computeForceWithDerivative(contacts, fm, Km, Dm, verbose); // TODO: not implemented yet
}

void ForceGeneralSDFContactWithSuction::test_derivatives() {
    std::cerr << "**************************** General-SDF Collision With Suction Derivatives ***************************" << std::endl;
    srand(time(0));
    // srand(1000);

    dtype eps = 1e-7;
    // generate random xi1, E1, E2, phi1, phi2
    Vector3 xi1 = Vector3::Random();
    Eigen::Quaternion<dtype> quat_1(Vector4::Random());
    quat_1.normalize();
    Matrix4 E1 = Matrix4::Identity();
    E1.topLeftCorner(3, 3) = quat_1.toRotationMatrix();
    E1.topRightCorner(3, 1) = Vector3::Random();
    Eigen::Quaternion<dtype> quat_2(Vector4::Random());
    quat_2.normalize();
    Matrix4 E2 = Matrix4::Identity();
    E2.topLeftCorner(3, 3) = quat_2.toRotationMatrix();
    E2.topRightCorner(3, 1) = Vector3::Random();
    Vector6 phi1 = Vector6::Random();
    Vector6 phi2 = Vector6::Random();
    // Matrix4 E1 = Matrix4::Identity();
    // Matrix4 E2 = Matrix4::Identity();
    // Vector6 phi1 = Vector6::Random();
    // Vector6 phi2 = Vector6::Random();
    Vector3 xw1 = E1.topLeftCorner(3, 3) * xi1 + E1.topRightCorner(3, 1);

    this->_contact_body->_E_0i = E1;
    this->_contact_body->_phi = phi1;
    this->_SDF_body->_E_0i = E2;
    this->_SDF_body->_phi = phi2;

    std::vector<Contact> contacts, suction_contacts;

    enable_suction();
    collision_detection(contacts, suction_contacts);

    VectorX fm;
    MatrixX Km, Dm;
    ForceGeneralSDFContact::computeForceWithDerivative(contacts, fm, Km, Dm, false);
    computeForceWithDerivative(suction_contacts, fm, Km, Dm, true);

    // Km
    MatrixX Km_fd = MatrixX::Zero(12, 12);
    for (int i = 0;i < 6;i++) {
        Vector6 dq = Vector6::Zero();
        dq[i] = eps;
        Matrix4 E1_pos = E1 * math::exp(dq);
        this->_contact_body->_E_0i = E1_pos;
        VectorX fm_pos;
        ForceGeneralSDFContact::computeForce(contacts, fm_pos);
        computeForce(suction_contacts, fm_pos);
        Km_fd.col(i) = (fm_pos - fm) / eps;
    }
    this->_contact_body->_E_0i = E1;
    
    print_error("Km(q1) ", Km.leftCols(6), Km_fd.leftCols(6));

    for (int i = 0;i < 6;i++) {
        Vector6 dq = Vector6::Zero();
        dq[i] = eps;
        Matrix4 E2_pos = E2 * math::exp(dq);
        this->_SDF_body->_E_0i = E2_pos;
        VectorX fm_pos;
        ForceGeneralSDFContact::computeForce(contacts, fm_pos);
        computeForce(suction_contacts, fm_pos);
        Km_fd.col(6 + i) = (fm_pos - fm) / eps;
    }
    this->_SDF_body->_E_0i = E2;
    print_error("Km(q2) ", Km.rightCols(6), Km_fd.rightCols(6));

    // Dm
    MatrixX Dm_fd = MatrixX::Zero(12, 12);
    for (int i = 0;i < 6;i++) {
        Vector6 phi1_pos = phi1;
        phi1_pos[i] += eps;
        this->_contact_body->_phi = phi1_pos;
        VectorX fm_pos;
        ForceGeneralSDFContact::computeForce(contacts, fm_pos);
        computeForce(suction_contacts, fm_pos);
        Dm_fd.col(i) = (fm_pos - fm) / eps;
    }
    this->_contact_body->_phi = phi1;
    print_error("Dm(phi1) ", Dm.leftCols(6), Dm_fd.leftCols(6));

    for (int i = 0;i < 6;i++) {
        Vector6 phi2_pos = phi2;
        phi2_pos[i] += eps;
        this->_SDF_body->_phi = phi2_pos;
        VectorX fm_pos;
        ForceGeneralSDFContact::computeForce(contacts, fm_pos);
        computeForce(suction_contacts, fm_pos);
        Dm_fd.col(6 + i) = (fm_pos - fm) / eps;
    }
    this->_SDF_body->_phi = phi2;
    print_error("Dm(phi2) ", Dm.rightCols(6), Dm_fd.rightCols(6));
}

}