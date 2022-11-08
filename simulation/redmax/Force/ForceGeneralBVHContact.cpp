#include "Force/ForceGeneralBVHContact.h"
#include "CollisionDetection/Contact.h"
#include "CollisionDetection/CollisionDetection.h"
#include "Body/Body.h"
#include "Body/BodyBVHObj.h"
#include "Body/BodyCuboid.h"
#include "Simulation.h"
#include "Joint/Joint.h"

namespace redmax {

ForceGeneralBVHContact::ForceGeneralBVHContact(
    Simulation* sim,
    Body* contact_body, Body* BVH_body,
    dtype kn, dtype kt, dtype mu, dtype damping) : Force(sim) {
    _contact_body = contact_body;
    _BVH_body = dynamic_cast<BodyBVHObj*>(BVH_body);
    if (_BVH_body == nullptr) {
        throw_error("The second body in Contact should be BVH body.");
    }

    _kn = kn;
    _kt = kt;
    _mu = mu;
    _damping = damping;
    _scale = 1.;
}

void ForceGeneralBVHContact::set_stiffness(dtype kn, dtype kt) {
    _kn = kn;
    _kt = kt;
}

void ForceGeneralBVHContact::set_friction(dtype mu) {
    _mu = mu;
}

void ForceGeneralBVHContact::set_damping(dtype damping) {
    _damping = damping;
}

void ForceGeneralBVHContact::set_scale(dtype scale) {
    _scale = scale;
}

void ForceGeneralBVHContact::computeForce(VectorX& fm, VectorX& fr, bool verbose) {
    auto t0 = clock();

    // detect the contact points between two bodies
    std::vector<Contact> contacts;
    contacts.clear();

    collision_detection_general_BVH(_contact_body, _BVH_body, contacts);

    auto t1 = clock();

    VectorX fm_sub;
    computeForce(contacts, fm_sub, verbose);

    fm.segment(_contact_body->_index[0], 6) += fm_sub.head(6);
    fm.segment(_BVH_body->_index[0], 6) += fm_sub.tail(6);

    auto t2 = clock();

    _time.add("collision detection", t1 - t0);
    _time.add("compute force", t2 - t1);
}

void ForceGeneralBVHContact::computeForceWithDerivative(
    VectorX& fm, VectorX& fr, 
    MatrixX& Km, MatrixX& Dm, 
    MatrixX& Kr, MatrixX& Dr,
    bool verbose) {
    
    auto t0 = clock();

    // detect the contact points between two bodies
    std::vector<Contact> contacts;
    contacts.clear();
    collision_detection_general_BVH(_contact_body, _BVH_body, contacts);

    auto t1 = clock();

    VectorX fm_sub;
    MatrixX Km_sub, Dm_sub;
    computeForceWithDerivative(contacts, fm_sub, Km_sub, Dm_sub, verbose);

    fm.segment(_contact_body->_index[0], 6) += fm_sub.head(6);
    fm.segment(_BVH_body->_index[0], 6) += fm_sub.tail(6);

    Km.block(_contact_body->_index[0], _contact_body->_index[0], 6, 6) += Km_sub.topLeftCorner(6, 6);
    Km.block(_contact_body->_index[0], _BVH_body->_index[0], 6, 6) += Km_sub.topRightCorner(6, 6);
    Km.block(_BVH_body->_index[0], _contact_body->_index[0], 6, 6) += Km_sub.bottomLeftCorner(6, 6);
    Km.block(_BVH_body->_index[0], _BVH_body->_index[0], 6, 6) += Km_sub.bottomRightCorner(6, 6);
    
    Dm.block(_contact_body->_index[0], _contact_body->_index[0], 6, 6) += Dm_sub.topLeftCorner(6, 6);
    Dm.block(_contact_body->_index[0], _BVH_body->_index[0], 6, 6) += Dm_sub.topRightCorner(6, 6);
    Dm.block(_BVH_body->_index[0], _contact_body->_index[0], 6, 6) += Dm_sub.bottomLeftCorner(6, 6);
    Dm.block(_BVH_body->_index[0], _BVH_body->_index[0], 6, 6) += Dm_sub.bottomRightCorner(6, 6);

    auto t2 = clock();

    _time.add("collision detection", t1 - t0);
    _time.add("compute force", t2 - t1);
}

void ForceGeneralBVHContact::computeForceWithDerivative(
    VectorX& fm, VectorX& fr, 
    MatrixX& Km, MatrixX& Dm, 
    MatrixX& Kr, MatrixX& Dr,
    MatrixX& dfm_dp, MatrixX& dfr_dp,
    bool verbose) {

    auto t0 = clock();

    // detect the contact points between two bodies
    std::vector<Contact> contacts;
    contacts.clear();
    collision_detection_general_BVH(_contact_body, _BVH_body, contacts);

    auto t1 = clock();

    VectorX fm_sub;
    MatrixX Km_sub, Dm_sub;
    computeForceWithDerivative(contacts, fm_sub, Km_sub, Dm_sub, dfm_dp, verbose);

    fm.segment(_contact_body->_index[0], 6) += fm_sub.head(6);
    fm.segment(_BVH_body->_index[0], 6) += fm_sub.tail(6);

    Km.block(_contact_body->_index[0], _contact_body->_index[0], 6, 6) += Km_sub.topLeftCorner(6, 6);
    Km.block(_contact_body->_index[0], _BVH_body->_index[0], 6, 6) += Km_sub.topRightCorner(6, 6);
    Km.block(_BVH_body->_index[0], _contact_body->_index[0], 6, 6) += Km_sub.bottomLeftCorner(6, 6);
    Km.block(_BVH_body->_index[0], _BVH_body->_index[0], 6, 6) += Km_sub.bottomRightCorner(6, 6);
    
    Dm.block(_contact_body->_index[0], _contact_body->_index[0], 6, 6) += Dm_sub.topLeftCorner(6, 6);
    Dm.block(_contact_body->_index[0], _BVH_body->_index[0], 6, 6) += Dm_sub.topRightCorner(6, 6);
    Dm.block(_BVH_body->_index[0], _contact_body->_index[0], 6, 6) += Dm_sub.bottomLeftCorner(6, 6);
    Dm.block(_BVH_body->_index[0], _BVH_body->_index[0], 6, 6) += Dm_sub.bottomRightCorner(6, 6);

    auto t2 = clock();

    _time.add("collision detection", t1 - t0);
    _time.add("compute force", t2 - t1);
}

// new damping model (proportional to ddot * d) and friction related to damping force as well
void ForceGeneralBVHContact::computeForce(std::vector<Contact> &contacts, VectorX& fm, bool verbose) {
    Matrix3 R1 = _contact_body->_E_0i.topLeftCorner(3, 3);
    Vector6 phi1 = _contact_body->_phi;
    Matrix3 R2 = _BVH_body->_E_0i.topLeftCorner(3, 3);
    Vector6 phi2 = _BVH_body->_phi;
    
    fm = VectorX::Zero(12);

    dtype kn = _kn * _scale * _contact_body->_contact_scale;
    dtype kt = _kt * _scale * _contact_body->_contact_scale;
    for (int i = 0;i < contacts.size();i++) {
        Vector3 xi1 = contacts[i]._xi;
        Vector3 xw1 = contacts[i]._xw;
        Vector3 xw1_dot = R1 * (math::skew(xi1).transpose() * phi1.head(3) + phi1.tail(3));
        dtype d = contacts[i]._d;
        Vector3 dd_dx = contacts[i]._normal;
        
        dtype ddot;
        Vector3 n, tdot, xi2;
        _BVH_body->collision_parallel(xw1, xw1_dot, d, dd_dx, n, ddot, tdot, xi2);

        if (verbose) {
            std::cerr << "General BVH Contact: d = " << d << ", contact body = " << _contact_body->_name << ", BVH body = " << _BVH_body->_name << ", xi1 = " << xi1.transpose() << std::endl;
            std::cerr << "normal = " << n.transpose() << std::endl;
            std::cerr << "normal.norm = " << n.norm() << std::endl;
            std::cerr << "d = " << d << std::endl;
            std::cerr << "ddot = " << ddot << std::endl;
            std::cerr << "xw1 = " << xw1.transpose() << std::endl;
            std::cerr << "centor = " << _BVH_body->_E_0i.topRightCorner(3, 1).transpose() << std::endl;
        }
        
        MatrixX Gamma1 = math::gamma(xi1);
        MatrixX Gamma2 = math::gamma(xi2);

        MatrixX GTRT1 = Gamma1.transpose() * R1.transpose();
        Vector6 ni1 = GTRT1 * n;
        MatrixX GTRT2 = Gamma2.transpose() * R2.transpose();
        Vector6 ni2 = GTRT2 * n;

        // contact force
        Vector6 fc1 = (-kn * d + _damping * ddot * d) * ni1;
        fm.head(6) += fc1;
        Vector6 fc2 = (kn * d - _damping * ddot * d) * ni2;
        fm.tail(6) += fc2;

        if (verbose && dynamic_cast<BodyCuboid*>(const_cast<Body*>(_contact_body)) != nullptr) {
            std::cerr << "contact_body = " << _contact_body->_name << ", BVH body = " << _BVH_body->_name << ", xi1 = " << xi1.transpose() << ", fc1 = " << fc1.transpose() << std::endl;
        }

        if (verbose) {
            std::cerr << "n = " << n.transpose() << ", fc2 = " << fc2.transpose() << std::endl;
        }
        // frictional force
        dtype fc_norm = fc1.norm();
        dtype tdot_norm = tdot.norm();
        if (_mu > constants::eps) {
            // general contact body
            if (_mu * fc_norm >= kt * tdot_norm - constants::eps) {  // static frictional force
                Vector6 fs = -kt * GTRT1 * tdot;
                fm.head(6) += fs;
            } else {                                // dynamic frictional force
                Vector6 fd = -_mu * fc_norm / tdot_norm * GTRT1 * tdot;
                fm.head(6) += fd;
            }
            // BVH contact body
            if (_mu * fc_norm >= kt * tdot_norm - constants::eps) {  // static frictional force
                Vector6 fs = kt * GTRT2 * tdot;
                fm.tail(6) += fs;
            } else {                                // dynanmic frictional force
                Vector6 fd = _mu * fc_norm / tdot_norm * GTRT2 * tdot;
                fm.tail(6) += fd;
            }
        }
    }
}

void ForceGeneralBVHContact::computeForceWithDerivative(
    std::vector<Contact> &contacts,
    VectorX& fm, MatrixX& Km, MatrixX& Dm,
    bool verbose) {
    
    fm = VectorX::Zero(12);
    Km = MatrixX::Zero(12, 12);
    Dm = MatrixX::Zero(12, 12);

    Matrix3 R1 = _contact_body->_E_0i.topLeftCorner(3, 3);
    Vector6 phi1 = _contact_body->_phi;
    Matrix3 R2 = _BVH_body->_E_0i.topLeftCorner(3, 3);
    Vector6 phi2 = _BVH_body->_phi;

    dtype kn = _kn * _scale * _contact_body->_contact_scale;
    dtype kt = _kt * _scale * _contact_body->_contact_scale;
    for (int i = 0;i < contacts.size();i++) {
        Vector3 xi1 = contacts[i]._xi;
        Vector3 xw1 = contacts[i]._xw;
        Vector3 xw1_dot = R1 * (math::skew(xi1).transpose() * phi1.head(3) + phi1.tail(3));
        dtype d = contacts[i]._d;
        Vector3 dd_dx = contacts[i]._normal;
        
        Matrix36 dxw1_dq1;
        dxw1_dq1.leftCols(3) = -R1 * math::skew(xi1);
        dxw1_dq1.rightCols(3) = R1;
        
        Matrix36 dxw1dot_dq1;
        dxw1dot_dq1.leftCols(3) = -R1 * math::skew(math::skew(xi1).transpose() * phi1.head(3) + phi1.tail(3));
        dxw1dot_dq1.rightCols(3).setZero();

        Matrix36 dxw1dot_dphi1;
        dxw1dot_dphi1.leftCols(3) = R1 * math::skew(xi1).transpose();
        dxw1dot_dphi1.rightCols(3) = R1;

        dtype ddot;
        Vector3 n, tdot, xi2;
        RowVector3 dd_dxw1, dddot_dxw1, dddot_dxw1dot;
        RowVector6 dd_dq2, dddot_dq2, dddot_dphi2;
        Matrix3 dn_dxw1, dtdot_dxw1, dtdot_dxw1dot, dxi2_dxw1;
        Matrix36 dn_dq2, dtdot_dq2, dtdot_dphi2, dxi2_dq2;
        _BVH_body->collision_parallel(xw1, xw1_dot, d, dd_dx, n, ddot, tdot, xi2,
                                    dd_dxw1, dd_dq2,
                                    dn_dxw1, dn_dq2,
                                    dddot_dxw1, dddot_dxw1dot,
                                    dddot_dq2, dddot_dphi2,
                                    dtdot_dxw1, dtdot_dxw1dot,
                                    dtdot_dq2, dtdot_dphi2,
                                    dxi2_dxw1, dxi2_dq2);

        Matrix36 dxi2_dq1 = dxi2_dxw1 * dxw1_dq1;

        MatrixX Gamma1 = math::gamma(xi1);
        MatrixX Gamma2 = math::gamma(xi2);
        
        // derivatives
        JacobianMatrixVector dGamma2_dxi2(3, 6, 3);
        dGamma2_dxi2.setZero();
        dGamma2_dxi2(0)(1, 2) = 1; dGamma2_dxi2(0)(2, 1) = -1;
        dGamma2_dxi2(1)(0, 2) = -1; dGamma2_dxi2(1)(2, 0) = 1;
        dGamma2_dxi2(2)(0, 1) = 1; dGamma2_dxi2(2)(1, 0) = -1;

        JacobianMatrixVector dGamma2_dq1(3, 6, 6);
        dGamma2_dq1.setZero();
        for (int j = 0;j < 6;j++)
            for (int k = 0;k < 3;k++)
                dGamma2_dq1(j) += dGamma2_dxi2(k) * dxi2_dq1(k, j);
        
        JacobianMatrixVector dGamma2_dq2(3, 6, 6);
        dGamma2_dq2.setZero();
        for (int j = 0;j < 6;j++)
            for (int k = 0;k < 3;k++)
                dGamma2_dq2(j) += dGamma2_dxi2(k) * dxi2_dq2(k, j);

        MatrixX GTRT1 = Gamma1.transpose() * R1.transpose();
        Vector6 ni1 = GTRT1 * n;
        MatrixX GTRT2 = Gamma2.transpose() * R2.transpose();
        Vector6 ni2 = GTRT2 * n;

        // contact force
        // contact force on general contact body
        Vector6 fc1 = (-kn * d + _damping * ddot * d) * ni1;

        fm.head(6) += fc1;

        // derivatives
        RowVector6 dd_dq1 = dd_dxw1 * dxw1_dq1;
        RowVector6 dddot_dq1 = dddot_dxw1 * dxw1_dq1 + dddot_dxw1dot * dxw1dot_dq1;
        Matrix36 dn_dq1 = dn_dxw1 * dxw1_dq1;
        RowVector6 dddot_dphi1 = dddot_dxw1dot * dxw1dot_dphi1;
        Matrix36 dtdot_dq1 = dtdot_dxw1 * dxw1_dq1 + dtdot_dxw1dot * dxw1dot_dq1;
        Matrix36 dtdot_dphi1 = dtdot_dxw1dot * dxw1dot_dphi1;

        Matrix6 dfc1_dq1 = - ni1 * (kn * dd_dq1 - _damping * dddot_dq1 * d - _damping * ddot * dd_dq1) -
                            (kn * d - _damping * ddot * d) * GTRT1 * dn_dq1;
        dfc1_dq1.leftCols(3) -= (kn * d - _damping * ddot * d) * Gamma1.transpose() * math::skew(R1.transpose() * n);
        Km.topLeftCorner(6, 6) += dfc1_dq1;

        Matrix6 dfc1_dq2 = - ni1 * (kn * dd_dq2 - _damping * dddot_dq2 * d - _damping * ddot * dd_dq2) - 
                            (kn * d - _damping * ddot * d) * GTRT1 * dn_dq2;
        Km.topRightCorner(6, 6) += dfc1_dq2;

        Matrix6 dfc1_dphi1 = ni1 * _damping * dddot_dphi1 * d;
        Dm.topLeftCorner(6, 6) += dfc1_dphi1;

        Matrix6 dfc1_dphi2 = ni1 * _damping * dddot_dphi2 * d;
        Dm.topRightCorner(6, 6) += dfc1_dphi2;

        // contact force on BVH contact body
        Vector6 fc2 = (kn * d - _damping * ddot * d) * ni2;

        fm.tail(6) += fc2;

        // derivatives
        Matrix6 tmp1;
        tmp1.setZero();
        for (int j = 0;j < 6;j++)
            tmp1.col(j) = dGamma2_dq1(j).transpose() * (R2.transpose() * n);
        Km.bottomLeftCorner(6, 6) += ni2 * (kn * dd_dq1 - _damping * dddot_dq1 * d - _damping * ddot * dd_dq1) + 
                                        (kn * d - _damping * ddot * d) * (tmp1 + GTRT2 * dn_dq1);
        Matrix6 tmp2;
        tmp2.setZero();
        for (int j = 0;j < 6;j++)
            tmp2.col(j) = dGamma2_dq2(j).transpose() * (R2.transpose() * n);
        Km.bottomRightCorner(6, 6) += ni2 * (kn * dd_dq2 - _damping * dddot_dq2 * d - _damping * ddot * dd_dq2) + 
                                        (kn * d - _damping * ddot * d) * (tmp2 + GTRT2 * dn_dq2);
        Km.block(6, 6, 6, 3) +=
            (kn * d - _damping * ddot * d) * Gamma2.transpose() * math::skew(R2.transpose() * n);

        Dm.bottomLeftCorner(6, 6) -= ni2 * _damping * dddot_dphi1 * d;
        Dm.bottomRightCorner(6, 6) -= ni2 * _damping * dddot_dphi2 * d;

        // frictional force
        if (_mu > constants::eps) {
            dtype fc_norm = fc1.norm();
            dtype tdot_norm = tdot.norm();
            RowVector6 dfc_norm_dq1 = (1. / fc_norm) * (fc1.transpose() * dfc1_dq1);
            RowVector6 dfc_norm_dq2 = (1. / fc_norm) * (fc1.transpose() * dfc1_dq2);
            RowVector6 dfc_norm_dphi1 = (1. / fc_norm) * (fc1.transpose() * dfc1_dphi1);
            RowVector6 dfc_norm_dphi2 = (1. / fc_norm) * (fc1.transpose() * dfc1_dphi2);
            RowVector6 dtdot_norm_dq1 = (1. / tdot_norm) * (tdot.transpose() * dtdot_dq1);
            RowVector6 dtdot_norm_dq2 = (1. / tdot_norm) * (tdot.transpose() * dtdot_dq2);
            RowVector6 dtdot_norm_dphi1 = (1. / tdot_norm) * (tdot.transpose() * dtdot_dphi1);
            RowVector6 dtdot_norm_dphi2 = (1. / tdot_norm) * (tdot.transpose() * dtdot_dphi2);

            // frictional force on the general contact body
            if (_mu * fc_norm >= kt * tdot_norm - constants::eps) {  // static frictional force
                // std::cerr << "static" << std::endl;
                Vector6 fs = -kt * GTRT1 * tdot;
                fm.head(6) += fs;

                // derivatives
                Km.topLeftCorner(6, 6) -= kt * GTRT1 * dtdot_dq1;
                Km.block(0, 0, 6, 3) -= kt * Gamma1.transpose() * math::skew(R1.transpose() * tdot);
                Km.topRightCorner(6, 6) -= kt * GTRT1 * dtdot_dq2;

                Dm.topLeftCorner(6, 6) -= kt * GTRT1 * dtdot_dphi1;
                Dm.topRightCorner(6, 6) -= kt * GTRT1 * dtdot_dphi2;
            } else {                                // dynamic frictional force
                // std::cerr << "dynamic" << std::endl;
                Vector6 fd = -_mu * fc_norm / tdot_norm * GTRT1 * tdot;
                fm.head(6) += fd;

                // derivatives
                Km.topLeftCorner(6, 6) -=
                    _mu / tdot_norm * GTRT1 * (tdot * dfc_norm_dq1 - fc_norm / tdot_norm * tdot * dtdot_norm_dq1 + fc_norm * dtdot_dq1);
                Km.block(0, 0, 6, 3) -=
                    _mu * Gamma1.transpose() * math::skew(fc_norm / tdot_norm * R1.transpose() * tdot);
                Km.topRightCorner(6, 6) -=
                    _mu / tdot_norm * GTRT1 * (tdot * dfc_norm_dq2 - fc_norm / tdot_norm * tdot * dtdot_norm_dq2 + fc_norm * dtdot_dq2);

                Dm.topLeftCorner(6, 6) -=
                    _mu / tdot_norm * GTRT1 * (tdot * dfc_norm_dphi1 - fc_norm / tdot_norm * tdot * dtdot_norm_dphi1 + fc_norm * dtdot_dphi1);
                Dm.topRightCorner(6, 6) -=
                    _mu / tdot_norm * GTRT1 * (tdot * dfc_norm_dphi2 - fc_norm / tdot_norm * tdot * dtdot_norm_dphi2 + fc_norm * dtdot_dphi2);
            }

            // frictional force on the BVH contact body
            if (_mu * fc_norm >= kt * tdot_norm - constants::eps) {  // static frictional force
                Vector6 fs = kt * GTRT2 * tdot;
                fm.tail(6) += fs;

                // derivatives
                Matrix6 tmp1;
                tmp1.setZero();
                for (int j = 0;j < 6;j++)
                    tmp1.col(j) = dGamma2_dq1(j).transpose() * (R2.transpose() * tdot);
                Km.bottomLeftCorner(6, 6) += kt * GTRT2 * dtdot_dq1 + kt * tmp1;
                
                Matrix6 tmp2;
                tmp2.setZero();
                for (int j = 0;j < 6;j++)
                    tmp2.col(j) = dGamma2_dq2(j).transpose() * (R2.transpose() * tdot);
                Km.bottomRightCorner(6, 6) += kt * GTRT2 * dtdot_dq2 + kt * tmp2;
                Km.block(6, 6, 6, 3) += kt * Gamma2.transpose() * math::skew(R2.transpose() * tdot);
                
                Dm.bottomLeftCorner(6, 6) += kt * GTRT2 * dtdot_dphi1;
                Dm.bottomRightCorner(6, 6) += kt * GTRT2 * dtdot_dphi2;
            } else {                                // dynanmic frictional force
                Vector6 fd = _mu * fc_norm / tdot_norm * GTRT2 * tdot;
                fm.tail(6) += fd;

                // derivatives
                Matrix6 tmp1;
                tmp1.setZero();
                for (int j = 0;j < 6;j++)
                    tmp1.col(j) = dGamma2_dq1(j).transpose() * (R2.transpose() * tdot);
                Km.bottomLeftCorner(6, 6) += 
                    _mu / tdot_norm * GTRT2 * (tdot * dfc_norm_dq1 - tdot * fc_norm / tdot_norm * dtdot_norm_dq1 + fc_norm * dtdot_dq1) +
                    _mu * fc_norm / tdot_norm * tmp1;
                
                Matrix6 tmp2;
                tmp2.setZero();
                for (int j = 0;j < 6;j++)
                    tmp2.col(j) = dGamma2_dq2(j).transpose() * (R2.transpose() * tdot);
                Km.bottomRightCorner(6, 6) += 
                    _mu / tdot_norm * GTRT2 * (tdot * dfc_norm_dq2 - tdot * fc_norm / tdot_norm * dtdot_norm_dq2 + fc_norm * dtdot_dq2) +
                    _mu * fc_norm / tdot_norm * tmp2;
                Km.block(6, 6, 6, 3) += 
                    _mu * Gamma2.transpose() * math::skew(fc_norm / tdot_norm * R2.transpose() * tdot);

                Dm.bottomLeftCorner(6, 6) += 
                    _mu / tdot_norm * GTRT2 * (tdot * dfc_norm_dphi1 - tdot * fc_norm / tdot_norm * dtdot_norm_dphi1 + fc_norm * dtdot_dphi1);
                Dm.bottomRightCorner(6, 6) += 
                    _mu / tdot_norm * GTRT2 * (tdot * dfc_norm_dphi2 - tdot * fc_norm / tdot_norm * dtdot_norm_dphi2 + fc_norm * dtdot_dphi2);
            }
        }
    }
}

void ForceGeneralBVHContact::computeForceWithDerivative(
    std::vector<Contact> &contacts,
    VectorX& fm, MatrixX& Km, MatrixX& Dm, MatrixX& dfm_dp, 
    bool verbose) {
    
    fm = VectorX::Zero(12);
    Km = MatrixX::Zero(12, 12);
    Dm = MatrixX::Zero(12, 12);

    Matrix3 R1 = _contact_body->_E_0i.topLeftCorner(3, 3);
    Vector6 phi1 = _contact_body->_phi;
    Vector3 w1_dot = phi1.head(3);
    Vector3 v1_dot = phi1.tail(3);
    Matrix3 R2 = _BVH_body->_E_0i.topLeftCorner(3, 3);
    Vector6 phi2 = _BVH_body->_phi;

    // tools for design params 3
    Matrix3 dxw1_dp3 = R1, dxw1dot_dp3 = MatrixX::Zero(3, 3);
    for (int k = 0;k < 3;k++) {
        dxw1dot_dp3.col(k) = R1 * math::skew(-Vector3::Unit(k)) * w1_dot;
    }

    if (verbose)
        std::cerr << "num general contacts = " << contacts.size() << std::endl;

    dtype kn = _kn * _scale * _contact_body->_contact_scale;
    dtype kt = _kt * _scale * _contact_body->_contact_scale;
    for (int i = 0;i < contacts.size();i++) {
        Vector3 xi1 = contacts[i]._xi;
        Vector3 xw1 = contacts[i]._xw;
        Matrix3 xi1_brac = math::skew(xi1);
        Vector3 xi1_dot = xi1_brac.transpose() * phi1.head(3) + phi1.tail(3);
        Vector3 xw1_dot = R1 * xi1_dot;
        int contact_id = contacts[i]._id;
        dtype d = contacts[i]._d;
        Vector3 dd_dx = contacts[i]._normal;

        Matrix36 dxw1_dq1;
        dxw1_dq1.leftCols(3) = -R1 * math::skew(xi1);
        dxw1_dq1.rightCols(3) = R1;
        
        Matrix36 dxw1dot_dq1;
        dxw1dot_dq1.leftCols(3) = -R1 * math::skew(math::skew(xi1).transpose() * phi1.head(3) + phi1.tail(3));
        dxw1dot_dq1.rightCols(3).setZero();

        Matrix36 dxw1dot_dphi1;
        dxw1dot_dphi1.leftCols(3) = R1 * math::skew(xi1).transpose();
        dxw1dot_dphi1.rightCols(3) = R1;

        dtype ddot;
        Vector3 n, tdot, xi2;
        RowVector3 dd_dxw1, dddot_dxw1, dddot_dxw1dot;
        RowVector6 dd_dq2, dddot_dq2, dddot_dphi2;
        Matrix3 dn_dxw1, dtdot_dxw1, dtdot_dxw1dot, dxi2_dxw1;
        Matrix36 dn_dq2, dtdot_dq2, dtdot_dphi2, dxi2_dq2;
        _BVH_body->collision_parallel(xw1, xw1_dot, d, dd_dx, n, ddot, tdot, xi2,
                                    dd_dxw1, dd_dq2,
                                    dn_dxw1, dn_dq2,
                                    dddot_dxw1, dddot_dxw1dot,
                                    dddot_dq2, dddot_dphi2,
                                    dtdot_dxw1, dtdot_dxw1dot,
                                    dtdot_dq2, dtdot_dphi2,
                                    dxi2_dxw1, dxi2_dq2);

        Matrix36 dxi2_dq1 = dxi2_dxw1 * dxw1_dq1;

        MatrixX Gamma1 = math::gamma(xi1);
        MatrixX Gamma2 = math::gamma(xi2);
        
        // derivatives
        JacobianMatrixVector dGamma2_dxi2(3, 6, 3);
        dGamma2_dxi2.setZero();
        dGamma2_dxi2(0)(1, 2) = 1; dGamma2_dxi2(0)(2, 1) = -1;
        dGamma2_dxi2(1)(0, 2) = -1; dGamma2_dxi2(1)(2, 0) = 1;
        dGamma2_dxi2(2)(0, 1) = 1; dGamma2_dxi2(2)(1, 0) = -1;

        JacobianMatrixVector dGamma2_dq1(3, 6, 6);
        dGamma2_dq1.setZero();
        for (int j = 0;j < 6;j++)
            for (int k = 0;k < 3;k++)
                dGamma2_dq1(j) += dGamma2_dxi2(k) * dxi2_dq1(k, j);
        
        JacobianMatrixVector dGamma2_dq2(3, 6, 6);
        dGamma2_dq2.setZero();
        for (int j = 0;j < 6;j++)
            for (int k = 0;k < 3;k++)
                dGamma2_dq2(j) += dGamma2_dxi2(k) * dxi2_dq2(k, j);

        MatrixX GTRT1 = Gamma1.transpose() * R1.transpose();
        Vector6 ni1 = GTRT1 * n;
        Vector3 n1 = R1.transpose() * n;
        MatrixX GTRT2 = Gamma2.transpose() * R2.transpose();
        Vector6 ni2 = GTRT2 * n;
        Vector3 n2 = R2.transpose() * n;

        // contact force
        // contact force on general contact body
        Vector6 fc1 = (-kn * d + _damping * ddot * d) * ni1;

        fm.head(6) += fc1;

        // derivatives
        // tools
        RowVector6 dd_dq1 = dd_dxw1 * dxw1_dq1;
        RowVector6 dddot_dq1 = dddot_dxw1 * dxw1_dq1 + dddot_dxw1dot * dxw1dot_dq1;
        Matrix36 dn_dq1 = dn_dxw1 * dxw1_dq1;
        RowVector6 dddot_dphi1 = dddot_dxw1dot * dxw1dot_dphi1;
        Matrix36 dtdot_dq1 = dtdot_dxw1 * dxw1_dq1 + dtdot_dxw1dot * dxw1dot_dq1;
        Matrix36 dtdot_dphi1 = dtdot_dxw1dot * dxw1dot_dphi1;
        // tools for design derivatives
        MatrixX dxw1_dp1 = MatrixX::Zero(3, 12), dxw1_dp2 = MatrixX::Zero(3, 12);
        MatrixX dxw1dot_dp1 = MatrixX::Zero(3, 12), dxw1dot_dp2 = MatrixX::Zero(3, 12);
        RowVectorX dd_dp1 = RowVectorX::Zero(_sim->_ndof_p1), dd_dp2 = RowVectorX::Zero(12), dd_dp3 = RowVectorX::Zero(3);
        RowVectorX dddot_dp1 = RowVectorX::Zero(_sim->_ndof_p1), dddot_dp2 = RowVectorX::Zero(12), dddot_dp3 = RowVectorX::Zero(3);
        MatrixX dn_dp1 = MatrixX::Zero(3, _sim->_ndof_p1), dn_dp2 = MatrixX::Zero(3, 12), dn_dp3 = MatrixX::Zero(3, 3);
        MatrixX dxi2_dp1 = MatrixX::Zero(3, _sim->_ndof_p1), dxi2_dp2 = MatrixX::Zero(3, 12), dxi2_dp3 = MatrixX::Zero(3, 3);
        MatrixX dtdot_dp1 = MatrixX::Zero(3, _sim->_ndof_p1), dtdot_dp2 = MatrixX::Zero(3, 12), dtdot_dp3 = MatrixX::Zero(3, 3);

        // dxw1_dp1 and dxw1dot_dp1
        // design params 1
        for (auto ancestor = _contact_body->_joint; ancestor != nullptr; ancestor = ancestor->_parent) {
            if (ancestor->_design_params_1._active) {
                for (int k = 0;k < ancestor->_design_params_1._ndof;k++) {
                    int idx = ancestor->_design_params_1._param_index(k);
                    Matrix3 dR_dp1 = _contact_body->_dE0i_dp1(idx).topLeftCorner(3, 3);
                    Vector3 dp_dp1 = _contact_body->_dE0i_dp1(idx).topRightCorner(3, 1);
                    Vector6 dphi_dp1 = _contact_body->_dphi_dp1.col(idx);
                    dxw1_dp1.col(k) = dR_dp1 * xi1 + dp_dp1;
                    dxw1dot_dp1.col(k) = dR_dp1 * xi1_dot + R1 * (xi1_brac.transpose() * dphi_dp1.head(3) + dphi_dp1.tail(3));
                }
                int idx = ancestor->_design_params_1._param_index(0);
                int ndof = ancestor->_design_params_1._ndof;
                dd_dp1.segment(idx, ndof) = dd_dxw1 * dxw1_dp1;
                dddot_dp1.segment(idx, ndof) = dddot_dxw1 * dxw1_dp1 + dddot_dxw1dot * dxw1dot_dp1;
                dn_dp1.middleCols(idx, ndof) = dn_dxw1 * dxw1_dp1;
                dxi2_dp1.middleCols(idx, ndof) = dxi2_dxw1 * dxw1_dp1;
                dtdot_dp1.middleCols(idx, ndof) = dtdot_dxw1 * dxw1_dp1 + dtdot_dxw1dot * dxw1dot_dp1;
            }
        }
        // design params 2
        if (_contact_body->_design_params_2._active) {
            for (int k = 0;k < _contact_body->_design_params_2._ndof;k++) {
                Matrix3 dR_dp2 = _contact_body->_dE0i_dp2(k).topLeftCorner(3, 3);
                Vector3 dp_dp2 = _contact_body->_dE0i_dp2(k).topRightCorner(3, 1);
                Vector6 dphi_dp2 = _contact_body->_dphi_dp2.col(k);
                dxw1_dp2.col(k) = dR_dp2 * xi1 + dp_dp2;
                dxw1dot_dp2.col(k) = dR_dp2 * xi1_dot + R1 * (xi1_brac.transpose() * dphi_dp2.head(3) + dphi_dp2.tail(3));
            }
            dd_dp2 = dd_dxw1 * dxw1_dp2;
            dddot_dp2 = dddot_dxw1 * dxw1_dp2 + dddot_dxw1dot * dxw1dot_dp2;
            dn_dp2 = dn_dxw1 * dxw1_dp2;
            dxi2_dp2 = dxi2_dxw1 * dxw1_dp2;
            dtdot_dp2 = dtdot_dxw1 * dxw1_dp2 + dtdot_dxw1dot * dxw1dot_dp2;
        }
        // design params 3
        if (_contact_body->_design_params_3._active) {
            dd_dp3 = dd_dxw1 * dxw1_dp3;
            dddot_dp3 = dddot_dxw1 * dxw1_dp3 + dddot_dxw1dot * dxw1dot_dp3;
            dn_dp3 = dn_dxw1 * dxw1_dp3;
            dxi2_dp3 = dxi2_dxw1 * dxw1_dp3;
            dtdot_dp3 = dtdot_dxw1 * dxw1_dp3 + dtdot_dxw1dot * dxw1dot_dp3;
        }
        JacobianMatrixVector dG1_dp3(3, 6, 3), dG2_dp1(3, 6, _sim->_ndof_p1), dG2_dp2(3, 6, 12), dG2_dp3(3, 6, 3);
        for (int k = 0;k < 3;k++) {
            dG1_dp3(k).leftCols(3) = math::skew(Vector3::Unit(k)).transpose();
        }
        for (auto ancestor = _contact_body->_joint; ancestor != nullptr; ancestor = ancestor->_parent) {
            if (ancestor->_design_params_1._active) {
                for (int k = 0;k < ancestor->_design_params_1._ndof;k++) {
                    int idx = ancestor->_design_params_1._param_index(k);
                    for (int l = 0;l < 3;l++)   {
                        dG2_dp1(idx) += dG1_dp3(l) * dxi2_dp1(l, idx);
                    }
                }
            }
        }
        if (_contact_body->_design_params_2._active) {
            for (int k = 0;k < _contact_body->_design_params_2._ndof;k++) {
                for (int l = 0;l < 3;l++) {
                    dG2_dp2(k) += dG1_dp3(l) * dxi2_dp2(l, k);
                }
            }
        }
        if (_contact_body->_design_params_3._active) {
            for (int k = 0;k < 3;k++)
                for (int l = 0;l < 3;l++) {
                    dG2_dp3(k) += dG1_dp3(l) * dxi2_dp3(l, k);
                }
        }

        // dfc1_dq1
        Matrix6 dfc1_dq1 = - ni1 * (kn * dd_dq1 - _damping * dddot_dq1 * d - _damping * ddot * dd_dq1) -
                            (kn * d - _damping * ddot * d) * GTRT1 * dn_dq1;
        dfc1_dq1.leftCols(3) -= (kn * d - _damping * ddot * d) * Gamma1.transpose() * math::skew(R1.transpose() * n);
        Km.topLeftCorner(6, 6) += dfc1_dq1;

        // dfc1_dq2
        Matrix6 dfc1_dq2 = - ni1 * (kn * dd_dq2 - _damping * dddot_dq2 * d - _damping * ddot * dd_dq2) - 
                            (kn * d - _damping * ddot * d) * GTRT1 * dn_dq2;
        Km.topRightCorner(6, 6) += dfc1_dq2;

        // dfc1_dphi1
        Matrix6 dfc1_dphi1 = ni1 * _damping * dddot_dphi1 * d;
        Dm.topLeftCorner(6, 6) += dfc1_dphi1;

        // dfc1_dphi2
        Matrix6 dfc1_dphi2 = ni1 * _damping * dddot_dphi2 * d;
        Dm.topRightCorner(6, 6) += dfc1_dphi2;

        // design derivatives
        // dfc1_dp
        MatrixX dfc1_dp1 = MatrixX::Zero(6, _sim->_ndof_p1), dfc1_dp2 = MatrixX::Zero(6, 12), dfc1_dp3 = MatrixX::Zero(6, 3), dfc1_dp6 = VectorX::Zero(6);
        // design params 1
        for (auto ancestor = _contact_body->_joint; ancestor != nullptr; ancestor = ancestor->_parent) {
            if (ancestor->_design_params_1._active) {
                for (int k = 0;k < ancestor->_design_params_1._ndof;k++) {
                    int idx = ancestor->_design_params_1._param_index(k);
                    Matrix3 dR_dp1 = _contact_body->_dE0i_dp1(idx).topLeftCorner(3, 3);
                    dfc1_dp1.col(idx) = 
                        -(kn * dd_dp1(idx) - _damping * dddot_dp1(idx) * d - _damping * ddot * dd_dp1(idx)) * ni1
                        -(kn * d - _damping * ddot * d) * (Gamma1.transpose() * dR_dp1.transpose() * n + GTRT1 * dn_dp1.col(idx));
                    dfm_dp.block(_contact_body->_index[0], idx, 6, 1) += dfc1_dp1.col(idx);
                }
            }
        }
        // design params 2
        if (_contact_body->_design_params_2._active) {
            for (int k = 0;k < _contact_body->_design_params_2._ndof;k++) {
                int idx = _contact_body->_design_params_2._param_index(k) + _sim->_ndof_p1;
                Matrix3 dR_dp2 = _contact_body->_dE0i_dp2(k).topLeftCorner(3, 3);
                dfc1_dp2.col(k) = 
                    -(kn * dd_dp2(k) - _damping * dddot_dp2(k) * d - _damping * ddot * dd_dp2(k)) * ni1
                    -(kn * d - _damping * ddot * d) * (Gamma1.transpose() * dR_dp2.transpose() * n + GTRT1 * dn_dp2.col(k));
                dfm_dp.block(_contact_body->_index[0], idx, 6, 1) += dfc1_dp2.col(k);
            }
        }
        // design params 3
        if (_contact_body->_design_params_3._active) {
            for (int k = 0;k < 3;k++) {
                int idx = _contact_body->_design_params_3._param_index(contact_id * 3 + k) + _sim->_ndof_p1 + _sim->_ndof_p2;
                dfc1_dp3.col(k) = 
                    -(kn * dd_dp3(k) - _damping * dddot_dp3(k) * d - _damping * ddot * dd_dp3(k)) * ni1
                    -(kn * d - _damping * ddot * d) * (dG1_dp3(k).transpose() * n1 + GTRT1 * dn_dp3.col(k));
                dfm_dp.block(_contact_body->_index[0], idx, 6, 1) += dfc1_dp3.col(k);
            }
        }
        // design params 6
        int p6_offset = _sim->_ndof_p1 + _sim->_ndof_p2 + _sim->_ndof_p3 + _sim->_ndof_p4 + _sim->_ndof_p5;
        if (_contact_body->_design_params_6._active) {
            int idx = _contact_body->_design_params_6._param_index(0) + p6_offset;
            dfc1_dp6 = -kn * d * ni1 / _contact_body->_contact_scale;
            dfm_dp.block(_contact_body->_index[0], idx, 6, 1) += dfc1_dp6;
        }

        // contact force on BVH contact body
        Vector6 fc2 = (kn * d - _damping * ddot * d) * ni2;

        fm.tail(6) += fc2;

        // derivatives
        Matrix6 tmp1;
        tmp1.setZero();
        for (int j = 0;j < 6;j++)
            tmp1.col(j) = dGamma2_dq1(j).transpose() * (R2.transpose() * n);
        Km.bottomLeftCorner(6, 6) += ni2 * (kn * dd_dq1 - _damping * dddot_dq1 * d - _damping * ddot * dd_dq1) + 
                                        (kn * d - _damping * ddot * d) * (tmp1 + GTRT2 * dn_dq1);
        Matrix6 tmp2;
        tmp2.setZero();
        for (int j = 0;j < 6;j++)
            tmp2.col(j) = dGamma2_dq2(j).transpose() * (R2.transpose() * n);
        Km.bottomRightCorner(6, 6) += ni2 * (kn * dd_dq2 - _damping * dddot_dq2 * d - _damping * ddot * dd_dq2) + 
                                        (kn * d - _damping * ddot * d) * (tmp2 + GTRT2 * dn_dq2);
        Km.block(6, 6, 6, 3) +=
            (kn * d - _damping * ddot * d) * Gamma2.transpose() * math::skew(R2.transpose() * n);

        Dm.bottomLeftCorner(6, 6) -= ni2 * _damping * dddot_dphi1 * d;
        Dm.bottomRightCorner(6, 6) -= ni2 * _damping * dddot_dphi2 * d;

        // design derivatives
        // dfc2_dp
        // design params 1
        for (auto ancestor = _contact_body->_joint; ancestor != nullptr; ancestor = ancestor->_parent) {
            if (ancestor->_design_params_1._active) {
                int idx_start = ancestor->_design_params_1._param_index(0);
                int ndof = ancestor->_design_params_1._ndof;
                for (int k = 0;k < ndof;k++) {
                    int idx = idx_start + k;
                    dfm_dp.block(_BVH_body->_index[0], idx, 6, 1) += 
                        (kn * dd_dp1(idx) - _damping * dddot_dp1(idx) * d - _damping * ddot * dd_dp1(idx)) * ni2
                        + (kn * d - _damping * ddot * d) * (dG2_dp1(idx).transpose() * n2 + GTRT2 * dn_dp1.col(idx));
                }
            }
        }
        // design params 2
        if (_contact_body->_design_params_2._active) {
            int idx_start = _contact_body->_design_params_2._param_index(0);
            int ndof = _contact_body->_design_params_2._ndof;
            for (int k = 0;k < ndof;k++) {
                int idx = idx_start + k + _sim->_ndof_p1;
                dfm_dp.block(_BVH_body->_index[0], idx, 6, 1) += 
                    (kn * dd_dp2(k) - _damping * dddot_dp2(k) * d - _damping * ddot * dd_dp2(k)) * ni2
                    + (kn * d - _damping * ddot * d) * (dG2_dp2(k).transpose() * n2 + GTRT2 * dn_dp2.col(k));
            }
        }
        // design params 3
        if (_contact_body->_design_params_3._active) {
            int idx_start = _contact_body->_design_params_3._param_index(contact_id * 3);
            for (int k = 0;k < 3;k++) {
                int idx = idx_start + k + _sim->_ndof_p1 + _sim->_ndof_p2;
                dfm_dp.block(_BVH_body->_index[0], idx, 6, 1) += 
                    (kn * dd_dp3(k) - _damping * dddot_dp3(k) * d - _damping * ddot * dd_dp3(k)) * ni2
                    + (kn * d - _damping * ddot * d) * (dG2_dp3(k).transpose() * n2 + GTRT2 * dn_dp3.col(k));
            }
        }
        // design params 6
        if (_contact_body->_design_params_6._active) {
            int idx = _contact_body->_design_params_6._param_index(0) + p6_offset;
            dfm_dp.block(_BVH_body->_index[0], idx, 6, 1) += kn * d * ni2 / _contact_body->_contact_scale;
        }

        // frictional force
        if (_mu > constants::eps) {
            dtype fc_norm = fc1.norm();
            dtype tdot_norm = tdot.norm();

            // frictional force on the general contact body
            if (_mu * fc_norm >= kt * tdot_norm - constants::eps) {  // static frictional force
                // general contact body
                Vector6 fs1 = -kt * GTRT1 * tdot;
                fm.head(6) += fs1;

                // derivatives
                Km.topLeftCorner(6, 6) -= kt * GTRT1 * dtdot_dq1;
                Km.block(0, 0, 6, 3) -= kt * Gamma1.transpose() * math::skew(R1.transpose() * tdot);
                Km.topRightCorner(6, 6) -= kt * GTRT1 * dtdot_dq2;

                Dm.topLeftCorner(6, 6) -= kt * GTRT1 * dtdot_dphi1;
                Dm.topRightCorner(6, 6) -= kt * GTRT1 * dtdot_dphi2;

                // design derivatives
                // design params 1
                for (auto ancestor = _contact_body->_joint; ancestor != nullptr; ancestor = ancestor->_parent) {
                    if (ancestor->_design_params_1._active) {
                        for (int k = 0;k < ancestor->_design_params_1._ndof;k++) {
                            int idx = ancestor->_design_params_1._param_index(k);
                            Matrix3 dR_dp1 = _contact_body->_dE0i_dp1(idx).topLeftCorner(3, 3);
                            dfm_dp.block(_contact_body->_index[0], idx, 6, 1) += 
                                -kt * Gamma1.transpose() * (dR_dp1.transpose() * tdot + R1.transpose() * dtdot_dp1.col(idx));
                        }
                    }
                }
                // design params 2
                if (_contact_body->_design_params_2._active) {
                    for (int k = 0;k < _contact_body->_design_params_2._ndof;k++) {
                        int idx = _contact_body->_design_params_2._param_index(k) + _sim->_ndof_p1;
                        Matrix3 dR_dp2 = _contact_body->_dE0i_dp2(k).topLeftCorner(3, 3);
                        dfm_dp.block(_contact_body->_index[0], idx, 6, 1) += 
                            -kt * Gamma1.transpose() * (dR_dp2.transpose() * tdot + R1.transpose() * dtdot_dp2.col(k));
                    }
                }
                // design params 3
                if (_contact_body->_design_params_3._active) {
                    for (int k = 0;k < 3;k++) {
                        int idx = _contact_body->_design_params_3._param_index(contact_id * 3 + k) + _sim->_ndof_p1 + _sim->_ndof_p2;
                        dfm_dp.block(_contact_body->_index[0], idx, 6, 1) += 
                            - kt * dG1_dp3(k).transpose() * R1.transpose() * tdot - kt * GTRT1 * dtdot_dp3.col(k);
                    }
                }
                // design params 6
                if (_contact_body->_design_params_6._active) {
                    int idx = _contact_body->_design_params_6._param_index(0) + p6_offset;
                    dfm_dp.block(_contact_body->_index[0], idx, 6, 1) += fs1 / _contact_body->_contact_scale;
                }

                // BVH contact body
                Vector6 fs2 = kt * GTRT2 * tdot;
                fm.tail(6) += fs2;

                // derivatives
                Matrix6 tmp1;
                tmp1.setZero();
                for (int j = 0;j < 6;j++)
                    tmp1.col(j) = dGamma2_dq1(j).transpose() * (R2.transpose() * tdot);
                Km.bottomLeftCorner(6, 6) += kt * GTRT2 * dtdot_dq1 + kt * tmp1;
                
                Matrix6 tmp2;
                tmp2.setZero();
                for (int j = 0;j < 6;j++)
                    tmp2.col(j) = dGamma2_dq2(j).transpose() * (R2.transpose() * tdot);
                Km.bottomRightCorner(6, 6) += kt * GTRT2 * dtdot_dq2 + kt * tmp2;
                Km.block(6, 6, 6, 3) += kt * Gamma2.transpose() * math::skew(R2.transpose() * tdot);
                
                Dm.bottomLeftCorner(6, 6) += kt * GTRT2 * dtdot_dphi1;
                Dm.bottomRightCorner(6, 6) += kt * GTRT2 * dtdot_dphi2;

                // design derivatives
                // design params 1
                for (auto ancestor = _contact_body->_joint; ancestor != nullptr; ancestor = ancestor->_parent) {
                    if (ancestor->_design_params_1._active) {
                        for (int k = 0;k < ancestor->_design_params_1._ndof;k++) {
                            int idx = ancestor->_design_params_1._param_index(k);
                            dfm_dp.block(_BVH_body->_index[0], idx, 6, 1) += 
                                kt * (dG2_dp1(idx).transpose() * R2.transpose() * tdot + GTRT2 * dtdot_dp1.col(idx));
                        }
                    }
                }
                // design params 2
                if (_contact_body->_design_params_2._active) {
                    for (int k = 0;k < _contact_body->_design_params_2._ndof;k++) {
                        int idx = _contact_body->_design_params_2._param_index(k) + _sim->_ndof_p1;
                        dfm_dp.block(_BVH_body->_index[0], idx, 6, 1) += 
                            kt * (dG2_dp2(k).transpose() * R2.transpose() * tdot + GTRT2 * dtdot_dp2.col(k));
                    }
                }
                // design params 3
                if (_contact_body->_design_params_3._active) {
                    for (int k = 0;k < 3;k++) {
                        int idx = _contact_body->_design_params_3._param_index(contact_id * 3 + k) + _sim->_ndof_p1 + _sim->_ndof_p2;
                        dfm_dp.block(_BVH_body->_index[0], idx, 6, 1) += 
                            kt * (dG2_dp3(k).transpose() * R2.transpose() * tdot + GTRT2 * dtdot_dp3.col(k));
                    }
                }
                // design params 6
                if (_contact_body->_design_params_6._active) {
                    int idx = _contact_body->_design_params_6._param_index(0) + p6_offset;
                    dfm_dp.block(_BVH_body->_index[0], idx, 6, 1) += fs2 / _contact_body->_contact_scale;
                }
            } else {                                // dynamic frictional force
                RowVector6 dfc_norm_dq1 = (1. / fc_norm) * (fc1.transpose() * dfc1_dq1);
                RowVector6 dfc_norm_dq2 = (1. / fc_norm) * (fc1.transpose() * dfc1_dq2);
                RowVector6 dfc_norm_dphi1 = (1. / fc_norm) * (fc1.transpose() * dfc1_dphi1);
                RowVector6 dfc_norm_dphi2 = (1. / fc_norm) * (fc1.transpose() * dfc1_dphi2);
                RowVector6 dtdot_norm_dq1 = (1. / tdot_norm) * (tdot.transpose() * dtdot_dq1);
                RowVector6 dtdot_norm_dq2 = (1. / tdot_norm) * (tdot.transpose() * dtdot_dq2);
                RowVector6 dtdot_norm_dphi1 = (1. / tdot_norm) * (tdot.transpose() * dtdot_dphi1);
                RowVector6 dtdot_norm_dphi2 = (1. / tdot_norm) * (tdot.transpose() * dtdot_dphi2);
                Vector3 a = tdot / tdot_norm;
                Matrix3 A = (Matrix3::Identity() - a * a.transpose()) / tdot_norm;
                RowVectorX dfc_norm_dp1 = RowVectorX::Zero(_sim->_ndof_p1), dfc_norm_dp2 = RowVectorX::Zero(12), dfc_norm_dp3 = RowVectorX::Zero(3);
                MatrixX da_dp1 = MatrixX::Zero(3, _sim->_ndof_p1), da_dp2 = MatrixX::Zero(3, 12), da_dp3 = MatrixX::Zero(3, 3);
                for (auto ancestor = _contact_body->_joint; ancestor != nullptr; ancestor = ancestor->_parent) {
                    if (ancestor->_design_params_1._active) {
                        int idx_start = ancestor->_design_params_1._param_index(0);
                        int ndof = ancestor->_design_params_1._ndof;
                        dfc_norm_dp1.segment(idx_start, ndof) = fc1.transpose() / fc_norm * dfc1_dp1.middleCols(idx_start, ndof);
                        da_dp1.middleCols(idx_start, ndof) = A * dtdot_dp1.middleCols(idx_start, ndof);
                    }
                }
                if (_contact_body->_design_params_2._active) {
                    dfc_norm_dp2 = fc1.transpose() / fc_norm * dfc1_dp2;
                    da_dp2 = A * dtdot_dp2;
                }
                if (_contact_body->_design_params_3._active) {
                    dfc_norm_dp3 = fc1.transpose() / fc_norm * dfc1_dp3;
                    da_dp3 = A * dtdot_dp3;
                }

                // general contact body
                Vector6 fd1 = -_mu * fc_norm / tdot_norm * GTRT1 * tdot;
                fm.head(6) += fd1;

                // derivatives
                Km.topLeftCorner(6, 6) -=
                    _mu / tdot_norm * GTRT1 * (tdot * dfc_norm_dq1 - fc_norm / tdot_norm * tdot * dtdot_norm_dq1 + fc_norm * dtdot_dq1);
                Km.block(0, 0, 6, 3) -=
                    _mu * Gamma1.transpose() * math::skew(fc_norm / tdot_norm * R1.transpose() * tdot);
                Km.topRightCorner(6, 6) -=
                    _mu / tdot_norm * GTRT1 * (tdot * dfc_norm_dq2 - fc_norm / tdot_norm * tdot * dtdot_norm_dq2 + fc_norm * dtdot_dq2);

                Dm.topLeftCorner(6, 6) -=
                    _mu / tdot_norm * GTRT1 * (tdot * dfc_norm_dphi1 - fc_norm / tdot_norm * tdot * dtdot_norm_dphi1 + fc_norm * dtdot_dphi1);
                Dm.topRightCorner(6, 6) -=
                    _mu / tdot_norm * GTRT1 * (tdot * dfc_norm_dphi2 - fc_norm / tdot_norm * tdot * dtdot_norm_dphi2 + fc_norm * dtdot_dphi2);

                // design derivatives
                // design params 1
                for (auto ancestor = _contact_body->_joint; ancestor != nullptr; ancestor = ancestor->_parent) {
                    if (ancestor->_design_params_1._active) {
                        for (int k = 0;k < ancestor->_design_params_1._ndof;k++) {
                            int idx = ancestor->_design_params_1._param_index(k);
                            Matrix3 dR_dp1 = _contact_body->_dE0i_dp1(idx).topLeftCorner(3, 3);
                            dfm_dp.block(_contact_body->_index[0], idx, 6, 1) += 
                                -_mu * (Gamma1.transpose() * dR_dp1.transpose() * fc_norm * a
                                    + GTRT1 * (dfc_norm_dp1(idx) * a + fc_norm * da_dp1.col(idx)));
                        }
                    }
                }
                // design params 2
                if (_contact_body->_design_params_2._active) {
                    for (int k = 0;k < _contact_body->_design_params_2._ndof;k++) {
                        int idx = _contact_body->_design_params_2._param_index(k) + _sim->_ndof_p1;
                        Matrix3 dR_dp2 = _contact_body->_dE0i_dp2(k).topLeftCorner(3, 3);
                        dfm_dp.block(_contact_body->_index[0], idx, 6, 1) += 
                            -_mu * (Gamma1.transpose() * dR_dp2.transpose() * fc_norm * a
                                + GTRT1 * (dfc_norm_dp2(k) * a + fc_norm * da_dp2.col(k)));
                    }
                }
                // design params 3
                if (_contact_body->_design_params_3._active) {
                    for (int k = 0;k < 3;k++) {
                        int idx = _contact_body->_design_params_3._param_index(contact_id * 3 + k) + _sim->_ndof_p1 + _sim->_ndof_p2;
                        dfm_dp.block(_contact_body->_index[0], idx, 6, 1) += 
                            -_mu * (dG1_dp3(k).transpose() * R1.transpose() * fc_norm * a)
                            -_mu * GTRT1 * (dfc_norm_dp3(k) * a + fc_norm * da_dp3.col(k));
                    }
                }
                // design params 6
                if (_contact_body->_design_params_6._active) {
                    int idx = _contact_body->_design_params_6._param_index(0) + p6_offset;
                    dfm_dp.block(_contact_body->_index[0], idx, 6, 1) += fd1 / _contact_body->_contact_scale;
                }

                // BVH contact body
                Vector6 fd2 = _mu * fc_norm / tdot_norm * GTRT2 * tdot;
                fm.tail(6) += fd2;

                // derivatives
                Matrix6 tmp1;
                tmp1.setZero();
                for (int j = 0;j < 6;j++)
                    tmp1.col(j) = dGamma2_dq1(j).transpose() * (R2.transpose() * tdot);
                Km.bottomLeftCorner(6, 6) += 
                    _mu / tdot_norm * GTRT2 * (tdot * dfc_norm_dq1 - tdot * fc_norm / tdot_norm * dtdot_norm_dq1 + fc_norm * dtdot_dq1) +
                    _mu * fc_norm / tdot_norm * tmp1;
                
                Matrix6 tmp2;
                tmp2.setZero();
                for (int j = 0;j < 6;j++)
                    tmp2.col(j) = dGamma2_dq2(j).transpose() * (R2.transpose() * tdot);
                Km.bottomRightCorner(6, 6) += 
                    _mu / tdot_norm * GTRT2 * (tdot * dfc_norm_dq2 - tdot * fc_norm / tdot_norm * dtdot_norm_dq2 + fc_norm * dtdot_dq2) +
                    _mu * fc_norm / tdot_norm * tmp2;
                Km.block(6, 6, 6, 3) += 
                    _mu * Gamma2.transpose() * math::skew(fc_norm / tdot_norm * R2.transpose() * tdot);

                Dm.bottomLeftCorner(6, 6) += 
                    _mu / tdot_norm * GTRT2 * (tdot * dfc_norm_dphi1 - tdot * fc_norm / tdot_norm * dtdot_norm_dphi1 + fc_norm * dtdot_dphi1);
                Dm.bottomRightCorner(6, 6) += 
                    _mu / tdot_norm * GTRT2 * (tdot * dfc_norm_dphi2 - tdot * fc_norm / tdot_norm * dtdot_norm_dphi2 + fc_norm * dtdot_dphi2);
                
                // design derivatives
                // design params 1
                for (auto ancestor = _contact_body->_joint; ancestor != nullptr; ancestor = ancestor->_parent) {
                    if (ancestor->_design_params_1._active) {
                        for (int k = 0;k < ancestor->_design_params_1._ndof;k++) {
                            int idx = ancestor->_design_params_1._param_index(k);
                            dfm_dp.block(_BVH_body->_index[0], idx, 6, 1) += 
                                _mu * (dG2_dp1(idx).transpose() * R2.transpose() * fc_norm * a
                                    + GTRT2 * (dfc_norm_dp1(idx) * a + fc_norm * da_dp1.col(idx)));
                        }
                    }
                }
                // design params 2
                if (_contact_body->_design_params_2._active) {
                    for (int k = 0;k < _contact_body->_design_params_2._ndof;k++) {
                        int idx = _contact_body->_design_params_2._param_index(k) + _sim->_ndof_p1;
                        dfm_dp.block(_BVH_body->_index[0], idx, 6, 1) += 
                                _mu * (dG2_dp2(k).transpose() * R2.transpose() * fc_norm * a
                                    + GTRT2 * (dfc_norm_dp2(k) * a + fc_norm * da_dp2.col(k)));
                    }
                }
                // design params 3
                if (_contact_body->_design_params_3._active) {
                    for (int k = 0;k < 3;k++) {
                        int idx = _contact_body->_design_params_3._param_index(contact_id * 3 + k) + _sim->_ndof_p1 + _sim->_ndof_p2;
                        dfm_dp.block(_BVH_body->_index[0], idx, 6, 1) += 
                                _mu * (dG2_dp3(k).transpose() * R2.transpose() * fc_norm * a
                                    + GTRT2 * (dfc_norm_dp3(k) * a + fc_norm * da_dp3.col(k)));
                    }
                }
                // design params 6
                if (_contact_body->_design_params_6._active) {
                    int idx = _contact_body->_design_params_6._param_index(0) + p6_offset;
                    dfm_dp.block(_BVH_body->_index[0], idx, 6, 1) += fd2 / _contact_body->_contact_scale;
                }
            }
        }
    }
}

void ForceGeneralBVHContact::test_derivatives() {
    std::cerr << "**************************** General-BVH Collision Derivatives ***************************" << std::endl;
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
    this->_BVH_body->_E_0i = E2;
    this->_BVH_body->_phi = phi2;

    std::vector<Contact> contacts;
    // contacts.push_back(Contact(xi1, xw1, 0.0, Vector3::Zero()));
    collision_detection_general_BVH(this->_contact_body, this->_BVH_body, contacts);
    std::cerr << "# Contacts: " << contacts.size() << std::endl;

    VectorX fm;
    MatrixX Km, Dm;
    computeForceWithDerivative(contacts, fm, Km, Dm, true);

    std::vector<Contact> contacts_copy = contacts;

    // Km
    MatrixX Km_fd = MatrixX::Zero(12, 12);
    for (int i = 0;i < 6;i++) {
        Vector6 dq = Vector6::Zero();
        dq[i] = eps;
        Matrix4 E1_pos = E1 * math::exp(dq);
        this->_contact_body->_E_0i = E1_pos;

        VectorX d;
        Matrix3X dd_dx;
        Matrix3X xw1_pos;
        xw1_pos.resize(3, contacts.size());
        for (int j = 0; j < contacts.size(); ++j) {
            xw1_pos.col(j) = E1_pos.topLeftCorner(3, 3) * contacts[j]._xi + E1_pos.topRightCorner(3, 1);
        }
        this->_BVH_body->distance_normal_parallel(xw1_pos, d, dd_dx);
        for (int j = 0; j < contacts.size(); ++j) {
            contacts[j]._xw = xw1_pos.col(j);
            contacts[j]._d = d[j];
            contacts[j]._normal = dd_dx.col(j);
        }
        
        VectorX fm_pos;
        computeForce(contacts, fm_pos);
        Km_fd.col(i) = (fm_pos - fm) / eps;
    }
    this->_contact_body->_E_0i = E1;

    contacts = contacts_copy;
    
    print_error("Km(q1) ", Km.leftCols(6), Km_fd.leftCols(6));

    for (int i = 0;i < 6;i++) {
        Vector6 dq = Vector6::Zero();
        dq[i] = eps;
        Matrix4 E2_pos = E2 * math::exp(dq);
        this->_BVH_body->_E_0i = E2_pos;

        VectorX d;
        Matrix3X dd_dx;
        Matrix3X xw1_;
        xw1_.resize(3, contacts.size());
        for (int j = 0; j < contacts.size(); ++j) {
            xw1_.col(j) = contacts[j]._xw;
        }
        this->_BVH_body->distance_normal_parallel(xw1_, d, dd_dx);
        for (int j = 0; j < contacts.size(); ++j) {
            contacts[j]._d = d[j];
            contacts[j]._normal = dd_dx.col(j);
        }

        VectorX fm_pos;
        computeForce(contacts, fm_pos);
        Km_fd.col(6 + i) = (fm_pos - fm) / eps;
    }
    this->_BVH_body->_E_0i = E2;
    print_error("Km(q2) ", Km.rightCols(6), Km_fd.rightCols(6));

    contacts = contacts_copy;

    // Dm
    MatrixX Dm_fd = MatrixX::Zero(12, 12);
    for (int i = 0;i < 6;i++) {
        Vector6 phi1_pos = phi1;
        phi1_pos[i] += eps;
        this->_contact_body->_phi = phi1_pos;
        VectorX fm_pos;
        computeForce(contacts, fm_pos);
        Dm_fd.col(i) = (fm_pos - fm) / eps;
    }
    this->_contact_body->_phi = phi1;
    print_error("Dm(phi1) ", Dm.leftCols(6), Dm_fd.leftCols(6));

    for (int i = 0;i < 6;i++) {
        Vector6 phi2_pos = phi2;
        phi2_pos[i] += eps;
        this->_BVH_body->_phi = phi2_pos;
        VectorX fm_pos;
        computeForce(contacts, fm_pos);
        Dm_fd.col(6 + i) = (fm_pos - fm) / eps;
    }
    this->_BVH_body->_phi = phi2;
    print_error("Dm(phi2) ", Dm.rightCols(6), Dm_fd.rightCols(6));
}

void ForceGeneralBVHContact::test_derivatives_runtime() {
    // detect the contact points between two bodies
    std::vector<Contact> contacts;
    contacts.clear();
    collision_detection_general_BVH(_contact_body, _BVH_body, contacts);
    std::vector<Contact> contacts_copy = contacts;
    
    // test derivative for BVH body
    for (int i = 0;i < contacts.size();i++) {
        Vector3 xi1 = contacts[i]._xi;
        Vector3 xw1 = _contact_body->_E_0i.topLeftCorner(3, 3) * xi1 + _contact_body->_E_0i.topRightCorner(3, 1);
        Vector3 xw1_dot = _contact_body->_E_0i.topLeftCorner(3, 3) * (math::skew(xi1).transpose() * _contact_body->_phi.head(3) + _contact_body->_phi.tail(3));
        _BVH_body->test_collision_derivatives_runtime(xw1, xw1_dot);
    }

    VectorX fm;
    MatrixX Km, Dm;
    computeForceWithDerivative(contacts, fm, Km, Dm);

    dtype eps = 1e-7;
    // generate random xi1, E1, E2, phi1, phi2
    MatrixX E1 = this->_contact_body->_E_0i;
    VectorX phi1 = this->_contact_body->_phi;
    MatrixX E2 = this->_BVH_body->_E_0i;
    VectorX phi2 = this->_BVH_body->_phi;

    // Km
    MatrixX Km_fd = MatrixX::Zero(12, 12);
    for (int i = 0;i < 6;i++) {
        Vector6 dq = Vector6::Zero();
        dq[i] = eps;
        Matrix4 E1_pos = E1 * math::exp(dq);
        this->_contact_body->_E_0i = E1_pos;

        VectorX d;
        Matrix3X dd_dx;
        Matrix3X xw1_pos;
        xw1_pos.resize(3, contacts.size());
        for (int j = 0; j < contacts.size(); ++j) {
            xw1_pos.col(j) = E1_pos.topLeftCorner(3, 3) * contacts[j]._xi + E1_pos.topRightCorner(3, 1);
        }
        this->_BVH_body->distance_normal_parallel(xw1_pos, d, dd_dx);
        for (int j = 0; j < contacts.size(); ++j) {
            contacts[j]._xw = xw1_pos.col(j);
            contacts[j]._d = d[j];
            contacts[j]._normal = dd_dx.col(j);
        }

        VectorX fm_pos;
        computeForce(contacts, fm_pos);
        Km_fd.col(i) = (fm_pos - fm) / eps;
    }
    this->_contact_body->_E_0i = E1;
    print_error("General-BVH Collision Derivatives: Km(q1) ", Km.leftCols(6), Km_fd.leftCols(6));

    contacts = contacts_copy;

    for (int i = 0;i < 6;i++) {
        Vector6 dq = Vector6::Zero();
        dq[i] = eps;
        Matrix4 E2_pos = E2 * math::exp(dq);
        this->_BVH_body->_E_0i = E2_pos;

        VectorX d;
        Matrix3X dd_dx;
        Matrix3X xw1_;
        xw1_.resize(3, contacts.size());
        for (int j = 0; j < contacts.size(); ++j) {
            xw1_.col(j) = contacts[j]._xw;
        }
        this->_BVH_body->distance_normal_parallel(xw1_, d, dd_dx);
        for (int j = 0; j < contacts.size(); ++j) {
            contacts[j]._d = d[j];
            contacts[j]._normal = dd_dx.col(j);
        }
        
        VectorX fm_pos;
        computeForce(contacts, fm_pos);
        Km_fd.col(6 + i) = (fm_pos - fm) / eps;
    }
    this->_BVH_body->_E_0i = E2;
    print_error("General-BVH Collision Derivatives: Km(q2) ", Km.rightCols(6), Km_fd.rightCols(6));

    contacts = contacts_copy;

    // Dm
    MatrixX Dm_fd = MatrixX::Zero(12, 12);
    for (int i = 0;i < 6;i++) {
        Vector6 phi1_pos = phi1;
        phi1_pos[i] += eps;
        this->_contact_body->_phi = phi1_pos;
        VectorX fm_pos;
        computeForce(contacts, fm_pos);
        Dm_fd.col(i) = (fm_pos - fm) / eps;
    }
    this->_contact_body->_phi = phi1;
    print_error("General-BVH Collision Derivatives: Dm(phi1) ", Dm.leftCols(6), Dm_fd.leftCols(6));

    for (int i = 0;i < 6;i++) {
        Vector6 phi2_pos = phi2;
        phi2_pos[i] += eps;
        this->_BVH_body->_phi = phi2_pos;
        VectorX fm_pos;
        computeForce(contacts, fm_pos);
        Dm_fd.col(6 + i) = (fm_pos - fm) / eps;
    }
    this->_BVH_body->_phi = phi2;
    print_error("General-BVH Collision Derivatives: Dm(phi2) ", Dm.rightCols(6), Dm_fd.rightCols(6));

    test_design_derivatives_runtime();
}

void ForceGeneralBVHContact::test_design_derivatives_runtime() {
    // detect the contact points between two bodies
    std::vector<Contact> contacts;
    contacts.clear();

    collision_detection_general_BVH(_contact_body, _BVH_body, contacts);

    if (contacts.size() > 0) {
        VectorX fm = VectorX::Zero(_sim->_ndof_m), fm_sub;
        MatrixX Km_sub, Dm_sub;
        MatrixX dfm_dp = MatrixX::Zero(_sim->_ndof_m, _sim->_ndof_p);
        computeForceWithDerivative(contacts, fm_sub, Km_sub, Dm_sub, dfm_dp, false);
        fm.segment(_contact_body->_index[0], 6) += fm_sub.head(6);
        fm.segment(_BVH_body->_index[0], 6) += fm_sub.tail(6);

        VectorX design_params = _sim->get_design_params();

        dtype eps = 1e-7;
        for (int ii = 0;ii < 1;ii++) {
            MatrixX dfm_dp_fd = MatrixX::Zero(_sim->_ndof_m, _sim->_ndof_p);

            for (int i = 0;i < _sim->_ndof_p;i++) {
                VectorX design_params_pos = design_params;
                design_params_pos(i) += eps;
                std::vector<Contact> contacts_pos = contacts;

                _sim->set_design_params(design_params_pos);
                _sim->update_robot(false);
                std::vector<Vector3> all_contact_points = _contact_body->get_contact_points();
                for (int j = 0;j < contacts.size();j++)
                    contacts_pos[j]._xi = all_contact_points[contacts_pos[j]._id];

                VectorX fm_pos = VectorX::Zero(_sim->_ndof_m), fm_sub_pos;
                MatrixX Km_pos, Dm_pos;
                MatrixX dfm_dp_pos = MatrixX::Zero(_sim->_ndof_m, _sim->_ndof_p);
                computeForceWithDerivative(contacts_pos, fm_sub_pos, Km_pos, Dm_pos, dfm_dp_pos, false);
                fm_pos.segment(_contact_body->_index[0], 6) += fm_sub_pos.head(6);
                fm_pos.segment(_BVH_body->_index[0], 6) += fm_sub_pos.tail(6);
                dfm_dp_fd.col(i) = (fm_pos - fm) / eps;
            }
            // std::cerr << "eps = " << eps << std::endl;
            print_error("General-BVH Collision Derivatives: dfm_dp1", dfm_dp.leftCols(_sim->_ndof_p1), dfm_dp_fd.leftCols(_sim->_ndof_p1));
            print_error("General-BVH Collision Derivatives: dfm_dp2", dfm_dp.middleCols(_sim->_ndof_p1, _sim->_ndof_p2), dfm_dp_fd.middleCols(_sim->_ndof_p1, _sim->_ndof_p2));
            print_error("General-BVH Collision Derivatives: dfm_dp3", dfm_dp.middleCols(_sim->_ndof_p1 + _sim->_ndof_p2, _sim->_ndof_p3), dfm_dp_fd.middleCols(_sim->_ndof_p1 + _sim->_ndof_p2, _sim->_ndof_p3));
            print_error("General-BVH Collision Derivatives: dfm_dp6", dfm_dp.middleCols(_sim->_ndof_p1 + _sim->_ndof_p2 + _sim->_ndof_p3 + _sim->_ndof_p4 + _sim->_ndof_p5, _sim->_ndof_p6), dfm_dp_fd.middleCols(_sim->_ndof_p1 + _sim->_ndof_p2 + _sim->_ndof_p3 + _sim->_ndof_p4 + _sim->_ndof_p5, _sim->_ndof_p6));
            eps /= 10.;
        }

        _sim->set_design_params(design_params);
        _sim->update_robot(true);
    }
}

}