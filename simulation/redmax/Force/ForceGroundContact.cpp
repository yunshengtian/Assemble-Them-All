#include "ForceGroundContact.h"
#include "Body/Body.h"
#include "Simulation.h"
#include "Joint/Joint.h"
#include "CollisionDetection/CollisionDetection.h"

namespace redmax {

ForceGroundContact::ForceGroundContact(
    Simulation* sim, 
    Body* contact_body, Matrix4 E_g, 
    dtype kn, dtype kt, dtype mu, dtype damping) : Force(sim) {
        
    _contact_body = contact_body;
    _E_g = E_g;
    _kn = kn;
    _kt = kt;
    _mu = mu;
    _damping = damping;
}

void ForceGroundContact::set_transform(Matrix4 E_g) {
    _E_g = E_g;
}

void ForceGroundContact::set_stiffness(dtype kn, dtype kt) {
    _kn = kn;
    _kt = kt;
}

void ForceGroundContact::set_friction(dtype mu) {
    _mu = mu;
}

void ForceGroundContact::set_damping(dtype damping) {
    _damping = damping;
}

void ForceGroundContact::computeForce(VectorX& fm, VectorX& fr, bool verbose) {
    auto t0 = clock();

    // detect the contact points between two bodies
    std::vector<Contact> contacts;
    contacts.clear();

    collision_detection_ground_object(_E_g, _contact_body, contacts);

    auto t1 = clock();

    computeForce(contacts, fm, verbose);

    auto t2 = clock();

    _time.add("collision detection", t1 - t0);
    _time.add("compute force", t2 - t1);
}

void ForceGroundContact::computeForceWithDerivative(
    VectorX& fm, VectorX& fr, 
    MatrixX& Km, MatrixX& Dm, 
    MatrixX& Kr, MatrixX& Dr,
    bool verbose) {
    
    auto t0 = clock();

    // detect the contact points between two bodies
    std::vector<Contact> contacts;
    contacts.clear();
    collision_detection_ground_object(_E_g, _contact_body, contacts);

    auto t1 = clock();

    computeForceWithDerivative(contacts, fm, Km, Dm, verbose);

    auto t2 = clock();

    _time.add("collision detection", t1 - t0);
    _time.add("compute force", t2 - t1);
}

void ForceGroundContact::computeForceWithDerivative(
    VectorX& fm, VectorX& fr, 
    MatrixX& Km, MatrixX& Dm, 
    MatrixX& Kr, MatrixX& Dr,
    MatrixX& dfm_dp, MatrixX& dfr_dp,
    bool verbose) {
    
    auto t0 = clock();

    // detect the contact points between two bodies
    std::vector<Contact> contacts;
    contacts.clear();
    collision_detection_ground_object(_E_g, _contact_body, contacts);

    auto t1 = clock();

    computeForceWithDerivative(contacts, fm, Km, Dm, dfm_dp, verbose);

    auto t2 = clock();

    _time.add("collision detection", t1 - t0);
    _time.add("compute force", t2 - t1);
}

void ForceGroundContact::computeForce(std::vector<Contact> &contacts, VectorX& fm, bool verbose) {
    Vector3 xg = _E_g.topRightCorner(3, 1);
    Vector3 ng = _E_g.block(0, 2, 3, 1);
    Matrix3 N = ng * ng.transpose();
    Matrix3 T = Matrix3::Identity() - N;
    Matrix3 R = (_contact_body->_E_0i).topLeftCorner(3, 3);
    Vector3 p = (_contact_body->_E_0i).topRightCorner(3, 1);
    Vector6 phi = _contact_body->_phi;

    for (int i = 0;i < contacts.size();i++) {
        Vector3 xi = contacts[i]._xi;
        Vector3 xw = R * xi + p;
        dtype d = ng.transpose() * (xw - xg);
        MatrixX xibrac = math::skew(xi);

        // Contact force
        MatrixX G = math::gamma(xi);
        MatrixX J = R * G;
        Vector3 vwi = J * phi;
        Vector3 fc = -_kn * ng * d - _damping * N * vwi;

        fm.segment(_contact_body->_index[0], 6).noalias() += J.transpose() * fc;
        
        // Friction force
        if (_mu < constants::eps)
            continue;
        
        Vector3 a = T * vwi;
        dtype anorm = a.norm();
        if (_mu * fabs(_kn * d) >= _kt * anorm - constants::eps) {
            // Static friction force
            Vector3 fs = -_kt * a;
            // std::cerr << "fs = " << fs.transpose() << std::endl;
            fm.segment(_contact_body->_index[0], 6).noalias() += J.transpose() * fs;
        } else {
            dtype mukn = _mu * _kn;
            Vector3 t = a / anorm;
            Vector3 fd = mukn * d * t;
            // std::cerr << "fd = " << fd.transpose() << std::endl;
            fm.segment(_contact_body->_index[0], 6).noalias() += J.transpose() * fd;
        }
    }   
}

void ForceGroundContact::computeForceWithDerivative(std::vector<Contact> &contacts, 
    VectorX& fm, MatrixX& Km, MatrixX& Dm, bool verbose) {
        
    Vector3 xg = _E_g.topRightCorner(3, 1);
    Vector3 ng = _E_g.block(0, 2, 3, 1);
    Matrix3 N = ng * ng.transpose();
    Matrix3 I = Matrix3::Identity();
    Matrix3 Z = Matrix3::Zero();
    Matrix3 T = I - N;
    Matrix3 R = (_contact_body->_E_0i).topLeftCorner(3, 3);
    Vector3 p = (_contact_body->_E_0i).topRightCorner(3, 1);
    Vector6 phi = _contact_body->_phi;
    Matrix3 e1b = math::skew(Vector3::UnitX());
    Matrix3 e2b = math::skew(Vector3::UnitY());
    Matrix3 e3b = math::skew(Vector3::UnitZ());
    Matrix3 RNR = R.transpose() * N * R;
    Matrix3 pxgtmp = math::skew(R.transpose() * N * (p - xg));
    
    for (int i = 0;i < contacts.size();i++) {
        Vector3 xi = contacts[i]._xi;
        Vector3 xw = R * xi + p;
        dtype d = ng.transpose() * (xw - xg);
        Vector3 RNRxi = RNR * xi;
        MatrixX xibrac = math::skew(xi);
        MatrixX G = math::gamma(xi);
        Vector3 Gphi = G * phi;
        Vector3 vwi = R * Gphi;
        Vector3 RNRGphi = RNR * Gphi;
        Matrix3 Gphibrac = math::skew(Gphi);

        // Contact force
        
        MatrixX J = R * G;

        Vector3 fc = -_kn * ng * d - _damping * N * vwi;
        fm.segment(_contact_body->_index[0], 6).noalias() += J.transpose() * fc;

        Matrix3 tmp1 = Matrix3::Zero();
        tmp1 << -e1b * RNRxi, -e2b * RNRxi, -e3b * RNRxi;
        tmp1.noalias() += - RNR * xibrac + pxgtmp;

        Matrix3 tmp2 = Matrix3::Zero();
        tmp2 << -e1b * RNRGphi, -e2b * RNRGphi, -e3b * RNRGphi;
        tmp2.noalias() += -RNR * Gphibrac;
        Km.block(_contact_body->_index[0], _contact_body->_index[0], 6, 6) -= 
            _kn * G.transpose() * (Matrix<dtype, 3, 6>() << tmp1, RNR).finished()
            + _damping * G.transpose() * (Matrix<dtype, 3, 6>() << tmp2, Z).finished();
        Dm.block(_contact_body->_index[0], _contact_body->_index[0], 6, 6) -=
            _damping * G.transpose() * RNR * G;

        // Friction force
        if (_mu < constants::eps)
            continue;
        
        Vector3 xwdot = J * phi;
        Vector3 a = T * xwdot;
        dtype anorm = a.norm();
        if (_mu * fabs(_kn * d) >= _kt * anorm - constants::eps) {
            // Static friction force
            Vector3 fs = -_kt * a;
            fm.segment(_contact_body->_index[0], 6).noalias() += J.transpose() * fs;
            
            Matrix6 D = -_kt * J.transpose() * T * J;
            Matrix3 B = R.transpose() * T * R;
            MatrixX tmp3(3, 6);
            tmp3 << (B * e1b - e1b * B) * Gphi, (B * e2b - e2b * B) * Gphi, (B * e3b - e3b * B) * Gphi, Z;
            Matrix6 K = -_kt * G.transpose() * tmp3;

            Dm.block(_contact_body->_index[0], _contact_body->_index[0], 6, 6).noalias() += D;
            Km.block(_contact_body->_index[0], _contact_body->_index[0], 6, 6).noalias() += K;
        } else {
            // Dynamic friction force
            dtype mukn = _mu * _kn;
            Vector3 t = a / anorm;
            Vector3 fd = mukn * d * t;
            fm.segment(_contact_body->_index[0], 6).noalias() += J.transpose() * fd;

            Matrix3 A = (I - t * t.transpose()) / anorm;
            Matrix6 D = mukn * J.transpose() * d * A * T * J;
            Vector3 Rt = R.transpose() * t;
            MatrixX K1(3, 6);
            K1 << e1b * Rt, e2b * Rt, e3b * Rt, Z;
            K1 *= -d;
            MatrixX K2 = R.transpose() * t * ng.transpose() * J;
            MatrixX tmp3(3, 6);
            tmp3 << math::skew(G * phi), Z;
            MatrixX K3 = -d * R.transpose() * A * T * R * tmp3;
            Matrix6 K = mukn * G.transpose() * (K1 + K2 + K3);

            Dm.block(_contact_body->_index[0], _contact_body->_index[0], 6, 6).noalias() += D;
            Km.block(_contact_body->_index[0], _contact_body->_index[0], 6, 6).noalias() += K;
        }
    }
}

void ForceGroundContact::computeForceWithDerivative(
    std::vector<Contact> &contacts, 
    VectorX& fm, MatrixX& Km, MatrixX& Dm, MatrixX& dfm_dp, 
    bool verbose) {

    Vector3 xg = _E_g.topRightCorner(3, 1);
    Vector3 ng = _E_g.block(0, 2, 3, 1);
    Matrix3 N = ng * ng.transpose();
    Matrix3 I = Matrix3::Identity();
    Matrix3 Z = Matrix3::Zero();
    Matrix3 T = I - N;
    Matrix3 R = (_contact_body->_E_0i).topLeftCorner(3, 3);
    Vector3 p = (_contact_body->_E_0i).topRightCorner(3, 1);
    Vector6 phi = _contact_body->_phi;
    Matrix3 e1b = math::skew(Vector3::UnitX());
    Matrix3 e2b = math::skew(Vector3::UnitY());
    Matrix3 e3b = math::skew(Vector3::UnitZ());
    Matrix3 RNR = R.transpose() * N * R;
    Matrix3 pxgtmp = math::skew(R.transpose() * N * (p - xg));
    
    if (verbose)
        std::cerr << "num ground contacts = " << contacts.size() << std::endl;

    for (int i = 0;i < contacts.size();i++) {
        Vector3 xi = contacts[i]._xi;
        Vector3 xw = R * xi + p;
        dtype d = ng.transpose() * (xw - xg);
        Vector3 RNRxi = RNR * xi;
        MatrixX xibrac = math::skew(xi);
        MatrixX G = math::gamma(xi);
        Vector3 Gphi = G * phi;
        Matrix36 RG = R * G;
        Vector3 vwi = R * Gphi;
        Vector3 RNRGphi = RNR * Gphi;
        Matrix3 Gphibrac = math::skew(Gphi);
        Vector3 Nvw = N * vwi;
        int contact_id = contacts[i]._id;

        // Contact force
        
        MatrixX J = R * G;

        Vector3 fc = -_kn * ng * d - _damping * Nvw;
        fm.segment(_contact_body->_index[0], 6).noalias() += J.transpose() * fc;

        // derivatives
        Matrix3 tmp1 = Matrix3::Zero();
        tmp1 << -e1b * RNRxi, -e2b * RNRxi, -e3b * RNRxi;
        tmp1.noalias() += - RNR * xibrac + pxgtmp;

        Matrix3 tmp2 = Matrix3::Zero();
        tmp2 << -e1b * RNRGphi, -e2b * RNRGphi, -e3b * RNRGphi;
        tmp2.noalias() += -RNR * Gphibrac;
        Km.block(_contact_body->_index[0], _contact_body->_index[0], 6, 6) -= 
            _kn * G.transpose() * (Matrix<dtype, 3, 6>() << tmp1, RNR).finished()
            + _damping * G.transpose() * (Matrix<dtype, 3, 6>() << tmp2, Z).finished();
        Dm.block(_contact_body->_index[0], _contact_body->_index[0], 6, 6) -=
            _damping * G.transpose() * RNR * G;

        // design derivatives
        // tools
        JacobianMatrixVector dG_dp3(3, 6, 3);
        for (int k = 0;k < 3;k++) {
            dG_dp3(k).leftCols(3) = math::skew(Vector3::Unit(k)).transpose();
        }

        RowVectorX dd_dp1 = RowVectorX::Zero(_sim->_ndof_p1), dd_dp2 = RowVectorX::Zero(12), dd_dp3 = RowVector3::Zero(3);
        MatrixX dvw_dp1 = MatrixX::Zero(3, _sim->_ndof_p1), dvw_dp2 = MatrixX::Zero(3, 12), dvw_dp3 = MatrixX::Zero(3, 3);
        // design params 1
        for (auto ancestor = _contact_body->_joint; ancestor != nullptr; ancestor = ancestor->_parent) {
            if (ancestor->_design_params_1._active) {
                for (int k = 0;k < ancestor->_design_params_1._ndof;k++) {
                    int idx = ancestor->_design_params_1._param_index(k);
                    Matrix3 dR_dp1 = _contact_body->_dE0i_dp1(idx).topLeftCorner(3, 3);
                    Vector3 dp_dp1 = _contact_body->_dE0i_dp1(idx).topRightCorner(3, 1);
                    dd_dp1(idx) = ng.transpose() * (dR_dp1 * xi + dp_dp1);
                    dvw_dp1.col(idx) = dR_dp1 * Gphi + RG * _contact_body->_dphi_dp1.col(idx);
                }
            }
        }
        // design params 2
        if (_contact_body->_design_params_2._active) {
            for (int k = 0;k < _contact_body->_design_params_2._ndof;k++) {
                Matrix3 dR_dp2 = _contact_body->_dE0i_dp2(k).topLeftCorner(3, 3);
                Vector3 dp_dp2 = _contact_body->_dE0i_dp2(k).topRightCorner(3, 1);
                dd_dp2(k) = ng.transpose() * (dR_dp2 * xi + dp_dp2);
                dvw_dp2.col(k) = dR_dp2 * Gphi;
            }
            dvw_dp2 += RG * _contact_body->_dphi_dp2;
        }
        // design params 3
        if (_contact_body->_design_params_3._active) {
            dd_dp3 = ng.transpose() * R;
            for (int k = 0;k < 3;k++)
                dvw_dp3.col(k) = R * dG_dp3(k) * phi;
        }

        // dfc_dp
        // design params 1
        for (auto ancestor = _contact_body->_joint; ancestor != nullptr; ancestor = ancestor->_parent) {
            if (ancestor->_design_params_1._active) {
                for (int k = 0;k < ancestor->_design_params_1._ndof;k++) {
                    int idx = ancestor->_design_params_1._param_index(k);
                    Matrix3 dR_dp1 = _contact_body->_dE0i_dp1(idx).topLeftCorner(3, 3);
                    dfm_dp.block(_contact_body->_index[0], idx, 6, 1) += G.transpose() * dR_dp1.transpose() * fc - J.transpose() * (_kn * ng * dd_dp1(idx) + _damping * N * dvw_dp1.col(idx));
                }
            }
        }
        // design params 2
        if (_contact_body->_design_params_2._active) {
            for (int k = 0;k < _contact_body->_design_params_2._ndof;k++) {
                int idx = _contact_body->_design_params_2._param_index(k) + _sim->_ndof_p1;
                Matrix3 dR_dp2 = _contact_body->_dE0i_dp2(k).topLeftCorner(3, 3);
                dfm_dp.block(_contact_body->_index[0], idx, 6, 1) += G.transpose() * dR_dp2.transpose() * fc - J.transpose() * (_kn * ng * dd_dp2(k) + _damping * N * dvw_dp2.col(k));
            }
        }
        // design params 3
        if (_contact_body->_design_params_3._active) {
            for (int k = 0;k < 3;k++) {
                int idx = _contact_body->_design_params_3._param_index(contact_id * 3 + k) + _sim->_ndof_p1 + _sim->_ndof_p2;
                dfm_dp.block(_contact_body->_index[0], idx, 6, 1) += dG_dp3(k).transpose() * R.transpose() * fc;
            }
            int idx = _contact_body->_design_params_3._param_index(contact_id * 3) + _sim->_ndof_p1 + _sim->_ndof_p2;
            dfm_dp.block(_contact_body->_index[0], idx, 6, 3) -= J.transpose() * (_kn * ng * dd_dp3 + _damping * N * dvw_dp3);
        }

        // Friction force
        if (_mu < constants::eps)
            continue;
        
        Vector3 xwdot = J * phi;
        Vector3 a = T * xwdot;
        dtype anorm = a.norm();
        if (_mu * fabs(_kn * d) >= _kt * anorm - constants::eps) {
            // Static friction force
            Vector3 fs = -_kt * a;
            fm.segment(_contact_body->_index[0], 6).noalias() += J.transpose() * fs;
            
            Matrix6 D = -_kt * J.transpose() * T * J;
            Matrix3 B = R.transpose() * T * R;
            MatrixX tmp3(3, 6);
            tmp3 << (B * e1b - e1b * B) * Gphi, (B * e2b - e2b * B) * Gphi, (B * e3b - e3b * B) * Gphi, Z;
            Matrix6 K = -_kt * G.transpose() * tmp3;

            Dm.block(_contact_body->_index[0], _contact_body->_index[0], 6, 6).noalias() += D;
            Km.block(_contact_body->_index[0], _contact_body->_index[0], 6, 6).noalias() += K;

            // design derivatives
            // design params 1
            MatrixX JTT = J.transpose() * T;
            for (auto ancestor = _contact_body->_joint; ancestor != nullptr; ancestor = ancestor->_parent) {
                if (ancestor->_design_params_1._active) {
                    for (int k = 0;k < ancestor->_design_params_1._ndof;k++) {
                        int idx = ancestor->_design_params_1._param_index(k);
                        Matrix3 dR_dp1 = _contact_body->_dE0i_dp1(idx).topLeftCorner(3, 3);
                        dfm_dp.block(_contact_body->_index[0], idx, 6, 1) -= _kt * (G.transpose() * dR_dp1.transpose() * a + JTT * dvw_dp1.col(idx));
                    }
                }
            }
            // design params 2
            if (_contact_body->_design_params_2._active) {
                for (int k = 0;k < _contact_body->_design_params_2._ndof;k++) {
                    int idx = _contact_body->_design_params_2._param_index(k) + _sim->_ndof_p1;
                    Matrix3 dR_dp2 = _contact_body->_dE0i_dp2(k).topLeftCorner(3, 3);
                    dfm_dp.block(_contact_body->_index[0], idx, 6, 1) -= _kt * (G.transpose() * dR_dp2.transpose() * a + JTT * dvw_dp2.col(k));
                }
            }
            // design params 3
            if (_contact_body->_design_params_3._active) {
                for (int k = 0;k < 3;k++) {
                    int idx = _contact_body->_design_params_3._param_index(contact_id * 3 + k) + _sim->_ndof_p1 + _sim->_ndof_p2;
                    dfm_dp.block(_contact_body->_index[0], idx, 6, 1) -= _kt * (dG_dp3(k).transpose() * R.transpose() * a + JTT * dvw_dp3.col(k));
                }
            }
        } else {
            // Dynamic friction force
            dtype mukn = _mu * _kn;
            Vector3 t = a / anorm;
            Vector3 fd = mukn * d * t;
            fm.segment(_contact_body->_index[0], 6).noalias() += J.transpose() * fd;

            Matrix3 A = (I - t * t.transpose()) / anorm;
            Matrix6 D = mukn * J.transpose() * d * A * T * J;
            Vector3 Rt = R.transpose() * t;
            MatrixX K1(3, 6);
            K1 << e1b * Rt, e2b * Rt, e3b * Rt, Z;
            K1 *= -d;
            MatrixX K2 = R.transpose() * t * ng.transpose() * J;
            MatrixX tmp3(3, 6);
            tmp3 << math::skew(G * phi), Z;
            MatrixX K3 = -d * R.transpose() * A * T * R * tmp3;
            Matrix6 K = mukn * G.transpose() * (K1 + K2 + K3);

            Dm.block(_contact_body->_index[0], _contact_body->_index[0], 6, 6).noalias() += D;
            Km.block(_contact_body->_index[0], _contact_body->_index[0], 6, 6).noalias() += K;

            // design derivatives
            // tool: dt_dp
            MatrixX dt_dp1 = MatrixX::Zero(3, _sim->_ndof_p1), dt_dp2 = MatrixX::Zero(3, _sim->_ndof_p2), dt_dp3 = MatrixX::Zero(3, 3);
            // design params 1
            Matrix3 AT = A * T;
            for (auto ancestor = _contact_body->_joint; ancestor != nullptr; ancestor = ancestor->_parent) {
                if (ancestor->_design_params_1._active) {
                    int idx = ancestor->_design_params_1._param_index(0);
                    int ndof = ancestor->_design_params_1._ndof;
                    dt_dp1.middleCols(idx, ndof) = AT * dvw_dp1.middleCols(idx, ndof);
                }
            }
            // design params 2
            if (_contact_body->_design_params_2._active) {
                dt_dp2 = AT * dvw_dp2;
            }
            // design params 3
            if (_contact_body->_design_params_3._active) {
                dt_dp3 = AT * dvw_dp3;
            }

            // dfd_dp
            // design params 1
            MatrixX JTt = J.transpose() * t;
            MatrixX RTt = R.transpose() * t;
            for (auto ancestor = _contact_body->_joint; ancestor != nullptr; ancestor = ancestor->_parent) {
                if (ancestor->_design_params_1._active) {
                    for (int k = 0;k < ancestor->_design_params_1._ndof;k++) {
                        int idx = ancestor->_design_params_1._param_index(k);
                        Matrix3 dR_dp1 = _contact_body->_dE0i_dp1(idx).topLeftCorner(3, 3);
                        dfm_dp.block(_contact_body->_index[0], idx, 6, 1) += mukn * (dd_dp1(idx) * JTt + d * G.transpose() * dR_dp1.transpose() * t + d * J.transpose() * dt_dp1.col(idx));
                    }
                }
            }
            // design params 2
            if (_contact_body->_design_params_2._active) {
                for (int k = 0;k < _contact_body->_design_params_2._ndof;k++) {
                    int idx = _contact_body->_design_params_2._param_index(k) + _sim->_ndof_p1;
                    Matrix3 dR_dp2 = _contact_body->_dE0i_dp2(k).topLeftCorner(3, 3);
                    dfm_dp.block(_contact_body->_index[0], idx, 6, 1) += mukn * (dd_dp2(k) * JTt + d * G.transpose() * dR_dp2.transpose() * t + d * J.transpose() * dt_dp2.col(k));
                }
            }
            // design params 3
            if (_contact_body->_design_params_3._active) {
                for (int k = 0;k < 3;k++) {
                    int idx = _contact_body->_design_params_3._param_index(contact_id * 3 + k) + _sim->_ndof_p1 + _sim->_ndof_p2;
                    dfm_dp.block(_contact_body->_index[0], idx, 6, 1) += mukn * (dd_dp3(k) * JTt + d * dG_dp3(k).transpose() * RTt + d * J.transpose() * dt_dp3.col(k));
                }
            }

        }
    }
}

// void ForceGroundContact::computeForce(std::vector<Contact> &contacts, VectorX& fm) {
//     Vector3 xg = _E_g.topRightCorner(3, 1);
//     Vector3 ng = _E_g.block(0, 2, 3, 1);
//     Matrix3 N = ng * ng.transpose();
//     Matrix3 T = Matrix3::Identity() - N;
//     Matrix3 R = (_contact_body->_E_0i).topLeftCorner(3, 3);
//     Vector6 phi = _contact_body->_phi;

//     for (int i = 0;i < contacts.size();i++) {
//         Vector3 xi = contacts[i]._xi;
//         // Vector3 xw = contacts[i]._xw;
//         dtype d = contacts[i]._d;
//         MatrixX xibrac = math::skew(xi);

//         // Contact force
//         MatrixX G = math::gamma(xi);
//         Vector3 vwi = R * (G * phi);
//         Vector3 fc = -_kn * ng * d + _damping * N * vwi * d;
        
//         fm.segment(_contact_body->_index[0], 6).noalias() += G.transpose() * (R.transpose() * fc);
        
//         // Friction force
//         if (_mu < constants::eps)
//             continue;
        
//         Vector3 a = T * vwi;
//         dtype anorm = a.norm();
//         if (_mu * fabs(_kn * d) >= _kt * anorm - constants::eps) {
//             // Static friction force
//             Vector3 fs = -_kt * a;
//             // std::cerr << "fs = " << fs.transpose() << std::endl;
//             fm.segment(_contact_body->_index[0], 6).noalias() += G.transpose() * (R.transpose() * fs);
//         } else {
//             dtype mukn = _mu * _kn;
//             Vector3 t = a / anorm;
//             Vector3 fd = mukn * d * t;
//             // std::cerr << "fd = " << fd.transpose() << std::endl;
//             fm.segment(_contact_body->_index[0], 6).noalias() += G.transpose() * (R.transpose() * fd);
//         }
//     }   
// }

// void ForceGroundContact::computeForceWithDerivative(std::vector<Contact> &contacts, VectorX& fm, MatrixX& Km, MatrixX& Dm) {
//     Vector3 xg = _E_g.topRightCorner(3, 1);
//     Vector3 ng = _E_g.block(0, 2, 3, 1);
//     Matrix3 N = ng * ng.transpose();
//     Matrix3 I = Matrix3::Identity();
//     Matrix3 Z = Matrix3::Zero();
//     Matrix3 T = I - N;
//     Matrix3 R = (_contact_body->_E_0i).topLeftCorner(3, 3);
//     Vector3 p = (_contact_body->_E_0i).topRightCorner(3, 1);
//     Vector6 phi = _contact_body->_phi;
//     Matrix3 e1b = math::skew(Vector3::UnitX());
//     Matrix3 e2b = math::skew(Vector3::UnitY());
//     Matrix3 e3b = math::skew(Vector3::UnitZ());
//     Matrix3 RNR = R.transpose() * N * R;
//     Matrix3 pxgtmp = math::skew(R.transpose() * N * (p - xg));
    
//     for (int i = 0;i < contacts.size();i++) {
//         Vector3 xi = contacts[i]._xi;
//         Vector3 xw = contacts[i]._xw;
//         dtype d = contacts[i]._d;
//         Vector3 RNRxi = RNR * xi;
//         MatrixX xibrac = math::skew(xi);
//         MatrixX G = math::gamma(xi);
//         Vector3 Gphi = G * phi;
//         Vector3 vwi = R * Gphi;
//         Vector3 RNRGphi = RNR * Gphi;
//         Matrix3 Gphibrac = math::skew(Gphi);

//         // Contact force
        
//         MatrixX J = R * G;

//         Vector3 fc = -_kn * ng * d + _damping * N * vwi * d;
//         fm.segment(_contact_body->_index[0], 6).noalias() += J.transpose() * fc;

//         Matrix3 tmp1 = Matrix3::Zero();
//         tmp1 << -e1b * RNRxi, -e2b * RNRxi, -e3b * RNRxi;
//         tmp1.noalias() += - RNR * xibrac + pxgtmp;

//         Matrix3 tmp2 = Matrix3::Zero();
//         tmp2 << -e1b * RNRGphi, -e2b * RNRGphi, -e3b * RNRGphi;
//         tmp2.noalias() += -RNR * Gphibrac;
//         Km.block(_contact_body->_index[0], _contact_body->_index[0], 6, 6) -= 
//             _kn * G.transpose() * (Matrix<dtype, 3, 6>() << tmp1, RNR).finished()
//             - _damping * G.transpose() * (Matrix<dtype, 3, 6>() << tmp2, Z).finished() * d
//             - _damping * J.transpose() * N * vwi * ng.transpose() * (Matrix<dtype, 3, 6>() << -R * xibrac, R).finished();
//         Dm.block(_contact_body->_index[0], _contact_body->_index[0], 6, 6) +=
//             _damping * G.transpose() * RNR * G * d;

//         // Friction force
//         if (_mu < constants::eps)
//             continue;
        
//         Vector3 xwdot = J * phi;
//         Vector3 a = T * xwdot;
//         dtype anorm = a.norm();
//         if (_mu * fabs(_kn * d) >= _kt * anorm - constants::eps) {
//             // Static friction force
//             Vector3 fs = -_kt * a;
//             fm.segment(_contact_body->_index[0], 6).noalias() += J.transpose() * fs;
            
//             Matrix6 D = -_kt * J.transpose() * T * J;
//             Matrix3 B = R.transpose() * T * R;
//             MatrixX tmp3(3, 6);
//             tmp3 << (B * e1b - e1b * B) * Gphi, (B * e2b - e2b * B) * Gphi, (B * e3b - e3b * B) * Gphi, Z;
//             Matrix6 K = -_kt * G.transpose() * tmp3;

//             Dm.block(_contact_body->_index[0], _contact_body->_index[0], 6, 6).noalias() += D;
//             Km.block(_contact_body->_index[0], _contact_body->_index[0], 6, 6).noalias() += K;
//         } else {
//             // Dynamic friction force
//             dtype mukn = _mu * _kn;
//             Vector3 t = a / anorm;
//             Vector3 fd = mukn * d * t;
//             fm.segment(_contact_body->_index[0], 6).noalias() += J.transpose() * fd;

//             Matrix3 A = (I - t * t.transpose()) / anorm;
//             Matrix6 D = mukn * J.transpose() * d * A * T * J;
//             Vector3 Rt = R.transpose() * t;
//             MatrixX K1(3, 6);
//             K1 << e1b * Rt, e2b * Rt, e3b * Rt, Z;
//             K1 *= -d;
//             MatrixX K2 = R.transpose() * t * ng.transpose() * J;
//             MatrixX tmp3(3, 6);
//             tmp3 << math::skew(G * phi), Z;
//             MatrixX K3 = -d * R.transpose() * A * T * R * tmp3;
//             Matrix6 K = mukn * G.transpose() * (K1 + K2 + K3);

//             Dm.block(_contact_body->_index[0], _contact_body->_index[0], 6, 6).noalias() += D;
//             Km.block(_contact_body->_index[0], _contact_body->_index[0], 6, 6).noalias() += K;
//         }
//     }
// }

void ForceGroundContact::test_derivatives_runtime() {
    test_design_derivatives_runtime();
}

void ForceGroundContact::test_design_derivatives_runtime() {
    // detect the contact points between two bodies
    std::vector<Contact> contacts;
    contacts.clear();
    collision_detection_ground_object(_E_g, _contact_body, contacts);

    if (contacts.size() > 0) {
        VectorX fm = VectorX::Zero(_sim->_ndof_m);
        MatrixX Km = MatrixX::Zero(_sim->_ndof_m, _sim->_ndof_m), Dm = MatrixX::Zero(_sim->_ndof_m, _sim->_ndof_m);
        MatrixX dfm_dp = MatrixX::Zero(_sim->_ndof_m, _sim->_ndof_p);
        computeForceWithDerivative(contacts, fm, Km, Dm, dfm_dp, false);

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

                VectorX fm_pos = VectorX::Zero(_sim->_ndof_m);
                MatrixX Km_pos = MatrixX::Zero(_sim->_ndof_m, _sim->_ndof_m), Dm_pos = MatrixX::Zero(_sim->_ndof_m, _sim->_ndof_m);
                MatrixX dfm_dp_pos = MatrixX::Zero(_sim->_ndof_m, _sim->_ndof_p);
                computeForceWithDerivative(contacts_pos, fm_pos, Km_pos, Dm_pos, dfm_dp_pos, false);
                dfm_dp_fd.col(i) = (fm_pos - fm) / eps;
            }
            // std::cerr << "eps = " << eps << std::endl;
            print_error("ForceGroundContact: dfm_dp1", dfm_dp.leftCols(_sim->_ndof_p1), dfm_dp_fd.leftCols(_sim->_ndof_p1));
            print_error("ForceGroundContact: dfm_dp2", dfm_dp.middleCols(_sim->_ndof_p1, _sim->_ndof_p2), dfm_dp_fd.middleCols(_sim->_ndof_p1, _sim->_ndof_p2));
            print_error("ForceGroundContact: dfm_dp3", dfm_dp.middleCols(_sim->_ndof_p1 + _sim->_ndof_p2, _sim->_ndof_p3), dfm_dp_fd.middleCols(_sim->_ndof_p1 + _sim->_ndof_p2, _sim->_ndof_p3));
            eps /= 10.;
        }

        _sim->set_design_params(design_params);
        _sim->update_robot(true);
    }
}

}