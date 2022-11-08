#include "Body/Body.h"
#include "Utils.h"
#include "Joint/Joint.h"
#include "Simulation.h"
#include "Robot.h"

namespace redmax {

const Eigen::Matrix4f BodyAnimator::AnimatedModelMatrix(const float t) {
    Matrix4f model_matrix = _body->_E_0i.cast<float>();
    if (_body->_sim->_options->_unit == "cm-g")
        model_matrix.topRightCorner(3, 1) /= 10.; // scale for better visualization
    else
        model_matrix.topRightCorner(3, 1) *= 10.; // scale for better visualization
    return model_matrix;
}

Body::Body(Simulation *sim, Joint *joint, dtype density): 
    _sim(sim),
    _design_params_2(this, DesignParameterType::TYPE2, false),
    _design_params_3(this, DesignParameterType::TYPE3, false),
    _design_params_4(this, DesignParameterType::TYPE4, false),
    _design_params_6(this, DesignParameterType::TYPE6, false) {
    
    _joint = joint;
    _joint->_body = this;
    if (_joint->_parent == nullptr)
        _parent = nullptr;
    else
        _parent = _joint->_parent->_body;

    _density = density;

    _contact_points.clear();

    _animator = nullptr;

    _color << 0.25f, 0.148f, 0.06475f;

    _use_texture = false;

    _contact_scale = 1.;

    _f_external.setZero();
}

Body::Body(Simulation *sim, Joint *joint, Matrix3 R_ji, Vector3 p_ji, dtype density): 
    _sim(sim),
    _design_params_2(this, DesignParameterType::TYPE2, false),
    _design_params_3(this, DesignParameterType::TYPE3, false),
    _design_params_4(this, DesignParameterType::TYPE4, false),
    _design_params_6(this, DesignParameterType::TYPE6, false) {

    _joint = joint;
    _joint->_body = this;
    if (_joint->_parent == nullptr)
        _parent = nullptr;
    else
        _parent = _joint->_parent->_body;

    _density = density;

    _E_ji = math::SE(R_ji, p_ji);
    _E_ij = math::Einv(_E_ji);
    _A_ji = math::Ad(_E_ji);
    _A_ij = math::Ad(_E_ij);
    
    _contact_points.clear();

    _animator = nullptr;

    _color << 0.25f, 0.148f, 0.06475f;

    _use_texture = false;

    _contact_scale = 1.;

    _f_external.setZero();
}

Body::~Body() {
    if (_animator) {
        delete _animator;
    }
}

void Body::set_transform(Matrix3 R_ji, Vector3 p_ji) {
    _E_ji = math::SE(R_ji, p_ji);
    _E_ij = math::Einv(_E_ji);
    _A_ji = math::Ad(_E_ji);
    _A_ij = math::Ad(_E_ij);
}

void Body::set_contacts(vector<Vector3>& contacts) {
    _contact_points = contacts;
}

void Body::init() {
    if (_joint->_ndof > 0) {
        _dAip_dq = JacobianMatrixVector(6, 6, _joint->_ndof);
        _dAipdot_dq = JacobianMatrixVector(6, 6, _joint->_ndof);
    }

    if (_sim->_robot->_ndof_p1 > 0) {
        _dE0i_dp1 = JacobianMatrixVector(4, 4, _sim->_robot->_ndof_p1);
        _dphi_dp1 = MatrixX::Zero(6, _sim->_robot->_ndof_p1);
    }

    if (_design_params_2._active && _design_params_2._ndof > 0) {
        _dAij_dp2 = JacobianMatrixVector(6, 6, _design_params_2._ndof);
        _dE0i_dp2 = JacobianMatrixVector(4, 4, _design_params_2._ndof);
        _dphi_dp2 = MatrixX::Zero(6, _design_params_2._ndof);
    }
}

void Body::set_color(Vector3 color) {
    _color = color.cast<float>();
    _use_texture = false;
}

void Body::set_texture(std::string texture_path) {
    _texture_path = texture_path;
    _use_texture = true;
}

Vector6 Body::get_external_force() const {
    return _f_external;
}

void Body::set_external_force(Vector6 force) {
    _f_external = force;
}

void Body::reset_external_force() {
    _f_external.setZero();
}

void Body::update(bool design_gradient) {
    // _E_0i
    _E_0i = _joint->_E_0j * _E_ji;
    _E_i0 = math::Einv(_E_0i);

    // _E_ip = _E_i0 * _E_0p
    if (_parent != nullptr) {
        _E_ip = _E_i0 * _parent->_E_0i;
    } else {
        _E_ip = SE3::Identity();
    }
    _E_pi = math::Einv(_E_ip);

    // _Ad_ip
    _A_ip = math::Ad(_E_ip);
    _A_pi = math::Ad(_E_pi);

    // _phi, based on noted eq. (89)
    _phi = _A_ij * _joint->_phi;

    /********************************* Derivatives ************************************/
    if (_parent != nullptr) {
        // derivatives: use the following formulas to derive the gradients:
        // Ad_BiBp = Ad0_BiJi * Ad(inv(Q_Ji)) * Ad0_Ji'Jp * Ad0_JpBp
        // Ad_dot_BiBp = -Ad0_BiJi * inv(Ad_Q) * Ad_Q_dot * inv(Ad_Q) * Ad0_Ji'Jp * Ad0_JpBp
        Matrix6 Aleft = -math::Ad(_E_ij * _joint->_Q_inv);
        Matrix4 E_jp_0 = _joint->_E_jp_0 * _parent->_E_ji; // transformation from Ji' to Bp
        Matrix6 Aright = math::Ad(_joint->_Q_inv * E_jp_0);

        _A_ip_dot = Aleft * _joint->_Adot * Aright;
        _dAip_dq.setZero();
        _dAipdot_dq.setZero();
        for (int k = 0;k < _joint->_ndof;k++) {
            _dAip_dq(k) = Aleft * _joint->_dA_dq(k) * Aright;
            Matrix6 tmp1 = _joint->_dA_dq(k) * _joint->_A_inv * _joint->_Adot;
            Matrix6 tmp2 = _joint->_Adot * _joint->_A_inv * _joint->_dA_dq(k);
            _dAipdot_dq(k) = Aleft * (_joint->_dAdot_dq(k) - tmp1 - tmp2) * Aright;
        }
    }

    if (design_gradient) {
        update_design_derivatives();
    }
}

void Body::update_design_derivatives() {
    // dE0i_dp1, dphi_dp1
    for (auto ancestor = _joint; ancestor != nullptr; ancestor = ancestor->_parent) {
        if (ancestor->_design_params_1._active) {
            for (int k = 0;k < ancestor->_design_params_1._ndof;k++) {
                int idx = ancestor->_design_params_1._param_index(k);
                _dE0i_dp1(idx) = _joint->_dE0j_dp1(idx) * _E_ji;
                _dphi_dp1.col(idx) = _A_ij * _joint->_dphi_dp1.col(idx);
            }
        }
    }

    if (_design_params_2._active) {
        // dAij_dp2
        Matrix3 RT = _E_ji.topLeftCorner(3, 3).transpose();
        Vector3 p = _E_ji.topRightCorner(3, 1);
        Matrix3 tmp = math::skew(-RT * p);
        for (int k = 0;k < 3;k++)
            for (int l = 0;l < 3;l++) {
                int idx = k * 3 + l;
                _dAij_dp2(idx)(l, k) = 1.;
                _dAij_dp2(idx)(l + 3, k + 3) = 1.;
                _dAij_dp2(idx).bottomLeftCorner(3, 3) = math::skew(-Vector3::Unit(l) * p(k)) * RT + (tmp * Vector3::Unit(l)) * RowVector3::Unit(k);
            }
        for (int k = 0;k < 3;k++) {
            int idx = 9 + k;
            Vector3 ek = Vector3::Unit(k);
            _dAij_dp2(idx).bottomLeftCorner(3, 3) = math::skew(-RT * ek) * RT;
        }

        // dE0i_dp2
        for (int k = 0;k < 3;k++) 
            for (int l = 0;l < 4;l++) {
                int idx;
                if (l < 3)
                    idx = k * 3 + l;
                else
                    idx = 9 + k;
                _dE0i_dp2(idx).setZero();
                _dE0i_dp2(idx).col(l) = _joint->_E_0j.col(k);
            }

        // dphi_dp2
        for (int k = 0;k < _design_params_2._ndof;k++) {
            _dphi_dp2.col(k) = _dAij_dp2(k) * _joint->_phi;
        }
    }
}

// activate design parameters
void Body::activate_design_parameters_type_2(bool active) {
    _design_params_2.activate(active);
}

void Body::activate_design_parameters_type_3(bool active) {
    int ndof = _contact_points.size() * 3;
    _design_params_3.activate(active, ndof);
}

void Body::activate_design_parameters_type_4(bool active) {
    _design_params_4.activate(active);
}

void Body::activate_design_parameters_type_6(bool active) {
    _design_params_6.activate(active);
}

void Body::computeMaximalForce(VectorX& fm) {
    Matrix3 R_0i = _E_0i.topLeftCorner(3, 3);

    // Coriolis Force
    Matrix6 ad_trans = math::ad(_phi).transpose();
    Vector6 f_coriolis = ad_trans * (_Inertia.cwiseProduct(_phi));

    // gravity
    Vector6 f_gravity;
    f_gravity.setZero();
    f_gravity.tail(3) = R_0i.transpose() * (_mass * _sim->_options->_gravity);

    // external force
    Vector6 f_external;
    f_external.setZero();
    // f_external.head(3) = R_0i.transpose() * _f_external.head(3);
    // f_external.tail(3) = R_0i.transpose() * _f_external.tail(3);
    f_external.head(3) = _Inertia.head(3).cwiseProduct(R_0i.transpose() * _f_external.head(3));
    f_external.tail(3) = _Inertia.tail(3).cwiseProduct(R_0i.transpose() * _f_external.tail(3));

    fm.segment(_index[0], 6) += f_coriolis + f_gravity + f_external;
}

// Km = dfm_dq
// Dm = dfm_dphi
void Body::computeMaximalForceWithDerivative(VectorX& fm, MatrixX& Km, MatrixX& Dm) {
    Matrix3 R_0i = _E_0i.topLeftCorner(3, 3);

    // Coriolis Force
    Matrix6 ad_trans = math::ad(_phi).transpose();
    Vector6 f_coriolis = ad_trans * (_Inertia.cwiseProduct(_phi));

    // gravity
    Vector6 f_gravity;
    f_gravity.setZero();
    f_gravity.tail(3) = R_0i.transpose() * (_mass * _sim->_options->_gravity);

    // external force
    Vector6 f_external;
    f_external.setZero();
    // f_external.head(3) = R_0i.transpose() * _f_external.head(3);
    // f_external.tail(3) = R_0i.transpose() * _f_external.tail(3);
    f_external.head(3) = _Inertia.head(3).cwiseProduct(R_0i.transpose() * _f_external.head(3));
    f_external.tail(3) = _Inertia.tail(3).cwiseProduct(R_0i.transpose() * _f_external.tail(3));

    fm.segment(_index[0], 6) += f_coriolis + f_gravity + f_external;

    /********************************* Derivatives ************************************/
    // Km comes from gravity and external force
    // Km.block(_index[0], _index[0], 3, 3) += math::skew(f_external.head(3));
    // Km.block(_index[3], _index[0], 3, 3) += math::skew(f_gravity.tail(3) + f_external.tail(3));
    Km.block(_index[3], _index[0], 3, 3) += math::skew(f_gravity.tail(3));
    MatrixX df_dq = MatrixX::Zero(6, 3);
    df_dq.topRows(3) = math::skew(R_0i.transpose() * _f_external.head(3));
    df_dq.bottomRows(3) = math::skew(R_0i.transpose() * _f_external.tail(3));
    for (int i = 0;i < 6;++i)
        df_dq.row(i) = df_dq.row(i) * _Inertia(i);
    Km.block(_index[0], _index[0], 6, 3) += df_dq;

    // Dm comes from Coriolis
    Matrix3 e1 = math::skew(Vector3::UnitX());
    Matrix3 e2 = math::skew(Vector3::UnitY());
    Matrix3 e3 = math::skew(Vector3::UnitZ());
    Vector3 Iw = _Inertia.head(3).cwiseProduct(_phi.head(3));
    Vector3 mv = _Inertia.tail(3).cwiseProduct(_phi.tail(3));

    // Dm += ad' * I
    for (int k = 0;k < 6;k++)
        Dm.block(_index[0], _index[0] + k, 6, 1) += ad_trans.col(k) * _Inertia(k);

    // Rest part in Dm (eq (1.39))
    Dm.block(_index[0], _index[0], 3, 1) -= e1 * Iw;
    Dm.block(_index[0], _index[1], 3, 1) -= e2 * Iw;
    Dm.block(_index[0], _index[2], 3, 1) -= e3 * Iw;
    Dm.block(_index[0], _index[3], 3, 1) -= e1 * mv;
    Dm.block(_index[0], _index[4], 3, 1) -= e2 * mv;
    Dm.block(_index[0], _index[5], 3, 1) -= e3 * mv;
    Dm.block(_index[3], _index[0], 3, 1) -= e1 * mv;
    Dm.block(_index[3], _index[1], 3, 1) -= e2 * mv;
    Dm.block(_index[3], _index[2], 3, 1) -= e3 * mv;
}

void Body::computeMaximalForceWithDerivative(
    VectorX& fm, 
    MatrixX& Km, MatrixX& Dm,
    MatrixX& dfm_dp) {
    
    Matrix3 R_0i = _E_0i.topLeftCorner(3, 3);

    // Coriolis Force
    Matrix6 ad_trans = math::ad(_phi).transpose();
    Vector6 Mmphi = _Inertia.cwiseProduct(_phi);
    Vector6 f_coriolis = ad_trans * Mmphi;
    Matrix6 adTMm;
    for (int i = 0;i <6;i ++)
        adTMm.col(i) = ad_trans.col(i) * _Inertia(i);

    // gravity
    Vector6 f_gravity;
    f_gravity.setZero();
    f_gravity.tail(3) = R_0i.transpose() * (_mass * _sim->_options->_gravity);

    // external force
    Vector6 f_external;
    f_external.setZero();
    // f_external.head(3) = R_0i.transpose() * _f_external.head(3);
    // f_external.tail(3) = R_0i.transpose() * _f_external.tail(3);
    f_external.head(3) = _Inertia.head(3).cwiseProduct(R_0i.transpose() * _f_external.head(3));
    f_external.tail(3) = _Inertia.tail(3).cwiseProduct(R_0i.transpose() * _f_external.tail(3));

    fm.segment(_index[0], 6) += f_coriolis + f_gravity + f_external;

    /********************************* Derivatives ************************************/
    // Km comes from gravity and external force
    // Km.block(_index[0], _index[0], 3, 3) += math::skew(f_external.head(3));
    // Km.block(_index[3], _index[0], 3, 3) += math::skew(f_gravity.tail(3) + f_external.tail(3));
    Km.block(_index[3], _index[0], 3, 3) += math::skew(f_gravity.tail(3));
    MatrixX df_dq = MatrixX::Zero(6, 3);
    df_dq.topRows(3) = math::skew(R_0i.transpose() * _f_external.head(3));
    df_dq.bottomRows(3) = math::skew(R_0i.transpose() * _f_external.tail(3));
    for (int i = 0;i < 6;++i)
        df_dq.row(i) = df_dq.row(i) * _Inertia(i);
    Km.block(_index[0], _index[0], 6, 3) += df_dq;

    // Dm comes from Coriolis
    Matrix3 e1 = math::skew(Vector3::UnitX());
    Matrix3 e2 = math::skew(Vector3::UnitY());
    Matrix3 e3 = math::skew(Vector3::UnitZ());
    Vector3 Iw = _Inertia.head(3).cwiseProduct(_phi.head(3));
    Vector3 mv = _Inertia.tail(3).cwiseProduct(_phi.tail(3));

    // Dm += ad' * I
    for (int k = 0;k < 6;k++)
        Dm.block(_index[0], _index[0] + k, 6, 1) += ad_trans.col(k) * _Inertia(k);

    // Rest part in Dm (eq (1.39))
    Dm.block(_index[0], _index[0], 3, 1) -= e1 * Iw;
    Dm.block(_index[0], _index[1], 3, 1) -= e2 * Iw;
    Dm.block(_index[0], _index[2], 3, 1) -= e3 * Iw;
    Dm.block(_index[0], _index[3], 3, 1) -= e1 * mv;
    Dm.block(_index[0], _index[4], 3, 1) -= e2 * mv;
    Dm.block(_index[0], _index[5], 3, 1) -= e3 * mv;
    Dm.block(_index[3], _index[0], 3, 1) -= e1 * mv;
    Dm.block(_index[3], _index[1], 3, 1) -= e2 * mv;
    Dm.block(_index[3], _index[2], 3, 1) -= e3 * mv;

    // design derivatives

    // f_gravity + f_external
    // design params 1
    for (auto ancestor = _joint;ancestor != nullptr;ancestor = ancestor->_parent) {
        if (ancestor->_design_params_1._active) {
            for (int k = 0;k < ancestor->_design_params_1._ndof;k++) {
                int idx = ancestor->_design_params_1._param_index(k);
                dfm_dp.block(_index[0], idx, 3, 1) += _dE0i_dp1(idx).topLeftCorner(3, 3).transpose() * _f_external.head(3);
                dfm_dp.block(_index[0] + 3, idx, 3, 1) += _dE0i_dp1(idx).topLeftCorner(3, 3).transpose() * (_mass * _sim->_options->_gravity + _f_external.tail(3));
            }
        }
    }
    // design params 2
    if (_design_params_2._active) {
        for (int k = 0;k < _design_params_2._ndof;k++) {
            int idx = _design_params_2._param_index(k);
            dfm_dp.block(_index[0], idx + _sim->_ndof_p1, 3, 1) += _dE0i_dp2(k).topLeftCorner(3, 3).transpose() * _f_external.head(3);
            dfm_dp.block(_index[0] + 3, idx + _sim->_ndof_p1, 3, 1) += _dE0i_dp2(k).topLeftCorner(3, 3).transpose() * (_mass * _sim->_options->_gravity + _f_external.tail(3));
        }
    }
    // design params 4
    int p4_offset = _sim->_ndof_p1 + _sim->_ndof_p2 + _sim->_ndof_p3;
    if (_design_params_4._active) {
        int idx = _design_params_4._param_index[0]; // mass
        dfm_dp.block(_index[0] + 3, idx + p4_offset, 3, 1) += _E_0i.topLeftCorner(3, 3).transpose() * _sim->_options->_gravity;
    }

    // f_coriolis
    // da(phi)'_dphi
    JacobianMatrixVector daT_dphi(6, 6, 6);
    for (int i = 0;i < 3;i++) {
        Matrix3 tmp = -math::skew(Vector3::Unit(i));
        daT_dphi(i).topLeftCorner(3, 3) = tmp;
        daT_dphi(i).bottomRightCorner(3, 3) = tmp;
        daT_dphi(i + 3).topRightCorner(3, 3) = tmp;
    }

    // design params 1
    for (auto ancestor = _joint;ancestor != nullptr;ancestor = ancestor->_parent) {
        if (ancestor->_design_params_1._active) {
            for (int k = 0;k < ancestor->_design_params_1._ndof;k++) {
                int idx = ancestor->_design_params_1._param_index(k);
                Matrix6 daT_dp1 = Matrix6::Zero();
                for (int j = 0;j < 6;j++)
                    daT_dp1 += daT_dphi(j) * _dphi_dp1(j, idx);
                dfm_dp.block(_index[0], idx, 6, 1) += daT_dp1 * Mmphi + adTMm * _dphi_dp1.col(idx);
            }
        }
    }

    // design params 2
    if (_design_params_2._active) {
        for (int k = 0;k < _design_params_2._ndof;k++) {
            int idx = _design_params_2._param_index[k];
            Matrix6 daT_dp2 = Matrix6::Zero();
            for (int j = 0;j < 6;j++)
                daT_dp2 += daT_dphi(j) * _dphi_dp2(j, k);
            dfm_dp.block(_index[0], idx + _sim->_ndof_p1, 6, 1) += daT_dp2 * Mmphi + adTMm * _dphi_dp2.col(k);
        }
    }

    // design params 4
    if (_design_params_4._active) {
        int idx = _design_params_4._param_index[0];
        // mass
        dfm_dp.block(_index[0], idx + p4_offset, 6, 1) += ad_trans.rightCols(3) * _phi.tail(3);
        // inertia
        for (int i = 0;i < 3;i++) {
            dfm_dp.block(_index[0], idx + 1 + i + p4_offset, 6, 1) += ad_trans.col(i) * _phi(i);
        }
    }
}

void Body::test_derivatives_runtime() {
    if (_joint->_ndof > 0 && _parent != nullptr) {
        // test _dAip_dq, _dAipdot_dq
        _joint->update();
        auto dAip_dq = _dAip_dq;
        auto dAipdot_dq = _dAipdot_dq;
        auto q = _joint->_q;
        auto qdot = _joint->_qdot;
        auto E_ip = _E_ip;
        auto A_ip = _A_ip;
        auto A_ip_dot = _A_ip_dot;

        std::string body_str = "Body " + std::to_string(_joint->_id);

        dtype h = 1e-8;
        JacobianMatrixVector dAip_dq_fd(6, 6, _joint->_ndof);
        JacobianMatrixVector dAipdot_dq_fd(6, 6, _joint->_ndof);
        for (int k = 0;k < _joint->_ndof;k++) {
            _joint->_q(k) += h;
            _joint->update();
            auto A_ip_pos = _A_ip;
            auto A_ip_dot_pos = _A_ip_dot;
            dAip_dq_fd(k) = (A_ip_pos - A_ip) / h;
            dAipdot_dq_fd(k) = (A_ip_dot_pos - A_ip_dot) / h;
            _joint->_q(k) -= h;
        }
        print_error(body_str + ": dAip_dq", dAip_dq, dAip_dq_fd);
        print_error(body_str + ": dAipdot_dq", dAipdot_dq, dAipdot_dq_fd);

        _joint->_q = q;
        _joint->_qdot = qdot;
        _joint->update();
    }

    {
        // test Km
        dtype h = 1e-8;
        std::string body_str = "Body " + std::to_string(_joint->_id);

        MatrixX Km = MatrixX::Zero(_sim->_ndof_m, _sim->_ndof_m);
        MatrixX Dm = MatrixX::Zero(_sim->_ndof_m, _sim->_ndof_m);
        VectorX fm = VectorX::Zero(_sim->_ndof_m);
        computeMaximalForceWithDerivative(fm, Km, Dm);

        MatrixX Km_fd = MatrixX::Zero(_sim->_ndof_m, _sim->_ndof_m);
        
        Matrix4 E = _E_0i;
        for (int k = 0;k < 6;++k) {
            Vector6 dq = Vector6::Zero();
            dq[k] = h;
            Matrix4 E_pos = E * math::exp(dq);
            _E_0i = E_pos;
            VectorX fm_pos = VectorX::Zero(_sim->_ndof_m);
            computeMaximalForce(fm_pos);
            Km_fd.col(_index[k]) = (fm_pos - fm) / h;
            _E_0i = E;
        }

        print_error(body_str + ": Maximal Force Km", Km, Km_fd);

        MatrixX Dm_fd = MatrixX::Zero(_sim->_ndof_m, _sim->_ndof_m);
        Vector6 phi = _phi;
        for (int k = 0;k < 6;++k) {
            Vector6 phi_pos = phi;
            phi_pos[k] += h;
            _phi = phi_pos;
            VectorX fm_pos = VectorX::Zero(_sim->_ndof_m);
            computeMaximalForce(fm_pos);
            Dm_fd.col(_index[k]) = (fm_pos - fm) / h;
            _phi = phi;
        }

        print_error(body_str + ": Maximal Force Dm", Dm, Dm_fd);

        _joint->update();
    }
}

Vector3 Body::position_in_world(Vector3 pos) const {
    Matrix3 R = _E_0i.topLeftCorner(3, 3);
    Vector3 p;
    p[0] = _E_0i(0, 3);
    p[1] = _E_0i(1, 3);
    p[2] = _E_0i(2, 3);
    return R * pos + p;
}

Vector3 Body::velocity_in_world(Vector3 pos) const {
    Matrix3 R = _E_0i.topLeftCorner(3, 3);
    Vector3 vel = R * (math::skew(pos).transpose() * _phi.head(3) + _phi.tail(3));
    return vel;
}

void Body::add_contact_body(std::string body_name) {
    _contact_body_history.insert(body_name);
}

std::vector<std::string> Body::get_contact_bodies() const {
    std::vector<std::string> contact_bodies;
    for (auto itr = _contact_body_history.begin(); itr != _contact_body_history.end(); itr++) {
        contact_bodies.push_back(*itr);
    }
    return contact_bodies;
}

void Body::clear_contact_bodies() {
    _contact_body_history.clear();
}

std::pair<Vector3, Vector3> Body::get_AABB() {
    std::pair<Vector3, Vector3> aabb;
    double p[2][3];
    p[0][0] = _bounding_box.first[0], p[0][1] = _bounding_box.first[1], p[0][2] = _bounding_box.first[2];
    p[1][0] = _bounding_box.second[0], p[1][1] = _bounding_box.second[1], p[1][2] = _bounding_box.second[2];
    for (int x = 0;x < 2;++x)
        for (int y = 0;y < 2;++y)
            for (int z = 0;z < 2;++z) {
                Vector3 vert(p[x][0], p[y][1], p[z][2]);
                Vector3 vert_world = position_in_world(vert);
                if (x == 0 && y == 0 && z == 0) {
                    aabb.first = vert_world;
                    aabb.second = vert_world;
                } else {
                    for (int axis = 0;axis < 3;++axis) {
                        aabb.first[axis] = min(aabb.first[axis], vert_world[axis]);
                        aabb.second[axis] = max(aabb.second[axis], vert_world[axis]);
                    }
                }
            }
    return aabb;
}

}