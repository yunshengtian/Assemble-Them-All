#include "Joint/Joint.h"
#include "Body/Body.h"
#include "Simulation.h"
#include "Robot.h"

namespace redmax {

Joint::Joint(Simulation *sim, int id, Joint* parent, Matrix3 R, Vector3 p, int ndof, Frame frame):
    _design_params_1(this, DesignParameterType::TYPE1, false),
    _design_params_5(this, DesignParameterType::TYPE5, false) {
    
    _sim = sim;

    _id = id;

    _parent = parent;
    if (_parent != nullptr) {
        _parent->_children.push_back(this);
    }

    _body = nullptr;
    
    if (frame == Frame::WORLD) {
        // E_pj_0 = E_p0_0 * E_0j_0
        if (_parent != nullptr) {
            _E_pj_0 = _parent->_E_j0_0 * math::SE(R, p);
        } else {
            _E_pj_0 = math::SE(R, p);
        }
    } else if (frame == Frame::LOCAL) {
        _E_pj_0 = math::SE(R, p);
    }

    _E_jp_0 = math::Einv(_E_pj_0);
    if (_parent == nullptr) {
        _E_j0_0 = _E_jp_0;
    } else {
        _E_j0_0 = _E_jp_0 * _parent->_E_j0_0;
    }

    _ndof = ndof;
    _children.clear();

    _Kr = _Dr = 0;
    _joint_limit_lower = (dtype)INT_MIN;
    _joint_limit_upper = (dtype)INT_MAX;
    _joint_limit_stiffness = 0.;

    if (_ndof > 0) {
        _q_rest = VectorX::Zero(_ndof);
        _q = VectorX::Zero(_ndof);
        _qdot = VectorX::Zero(_ndof);
        _S_j = MatrixX::Zero(6, _ndof);
        _S_j_dot = MatrixX::Zero(6, _ndof);
        _dSj_dq = JacobianMatrixVector(6, _ndof, _ndof);
        _dSjdot_dq = JacobianMatrixVector(6, _ndof, _ndof);
        _dA_dq = JacobianMatrixVector(6, 6, _ndof);
        _dAdot_dq = JacobianMatrixVector(6, 6, _ndof);
        _dQ_dq = JacobianMatrixVector(4, 4, _ndof);
        _dEpj_dq = JacobianMatrixVector(4, 4, _ndof);
    }
}

void Joint::init() {
    if (_sim->_robot->_ndof_p1 > 0) {
        _dAjp0_dp1 = JacobianMatrixVector(6, 6, 12);
        _dE0j_dp1 = JacobianMatrixVector(4, 4, _sim->_robot->_ndof_p1);
        _dphi_dp1 = MatrixX::Zero(6, _sim->_robot->_ndof_p1);
    }
}

void Joint::set_damping(dtype damping) {
    _Dr = damping;
    if (_design_params_5._active) {
        _design_params_5._params[1] = _Dr;
    }
}

void Joint::set_stiffness(dtype stiffness, VectorX q_rest) {
    assert(q_rest.size() == _ndof);
    _Kr = stiffness;
    _q_rest = q_rest;
    _q = q_rest;
    if (_design_params_5._active) {
        _design_params_5._params[0] = _Kr;
    }
}

void Joint::set_joint_limit(dtype lower, dtype upper, dtype joint_limit_stiffness) {
    _joint_limit_lower = lower;
    _joint_limit_upper = upper;
    _joint_limit_stiffness = joint_limit_stiffness;
}

void Joint::update(bool design_gradient) {
    // Q_inv, A_inv
    _Q_inv = math::Einv(_Q);
    _A_inv = math::Ad(_Q_inv);

    // _E_pj = _E_pj_0 * _Q
    _E_pj = _E_pj_0 * _Q;
    _E_jp = math::Einv(_E_pj);

    // dEpj_dq
    for (int i = 0;i < _ndof;i++)
        _dEpj_dq(i) = _E_pj_0 * _dQ_dq(i);

    // _A_jp
    _A_jp = math::Ad(_E_jp);

    // _E_0j
    if (_parent == nullptr) {
        _E_0j = _E_pj;
    } else {
        _E_0j = _parent->_E_0j * _E_pj;
    }

    // _phi_j
    if (_ndof > 0) {
        _phi = _S_j * _qdot;
    } else {
        _phi.setZero();
    }
    if (_parent != nullptr) {
        _phi += _A_jp * _parent->_phi;
    }

    if (design_gradient) {
        update_design_derivatives();
    }

    // update its body
    if (_body != nullptr) {
        _body->update(design_gradient);
    }

    // update its children
    for (auto child : _children) {
        child->update(design_gradient);
    }
}

// update design derivative related jacobians
void Joint::update_design_derivatives() {
    // dAjp0_dp1
    if (_design_params_1._active) {
        Matrix3 RT = _E_pj_0.topLeftCorner(3, 3).transpose();
        Vector3 p = _E_pj_0.topRightCorner(3, 1);
        Matrix3 tmp = math::skew(-RT * p);
        for (int k = 0;k < 3;k++)
            for (int l = 0;l < 3;l++) {
                int idx = k * 3 + l;
                _dAjp0_dp1(idx)(l, k) = 1.;
                _dAjp0_dp1(idx)(l + 3, k + 3) = 1.;
                _dAjp0_dp1(idx).bottomLeftCorner(3, 3) = math::skew(-Vector3::Unit(l) * p(k)) * RT + (tmp * Vector3::Unit(l)) * RowVector3::Unit(k);
            }
        for (int k = 0;k < 3;k++) {
            Vector3 ek = Vector3::Unit(k);
            _dAjp0_dp1(9 + k).bottomLeftCorner(3, 3) = math::skew(-RT * ek) * RT;
        }
    }

    // dE0j_dp1
    if (_sim->_robot->_ndof_p1 > 0) {
        if (_design_params_1._active) {
            // the type 1 design parameters of the current joint.
            // E_0p * dEpj0_dpk * Q
            for (int k = 0;k < 3;k++)
                for (int l = 0;l < 4;l++) {
                    int idx;
                    if (l < 3)
                        idx = _design_params_1._param_index(k * 3 + l);
                    else
                        idx = _design_params_1._param_index(9 + k);
                    Matrix4 dEpj0_dp1 = Matrix4::Zero();
                    dEpj0_dp1(k, l) = 1.;
                    if (_parent == nullptr) {
                        _dE0j_dp1(idx) = dEpj0_dp1 * _Q;
                    } else {
                        _dE0j_dp1(idx) = _parent->_E_0j * dEpj0_dp1 * _Q;
                    }
                }
        }
        // the type 1 design parameters of ancestor's joints
        for (auto ancestor = _parent; ancestor != nullptr; ancestor = ancestor->_parent) {
            if (ancestor->_design_params_1._active) {
                for (int k = 0;k < ancestor->_design_params_1._ndof;k++) {
                    int idx = ancestor->_design_params_1._param_index(k);
                    _dE0j_dp1(idx) = _parent->_dE0j_dp1(idx) * _E_pj_0 * _Q;
                }
            }
        }
    }

    // dphi_dp1
    if (_sim->_robot->_ndof_p1 > 0) {
        if (_parent != nullptr) {
            if (_design_params_1._active) {
                // the type 1 design parameters of the current joint.
                for (int k = 0;k < _design_params_1._ndof;k++) {
                    int idx = _design_params_1._param_index(k);
                    _dphi_dp1.col(idx) = _A_inv * (_dAjp0_dp1(k) * _parent->_phi);
                }
            }
            // the type 1 design parameters of ancestor's joints
            for (auto ancestor = _parent; ancestor != nullptr; ancestor = ancestor->_parent) {
                if (ancestor->_design_params_1._active) {
                    for (int k = 0;k < ancestor->_design_params_1._ndof;k++) {
                        int idx = ancestor->_design_params_1._param_index(k);
                        _dphi_dp1.col(idx) = _A_jp * _parent->_dphi_dp1.col(idx);
                    }
                }
            }
        }
    }
}

// activate design parameters
void Joint::activate_design_parameters_type_1(bool active) {
    _design_params_1.activate(active);
}

void Joint::activate_design_parameters_type_5(bool active) {
    _design_params_5.activate(active);
}

void Joint::computeJointForce(VectorX& fr) {
    if (_ndof > 0) {
        // stiffness and damping force
        fr.segment(_index[0], _ndof) += - _Kr * (_q - _q_rest) - _Dr * _qdot;
        // joint limit force
        for (int i = 0;i < _ndof;++i) {
            if (_q[i] < _joint_limit_lower)
                fr(_index[i]) += _joint_limit_stiffness * (_joint_limit_lower - _q[i]);
            if (_q[i] > _joint_limit_upper)
                fr(_index[i]) += _joint_limit_stiffness * (_joint_limit_upper - _q[i]);
        }
    }
}

void Joint::computeJointForceWithDerivative(VectorX& fr, MatrixX& Kr, MatrixX& Dr) {
    if (_ndof > 0) {
        // stiffness and damping force
        fr.segment(_index[0], _ndof) += - _Kr * (_q - _q_rest) - _Dr * _qdot;
        // joint limit force
        for (int i = 0;i < _ndof;++i) {
            if (_q[i] < _joint_limit_lower)
                fr(_index[i]) += _joint_limit_stiffness * (_joint_limit_lower - _q[i]);
            if (_q[i] > _joint_limit_upper)
                fr(_index[i]) += _joint_limit_stiffness * (_joint_limit_upper - _q[i]);
        }

        /********************************* Derivatives ************************************/
        // stiffness and damping force
        for (int k = 0;k < _ndof;k++) {
            Kr(_index[k], _index[k]) -= _Kr;
            Dr(_index[k], _index[k]) -= _Dr;
        }
        // joint limit force
        for (int i = 0;i < _ndof;++i) {
            if (_q[i] < _joint_limit_lower || _q[i] > _joint_limit_upper)
                Kr(_index[i], _index[i]) += -_joint_limit_stiffness;
        }
    }
}

void Joint::computeJointForceWithDerivative(VectorX& fr, MatrixX& Kr, MatrixX& Dr, MatrixX& dfr_dp) {
    computeJointForceWithDerivative(fr, Kr, Dr);
}

Vector3 Joint::position_in_world(Vector3& pos) {
    Matrix3 R = _E_0j.topLeftCorner(3, 3);
    Vector3 p = _E_0j.topRightCorner(3, 1);
    return R * pos + p;
}

Vector3 Joint::velocity_in_world(Vector3& pos) {
    Matrix3 R = _E_0j.topLeftCorner(3, 3);
    Vector3 vel = R * (math::skew(pos).transpose() * _phi.head(3) + _phi.tail(3));
    return vel;
}

void Joint::test_derivatives() {
    if (_ndof > 0) {
        // test _dSj_dq, _dSjdot_dq, _dA_dq, _dAdot_dq, _dQ_dq
        // srand(time(0));
        srand(1000);
        VectorX q = VectorX::Random(_ndof) * 1.;
        VectorX qdot = VectorX::Random(_ndof) * 1.;
        
        _q = q;
        _qdot = qdot;

        update();
        auto Q = _Q;
        auto S_j = _S_j;
        auto S_j_dot = _S_j_dot;
        auto A = _A;
        auto Adot = _Adot;
        auto dQ_dq = _dQ_dq;
        auto dSj_dq = _dSj_dq;
        auto dSjdot_dq = _dSjdot_dq;
        auto dA_dq = _dA_dq;
        auto dAdot_dq = _dAdot_dq;
        std::string joint_str = "Joint " + std::to_string(_id);
        dtype h = 1e-7;
        JacobianMatrixVector dQ_dq_fd = JacobianMatrixVector(4, 4, _ndof);
        JacobianMatrixVector dSj_dq_fd = JacobianMatrixVector(6, _ndof, _ndof);
        JacobianMatrixVector dSjdot_dq_fd = JacobianMatrixVector(6, _ndof, _ndof);
        JacobianMatrixVector dA_dq_fd = JacobianMatrixVector(6, 6, _ndof);
        JacobianMatrixVector dAdot_dq_fd = JacobianMatrixVector(6, 6, _ndof);
        for (int k = 0;k < _ndof;k++) {
            _q(k) += h;
            update();
            auto Q_pos = _Q;
            auto S_j_pos = _S_j;
            auto S_j_dot_pos = _S_j_dot;
            auto A_pos = _A;
            auto Adot_pos = _Adot;
            dQ_dq_fd(k) = (Q_pos - Q) / h;
            dSj_dq_fd(k) = (S_j_pos - S_j) / h;
            dSjdot_dq_fd(k) = (S_j_dot_pos - S_j_dot) / h;
            dA_dq_fd(k) = (A_pos - A) / h;
            dAdot_dq_fd(k) = (Adot_pos - Adot) / h;
            _q(k) -= h;
        }
        print_error(joint_str + ": dQ_dq", dQ_dq, dQ_dq_fd);
        print_error(joint_str + ": dSj_dq", dSj_dq, dSj_dq_fd);
        print_error(joint_str + ": dSjdot_dq", dSjdot_dq, dSjdot_dq_fd);
        print_error(joint_str + ": dA_dq", dA_dq, dA_dq_fd);
        print_error(joint_str + ": dAdot_dq", dAdot_dq, dAdot_dq_fd);
        
        Matrix6 Adot_fd = Matrix6::Zero();
        MatrixX Sdot_fd = _S_j; Sdot_fd.setZero();
        {
            _q += _qdot * h;
            update();
            auto A_pos = _A;
            auto S_pos = _S_j;
            Adot_fd = (A_pos - A) / h;
            Sdot_fd = (S_pos - S_j) / h;
        }
        print_error(joint_str + ": Adot", Adot, Adot_fd);
        print_error(joint_str + ": Sdot", S_j_dot, Sdot_fd);
    }
}

void Joint::test_derivatives_runtime() {
    if (_ndof > 0) {
        // test _dSj_dq, _dSjdot_dq, _dA_dq, _dAdot_dq
        update();
        auto S_j = _S_j;
        auto S_j_dot = _S_j_dot;
        auto A = _A;
        auto Adot = _Adot;
        auto dSj_dq = _dSj_dq;
        auto dSjdot_dq = _dSjdot_dq;
        auto dA_dq = _dA_dq;
        auto dAdot_dq = _dAdot_dq;
        auto q = _q;
        auto qdot = _qdot;
        std::string joint_str = "Joint " + std::to_string(_id);
        dtype h = 1e-8;
        JacobianMatrixVector dSj_dq_fd = JacobianMatrixVector(6, _ndof, _ndof);
        JacobianMatrixVector dSjdot_dq_fd = JacobianMatrixVector(6, _ndof, _ndof);
        JacobianMatrixVector dA_dq_fd = JacobianMatrixVector(6, 6, _ndof);
        JacobianMatrixVector dAdot_dq_fd = JacobianMatrixVector(6, 6, _ndof);
        for (int k = 0;k < _ndof;k++) {
            _q(k) += h;
            update();
            auto S_j_pos = _S_j;
            auto S_j_dot_pos = _S_j_dot;
            auto A_pos = _A;
            auto Adot_pos = _Adot;
            dSj_dq_fd(k) = (S_j_pos - S_j) / h;
            dSjdot_dq_fd(k) = (S_j_dot_pos - S_j_dot) / h;
            dA_dq_fd(k) = (A_pos - A) / h;
            dAdot_dq_fd(k) = (Adot_pos - Adot) / h;
            _q(k) -= h;
        }
        // if (_id == 0) {
        //     // std::cerr << "====================== Analytical ======================" << std::endl;
        //     // std::cerr << "S_j = " << S_j.transpose() << std::endl;
        //     // std::cerr << "S_j_dot = " << S_j_dot.transpose() << std::endl;
        //     // std::cerr << "A = " << std::endl << A << std::endl;
        //     // std::cerr << "Adot = " << std::endl <<  Adot << std::endl;
        //     // std::cerr << "dSj_dq = " << std::endl << dSj_dq << std::endl;
        //     // std::cerr << "dSjdot_dq = " << std::endl << dSjdot_dq << std::endl;
        //     std::cerr << "dA_dq = " << std::endl << dA_dq << std::endl;
        //     // std::cerr << "dAdot_dq = " << std::endl << dAdot_dq << std::endl;
            
        //     // std::cerr << "====================== Finite Difference ======================" << std::endl;
        //     // std::cerr << "dSj_dq_fd = " << std::endl << dSj_dq_fd << std::endl;
        //     // std::cerr << "dSjdot_dq_fd = " << std::endl << dSjdot_dq_fd << std::endl;
        //     std::cerr << "dA_dq_fd = " << std::endl << dA_dq_fd << std::endl;
        //     // std::cerr << "dAdot_dq_fd = " << std::endl << dAdot_dq_fd << std::endl;
        // }

        // std::cerr << "====================== Error ======================" << std::endl;
        print_error(joint_str + ": dSj_dq", dSj_dq, dSj_dq_fd);
        print_error(joint_str + ": dSjdot_dq", dSjdot_dq, dSjdot_dq_fd);
        print_error(joint_str + ": dA_dq", dA_dq, dA_dq_fd);
        print_error(joint_str + ": dAdot_dq", dAdot_dq, dAdot_dq_fd);

        // restore
        _q = q;
        _qdot = qdot;
        update();
    }

}

}