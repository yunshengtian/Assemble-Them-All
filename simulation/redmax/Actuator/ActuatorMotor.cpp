#include "Actuator/ActuatorMotor.h"
#include "Joint/Joint.h"

namespace redmax {

ActuatorMotor::ActuatorMotor(std::string name, Joint* joint, ControlMode control_mode, VectorX ctrl_min, VectorX ctrl_max, VectorX ctrl_P, VectorX ctrl_D)
    : Actuator(name, joint->_ndof, ctrl_min, ctrl_max) {
    
    _joint = joint;
    _control_mode = control_mode;
    _ctrl_P = ctrl_P;
    _ctrl_D = ctrl_D;
}

ActuatorMotor::ActuatorMotor(std::string name, Joint* joint, ControlMode control_mode, dtype ctrl_min, dtype ctrl_max, dtype ctrl_P, dtype ctrl_D)
    : Actuator(name, joint->_ndof, ctrl_min, ctrl_max) {
    
    _joint = joint;
    _control_mode = control_mode;
    _ctrl_P = VectorX::Constant(_ndof, ctrl_P);
    _ctrl_D = VectorX::Constant(_ndof, ctrl_D);
}

void ActuatorMotor::update_dofs(const VectorX& dofs, const VectorX& dofs_vel) {
    for (int i = 0;i < _ndof;i++) {
        _dofs[i] = dofs[_joint->_index[i]];
        _dofs_vel[i] = dofs_vel[_joint->_index[i]];
    }
}

void ActuatorMotor::computeForce(VectorX& fm, VectorX& fr) {
    if (_control_mode == ControlMode::FORCE) {
        VectorX u = _u.cwiseMin(VectorX::Ones(_joint->_ndof)).cwiseMax(VectorX::Ones(_joint->_ndof) * -1.);
        fr.segment(_joint->_index[0], _joint->_ndof) += map_value(u, VectorX::Constant(_joint->_ndof, -1), VectorX::Ones(_joint->_ndof, 1), _ctrl_min, _ctrl_max);
        // fr.segment(_joint->_index[0], _joint->_ndof) += ((u + VectorX::Ones(_joint->_ndof)) / 2.).cwiseProduct(_ctrl_max - _ctrl_min) + _ctrl_min;
    } else {
        _pos_error = _u - _dofs;
        _vel_error = - _dofs_vel;
        fr.segment(_joint->_index[0], _joint->_ndof) += (_ctrl_P.cwiseProduct(_pos_error) + _ctrl_D.cwiseProduct(_vel_error)).cwiseMin(_ctrl_max).cwiseMax(_ctrl_min);
    }
}

void ActuatorMotor::computeForceWithDerivative(VectorX& fm, VectorX& fr, MatrixX& Km, MatrixX& Dm, MatrixX& Kr, MatrixX& Dr) {
    // VectorX u = _u.cwiseMin(VectorX::Ones(_joint->_ndof)).cwiseMax(VectorX::Ones(_joint->_ndof) * -1.);
    // fr.segment(_joint->_index[0], _joint->_ndof) += ((u + VectorX::Ones(_joint->_ndof)) / 2.).cwiseProduct(_ctrl_max - _ctrl_min) + _ctrl_min;
    computeForce(fm, fr);
}

void ActuatorMotor::compute_dfdu(MatrixX& dfm_du, MatrixX& dfr_du) {
    if (_control_mode == ControlMode::FORCE) {
        for (int i = 0;i < _joint->_ndof;i++)
            if (_u[i] >= -1. && _u[i] <= 1.) {
                dfr_du(_joint->_index[i], _index[i]) += (_ctrl_max[i] - _ctrl_min[i]) / 2.;
            }
    } else { // TODO: to be tested
        for (int i = 0;i < _joint->_ndof;i++) {
            dtype f = _ctrl_P[i] * _pos_error[i] + _ctrl_D[i] * _vel_error[i];
            if (f >= _ctrl_min[i] && f <= _ctrl_max[i])
                dfr_du(_joint->_index[i], _index[i]) += _ctrl_P[i];
        }
    }
}

}