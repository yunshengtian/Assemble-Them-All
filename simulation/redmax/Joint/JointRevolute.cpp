#include "Joint/JointRevolute.h"
#include "Utils.h"

namespace redmax {

JointRevolute::JointRevolute(Simulation *sim, int id, Vector3 axis, Joint* parent, Matrix3 R_pj_0, Vector3 p_pj_0, Joint::Frame frame) 
    : Joint(sim, id, parent, R_pj_0, p_pj_0, 1, frame) {
    _axis = axis / axis.norm();
    if (frame == Joint::Frame::WORLD) {
        _axis = _E_j0_0.topLeftCorner(3, 3) * _axis;
    }

    _S_j.setZero();
    _S_j.block(0, 0, 3, 1) = _axis;
}

// update _Q_pj matrix
/*  notes section 3.2.1
    R(q) = exp([axis * q])
    Q(q) = [R 0]
           [0 1]
    A(q) = [R 0]
           [0 R]
    dR_dq = R * [axis]
    dA_dq = [dR_dq     0]
            [0     dR_dq]
    Rdot = dR_dq * qdot
    Adot = [Rdot    0]
           [0    Rdot]
    dRdot_dq = R * [axis] * [axis] * qdot
    dAdot_dq = [dRdot_dq        0]
               [0        dRdot_dq]
    dS_dq = 0
    dSdot_dq = 0
    dQ_dq = [dR_dq 0]
            [0     0]
*/
void JointRevolute::update(bool design_gradient) {
    Matrix3 R = Eigen::AngleAxis<dtype>(_q(0), _axis).matrix();
    _Q = math::SE(R, Vector3::Zero());
    
    _A = math::Ad(_Q);

    Matrix3 axis_brac = math::skew(_axis);
    Matrix3 dR_dq = R * axis_brac;
    
    _dA_dq(0).block(0, 0, 3, 3) = dR_dq;
    _dA_dq(0).block(3, 3, 3, 3) = dR_dq;

    Matrix3 Rdot = dR_dq * _qdot(0);

    _Adot.setZero();
    _Adot.block(0, 0, 3, 3) = Rdot;
    _Adot.block(3, 3, 3, 3) = Rdot;

    Matrix3 dRdot_dq = dR_dq * axis_brac * _qdot(0);

    _dAdot_dq.setZero();
    _dAdot_dq(0).block(0, 0, 3, 3) = dRdot_dq;
    _dAdot_dq(0).block(3, 3, 3, 3) = dRdot_dq;
    
    _dSj_dq.setZero();
    _dSjdot_dq.setZero();

    _dQ_dq.setZero();
    _dQ_dq(0).topLeftCorner(3, 3) = dR_dq;

    Joint::update(design_gradient);
}

}