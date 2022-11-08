#include "Joint/JointFree3DEuler.h"
#include "Utils.h"

namespace redmax {

JointFree3DEuler::JointFree3DEuler(Simulation *sim, int id, Joint* parent, Matrix3 R_pj_0, Vector3 p_pj_0, 
    JointSphericalEuler::Chart chart, Joint::Frame frame) 
    : Joint(sim, id, parent, R_pj_0, p_pj_0, 6, frame) {
    
    _joint_translational = new JointTranslational(sim, -1, nullptr, R_pj_0, p_pj_0, frame);
    _joint_spherical_euler = new JointSphericalEuler(sim, -1, nullptr, R_pj_0, p_pj_0, chart, frame);
}

void JointFree3DEuler::update(bool design_gradient) {
    Vector3 p = _q.head(3);
    Vector3 r = _q.tail(3);
    Vector3 pdot = _qdot.head(3);
    Vector3 rdot = _qdot.tail(3);
    
    _joint_translational->_q = p;
    _joint_translational->_qdot = pdot;
    _joint_translational->update();

    _joint_spherical_euler->_q = r;
    _joint_spherical_euler->_qdot = rdot;
    _joint_spherical_euler->update();

    Matrix3 R = _joint_spherical_euler->_Q.topLeftCorner(3, 3);
    Matrix3 Rdot = _joint_spherical_euler->_Adot.topLeftCorner(3, 3);
    Matrix3 T = _joint_spherical_euler->_S_j.topRows(3);
    Matrix3 Tdot = _joint_spherical_euler->_S_j_dot.topRows(3);

    JacobianMatrixVector dR_dr(3, 3, 3);
    JacobianMatrixVector dRdot_dr(3, 3, 3);
    JacobianMatrixVector dT_dr(3, 3, 3);
    JacobianMatrixVector dTdot_dr(3, 3, 3);
    for (int i = 0;i < 3;i++) {
        dR_dr(i) = _joint_spherical_euler->_dA_dq(i).topLeftCorner(3, 3);
        dRdot_dr(i) = _joint_spherical_euler->_dAdot_dq(i).topLeftCorner(3, 3);
        dT_dr(i) = _joint_spherical_euler->_dSj_dq(i).topRows(3);
        dTdot_dr(i) = _joint_spherical_euler->_dSjdot_dq(i).topRows(3);
    }

    _Q.setIdentity();
    _Q.topLeftCorner(3, 3) = R;
    _Q.topRightCorner(3, 1) = p;

    _dQ_dq.setZero();
    for (int i = 0;i < 3;i++) {
        _dQ_dq(i)(i, 3) = 1.;
    }
    for (int i = 3;i < 6;i++) {
        _dQ_dq(i).topLeftCorner(3, 3) = dR_dr(i - 3);
    }

    _A = math::Ad(_Q);
    
    _dA_dq.setZero();
    for (int i = 0;i < 3;i++) {
        _dA_dq(i).bottomLeftCorner(3, 3) = math::skew(Vector3::Unit(i)) * R;
    }
    for (int i = 3;i < 6;i++) {
        _dA_dq(i).topLeftCorner(3, 3) = dR_dr(i - 3);
        _dA_dq(i).bottomLeftCorner(3, 3) = math::skew(p) * dR_dr(i - 3);
        _dA_dq(i).bottomRightCorner(3, 3) = dR_dr(i - 3);
    }

    _Adot.setZero();
    _Adot.topLeftCorner(3, 3) = Rdot; 
    _Adot.bottomLeftCorner(3, 3) = math::skew(pdot) * R + math::skew(p) * Rdot;
    _Adot.bottomRightCorner(3, 3) = Rdot;

    _dAdot_dq.setZero();
    for (int i = 0;i < 3;i++) {
        _dAdot_dq(i).bottomLeftCorner(3, 3) = math::skew(Vector3::Unit(i)) * Rdot;
    }
    for (int i = 3;i < 6;i++) {
        _dAdot_dq(i).topLeftCorner(3, 3) = dRdot_dr(i - 3);
        _dAdot_dq(i).bottomLeftCorner(3, 3) = math::skew(pdot) * dR_dr(i - 3) + math::skew(p) * dRdot_dr(i - 3);
        _dAdot_dq(i).bottomRightCorner(3, 3) = dRdot_dr(i - 3);
    }

    _S_j.setZero();
    _S_j.topRightCorner(3, 3) = T;
    _S_j.bottomLeftCorner(3, 3) = R.transpose();

    _dSj_dq.setZero();
    for (int i = 3;i < 6;i++) {
        _dSj_dq(i).topRightCorner(3, 3) = dT_dr(i - 3);
        _dSj_dq(i).bottomLeftCorner(3, 3) = dR_dr(i - 3).transpose();
    }

    _S_j_dot.setZero();
    _S_j_dot.topRightCorner(3, 3) = Tdot;
    _S_j_dot.bottomLeftCorner(3, 3) = Rdot.transpose();

    _dSjdot_dq.setZero();
    for (int i = 3;i < 6;i++) {
        _dSjdot_dq(i).topRightCorner(3, 3) = dTdot_dr(i - 3);
        _dSjdot_dq(i).bottomLeftCorner(3, 3) = dRdot_dr(i - 3).transpose();
    }

    Joint::update(design_gradient);
}

}