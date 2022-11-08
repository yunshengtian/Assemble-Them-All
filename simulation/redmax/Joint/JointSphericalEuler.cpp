#include "Joint/JointSphericalEuler.h"
#include "Utils.h"

namespace redmax {

JointSphericalEuler::JointSphericalEuler(
    Simulation *sim, int id, Joint * parent, Matrix3 R_pj_0, Vector3 p_pj_0,
    Chart chart, Joint::Frame frame) 
    : Joint(sim, id, parent, R_pj_0, p_pj_0, 3, frame) {
    _chart = chart;
}

/* notes section 2.4.5 and 3.2.5
*/
void JointSphericalEuler::update(bool design_gradient) {
    dtype qdot1 = _qdot(0), qdot2 = _qdot(1), qdot3 = _qdot(2);
    dtype c1 = cos(_q(0)), s1 = sin(_q(0));
    dtype c2 = cos(_q(1)), s2 = sin(_q(1));
    dtype c3 = cos(_q(2)), s3 = sin(_q(2));

    Matrix3 R;
    Matrix3 dR_dq0, dR_dq1, dR_dq2;
    Matrix3 dRdot_dq0, dRdot_dq1, dRdot_dq2;

    if (_chart == Chart::XYZ) {
        R <<               c2 * c3,               -c2 * s3,     s2, 
            c1 * s3 + c3 * s1 * s2, c1 * c3 - s1 * s2 * s3, -c2 * s1,
            s1 * s3 - c1 * c3 * s2, c3 * s1 + c1 * s2 * s3, c1 * c2;
        dR_dq0 <<                      0,                       0,        0,
                  c1 * c3 * s2 - s1 * s3, -c3 * s1 - c1 * s2 * s3, -c1 * c2,
                  c1 * s3 + c3 * s1 * s2,  c1 * c3 - s1 * s2 * s3, -c2 * s1;
        dR_dq1 <<      -c3 * s2,   s2 * s3, c2,
                   c2 * c3 * s1, -c2 * s1 * s3,  s1 * s2,
                  -c1 * c2 * c3,  c1 * c2 * s3, -c1 * s2;
        dR_dq2 <<               -c2 * s3,                -c2 * c3, 0,
                  c1 * c3 - s1 * s2 * s3, -c1 * s3 - c3 * s1 * s2, 0,
                  c3 * s1 + c1 * s2 * s3,  c1 * c3 * s2 - s1 * s3, 0;
        dRdot_dq0 << 0, 0, 0,
                     qdot2 * c1 * c2 * c3 - qdot3 * (c3 * s1 + c1 * s2 * s3) - qdot1 * (c1 * s3 + c3 * s1 * s2),
                     qdot3 * (s1 * s3 - c1 * c3 * s2) - qdot1 * (c1 * c3 - s1 * s2 * s3) - qdot2 * c1 * c2 * s3,
                     qdot1 * c2 * s1 + qdot2 * c1 * s2,
                     qdot3 * (c1 * c3 - s1 * s2 * s3) - qdot1 * (s1 * s3 - c1 * c3 * s2) + qdot2 * c2 * c3 * s1,
                     -qdot1 * (c3 * s1 + c1 * s2 * s3) - qdot3 * (c1 * s3 + c3 * s1 * s2) - qdot2 * c2 * s1 * s3,
                     qdot2 * s1 * s2 - qdot1 * c1 * c2;
        dRdot_dq1 << qdot3 * s2 * s3 - qdot2 * c2 * c3, qdot2 * c2 * s3 + qdot3 * c3 * s2, -qdot2 * s2,
                     qdot1 * c1 * c2 * c3 - qdot2 * c3 * s1 * s2 - qdot3 * c2 * s1 * s3,
                     qdot2 * s1 * s2 * s3 - qdot3 * c2 * c3 * s1 - qdot1 * c1 * c2 * s3,
                     qdot1 * c1 * s2 + qdot2 * c2 * s1,
                     qdot1 * c2 * c3 * s1 + qdot2 * c1 * c3 * s2 + qdot3 * c1 * c2 * s3,
                     qdot3 * c1 * c2 * c3 - qdot1 * c2 * s1 * s3 - qdot2 * c1 * s2 * s3,
                     qdot1 * s1 * s2 - qdot2 * c1 * c2;
        dRdot_dq2 << qdot2 * s2 * s3 - qdot3 * c2 * c3, qdot2 * c3 * s2 + qdot3 * c2 * s3, 0,
                     -qdot1 * (c3 * s1 + c1 * s2 * s3) - qdot3 * (c1 * s3 + c3 * s1 * s2) - qdot2 * c2 * s1 * s3,
                     qdot1 * (s1 * s3 - c1 * c3 * s2) - qdot3 * (c1 * c3 - s1 * s2 * s3) - qdot2 * c2 * c3 * s1, 0,
                     qdot1 * (c1 * c3 - s1 * s2 * s3) - qdot3 * (s1 * s3 - c1 * c3 * s2) + qdot2 * c1 * c2 * s3, 
                     qdot2 * c1 * c2 * c3 - qdot3 * (c3 * s1 + c1 * s2 * s3) - qdot1 * (c1 * s3 + c3 * s1 * s2), 0;
        _S_j.setZero();
        _S_j(0, 0) = c2 * c3; _S_j(0, 1) = s3;
        _S_j(1, 0) = -c2 * s3; _S_j(1, 1) = c3;
        _S_j(2, 0) = s2; _S_j(2, 2) = 1;

        _S_j_dot.setZero();
        _S_j_dot(0, 0) = -qdot2 * c3 * s2 - qdot3 * c2 * s3; _S_j_dot(0, 1) = qdot3 * c3;
        _S_j_dot(1, 0) = qdot2 * s2 * s3 - qdot3 * c2 * c3; _S_j_dot(1, 1) = -qdot3 * s3;
        _S_j_dot(2, 0) = qdot2 * c2;

        _dSj_dq.setZero();
        _dSj_dq(1)(0, 0) = -c3 * s2; _dSj_dq(1)(1, 0) = s2 * s3; _dSj_dq(1)(2, 0) = c2;
        _dSj_dq(2)(0, 0) = -c2 * s3; _dSj_dq(2)(0, 1) = c3; 
        _dSj_dq(2)(1, 0) = -c2 * c3; _dSj_dq(2)(1, 1) = -s3;

        _dSjdot_dq.setZero();
        _dSjdot_dq(1)(0, 0) = qdot3 * s2 * s3 - qdot2 * c2 * c3;
        _dSjdot_dq(1)(1, 0) = qdot2 * c2 * s3 + qdot3 * c3 * s2;
        _dSjdot_dq(1)(2, 0) = -qdot2 * s2;
        _dSjdot_dq(2)(0, 0) = qdot2 * s2 * s3 - qdot3 * c2 * c3;
        _dSjdot_dq(2)(0, 1) = -qdot3 * s3;
        _dSjdot_dq(2)(1, 0) = qdot2 * c3 * s2 + qdot3 * c2 * s3;
        _dSjdot_dq(2)(1, 1) = -qdot3 * c3;
    }

    if (_S_j.topRows(3).determinant() < 1e-5) {
        std::cerr << "gimbal lock" << std::endl;
    }

    _Q.setIdentity();
    _Q.topLeftCorner(3, 3) = R;

    _dQ_dq.setZero();
    _dQ_dq(0).topLeftCorner(3, 3) = dR_dq0;
    _dQ_dq(1).topLeftCorner(3, 3) = dR_dq1;
    _dQ_dq(2).topLeftCorner(3, 3) = dR_dq2;

    _A = math::Ad(_Q);
    
    _dA_dq.setZero();
    _dA_dq(0).topLeftCorner(3, 3) = dR_dq0; _dA_dq(0).bottomRightCorner(3, 3) = dR_dq0;
    _dA_dq(1).topLeftCorner(3, 3) = dR_dq1; _dA_dq(1).bottomRightCorner(3, 3) = dR_dq1;
    _dA_dq(2).topLeftCorner(3, 3) = dR_dq2; _dA_dq(2).bottomRightCorner(3, 3) = dR_dq2;

    Matrix3 Rdot = dR_dq0 * _qdot(0) + dR_dq1 * _qdot(1) + dR_dq2 * _qdot(2);
    _Adot.setZero();
    _Adot.topLeftCorner(3, 3) = Rdot; _Adot.bottomRightCorner(3, 3) = Rdot;

    _dAdot_dq.setZero();
    _dAdot_dq(0).topLeftCorner(3, 3) = dRdot_dq0; _dAdot_dq(0).bottomRightCorner(3, 3) = dRdot_dq0;
    _dAdot_dq(1).topLeftCorner(3, 3) = dRdot_dq1; _dAdot_dq(1).bottomRightCorner(3, 3) = dRdot_dq1;
    _dAdot_dq(2).topLeftCorner(3, 3) = dRdot_dq2; _dAdot_dq(2).bottomRightCorner(3, 3) = dRdot_dq2;

    Joint::update(design_gradient);
}

}