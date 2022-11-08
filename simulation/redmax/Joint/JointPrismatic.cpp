#include "Joint/JointPrismatic.h"
#include "Utils.h"

namespace redmax {

// update _Q_pj matrix
/*  notes section 3.2.1
    Q(q) = [I a*q]
           [0   I]
    A(q) = [I     0]
           [[a]*q I]
    dA_dq = [0    0]
            [[a]  0]
    Adot = [0         0]
           [[a]*qdot  0]
    dAdot_dq = 0
    dS_dq = 0
    dSdot_dq = 0
    dQ_dq = [0 a]
            [0 0]
*/
void JointPrismatic::update(bool design_gradient) {
    _Q.setIdentity();
    _Q.topRightCorner(3, 1) = _axis * _q(0);
    
    _A = math::Ad(_Q);

    Matrix3 axis_brac = math::skew(_axis);

    _dA_dq(0).setZero();
    _dA_dq(0).bottomLeftCorner(3, 3).noalias() = axis_brac;

    _Adot.setZero();
    _Adot.bottomLeftCorner(3, 3).noalias() = axis_brac * _qdot(0);
    
    _dAdot_dq.setZero();
    
    _dSj_dq.setZero();
    _dSjdot_dq.setZero();

    _dQ_dq.setZero();
    _dQ_dq(0).topRightCorner(3, 1) = _axis;
    
    Joint::update(design_gradient);
}

}