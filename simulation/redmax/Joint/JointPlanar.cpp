#include "Joint/JointPlanar.h"
#include "Utils.h"

namespace redmax {

// update _Q_pj matrix
void JointPlanar::update(bool design_gradient) {
    _Q.setIdentity();
    _Q.topRightCorner(3, 1) = _axis0 * _q(0) + _axis1 * _q(1);
    
    _A = math::Ad(_Q);

    Matrix3 axis0_brac = math::skew(_axis0); 
    Matrix3 axis1_brac = math::skew(_axis1); 

    _dA_dq.setZero();
    _dA_dq(0).bottomLeftCorner(3, 3).noalias() = axis0_brac;
    _dA_dq(1).bottomLeftCorner(3, 3).noalias() = axis1_brac;

    _Adot.setZero();
    _Adot.bottomLeftCorner(3, 3).noalias() = axis0_brac * _qdot(0) + axis1_brac * _qdot(1);
    
    _dAdot_dq.setZero();
    
    _dSj_dq.setZero();
    _dSjdot_dq.setZero();

    _dQ_dq.setZero();
    _dQ_dq(0).topRightCorner(3, 1) = _axis0;
    _dQ_dq(1).topRightCorner(3, 1) = _axis1;
    
    Joint::update(design_gradient);
}

}