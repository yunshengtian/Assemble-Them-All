#include "Joint/JointFixed.h"

namespace redmax {

// update _Q
void JointFixed::update(bool design_gradient) {
    _Q.setIdentity();
    _A = math::Ad(_Q);
    _dA_dq.setZero();
    _Adot.setZero();
    
    Joint::update(design_gradient);
}

}