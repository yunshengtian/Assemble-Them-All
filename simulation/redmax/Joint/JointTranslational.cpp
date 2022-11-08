#include "Joint/JointTranslational.h"
#include "Utils.h"

namespace redmax {

JointTranslational::JointTranslational(Simulation *sim, int id, Joint* parent, Matrix3 R_pj_0, Vector3 p_pj_0, Joint::Frame frame) 
    : Joint(sim, id, parent, R_pj_0, p_pj_0, 3, frame) {}

void JointTranslational::update(bool design_gradient) {
    _Q.setIdentity();
    _Q.topRightCorner(3, 1) = _q;
    
    _A = math::Ad(_Q);

    _S_j.setZero();
    _S_j.bottomRows(3).setIdentity();

    _Adot.setZero();
    _Adot.bottomLeftCorner(3, 3) = math::skew(_qdot);

    _S_j_dot.setZero();

    _dSj_dq.setZero();

    _dSjdot_dq.setZero();

    _dA_dq.setZero();
    for (int i = 0;i < 3;i++) {
        _dA_dq(i).bottomLeftCorner(3, 3) = math::skew(Vector3::Unit(i));
    }

    _dAdot_dq.setZero();

    _dQ_dq.setZero();
    for (int i = 0;i < 3;i++) {
        _dQ_dq(i)(i, 3) = 1;
    }

    Joint::update(design_gradient);
}

}