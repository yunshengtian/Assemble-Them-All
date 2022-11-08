#pragma once
#include "Common.h"
#include "Joint/Joint.h"

namespace redmax {

// JointPlanar is a planar joint
class JointPlanar : public Joint {
public:
    Vector3 _axis0, _axis1;

    JointPlanar(Simulation *sim, int id, Vector3 axis0, Vector3 axis1, Joint *parent, Matrix3 R_pj_0, Vector3 p_pj_0,
        Joint::Frame frame = Joint::Frame::LOCAL) : Joint(sim, id, parent, R_pj_0, p_pj_0, 2, frame) {
        _axis0 = axis0 / axis0.norm();
        _axis1 = axis1 / axis1.norm();
        if (frame == Joint::Frame::WORLD) {
            _axis0 = _E_j0_0.topLeftCorner(3, 3) * _axis0;
            _axis1 = _E_j0_0.topLeftCorner(3, 3) * _axis1;
        }

        _S_j.setZero();
        _S_j.block(3, 0, 3, 1) = _axis0;
        _S_j.block(3, 1, 3, 1) = _axis1;
    }

    virtual void update(bool design_gradient = false);
};

}