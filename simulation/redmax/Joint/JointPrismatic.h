#pragma once
#include "Joint/Joint.h"
#include "Common.h"

namespace redmax {

class JointPrismatic : public Joint {
public:
    Vector3 _axis;

    JointPrismatic(Simulation *sim, int id, Vector3 axis, Joint* parent, Matrix3 R_pj_0, Vector3 p_pj_0, Joint::Frame frame = Joint::Frame::LOCAL) 
        : Joint(sim, id, parent, R_pj_0, p_pj_0, 1, frame) {
        _axis = axis / axis.norm();
        if (frame == Joint::Frame::WORLD) {
            _axis = _E_j0_0.topLeftCorner(3, 3) * _axis;
        }

        _S_j.setZero();
        _S_j.block(3, 0, 3, 1) = _axis;
    }
    
    virtual void update(bool design_gradient = false);
};

}