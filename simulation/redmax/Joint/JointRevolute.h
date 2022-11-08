#pragma once
#include "Joint/Joint.h"
#include "Common.h"

namespace redmax {

class JointRevolute : public Joint {
public:
    Vector3 _axis;

    JointRevolute(Simulation *sim, int id, Vector3 axis, Joint* parent, Matrix3 R_pj_0, Vector3 p_pj_0, Joint::Frame frame = Joint::Frame::LOCAL);
    
    virtual void update(bool design_gradient = false);
};

}