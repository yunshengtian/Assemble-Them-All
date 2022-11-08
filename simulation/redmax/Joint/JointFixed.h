#pragma once
#include "Common.h"
#include "Joint/Joint.h"

namespace redmax {

class JointFixed : public Joint {
public:
    JointFixed(Simulation *sim, int id, Joint* parent, Matrix3 R_pj_0, Vector3 p_pj_0, Joint::Frame frame = Joint::Frame::LOCAL) 
        : Joint(sim, id, parent, R_pj_0, p_pj_0, 0, frame) {}

    virtual void update(bool design_gradient = false);
};

}