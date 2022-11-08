#pragma once
#include "Common.h"
#include "Joint/Joint.h"

namespace redmax {

// JointTranslational is a 3D translational joint
class JointTranslational : public Joint {
public:
    JointTranslational(Simulation *sim, int id, Joint *parent, Matrix3 R_pj_0, Vector3 p_pj_0,
        Joint::Frame frame = Joint::Frame::LOCAL);

    virtual void update(bool design_gradient = false);
};

}