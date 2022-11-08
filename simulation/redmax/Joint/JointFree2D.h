#pragma once
#include "Common.h"
#include "Joint/Joint.h"

namespace redmax {

// JointFree2D is a planar joint on X-Y plane and a revolute joint around Z axis.
class JointFree2D : public Joint {
public:
    JointFree2D(Simulation *sim, int id, Joint* parent, Matrix3 R_pj_0, Vector3 p_pj_0, Joint::Frame frame = Joint::Frame::LOCAL);

    virtual void update(bool design_gradient = false);
};

}
