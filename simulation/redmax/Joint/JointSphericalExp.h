#pragma once
#include "Common.h"
#include "Joint/Joint.h"

namespace redmax {

// JointSphericalExp is a spherical joint represented by exponential coordinates
class JointSphericalExp : public Joint {
public:
    JointSphericalExp(Simulation *sim, int id, Joint * parent, Matrix3 R_pj_0, Vector3 p_pj_0,
        Joint::Frame frame = Joint::Frame::LOCAL);

    virtual void update(bool design_gradient = false);

    void inner_update();
    
    bool reparam();
};

}