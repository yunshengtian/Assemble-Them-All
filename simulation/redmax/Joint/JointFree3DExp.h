#pragma once
#include "Common.h"
#include "Joint/Joint.h"
#include "Joint/JointSphericalExp.h"
#include "Joint/JointTranslational.h"

namespace redmax {

// JointFree3D is a 3D free joint represented by euler angles
class JointFree3DExp : public Joint {
public:
    JointFree3DExp(Simulation *sim, int id, Joint * parent, Matrix3 R_pj_0, Vector3 p_pj_0,
        Joint::Frame frame = Joint::Frame::LOCAL);

    virtual void update(bool design_gradient = false);

    bool reparam();
    
    JointTranslational* _joint_translational;
    JointSphericalExp* _joint_spherical_exp;
};

}