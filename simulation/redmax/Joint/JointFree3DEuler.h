#pragma once
#include "Common.h"
#include "Joint/Joint.h"
#include "Joint/JointSphericalEuler.h"
#include "Joint/JointTranslational.h"

namespace redmax {

// JointFree3D is a 3D free joint represented by euler angles
class JointFree3DEuler : public Joint {
public:
    JointFree3DEuler(Simulation *sim, int id, Joint * parent, Matrix3 R_pj_0, Vector3 p_pj_0,
        JointSphericalEuler::Chart chart = JointSphericalEuler::Chart::XYZ, 
        Joint::Frame frame = Joint::Frame::LOCAL);

    virtual void update(bool design_gradient = false);

    JointTranslational* _joint_translational;
    JointSphericalEuler* _joint_spherical_euler;
};

}