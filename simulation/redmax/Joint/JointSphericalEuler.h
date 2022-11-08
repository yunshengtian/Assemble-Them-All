#pragma once
#include "Common.h"
#include "Joint/Joint.h"

namespace redmax {

// JointSphericalEuler is a spherical joint represented by euler angles
class JointSphericalEuler : public Joint {
public:
    enum Chart {
        XYX, XZX, YZY, YXY, ZXZ, ZYZ,
        XYZ, XZY, YZX, YXZ, ZXY, ZYX
    };

    JointSphericalEuler(Simulation *sim, int id, Joint * parent, Matrix3 R_pj_0, Vector3 p_pj_0,
        Chart chart = Chart::XYZ, Joint::Frame frame = Joint::Frame::LOCAL);

    virtual void update(bool design_gradient = false);

    Chart _chart;
};

}