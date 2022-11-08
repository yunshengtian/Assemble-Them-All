#pragma once
#include "DesignParameters.h"

namespace redmax {

class Joint;

class JointDesignParameters : public DesignParameters {
public:
    JointDesignParameters(Joint* joint, DesignParameterType type, bool active, int ndof = 0);

    Joint* _joint;

    void update_params();
    void update_params(const VectorX params);
};

}