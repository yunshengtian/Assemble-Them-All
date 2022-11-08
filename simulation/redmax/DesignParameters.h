#pragma once

#include "Common.h"

namespace redmax {

enum DesignParameterType {
    TYPE1 = 1,  // joint transformation (kinematic)
    TYPE2 = 2,  // body transformation (kinematic)
    TYPE3 = 3,  // body contact points
    TYPE4 = 4,  // body mass properties
    TYPE5 = 5,  // joint coefficients
    TYPE6 = 6,  // body contact coefficients
};

class DesignParameters {
public:
    DesignParameters(DesignParameterType type, bool active, int ndof = 0);

    DesignParameterType _type;
    bool _active;
    int _ndof;
    VectorX _params;
    VectorXi _param_index;

    void activate(bool active = true, int ndof = 0);

    virtual void update_params() = 0;
    virtual void update_params(const VectorX params) = 0;
    void get_params(VectorX &all_params);
    VectorX get_params();
};

}