#pragma once
#include "DesignParameters.h"

namespace redmax {

class Body;

class BodyDesignParameters : public DesignParameters {
public:
    BodyDesignParameters();
    BodyDesignParameters(Body* body, DesignParameterType type, bool active, int ndof = 0);

    Body* _body;

    void update_params();
    void update_params(const VectorX params);
};

}