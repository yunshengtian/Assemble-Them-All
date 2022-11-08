#include "DesignParameters.h"

namespace redmax {

DesignParameters::DesignParameters(DesignParameterType type, bool active, int ndof) {
    _type = type;
    _active = active;
    _ndof = 0;
    if (_active) {
        activate(_active, ndof);
    }
}

void DesignParameters::activate(bool active, int ndof) {
    _active = active;
    if (_active) {
        if (_type == DesignParameterType::TYPE1) {
            _ndof = 12; // joint transformation (9 for rotation, 3 for translation)
        } else if (_type == DesignParameterType::TYPE2) {
            _ndof = 12; // body transformation (9 for rotation, 3 for translation)
        } else if (_type == DesignParameterType::TYPE3) {
            _ndof = ndof; // contact points ndof should be passed in
        } else if (_type == DesignParameterType::TYPE4) {
            _ndof = 4; // body mass properties (1 for mass, 3 for inertia)
        } else if (_type == DesignParameterType::TYPE5) {
            _ndof = 2; // joint coefficients (1 for stiffness, 1 for damping)
        } else if (_type == DesignParameterType::TYPE6) {
            _ndof = 1; // body contact coefficients
        } else {
            throw_error("[Error] Undefined design parameter type");
        }
        _params = VectorX::Zero(_ndof);
        _param_index = VectorXi::Zero(_ndof);

        update_params();
    } else {
        _ndof = 0;
    }
}

void DesignParameters::get_params(VectorX &all_params) {
    if (_active) {
        for (int i = 0;i < _param_index.size();i++) {
            all_params[_param_index[i]] = _params[i];
        }
    }
}

VectorX DesignParameters::get_params() {
    return _params;
}

}