#include "Joint/JointDesignParameters.h"
#include "Joint/Joint.h"
#include "Utils.h"

namespace redmax {

JointDesignParameters::JointDesignParameters(
    Joint* joint, DesignParameterType type, bool active, int ndof)
    : DesignParameters(type, active, ndof) {
    
    if (_type != DesignParameterType::TYPE1 && _type != DesignParameterType::TYPE5) {
        throw_error("[Error] Joint can only have type I, V design parameters");
    }

    _joint = joint;
}

void JointDesignParameters::update_params() {
    if (_type == DesignParameterType::TYPE1) {
        _params = math::flatten_E(_joint->_E_pj_0);
    } else if (_type == DesignParameterType::TYPE5) {
        _params[0] = _joint->_Kr;
        _params[1] = _joint->_Dr;
    } else {
        throw_error("[Error] Joint can only have type I, V design parameters");
    }
}

void JointDesignParameters::update_params(const VectorX params) {
    if (!_active)
        return;
        
    for (int i = 0;i < _ndof;i++) {
        _params[i] = params[_param_index[i]];
    }

    if (_type == DesignParameterType::TYPE1) {
        _joint->_E_pj_0 = math::compose_E(_params);
        _joint->_E_jp_0 = math::Einv(_joint->_E_pj_0);
    } else if (_type == DesignParameterType::TYPE5) {
        _joint->_Kr = _params[0];
        _joint->_Dr = _params[1];
    } else {
        throw_error("[Error] Joint can only have type I, V design parameters");
    }
}

}