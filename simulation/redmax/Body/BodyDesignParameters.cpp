#include "Body/BodyDesignParameters.h"
#include "Body/Body.h"
#include "Utils.h"

namespace redmax {

BodyDesignParameters::BodyDesignParameters(
    Body* body, DesignParameterType type, bool active, int ndof)
    : DesignParameters(type, active, ndof) {
    
    if (_type != DesignParameterType::TYPE2 
        && _type != DesignParameterType::TYPE3
        && _type != DesignParameterType::TYPE4
        && _type != DesignParameterType::TYPE6) {
        throw_error("[Error] Body can only have type II, III, IV design parameters");
    }

    _body = body;
}

void BodyDesignParameters::update_params() {
    if (_type == DesignParameterType::TYPE2) {
        _params = math::flatten_E(_body->_E_ji);
    } else if (_type == DesignParameterType::TYPE3) {
        std::vector<Vector3> contacts = _body->get_contact_points();
        for (int i = 0;i < contacts.size();i++) {
            _params.segment(i * 3, 3) = contacts[i];
        }
    } else if (_type == DesignParameterType::TYPE4) {
        _params[0] = _body->_mass;
        _params.tail(3) = _body->_Inertia.head(3);
    } else if (_type == DesignParameterType::TYPE6) {
        _params[0] = _body->_contact_scale;
    } else {
        throw_error("[Error] Joint can only have type I, V design parameters");
    }
}

void BodyDesignParameters::update_params(const VectorX params) {
    if (!_active)
        return;

    for (int i = 0;i < _ndof;i++) {
        _params[i] = params[_param_index[i]];
    }

    if (_type == DesignParameterType::TYPE2) {
        _body->_E_ji = math::compose_E(_params);
        _body->_E_ij = math::Einv(_body->_E_ji);
        _body->_A_ji = math::Ad(_body->_E_ji);
        _body->_A_ij = math::Ad(_body->_E_ij);
    } else if (_type == DesignParameterType::TYPE3) {
        for (int i = 0;i < _body->_contact_points.size();i++) {
            for (int j = 0;j < 3;j++) 
                _body->_contact_points[i](j) = _params[i * 3 + j];
        }
    } else if (_type == DesignParameterType::TYPE4) {
        _body->_mass = _params[0];
        _body->_Inertia.head(3) = _params.tail(3);
        _body->_Inertia.tail(3).setConstant(_params[0]);
        // std::cerr << "body " << _body->_name << ", I = " << _body->_Inertia.transpose() << std::endl;
    } else if (_type == DesignParameterType::TYPE6) {
        _body->_contact_scale = _params(0);
    } else {
        throw_error("[Error] Joint can only have type I, V design parameters");
    }
}

}