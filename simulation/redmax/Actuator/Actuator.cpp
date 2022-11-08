#include "Actuator/Actuator.h"

namespace redmax {

Actuator::Actuator(std::string name, int ndof, VectorX ctrl_min, VectorX ctrl_max) {
    _name = name;
    _ndof = ndof;
    _ctrl_min = ctrl_min;
    _ctrl_max = ctrl_max;
    _index.clear();
    _u = VectorX::Zero(_ndof);
    _dofs = VectorX::Zero(_ndof);
    _dofs_vel = VectorX::Zero(_ndof);
}

Actuator::Actuator(std::string name, int ndof, dtype ctrl_min, dtype ctrl_max) {
    _name = name;
    _ndof = ndof;
    _ctrl_min = VectorX(ndof);
    _ctrl_max = VectorX(ndof);
    for (int i = 0;i < ndof;i++) {
        _ctrl_min[i] = ctrl_min;
        _ctrl_max[i] = ctrl_max;
    }
    _index.clear();
    _u = VectorX::Zero(_ndof);
    _dofs = VectorX::Zero(_ndof);
    _dofs_vel = VectorX::Zero(_ndof);
}

void Actuator::get_ctrl_range(VectorX& ctrl_min, VectorX& ctrl_max) {
    for (int i = 0;i < _ndof;i++) {
        ctrl_min[_index[i]] = _ctrl_min[i];
        ctrl_max[_index[i]] = _ctrl_max[i];
    }
}

void Actuator::set_u(const VectorX& u) {
    for (int i = 0;i < _ndof;i++) {
        // _u[i] = min(max(u[_index[i]], -1.), 1.);
        // _u[i] = min(max(u[_index[i]], _ctrl_min[i]), _ctrl_max[i]);
        _u[i] = u[_index[i]];
    }
}

}