#include "Force/ForceCable.h"
#include "Body/Body.h"
#include "Simulation.h"
#include "Joint/Joint.h"

namespace redmax {

ForceCable::ForceCable(Simulation* sim, dtype stiffness, dtype damping) : ForceSpringMultiPointGeneric(sim) {
    _stiffness = stiffness;
    _damping = damping;
    _l_rest = 0.;
}

void ForceCable::set_stiffness(dtype stiffness) { 
    _stiffness = stiffness;
}

void ForceCable::set_damping(dtype damping) {
    _damping = damping;
}

void ForceCable::init() {
    _l_rest = 0.;
    int n = _bodies.size();
    for (int k = 0;k < n - 1;k++) {
        const Matrix3& R0 = _bodies[k]->_E_0i.topLeftCorner(3, 3);
        const Vector3& p0 = _bodies[k]->_E_0i.topRightCorner(3, 1);
        const Matrix3& R1 = _bodies[k + 1]->_E_0i.topLeftCorner(3, 3);
        const Vector3& p1 = _bodies[k + 1]->_E_0i.topRightCorner(3, 1);
        Vector3 xw0 = R0 * _xls[k] + p0;
        Vector3 xw1 = R1 * _xls[k + 1] + p1;
        _l_rest += (xw0 - xw1).norm();
    }
}

void ForceCable::computeSpringForce(dtype l, dtype ldot, dtype& f) {
    dtype strain = (l - _l_rest) / _l_rest;
    dtype dstrain = ldot / _l_rest;
    if (strain > 0) {
        f = _stiffness * strain + _damping * dstrain;
    } else {
        f = 0.;
    }
}

void ForceCable::computeSpringForceWithDerivative(dtype l, dtype ldot, dtype& f, dtype& df_dl, dtype& df_dldot, bool verbose) {
    dtype strain = (l - _l_rest) / _l_rest;
    dtype dstrain = ldot / _l_rest;
    if (strain > 0) {
        f = _stiffness * strain + _damping * dstrain;
        df_dl = _stiffness / _l_rest;
        df_dldot = _damping / _l_rest;
    } else {
        f = 0.;
        df_dl = df_dldot = 0.;
    }
}

}