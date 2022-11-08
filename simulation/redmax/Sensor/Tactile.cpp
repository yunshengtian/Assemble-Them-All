#include "Body/Body.h"
#include "Body/BodyPrimitiveShape.h"
#include "Robot.h"
#include "Sensor/Tactile.h"

namespace redmax {

void Tactile::compute_tactile_values() {
    Matrix4 E_wi = _body->_E_0i;
    Matrix3 R_wi = E_wi.topLeftCorner(3, 3);
    Vector3 p_wi = E_wi.topRightCorner(3, 1);

    for (int i = 0;i < _depth.size();i++)
        _depth[i] = 0.;

    for (auto body : _robot->_bodies) {
        if (dynamic_cast<BodyPrimitiveShape*>(const_cast<Body*>(body)) != nullptr) {
            BodyPrimitiveShape* primitive_body = dynamic_cast<BodyPrimitiveShape*>(const_cast<Body*>(body));
            for (int i = 0;i < _pos_i.size();i++) {
                Vector3 pos_w = R_wi * _pos_i[i] + p_wi;
                dtype d = primitive_body->distance(pos_w);
                if (d < 0.) {
                    _depth[i] = max(_depth[i], -d);
                }
            }
        }
    }
}

};