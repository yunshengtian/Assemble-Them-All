#include "Sensor/TactileRectArray.h"
#include "Body/Body.h"

namespace redmax {

TactileRectArray::TactileRectArray(
    Robot* robot, Body* body, string name, Vector3 rect_pos0, Vector3 rect_pos1, 
    Vector3 axis0, Vector3 axis1, Vector2i resolution) : Tactile(robot, body, name) {
    
    _rect_pos0 = rect_pos0;
    _rect_pos1 = rect_pos1;
    _axis0 = axis0;
    _axis1 = axis1;
    _resolution = resolution;

    init();
}

void TactileRectArray::init() {
    // check whether input info is compatible
    dtype length_dir0 = (_rect_pos1 - _rect_pos0).dot(_axis0);
    dtype length_dir1 = (_rect_pos1 - _rect_pos0).dot(_axis1);
    if ((_rect_pos0 + length_dir0 * _axis0 + length_dir1 * _axis1 - _rect_pos1).norm() > 1e-5) {
        throw_error("Tactile info for " + _body->_name + " is incompatible");
    }

    // pre-compute the tactile sensor point locations
    Vector3 step_axis0 = length_dir0 / (_resolution(0) - 1) * _axis0;
    Vector3 step_axis1 = length_dir1 / (_resolution(1) - 1) * _axis1;

    _pos_i.clear(); _image_pos.clear(); _depth.clear();
    for (int i = 0;i < _resolution(0);i++)
        for (int j = 0;j < _resolution(1);j++) {
            _pos_i.push_back(_rect_pos0 + step_axis0 * i + step_axis1 * j);
            _image_pos.push_back(Vector2i(i, j));
            _depth.push_back(0.);
        }
}

};