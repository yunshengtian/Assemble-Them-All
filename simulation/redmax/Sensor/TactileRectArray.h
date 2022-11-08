#pragma once
#include "Sensor/Tactile.h"

namespace redmax {

class TactileRectArray : public Tactile {
public:
    TactileRectArray(Robot* robot, Body* body, string name, Vector3 rect_pos0, Vector3 rect_pos1, Vector3 axis0, Vector3 axis1, Vector2i resolution);

    void init();

    Vector3 _rect_pos0, _rect_pos1;
    Vector3 _axis0, _axis1;
    Vector2i _resolution;
};

};