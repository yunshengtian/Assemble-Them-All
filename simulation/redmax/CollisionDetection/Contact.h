#pragma once
#include "Utils.h"

namespace redmax {

class Contact {
public:
    Vector3 _xi;     // contact point in body frame.
    Vector3 _xw;     // contact point in world frame.
    dtype _d;        // penetration depth.
    Vector3 _normal; // contact normal
    int _id;         // id of contact point in body

    Contact(Vector3 xi, Vector3 xw, dtype d, Vector3 normal, int id = 0) {
        _xi = xi;
        _xw = xw;
        _d = d;
        _normal = normal;
        _id = id;
    }
};
}