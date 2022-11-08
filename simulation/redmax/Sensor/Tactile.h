#pragma once
#include "Common.h"
#include "Utils.h"

namespace redmax {

class Body;
class Robot;

class Tactile {
public:
    Tactile(Robot* robot, Body* body, string name) : _robot(robot), _body(body), _name(name) {}

    virtual void init() = 0; // init class member variables

    void compute_tactile_values();    

    Robot* _robot;
    Body* _body;
    string _name;

    std::vector<Vector3> _pos_i; // positions of each tactile sensor point in the associated body frame
    std::vector<Vector2i> _image_pos; // positions of each tactile sensor point in the tactile sensor image
    std::vector<dtype> _depth;
};

};