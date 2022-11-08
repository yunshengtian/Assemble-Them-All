#include "VirtualObject/VirtualObject.h"

namespace redmax {
const Eigen::Matrix4f VirtualObjectAnimator::AnimatedModelMatrix(const float t) {
    Matrix4f model_matrix = (_virtual_object->get_transform_matrix()).cast<float>();
    model_matrix.topRightCorner(3, 1) /= 10.; // scale for better visualization
    return model_matrix;
}

VirtualObject::VirtualObject(std::string name) : _name(name) {}

}