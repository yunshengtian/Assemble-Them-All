#pragma once
#include "Common.h"
#include "Utils.h"

namespace redmax {

class VirtualObject;

class VirtualObjectAnimator : public opengl_viewer::Animator {
public:
    VirtualObjectAnimator(VirtualObject* virtual_object) : _virtual_object(virtual_object) {}

    const Eigen::Matrix4f AnimatedModelMatrix(const float t);

    VirtualObject* _virtual_object;
};

class VirtualObject {
public:
    VirtualObject(std::string name);

    virtual void update_data(VectorX data) {}

    virtual Matrix4 get_transform_matrix() = 0;

    virtual void get_rendering_objects(
        std::vector<Matrix3Xf>& vertex_list, 
        std::vector<Matrix3Xi>& face_list,
        std::vector<opengl_viewer::Option>& option_list,
        std::vector<opengl_viewer::Animator*>& animator_list) {};

    const std::string _name;

protected:
    VirtualObjectAnimator* _animator;
};

}