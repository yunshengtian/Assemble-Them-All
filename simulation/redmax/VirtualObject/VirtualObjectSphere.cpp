#include "VirtualObject/VirtualObjectSphere.h"

namespace redmax {

VirtualObjectSphere::VirtualObjectSphere(std::string name, Vector3 pos_world, dtype radius, Vector3 color) 
    : VirtualObject(name), _pos_world(pos_world), _radius(radius), _color(color) {}

void VirtualObjectSphere::update_data(VectorX data) {
    if (data.size() != 3 && data.size() != 6) {
        throw_error("Update virtual sphere with wrong data size." );
    }

    if (data.size() == 3) {
        _pos_world = data.head(3);
    } else {
        _pos_world = data.head(3);
        _color = data.tail(3);
    }
}

Matrix4 VirtualObjectSphere::get_transform_matrix() {
    Matrix4 E = Matrix4::Identity();
    E.topRightCorner(3, 1) = _pos_world;
    return E;
}

// rendering
void VirtualObjectSphere::get_rendering_objects(
    std::vector<Matrix3Xf>& vertex_list, 
    std::vector<Matrix3Xi>& face_list,
    std::vector<opengl_viewer::Option>& option_list,
    std::vector<opengl_viewer::Animator*>& animator_list) {
    
    Matrix3Xf vertex;
    Matrix3Xi face;
    Matrix2Xf uv;
    opengl_viewer::Option object_option;

    opengl_viewer::ReadFromObjFile(
        std::string(GRAPHICS_CODEBASE_SOURCE_DIR) + "/resources/meshes/sphere.obj",
        vertex, face, uv);

    for (int i = 0;i < 3;i++)
        vertex.row(i) *= (float)_radius / 10.;

    object_option.SetBoolOption("smooth normal", false);
    object_option.SetVectorOption("ambient", _color(0), _color(1), _color(2));
    object_option.SetVectorOption("diffuse", 0.4f, 0.2368f, 0.1036f);
    object_option.SetVectorOption("specular", 0.774597f, 0.458561f, 0.200621f);
    object_option.SetFloatOption("shininess", 76.8f);

    _animator = new VirtualObjectAnimator(this);

    vertex_list.push_back(vertex);
    face_list.push_back(face);
    option_list.push_back(object_option);
    animator_list.push_back(_animator);
}

}