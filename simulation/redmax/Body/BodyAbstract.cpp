#include "Body/BodyAbstract.h"
#include "tiny_obj_loader.h"
#include "Simulation.h"

namespace redmax {

// BodyAbstract::BodyAbstract(
//     Simulation* sim, Joint* joint, 
//     Matrix3 R_ji, Vector3 p_ji,
//     dtype mass, Vector3 Inertia,
//     std::string rendering_mesh_filename)
//     : Body(sim, joint, R_ji, p_ji, 0.0) {

//     if (rendering_mesh_filename != "") {
//         _rendering_mesh_exists = true;
//         load_visual_mesh(rendering_mesh_filename);
//     } else {
//         _rendering_mesh_exists = false;
//     }
    
//     _mass = mass;
//     _Inertia.head(3) = Inertia;
//     _Inertia.tail(3).setConstant(_mass);
// }

BodyAbstract::BodyAbstract(Simulation* sim, Joint* joint, 
    Matrix3 R_ji, Vector3 p_ji,
    dtype mass, VectorX Inertia,
    std::string visual_mesh_filename, Matrix3 visual_frame_R, Vector3 visual_frame_pos,  
    std::vector<Vector3> collision_points, Matrix3 collision_frame_R, Vector3 collision_frame_pos,
    Vector3 scale)
    : Body(sim, joint, 0.0) {
    
    _mass = mass;

    if (visual_mesh_filename != "") {
        _rendering_mesh_exists = true;
        load_visual_mesh(visual_mesh_filename, scale);
    } else {
        _rendering_mesh_exists = false;
    }

    if (Inertia.size() == 6) { // xx, xy, xz, yy, yz, zz // DEBUG: TO BE TESTED
        // get the principal axes for inertia tensor by eigenvalue decomposition
        // https://en.wikipedia.org/wiki/Moment_of_inertia#Principal_axes

        Matrix3 I;
        I(0, 0) = Inertia(0), I(0, 1) = I(1, 0) = Inertia(1), I(0, 2) = I(2, 0) = Inertia(2);
        I(1, 1) = Inertia(3), I(1, 2) = I(2, 1) = Inertia(4);
        I(2, 2) = Inertia(5);
        Eigen::SelfAdjointEigenSolver<Matrix3> eigensolver(I);
        Vector3 eigen_values = eigensolver.eigenvalues();
        Matrix3 eigen_vector3 = eigensolver.eigenvectors();
        _Inertia.head(3) = eigen_values;
        _Inertia.tail(3).setConstant(_mass);

        Matrix3 R_oldi_i = eigen_vector3;
        if (R_oldi_i.determinant() < 0.0) {
            R_oldi_i = -R_oldi_i;
        }

        // update the new body frame
        R_ji = R_ji * R_oldi_i;
        
        // update the visual mesh R
        visual_frame_R = R_oldi_i.transpose() * visual_frame_R;

        // update the collision mesh R
        collision_frame_R = R_oldi_i.transpose() * collision_frame_R;
    } else {
        _Inertia.head(3) = Inertia;
        _Inertia.tail(3).setConstant(_mass);
    }

    // transform the visual mesh into body frame
    if (_rendering_mesh_exists) {
        _V = (visual_frame_R * _V).colwise() + visual_frame_pos;
    }

    // transform the collision points into body frame
    for (int i = 0;i < collision_points.size();i++)
        collision_points[i] = collision_frame_R * collision_points[i].cwiseProduct(scale) + collision_frame_pos;
    
    set_contacts(collision_points);
    set_transform(R_ji, p_ji);
}


void BodyAbstract::load_visual_mesh(std::string filename, Vector3 scale) {
    std::vector<tinyobj::shape_t> obj_shape;
    std::vector<tinyobj::material_t> obj_material;
    tinyobj::attrib_t attrib;
    std::string err;
    tinyobj::LoadObj(&attrib, &obj_shape, &obj_material, &err, filename.c_str());

    int num_vertices = (int)attrib.vertices.size() / 3;
    _V.resize(3, num_vertices);
    for (int i = 0;i < num_vertices;i++) {
        _V.col(i) = Vector3(attrib.vertices[i * 3], 
            attrib.vertices[i * 3 + 1],
            attrib.vertices[i * 3 + 2]).cwiseProduct(scale);
    }
    
    int num_elements = (int)obj_shape[0].mesh.indices.size() / 3;
    _F.resize(3, num_elements);
    for (int i = 0;i < num_elements;i++) {
        _F.col(i) = Vector3i(obj_shape[0].mesh.indices[i * 3].vertex_index,
            obj_shape[0].mesh.indices[i * 3 + 1].vertex_index,
            obj_shape[0].mesh.indices[i * 3 + 2].vertex_index);
    }
}

Matrix3X BodyAbstract::get_vertices() const {
    return _V;
}

Matrix3Xi BodyAbstract::get_faces() const {
    return _F;
}

void BodyAbstract::set_rendering_mesh_vertices(const Matrix3X V) {
    _V = V;
}

void BodyAbstract::set_rendering_mesh(const Matrix3X V, const Matrix3Xi F) {
    _V = V;
    _F = F;
}
    
void BodyAbstract::get_rendering_objects(
    std::vector<Matrix3Xf>& vertex_list, 
    std::vector<Matrix3Xi>& face_list,
    std::vector<opengl_viewer::Option>& option_list,
    std::vector<opengl_viewer::Animator*>& animator_list) {

    _animator = new BodyAnimator(this);

    if (_rendering_mesh_exists) {
        opengl_viewer::Option object_option;

        object_option.SetBoolOption("smooth normal", false);
        object_option.SetVectorOption("ambient", _color(0), _color(1), _color(2));
        object_option.SetVectorOption("diffuse", 0.4f, 0.2368f, 0.1036f);
        object_option.SetVectorOption("specular", 0.774597f, 0.658561f, 0.400621f);
        object_option.SetFloatOption("shininess", 76.8f);

        Matrix3Xf vertex = _V.cast<float>();
        if (_sim->_options->_unit == "cm-g") 
            vertex = vertex / 10.;
        else
            vertex = vertex * 10.;

        vertex_list.push_back(vertex);
        face_list.push_back(_F);
        option_list.push_back(object_option);
    } else {
        opengl_viewer::Option object_option;

        object_option.SetVectorOption("ambient", _color(0), _color(1), _color(2));
        object_option.SetVectorOption("diffuse", 0.4f, 0.2368f, 0.1036f);
        object_option.SetVectorOption("specular", 0.474597f, 0.358561f, 0.200621f);
        object_option.SetFloatOption("shininess", 46.8f);

        Matrix3Xf cube_vertex;
        Matrix3Xi cube_face;
        Matrix2Xf cube_uv;

        opengl_viewer::ReadFromObjFile(
            std::string(GRAPHICS_CODEBASE_SOURCE_DIR) + "/resources/meshes/cube.obj",
            cube_vertex, cube_face, cube_uv);
        
        Vector3 length;
        length(0) = sqrt((_Inertia(2) + _Inertia(1) - _Inertia(0)) * 6 / _mass);
        length(1) = sqrt((_Inertia(2) + _Inertia(0) - _Inertia(1)) * 6 / _mass);
        length(2) = sqrt((_Inertia(1) + _Inertia(0) - _Inertia(2)) * 6 / _mass);
        
        for (int i = 0;i < 3;i++)
            cube_vertex.row(i) *= (float)length(i);
        if (_sim->_options->_unit == "cm-g") 
            cube_vertex /= 10.;
        else
            cube_vertex *= 10.;
        
        vertex_list.push_back(cube_vertex);
        face_list.push_back(cube_face);
        option_list.push_back(object_option);
    }

    animator_list.push_back(_animator);
}

}