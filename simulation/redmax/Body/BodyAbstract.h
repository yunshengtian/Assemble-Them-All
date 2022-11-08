#pragma once
#include "Body/Body.h"

namespace redmax {

/**
 * BodyAbstract is an abstract body.
 * The mass, inertia and contacts are directly specified instead of being inferred from object.
 * the rendering_mesh is only used for rendering but not related to any simulation dynamics
 * */
class BodyAbstract : public Body {
public:
    bool _rendering_mesh_exists;
    Matrix3X _V;            // vertices
    Matrix3Xi _F;            // face elements

    // BodyAbstract(Simulation* sim, Joint* joint, 
    //     Matrix3 R_ji, Vector3 p_ji,
    //     dtype mass, Vector3 Inertia,
    //     std::string rendering_mesh_filename = "");
    
    BodyAbstract(Simulation* sim, Joint* joint, 
        Matrix3 R_ji, Vector3 p_ji,
        dtype mass, VectorX Inertia,
        std::string visual_mesh_filename = "", Matrix3 visual_frame_R = Matrix3::Identity(), Vector3 visual_frame_pos = Vector3::Zero(), 
        std::vector<Vector3> collision_points = std::vector<Vector3>(), Matrix3 collision_frame_R = Matrix3::Identity(), Vector3 collision_frame_pos = Vector3::Zero(),
        Vector3 scale = Vector3::Ones());

    void load_visual_mesh(std::string filename, Vector3 scale);

    // mesh
    Matrix3X get_vertices() const;
    Matrix3Xi get_faces() const;

    void set_rendering_mesh_vertices(const Matrix3X V);
    void set_rendering_mesh(const Matrix3X V, const Matrix3Xi F);
    
    // rendering
    void get_rendering_objects(
        std::vector<Matrix3Xf>& vertex_list, 
        std::vector<Matrix3Xi>& face_list,
        std::vector<opengl_viewer::Option>& option_list,
        std::vector<opengl_viewer::Animator*>& animator_list);
};

}