#pragma once
#include "Body/BodyPrimitiveShape.h"

namespace redmax {

class BodyCylinder : public BodyPrimitiveShape {
public:
    dtype _height;
    dtype _radius;

    Matrix3X _V;
    Matrix3Xi _F;

    BodyCylinder(Simulation *sim, Joint *joint, dtype height, dtype radius, 
                Matrix3 R_ji, Vector3 p_ji, dtype density);

    void computeMassMatrix();

    // contact candidate points
    void precompute_contact_points();
    std::vector<Vector3> get_general_contact_points();

    // mesh
    Matrix3X get_vertices() const;
    Matrix3Xi get_faces() const;

    // rendering
    void get_rendering_objects(
        std::vector<Matrix3Xf>& vertex_list, 
        std::vector<Matrix3Xi>& face_list,
        std::vector<opengl_viewer::Option>& option_list,
        std::vector<opengl_viewer::Animator*>& animator_list);

    /**
     * analytical distance function
     * @param:
     * xw: the location of the query position in world frame
     * @return: contact distance of xw
    **/
    dtype distance(Vector3 xw);

    /**
     * analytical distance function return dist, normal, x2i, ddot, tdot
     * @return
     * dist: the distance of xw
     * normal: the normalized contact normal
     * xi2: the contact position on this body
     * ddot: the magnitude of the velocity on normal direction
     * tdot: the tangential velocity
     **/
    void collision(Vector3 xw, Vector3 xw_dot, /* input */
                    dtype &d, Vector3 &n,  /* output */
                    dtype &ddot, Vector3 &tdot,
                    Vector3 &xi2);

    /**
     * analytical distance function return dist, normal, x2i, ddot, tdot and derivatives
     * @return
     * dist: the distance of xw
     * normal: the normalized contact normal
     * xi2: the contact position on this body
     * ddot: the magnitude of the velocity on normal direction
     * tdot: the tangential velocity
     * derivatives
     **/
    void collision(Vector3 xw, Vector3 xw_dot, /* input */
                    dtype &d, Vector3 &n,  /* output */
                    dtype &ddot, Vector3 &tdot,
                    Vector3 &xi2,
                    RowVector3 &dd_dxw, RowVector6 &dd_dq2, /* derivatives for d */ 
                    Matrix3 &dn_dxw, Matrix36 &dn_dq2, /* derivatives for n */
                    RowVector3 &dddot_dxw, RowVector3 &dddot_dxwdot, /* derivatives for ddot */
                    RowVector6 &dddot_dq2, RowVector6 &dddot_dphi2,
                    Matrix3 &dtdot_dxw, Matrix3 &dtdot_dxwdot, /* derivatives for tdot */
                    Matrix36 &dtdot_dq2, Matrix36 &dtdot_dphi2,
                    Matrix3 &dxi2_dxw, Matrix36 &dxi2_dq2 /* derivatives for xi2 */);


    void test_collision_derivatives();
    void test_collision_derivatives_runtime(Vector3 xw, Vector3 xw_dot);

private:
    void load_mesh();
};

}