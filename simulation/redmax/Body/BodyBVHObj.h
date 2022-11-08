#pragma once
#include "Body/BodyMeshObj.h"
#include "Body/BVH/BVHEngine.h"
#include "array3.h"
#include "vec.h"

namespace redmax {

/**
 * BodyBVHObj define a body loaded from an .obj file (triangle surface mesh) and compute a BVH for it.
 **/
class BodyBVHObj : public BodyMeshObj {
public:
    Matrix3X _BVH_V;                  // BVH mesh vertices
    Matrix3Xi _BVH_F;                 // BVH mesh face elements
    std::shared_ptr<BVHEngine> _BVH;  // BVH engine
    Vector3 _min_box, _max_box;       // bounding box of SDF grid
    dtype _dx, _dy, _dz;              // delta of SDF grid
    dtype _col_th;                    // collision threshold
    sdfgen::Array3<dtype> _SDF;       // voxel grid of SDF
    bool _use_sdf_gradient = false;

    BodyBVHObj(Simulation* sim, Joint* joint,
                    std::string filename_mesh, 
                    std::string filename_BVH_mesh,
                    Matrix3 R, Vector3 p, 
                    dtype dx, int res, dtype col_th,
                    TransformType transform_type = BODY_TO_JOINT,
                    dtype density = (dtype)1.0,
                    Vector3 scale = Vector3::Ones(),
                    bool adaptive_sample = false,
                    bool load_sdf = false,
                    bool save_sdf = false);
    void load_BVH_mesh(std::string filename_BVH_mesh, Vector3 scale);

    // mesh
    Matrix3X get_vertices() const;
    Matrix3Xi get_faces() const;

    /**
     * analytical distance function
     * @param:
     * xw: the location of the query position in world frame
     * @return: contact distance of xw
    **/
    void distance(Vector3 xw, dtype &d);
    void distance_local(Vector3 x, dtype &d);
    void distance_parallel(Matrix3X xw, VectorX &d);

    void distance_normal(Vector3 xw, dtype &d, Vector3 &dd_dx);
    void distance_normal_local(Vector3 x, dtype &d, Vector3 &dd_dx);
    void distance_normal_parallel(Matrix3X xw, VectorX &d, Matrix3X &dd_dx);

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
    void collision_parallel(Vector3 xw, Vector3 xw_dot, dtype d, Vector3 dd_dx, /* input */
                    Vector3 &n,  /* output */
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
    void collision_parallel(Vector3 xw, Vector3 xw_dot, dtype d, Vector3 dd_dx, /* input */
                    Vector3 &n,  /* output */
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

    void precompute_bounding_box();

    void precompute_SDF(dtype dx, int res);
    bool load_SDF();
    void save_SDF();
    dtype query_SDF(Vector3 x);
    Vector3 query_dSDF(Vector3 x);
    Matrix3 query_ddSDF(Vector3 x);
};

}
