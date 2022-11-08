#pragma once
#include "Common.h"
#include "Utils.h"
#include "BodyDesignParameters.h"
#include "CollisionDetection/Contact.h"
#include <set>

namespace redmax {

class Joint;
class Simulation;
class Body;

class BodyAnimator : public opengl_viewer::Animator {
public:
    BodyAnimator(Body* body) : _body(body) {}

    const Eigen::Matrix4f AnimatedModelMatrix(const float t);

    Body* _body;
};

class Body {
public:
    // simulation
    Simulation* _sim;

    // name
    std::string _name = "";

    // structure
    Joint* _joint;
    Body* _parent;

    // contacts
    std::vector<Vector3> _contact_points;
    std::set<std::string> _contact_body_history;

    std::pair<Vector3, Vector3> _bounding_box; // bounding box of the body

    // index
    vector<int> _index;                 // index for phi in maximal coordinates

    // data
    dtype _density;
    dtype _mass;
    Vector6 _Inertia;                   // diagonal of mass matrix M_i in body frame
    dtype _contact_scale;               // contact coefficients scale
    
    // constants
    SE3 _E_ij, _E_ji;                   // transformation between body i and its joint j, constant
    Matrix6 _A_ij, _A_ji;             // ajoint matrix between body i and its joint j, constant

    // variables
    SE3 _E_0i, _E_i0;                   // transformation between body i and world frame
    SE3 _E_ip, _E_pi;                   // transformation between body i and its parant
    Matrix6 _A_ip, _A_pi;             // adjoint matrix between body i and its parent body
    Matrix6 _A_ip_dot;                 // time derivative of _Ad_ip, updated based on notes eq. (97)

    // phi
    se3 _phi;                           // spatial velocity
    Vector6 _phi_dot;                   // TODO: leave to implement

    // external force
    Vector6 _f_external;

    // design parameters
    BodyDesignParameters _design_params_2, _design_params_3, _design_params_4, _design_params_6;

    // derivatives
    JacobianMatrixVector _dAip_dq, _dAipdot_dq; // represented in body frame.

    // design derivatives
    JacobianMatrixVector _dAij_dp2;
    JacobianMatrixVector _dE0i_dp1, _dE0i_dp2;
    MatrixX _dphi_dp1, _dphi_dp2;

    // rendering
    Vector3f _color;
    std::string _texture_path;
    bool _use_texture;

    Body(Simulation* sim, Joint* joint, dtype density);

    Body(Simulation* sim, Joint* joint, Matrix3 R_ji, Vector3 p_ji, dtype density);

    ~Body();

    // init body
    void init();

    void set_transform(Matrix3 R_ji, Vector3 p_ji);

    void set_contacts(vector<Vector3>& contacts);

    // set rendering color
    void set_color(Vector3 color);
    void set_texture(std::string texture_path);

    // mesh
    virtual Matrix3X get_vertices() const = 0;
    virtual Matrix3Xi get_faces() const = 0;

    // external force
    Vector6 get_external_force() const;
    void set_external_force(Vector6 force);
    void reset_external_force();

    // update
    void update(bool design_gradient = false);

    // activate design parameters
    void activate_design_parameters_type_2(bool active = true);
    void activate_design_parameters_type_3(bool active = true);
    void activate_design_parameters_type_4(bool active = true);
    void activate_design_parameters_type_6(bool active = true);

    // compute maximal force
    void computeMaximalForce(VectorX& fm);
    void computeMaximalForceWithDerivative(VectorX& fm, MatrixX& Km, MatrixX& Dm);
    void computeMaximalForceWithDerivative(
        VectorX& fm, 
        MatrixX& Km, MatrixX& Dm,
        MatrixX& dfm_dp);

    void test_derivatives_runtime();

    // contact points
    std::vector<Vector3> get_contact_points() const { return _contact_points; };

    // contacts
    void add_contact_body(std::string body_name);
    std::vector<std::string> get_contact_bodies() const;
    void clear_contact_bodies();
    
    // coordinates transform
    Vector3 position_in_world(Vector3 pos) const;
    // se3 spatial_velocity_in_world(Vector3 pos);
    Vector3 velocity_in_world(Vector3 pos) const;

    // rendering
    virtual void get_rendering_objects(
        std::vector<Matrix3Xf>& vertex_list, 
        std::vector<Matrix3Xi>& face_list,
        std::vector<opengl_viewer::Option>& option_list,
        std::vector<opengl_viewer::Animator*>& animator_list) {};

    std::pair<Vector3, Vector3> get_AABB();
    
    // // contact related
    // virtual dtype distance(Vector3 xw) { return (dtype)0.; };
    
    // /**
    //  * analytical distance function return dist, normal, x2i, ddot, tdot
    //  * @return
    //  * dist: the distance of xw
    //  * normal: the normalized contact normal
    //  * xi2: the contact position on this body
    //  * ddot: the magnitude of the velocity on normal direction
    //  * tdot: the tangential velocity
    //  **/
    // virtual void collision(Vector3 xw, Vector3 xw_dot, /* input */
    //                 dtype &d, Vector3 &n,  /* output */
    //                 dtype &ddot, Vector3 &tdot,
    //                 Vector3 &xi2) {};
    
    // /**
    //  * analytical distance function return dist, normal, x2i, ddot, tdot and derivatives
    //  * @return
    //  * dist: the distance of xw
    //  * normal: the normalized contact normal
    //  * xi2: the contact position on this body
    //  * ddot: the magnitude of the velocity on normal direction
    //  * tdot: the tangential velocity
    //  * derivatives
    //  **/
    // virtual void collision(Vector3 xw, Vector3 xw_dot, /* input */
    //                 dtype &d, Vector3 &n,  /* output */
    //                 dtype &ddot, Vector3 &tdot,
    //                 Vector3 &xi2,
    //                 RowVector3 &dd_dxw, RowVector6 &dd_dq2, /* derivatives for d */ 
    //                 Matrix3 &dn_dxw, Matrix36 &dn_dq2, /* derivatives for n */
    //                 RowVector3 &dddot_dxw, RowVector3 &dddot_dxwdot, /* derivatives for ddot */
    //                 RowVector6 &dddot_dq2, RowVector6 &dddot_dphi2,
    //                 Matrix3 &dtdot_dxw, Matrix3 &dtdot_dxwdot, /* derivatives for tdot */
    //                 Matrix36 &dtdot_dq2, Matrix36 &dtdot_dphi2,
    //                 Matrix3 &dxi2_dxw, Matrix36 &dxi2_dq2 /* derivatives for xi2 */) {}

protected:
    BodyAnimator* _animator;

    virtual void update_design_derivatives();
};

}