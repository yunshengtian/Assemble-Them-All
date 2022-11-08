#include "Body/BodySphere.h"
#include "Simulation.h"

namespace redmax {

BodySphere::BodySphere(
    Simulation *sim, Joint *joint, dtype radius,
    Matrix3 R_ji, Vector3 p_ji,
    dtype density)
    : BodyPrimitiveShape(sim, joint, R_ji, p_ji, density) {
    
    _radius = radius;

    computeMassMatrix();
}

void BodySphere::computeMassMatrix() {
    _mass = 4. / 3. * constants::pi * _radius * _radius * _radius * _density;
    _Inertia.head(3) = Vector3::Constant(0.4 * _mass * _radius * _radius);
    _Inertia.tail(3) = Vector3::Constant(_mass);
}

Matrix3X BodySphere::get_vertices() const {
    std::cerr << "[Error] get_vertices not implemented" << std::endl;
    throw "error";
}

Matrix3Xi BodySphere::get_faces() const {
    std::cerr << "[Error] get_faces not implemented" << std::endl;
    throw "error";
}

void BodySphere::get_rendering_objects(
    std::vector<Matrix3Xf>& vertex_list, 
    std::vector<Matrix3Xi>& face_list,
    std::vector<opengl_viewer::Option>& option_list,
    std::vector<opengl_viewer::Animator*>& animator_list) {
    
    Matrix3Xf sphere_vertex;
    Matrix3Xi sphere_face;
    Matrix2Xf sphere_uv;
    opengl_viewer::Option object_option;

    opengl_viewer::ReadFromObjFile(
        std::string(GRAPHICS_CODEBASE_SOURCE_DIR) + "/resources/meshes/sphere.obj",
        sphere_vertex, sphere_face, sphere_uv);
    
    sphere_vertex *= (float)_radius;
    if (_sim->_options->_unit == "cm-g") 
        sphere_vertex /= 10.;
    else
        sphere_vertex *= 10.;
        
    if (!_use_texture) {
        object_option.SetVectorOption("ambient", _color(0), _color(1), _color(2));
        object_option.SetVectorOption("diffuse", 0.4f, 0.2368f, 0.1036f);
        object_option.SetVectorOption("specular", 0.474597f, 0.358561f, 0.200621f);
        object_option.SetFloatOption("shininess", 46.8f);
    } else {
        opengl_viewer::Image checker_texture;
        checker_texture.Initialize(std::string(GRAPHICS_CODEBASE_SOURCE_DIR) + "//" + _texture_path);

        object_option.SetBoolOption("smooth normal", false);

        object_option.SetVectorOption("ambient", 0.7f, 0.7f, 0.7f);
        object_option.SetVectorOption("diffuse", 1.0f, 1.0f, 1.0f);
        object_option.SetVectorOption("specular", 1.0f, 1.0f, 1.0f);
        object_option.SetFloatOption("shininess", 1.5f);
        
        object_option.SetMatrixOption("uv", sphere_uv);
        object_option.SetMatrixOption("texture", checker_texture.rgb_data());
        object_option.SetIntOption("texture row num", checker_texture.row_num());
        object_option.SetIntOption("texture col num", checker_texture.col_num());
        object_option.SetStringOption("texture mag filter", "nearest");
    }

    _animator = new BodyAnimator(this);

    vertex_list.push_back(sphere_vertex);
    face_list.push_back(sphere_face);
    option_list.push_back(object_option);
    animator_list.push_back(_animator);
}

dtype BodySphere::distance(Vector3 xw) {
    Vector3 p = _E_0i.topRightCorner(3, 1);
    dtype d = (xw - p).norm() - _radius;
    return d;
}

// void BodySphere::collision(
//     Vector3 xw, Vector3 xw_dot, /* input */
//     dtype &d, Vector3 &n,  /* output */
//     dtype &ddot, Vector3 &tdot,
//     Vector3 &xi2) {
    
//     Matrix3 R2 = _E_0i.topLeftCorner(3, 3);
//     Vector3 p2 = _E_0i.topRightCorner(3, 1);
//     Vector3 v2_dot = _phi.tail(3);
//     Vector3 p2_dot = R2 * v2_dot;
//     Vector3 vec = xw - p2;
//     Vector3 vec_dot = xw_dot - p2_dot;
//     dtype l = (xw - p2).norm();

//     d = l - _radius;
//     n = vec / l;
//     ddot = vec.dot(vec_dot) / l;
//     tdot = vec_dot - ddot * n;
//     xi2 = _radius / l * R2.transpose() * vec;

//     // std::cerr << "p2_dot = " << p2_dot.transpose() << std::endl;
// }

// void BodySphere::collision(
//     Vector3 xw, Vector3 xw_dot, /* input */
//     dtype &d, Vector3 &n,  /* output */
//     dtype &ddot, Vector3 &tdot,
//     Vector3 &xi2,
//     RowVector3 &dd_dxw, RowVector6 &dd_dq2, /* derivatives for d */ 
//     Matrix3 &dn_dxw, Matrix36 &dn_dq2, /* derivatives for n */
//     RowVector3 &dddot_dxw, RowVector3 &dddot_dxwdot, /* derivatives for ddot */
//     RowVector6 &dddot_dq2, RowVector6 &dddot_dphi2,
//     Matrix3 &dtdot_dxw, Matrix3 &dtdot_dxwdot, /* derivatives for tdot */
//     Matrix36 &dtdot_dq2, Matrix36 &dtdot_dphi2,
//     Matrix3 &dxi2_dxw, Matrix36 &dxi2_dq2 /* derivatives for xi2 */) {
    
//     Matrix3 I = Matrix3::Identity();
//     Matrix3 R2 = _E_0i.topLeftCorner(3, 3);
//     Vector3 p2 = _E_0i.topRightCorner(3, 1);
//     Vector3 v2_dot = _phi.tail(3);
//     Vector3 p2_dot = R2 * v2_dot;
//     Vector3 vec = xw - p2;
//     Vector3 vec_dot = xw_dot - p2_dot;
//     dtype l = vec.norm();

//     /**************** values ****************/
//     d = l - _radius;
//     n = vec / l;
//     ddot = vec.dot(vec_dot) / l;
//     tdot = vec_dot - ddot * n;
//     xi2 = _radius / l * R2.transpose() * vec;

//     /**************** derivatives ****************/
//     // d
//     dd_dxw = vec.transpose() / l;
//     dd_dq2.head(3).setZero();
//     dd_dq2.tail(3) = -vec.transpose() * R2 / l;

//     // n
//     dn_dxw = I / l - vec * vec.transpose() / pow(l, 3);
//     dn_dq2.leftCols(3).setZero();
//     dn_dq2.rightCols(3) = -R2 / l + vec * (vec.transpose() * R2) / pow(l, 3);

//     // ddot
//     dddot_dxw = vec_dot.transpose() / l - vec.transpose() * vec_dot * vec.transpose() / pow(l, 3);
//     dddot_dxwdot = vec.transpose() / l;
//     dddot_dq2.head(3) = vec.transpose() * R2 * math::skew(v2_dot) / l;
//     dddot_dq2.tail(3) = (-vec_dot.transpose() / l + vec.transpose() * vec_dot * vec.transpose() / pow(l, 3)) * R2;
//     dddot_dphi2.head(3).setZero();
//     dddot_dphi2.tail(3) = -vec.transpose() * R2 / l;

//     // tdot
//     dtdot_dxw = -n * dddot_dxw - ddot * dn_dxw;
//     dtdot_dxwdot = I - n * dddot_dxwdot;
//     dtdot_dq2.leftCols(3) = R2 * math::skew(v2_dot) - n * dddot_dq2.head(3) - ddot * dn_dq2.leftCols(3);
//     dtdot_dq2.rightCols(3) = -n * dddot_dq2.tail(3) - ddot * dn_dq2.rightCols(3);
//     dtdot_dphi2.leftCols(3).setZero();
//     dtdot_dphi2.rightCols(3) = -R2 - n * dddot_dphi2.tail(3);

//     // xi2
//     dxi2_dxw = _radius * R2.transpose() * (I / l - vec * vec.transpose() / pow(l, 3));
//     dxi2_dq2.leftCols(3) = math::skew(xi2);
//     dxi2_dq2.rightCols(3) = -_radius / l * (I - (R2.transpose() * vec) * (vec.transpose() * R2) / pow(l, 2));
// }

void BodySphere::collision(
    Vector3 xw, Vector3 xw_dot, /* input */
    dtype &d, Vector3 &n,  /* output */
    dtype &ddot, Vector3 &tdot,
    Vector3 &xi2) {
    
    Matrix3 I = Matrix3::Identity();
    Matrix3 R2 = _E_0i.topLeftCorner(3, 3);
    Vector3 p2 = _E_0i.topRightCorner(3, 1);
    Vector3 w2_dot = _phi.head(3);
    Vector3 v2_dot = _phi.tail(3);
    Vector3 p2_dot = R2 * v2_dot;
    Vector3 vec = xw - p2;
    Vector3 vec_dot = xw_dot - p2_dot;
    dtype l = vec.norm();

    d = l - _radius;
    n = vec / l;
    ddot = vec.dot(vec_dot) / l;
    xi2 = _radius / l * R2.transpose() * vec;

    Vector3 xw2_dot = R2 * (w2_dot.cross(xi2) + v2_dot);
    Vector3 vw = xw_dot - xw2_dot;
    tdot = (I - n * n.transpose()) * vw;
}

// refer to https://www.overleaf.com/project/5cd0a78e0e41fe0f8632a42a
// section 1.4: General Penalty-based Contact
void BodySphere::collision(
    Vector3 xw, Vector3 xw_dot, /* input */
    dtype &d, Vector3 &n,  /* output */
    dtype &ddot, Vector3 &tdot,
    Vector3 &xi2,
    RowVector3 &dd_dxw, RowVector6 &dd_dq2, /* derivatives for d */ 
    Matrix3 &dn_dxw, Matrix36 &dn_dq2, /* derivatives for n */
    RowVector3 &dddot_dxw, RowVector3 &dddot_dxwdot, /* derivatives for ddot */
    RowVector6 &dddot_dq2, RowVector6 &dddot_dphi2,
    Matrix3 &dtdot_dxw, Matrix3 &dtdot_dxwdot, /* derivatives for tdot */
    Matrix36 &dtdot_dq2, Matrix36 &dtdot_dphi2,
    Matrix3 &dxi2_dxw, Matrix36 &dxi2_dq2/* derivatives for xi2 */) {
    
    Matrix3 I = Matrix3::Identity();
    Matrix3 R2 = _E_0i.topLeftCorner(3, 3);
    Vector3 p2 = _E_0i.topRightCorner(3, 1);
    Vector3 w2_dot = _phi.head(3);
    Vector3 v2_dot = _phi.tail(3);
    Vector3 p2_dot = R2 * v2_dot;
    Vector3 vec = xw - p2;
    Vector3 vec_dot = xw_dot - p2_dot;
    dtype l = vec.norm();

    /**************** values ****************/
    d = l - _radius;
    n = vec / l;
    ddot = vec.dot(vec_dot) / l;
    xi2 = _radius / l * R2.transpose() * vec;
    Vector3 xw2_dot = R2 * (w2_dot.cross(xi2) + v2_dot);
    Vector3 vw = xw_dot - xw2_dot;
    tdot = (I - n * n.transpose()) * vw;

    /**************** derivatives ****************/

    // l
    RowVector3 dl_dxw = vec.transpose() / l;
    RowVector6 dl_dq2;
    dl_dq2.head(3).setZero();
    dl_dq2.tail(3) = -vec.transpose() * R2 / l;

    // d
    dd_dxw = dl_dxw;
    dd_dq2 = dl_dq2;

    // n
    dn_dxw = I / l - vec * vec.transpose() / pow(l, 3);
    dn_dq2.leftCols(3).setZero();
    dn_dq2.rightCols(3) = -R2 / l + vec * (vec.transpose() * R2) / pow(l, 3);

    // ddot
    dddot_dxw = vec_dot.transpose() / l - vec.transpose() * vec_dot * vec.transpose() / pow(l, 3);
    dddot_dxwdot = vec.transpose() / l;
    dddot_dq2.head(3) = vec.transpose() * R2 * math::skew(v2_dot) / l;
    dddot_dq2.tail(3) = (-vec_dot.transpose() / l + vec.transpose() * vec_dot * vec.transpose() / pow(l, 3)) * R2;
    dddot_dphi2.head(3).setZero();
    dddot_dphi2.tail(3) = -vec.transpose() * R2 / l;

    // xi2
    dxi2_dxw = _radius * R2.transpose() * (I / l - vec * vec.transpose() / pow(l, 3));
    dxi2_dq2.leftCols(3) = math::skew(xi2);
    dxi2_dq2.rightCols(3) = -_radius / l * (I - (R2.transpose() * vec) * (vec.transpose() * R2) / pow(l, 2));

    // xw2_dot
    Matrix3 dxw2dot_dxw = R2 * math::skew(w2_dot) * dxi2_dxw;
    Matrix3 dxw2dot_dxwdot = Matrix3::Zero();
    Matrix36 dxw2dot_dq2;
    dxw2dot_dq2.leftCols(3) = -R2 * math::skew(w2_dot.cross(xi2) + v2_dot) + R2 * math::skew(w2_dot) * dxi2_dq2.leftCols(3);
    dxw2dot_dq2.rightCols(3) = R2 * math::skew(w2_dot) * dxi2_dq2.rightCols(3);
    Matrix36 dxw2dot_dphi2;
    Vector3 e1 = Vector3::UnitX(), e2 = Vector3::UnitY(), e3 = Vector3::UnitZ();
    dxw2dot_dphi2.col(0) = R2 * math::skew(e1) * xi2;
    dxw2dot_dphi2.col(1) = R2 * math::skew(e2) * xi2;
    dxw2dot_dphi2.col(2) = R2 * math::skew(e3) * xi2;
    dxw2dot_dphi2.rightCols(3) = R2;

    // tdot
    dtdot_dxw = -dxw2dot_dxw - n.transpose() * vw * dn_dxw - n * (vw.transpose() * dn_dxw - n.transpose() * dxw2dot_dxw);
    dtdot_dxwdot = (I - n * n.transpose()) * (I - dxw2dot_dxwdot);
    dtdot_dq2.leftCols(3) = (I - n * n.transpose()) * (-dxw2dot_dq2.leftCols(3));
    dtdot_dq2.rightCols(3) = -dxw2dot_dq2.rightCols(3) - n.transpose() * vw * dn_dq2.rightCols(3) - n * (vw.transpose() * dn_dq2.rightCols(3) - n.transpose() * dxw2dot_dq2.rightCols(3));;
    dtdot_dphi2 = (I - n * n.transpose()) * (- dxw2dot_dphi2);
}

void BodySphere::test_collision_derivatives() {
    srand(1000);
    // srand(time(0));
    dtype eps = 1e-7;
    // generate random xw, xw_dot, E_2, phi_2
    Vector3 xw1 = Vector3::Random() * 10.;
    Vector3 xw1_dot = Vector3::Random() * 10.;
    Eigen::Quaternion<dtype> quat_2(Vector4::Random());
    quat_2.normalize();
    Matrix4 E2 = Matrix4::Identity();
    E2.topLeftCorner(3, 3) = quat_2.toRotationMatrix();
    E2.topRightCorner(3, 1) = Vector3::Random() * 10.;
    Vector6 phi2 = Vector6::Random() * 10.;
    _radius = 2.;

    this->_E_0i = E2;
    this->_phi = phi2;
    dtype d, ddot;
    Vector3 n, tdot, xi2;
    RowVector3 dd_dxw1, dddot_dxw1, dddot_dxw1dot;
    RowVector6 dd_dq2, dddot_dq2, dddot_dphi2;
    Matrix3 dn_dxw1, dtdot_dxw1, dtdot_dxw1dot, dxi2_dxw1;
    Matrix36 dn_dq2, dtdot_dq2, dtdot_dphi2, dxi2_dq2;
    collision(xw1, xw1_dot, d, n, ddot, tdot, xi2,
                dd_dxw1, dd_dq2,
                dn_dxw1, dn_dq2,
                dddot_dxw1, dddot_dxw1dot,
                dddot_dq2, dddot_dphi2,
                dtdot_dxw1, dtdot_dxw1dot,
                dtdot_dq2, dtdot_dphi2,
                dxi2_dxw1, dxi2_dq2);

    // test time derivatives
    {
        dtype d_ori, ddot_ori;
        Vector3 n_ori, tdot_ori, xi2_ori;
        collision(xw1, xw1_dot, d_ori, n_ori, ddot_ori, tdot_ori, xi2_ori);

        dtype ddot_fd;

        Vector3 xw1_pos = xw1 + eps * xw1_dot;
        Vector6 dq = eps * phi2;
        Matrix4 E2_pos = E2 * math::exp(dq);
        this->_E_0i = E2_pos;
        dtype d_pos, ddot_pos;
        Vector3 n_pos, tdot_pos, xi2_pos;
        collision(xw1_pos, xw1_dot, d_pos, n_pos, ddot_pos, tdot_pos, xi2_pos);
        ddot_fd = (d_pos - d) / eps;
        this->_E_0i = E2;

        print_error("Body Sphere Collision Derivatives: ddot", ddot_ori, ddot_fd);
    }

    // test dxw1 related
    RowVector3 dd_dxw1_fd, dddot_dxw1_fd;
    Matrix3 dn_dxw1_fd, dtdot_dxw1_fd, dxi2_dxw1_fd;
    for (int i = 0;i < 3;i++) {
        Vector3 xw1_pos = xw1;
        xw1_pos[i] += eps;
        dtype d_pos, ddot_pos;
        Vector3 n_pos, tdot_pos, xi2_pos;
        collision(xw1_pos, xw1_dot, d_pos, n_pos, ddot_pos, tdot_pos, xi2_pos);
        dd_dxw1_fd[i] = (d_pos - d) / eps;
        dddot_dxw1_fd[i] = (ddot_pos - ddot) / eps;
        dn_dxw1_fd.col(i) = (n_pos - n) / eps;
        dtdot_dxw1_fd.col(i) = (tdot_pos - tdot) / eps;
        dxi2_dxw1_fd.col(i) = (xi2_pos - xi2) / eps;
    }
    print_error("Body Sphere Collision Derivatives: dd_dxw1", dd_dxw1, dd_dxw1_fd);
    print_error("Body Sphere Collision Derivatives: dddot_dxw1", dddot_dxw1, dddot_dxw1_fd);
    print_error("Body Sphere Collision Derivatives: dn_dxw1", dn_dxw1, dn_dxw1_fd);
    print_error("Body Sphere Collision Derivatives: dtdot_dxw1", dtdot_dxw1, dtdot_dxw1_fd);
    print_error("Body Sphere Collision Derivatives: dxi2_dxw1", dxi2_dxw1, dxi2_dxw1_fd);

    // test dxw1dot related
    RowVector3 dddot_dxw1dot_fd;
    Matrix3 dtdot_dxw1dot_fd;
    for (int i = 0;i < 3;i++) {
        Vector3 xw1dot_pos = xw1_dot;
        xw1dot_pos[i] += eps;
        dtype d_pos, ddot_pos;
        Vector3 n_pos, tdot_pos, xi2_pos;
        collision(xw1, xw1dot_pos, d_pos, n_pos, ddot_pos, tdot_pos, xi2_pos);
        dddot_dxw1dot_fd[i] = (ddot_pos - ddot) / eps;
        dtdot_dxw1dot_fd.col(i) = (tdot_pos - tdot) / eps;
    }
    print_error("Body Sphere Collision Derivatives: dddot_dxw1dot", dddot_dxw1dot, dddot_dxw1dot_fd);
    print_error("Body Sphere Collision Derivatives: dtdot_dxw1dot", dtdot_dxw1dot, dtdot_dxw1dot_fd);

    // test dq2 related
    RowVector6 dd_dq2_fd, dddot_dq2_fd;
    Matrix36 dn_dq2_fd, dtdot_dq2_fd, dxi2_dq2_fd;
    for (int i = 0;i < 6;i++) {
        Vector6 dq = Vector6::Zero();
        dq[i] = eps;
        Matrix4 E2_pos = E2 * math::exp(dq);
        dtype d_pos, ddot_pos;
        Vector3 n_pos, tdot_pos, xi2_pos;
        this->_E_0i = E2_pos;
        collision(xw1, xw1_dot, d_pos, n_pos, ddot_pos, tdot_pos, xi2_pos);
        dd_dq2_fd[i] = (d_pos - d) / eps;
        dddot_dq2_fd[i] = (ddot_pos - ddot) / eps;
        dn_dq2_fd.col(i) = (n_pos - n) / eps;
        dtdot_dq2_fd.col(i) = (tdot_pos - tdot) / eps;
        dxi2_dq2_fd.col(i) = (xi2_pos - xi2) / eps;
    }
    print_error("Body Sphere Collision Derivatives: dd_dq2", dd_dq2, dd_dq2_fd);
    print_error("Body Sphere Collision Derivatives: dddot_dq2", dddot_dq2, dddot_dq2_fd);
    print_error("Body Sphere Collision Derivatives: dn_dw2", dn_dq2.leftCols(3), dn_dq2_fd.leftCols(3));
    print_error("Body Sphere Collision Derivatives: dn_dv2", dn_dq2.rightCols(3), dn_dq2_fd.rightCols(3));
    print_error("Body Sphere Collision Derivatives: dtdot_dw2", dtdot_dq2.leftCols(3), dtdot_dq2_fd.leftCols(3));
    print_error("Body Sphere Collision Derivatives: dtdot_dv2", dtdot_dq2.rightCols(3), dtdot_dq2_fd.rightCols(3));
    print_error("Body Sphere Collision Derivatives: dxi2_dq2", dxi2_dq2, dxi2_dq2_fd);
    this->_E_0i = E2;

    // test dphi2 related
    RowVector6 dddot_dphi2_fd;
    Matrix36 dtdot_dphi2_fd;
    for (int i = 0;i < 6;i++) {
        Vector6 phi2_pos = phi2;
        phi2_pos[i] += eps;
        this->_phi = phi2_pos;
        dtype d_pos, ddot_pos;
        Vector3 n_pos, tdot_pos, xi2_pos;
        collision(xw1, xw1_dot, d_pos, n_pos, ddot_pos, tdot_pos, xi2_pos);
        dddot_dphi2_fd[i] = (ddot_pos - ddot) / eps;
        dtdot_dphi2_fd.col(i) = (tdot_pos - tdot) / eps;
    }
    print_error("Body Sphere Collision Derivatives: dddot_dw2dot", dddot_dphi2.head(3), dddot_dphi2_fd.head(3));
    print_error("Body Sphere Collision Derivatives: dddot_dv2dot", dddot_dphi2.tail(3), dddot_dphi2_fd.tail(3));
    print_error("Body Sphere Collision Derivatives: dtdot_dw2dot", dtdot_dphi2.leftCols(3), dtdot_dphi2_fd.leftCols(3));
    print_error("Body Sphere Collision Derivatives: dtdot_dv2dot", dtdot_dphi2.rightCols(3), dtdot_dphi2_fd.rightCols(3));
    this->_phi = phi2;
}
}