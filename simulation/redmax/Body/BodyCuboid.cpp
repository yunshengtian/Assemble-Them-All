#include "BodyCuboid.h"
#include "Simulation.h"

namespace redmax {

BodyCuboid::BodyCuboid(Simulation *sim, Joint *joint, Vector3 length, 
            Matrix3 R_ji, Vector3 p_ji, dtype density) 
            : BodyPrimitiveShape(sim, joint, R_ji, p_ji, density) {
    _length = length;
    precompute_contact_points();
    computeMassMatrix();
}

void BodyCuboid::computeMassMatrix() {
    _mass = _length.prod() * _density;
    _Inertia(0) = _mass / 12.0 * (_length(1) * _length(1) + _length(2) * _length(2));
    _Inertia(1) = _mass / 12.0 * (_length(0) * _length(0) + _length(2) * _length(2));
    _Inertia(2) = _mass / 12.0 * (_length(0) * _length(0) + _length(1) * _length(1));
    _Inertia.tail(3).noalias() = Vector3::Constant(_mass);
}

void BodyCuboid::precompute_contact_points() {
    _contact_points.clear();
    for (int i = -1;i <= 1;i += 2)
        for (int j = -1;j <= 1;j += 2) 
            for (int k = -1;k <= 1;k += 2)
                _contact_points.push_back(Vector3(i * _length[0] * 0.5, j * _length[1] * 0.5, k * _length[2] * 0.5));
}

// std::vector<Vector3> BodyCuboid::get_general_contact_points() {
//     std::vector<Vector3> contact_points;
//     contact_points.clear();
//     for (int i = 0;i < 4;i++)
//         for (int j = 0;j < 4;j++)
//             for (int k = 0;k < 4;k++) {
//                 Vector3 pos;
//                 pos(0) = i / 3 * _length[0] - _length[0] / 2.;
//                 pos(1) = j / 3 * _length[1] - _length[1] / 2.;
//                 pos(2) = k / 3 * _length[2] - _length[2] / 2.;
//                 contact_points.push_back(pos);
//             }
//     return contact_points;
// }

// std::vector<Vector3> BodyCuboid::get_general_contact_points() {
//     std::vector<Vector3> contact_points;
//     contact_points.clear();
//     for (int i = 0;i < 4;i++)
//         for (int j = 0;j < 4;j++)
//             for (int k = 0;k < 4;k++) {
//                 Vector3 pos;
//                 pos(0) = (double)i / 3 * _length[0] - _length[0] / 2.;
//                 pos(1) = (double)j / 3 * _length[1] - _length[1] / 2.;
//                 pos(2) = (double)k / 3 * _length[2] - _length[2] / 2.;
//                 contact_points.push_back(pos);
//             }
//     return contact_points;
// }

std::vector<Vector3> BodyCuboid::get_general_contact_points() {
    std::vector<Vector3> contact_points;
    contact_points.clear();
    for (int i = 0;i < 6;i++)
        for (int j = 0;j < 6;j++)
            for (int k = 0;k < 6;k++) {
                Vector3 pos;
                pos(0) = (double)i / 5 * _length[0] - _length[0] / 2.;
                pos(1) = (double)j / 5 * _length[1] - _length[1] / 2.;
                pos(2) = (double)k / 5 * _length[2] - _length[2] / 2.;
                contact_points.push_back(pos);
            }
    return contact_points;
}

Matrix3X BodyCuboid::get_vertices() const {
    Matrix3X vertices(3, 8);
    vertices << 0, 0, 0, 1, 0, 1, 1, 1,
                0, 0, 1, 0, 1, 0, 1, 1,
                0, 1, 0, 0, 1, 1, 0, 1;
    vertices.array() -= 0.5;
    for (int i = 0; i < vertices.cols(); ++i) {
        vertices.col(i) = (vertices.col(i).cwiseProduct(_length)).eval();
    }
    return vertices;
}

Matrix3Xi BodyCuboid::get_faces() const {
    Matrix3Xi faces(3, 12);
    faces << 1, 4, 1, 2, 1, 3, 2, 6, 4, 7, 3, 5,
             3, 3, 4, 4, 2, 2, 6, 8, 7, 8, 5, 8,
             4, 7, 2, 6, 3, 5, 5, 5, 6, 6, 7, 7;
    faces.array() -= 1;
    return faces;
}

// // rendering
// void BodyCuboid::get_rendering_objects(
//     std::vector<Matrix3Xf>& vertex_list, 
//     std::vector<Matrix3Xi>& face_list,
//     std::vector<opengl_viewer::Option>& option_list,
//     std::vector<opengl_viewer::Animator*>& animator_list) {

//     Matrix3Xf cube_vertex;
//     Matrix3Xi cube_face;
//     Matrix2Xf cube_uv;
//     opengl_viewer::Option object_option;

//     opengl_viewer::ReadFromObjFile(
//         std::string(GRAPHICS_CODEBASE_SOURCE_DIR) + "/resources/meshes/cube.obj",
//         cube_vertex, cube_face, cube_uv);

//     for (int i = 0;i < 3;i++)
//         cube_vertex.row(i) *= (float)_length(i) / 10.;

//     // object_option.SetBoolOption("smooth normal", false);
//     // object_option.SetVectorOption("ambient", 0.2f, 0.2f, 0.2f);
//     // object_option.SetVectorOption("diffuse", 0.5f, 0.5f, 0.5f);
//     // object_option.SetVectorOption("specular", 1.f, 1.f, 1.f);
//     // object_option.SetFloatOption("shininess", 2.0f);
//     object_option.SetBoolOption("smooth normal", false);
//     // object_option.SetVectorOption("ambient", 0.25f, 0.148f, 0.06475f);
//     object_option.SetVectorOption("ambient", _color(0), _color(1), _color(2));
//     object_option.SetVectorOption("diffuse", 0.4f, 0.2368f, 0.1036f);
//     object_option.SetVectorOption("specular", 0.774597f, 0.458561f, 0.200621f);
//     object_option.SetFloatOption("shininess", 76.8f);

//     _animator = new BodyAnimator(this);

//     vertex_list.push_back(cube_vertex);
//     face_list.push_back(cube_face);
//     option_list.push_back(object_option);
//     animator_list.push_back(_animator);
// }

// rendering
void BodyCuboid::get_rendering_objects(
    std::vector<Matrix3Xf>& vertex_list, 
    std::vector<Matrix3Xi>& face_list,
    std::vector<opengl_viewer::Option>& option_list,
    std::vector<opengl_viewer::Animator*>& animator_list) {

    Matrix3Xf cube_vertex;
    Matrix3Xi cube_face;
    Matrix2Xf cube_uv;
    opengl_viewer::Option object_option;

    opengl_viewer::ReadFromObjFile(
        std::string(GRAPHICS_CODEBASE_SOURCE_DIR) + "/resources/meshes/cube.obj",
        cube_vertex, cube_face, cube_uv);

    // set cube uv
    Eigen::Vector2f texture_bottom_left;
    Eigen::Vector2f texture_top_right;
    Eigen::Vector2f texture_pos;
    for (int i = 0;i < cube_vertex.cols();i++) {
        if (cube_vertex(0, i) <= -0.5 + math::eps) { // middle (1)
            texture_bottom_left = Eigen::Vector2f(0.25, 2. / 3.);
            texture_top_right = Eigen::Vector2f(0., 1. / 3.);
            texture_pos = Eigen::Vector2f(cube_vertex(1, i), cube_vertex(2, i)) + Eigen::Vector2f(0.5, 0.5);
        } else if (cube_vertex(1, i) <= -0.5 + math::eps) { // middle (2)
            texture_bottom_left = Eigen::Vector2f(0.25, 2. / 3.);
            texture_top_right = Eigen::Vector2f(0.5, 1. / 3.);
            texture_pos = Eigen::Vector2f(cube_vertex(0, i), cube_vertex(2, i)) + Eigen::Vector2f(0.5, 0.5);
        } else if (cube_vertex(0, i) >= 0.5 - math::eps) { // middle (3)
            texture_bottom_left = Eigen::Vector2f(0.5, 2. / 3.);
            texture_top_right = Eigen::Vector2f(0.75, 1. / 3.);
            texture_pos = Eigen::Vector2f(cube_vertex(1, i), cube_vertex(2, i)) + Eigen::Vector2f(0.5, 0.5);
        } else if (cube_vertex(1, i) >= 0.5 - math::eps) { // middle (4)
            texture_bottom_left = Eigen::Vector2f(1., 2. / 3.);
            texture_top_right = Eigen::Vector2f(0.75, 1. / 3.);
            texture_pos = Eigen::Vector2f(cube_vertex(0, i), cube_vertex(2, i)) + Eigen::Vector2f(0.5, 0.5);
        } else if (cube_vertex(2, i) >= 0.5 - math::eps) { // top
            texture_bottom_left = Eigen::Vector2f(0.25, 1. / 3.);
            texture_top_right = Eigen::Vector2f(0.5 - 0.01, 0. + 0.01);
            texture_pos = Eigen::Vector2f(cube_vertex(0, i), cube_vertex(1, i)) + Eigen::Vector2f(0.5, 0.5);
        } else { // bottom
            texture_bottom_left = Eigen::Vector2f(0.25, 2. / 3.);
            texture_top_right = Eigen::Vector2f(0.5 - 0.01, 1. - 0.01);
            texture_pos = Eigen::Vector2f(cube_vertex(0, i), cube_vertex(1, i)) + Eigen::Vector2f(0.5, 0.5);
        }
        cube_uv.col(i) = texture_bottom_left + (texture_top_right - texture_bottom_left).cwiseProduct(texture_pos);
    }

    for (int i = 0;i < 3;i++)
        cube_vertex.row(i) *= (float)_length(i);

    if (_sim->_options->_unit == "cm-g") 
        cube_vertex /= 10.;
    else
        cube_vertex *= 10.;

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
        
        object_option.SetMatrixOption("uv", cube_uv);
        object_option.SetMatrixOption("texture", checker_texture.rgb_data());
        object_option.SetIntOption("texture row num", checker_texture.row_num());
        object_option.SetIntOption("texture col num", checker_texture.col_num());
        object_option.SetStringOption("texture mag filter", "nearest");
    }

    _animator = new BodyAnimator(this);

    vertex_list.push_back(cube_vertex);
    face_list.push_back(cube_face);
    option_list.push_back(object_option);
    animator_list.push_back(_animator);
}

dtype BodyCuboid::distance(Vector3 xw) {
    Vector3 x = _E_0i.topLeftCorner(3, 3).transpose() * (xw - _E_0i.topRightCorner(3, 1));
    Vector3 s = _length / 2.;
    dtype d = -99999999.;
    for (int i = 0;i < 3;i++) {
        d = max(d, max(x[i] - s[i], -x[i] - s[i]));
    }

    return d;
}

void BodyCuboid::collision(
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
    Vector3 s = _length / 2.;

    Vector3 x = R2.transpose() * (xw - p2);

    d = -9999999.;
    Vector3 e;
    for (int i = 0;i < 3;i++) {
        if (x[i] - s[i] > d) {
            d = x[i] - s[i];
            e = Vector3::Unit(i);
        }
        if (-x[i] - s[i] > d) {
            d = -x[i] - s[i];
            e = -Vector3::Unit(i);
        }
    }
    
    n = R2 * e;

    xi2 = x - d * e;

    ddot = e.transpose() * (math::skew(w2_dot).transpose() * x + R2.transpose() * xw_dot - v2_dot);

    Vector3 xw2_dot = R2 * (w2_dot.cross(xi2) + v2_dot);
    Vector3 vw = xw_dot - xw2_dot;
    tdot = (I - n * n.transpose()) * vw;
}

// refer to https://www.overleaf.com/project/5cd0a78e0e41fe0f8632a42a
// section 1.4: General Penalty-based Contact
void BodyCuboid::collision(
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
    Vector3 s = _length / 2.;

    /**************** values ****************/
    Vector3 x = R2.transpose() * (xw - p2);
    d = -9999999.;
    Vector3 e;
    for (int i = 0;i < 3;i++) {
        if (x[i] - s[i] > d) {
            d = x[i] - s[i];
            e = Vector3::Unit(i);
        }
        if (-x[i] - s[i] > d) {
            d = -x[i] - s[i];
            e = -Vector3::Unit(i);
        }
    }
    n = R2 * e;
    xi2 = x - d * e;
    ddot = e.transpose() * (math::skew(w2_dot).transpose() * x + R2.transpose() * xw_dot - v2_dot);
    Vector3 xw2_dot = R2 * (w2_dot.cross(xi2) + v2_dot);
    Vector3 vw = xw_dot - xw2_dot;
    tdot = (I - n * n.transpose()) * vw;

    /**************** derivatives ****************/
    
    // x
    Matrix3 dx_dxw = R2.transpose();
    Matrix3 dx_dw2 = math::skew(x);
    Matrix3 dx_dv2 = -I;

    // d
    dd_dxw = e.transpose() * dx_dxw;
    dd_dq2.head(3) = e.transpose() * dx_dw2;
    dd_dq2.tail(3) = -e.transpose();

    // n
    dn_dxw.setZero();
    dn_dq2.leftCols(3) = -R2 * math::skew(e);
    dn_dq2.rightCols(3).setZero();

    // xi2
    dxi2_dxw = dx_dxw - e * dd_dxw;
    dxi2_dq2.leftCols(3) = dx_dw2 - e * dd_dq2.leftCols(3);
    dxi2_dq2.rightCols(3) = -I + e * e.transpose();

    // ddot
    dddot_dxw = e.transpose() * math::skew(w2_dot).transpose() * R2.transpose();
    dddot_dxwdot = e.transpose() * R2.transpose();
    dddot_dq2.leftCols(3) = e.transpose() * (math::skew(w2_dot).transpose() * math::skew(x) 
                + math::skew(R2.transpose() * xw_dot));
    dddot_dq2.rightCols(3) = -e.transpose() * math::skew(w2_dot).transpose();
    dddot_dphi2.col(0) = e.transpose() * math::skew(Vector3::UnitX()).transpose() * x;
    dddot_dphi2.col(1) = e.transpose() * math::skew(Vector3::UnitY()).transpose() * x;
    dddot_dphi2.col(2) = e.transpose() * math::skew(Vector3::UnitZ()).transpose() * x;
    dddot_dphi2.rightCols(3) = -e.transpose();

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
    dtdot_dxw = -(I - n * n.transpose()) * dxw2dot_dxw;
    dtdot_dxwdot = (I - n * n.transpose()) * (I - dxw2dot_dxwdot);
    dtdot_dq2.leftCols(3) = -dxw2dot_dq2.leftCols(3) - n.transpose() * vw * dn_dq2.leftCols(3) - n * (vw.transpose() * dn_dq2.leftCols(3) - n.transpose() * dxw2dot_dq2.leftCols(3));
    dtdot_dq2.rightCols(3) = (I - n * n.transpose()) * (-dxw2dot_dq2.rightCols(3));
    dtdot_dphi2 = (I - n * n.transpose()) * (- dxw2dot_dphi2);
}

void BodyCuboid::test_collision_derivatives() {
    // srand(1000);
    srand(time(0));
    dtype eps = 1e-7;
    // generate random xw, xw_dot, E_2, phi_2
    Eigen::Quaternion<dtype> quat_2(Vector4::Random());
    quat_2.normalize();
    Matrix4 E2 = Matrix4::Identity();
    E2.topLeftCorner(3, 3) = quat_2.toRotationMatrix();
    E2.topRightCorner(3, 1) = Vector3::Random() * 10.;
    Vector6 phi2 = Vector6::Random() * 10.;
    _length = Vector3(1.0, 2.0, 0.5);
    Vector3 x = Vector3::Random().cwiseProduct(_length / 2.);
    Vector3 xw1 = E2.topLeftCorner(3, 3) * x + E2.topRightCorner(3, 1);
    Vector3 xw1_dot = Vector3::Random() * 10.;

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

        print_error("ddot", ddot_ori, ddot_fd);
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
    print_error("dd_dxw1", dd_dxw1, dd_dxw1_fd);
    print_error("dddot_dxw1", dddot_dxw1, dddot_dxw1_fd);
    print_error("dn_dxw1", dn_dxw1, dn_dxw1_fd);
    print_error("dtdot_dxw1", dtdot_dxw1, dtdot_dxw1_fd);
    print_error("dxi2_dxw1", dxi2_dxw1, dxi2_dxw1_fd);

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
    print_error("dddot_dxw1dot", dddot_dxw1dot, dddot_dxw1dot_fd);
    print_error("dtdot_dxw1dot", dtdot_dxw1dot, dtdot_dxw1dot_fd);

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
    print_error("dd_dq2", dd_dq2, dd_dq2_fd);
    print_error("dddot_dq2", dddot_dq2, dddot_dq2_fd);
    print_error("dn_dw2", dn_dq2.leftCols(3), dn_dq2_fd.leftCols(3));
    print_error("dn_dv2", dn_dq2.rightCols(3), dn_dq2_fd.rightCols(3));
    print_error("dtdot_dw2", dtdot_dq2.leftCols(3), dtdot_dq2_fd.leftCols(3));
    print_error("dtdot_dv2", dtdot_dq2.rightCols(3), dtdot_dq2_fd.rightCols(3));
    print_error("dxi2_dq2", dxi2_dq2, dxi2_dq2_fd);
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
    print_error("dddot_dw2dot", dddot_dphi2.head(3), dddot_dphi2_fd.head(3));
    print_error("dddot_dv2dot", dddot_dphi2.tail(3), dddot_dphi2_fd.tail(3));
    print_error("dtdot_dw2dot", dtdot_dphi2.leftCols(3), dtdot_dphi2_fd.leftCols(3));
    print_error("dtdot_dv2dot", dtdot_dphi2.rightCols(3), dtdot_dphi2_fd.rightCols(3));
    this->_phi = phi2;
}

void BodyCuboid::test_collision_derivatives_runtime(Vector3 xw, Vector3 xw_dot) {
    dtype eps = 1e-7;

    Matrix4 E2 = this->_E_0i;
    Vector6 phi2 = this->_phi;
    Vector3 xw1 = xw;
    Vector3 xw1_dot = xw_dot;
    Vector3 x = E2.topLeftCorner(3, 3).transpose() * (xw1 - E2.topRightCorner(3, 1));

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
        // dtype d_ori, ddot_ori;
        // Vector3 n_ori, tdot_ori, xi2_ori;
        // collision(xw1, xw1_dot, d_ori, n_ori, ddot_ori, tdot_ori, xi2_ori);

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

        print_error("Body Cuboid Collision Derivatives: ddot", ddot, ddot_fd);
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
    print_error("Body Cuboid Collision Derivatives: dd_dxw1", dd_dxw1, dd_dxw1_fd);
    print_error("Body Cuboid Collision Derivatives: dddot_dxw1", dddot_dxw1, dddot_dxw1_fd);
    print_error("Body Cuboid Collision Derivatives: dn_dxw1", dn_dxw1, dn_dxw1_fd);
    print_error("Body Cuboid Collision Derivatives: dtdot_dxw1", dtdot_dxw1, dtdot_dxw1_fd);
    print_error("Body Cuboid Collision Derivatives: dxi2_dxw1", dxi2_dxw1, dxi2_dxw1_fd);

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
    print_error("Body Cuboid Collision Derivatives: dddot_dxw1dot", dddot_dxw1dot, dddot_dxw1dot_fd);
    print_error("Body Cuboid Collision Derivatives: dtdot_dxw1dot", dtdot_dxw1dot, dtdot_dxw1dot_fd);

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
    print_error("Body Cuboid Collision Derivatives: dd_dq2", dd_dq2, dd_dq2_fd);
    print_error("Body Cuboid Collision Derivatives: dddot_dq2", dddot_dq2, dddot_dq2_fd);
    print_error("Body Cuboid Collision Derivatives: dn_dw2", dn_dq2.leftCols(3), dn_dq2_fd.leftCols(3));
    print_error("Body Cuboid Collision Derivatives: dn_dv2", dn_dq2.rightCols(3), dn_dq2_fd.rightCols(3));
    print_error("Body Cuboid Collision Derivatives: dtdot_dw2", dtdot_dq2.leftCols(3), dtdot_dq2_fd.leftCols(3));
    print_error("Body Cuboid Collision Derivatives: dtdot_dv2", dtdot_dq2.rightCols(3), dtdot_dq2_fd.rightCols(3));
    print_error("Body Cuboid Collision Derivatives: dxi2_dq2", dxi2_dq2, dxi2_dq2_fd);
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
    print_error("Body Cuboid Collision Derivatives: dddot_dw2dot", dddot_dphi2.head(3), dddot_dphi2_fd.head(3));
    print_error("Body Cuboid Collision Derivatives: dddot_dv2dot", dddot_dphi2.tail(3), dddot_dphi2_fd.tail(3));
    print_error("Body Cuboid Collision Derivatives: dtdot_dw2dot", dtdot_dphi2.leftCols(3), dtdot_dphi2_fd.leftCols(3));
    print_error("Body Cuboid Collision Derivatives: dtdot_dv2dot", dtdot_dphi2.rightCols(3), dtdot_dphi2_fd.rightCols(3));
    this->_phi = phi2;
}

}