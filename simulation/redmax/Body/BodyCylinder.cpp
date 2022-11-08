#include "BodyCylinder.h"
#include "Simulation.h"
#include "tiny_obj_loader.h"

namespace redmax {

BodyCylinder::BodyCylinder(Simulation *sim, Joint *joint, dtype height, dtype radius, 
            Matrix3 R_ji, Vector3 p_ji, dtype density) 
            : BodyPrimitiveShape(sim, joint, R_ji, p_ji, density) {
    _height = height;
    _radius = radius;

    load_mesh();
    _V.topRows(2) *= (float)_radius / 0.5f;
    _V.bottomRows(1) *= (float)_height;
    
    precompute_contact_points();
    computeMassMatrix();
}

void BodyCylinder::load_mesh() {
    std::string filename = std::string(GRAPHICS_CODEBASE_SOURCE_DIR) + "/resources/meshes/cylinder.obj";
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
            attrib.vertices[i * 3 + 2]);
    }
    
    int num_elements = (int)obj_shape[0].mesh.indices.size() / 3;
    _F.resize(3, num_elements);
    for (int i = 0;i < num_elements;i++) {
        _F.col(i) = Vector3i(obj_shape[0].mesh.indices[i * 3].vertex_index,
            obj_shape[0].mesh.indices[i * 3 + 1].vertex_index,
            obj_shape[0].mesh.indices[i * 3 + 2].vertex_index);
    }
}

void BodyCylinder::computeMassMatrix() {
    _mass = constants::pi * _radius * _radius * _height * _density;
    _Inertia(0) = _mass / 4.0 * (_radius * _radius) + _mass / 12.0 * (_height * _height);
    _Inertia(1) = _mass / 4.0 * (_radius * _radius) + _mass / 12.0 * (_height * _height);
    _Inertia(2) = _mass / 2.0 * (_radius * _radius);
    _Inertia.tail(3) = Vector3::Constant(_mass);
}

void BodyCylinder::precompute_contact_points() {
    srand(1000);
    _contact_points.clear();
    int sample_rate = 1;

    for (int i = 0;i < _V.cols();i++) {
        int p = rand() % sample_rate;
        if (p == 0) {
            _contact_points.push_back(_V.col(i));
        }
    }
}

Matrix3X BodyCylinder::get_vertices() const {
    return _V;
}

Matrix3Xi BodyCylinder::get_faces() const {
    return _F;
}

// rendering
void BodyCylinder::get_rendering_objects(
    std::vector<Matrix3Xf>& vertex_list, 
    std::vector<Matrix3Xi>& face_list,
    std::vector<opengl_viewer::Option>& option_list,
    std::vector<opengl_viewer::Animator*>& animator_list) {

    Matrix3Xf vertex = _V.cast<float>();

    if (_sim->_options->_unit == "cm-g") 
        vertex /= 10.;
    else
        vertex *= 10.;
        
    opengl_viewer::Option object_option;

    object_option.SetBoolOption("smooth normal", false);
    // object_option.SetVectorOption("ambient", 0.25f, 0.148f, 0.06475f);
    object_option.SetVectorOption("ambient", _color(0), _color(1), _color(2));
    object_option.SetVectorOption("diffuse", 0.4f, 0.2368f, 0.1036f);
    object_option.SetVectorOption("specular", 0.774597f, 0.658561f, 0.400621f);
    object_option.SetFloatOption("shininess", 76.8f);

    _animator = new BodyAnimator(this);

    vertex_list.push_back(vertex);
    face_list.push_back(_F);
    option_list.push_back(object_option);
    animator_list.push_back(_animator);
}

dtype BodyCylinder::distance(Vector3 xw) {
    Vector3 x = _E_0i.topLeftCorner(3, 3).transpose() * (xw - _E_0i.topRightCorner(3, 1));
    dtype dxy = std::sqrt(x[0] * x[0] + x[1] * x[1]) - _radius;
    dtype dz = max(x[2] - _height / 2., -x[2] - _height / 2.);
    dtype d = max(dxy, dz);
    return d;
}

void BodyCylinder::collision(
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
    Vector3 x = R2.transpose() * (xw - p2);

    dtype dxy = std::sqrt(x[0] * x[0] + x[1] * x[1]) - _radius;
    dtype dz = max(x[2] - _height / 2., -x[2] - _height / 2.);
    Vector3 e;
    if (dxy > dz) {
        d = dxy;
        e = Vector3(x[0], x[1], 0);
        e /= e.norm();
    } else {
        d = dz;
        if (x[2] > 0) e = Vector3::Unit(2);
        else e = -Vector3::Unit(2);
    }
    n = R2 * e;
    xi2 = x - d * e;
    Vector3 xdot = math::skew(w2_dot).transpose() * x + R2.transpose() * xw_dot - v2_dot;
    ddot = e.transpose() * xdot;
    Vector3 xw2_dot = R2 * (w2_dot.cross(xi2) + v2_dot);
    Vector3 vw = xw_dot - xw2_dot;
    tdot = (I - n * n.transpose()) * vw;
}

// refer to https://www.overleaf.com/project/5cd0a78e0e41fe0f8632a42a
// section 1.4: General Penalty-based Contact
void BodyCylinder::collision(
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
    Vector3 x = R2.transpose() * (xw - p2);

    /**************** values ****************/
    dtype dxy = std::sqrt(x[0] * x[0] + x[1] * x[1]) - _radius;
    dtype dz = max(x[2] - _height / 2., -x[2] - _height / 2.);
    Vector3 e;
    if (dxy > dz) {
        d = dxy;
        e = Vector3(x[0], x[1], 0);
        e /= e.norm();
    } else {
        d = dz;
        if (x[2] > 0) e = Vector3::Unit(2);
        else e = -Vector3::Unit(2);
    }
    n = R2 * e;
    xi2 = x - d * e;
    Vector3 xdot = math::skew(w2_dot).transpose() * x + R2.transpose() * xw_dot - v2_dot;
    ddot = e.transpose() * xdot;
    Vector3 xw2_dot = R2 * (w2_dot.cross(xi2) + v2_dot);
    Vector3 vw = xw_dot - xw2_dot;
    tdot = (I - n * n.transpose()) * vw;

    /**************** derivatives ****************/

    Vector3 e1 = Vector3::UnitX(), e2 = Vector3::UnitY(), e3 = Vector3::UnitZ();
    
    // x
    Matrix3 dx_dxw = R2.transpose();
    Matrix3 dx_dw2 = math::skew(x);
    Matrix3 dx_dv2 = -I;

    // d
    Vector3 dd_dx = e;
    dd_dxw = dd_dx.transpose() * dx_dxw;
    dd_dq2.head(3) = dd_dx.transpose() * dx_dw2;
    dd_dq2.tail(3) = -dd_dx.transpose();

    // n
    Matrix3 de_dx = Matrix3::Zero();
    if (dxy > dz) {
        de_dx.col(0) = Vector3(x[1] * x[1], -x[0] * x[1], 0);
        de_dx.col(1) = Vector3(-x[0] * x[1], x[0] * x[0], 0);
        de_dx /= std::pow(std::sqrt(x[0] * x[0] + x[1] * x[1]), 3);
    }
    Matrix3 de_dxw = de_dx * dx_dxw;
    Matrix3 de_dw2 = de_dx * dx_dw2;
    Matrix3 de_dv2 = de_dx * dx_dv2;

    dn_dxw = R2 * de_dxw;
    dn_dq2.leftCols(3) = R2 * de_dw2;
    dn_dq2.col(0) += R2 * math::skew(e1) * e;
    dn_dq2.col(1) += R2 * math::skew(e2) * e;
    dn_dq2.col(2) += R2 * math::skew(e3) * e;
    dn_dq2.rightCols(3) = R2 * de_dv2;

    // xi2
    dxi2_dxw = dx_dxw - e * dd_dxw - d * de_dxw;
    dxi2_dq2.leftCols(3) = dx_dw2 - e * dd_dq2.head(3) - d * de_dw2;
    dxi2_dq2.rightCols(3) = dx_dv2 - e * dd_dq2.tail(3) - d * de_dv2;

    // ddot
    Matrix3 ddd_dx_dx = de_dx;
    dddot_dxw = (dd_dx.transpose() * math::skew(w2_dot).transpose() + xdot.transpose() * ddd_dx_dx) * R2.transpose();
    dddot_dxwdot = dd_dx.transpose() * R2.transpose();
    dddot_dq2.leftCols(3) = dd_dx.transpose() * (-math::skew(w2_dot) * math::skew(x) + math::skew(R2.transpose() * xw_dot)) + xdot.transpose() * ddd_dx_dx * dx_dw2;
    dddot_dq2.rightCols(3) = dd_dx.transpose() * math::skew(w2_dot) + xdot.transpose() * ddd_dx_dx * dx_dv2;
    dddot_dphi2.leftCols(3) = dd_dx.transpose() * math::skew(x);
    dddot_dphi2.rightCols(3) = dd_dx.transpose() * (-I);

    // xw2_dot
    Matrix3 dxw2dot_dxw = R2 * math::skew(w2_dot) * dxi2_dxw;
    Matrix3 dxw2dot_dxwdot = Matrix3::Zero();
    Matrix36 dxw2dot_dq2;
    dxw2dot_dq2.leftCols(3) = -R2 * math::skew(w2_dot.cross(xi2) + v2_dot) + R2 * math::skew(w2_dot) * dxi2_dq2.leftCols(3);
    dxw2dot_dq2.rightCols(3) = R2 * math::skew(w2_dot) * dxi2_dq2.rightCols(3);
    Matrix36 dxw2dot_dphi2;    
    dxw2dot_dphi2.col(0) = R2 * math::skew(e1) * xi2;
    dxw2dot_dphi2.col(1) = R2 * math::skew(e2) * xi2;
    dxw2dot_dphi2.col(2) = R2 * math::skew(e3) * xi2;
    dxw2dot_dphi2.rightCols(3) = R2;

    // tdot
    dtdot_dxw = -dxw2dot_dxw - n.transpose() * vw * dn_dxw - n * vw.transpose() * dn_dxw + n * n.transpose() * dxw2dot_dxw;
    dtdot_dxwdot = (I - n * n.transpose()) * (I - dxw2dot_dxwdot);
    dtdot_dq2.leftCols(3) = -dxw2dot_dq2.leftCols(3) - n.transpose() * vw * dn_dq2.leftCols(3) - n * (vw.transpose() * dn_dq2.leftCols(3) - n.transpose() * dxw2dot_dq2.leftCols(3));
    dtdot_dq2.rightCols(3) = -dxw2dot_dq2.rightCols(3) - n.transpose() * vw * dn_dq2.rightCols(3) - n * (vw.transpose() * dn_dq2.rightCols(3) - n.transpose() * dxw2dot_dq2.rightCols(3));
    dtdot_dphi2 = (I - n * n.transpose()) * (- dxw2dot_dphi2);
}

void BodyCylinder::test_collision_derivatives() {
    // srand(1000);
    srand(time(0));
    dtype eps = 1e-7;
    // generate random xw, xw_dot, E_2, phi_2
    Eigen::Quaternion<dtype> quat_2(Vector4::Random());
    quat_2.normalize();
    Matrix4 E2 = Matrix4::Identity();
    E2.topLeftCorner(3, 3) = quat_2.toRotationMatrix();
    E2.topRightCorner(3, 1) = Vector3::Random();
    Vector6 phi2 = Vector6::Random();

    _radius = 0.5;
    _height = 1.0;

    Vector3 x = Vector3::Random();
    Vector3 xw1 = E2.topLeftCorner(3, 3) * x + E2.topRightCorner(3, 1);
    Vector3 xw1_dot = Vector3::Random();

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

void BodyCylinder::test_collision_derivatives_runtime(Vector3 xw, Vector3 xw_dot) {
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

        print_error("Body Cylinder Collision Derivatives: ddot", ddot, ddot_fd);
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
    print_error("Body Cylinder Collision Derivatives: dd_dxw1", dd_dxw1, dd_dxw1_fd);
    print_error("Body Cylinder Collision Derivatives: dddot_dxw1", dddot_dxw1, dddot_dxw1_fd);
    print_error("Body Cylinder Collision Derivatives: dn_dxw1", dn_dxw1, dn_dxw1_fd);
    print_error("Body Cylinder Collision Derivatives: dtdot_dxw1", dtdot_dxw1, dtdot_dxw1_fd);
    print_error("Body Cylinder Collision Derivatives: dxi2_dxw1", dxi2_dxw1, dxi2_dxw1_fd);

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
    print_error("Body Cylinder Collision Derivatives: dddot_dxw1dot", dddot_dxw1dot, dddot_dxw1dot_fd);
    print_error("Body Cylinder Collision Derivatives: dtdot_dxw1dot", dtdot_dxw1dot, dtdot_dxw1dot_fd);

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
    print_error("Body Cylinder Collision Derivatives: dd_dq2", dd_dq2, dd_dq2_fd);
    print_error("Body Cylinder Collision Derivatives: dddot_dq2", dddot_dq2, dddot_dq2_fd);
    print_error("Body Cylinder Collision Derivatives: dn_dw2", dn_dq2.leftCols(3), dn_dq2_fd.leftCols(3));
    print_error("Body Cylinder Collision Derivatives: dn_dv2", dn_dq2.rightCols(3), dn_dq2_fd.rightCols(3));
    print_error("Body Cylinder Collision Derivatives: dtdot_dw2", dtdot_dq2.leftCols(3), dtdot_dq2_fd.leftCols(3));
    print_error("Body Cylinder Collision Derivatives: dtdot_dv2", dtdot_dq2.rightCols(3), dtdot_dq2_fd.rightCols(3));
    print_error("Body Cylinder Collision Derivatives: dxi2_dq2", dxi2_dq2, dxi2_dq2_fd);
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
    print_error("Body Cylinder Collision Derivatives: dddot_dw2dot", dddot_dphi2.head(3), dddot_dphi2_fd.head(3));
    print_error("Body Cylinder Collision Derivatives: dddot_dv2dot", dddot_dphi2.tail(3), dddot_dphi2_fd.tail(3));
    print_error("Body Cylinder Collision Derivatives: dtdot_dw2dot", dtdot_dphi2.leftCols(3), dtdot_dphi2_fd.leftCols(3));
    print_error("Body Cylinder Collision Derivatives: dtdot_dv2dot", dtdot_dphi2.rightCols(3), dtdot_dphi2_fd.rightCols(3));
    this->_phi = phi2;
}

}