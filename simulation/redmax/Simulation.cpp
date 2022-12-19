#include "Simulation.h"
#include "Robot.h"
#include "Body/Body.h"
#include "Body/BodyPrimitiveShape.h"
#include "Body/BodySDFObj.h"
#include "Body/BodyBVHObj.h"
#include "Joint/Joint.h"
#include "SimViewer.h"
#include "Force/Force.h"

namespace redmax {

Simulation::~Simulation() {
}

/******************************************* Forward Dynamics ************************************************/

// init simulation
void Simulation::init(bool verbose) {
    _robot->init(verbose);

    _ndof_r = _robot->_ndof_r;
    _ndof_m = _robot->_ndof_m;
    _ndof_u = _robot->_ndof_u;
    _ndof_var = _robot->_ndof_var;
    _ndof_p = _robot->_ndof_p;
    _ndof_p1 = _robot->_ndof_p1;
    _ndof_p2 = _robot->_ndof_p2;
    _ndof_p3 = _robot->_ndof_p3;
    _ndof_p4 = _robot->_ndof_p4;
    _ndof_p5 = _robot->_ndof_p5;
    _ndof_p6 = _robot->_ndof_p6;

    _dM_dp = JacobianMatrixVector(_ndof_r, _ndof_r, _ndof_p);
    _df_dp = MatrixX::Zero(_ndof_r, _ndof_p);

    // q_init and qdot_init
    _q_init = get_q();
    _qdot_init = get_qdot();
}

// init states
void Simulation::set_state_init(const VectorX q_init, const VectorX qdot_init) {
    set_q_init(q_init);
    set_qdot_init(qdot_init);
}

void Simulation::set_q_init(const VectorX q_init) {
    if (q_init.size() != _ndof_r) {
        std::cerr << "[Error] set_q_init: q_init.size() != _ndof_r." << std::endl;
        throw "error";
    }
    _q_init = q_init;
}

void Simulation::set_qdot_init(const VectorX qdot_init) {
    if (qdot_init.size() != _ndof_r) {
        std::cerr << "[Error] set_qdot_init: qdot_init.size() != _ndof_r." << std::endl;
        throw "error";
    }
    _qdot_init = qdot_init;
}

const VectorX Simulation::get_q_init() {
    return _q_init;
}

const VectorX Simulation::get_qdot_init() {
    return _qdot_init;
}

// states
void Simulation::set_state(const VectorX q, const VectorX qdot) {
    set_q(q);
    set_qdot(qdot);
}

void Simulation::set_q(const VectorX q) {
    _robot->set_q(q);
}

void Simulation::set_qdot(const VectorX qdot) {
    _robot->set_qdot(qdot);
}

const VectorX Simulation::get_q() {
    return _robot->get_q();
}

const VectorX Simulation::get_qdot() {
    return _robot->get_qdot();
}

void Simulation::reparam() {
    _robot->reparam();
}

const std::vector<VectorX> Simulation::get_q_his() {
    return _q_his;
}

const std::vector<VectorX> Simulation::get_qdot_his() {
    return _qdot_his;
}

void Simulation::set_state_his(std::vector<VectorX> q_his, std::vector<VectorX> qdot_his) {
    _backward_flag = false;
    _backward_design_params_flag = false;
    _backward_info.clear();

    _time_report.reset();
    _robot->reset_time_report();

    _q_his.clear();
    _qdot_his.clear();
    _h = _options->_h;

    for (int i = 0; i < q_his.size(); ++i) {
        _q_his.push_back(q_his[i]);
        _qdot_his.push_back(qdot_his[i]);
    }
}

const VectorX Simulation::get_variables() {
    return _robot->get_variables();
}

// control variables
void Simulation::set_u(const VectorX& u) {
    if (u.size() != _ndof_u) {
        std::cerr << "[Error] set_u: u.size() != _ndof_u." << std::endl;
        throw "error";
    }
    _robot->set_u(u);
}

void Simulation::get_ctrl_range(VectorX& ctrl_min, VectorX& ctrl_max) {
    _robot->get_ctrl_range(ctrl_min, ctrl_max);
}

void Simulation::print_ctrl_info() {
    _robot->print_ctrl_info();
}

// joint related
Joint* Simulation::get_joint(const string name) {
    auto it = _joint_map.find(name);
    if (it == _joint_map.end()) {
        std::string error_msg = "Get joint error: Joint " + name + " does not exist";
        throw_error(error_msg);
    }
    return _joint_map[name];
}

VectorX Simulation::get_joint_q(const string name, bool world_frame) {
    Joint* joint = get_joint(name);
    VectorX q = joint->_q;
    if (world_frame) {
        Matrix3 R2 = joint->_E_0j.topLeftCorner(3, 3);
        Vector3 p2 = joint->_E_0j.topRightCorner(3, 1);
        q.head(3) = (R2 * q.head(3) + p2).eval();
        q.tail(3) = (R2 * q.tail(3)).eval();
    }
    return q;
}

VectorX Simulation::get_joint_qdot(const string name, bool world_frame) {
    Joint* joint = get_joint(name);
    VectorX qdot = joint->_qdot;
    if (world_frame) {
        Matrix3 R2 = joint->_E_0j.topLeftCorner(3, 3);
        qdot.head(3) = (R2 * qdot.head(3)).eval();
        qdot.tail(3) = (R2 * qdot.tail(3)).eval();
    }
    return qdot;
}

void Simulation::set_joint_q(const string name, VectorX q, bool world_frame) {
    Joint* joint = get_joint(name);
    if (world_frame) {
        Matrix3 R2_inv = joint->_E_0j.topLeftCorner(3, 3).inverse();
        Vector3 p2_inv = -joint->_E_0j.topRightCorner(3, 1);
        joint->_q.head(3) = R2_inv * (q.head(3) + p2_inv);
        joint->_q.tail(3) = R2_inv * q.tail(3);
    } else {
        joint->_q.head(3) = q.head(3);
        joint->_q.tail(3) = q.tail(3);
    }
}

void Simulation::set_joint_qdot(const string name, VectorX qdot, bool world_frame) {
    Joint* joint = get_joint(name);
    if (world_frame) {
        Matrix3 R2_inv = joint->_E_0j.topLeftCorner(3, 3).inverse();
        joint->_qdot.head(3) = R2_inv * qdot.head(3);
        joint->_qdot.tail(3) = R2_inv * qdot.tail(3);
    } else {
        joint->_qdot.head(3) = qdot.head(3);
        joint->_qdot.tail(3) = qdot.tail(3);
    }
}

void Simulation::set_joint_state(const string name, VectorX q, VectorX qdot, bool world_frame) {
    Joint* joint = get_joint(name);
    set_joint_q(name, q, world_frame);
    set_joint_qdot(name, qdot, world_frame);
}

// body related
Body* Simulation::get_body(const string name) {
    auto it = _body_map.find(name);
    if (it == _body_map.end()) {
        std::string error_msg = "Get body error: Body " + name + " does not exist";
        throw_error(error_msg);
    }
    return _body_map[name];
}

dtype Simulation::get_body_mass(const string name) {
    return get_body(name)->_mass;
}

SE3 Simulation::get_body_E0i(const string name) {
    Body* body = get_body(name);
    return body->_E_0i;
}

SE3 Simulation::get_body_Ei0(const string name) {
    Body* body = get_body(name);
    return body->_E_i0;
}

Matrix3X Simulation::get_body_vertices(const string name, bool world_frame) {
    Body* body = get_body(name);
    Matrix3X vertices = body->get_vertices();
    if (world_frame) {
        Matrix3 R2 = body->_E_0i.topLeftCorner(3, 3);
        Vector3 p2 = body->_E_0i.topRightCorner(3, 1);
        for (int i = 0; i < vertices.cols(); ++i) {
            vertices.col(i) = (R2 * vertices.col(i) + p2).eval();
        }
    }
    return vertices;
}

Matrix3Xi Simulation::get_body_faces(const string name) {
    return get_body(name)->get_faces();
}

void Simulation::set_body_external_force(const string name, const VectorX& force) {
    get_body(name)->set_external_force(force);
}

dtype Simulation::get_body_distance(const string name_from, const string name_to) {
    Matrix3X vertices_from = get_body_vertices(name_from, true);
    Body* body_to = get_body(name_to);
    VectorX d;
    d.resize(vertices_from.cols(), 1);
    if (dynamic_cast<BodyPrimitiveShape*>(const_cast<Body*>(body_to)) != nullptr) {
        BodyPrimitiveShape* body_to_ = dynamic_cast<BodyPrimitiveShape*>(const_cast<Body*>(body_to));
        for (int i = 0; i < d.size(); ++i) {
            d(i) = body_to_->distance(vertices_from.col(i));
        }
    } else if (dynamic_cast<BodySDFObj*>(const_cast<Body*>(body_to)) != nullptr) {
        BodySDFObj* body_to_ = dynamic_cast<BodySDFObj*>(const_cast<Body*>(body_to));
        for (int i = 0; i < d.size(); ++i) {
            d(i) = body_to_->distance(vertices_from.col(i));
        }
    } else if (dynamic_cast<BodyBVHObj*>(const_cast<Body*>(body_to)) != nullptr) {
        BodyBVHObj* body_to_ = dynamic_cast<BodyBVHObj*>(const_cast<Body*>(body_to));
        body_to_->distance_parallel(vertices_from, d);
    }
    return d.minCoeff();
}

bool Simulation::body_in_contact(const string name_a, const string name_b, dtype eps) {
    dtype d_ab = get_body_distance(name_a, name_b);
    dtype d_ba = get_body_distance(name_b, name_a);
    dtype d_min = std::min(d_ab, d_ba);
    return d_min < eps;
}

std::vector<std::string> Simulation::get_contact_bodies(const string name) {
    Body* body = get_body(name);
    return body->get_contact_bodies();
}

void Simulation::clear_contact_bodies(const string name) {
    Body* body = get_body(name);
    body->clear_contact_bodies();
}

void Simulation::clear_all_contact_bodies() {
    for (auto iter = _body_map.begin(); iter != _body_map.end(); iter++) {
        Body* body = iter->second;
        body->clear_contact_bodies();
    }
}

void Simulation::clear_all_saved_sdfs() {
    for (auto iter = _body_map.begin(); iter != _body_map.end(); iter++) {
        Body* body = iter->second;
        if (dynamic_cast<BodySDFObj*>(const_cast<Body*>(body)) != nullptr) {
            BodySDFObj* body_ = dynamic_cast<BodySDFObj*>(const_cast<Body*>(body));
            body_->clear_saved_SDF();
        }
    }
}

// tactile data
std::vector<dtype> Simulation::get_tactile_depth(std::string name) {
    return _robot->get_tactile_depth(name);
}

std::vector<Vector2i> Simulation::get_tactile_image_pos(std::string name) {
    return _robot->get_tactile_image_pos(name);
}

// design parameters
void Simulation::set_design_params(const VectorX &design_params) {
    _robot->set_design_params(design_params);
}  

VectorX Simulation::get_design_params() {
    return _robot->get_design_params();
}

void Simulation::print_design_params_info() {
    if (_ndof_p > 0)
        _robot->print_design_params_info();
}

// set contact coefficient scale
void Simulation::set_contact_scale(dtype scale) {
    _robot->set_contact_scale(scale);
}

void Simulation::set_rendering_mesh_vertices(const std::vector<Matrix3X> &Vs) {
    _robot->set_rendering_mesh_vertices(Vs);
}

void Simulation::set_rendering_mesh(const std::vector<Matrix3X> &Vs, const std::vector<Matrix3Xi> &Fs) {
    _robot->set_rendering_mesh(Vs, Fs);
}

// virtual objects
void Simulation::update_virtual_object(std::string name, VectorX data) {
    _robot->update_virtual_object(name, data);
}

// update robot
void Simulation::update_robot(bool design_gradient) {
    _robot->update(design_gradient);
}

// compute the matrices M, f, so that M*qddot = f
void Simulation::computeMatrices(MatrixX& M, VectorX& f) {

    auto t_compute_matrices_start = clock();

    _robot->update();

    VectorX qdot = get_qdot();
    
    MatrixX J, Jdot; 
    _robot->computeJointJacobian(J, Jdot);

    VectorX Mm; 
    _robot->computeMaximalMassMatrix(Mm);

    auto t_compute_df_start = clock();

    VectorX fm, fr; 
    _robot->computeForce(fm, fr);

    _time_report._time_compute_df += clock() - t_compute_df_start;

    _time_report._time_compute_matrices += clock() - t_compute_matrices_start;

    auto t_compose_matrices_start = clock();

    MatrixX JT = J.transpose();
    
    MatrixX MmJ(_ndof_m, _ndof_r);
    for (int i = 0;i < MmJ.rows();i++) 
        MmJ.row(i).noalias() = J.row(i) * Mm(i);

    MatrixX MmJdot(_ndof_m, _ndof_r);
    for (int i = 0;i < MmJdot.rows();i++)
        MmJdot.row(i).noalias() = Jdot.row(i) * Mm(i);

    M = JT * MmJ;
    f = JT * (fm - MmJdot * qdot) + fr;

    _time_report._time_compose_matrices += clock() - t_compose_matrices_start;
}

// compute the matrices M, f, so that M*qddot = f;
// compute derivatives dM_dq, K = df_dq, D = df_dqdot
void Simulation::computeMatrices(MatrixX& M, VectorX& f, JacobianMatrixVector& dM_dq, MatrixX& K, MatrixX& D) {

    auto t_compute_matrices_start = clock();

    _robot->update();

    VectorX qdot = get_qdot();

    MatrixX J, Jdot; 
    JacobianMatrixVector dJ_dq, dJdot_dq; 
    _robot->computeJointJacobianWithDerivative(J, Jdot, dJ_dq, dJdot_dq);

    VectorX Mm; 
    _robot->computeMaximalMassMatrix(Mm);

    auto t_compute_df_start = clock();

    VectorX fm, fr;
    MatrixX Km, Dm, Kr, Dr;
    _robot->computeForceWithDerivative(fm, fr, Km, Dm, Kr, Dr);

    _time_report._time_compute_df += clock() - t_compute_df_start;
    
    _time_report._time_compute_matrices += clock() - t_compute_matrices_start;

    auto t_compose_matrices_start = clock();

    MatrixX JT = J.transpose();

    MatrixX MmJ(_ndof_m, _ndof_r);
    for (int i = 0;i < MmJ.rows();i++)
        MmJ.row(i).noalias() = J.row(i) * Mm(i);

    MatrixX MmJdot(_ndof_m, _ndof_r);
    for (int i = 0;i < MmJdot.rows();i++)
        MmJdot.row(i).noalias() = Jdot.row(i) * Mm(i);

    M = JT * MmJ;

    MatrixX JTMm = MmJ.transpose();

    // fqvv = -J'*Mm*Jdot*qdot
    MatrixX JTMmJdot = JTMm * Jdot;
    VectorX fqvv = -JTMmJdot * qdot;

    f = JT * fm + fqvv + fr;
    
    /********************************* Derivatives ************************************/
    // dM_dq = (dJ_dq)'*Mm*J + J'*Mm*dJ_dq
    dM_dq = JacobianMatrixVector(_ndof_r, _ndof_r, _ndof_r);
    for (int k = 0;k < _ndof_r;k++) {
        MatrixX tmp = dJ_dq(k).transpose() * MmJ;
        dM_dq(k) = tmp + tmp.transpose();
    }
    
    // dJdq(k) * qdot
    MatrixX dJdq_qdot(_ndof_m, _ndof_r);
    for (int k = 0;k < _ndof_r;k++) {
        dJdq_qdot.col(k) = dJ_dq(k) * qdot;
    }

    // Kqvv = dfqvv_dq = -(dJ_dq)'*Mm*Jdot*qdot-J'*Mm*dJdot_dq*qdot
    // Dqvv = dfqvv_dqdot = -J'*Mm*dJdot_dqdot*qdot-J'*Mm*Jdot = -J'*Mm*dJ_dq*qdot-J'*Mm*Jdot (using the fact dJdot_dqdot = dJ_dq)
    MatrixX Kqvv(_ndof_r, _ndof_r);
    MatrixX Dqvv = -JTMmJdot - JTMm * dJdq_qdot;
    MatrixX MmJdotqdot = MmJdot * qdot;
    for (int k = 0;k < _ndof_r;k++) {
        Kqvv.col(k) = -dJ_dq(k).transpose() * MmJdotqdot - JTMm * (dJdot_dq(k) * qdot);
    }

    // Km = d(J'*fm)_dq = dJ'_dq * fm + J'*(Km*J + Dm * dqmdot_dqr), notes eq (3.49)
    // Dm = d(J'*fm)_dqdot = J'*Dm*J
    MatrixX JTDm = JT * Dm;
    K = Kqvv + Kr + JT * Km * J + JTDm * dJdq_qdot;
    for (int k = 0;k < _ndof_r;k++) {
        K.col(k) += dJ_dq(k).transpose() * fm;
    }

    D = Dqvv + Dr + JTDm * J;

    _time_report._time_compose_matrices += clock() - t_compose_matrices_start;
}

// compute the matrices M, f, so that M*qddot = f;
// compute derivatives dM_dq, K = df_dq, D = df_dqdot
// compute derivatives w.r.t control
void Simulation::computeMatrices(MatrixX& M, VectorX& f, JacobianMatrixVector& dM_dq, MatrixX& K, MatrixX& D, MatrixX& df_du) {

    auto t_compute_matrices_start = clock();

    _robot->update();

    VectorX qdot = get_qdot();

    MatrixX J, Jdot; 
    JacobianMatrixVector dJ_dq, dJdot_dq; 
    _robot->computeJointJacobianWithDerivative(J, Jdot, dJ_dq, dJdot_dq);

    VectorX Mm; 
    _robot->computeMaximalMassMatrix(Mm);

    auto t_compute_df_start = clock();
    
    VectorX fm, fr;
    MatrixX Km, Dm, Kr, Dr;
    _robot->computeForceWithDerivative(fm, fr, Km, Dm, Kr, Dr);

    _time_report._time_compute_df += clock() - t_compute_df_start;

    _time_report._time_compute_matrices += clock() - t_compute_matrices_start;

    auto t_compose_matrices_start = clock();

    MatrixX JT = J.transpose();

    MatrixX MmJ(_ndof_m, _ndof_r);
    for (int i = 0;i < MmJ.rows();i++)
        MmJ.row(i).noalias() = J.row(i) * Mm(i);

    MatrixX MmJdot(_ndof_m, _ndof_r);
    for (int i = 0;i < MmJdot.rows();i++)
        MmJdot.row(i).noalias() = Jdot.row(i) * Mm(i);

    M = JT * MmJ;

    MatrixX JTMm = MmJ.transpose(); // TODO: need to check whether Mm is symmetric

    // fqvv = -J'*Mm*Jdot*qdot
    MatrixX JTMmJdot = JTMm * Jdot;
    VectorX fqvv = -JTMmJdot * qdot;

    f = JT * fm + fqvv + fr;
    
    /********************************* Derivatives ************************************/
    // dM_dq = (dJ_dq)'*Mm*J + J'*Mm*dJ_dq
    dM_dq = JacobianMatrixVector(_ndof_r, _ndof_r, _ndof_r);
    for (int k = 0;k < _ndof_r;k++) {
        MatrixX tmp = dJ_dq(k).transpose() * MmJ;
        dM_dq(k) = tmp + tmp.transpose();
    }

    // dJdq(k) * qdot
    MatrixX dJdq_qdot(_ndof_m, _ndof_r);
    for (int k = 0;k < _ndof_r;k++) {
        dJdq_qdot.col(k) = dJ_dq(k) * qdot;
    }

    // Kqvv = dfqvv_dq = -(dJ_dq)'*Mm*Jdot*qdot-J'*Mm*dJdot_dq*qdot
    // Dqvv = dfqvv_dqdot = -J'*Mm*dJdot_dqdot*qdot-J'*Mm*Jdot = -J'*Mm*dJ_dq*qdot-J'*Mm*Jdot (using the fact dJdot_dqdot = dJ_dq)
    MatrixX Kqvv(_ndof_r, _ndof_r);
    MatrixX Dqvv = -JTMmJdot - JTMm * dJdq_qdot;
    MatrixX MmJdotqdot = MmJdot * qdot;
    for (int k = 0;k < _ndof_r;k++) {
        Kqvv.col(k) = -dJ_dq(k).transpose() * MmJdotqdot - JTMm * (dJdot_dq(k) * qdot);
    }

    // Km = d(J'*fm)_dq = dJ'_dq * fm + J'*(Km*J + Dm * dqmdot_dqr), notes eq (3.49)
    // Dm = d(J'*fm)_dqdot = J'*Dm*J
    MatrixX JTDm = JT * Dm;
    K = Kqvv + Kr + JT * Km * J + JTDm * dJdq_qdot;
    for (int k = 0;k < _ndof_r;k++) {
        K.col(k) += dJ_dq(k).transpose() * fm;
    }

    D = Dqvv + Dr + JTDm * J;

    // df_du = JT * dfm_du + dfr_du
    MatrixX dfm_du, dfr_du;
    _robot->compute_dfdu(dfm_du, dfr_du);
    df_du = JT * dfm_du + dfr_du;

    _time_report._time_compose_matrices += clock() - t_compose_matrices_start;
}

// compute the matrices M, f, so that M*qddot = f;
// compute derivatives dM_dq, K = df_dq, D = df_dqdot
// compute derivatives w.r.t control
// compute derivatives w.r.t design parameters
void Simulation::computeMatrices(
        MatrixX& M, VectorX& f,
        JacobianMatrixVector& dM_dq, MatrixX& K, MatrixX& D,
        MatrixX& df_du,
        JacobianMatrixVector& dM_dp, MatrixX& df_dp) {

    auto t_compute_matrices_start = clock();

    _robot->update();

    VectorX qdot = get_qdot();

    auto t_compute_dJ_start = clock();

    MatrixX J, Jdot; 
    JacobianMatrixVector dJ_dq, dJdot_dq; 
    SparseJacobianMatrixVector dJ_dp1, dJ_dp2, dJdot_dp1, dJdot_dp2;
    // JacobianMatrixVector dJ_dp1, dJ_dp2, dJdot_dp1, dJdot_dp2;
    _robot->computeJointJacobianWithDerivative(J, Jdot, dJ_dq, dJdot_dq, dJ_dp1, dJ_dp2, dJdot_dp1, dJdot_dp2);
    
    _time_report._time_compute_dJ += clock() - t_compute_dJ_start;

    VectorX Mm; 
    _robot->computeMaximalMassMatrix(Mm);
    
    auto t_compute_df_start = clock();

    VectorX fm, fr;
    MatrixX Km, Dm, Kr, Dr;
    // SparseMatrixX dfm_dp, dfr_dp;
    MatrixX dfm_dp, dfr_dp;
    _robot->computeForceWithDerivative(fm, fr, Km, Dm, Kr, Dr, dfm_dp, dfr_dp);

    _time_report._time_compute_df += clock() - t_compute_df_start;

    _time_report._time_compute_matrices += clock() - t_compute_matrices_start;

    auto t_compose_matrices_start = clock();

    MatrixX JT = J.transpose();

    MatrixX MmJ(_ndof_m, _ndof_r);
    for (int i = 0;i < MmJ.rows();i++)
        MmJ.row(i).noalias() = J.row(i) * Mm(i);

    MatrixX MmJdot(_ndof_m, _ndof_r);
    for (int i = 0;i < MmJdot.rows();i++)
        MmJdot.row(i).noalias() = Jdot.row(i) * Mm(i);

    M = JT * MmJ;

    MatrixX JTMm = MmJ.transpose();

    // fqvv = -J'*Mm*Jdot*qdot
    MatrixX JTMmJdot = JTMm * Jdot;
    VectorX Jdotqdot = Jdot * qdot;
    VectorX fqvv = -JTMmJdot * qdot;

    f = JT * fm + fqvv + fr;
    
    /********************************* Derivatives ************************************/
    // dM_dq = (dJ_dq)'*Mm*J + J'*Mm*dJ_dq
    dM_dq = JacobianMatrixVector(_ndof_r, _ndof_r, _ndof_r);
    for (int k = 0;k < _ndof_r;k++) {
        MatrixX tmp = dJ_dq(k).transpose() * MmJ;
        dM_dq(k) = tmp + tmp.transpose();
    }

    // dJdq(k) * qdot
    MatrixX dJdq_qdot(_ndof_m, _ndof_r);
    for (int k = 0;k < _ndof_r;k++) {
        dJdq_qdot.col(k) = dJ_dq(k) * qdot;
    }

    // Kqvv = dfqvv_dq = -(dJ_dq)'*Mm*Jdot*qdot-J'*Mm*dJdot_dq*qdot
    // Dqvv = dfqvv_dqdot = -J'*Mm*dJdot_dqdot*qdot-J'*Mm*Jdot = -J'*Mm*dJ_dq*qdot-J'*Mm*Jdot (using the fact dJdot_dqdot = dJ_dq)
    MatrixX Kqvv(_ndof_r, _ndof_r);
    MatrixX Dqvv = -JTMmJdot - JTMm * dJdq_qdot;
    MatrixX MmJdotqdot = MmJdot * qdot;
    for (int k = 0;k < _ndof_r;k++) {
        Kqvv.col(k) = -dJ_dq(k).transpose() * MmJdotqdot - JTMm * (dJdot_dq(k) * qdot);
    }

    // Km = d(J'*fm)_dq = dJ'_dq * fm + J'*(Km*J + Dm * dqmdot_dqr), notes eq (3.49)
    // Dm = d(J'*fm)_dqdot = J'*Dm*J
    MatrixX JTDm = JT * Dm;
    K = Kqvv + Kr + JT * Km * J + JTDm * dJdq_qdot;
    for (int k = 0;k < _ndof_r;k++) {
        K.col(k) += dJ_dq(k).transpose() * fm;
    }

    D = Dqvv + Dr + JTDm * J;

    // df_du = JT * dfm_du + dfr_du
    MatrixX dfm_du, dfr_du;
    _robot->compute_dfdu(dfm_du, dfr_du);
    df_du = JT * dfm_du + dfr_du;

    // design derivatives

    auto t_dM_dp_start = clock();

    // dM_dp
    auto t_dM_dp1_start = clock();
    // design params 1
    for (int i = 0;i < _ndof_p1;i++) {
        MatrixX tmp = dJ_dp1(i).transpose() * MmJ;
        dM_dp(i) = tmp + tmp.transpose();
    }
    _time_report._time_dM_dp1 += clock() - t_dM_dp1_start;
    auto t_dM_dp2_start = clock();
    // design params 2
    for (int i = 0;i < _ndof_p2;i++) {
        MatrixX tmp = dJ_dp2(i).transpose() * MmJ;
        dM_dp(_ndof_p1 + i) = tmp + tmp.transpose();
    }
    _time_report._time_dM_dp2 += clock() - t_dM_dp2_start;
    auto t_dM_dp4_start = clock();
    // design params 4
    int p4_offset = _ndof_p1 + _ndof_p2 + _ndof_p3;
    for (auto body : _robot->_bodies) {
        if (body->_design_params_4._active) {
            int m_id = body->_index[0];
            int p_id = body->_design_params_4._param_index[0];
            // mass
            dM_dp(p4_offset + p_id) = JT.middleCols(m_id + 3, 3) * J.middleRows(m_id + 3, 3);
            // inertia
            for (int k = 0;k < 3;k++) {
                dM_dp(p4_offset + p_id + 1 + k) = JT.col(m_id + k) * J.row(m_id + k);
            }
        }
    }
    _time_report._time_dM_dp4 += clock() - t_dM_dp4_start;
    _time_report._time_dM_dp += clock() - t_dM_dp_start;

    auto t_df_dp_start = clock();

    // df_dp
    df_dp = dfr_dp;

    // fqvv component
    // design params 1
    for (int i = 0;i < _ndof_p1;i++) {
        df_dp.col(i) -= dJ_dp1(i).transpose() * MmJdotqdot + JTMm * dJdot_dp1(i) * qdot;
    }
    // design params 2
    for (int i = 0;i < _ndof_p2;i++) {
        df_dp.col(i + _ndof_p1) -= dJ_dp2(i).transpose() * MmJdotqdot + JTMm * dJdot_dp2(i) * qdot;
    }
    // design params 4
    for (auto body : _robot->_bodies) {
        if (body->_design_params_4._active) {
            int m_id = body->_index[0];
            int p_id = body->_design_params_4._param_index[0];
            // mass
            df_dp.col(p4_offset + p_id) -= JT.middleCols(m_id + 3, 3) * Jdotqdot.segment(m_id + 3, 3);
            // inertia
            for (int k = 0;k < 3;k++) {
                df_dp.col(p4_offset + p_id + 1 + k) -= JT.col(m_id + k) * Jdotqdot(m_id + k);
            }
        }
    }

    // fm component: dJdp * fm + JT * dfm_dp
    // dJdp * fm
    // design params 1
    for (int i = 0;i < _ndof_p1;i++) {
        df_dp.col(i) += dJ_dp1(i).transpose() * fm;
    }
    // design params 2
    for (int i = 0;i < _ndof_p2;i++) {
        df_dp.col(i + _ndof_p1) += dJ_dp2(i).transpose() * fm;
    }
    // JT * dfm_dp
    for (int i = 0;i < _ndof_p;i++) {
        df_dp.col(i) += JT * dfm_dp.col(i);
    }

    _time_report._time_df_dp += clock() - t_df_dp_start;

    _time_report._time_compose_matrices += clock() - t_compose_matrices_start;
}

void Simulation::computeVariablesWithDerivative(VectorX& variables, MatrixX& dvar_dq) {
    _robot->computeVariablesWithDerivative(variables, dvar_dq);
}

void Simulation::computeVariablesWithDerivative(VectorX& variables, MatrixX& dvar_dq, MatrixX& dvar_dp) {
    _robot->computeVariablesWithDerivative(variables, dvar_dq, dvar_dp);
}

void Simulation::test_derivatives_runtime() {
    _robot->test_derivatives_runtime();

    // test dM_dq, K, D
    _robot->update();

    auto q = get_q();
    auto qdot = get_qdot();

    MatrixX M;
    VectorX f;
    JacobianMatrixVector dM_dq;
    MatrixX K, D;
    computeMatrices(M, f, dM_dq, K, D);

    // printf("**************************** Simulation ****************************\n");
    dtype h = 1e-7;
    for (int ii = 0;ii < 1;ii++) {
        // printf("---------------------- eps = %.9lf ----------------------------\n", h);
        JacobianMatrixVector dM_dq_fd(_ndof_r, _ndof_r, _ndof_r);
        MatrixX K_fd = MatrixX::Zero(_ndof_r, _ndof_r);
        for (int k = 0;k < _ndof_r;k++) {
            auto q_pos = q;
            q_pos(k) += h;
            set_q(q_pos);

            _robot->update();

            MatrixX M_pos;
            VectorX f_pos;
            computeMatrices(M_pos, f_pos);

            dM_dq_fd(k) = (M_pos - M) / h;
            K_fd.col(k) = (f_pos - f) / h;
        }

        print_error("Simulation: dM_dq", dM_dq, dM_dq_fd);
        print_error("Simulation: K", K, K_fd);

        set_q(q);
        _robot->update();

        MatrixX D_fd = MatrixX::Zero(_ndof_r, _ndof_r);
        for (int k = 0;k < _ndof_r;k++) {
            auto qdot_pos = qdot;
            qdot_pos(k) += h;
            set_qdot(qdot_pos);

            _robot->update();

            MatrixX M_pos;
            VectorX f_pos;
            computeMatrices(M_pos, f_pos);

            D_fd.col(k) = (f_pos - f) / h;
        }

        print_error("Simulation: D", D, D_fd);

        set_qdot(qdot);
        _robot->update();

        h /= 10.;
    }

    // test_design_derivatives_runtime();
}

void Simulation::test_design_derivatives_runtime() {
    dtype eps = 1e-8;

    VectorX design_params = get_design_params();
    _robot->update(true);

    MatrixX M;
    VectorX f;
    JacobianMatrixVector dM_dq;
    MatrixX K, D;
    MatrixX df_du;
    computeMatrices(M, f, dM_dq, K, D, df_du, _dM_dp, _df_dp);

    if (_ndof_p > 0) {
        JacobianMatrixVector dM_dp_fd(_ndof_r, _ndof_r, _ndof_p);
        MatrixX df_dp_fd = MatrixX::Zero(_ndof_r, _ndof_p);
        for (int i = 0;i < _ndof_p;i++) {
            VectorX design_params_pos = design_params;
            design_params_pos(i) += eps;
            set_design_params(design_params_pos);
            _robot->update(false);
            MatrixX M_pos;
            VectorX f_pos;
            computeMatrices(M_pos, f_pos);
            dM_dp_fd(i) = (M_pos - M) / eps;
            df_dp_fd.col(i) = (f_pos - f) / eps;
        }
        print_error("Simulation: dM_dp", _dM_dp, dM_dp_fd);
        print_error("Simulation: df_dp", _df_dp, df_dp_fd);
    }

    set_design_params(design_params);
    _robot->update(true);
}

void Simulation::reset(bool backward_flag, bool backward_design_params_flag) {
    set_state(_q_init, _qdot_init);
    clear_all_contact_bodies();
    update_robot();
    _robot->reset_external_force();
    _robot->update_actuator_dofs(_q_init, _qdot_init);
    
    _q_his.clear();
    _qdot_his.clear();
    _q_his.push_back(_q_init);
    _qdot_his.push_back(_qdot_init);

    _backward_flag = backward_flag;
    _backward_design_params_flag = backward_design_params_flag;
    _backward_info.clear();

    _time_report.reset();
    _robot->reset_time_report();
}

void Simulation::forward(int num_steps, bool verbose, bool test_derivatives) {
    
    _verbose = verbose;

    if (_q_his.size() == 0) {
        std::cerr << "[Error] Please call simulation.reset() before simulation.forward()." << std::endl;
        throw "error";
    }

    for (int i = 0;i < num_steps;i++) {
        VectorX q = get_q();
        VectorX qdot = get_qdot();

        // if (verbose) {
        //     std::cerr << "q = " << q.transpose() << ", qdot = " << qdot.transpose() << std::endl;
        // }

        assert((q - _q_his[_q_his.size() - 1]).norm() < 1e-7 && (qdot - _qdot_his[_qdot_his.size() - 1]).norm() < 1e-7);

        VectorX q_next, qdot_next;
        if (_options->_integrator == "BDF1") {
            integration_BDF1(q, qdot, _options->_h, q_next, qdot_next);
        } else if (_options->_integrator == "BDF2") {
            if (_q_his.size() == 1) {
                integration_BDF1(q, qdot, _options->_h, q_next, qdot_next);
            } else{
                VectorX q_prev = _q_his[_q_his.size() - 2];
                VectorX qdot_prev = _qdot_his[_qdot_his.size() - 2];
                integration_BDF2(q_prev, qdot_prev, q, qdot, _options->_h, q_next, qdot_next);
            }
        } else if (_options->_integrator == "Euler") {
            integration_Euler(q, qdot, _options->_h, q_next, qdot_next);
        } else if (_options->_integrator == "Midpoint") {
            integration_Midpoint(q, qdot, _options->_h, q_next, qdot_next);
        } else if (_options->_integrator == "RK4") {
            integration_RK4(q, qdot, _options->_h, q_next, qdot_next);
        } else {
            std::cerr << "[Error] Integrator " << _options->_integrator << " has not been implemented." << std::endl;
            throw "error";
        }

        set_state(q_next, qdot_next);
        reparam();
        q_next = get_q();
        qdot_next = get_qdot();
        
        update_robot(_backward_flag);

        // update actuator's dofs
        _robot->update_actuator_dofs(q_next, qdot_next);

        _q_his.push_back(q_next);
        _qdot_his.push_back(qdot_next);

        if (_backward_flag) {
            if (_backward_design_params_flag) {
                VectorX variables;
                MatrixX dvar_dq, dvar_dp;
                computeVariablesWithDerivative(variables, dvar_dq, dvar_dp);
                _backward_info._dvar_dq.push_back(dvar_dq);
                _backward_info._dvar_dp.push_back(dvar_dp);
            } else {
                VectorX variables;
                MatrixX dvar_dq;
                computeVariablesWithDerivative(variables, dvar_dq);
                _backward_info._dvar_dq.push_back(dvar_dq);
            }
        }

        if (test_derivatives && !_newton_converge) {
            test_derivatives_runtime();
        }

        // VectorX fm = VectorX::Zero(_ndof_m);
        // VectorX fr = VectorX::Zero(_ndof_r);
        // for (auto force : _robot->_forces) {
        //     force->computeForce(fm, fr, true);
        // }
    }
}

bool Simulation::is_converged() const {
    return _newton_converge;
}

void Simulation::newton(VectorX& x,
                        Func func, 
                        Func_With_Derivatives func_with_derivatives) {
    
    auto t_newton_start = clock();

    if (_ndof_r == 0)
        return;

    dtype tol = 1e-8;
    int MaxIter_Newton = max(20 * _ndof_r, 100);
    // int MaxIter_Newton = max(100 * _ndof_r, 100);
    // int MaxIter_LS = 300;
    int MaxIter_LS = 20;
    // int MaxIter_LS_Fail_Strike = 1;
    int MaxIter_LS_Fail_Strike = 3;

    bool success_newton = false;
    dtype g_last;
    int fail_strike = 0;
    std::vector<dtype> g_his;
    // std::cerr << "newton's solver" << std::endl;
    for (int iter_newton = 0;iter_newton < MaxIter_Newton;iter_newton ++) {
        VectorX g;
        MatrixX H;
        (this->*func_with_derivatives)(x, g, H, false);

        // new_x = -inv(H) * g + x
        // VectorX dx = -H.inverse() * g;
        VectorX dx = H.partialPivLu().solve(-g);

        // line search new_x = x + alpha * dx
        VectorX g_new;
        dtype gnorm = g.norm();
        dtype alpha = 1.;
        bool success_ls = false;
        // for (int trial = 0;trial < MaxIter_LS;trial ++, alpha *= 0.95) {
        for (int trial = 0;trial < MaxIter_LS;trial ++, alpha *= 0.5) {
            auto t_eval_start = clock();

            (this->*func)(x + alpha * dx, g_new);

            if (g_new.norm() < gnorm) {
                success_ls = true;
                break;
            }
        }

        // std::cerr << "g_new = " << g_new.norm() << ", g = " << gnorm << std::endl;

        if (success_ls) {
            fail_strike = 0;
        } else {
            fail_strike += 1;
            if (fail_strike >= MaxIter_LS_Fail_Strike)
                break;
        }

        x = x + alpha * dx;

        if (g_new.norm() < tol) {
            success_newton = true;
            break;
        }

        g_last = g_new.norm();
        g_his.push_back(g_last);
    }

    if (!success_newton) {
        _newton_converge = false;
        if (_verbose) {
            std::cerr << "g = " << g_last << std::endl;
            printf("Newton method did not converge\n");
        }
    } else {
        _newton_converge = true;
    }

    _time_report._time_solver += clock() - t_newton_start;
}

void Simulation::evaluate_g_BDF1(const VectorX& q1, VectorX& g) {
    VectorX qdot1 = (q1 - _q0) / _h;

    set_state(q1, qdot1);
    update_robot();

    MatrixX M;
    VectorX f;
    computeMatrices(M, f);

    g = M * (q1 - _q0 - _h * _qdot0) - _h * _h * f;
}

void Simulation::evaluate_g_with_derivatives_BDF1(const VectorX& q1, VectorX& g, MatrixX& H, bool save_backward_info) {
    VectorX qdot1 = (q1 - _q0) / _h;

    set_state(q1, qdot1);
    update_robot(save_backward_info && _backward_design_params_flag);

    if (!save_backward_info) {
        MatrixX M;
        VectorX f;
        JacobianMatrixVector dM_dq;
        MatrixX K, D;
        computeMatrices(M, f, dM_dq, K, D);

        VectorX dq_tmp = q1 - _q0 - _h * _qdot0;
        g = M * dq_tmp - _h * _h * f;
        H = M - _h * _h * K - _h * D;
        for (int k = 0;k < _ndof_r;k++) {
            H.col(k) += dM_dq(k) * dq_tmp;
        }
    } else {
        if (_backward_design_params_flag) {
            MatrixX M;
            VectorX f;
            JacobianMatrixVector dM_dq;
            MatrixX K, D;
            MatrixX df_du;
            computeMatrices(M, f, dM_dq, K, D, df_du, _dM_dp, _df_dp);

            VectorX dq_tmp = q1 - _q0 - _h * _qdot0;
            g = M * dq_tmp - _h * _h * f;
            H = M - _h * _h * K - _h * D;
            for (int k = 0;k < _ndof_r;k++) {
                H.col(k) += dM_dq(k) * dq_tmp;
            }
            MatrixX dg_dp = -_h * _h * _df_dp;
            // for (int k = 0;k < _ndof_p;k++) {
            //     dg_dp.col(k) += dM_dp(k) * dq_tmp;
            // }
            for (int k = 0;k < _ndof_p1 + _ndof_p2;k++) {
                dg_dp.col(k) += _dM_dp(k) * dq_tmp;
            }
            for (int k = _ndof_p1 + _ndof_p2 + _ndof_p3;k < _ndof_p1 + _ndof_p2 + _ndof_p3 + _ndof_p4;k++) {
                dg_dp.col(k) += _dM_dp(k) * dq_tmp;
            }

            // save backward data
            _backward_info._M.push_back(M);
            _backward_info._D.push_back(D);
            _backward_info._dg_dp.push_back(dg_dp);
            _backward_info._dg_du.push_back(-_h * _h * df_du);
        } else {
            MatrixX M;
            VectorX f;
            JacobianMatrixVector dM_dq;
            MatrixX K, D;
            MatrixX df_du;
            computeMatrices(M, f, dM_dq, K, D, df_du);

            VectorX dq_tmp = q1 - _q0 - _h * _qdot0;
            g = M * dq_tmp - _h * _h * f;
            H = M - _h * _h * K - _h * D;
            for (int k = 0;k < _ndof_r;k++) {
                H.col(k) += dM_dq(k) * dq_tmp;
            }
            
            // save backward data
            _backward_info._M.push_back(M);
            _backward_info._D.push_back(D);
            _backward_info._dg_du.push_back(-_h * _h * df_du);
        }
    }
}

void Simulation::integration_Euler(
    const VectorX q0, const VectorX qdot0, const dtype h, 
        VectorX& q1, VectorX& qdot1) {
    _q0 = q0;
    _qdot0 = qdot0;
    _h = h;
    MatrixX M;
    VectorX f;

    set_state(q0, qdot0);
    update_robot();
    computeMatrices(M, f);

    q1 = q0 + h * qdot0;
    qdot1 = qdot0 + h * M.inverse() * f;
}

void Simulation::integration_Midpoint(
    const VectorX q0, const VectorX qdot0, const dtype h, 
        VectorX& q1, VectorX& qdot1) {
    _q0 = q0;
    _qdot0 = qdot0;
    _h = h;
    MatrixX M;
    VectorX f;
    VectorX q0_5, qdot0_5;

    set_state(q0, qdot0);
    update_robot();
    computeMatrices(M, f);

    q0_5 = q0 + h / 2. * qdot0;
    qdot0_5 = qdot0 + h / 2. * M.inverse() * f;

    set_state(q0_5, qdot0_5);
    update_robot();
    computeMatrices(M, f);

    q1 = q0 + h * qdot0_5;
    qdot1 = qdot0 + h * M.inverse() * f;
}

void Simulation::integration_RK4(
    const VectorX q0, const VectorX qdot0, const dtype h, 
        VectorX& q1, VectorX& qdot1) {
    _q0 = q0;
    _qdot0 = qdot0;
    _h = h;
    MatrixX M;
    VectorX f;
    VectorX kq1, kq2, kq3, kq4;
    VectorX kqdot1, kqdot2, kqdot3, kqdot4;

    kq1 = qdot0;
    set_state(q0, qdot0);
    update_robot();
    computeMatrices(M, f);
    kqdot1 = M.inverse() * f;

    kq2 = qdot0 + h / 2. * kqdot1;
    set_state(q0 + h / 2. * kq1, kq2);
    update_robot();
    computeMatrices(M, f);
    kqdot2 = M.inverse() * f;

    kq3 = qdot0 + h / 2. * kqdot2;
    set_state(q0 + h / 2. * kq2, kq3);
    update_robot();
    computeMatrices(M, f);
    kqdot3 = M.inverse() * f;

    kq4 = qdot0 + h * kqdot3;
    set_state(q0 + h * kq3, kq4);
    update_robot();
    computeMatrices(M, f);
    kqdot4 = M.inverse() * f;

    q1 = q0 + h / 6. * (kq1 + 2 * kq2 + 2 * kq3 + kq4);
    qdot1 = qdot0 + h / 6. * (kqdot1 + 2 * kqdot2 + 2 * kqdot3 + kqdot4);
}

void Simulation::integration_BDF1(
        const VectorX q0, const VectorX qdot0, const dtype h, 
        VectorX& q1, VectorX& qdot1) {
    _q0 = q0;
    _qdot0 = qdot0;
    _h = h;

    q1 = q0 + h * qdot0; // initial guess

    newton(q1, &Simulation::evaluate_g_BDF1, &Simulation::evaluate_g_with_derivatives_BDF1);

    qdot1 = (q1 - q0) / h;

    // save backward data: H_inv, M, D, dg_du
    if (_backward_flag) {
        auto t_save_backward_start = clock();
        VectorX g;
        MatrixX H;
        evaluate_g_with_derivatives_BDF1(q1, g, H, true);
        _backward_info._H_lu.push_back(H.partialPivLu());
        _time_report._time_save_backward += clock() - t_save_backward_start;
    }
}

void Simulation::evaluate_g_SDIRK2b(const VectorX& q1, VectorX& g) {
    dtype alpha = (2. - sqrt(2.)) / 2.;
    VectorX qdot1 = (q1 + (1. / alpha - 2.) * _q0 - (1. - alpha) / alpha * _q_alpha) / (alpha * _h);

    set_state(q1, qdot1);
    update_robot();

    MatrixX M;
    VectorX f;
    computeMatrices(M, f);

    g = M * (q1 - _q0 - (2. * alpha - 1.) * _h * _qdot0 - 2. * (1. - alpha) * _h * _qdot_alpha) - alpha * alpha * _h * _h * f;
}

void Simulation::evaluate_g_with_derivatives_SDIRK2b(const VectorX& q1, VectorX& g, MatrixX& H, bool save_backward_info) {
    dtype alpha = (2. - sqrt(2.)) / 2.;
    VectorX qdot1 = (q1 + (1. / alpha - 2.) * _q0 - (1. - alpha) / alpha * _q_alpha) / (alpha * _h);
    
    set_state(q1, qdot1);
    update_robot(save_backward_info && _backward_design_params_flag);

    if (!save_backward_info) {
        MatrixX M;
        VectorX f;
        JacobianMatrixVector dM_dq;
        MatrixX K, D;
        computeMatrices(M, f, dM_dq, K, D);

        VectorX dq_tmp = q1 - _q0 - (2. * alpha - 1.) * _h * _qdot0 - 2. * (1. - alpha) * _h * _qdot_alpha;
        g = M * dq_tmp - alpha * alpha * _h * _h * f;
        H = M - alpha * alpha * _h * _h * K - alpha * _h * D;
        for (int k = 0;k < _ndof_r;k++) {
            H.col(k) += dM_dq(k) * dq_tmp;
        }
    } else {
        if (_backward_design_params_flag) {
            MatrixX M;
            VectorX f;
            JacobianMatrixVector dM_dq;
            MatrixX K, D;
            MatrixX df_du;
            computeMatrices(M, f, dM_dq, K, D, df_du, _dM_dp, _df_dp);

            VectorX dq_tmp = q1 - _q0 - (2. * alpha - 1.) * _h * _qdot0 - 2. * (1. - alpha) * _h * _qdot_alpha;
            g = M * dq_tmp - alpha * alpha * _h * _h * f;
            H = M - alpha * alpha * _h * _h * K - alpha * _h * D;
            for (int k = 0;k < _ndof_r;k++) {
                H.col(k) += dM_dq(k) * dq_tmp;
            }
            MatrixX dg_dp = -alpha * alpha * _h * _h * _df_dp;
            // for (int k = 0;k < _ndof_p;k++) {
            //     dg_dp.col(k) += dM_dp(k) * dq_tmp;
            // }
            for (int k = 0;k < _ndof_p1 + _ndof_p2;k++) {
                dg_dp.col(k) += _dM_dp(k) * dq_tmp;
            }
            for (int k = _ndof_p1 + _ndof_p2 + _ndof_p3;k < _ndof_p1 + _ndof_p2 + _ndof_p3 + _ndof_p4;k++) {
                dg_dp.col(k) += _dM_dp(k) * dq_tmp;
            }

            // save backward_data
            _backward_info._M.push_back(M);
            _backward_info._D.push_back(D);
            _backward_info._dg_dp.push_back(dg_dp);
            _backward_info._dg_du.push_back(-alpha * alpha * _h * _h * df_du);
        } else {
            MatrixX M;
            VectorX f;
            JacobianMatrixVector dM_dq;
            MatrixX K, D;
            MatrixX df_du;
            computeMatrices(M, f, dM_dq, K, D, df_du);

            VectorX dq_tmp = q1 - _q0 - (2. * alpha - 1.) * _h * _qdot0 - 2. * (1. - alpha) * _h * _qdot_alpha;
            g = M * dq_tmp - alpha * alpha * _h * _h * f;
            H = M - alpha * alpha * _h * _h * K - alpha * _h * D;
            for (int k = 0;k < _ndof_r;k++) {
                H.col(k) += dM_dq(k) * dq_tmp;
            }

            // save backward_data
            _backward_info._M.push_back(M);
            _backward_info._D.push_back(D);
            _backward_info._dg_du.push_back(-alpha * alpha * _h * _h * df_du);
        }
    }
}

void Simulation::integration_SDIRK2(
        const VectorX q0, const VectorX qdot0, const dtype h, 
        VectorX& q1, VectorX& qdot1) {
    // step 1: 
    // solve for q_alpha
    // compute q_alpha and qdot_alpha by BDF1
    dtype alpha = (2. - sqrt(2.)) / 2.;
    VectorX q_alpha, qdot_alpha;
    
    integration_BDF1(q0, qdot0, alpha * h, q_alpha, qdot_alpha);

    // step 2:
    // compute q1 = q0 + (2alpha-1)qdot0 + 2(1-alpha)qdot_alpha + alpha^2*h^2*inv(M)*f(q1, qdot1)
    // q1_dot = (q1 + (1/a - 2) * q0 - (1 - a) / a * q_alpha) / (a * h)
    
    // save for newton's method
    _q0 = q0; _qdot0 = qdot0;
    _q_alpha = q_alpha; _qdot_alpha = qdot_alpha;
    _h = h;

    q1 = q_alpha + (1 - alpha) * h * qdot_alpha;

    newton(q1, &Simulation::evaluate_g_SDIRK2b, &Simulation::evaluate_g_with_derivatives_SDIRK2b);

    qdot1 = (q1 + (1. / alpha - 2.) * q0 - (1. - alpha) / alpha * q_alpha) / (alpha * h);

    // save backward data: H_inv, M, D, dg_du
    if (_backward_flag) {
        VectorX g;
        MatrixX H;
        auto t_save_backward_start = clock();
        evaluate_g_with_derivatives_SDIRK2b(q1, g, H, true);
        _backward_info._H_lu.push_back(H.partialPivLu());
        _time_report._time_save_backward += clock() - t_save_backward_start;
    }
}

void Simulation::evaluate_g_BDF2(const VectorX& q2, VectorX& g) {
    VectorX qdot2 = 3. / (2. * _h) * (q2 - 4. / 3. * _q1 + 1. / 3. * _q0);

    set_state(q2, qdot2);
    update_robot();

    MatrixX M;
    VectorX f;
    computeMatrices(M, f);

    g = M * (q2 - 4. / 3. * _q1 + 1. / 3. * _q0 - 8. / 9. * _h * _qdot1 + 2. / 9. * _h * _qdot0) - 4. / 9. * _h * _h * f;
}

void Simulation::evaluate_g_with_derivatives_BDF2(const VectorX& q2, VectorX& g, MatrixX& H, bool save_backward_info) {
    VectorX qdot2 = 3. / (2. * _h) * (q2 - 4. / 3. * _q1 + 1. / 3. * _q0);

    set_state(q2, qdot2);
    update_robot(save_backward_info && _backward_design_params_flag);

    if (!save_backward_info) {
        MatrixX M;
        VectorX f;
        JacobianMatrixVector dM_dq;
        MatrixX K, D;
        computeMatrices(M, f, dM_dq, K, D);

        VectorX dq_tmp = q2 - 4. / 3. * _q1 + 1. / 3. * _q0 - 8. / 9. * _h * _qdot1 + 2. / 9. * _h * _qdot0;
        g = M * dq_tmp - 4. / 9. * _h * _h * f;
        H = M - 4. / 9. * _h * _h * K - 2. / 3. * _h * D;
        for (int k = 0;k < _ndof_r;k++) {
            H.col(k) += dM_dq(k) * dq_tmp;
        }
    } else {
        if (_backward_design_params_flag) {
            MatrixX M;
            VectorX f;
            JacobianMatrixVector dM_dq;
            MatrixX K, D;
            MatrixX df_du;
            computeMatrices(M, f, dM_dq, K, D, df_du, _dM_dp, _df_dp);

            VectorX dq_tmp = q2 - 4. / 3. * _q1 + 1. / 3. * _q0 - 8. / 9. * _h * _qdot1 + 2. / 9. * _h * _qdot0;
            g = M * dq_tmp - 4. / 9. * _h * _h * f;
            H = M - 4. / 9. * _h * _h * K - 2. / 3. * _h * D;
            for (int k = 0;k < _ndof_r;k++) {
                H.col(k) += dM_dq(k) * dq_tmp;
            }
            MatrixX dg_dp = -4./9. * _h * _h * _df_dp;
            for (int k = 0;k < _ndof_p1 + _ndof_p2;k++) {
                dg_dp.col(k) += _dM_dp(k) * dq_tmp;
            }
            for (int k = _ndof_p1 + _ndof_p2 + _ndof_p3;k < _ndof_p1 + _ndof_p2 + _ndof_p3 + _ndof_p4;k++) {
                dg_dp.col(k) += _dM_dp(k) * dq_tmp;
            }

            // save backward data
            _backward_info._M.push_back(M);
            _backward_info._D.push_back(D);
            _backward_info._dg_dp.push_back(dg_dp);
            _backward_info._dg_du.push_back(-4./9. * _h * _h * df_du);
        } else {
            MatrixX M;
            VectorX f;
            JacobianMatrixVector dM_dq;
            MatrixX K, D;
            MatrixX df_du;
            computeMatrices(M, f, dM_dq, K, D, df_du);

            VectorX dq_tmp = q2 - 4. / 3. * _q1 + 1. / 3. * _q0 - 8. / 9. * _h * _qdot1 + 2. / 9. * _h * _qdot0;
            g = M * dq_tmp - 4. / 9. * _h * _h * f;
            H = M - 4. / 9. * _h * _h * K - 2. / 3. * _h * D;
            for (int k = 0;k < _ndof_r;k++) {
                H.col(k) += dM_dq(k) * dq_tmp;
            }

            // save backward data
            _backward_info._M.push_back(M);
            _backward_info._D.push_back(D);
            _backward_info._dg_du.push_back(-4./9. * _h * _h * df_du);
        }
    }
}

/*
u(k+1) = 4/3 * u(k) - 1/3 * u(k-1) + 2/3 * h * f(u(k+1))
*/
void Simulation::integration_BDF2(
        const VectorX q0, const VectorX qdot0, const VectorX q1, const VectorX qdot1, const dtype h, 
        VectorX& q2, VectorX& qdot2) {

    // save for newton's method
    _q0 = q0; _qdot0 = qdot0;
    _q1 = q1; _qdot1 = qdot1;
    _h = h;

    q2 = q1 + h * qdot1;

    newton(q2, &Simulation::evaluate_g_BDF2, &Simulation::evaluate_g_with_derivatives_BDF2);

    qdot2 = 3. / (2. * h) * (q2 - 4. / 3. * q1 + 1. / 3. * q0);
    
    // save backward data: H_inv, M, D, dg_du
    if (_backward_flag) {
        VectorX g;
        MatrixX H;
        auto t_save_backward_start = clock();
        evaluate_g_with_derivatives_BDF2(q2, g, H, true);
        _backward_info._H_lu.push_back(H.partialPivLu());
        _time_report._time_save_backward += clock() - t_save_backward_start;
    }
}

// backward related

// NOTE: set the variable indicators and set terminal derivatives in _backward_info
void Simulation::backward() {
    if (_q_his.size() <= 1) {
        std::cerr << "[Error] Please call simulation.forward() before simulation.backward()." << std::endl;
        throw "error";
    }
    
    int T = _q_his.size() - 1;

    // sanity checks
    if (_backward_info._flag_q0) {
        if (_backward_info._df_dq0.size() != _ndof_r) {
            throw_error("_backward_info._df_dq0.size != _ndof_r");
        }
    }
    if (_backward_info._flag_qdot0) {
        if (_backward_info._df_dqdot0.size() != _ndof_r) {
            throw_error("_backward_info._df_dqdot0.size != _ndof_r");
        }
    }
    if (_backward_info._flag_p) {
        if (_backward_info._df_dp.size() != _ndof_p) {
            throw_error("_backward_info._df_dp.size != _ndof_p");
        }
    }
    if (_backward_info._flag_u) {
        if (_backward_info._df_du.size() != _ndof_u * T) {
            throw_error("_backward_info._df_du.size != _ndof_u * T");
        }
    }
    if (_backward_info._df_dq.size() != _ndof_r * T) {
        throw_error("_backward_info._df_dq.size != _ndof_r * T");
    }
    if (_backward_info._df_dvar.size() != _ndof_var * T) {
        throw_error("_backward_info._df_dvar.size != _ndof_var * T");
    }

    auto t_backward_start = clock();

    if (_options->_integrator == "BDF1") {
        backward_BDF1();
    } else if (_options->_integrator == "BDF2") {
        backward_BDF2();
    }

    _time_report._time_backward += clock() - t_backward_start;
}

void Simulation::backward_BDF1() {
    int T = _q_his.size() - 1;

    // initialize
    VectorX& df_dq = _backward_info._df_dq;
    VectorX& df_dvar = _backward_info._df_dvar;
    
    // backward propagation
    VectorX z = VectorX::Zero(T * _ndof_r);
    dtype h = _options->_h;
    
    for (int i = T;i >= 1;i--) {
        int k = i - 1;
        // VectorX yk = df_dq.segment(k * _ndof_r, _ndof_r) + _dvar_dq_his[i].transpose() * df_dvar.segment(k * _ndof_var, _ndof_var);
        VectorX yk = df_dq.segment(k * _ndof_r, _ndof_r) + _backward_info._dvar_dq[k].transpose() * df_dvar.segment(k * _ndof_var, _ndof_var);
        if (k < T - 1) {
            // Add contributions from step k + 1
            MatrixX M = _backward_info._M[k + 1];
            MatrixX D = _backward_info._D[k + 1];
            MatrixX H = -2. * M + h * D;
            yk -= H.transpose() * z.segment((k + 1) * _ndof_r, _ndof_r);
        }
        if (k < T - 2) {
            // Add contributions from step k + 2
            MatrixX& H = _backward_info._M[k + 2];
            yk -= H.transpose() * z.segment((k + 2) * _ndof_r, _ndof_r);
        }
        // z.segment(k * _ndof_r, _ndof_r) = _backward_info._H_inv[k].transpose() * yk;
        z.segment(k * _ndof_r, _ndof_r) = _backward_info._H_lu[k].transpose().solve(yk);
    }

    if (_backward_info._flag_q0) {
        // df_dq0 = _df_dq0 - dg_dq0' * z
        _backward_results._df_dq0 = _backward_info._df_dq0;

        // dg(1)_dq0 = -M(1) + hD(1)
        if (T > 0) {
            MatrixX dg1_dq0 = -_backward_info._M[0] + h * _backward_info._D[0];
            _backward_results._df_dq0 -= dg1_dq0.transpose() * z.head(_ndof_r);
        }

        // dg(2)_dq0 = M(2)
        if (T > 1) {
            MatrixX& dg2_dq0 = _backward_info._M[1];
            _backward_results._df_dq0 -= dg2_dq0.transpose() * z.segment(_ndof_r, _ndof_r);
        }
    }

    if (_backward_info._flag_qdot0) {
        // df_dqdot0 = _df_dqdot0 - dg_dqdot0' * z
        _backward_results._df_dqdot0 = _backward_info._df_dqdot0;
        
        // dg(1)_dqdot0 = -hM(1)
        if (T > 0) {
            MatrixX dg1_dqdot0 = -h * _backward_info._M[0];
            _backward_results._df_dqdot0 -= dg1_dqdot0.transpose() * z.head(_ndof_r);
        } 
    }
    
    if (_backward_info._flag_p) { 
        // df_dp = _df_dp - dg_dp' * z
        _backward_results._df_dp = _backward_info._df_dp;
        for (int k = 0;k < T;k++)
            _backward_results._df_dp += 
                - _backward_info._dg_dp[k].transpose() * z.segment(k * _ndof_r, _ndof_r)
                + _backward_info._dvar_dp[k].transpose() * _backward_info._df_dvar.segment(k * _ndof_var, _ndof_var);
    }

    if (_backward_info._flag_u) { 
        // df_du = _df_du - dg_du' * z
        _backward_results._df_du = _backward_info._df_du;
        for (int k = 0;k < T;k++)
            _backward_results._df_du.segment(k * _ndof_u, _ndof_u) -= 
                _backward_info._dg_du[k].transpose() * z.segment(k * _ndof_r, _ndof_r);
    }
}

void Simulation::backward_BDF2() {
    int T = _q_his.size() - 1;

    dtype alpha = (2. - sqrt(2.)) / 2.;

    // initialize
    VectorX& df_dq = _backward_info._df_dq;
    VectorX& df_dvar = _backward_info._df_dvar;
    
    // backward propagation
    VectorX z = VectorX::Zero((T + 1) * _ndof_r);
    dtype h = _options->_h;
    for (int k = T;k >= 1;k--) {
        // VectorX yk = df_dq.segment((k - 1) * _ndof_r, _ndof_r) + _dvar_dq_his[k].transpose() * df_dvar.segment((k - 1) * _ndof_var, _ndof_var);
        VectorX yk = df_dq.segment((k - 1) * _ndof_r, _ndof_r) + _backward_info._dvar_dq[k - 1].transpose() * df_dvar.segment((k - 1) * _ndof_var, _ndof_var);
        if (k <= T - 1) {
            // Add contributions from step k + 1
            MatrixX M = _backward_info._M[k + 1];
            MatrixX D = _backward_info._D[k + 1];
            if (k > 1) { // BDF2 step    
                MatrixX H = -8. / 3. * M + 8. / 9. * h * D;
                yk -= H.transpose() * z.segment((k + 1) * _ndof_r, _ndof_r);
            } else {     // SDIRK2 step
                MatrixX H = -(8. / (9. * alpha) + 4. / 3.) * M + 8. / 9. * h * D;
                yk -= H.transpose() * z.segment((k + 1) * _ndof_r, _ndof_r);
            }
        }

        if (k <= T - 2) {
            // Add contributions from step k + 2
            MatrixX M = _backward_info._M[k + 2];
            MatrixX D = _backward_info._D[k + 2];
            if (k > 1) { // BDF2 step    
                MatrixX H = 22. / 9. * M - 2. / 9. * h * D;
                yk -= H.transpose() * z.segment((k + 2) * _ndof_r, _ndof_r);
            } else {     // SDIRK2 step
                MatrixX H = (2. / (9. * alpha) + 19. / 9.) * M - 2. / 9. * h * D;
                yk -= H.transpose() * z.segment((k + 2) * _ndof_r, _ndof_r);
            }
        }

        if (k <= T - 3) {
            // Add contributions from step k + 3
            MatrixX M = _backward_info._M[k + 3];
            MatrixX H = - 8. / 9. * M;
            yk -= H.transpose() * z.segment((k + 3) * _ndof_r, _ndof_r);
        }

        if (k <= T - 4) {
            // Add contributions from step k + 4
            MatrixX M = _backward_info._M[k + 4];
            MatrixX H = 1. / 9. * M;
            yk -= H.transpose() * z.segment((k + 4) * _ndof_r, _ndof_r);
        }

        z.segment(k * _ndof_r, _ndof_r) = _backward_info._H_lu[k].transpose().solve(yk);
    }

    // SDIRK2 (a) step
    VectorX ya = VectorX::Zero(_ndof_r); // qa is an internal state, so df_dya = 0

    if (T >= 1) {
        MatrixX M = _backward_info._M[1];
        MatrixX D = _backward_info._D[1];
        MatrixX H = 2. * (alpha - 1.) / alpha * M + h * (1. - alpha) * D;
        ya -= H.transpose() * z.segment(_ndof_r, _ndof_r);
    }

    if (T >= 2) {
        MatrixX M = _backward_info._M[2];
        MatrixX H = 8. * (1. - alpha) / (9. * alpha * alpha) * M;
        ya -= H.transpose() * z.segment(_ndof_r * 2, _ndof_r);
    }

    if (T >= 3) {
        MatrixX M = _backward_info._M[3];
        MatrixX H = 2. * (alpha - 1.) / (9. * alpha * alpha) * M;
        ya -= H.transpose() * z.segment(_ndof_r * 3, _ndof_r);
    }

    // z.head(_ndof_r) = _backward_info._H_inv[0].transpose() * ya;
    z.head(_ndof_r) = _backward_info._H_lu[0].transpose().solve(ya);

    // aggregate everything together
    if (_backward_info._flag_q0) {
        // df_dq0 = _df_dq0 - dg_dq0' * z
        _backward_results._df_dq0 = _backward_info._df_dq0;

        if (T >= 1) {
            // dga_dq0 = -M + a * h * D
            MatrixX dga_dq0 = -_backward_info._M[0] + alpha * h * _backward_info._D[0];
            _backward_results._df_dq0 -= dga_dq0.transpose() * z.head(_ndof_r);
            // dg1_dq0 = (2-3a) / a * M - a * h * (1/a + 2) * D
            MatrixX dg1_dq0 = (2. - 3. * alpha) / alpha * _backward_info._M[1] - alpha * h * (1. / alpha - 2.) * _backward_info._D[1];
            _backward_results._df_dq0 -= dg1_dq0.transpose() * z.segment(_ndof_r, _ndof_r);
        }

        if (T >= 2) {
            // dg2_dq0 = (3a^2 + 16a - 8) / (9a^2) * M - 2/9 * h * D
            MatrixX dg2_dq0 = (3. * alpha * alpha + 16. * alpha - 8.) / (9. * alpha * alpha) * _backward_info._M[2]
                                - 2. / 9. * h * _backward_info._D[2];
            _backward_results._df_dq0 -= dg2_dq0.transpose() * z.segment(_ndof_r * 2, _ndof_r);
        }

        if (T >= 3) {
            // dg3_dq0 = -(4a^2 + 4a - 2) / (9a^2) * M
            MatrixX dg3_dq0 = - (4. * alpha * alpha + 4. * alpha - 2.) / (9. * alpha * alpha) * _backward_info._M[3];
            _backward_results._df_dq0 -= dg3_dq0.transpose() * z.segment(_ndof_r * 3, _ndof_r);
        }

        if (T >= 4) {
            // dg4_dq0 = 1/9 * M
            MatrixX dg4_dq0 = _backward_info._M[4] / 9.;
            _backward_results._df_dq0 -= dg4_dq0.transpose() * z.segment(_ndof_r * 4, _ndof_r);
        }
    }

    if (_backward_info._flag_qdot0) {
        // df_dqdot0 = _df_dqdot0 - dg_dqdot0' * z
        _backward_results._df_dqdot0 = _backward_info._df_dqdot0;
        
        if (T >= 1) {
            // dga_dqdot0 = -a * h * M
            MatrixX dga_dqdot0 = -alpha * h * _backward_info._M[0];
            _backward_results._df_dqdot0 -= dga_dqdot0.transpose() * z.head(_ndof_r);
            // dg1_dqdot0 = -(2a - 1) * h * M
            MatrixX dg1_dqdot0 = -(2. * alpha - 1.) * h * _backward_info._M[1];
            _backward_results._df_dqdot0 -= dg1_dqdot0.transpose() * z.segment(_ndof_r, _ndof_r);
        }

        if (T >= 2) {
            // dg2_dqdot0 = 2/9 * h * M
            MatrixX dg2_dqdot0 = 2. / 9. * h * _backward_info._M[2];
            _backward_results._df_dqdot0 -= dg2_dqdot0.transpose() * z.segment(_ndof_r * 2, _ndof_r);
        }
    }

    if (_backward_info._flag_p) {
        // df_dp = _df_dp - dg_dp' * z
        _backward_results._df_dp = _backward_info._df_dp;
        for (int k = 0;k <= T;k++)
            _backward_results._df_dp -= 
                _backward_info._dg_dp[k].transpose() * z.segment(k * _ndof_r, _ndof_r);
        for (int k = 0;k < T;k++)
            _backward_results._df_dp += 
                _backward_info._dvar_dp[k].transpose() * _backward_info._df_dvar.segment(k * _ndof_var, _ndof_var);
    }

    if (_backward_info._flag_u) {
        // df_du = _df_du - dg_du' * z
        _backward_results._df_du = _backward_info._df_du;
        for (int k = 0;k < T;k++)
            _backward_results._df_du.segment(k * _ndof_u, _ndof_u) -= 
                _backward_info._dg_du[k + 1].transpose() * z.segment((k + 1) * _ndof_r, _ndof_r);
        _backward_results._df_du.head(_ndof_u) -= _backward_info._dg_du[0].transpose() * z.head(_ndof_r);
    }
}

// inverse dynamics

VectorX Simulation::inverse_dynamics_BDF1(const VectorX q0, const VectorX qdot0, const VectorX q1, bool reset) {

    // current simulation states
    VectorX q_curr, qdot_curr;
    std::vector<Vector6> f_external_curr;
    if (reset) {
        q_curr = get_q();
        qdot_curr = get_qdot();
        f_external_curr = _robot->get_external_force();
    }

    VectorX qdot1 = (q1 - q0) / _h;
    
    set_state(q1, qdot1);
    update_robot();

    MatrixX J, Jdot, MmJ(_ndof_m, _ndof_r), MmJdot(_ndof_m, _ndof_r);
    VectorX Mm, fm, fr;
    _robot->reset_external_force();
    _robot->computeJointJacobian(J, Jdot);    
    _robot->computeMaximalMassMatrix(Mm);
    _robot->computeForce(fm, fr);
    for (int i = 0;i < MmJ.rows();i++) 
        MmJ.row(i).noalias() = J.row(i) * Mm(i);
    for (int i = 0;i < MmJdot.rows();i++)
        MmJdot.row(i).noalias() = Jdot.row(i) * Mm(i);

    // compute body-frame external force
    VectorX f_maximal = fm - MmJdot * qdot1 + J * fr;
    VectorX f_external = MmJ * (qdot1 - qdot0) / _h - f_maximal;

    // compute world-frame external force
    for (auto body : _robot->_bodies) {
        Matrix3 R_0i = body->_E_0i.topLeftCorner(3, 3);
        Vector6 f_body_external = f_external.segment(body->_index[0], 6);
        f_external.segment(body->_index[0], 3) = R_0i * f_body_external.head(3);
        f_external.segment(body->_index[0] + 3, 3) = R_0i * f_body_external.tail(3);
    }

    if (reset) {
        set_state(q_curr, qdot_curr);
        update_robot();
        _robot->set_external_force(f_external_curr);
    }

    return f_external;
}

VectorX Simulation::inverse_dynamics_BDF2(const VectorX q0, const VectorX qdot0, const VectorX q1, const VectorX qdot1, const VectorX q2, bool reset) {

    // current simulation states
    VectorX q_curr, qdot_curr;
    std::vector<Vector6> f_external_curr;
    if (reset) {
        q_curr = get_q();
        qdot_curr = get_qdot();
        f_external_curr = _robot->get_external_force();
    }
    
    VectorX qdot2 = 3. / (2. * _h) * (q2 - 4. / 3. * q1 + 1. / 3. * q0);
    
    set_state(q2, qdot2);
    update_robot();

    MatrixX J, Jdot, MmJ(_ndof_m, _ndof_r), MmJdot(_ndof_m, _ndof_r);
    VectorX Mm, fm, fr;
    _robot->reset_external_force();
    _robot->computeJointJacobian(J, Jdot);    
    _robot->computeMaximalMassMatrix(Mm);
    _robot->computeForce(fm, fr);
    for (int i = 0;i < MmJ.rows();i++) 
        MmJ.row(i).noalias() = J.row(i) * Mm(i);
    for (int i = 0;i < MmJdot.rows();i++)
        MmJdot.row(i).noalias() = Jdot.row(i) * Mm(i);

    // compute body-frame external force
    VectorX f_maximal = fm - MmJdot * qdot1 + J * fr;
    VectorX f_external = MmJ * (qdot2 - 4. / 3. * qdot1 + 1. / 3. * qdot0) / (2. / 3. * _h) - f_maximal;

    // compute world-frame external force
    for (auto body : _robot->_bodies) {
        Matrix3 R_0i = body->_E_0i.topLeftCorner(3, 3);
        Vector6 f_body_external = f_external.segment(body->_index[0], 6);
        f_external.segment(body->_index[0], 3) = R_0i * f_body_external.head(3);
        f_external.segment(body->_index[0] + 3, 3) = R_0i * f_body_external.tail(3);
    }

    if (reset) {
        set_state(q_curr, qdot_curr);
        update_robot();
        _robot->set_external_force(f_external_curr);
    }

    return f_external;
}

std::vector<VectorX> Simulation::inverse_dynamics_his(bool forward, bool reset) {

    if (_options->_integrator == "BDF1") {
        return inverse_dynamics_his_BDF1(forward, reset);
    } else if (_options->_integrator == "BDF2") {
        return inverse_dynamics_his_BDF2(forward, reset);
    } else {
        std::cerr << "[Error] Inverse dynamics of integrator " << _options->_integrator << " has not been implemented." << std::endl;
        throw "error";
    }
}

std::vector<VectorX> Simulation::inverse_dynamics_his_BDF1(bool forward, bool reset) {
    
    if (_q_his.size() <= 1) {
        std::cerr << "[Error] Please call simulation.forward() before computing inverse dynamics of the history." << std::endl;
        throw "error";
    }

    // current simulation states
    VectorX q_curr, qdot_curr;
    std::vector<Vector6> f_external_curr;
    if (reset) {
        q_curr = get_q();
        qdot_curr = get_qdot();
        f_external_curr = _robot->get_external_force();
    }

    std::vector<VectorX> f_externals;
    auto his_size = _q_his.size();

    VectorX q0, qdot0, q1;
    if (forward) qdot0 = _qdot_his[0];
    // else qdot0 = -_qdot_his[his_size - 1]; // use reverse velocity to init
    else qdot0 = VectorX::Zero(_qdot_his[his_size - 1].size()); // use zero velocity to init

    for (int i = 0; i < his_size - 1; ++i) {

        if (forward) q0 = _q_his[i], q1 = _q_his[i + 1];
        else q0 = _q_his[his_size - 1 - i], q1 = _q_his[his_size - 2 - i];

        VectorX f_external = inverse_dynamics_BDF1(q0, qdot0, q1, false);
        f_externals.push_back(f_external);

        qdot0 = (q1 - q0) / _h;
    }

    if (reset) {
        set_state(q_curr, qdot_curr);
        update_robot();
        _robot->set_external_force(f_external_curr);
    }

    return f_externals;
}

std::vector<VectorX> Simulation::inverse_dynamics_his_BDF2(bool forward, bool reset) {

    if (_q_his.size() <= 1) {
        std::cerr << "[Error] Please call simulation.forward() before computing inverse dynamics of the history." << std::endl;
        throw "error";
    }

    // current simulation states
    VectorX q_curr, qdot_curr;
    std::vector<Vector6> f_external_curr;
    if (reset) {
        q_curr = get_q();
        qdot_curr = get_qdot();
        f_external_curr = _robot->get_external_force();
    }

    std::vector<VectorX> f_externals;
    auto his_size = _q_his.size();

    VectorX q0, qdot0, q1, qdot1, q2;
    if (forward) qdot0 = _qdot_his[0];
    // else qdot0 = -_qdot_his[his_size - 1]; // use reverse velocity to init
    else qdot0 = VectorX::Zero(_qdot_his[his_size - 1].size()); // use zero velocity to init

    for (int i = 0; i < his_size - 1; ++i) {

        VectorX f_external;

        if (i == 0) { // BDF1

            if (forward) q0 = _q_his[i], q1 = _q_his[i + 1];
            else q0 = _q_his[his_size - 1 - i], q1 = _q_his[his_size - 2 - i];
            f_external = inverse_dynamics_BDF1(q0, qdot0, q1, false);
            qdot1 = (q1 - q0) / _h;

        } else { // BDF2

            if (forward) q0 = _q_his[i - 1], q1 = _q_his[i], q2 = _q_his[i + 1];
            else q0 = _q_his[his_size - i], q1 = _q_his[his_size - 1 - i], q2 = _q_his[his_size - 2 - i];
            f_external = inverse_dynamics_BDF2(q0, qdot0, q1, qdot1, q2, false);
            qdot0 = qdot1;
            qdot1 = 3. / (2. * _h) * (q2 - 4. / 3. * q1 + 1. / 3. * q0);
        }

        f_externals.push_back(f_external);
    }

    if (reset) {
        set_state(q_curr, qdot_curr);
        update_robot();
        _robot->set_external_force(f_external_curr);
    }

    return f_externals;
}

// viewer related

void Simulation::get_rendering_objects(
    std::vector<Matrix3Xf>& vertex_list, 
    std::vector<Matrix3Xi>& face_list,
    std::vector<opengl_viewer::Option>& option_list,
    std::vector<opengl_viewer::Animator*>& animator_list) {
    
    _robot->get_rendering_objects(vertex_list, face_list, option_list, animator_list);
}

void Simulation::init_viewer() {
    if (!_viewer) {
        _viewer = std::make_shared<SimViewer>(this);
    }

    _viewer->initialize();

    _viewer_step = 0;
    set_state(_q_his[_viewer_step], _qdot_his[_viewer_step]);
    update_robot();
}

// advance the viewer steps, and return whether it completes one trajectory.
bool Simulation::advance_viewer_step(int num_steps) {
    _viewer_step += num_steps;

    bool done = false;

    if (_viewer_step >= _q_his.size()) {
        done = true;
        if (_viewer_options->_loop) {
            _viewer_step %= _q_his.size();
        } else {
            _viewer_step = _q_his.size() - 1;
        }
    }

    set_state(_q_his[_viewer_step], _qdot_his[_viewer_step]);
    update_robot();

    return done;
}

void Simulation::replay() {
    if (_q_his.size() == 0) {
        std::cerr << "[Error] Please call simulation.reset() before simulation.run_viewer()." << std::endl;
        throw "error";
    }
    
    // save the current q and qdot
    VectorX q = get_q();
    VectorX qdot = get_qdot();

    init_viewer();

    _viewer->run();

    // restore the q and qdot
    set_state(q, qdot);
    update_robot();
}

void Simulation::export_replay(std::string folder) {
    for (int i = 0;i < _q_his.size();i++) {
        set_state(_q_his[i], _qdot_his[i]);
        update_robot();
        std::string path = folder + "/" + to_string(i) + ".txt";
        FILE* fp = fopen(path.c_str(), "w");
        fprintf(fp, "%d\n", _robot->_bodies.size());
        for (auto body : _robot->_bodies) {
            for (int j = 0;j < 4;j++) {
                for (int k = 0;k < 4;k++) {
                    fprintf(fp, "%.6lf ", body->_E_0i(j, k));
                }
                fprintf(fp, "\n");
            }
        }
        fclose(fp);
    }
}

void Simulation::print_time_report() {
    std::cerr << "----------- Time Report -----------" << std::endl;
    std::cerr << "|Simulation                       |" << std::endl;
    std::cerr << "|---------------------------------|" << std::endl;
    std::cerr << "|" << std::setw(20) << std::left << "Solver" << "|" << std::setw(7) << std::right << _time_report._time_solver / 1000 << "(ms) |" << std::endl;
    std::cerr << "|" << std::setw(20) << std::left << "Compute Matrices" << "|" << std::setw(7) << std::right << _time_report._time_compute_matrices / 1000 << "(ms) |" << std::endl;
    std::cerr << "|" << std::setw(20) << std::left << "-- Compute dJ" << "|" << std::setw(7) << std::right << _time_report._time_compute_dJ / 1000 << "(ms) |" << std::endl;
    std::cerr << "|" << std::setw(20) << std::left << "-- Compute df" << "|" << std::setw(7) << std::right << _time_report._time_compute_df / 1000 << "(ms) |" << std::endl;
    std::cerr << "|" << std::setw(20) << std::left << "Compose Matrices" << "|" << std::setw(7) << std::right << _time_report._time_compose_matrices / 1000 << "(ms) |" << std::endl;
    std::cerr << "|" << std::setw(20) << std::left << "-- Compose dM_dp" << "|" << std::setw(7) << std::right << _time_report._time_dM_dp / 1000 << "(ms) |" << std::endl;
    std::cerr << "|" << std::setw(20) << std::left << "---- Compose dM_dp1" << "|" << std::setw(7) << std::right << _time_report._time_dM_dp1 / 1000 << "(ms) |" << std::endl;
    std::cerr << "|" << std::setw(20) << std::left << "---- Compose dM_dp2" << "|" << std::setw(7) << std::right << _time_report._time_dM_dp2 / 1000 << "(ms) |" << std::endl;
    std::cerr << "|" << std::setw(20) << std::left << "---- Compose dM_dp4" << "|" << std::setw(7) << std::right << _time_report._time_dM_dp4 / 1000 << "(ms) |" << std::endl;
    std::cerr << "|" << std::setw(20) << std::left << "-- Compose df_dp" << "|" << std::setw(7) << std::right << _time_report._time_df_dp / 1000 << "(ms) |" << std::endl;
    std::cerr << "|" << std::setw(20) << std::left << "Save Backward" << "|" << std::setw(7) << std::right << _time_report._time_save_backward / 1000 << "(ms) |" << std::endl;
    std::cerr << "|" << std::setw(20) << std::left << "Backward" << "|" << std::setw(7) << std::right << _time_report._time_backward / 1000 << "(ms) |" << std::endl;
    // _robot->print_time_report();
    std::cerr << "-----------------------------------" << std::endl;
}

}