#include "Robot.h"
#include "Joint/Joint.h"
#include "Body/Body.h"
#include "Force/Force.h"
#include "Actuator/Actuator.h"
#include "EndEffector/EndEffector.h"
#include "VirtualObject/VirtualObject.h"
#include "Joint/JointDesignParameters.h"
#include "Body/BodyDesignParameters.h"
#include "Body/BodyAbstract.h"
#include "Body/BodyPrimitiveShape.h"
#include "Body/BodySDFObj.h"
#include "Body/BodyBVHObj.h"
#include "Force/ForceGeneralPrimitiveContact.h"
#include "Force/ForceGeneralSDFContact.h"
#include "Force/ForceGeneralBVHContact.h"
#include "Sensor/Tactile.h"

namespace redmax {

void Robot::add_tactile(Tactile* tactile) {
    _tactiles.push_back(tactile);
}

void Robot::add_force(Force* force) {
    _forces.push_back(force);
}

void Robot::add_actuator(Actuator* actuator) {
    _actuators.push_back(actuator);
}

void Robot::add_end_effector(EndEffector* end_effector) {
    _end_effectors.push_back(end_effector);
}

void Robot::add_virtual_object(VirtualObject* virtual_object) {
    _virtual_objects.push_back(virtual_object);
}

// init robot
void Robot::construct_dfs_order(Joint* now) {
    _joints.push_back(now);
    for (auto child : now->_children)
        construct_dfs_order(child);
}

void Robot::init(bool verbose) {
    assert(_root_joints.size() > 0);

    // make array for joints and bodies in bfs/dfs order
    _joints.clear();

    // // bfs order
    // int l1 = 0, l2 = 0;
    // for (Joint* root_joint : _root_joints) {
    //     _joints.push_back(root_joint);
    //     l2 += 1;
    //     for (;l1 < l2;++l1) {
    //         for (auto joint : _joints[l1]->_children) {
    //             _joints.push_back(joint);
    //             ++l2;
    //         }
    //     }
    // }

    // dfs order
    for (Joint* root_joint : _root_joints) {
        construct_dfs_order(root_joint);
    }

    _bodies.clear();
    for (auto joint : _joints) {
        if (joint->_body != nullptr)
            _bodies.push_back(joint->_body);
    }

    if (verbose) {
        // print order for check
        cerr << "joint order : ";
        for (auto joint : _joints) {
            cerr << joint->_id << " ";
        }
        cerr << endl;
    }

    // count dofs, assign index for each joints' dof in reduced matrix and 
    // index for body's dof in maximal matrix
    _ndof_r = 0;
    for (auto joint : _joints) {
        joint->_index.clear();
        for (int i = 0;i < joint->_ndof;++i)
            joint->_index.push_back(_ndof_r++);
    }

    _ndof_m = 0;
    for (auto body : _bodies) {
        body->_index.clear();
        for (int i = 0;i < 6;++i)
            body->_index.push_back(_ndof_m++);
    }

    _ndof_u = 0;
    for (auto actuator : _actuators) {
        actuator->_index.clear();
        for (int i = 0;i < actuator->_ndof;i++) {
            actuator->_index.push_back(_ndof_u++);
        }
    }

    _ndof_var = 0;
    for (auto end_effector : _end_effectors) {
        end_effector->_index.clear();
        for (int i = 0;i < end_effector->_ndof;i++) {
            end_effector->_index.push_back(_ndof_var++);
        }
    }

    // count design dofs
    _ndof_p1 = _ndof_p2 = _ndof_p3 = _ndof_p4 = _ndof_p5 = _ndof_p6 = 0;
    for (auto joint : _joints) {
        if (joint->_design_params_1._active)
            for (int i = 0;i < joint->_design_params_1._ndof;i++)
                joint->_design_params_1._param_index(i) = _ndof_p1++;
        if (joint->_design_params_5._active)
            for (int i = 0;i < joint->_design_params_5._ndof;i++)
                joint->_design_params_5._param_index(i) = _ndof_p5++;
    }
    for (auto body : _bodies) {
        if (body->_design_params_2._active)
            for (int i = 0;i < body->_design_params_2._ndof;i++)
                body->_design_params_2._param_index(i) = _ndof_p2++;
        if (body->_design_params_3._active)
            for (int i = 0;i < body->_design_params_3._ndof;i++)
                body->_design_params_3._param_index(i) = _ndof_p3++;
        if (body->_design_params_4._active)
            for (int i = 0;i < body->_design_params_4._ndof;i++)
                body->_design_params_4._param_index(i) = _ndof_p4++;
        if (body->_design_params_6._active)
            for (int i = 0;i < body->_design_params_6._ndof;i++)
                body->_design_params_6._param_index(i) = _ndof_p6++;
    }
    _ndof_p = _ndof_p1 + _ndof_p2 + _ndof_p3 + _ndof_p4 + _ndof_p5 + _ndof_p6;

    // init joints and bodies and springs
    for (auto joint : _joints) {
        joint->init();
    }
    for (auto body : _bodies) {
        body->init();
    }

    // update joints and bodies to fill the data
    update();

    // init force after update
    for (auto force : _forces) {
        force->init();
    }

    // init dJ_dp1, dJ_dp2, dJdot_dp1, dJdot_dp2
    dJ_dp1 = JacobianMatrixVector(_ndof_m, _ndof_r, _ndof_p1);
    dJ_dp2 = JacobianMatrixVector(_ndof_m, _ndof_r, 12);
    dJdot_dp1 = JacobianMatrixVector(_ndof_m, _ndof_r, _ndof_p1);
    dJdot_dp2 = JacobianMatrixVector(_ndof_m, _ndof_r, 12);

    // init dfm_dp
    dfm_dp = MatrixX::Zero(_ndof_m, _ndof_p);
    dfr_dp = MatrixX::Zero(_ndof_r, _ndof_p);
}

void Robot::update(bool design_gradient) {
    for (Joint* joint : _root_joints)
        joint->update(design_gradient);
}

bool Robot::check_valid() {
    for (auto joint : _joints) {
        if (!joint->check_valid()) 
            return false;
    }
    return true;
}

void Robot::reparam() {
    for (auto joint : _joints) {
        joint->reparam();
    }
    
    update();
}

void Robot::set_q(const VectorX q) {
    assert(q.size() == _ndof_r);
    for (auto joint : _joints) {
        for (int i = 0;i < joint->_ndof;++i) {
            joint->_q[i] = q[joint->_index[i]];
        }
    }
}

void Robot::set_qdot(const VectorX qdot) {
    assert(qdot.size() == _ndof_r);
    for (auto joint : _joints) {
        for (int i = 0;i < joint->_ndof;++i) {
            joint->_qdot[i] = qdot[joint->_index[i]];
        }
    }
}

void Robot::set_state(const VectorX q, const VectorX qdot) {
    set_q(q);
    set_qdot(qdot);
}

VectorX Robot::get_q() {
    VectorX q = VectorX::Zero(_ndof_r);
    for (auto joint : _joints) {
        for (int i = 0;i < joint->_ndof;++i) {
            q[joint->_index[i]] = joint->_q[i];
        }
    }
    return q;
}

VectorX Robot::get_qdot() {
    VectorX qdot = VectorX::Zero(_ndof_r);
    for (auto joint : _joints) {
        for (int i = 0;i < joint->_ndof;++i) {
            qdot[joint->_index[i]] = joint->_qdot[i];
        }
    }
    return qdot;
}

VectorX Robot::get_phi() {
    VectorX phi = VectorX::Zero(_ndof_m);
    for (auto body : _bodies) {
        for (int i = 0;i < 6;i++) {
            phi[body->_index[i]] = body->_phi[i];
        }
    }
    return phi;
}

VectorX Robot::get_variables() {
    VectorX variables = VectorX::Zero(_ndof_var);
    for (auto end_effector : _end_effectors) {
        end_effector->computeVariables(variables);
    }
    return variables;
}

// control variables
void Robot::update_actuator_dofs(const VectorX& q_next, const VectorX& qdot_next) {
    for (auto actuator : _actuators) {
        actuator->update_dofs(q_next, qdot_next);
    }
}

VectorX Robot::get_u() {
    VectorX u = VectorX::Zero(_ndof_u);
    for (auto actuator : _actuators) {
        for (int i = 0;i < actuator->_ndof;i++)
            u[actuator->_index[i]] = actuator->_u[i];
    }
    return u;
}

void Robot::set_u(const VectorX& u) {
    assert(u.size() == _ndof_u);
    for (auto actuator : _actuators) {
        actuator->set_u(u);
    }
}

void Robot::get_ctrl_range(VectorX& ctrl_min, VectorX& ctrl_max) {
    ctrl_min = VectorX::Zero(_ndof_u);
    ctrl_max = VectorX::Zero(_ndof_u);
    for (auto actuator : _actuators) {
        actuator->get_ctrl_range(ctrl_min, ctrl_max);
    }
}

void Robot::print_ctrl_info() {
    std::cerr << " ------------- Control Info ------------ " << std::endl;
    std::cerr << "|             Name             |  #Dof  |" << std::endl;
    std::cerr << "|------------------------------|--------|" << std::endl;
    for (auto actuator : _actuators) {
        std::cerr << "|" << std::setw(30) << std::left << actuator->_name << "|" << std::setw(8) << std::right << actuator->_ndof << "|" << std::endl;
    }
    std::cerr << " --------------------------------------- " << std::endl;
}

// external force
std::vector<Vector6> Robot::get_external_force() const {
    std::vector<Vector6> forces;
    for (auto body : _bodies) {
        forces.push_back(body->get_external_force());
    }
    return forces;
}

void Robot::set_external_force(std::vector<Vector6> forces) {
    if (forces.size() != _bodies.size()) {
        throw_error("Size of external forces does not match size of bodies");
    }
    for (int i = 0; i < _bodies.size(); ++i) {
        _bodies[i]->set_external_force(forces[i]);
    }
}

void Robot::reset_external_force() {
    for (auto body : _bodies) {
        body->reset_external_force();
    }
}

// tactile data
std::vector<dtype> Robot::get_tactile_depth(std::string name) {
    for (auto tactile : _tactiles) {
        if (tactile->_name == name) {
            tactile->compute_tactile_values();
            return tactile->_depth;
        }
    }
    return std::vector<dtype>();
}

std::vector<Vector2i> Robot::get_tactile_image_pos(std::string name) {
    for (auto tactile : _tactiles) {
        if (tactile->_name == name) {
            return tactile->_image_pos;
        }
    }
    return std::vector<Vector2i>();
}

// design parameters
VectorX Robot::get_design_params() {
    VectorX design_params_1 = VectorX::Zero(_ndof_p1);
    VectorX design_params_5 = VectorX::Zero(_ndof_p5);
    for (auto joint : _joints) {
        joint->_design_params_1.get_params(design_params_1);
        joint->_design_params_5.get_params(design_params_5);
    }

    VectorX design_params_2 = VectorX::Zero(_ndof_p2);
    VectorX design_params_3 = VectorX::Zero(_ndof_p3);
    VectorX design_params_4 = VectorX::Zero(_ndof_p4);
    VectorX design_params_6 = VectorX::Zero(_ndof_p6);
    for (auto body : _bodies) {
        body->_design_params_2.get_params(design_params_2);
        body->_design_params_3.get_params(design_params_3);
        body->_design_params_4.get_params(design_params_4);
        body->_design_params_6.get_params(design_params_6);
    }

    VectorX design_params = VectorX::Zero(_ndof_p1 + _ndof_p2 + _ndof_p3 + _ndof_p4 + _ndof_p5 + _ndof_p6);
    design_params << design_params_1, design_params_2, design_params_3, design_params_4, design_params_5, design_params_6;

    return design_params;
}

void Robot::set_design_params(const VectorX &design_params) {
    if (_ndof_p1 > 0) {
        VectorX design_params_1 = design_params.segment(0, _ndof_p1);
        for (auto joint : _joints) {
            joint->_design_params_1.update_params(design_params_1);
        }
    }
    if (_ndof_p2 > 0) {
        VectorX design_params_2 = design_params.segment(_ndof_p1, _ndof_p2);
        for (auto body : _bodies) {
            body->_design_params_2.update_params(design_params_2);
        }
    }
    if (_ndof_p3 > 0) {
        VectorX design_params_3 = design_params.segment(_ndof_p1 + _ndof_p2, _ndof_p3);
        for (auto body : _bodies) {
            body->_design_params_3.update_params(design_params_3);
        }
    }
    if (_ndof_p4 > 0) {
        VectorX design_params_4 = design_params.segment(_ndof_p1 + _ndof_p2 + _ndof_p3, _ndof_p4);
        for (auto body : _bodies) {
            body->_design_params_4.update_params(design_params_4);
        }
    }
    if (_ndof_p5 > 0) {
        VectorX design_params_5 = design_params.segment(_ndof_p1 + _ndof_p2 + _ndof_p3 + _ndof_p4, _ndof_p5);
        for (auto joint : _joints) {
            joint->_design_params_5.update_params(design_params_5);
        }
    }
    if (_ndof_p6 > 0) {
        VectorX design_params_6 = design_params.segment(_ndof_p1 + _ndof_p2 + _ndof_p3 + _ndof_p4 + _ndof_p5, _ndof_p6);
        for (auto body : _bodies) {
            body->_design_params_6.update_params(design_params_6);
        }
    }
}

void Robot::print_design_params_info() {
    std::cerr << " ---------- Design Params Info ---------- " << std::endl;
    if (_ndof_p1 > 0) {
        std::cerr << "|            Type I            |  #Dof  |" << std::endl;
        std::cerr << "|------------------------------|--------|" << std::endl;
        for (auto joint : _joints) {
            if (joint->_design_params_1._active)
                std::cerr << "|" << std::setw(30) << std::left << joint->_name << "|" << std::setw(8) << std::right << joint->_design_params_1._ndof << "|" << std::endl;
        }
        std::cerr << "|---------------------------------------| " << std::endl;
    }
    if (_ndof_p2 > 0) {
        std::cerr << "|            Type II           |  #Dof  |" << std::endl;
        std::cerr << "|------------------------------|--------|" << std::endl;
        for (auto body : _bodies) {
            if (body->_design_params_2._active)
                std::cerr << "|" << std::setw(30) << std::left << body->_name << "|" << std::setw(8) << std::right << body->_design_params_2._ndof << "|" << std::endl;
        }
        std::cerr << "|---------------------------------------|" << std::endl;
    }
    if (_ndof_p3 > 0) {
        std::cerr << "|           Type III           |  #Dof  |" << std::endl;
        std::cerr << "|------------------------------|--------|" << std::endl;
        for (auto body : _bodies) {
            if (body->_design_params_3._active)
                std::cerr << "|" << std::setw(30) << std::left << body->_name << "|" << std::setw(8) << std::right << body->_design_params_3._ndof << "|" << std::endl;
        }
        std::cerr << "|---------------------------------------|" << std::endl;
    }
    if (_ndof_p4 > 0) {
        std::cerr << "|            Type IV           |  #Dof  |" << std::endl;
        std::cerr << "|------------------------------|--------|" << std::endl;
        for (auto body : _bodies) {
            if (body->_design_params_4._active)
                std::cerr << "|" << std::setw(30) << std::left << body->_name << "|" << std::setw(8) << std::right << body->_design_params_4._ndof << "|" << std::endl;
        }
        std::cerr << "|---------------------------------------|" << std::endl;
    }
    if (_ndof_p5 > 0) {
        std::cerr << "|            Type V            |  #Dof  |" << std::endl;
        std::cerr << "|------------------------------|--------|" << std::endl;
        for (auto joint : _joints) {
            if (joint->_design_params_5._active)
                std::cerr << "|" << std::setw(30) << std::left << joint->_name << "|" << std::setw(8) << std::right << joint->_design_params_5._ndof << "|" << std::endl;
        }
        std::cerr << " --------------------------------------- " << std::endl;
    }
    if (_ndof_p6 > 0) {
        std::cerr << "|            Type VI           |  #Dof  |" << std::endl;
        std::cerr << "|------------------------------|--------|" << std::endl;
        for (auto body : _bodies) {
            if (body->_design_params_6._active)
                std::cerr << "|" << std::setw(30) << std::left << body->_name << "|" << std::setw(8) << std::right << body->_design_params_6._ndof << "|" << std::endl;
        }
        std::cerr << "|---------------------------------------|" << std::endl;
    }
}

void Robot::set_contact_scale(dtype scale) {
    for (auto force : _forces) {
        if (dynamic_cast<ForceGeneralPrimitiveContact*>(force) != nullptr) {
            dynamic_cast<ForceGeneralPrimitiveContact*>(force)->set_scale(scale);
        } else if (dynamic_cast<ForceGeneralSDFContact*>(force) != nullptr) {
            dynamic_cast<ForceGeneralSDFContact*>(force)->set_scale(scale);
        } else if (dynamic_cast<ForceGeneralBVHContact*>(force) != nullptr) {
            dynamic_cast<ForceGeneralBVHContact*>(force)->set_scale(scale);
        }
    }
}

// rendering mesh for abstract bodies
void Robot::set_rendering_mesh_vertices(const std::vector<Matrix3X> &Vs) {
    int idx = 0;
    for (auto body : _bodies) {
        if (dynamic_cast<BodyAbstract*>(body) != nullptr) {
            dynamic_cast<BodyAbstract*>(body)->set_rendering_mesh_vertices(Vs[idx]);
            idx += 1;
        }
    }
}

void Robot::set_rendering_mesh(const std::vector<Matrix3X> &Vs, const std::vector<Matrix3Xi> &Fs) {
    int idx = 0;
    for (auto body : _bodies) {
        if (dynamic_cast<BodyAbstract*>(body) != nullptr) {
            dynamic_cast<BodyAbstract*>(body)->set_rendering_mesh(Vs[idx], Fs[idx]);
            idx += 1;
        }
    }
}

void Robot::update_virtual_object(std::string name, VectorX data) {
    for (auto virtual_object : _virtual_objects) {
        if (virtual_object->_name == name)
            virtual_object->update_data(data);
    }
}

void Robot::computeMaximalMassMatrix(VectorX& Mm) {
    // compute the mass matrix in maximal coordinates
    // M = diagonal(M_0, M_1, ..., M_n)
    Mm = VectorX::Zero(_ndof_m);
    for (auto body : _bodies) {
        for (int i = 0;i < 6;i++) {
            Mm(body->_index[i]) = body->_Inertia[i];
        }
    }
}

void Robot::computeForce(VectorX& fm, VectorX& fr) {
    // compute the force vector in maximal coordinates
    // fm = f_coriolis + f_gravity + f_external
    fm = VectorX::Zero(_ndof_m);
    for (auto body : _bodies) {
        body->computeMaximalForce(fm);
    }
    // compute the force vector in reduced coordinates
    fr = VectorX::Zero(_ndof_r);
    for (auto joint : _joints) 
        joint->computeJointForce(fr);
    for (auto force : _forces)
        force->computeForce(fm, fr);

    for (auto actuator : _actuators)
        actuator->computeForce(fm, fr);
}

void Robot::computeForceWithDerivative(VectorX& fm, VectorX& fr, MatrixX& Km, MatrixX& Dm, MatrixX& Kr, MatrixX& Dr, bool verbose) {
    // compute the force vector in maximal coordinates and their derivatives w.r.t phi
    // fm = f_coriolis + f_gravity + f_external
    fm = VectorX::Zero(_ndof_m);
    Km = MatrixX::Zero(_ndof_m, _ndof_m);
    Dm = MatrixX::Zero(_ndof_m, _ndof_m);
    for (auto body : _bodies) {
        body->computeMaximalForceWithDerivative(fm, Km, Dm);
    }
    fr = VectorX::Zero(_ndof_r);
    Kr = MatrixX::Zero(_ndof_r, _ndof_r);
    Dr = MatrixX::Zero(_ndof_r, _ndof_r);
    for (auto joint : _joints) 
        joint->computeJointForceWithDerivative(fr, Kr, Dr);

    for (auto force : _forces) 
        force->computeForceWithDerivative(fm, fr, Km, Dm, Kr, Dr, verbose);
    
    for (auto actuator : _actuators)
        actuator->computeForceWithDerivative(fm, fr, Km, Dm, Kr, Dr);
}

void Robot::computeForceWithDerivative(
    VectorX& fm, VectorX& fr, 
    MatrixX& Km, MatrixX& Dm, MatrixX& Kr, MatrixX& Dr, 
    MatrixX& dfm_dp, MatrixX& dfr_dp, bool verbose) {
    
    // compute the force vector in maximal coordinates and their derivatives w.r.t phi
    // fm = f_coriolis + f_gravity + f_external
    fm = VectorX::Zero(_ndof_m);
    Km = MatrixX::Zero(_ndof_m, _ndof_m);
    Dm = MatrixX::Zero(_ndof_m, _ndof_m);
    dfm_dp = MatrixX::Zero(_ndof_m, _ndof_p);
    for (auto body : _bodies) {
        body->computeMaximalForceWithDerivative(fm, Km, Dm, dfm_dp);
    }
    fr = VectorX::Zero(_ndof_r);
    Kr = MatrixX::Zero(_ndof_r, _ndof_r);
    Dr = MatrixX::Zero(_ndof_r, _ndof_r);
    dfr_dp = MatrixX::Zero(_ndof_r, _ndof_p);
    for (auto joint : _joints) {
        joint->computeJointForceWithDerivative(fr, Kr, Dr, dfr_dp);
    }

    for (auto force : _forces) {
        force->computeForceWithDerivative(fm, fr, Km, Dm, Kr, Dr, dfm_dp, dfr_dp, verbose);
    }
    
    for (auto actuator : _actuators) {
        actuator->computeForceWithDerivative(fm, fr, Km, Dm, Kr, Dr);
    }
}

void Robot::computeForceWithDerivative(
    VectorX& fm, VectorX& fr, 
    MatrixX& Km, MatrixX& Dm, MatrixX& Kr, MatrixX& Dr, 
    SparseMatrixX& dfm_dp_sparse, SparseMatrixX& dfr_dp_sparse, bool verbose) {
    
    // compute the force vector in maximal coordinates and their derivatives w.r.t phi
    // fm = f_coriolis + f_gravity + f_external
    fm = VectorX::Zero(_ndof_m);
    Km = MatrixX::Zero(_ndof_m, _ndof_m);
    Dm = MatrixX::Zero(_ndof_m, _ndof_m);
    dfm_dp.setZero();
    for (auto body : _bodies) {
        body->computeMaximalForceWithDerivative(fm, Km, Dm, dfm_dp);
    }
    fr = VectorX::Zero(_ndof_r);
    Kr = MatrixX::Zero(_ndof_r, _ndof_r);
    Dr = MatrixX::Zero(_ndof_r, _ndof_r);
    dfr_dp.setZero();
    for (auto joint : _joints) {
        joint->computeJointForceWithDerivative(fr, Kr, Dr, dfr_dp);
    }

    for (auto force : _forces) {
        force->computeForceWithDerivative(fm, fr, Km, Dm, Kr, Dr, dfm_dp, dfr_dp, verbose);
    }
    
    for (auto actuator : _actuators) {
        actuator->computeForceWithDerivative(fm, fr, Km, Dm, Kr, Dr);
    }

    // compose dfm_dp
    dfm_dp_sparse = dfm_dp.sparseView();
    
    // compose dfr_dp
    dfr_dp_sparse = dfr_dp.sparseView();
}


void Robot::computeForceWithDerivative(
    VectorX& fm, VectorX& fr, 
    MatrixX& Km, MatrixX& Dm, MatrixX& Kr, MatrixX& Dr, 
    MatrixX& dfm_dp, SparseMatrixX& dfr_dp_sparse, bool verbose) {
    
    // compute the force vector in maximal coordinates and their derivatives w.r.t phi
    // fm = f_coriolis + f_gravity + f_external
    fm = VectorX::Zero(_ndof_m);
    Km = MatrixX::Zero(_ndof_m, _ndof_m);
    Dm = MatrixX::Zero(_ndof_m, _ndof_m);
    dfm_dp = MatrixX::Zero(_ndof_m, _ndof_p);
    for (auto body : _bodies) {
        body->computeMaximalForceWithDerivative(fm, Km, Dm, dfm_dp);
    }
    fr = VectorX::Zero(_ndof_r);
    Kr = MatrixX::Zero(_ndof_r, _ndof_r);
    Dr = MatrixX::Zero(_ndof_r, _ndof_r);
    dfr_dp.setZero();
    for (auto joint : _joints) {
        joint->computeJointForceWithDerivative(fr, Kr, Dr, dfr_dp);
    }

    for (auto force : _forces) {
        force->computeForceWithDerivative(fm, fr, Km, Dm, Kr, Dr, dfm_dp, dfr_dp, verbose);
    }
    
    for (auto actuator : _actuators) {
        actuator->computeForceWithDerivative(fm, fr, Km, Dm, Kr, Dr);
    }

    // compose dfr_dp
    dfr_dp_sparse = dfr_dp.sparseView(0, 0);
}

void Robot::computeJointJacobian(MatrixX& J, MatrixX& Jdot) {
    // compute joint jacobian J_mr: phi = J_mr * qdot
    // compute its time derivative Jdot
    // follow pseudo-code Algorithm (2)
    J = MatrixX::Zero(_ndof_m, _ndof_r);
    Jdot = MatrixX::Zero(_ndof_m, _ndof_r);
    for (auto joint : _joints) {
        if (joint->_ndof > 0) {
            J.block(joint->_body->_index[0], joint->_index[0], 6, joint->_ndof).noalias() = joint->_body->_A_ij * joint->_S_j;
            Jdot.block(joint->_body->_index[0], joint->_index[0], 6, joint->_ndof).noalias() = joint->_body->_A_ij * joint->_S_j_dot;
        }
        Joint* parent = joint->_parent;
        // MatrixX tmp = joint->_A_jp;
        // MatrixX tmp1 = -joint->_A_inv * joint->_Adot * joint->_A_inv * math::Ad(joint->_E_jp_0);
        for (Joint* now = parent; now != nullptr; now = now->_parent) {
            if (now->_ndof > 0) {
                J.block(joint->_body->_index[0], now->_index[0], 6, now->_ndof).noalias() = joint->_body->_A_ip * J.block(parent->_body->_index[0], now->_index[0], 6, now->_ndof);
                Jdot.block(joint->_body->_index[0], now->_index[0], 6, now->_ndof).noalias() =
                    joint->_body->_A_ip_dot * J.block(parent->_body->_index[0], now->_index[0], 6, now->_ndof)
                    + joint->_body->_A_ip * Jdot.block(parent->_body->_index[0], now->_index[0], 6, now->_ndof);
                // J.block(joint->_body->_index[0], now->_index[0], 6, now->_ndof) = joint->_body->_A_ij * tmp * now->_S_j;
                // Jdot.block(joint->_body->_index[0], now->_index[0], 6, now->_ndof).noalias() =
                //     joint->_body->_A_ij * tmp1 * now->_S_j
                //     + joint->_body->_A_ij * tmp * now->_S_j_dot;
            }
            // tmp1 = tmp1 * now->_A_jp - tmp * now->_A_inv * now->_Adot * now->_A_inv * math::Ad(now->_E_jp_0);
            // tmp = tmp * now->_A_jp;
        }
    }
}

void Robot::computeJointJacobianWithDerivative(MatrixX& J, MatrixX& Jdot, JacobianMatrixVector& dJ_dq, JacobianMatrixVector& dJdot_dq) {
    // compute joint jacobian J_mr: phi = J_mr * qdot
    // compute its time derivative Jdot
    // compute their derivative w.r.t q
    // follow pseudo-code Slgorithm (11)
    J = MatrixX::Zero(_ndof_m, _ndof_r);
    Jdot = MatrixX::Zero(_ndof_m, _ndof_r);
    dJ_dq = JacobianMatrixVector(_ndof_m, _ndof_r, _ndof_r);
    dJdot_dq = JacobianMatrixVector(_ndof_m, _ndof_r, _ndof_r);
    for (auto joint : _joints) {
        if (joint->_ndof > 0) {
            J.block(joint->_body->_index[0], joint->_index[0], 6, joint->_ndof).noalias() = joint->_body->_A_ij * joint->_S_j;
            Jdot.block(joint->_body->_index[0], joint->_index[0], 6, joint->_ndof).noalias() = joint->_body->_A_ij * joint->_S_j_dot;
            for (int k = 0;k < joint->_ndof;k++) {
                dJ_dq(joint->_index[k]).block(joint->_body->_index[0], joint->_index[0], 6, joint->_ndof).noalias() = joint->_body->_A_ij * joint->_dSj_dq(k);
                dJdot_dq(joint->_index[k]).block(joint->_body->_index[0], joint->_index[0], 6, joint->_ndof).noalias() = joint->_body->_A_ij * joint->_dSjdot_dq(k);
            }
        }
        Joint* parent = joint->_parent;
        // MatrixX tmp = joint->_A_jp;
        // MatrixX tmp1 = -joint->_A_inv * joint->_Adot * joint->_A_inv * math::Ad(joint->_E_jp_0);
        for (Joint* now = parent; now != nullptr; now = now->_parent) {
            if (now->_ndof > 0) {
                J.block(joint->_body->_index[0], now->_index[0], 6, now->_ndof).noalias() = joint->_body->_A_ip * J.block(parent->_body->_index[0], now->_index[0], 6, now->_ndof);
                Jdot.block(joint->_body->_index[0], now->_index[0], 6, now->_ndof).noalias() =
                    joint->_body->_A_ip_dot * J.block(parent->_body->_index[0], now->_index[0], 6, now->_ndof) +
                    joint->_body->_A_ip * Jdot.block(parent->_body->_index[0], now->_index[0], 6, now->_ndof);
                // J.block(joint->_body->_index[0], now->_index[0], 6, now->_ndof) = joint->_body->_A_ij * tmp * now->_S_j;
                // Jdot.block(joint->_body->_index[0], now->_index[0], 6, now->_ndof).noalias() =
                //     joint->_body->_A_ij * tmp1 * now->_S_j
                //     + joint->_body->_A_ij * tmp * now->_S_j_dot;
                
                for (int k = 0;k < joint->_ndof;k++) {
                    dJ_dq(joint->_index[k]).block(joint->_body->_index[0], now->_index[0], 6, now->_ndof).noalias() = 
                        joint->_body->_dAip_dq(k) * J.block(parent->_body->_index[0], now->_index[0], 6, now->_ndof);

                    dJdot_dq(joint->_index[k]).block(joint->_body->_index[0], now->_index[0], 6, now->_ndof).noalias() = 
                        joint->_body->_dAipdot_dq(k) * J.block(parent->_body->_index[0], now->_index[0], 6, now->_ndof) + 
                        joint->_body->_dAip_dq(k) * Jdot.block(parent->_body->_index[0], now->_index[0], 6, now->_ndof);
                }

                for (Joint* nex = parent; nex != nullptr; nex = nex->_parent) {
                    for (int k = 0;k < nex->_ndof;k++) {
                        dJ_dq(nex->_index[k]).block(joint->_body->_index[0], now->_index[0], 6, now->_ndof).noalias() = 
                            joint->_body->_A_ip * dJ_dq(nex->_index[k]).block(parent->_body->_index[0], now->_index[0], 6, now->_ndof);
                        dJdot_dq(nex->_index[k]).block(joint->_body->_index[0], now->_index[0], 6, now->_ndof).noalias() = 
                            joint->_body->_A_ip_dot * dJ_dq(nex->_index[k]).block(parent->_body->_index[0], now->_index[0], 6, now->_ndof) + 
                            joint->_body->_A_ip * dJdot_dq(nex->_index[k]).block(parent->_body->_index[0], now->_index[0], 6, now->_ndof);
                    }
                    if (nex == now) // TODO: check correctness
                        break;
                }
            }
            // tmp1 = tmp1 * now->_A_jp - tmp * now->_A_inv * now->_Adot * now->_A_inv * math::Ad(now->_E_jp_0);
            // tmp = tmp * now->_A_jp;
        }
    }
}

void Robot::computeJointJacobianWithDerivative(
    MatrixX& J, MatrixX& Jdot, 
    JacobianMatrixVector& dJ_dq, JacobianMatrixVector& dJdot_dq,
    JacobianMatrixVector& dJ_dp1, JacobianMatrixVector& dJ_dp2,
    JacobianMatrixVector& dJdot_dp1, JacobianMatrixVector& dJdot_dp2) {

    computeJointJacobianWithDerivative(J, Jdot, dJ_dq, dJdot_dq);

    dJ_dp1 = JacobianMatrixVector(_ndof_m, _ndof_r, _ndof_p1);
    dJ_dp2 = JacobianMatrixVector(_ndof_m, _ndof_r, _ndof_p2);
    dJdot_dp1 = JacobianMatrixVector(_ndof_m, _ndof_r, _ndof_p1);
    dJdot_dp2 = JacobianMatrixVector(_ndof_m, _ndof_r, _ndof_p2);
    
    for (auto joint : _joints) {
        if (joint->_body == nullptr)
            continue;

        if (joint->_ndof > 0) {
            // design params 2
            if (joint->_body->_design_params_2._active) {
                for (int k = 0;k < joint->_body->_design_params_2._ndof;k++) {
                    int idx = joint->_body->_design_params_2._param_index(k);
                    dJ_dp2(idx).block(joint->_body->_index[0], joint->_index[0], 6, joint->_ndof).noalias() = joint->_body->_dAij_dp2(k) * joint->_S_j;
                    dJdot_dp2(idx).block(joint->_body->_index[0], joint->_index[0], 6, joint->_ndof).noalias() = joint->_body->_dAij_dp2(k) * joint->_S_j_dot;
                }
            }
        }
        Joint* Jp = joint->_parent;
        
        if (Jp == nullptr) continue;

        Body* Bp = joint->_body->_parent;
        Matrix6 A_BiJi_bar = math::Ad(joint->_body->_E_ij * joint->_Q_inv);
        Matrix6 A_JiBp = joint->_A_jp * Bp->_A_ji;
        Matrix6 Aleft = A_BiJi_bar * joint->_Adot * joint->_A_inv;
        Matrix6 Aright = joint->_A_inv * joint->_Adot * A_JiBp;
        for (Joint* now = Jp; now != nullptr; now = now->_parent) {
            if (now->_ndof > 0) {
                // design params 1
                for (Joint* middle = joint; middle != now; middle = middle->_parent) {
                    if (middle->_design_params_1._active) {
                        // dJ_dp1(Bi, j) = dAip_dp1 * J(Bp, j) + Aip * dJ_dp1(Bp, j)
                        // dJdot_dp1 = dAdotip_dp1 * J(Bp, j) + Adot_ip * dJ_dp1(Bp, j) + dAip_dp1 * Jdot(Bp, j) + Aip * dJdot_dp1(Bp, j)

                        // Aip * dJ_dp1(Bp, j)
                        // Adot_ip * dJ_dp1(Bp, j) + Aip * dJdot_dp1(Bp, j)
                        if (middle != joint) {
                            for (int k = 0;k < middle->_design_params_1._ndof;k++) {
                                int idx = middle->_design_params_1._param_index(k);
                                dJ_dp1(idx).block(joint->_body->_index[0], now->_index[0], 6, now->_ndof).noalias() += joint->_body->_A_ip * dJ_dp1(idx).block(Bp->_index[0], now->_index[0], 6, now->_ndof);
                                dJdot_dp1(idx).block(joint->_body->_index[0], now->_index[0], 6, now->_ndof).noalias() += 
                                    joint->_body->_A_ip_dot * dJ_dp1(idx).block(Bp->_index[0], now->_index[0], 6, now->_ndof)
                                    + joint->_body->_A_ip * dJdot_dp1(idx).block(Bp->_index[0], now->_index[0], 6, now->_ndof);
                            }
                        }

                        // dAip_dp1 * J(Bp, j)
                        // dAdotip_dp1 * J(Bp, j) + dAip_dp1 * Jdot(Bp, j)
                        if (middle == joint) {
                            for (int k = 0;k < middle->_design_params_1._ndof;k++) {
                                int idx = middle->_design_params_1._param_index(k);
                                Matrix6 dAip_dp1 = A_BiJi_bar * joint->_dAjp0_dp1(k) * Bp->_A_ji;
                                Matrix6 dAipdot_dp1 = -Aleft * joint->_dAjp0_dp1(k) * Bp->_A_ji;
                                
                                dJ_dp1(idx).block(joint->_body->_index[0], now->_index[0], 6, now->_ndof).noalias() += dAip_dp1 * J.block(Bp->_index[0], now->_index[0], 6, now->_ndof);
                                dJdot_dp1(idx).block(joint->_body->_index[0], now->_index[0], 6, now->_ndof).noalias() += 
                                    dAipdot_dp1 * J.block(Bp->_index[0], now->_index[0], 6, now->_ndof)
                                     + dAip_dp1 * Jdot.block(Bp->_index[0], now->_index[0], 6, now->_ndof);
                            }
                        }
                    }
                }

                // design params 2
                if (joint->_body->_design_params_2._active) {
                    // For type II of Bi: 
                    // dJ_dp2(Bi, j) = dAip_dp2 * J(Bp, j)
                    // dJdot_dp2(Bi, j) = dAipdot_dp2 * J(BP, j) + dAip_dp2 * Jdot(Bp, j)
                    for (int k = 0;k < joint->_body->_design_params_2._ndof;k++) {
                        int idx = joint->_body->_design_params_2._param_index(k);
                        Matrix6 dAip_dp2 = joint->_body->_dAij_dp2(k) * A_JiBp;
                        Matrix6 dAipdot_dp2 = -joint->_body->_dAij_dp2(k) * Aright;

                        dJ_dp2(idx).block(joint->_body->_index[0], now->_index[0], 6, now->_ndof).noalias() += dAip_dp2 * J.block(Bp->_index[0], now->_index[0], 6, now->_ndof);
                        dJdot_dp2(idx).block(joint->_body->_index[0], now->_index[0], 6, now->_ndof).noalias() += 
                            dAipdot_dp2 * J.block(Bp->_index[0], now->_index[0], 6, now->_ndof)
                            + dAip_dp2 * Jdot.block(Bp->_index[0], now->_index[0], 6, now->_ndof);
                    }
                }
            }
        }
    }
}

void Robot::computeJointJacobianWithDerivative(
    MatrixX& J, MatrixX& Jdot, 
    JacobianMatrixVector& dJ_dq, JacobianMatrixVector& dJdot_dq,
    SparseJacobianMatrixVector& dJ_dp1_sparse, SparseJacobianMatrixVector& dJ_dp2_sparse,
    SparseJacobianMatrixVector& dJdot_dp1_sparse, SparseJacobianMatrixVector& dJdot_dp2_sparse) {

    computeJointJacobianWithDerivative(J, Jdot, dJ_dq, dJdot_dq);

    for (auto joint : _joints) {   
        if (joint->_body == nullptr)
            continue;

        // special processing for design params 2
        if (joint->_ndof > 0) {
            // design params 2
            if (joint->_body->_design_params_2._active) {
                for (int k = 0;k < 12;k++) {
                    dJ_dp2(k).block(joint->_body->_index[0], joint->_index[0], 6, joint->_ndof).noalias() = joint->_body->_dAij_dp2(k) * joint->_S_j;
                    dJdot_dp2(k).block(joint->_body->_index[0], joint->_index[0], 6, joint->_ndof).noalias() = joint->_body->_dAij_dp2(k) * joint->_S_j_dot;
                }
            }
        }

        Joint* Jp = joint->_parent;
        
        if (Jp == nullptr) continue;

        Body* Bp = joint->_body->_parent;
        Matrix6 A_BiJi_bar = math::Ad(joint->_body->_E_ij * joint->_Q_inv);
        Matrix6 A_JiBp = joint->_A_jp * Bp->_A_ji;
        Matrix6 Aleft = A_BiJi_bar * joint->_Adot * joint->_A_inv;
        Matrix6 Aright = joint->_A_inv * joint->_Adot * A_JiBp;

        for (Joint* now = Jp; now != nullptr; now = now->_parent) {
            if (now->_ndof > 0) {
                // design params 1
                for (Joint* middle = joint; middle != now; middle = middle->_parent) {
                    if (middle->_design_params_1._active) {
                        // dJ_dp1(Bi, j) = dAip_dp1 * J(Bp, j) + Aip * dJ_dp1(Bp, j)
                        // dJdot_dp1 = dAdotip_dp1 * J(Bp, j) + Adot_ip * dJ_dp1(Bp, j) + dAip_dp1 * Jdot(Bp, j) + Aip * dJdot_dp1(Bp, j)

                        // Aip * dJ_dp1(Bp, j)
                        // Adot_ip * dJ_dp1(Bp, j) + Aip * dJdot_dp1(Bp, j)
                        if (middle != joint) {
                            for (int k = 0;k < 12;k++) {
                                int idx = middle->_design_params_1._param_index(k);
                                dJ_dp1(idx).block(joint->_body->_index[0], now->_index[0], 6, now->_ndof).noalias() = joint->_body->_A_ip * dJ_dp1(idx).block(Bp->_index[0], now->_index[0], 6, now->_ndof);
                                dJdot_dp1(idx).block(joint->_body->_index[0], now->_index[0], 6, now->_ndof).noalias() = 
                                    joint->_body->_A_ip_dot * dJ_dp1(idx).block(Bp->_index[0], now->_index[0], 6, now->_ndof)
                                    + joint->_body->_A_ip * dJdot_dp1(idx).block(Bp->_index[0], now->_index[0], 6, now->_ndof);
                            }
                        }
                        if (middle == joint) {
                            auto J_right = Bp->_A_ji * J.block(Bp->_index[0], now->_index[0], 6, now->_ndof);
                            auto Jdot_right = Bp->_A_ji * Jdot.block(Bp->_index[0], now->_index[0], 6, now->_ndof);
                            for (int k = 0;k < 12;k++) {
                                int idx = middle->_design_params_1._param_index(k);
                                Matrix6 dAip_dp1 = A_BiJi_bar * joint->_dAjp0_dp1(k);
                                
                                dJ_dp1(idx).block(joint->_body->_index[0], now->_index[0], 6, now->_ndof).noalias() = dAip_dp1 * J_right;
                                dJdot_dp1(idx).block(joint->_body->_index[0], now->_index[0], 6, now->_ndof).noalias() = 
                                    -Aleft * joint->_dAjp0_dp1(k) * J_right
                                     + dAip_dp1 * Jdot_right;
                            }
                        }
                    }
                }

                // design params 2
                if (joint->_body->_design_params_2._active) {
                    // For type II of Bi: 
                    // dJ_dp2(Bi, j) = dAip_dp2 * J(Bp, j)
                    // dJdot_dp2(Bi, j) = dAipdot_dp2 * J(BP, j) + dAip_dp2 * Jdot(Bp, j)
                    for (int k = 0;k < 12;k++) {
                        Matrix6 dAip_dp2 = joint->_body->_dAij_dp2(k) * A_JiBp;

                        dJ_dp2(k).block(joint->_body->_index[0], now->_index[0], 6, now->_ndof).noalias() = dAip_dp2 * J.block(Bp->_index[0], now->_index[0], 6, now->_ndof);
                        dJdot_dp2(k).block(joint->_body->_index[0], now->_index[0], 6, now->_ndof).noalias() = 
                            -joint->_body->_dAij_dp2(k) * Aright * J.block(Bp->_index[0], now->_index[0], 6, now->_ndof)
                            + dAip_dp2 * Jdot.block(Bp->_index[0], now->_index[0], 6, now->_ndof);
                    }
                }
            }
        }
    }

    // compose sparse matrices for dJ_dp1, dJ_dp2, dJdot_dp1, dJdot_dp2
    dJ_dp1_sparse = dJ_dp1.sparseView();
    dJdot_dp1_sparse = dJdot_dp1.sparseView();

    dJ_dp2_sparse = SparseJacobianMatrixVector(_ndof_m, _ndof_r, _ndof_p2);
    dJdot_dp2_sparse = SparseJacobianMatrixVector(_ndof_m, _ndof_r, _ndof_p2);
    for (int k = 0;k < 12;k++) {
        for (auto body : _bodies) {
            if (body->_design_params_2._active) {
                std::vector<TripletX> triplets;
                for (int i = 0;i < 6;i++)
                    for (int j = 0;j < _ndof_r;j++) {
                        if (fabs(dJ_dp2(k)(body->_index[i], j)) > 0) {
                            triplets.push_back(TripletX(body->_index[i], j, dJ_dp2(k)(body->_index[i], j)));
                        }
                    }
                dJ_dp2_sparse(body->_design_params_2._param_index(k)).setFromTriplets(triplets.begin(), triplets.end());
            
                triplets.clear();
                for (int i = 0;i < 6;i++)
                    for (int j = 0;j < _ndof_r;j++) {
                        if (fabs(dJdot_dp2(k)(body->_index[i], j)) > 0) {
                            triplets.push_back(TripletX(body->_index[i], j, dJdot_dp2(k)(body->_index[i], j)));
                        }
                    }
                dJdot_dp2_sparse(body->_design_params_2._param_index(k)).setFromTriplets(triplets.begin(), triplets.end());
            }
        }
    }
}

void Robot::compute_dfdu(MatrixX& dfm_du, MatrixX& dfr_du) {
    dfm_du = MatrixX::Zero(_ndof_m, _ndof_u);
    dfr_du = MatrixX::Zero(_ndof_r, _ndof_u);
    for (auto actuator : _actuators) {
        actuator->compute_dfdu(dfm_du, dfr_du);
    }
}

void Robot::computeVariables(VectorX& variables) {
    variables = VectorX::Zero(_ndof_var);
    for (auto end_effector : _end_effectors) {
        end_effector->computeVariables(variables);
    }
}

void Robot::computeVariablesWithDerivative(VectorX& variables, MatrixX& dvar_dq) {
    variables = VectorX::Zero(_ndof_var);
    dvar_dq = MatrixX::Zero(_ndof_var, _ndof_r);
    for (auto end_effector : _end_effectors) {
        end_effector->computeVariablesWithDerivative(variables, dvar_dq);
    }
}

void Robot::computeVariablesWithDerivative(VectorX& variables, MatrixX& dvar_dq, MatrixX& dvar_dp) {
    variables = VectorX::Zero(_ndof_var);
    dvar_dq = MatrixX::Zero(_ndof_var, _ndof_r);
    dvar_dp = MatrixX::Zero(_ndof_var, _ndof_p);
    for (auto end_effector : _end_effectors) {
        end_effector->computeVariablesWithDerivative(variables, dvar_dq, dvar_dp);
    }
}

void Robot::test_derivatives_runtime() {
    for (Joint* joint : _joints) {
        joint->test_derivatives_runtime();
    }

    for (Body* body : _bodies) {
        body->test_derivatives_runtime();
    }

    for (Force* force : _forces) {
        force->test_derivatives_runtime();
    }

    VectorX fm, fr;
    MatrixX Km, Dm, Kr, Dr;
    computeForceWithDerivative(fm, fr, Km, Dm, Kr, Dr, true);

    MatrixX J, Jdot;
    JacobianMatrixVector dJ_dq, dJdot_dq;
    computeJointJacobianWithDerivative(J, Jdot, dJ_dq, dJdot_dq);

    MatrixX dfm_du, dfr_du;
    compute_dfdu(dfm_du, dfr_du);

    VectorX variables;
    MatrixX dvar_dq;
    MatrixX dvar_dp;
    computeVariablesWithDerivative(variables, dvar_dq, dvar_dp);

    auto q = get_q();
    auto qdot = get_qdot();
    auto u = get_u();

    dtype h = 1e-7;

    for (int ii = 0;ii < 1;ii++) {
        // printf("---------------------- eps = %.9lf ----------------------------\n", h);
        VectorX phi = get_phi();
        print_error("Robot: J", phi, J * qdot);

        MatrixX Kr_fd = MatrixX::Zero(_ndof_r, _ndof_r);
        JacobianMatrixVector dJ_dq_fd(_ndof_m, _ndof_r, _ndof_r);
        JacobianMatrixVector dJdot_dq_fd(_ndof_m, _ndof_r, _ndof_r);
        MatrixX dvar_dq_fd = MatrixX::Zero(_ndof_var, _ndof_r);
        
        for (int k = 0;k < _ndof_r;k++) {
            auto q_pos = q;
            q_pos(k) += h;
            set_q(q_pos);
            update();

            VectorX fm_pos, fr_pos;
            computeForce(fm_pos, fr_pos);

            MatrixX J_pos, Jdot_pos;
            computeJointJacobian(J_pos, Jdot_pos);

            VectorX variables_pos;
            computeVariables(variables_pos);

            Kr_fd.col(k) = (fr_pos - fr) / h;
            dJ_dq_fd(k) = (J_pos - J) / h;
            dJdot_dq_fd(k) = (Jdot_pos - Jdot) / h;
            dvar_dq_fd.col(k) = (variables_pos - variables) / h;
        }

        print_error("Robot: Kr", Kr, Kr_fd);
        print_error("Robot: dJ_dq", dJ_dq, dJ_dq_fd);
        print_error("Robot: dJdot_dq", dJdot_dq, dJdot_dq_fd);
        print_error("Robot: dvar_dq", dvar_dq_fd, dvar_dq);

        set_q(q);
        update();
        
        MatrixX Dr_fd = MatrixX::Zero(_ndof_r, _ndof_r);

        for (int k = 0;k < _ndof_r;k++) {
            auto qdot_pos = qdot;
            qdot_pos(k) += h;
            set_qdot(qdot_pos);
            update();

            VectorX fm_pos, fr_pos;
            computeForce(fm_pos, fr_pos);

            Dr_fd.col(k) = (fr_pos - fr) / h;
        }

        print_error("Robot: Dr", Dr, Dr_fd);

        set_qdot(qdot);
        update();

        MatrixX Km_fd = MatrixX::Zero(_ndof_m, _ndof_m);
        
        for (int k = 0;k < _ndof_m;k++) {
            int body_idx = k / 6;
            Matrix4 E = _bodies[body_idx]->_E_0i;
            Vector6 dq = Vector6::Zero();
            dq[k % 6] = h;
            Matrix4 E_pos = E * math::exp(dq);
            _bodies[body_idx]->_E_0i = E_pos;
            VectorX fm_pos, fr_pos;
            computeForce(fm_pos, fr_pos);
            Km_fd.col(k) = (fm_pos - fm) / h;
            _bodies[body_idx]->_E_0i = E;
        }

        print_error("Robot: Km:", Km, Km_fd);

        MatrixX Dm_fd = MatrixX::Zero(_ndof_m, _ndof_m);
        
        for (int k = 0;k < _ndof_m;k++) {
            int body_idx = k / 6;
            Vector6 phi = _bodies[body_idx]->_phi;
            Vector6 phi_pos = phi;
            phi_pos[k % 6] += h;
            _bodies[body_idx]->_phi = phi_pos;
            VectorX fm_pos, fr_pos;
            computeForce(fm_pos, fr_pos);
            Dm_fd.col(k) = (fm_pos - fm) / h;
            _bodies[body_idx]->_phi = phi;
        }

        print_error("Robot: Dm:", Dm, Dm_fd);

        MatrixX dfm_du_fd = MatrixX::Zero(_ndof_m, _ndof_u);
        MatrixX dfr_du_fd = MatrixX::Zero(_ndof_r, _ndof_u);
        for (int k = 0;k < _ndof_u;k++) {
            VectorX u_pos = u;
            u_pos[k] += h;
            set_u(u_pos);
            VectorX fm_pos, fr_pos;
            computeForce(fm_pos, fr_pos);
            dfm_du_fd.col(k) = (fm_pos - fm) / h;
            dfr_du_fd.col(k) = (fr_pos - fr) / h;
        }

        set_u(u);

        print_error("Robot: dfm_du", dfm_du, dfm_du_fd);
        print_error("Robot: dfr_du", dfr_du, dfr_du_fd);

        h /= 10.;
    }

    // test design derivatives
    test_design_derivatives_runtime();
}

void Robot::test_design_derivatives_runtime() {
    dtype eps = 1e-7;
    // test joint design derivatives
    // test for type 1
    if (_ndof_p1 > 0) {
        VectorX design_params = get_design_params();
        VectorX design_params_1 = design_params.head(_ndof_p1);
        update(true);
        
        std::vector<JacobianMatrixVector> dAjp0_dp1, dE0j_dp1;
        std::vector<MatrixX> dphi_dp1;
        std::vector<MatrixX> Ajp0, E0j;
        std::vector<VectorX> phi;
        for (auto joint : _joints) {
            dAjp0_dp1.push_back(joint->_dAjp0_dp1);
            dE0j_dp1.push_back(joint->_dE0j_dp1);
            dphi_dp1.push_back(joint->_dphi_dp1);
            Ajp0.push_back(math::Ad(joint->_E_jp_0));
            E0j.push_back(joint->_E_0j);
            phi.push_back(joint->_phi);
        }

        std::vector<JacobianMatrixVector> dAjp0_dp1_fd, dE0j_dp1_fd;
        std::vector<MatrixX> dphi_dp1_fd;
        for (auto joint : _joints) {
            dAjp0_dp1_fd.push_back(JacobianMatrixVector(6, 6, _ndof_p1));
            dE0j_dp1_fd.push_back(JacobianMatrixVector(4, 4, _ndof_p1));
            dphi_dp1_fd.push_back(MatrixX::Zero(6, _ndof_p1));
        }
        for (int i = 0;i < design_params_1.size();i++) {
            VectorX design_params_pos = design_params;
            design_params_pos[i] += eps;
            set_design_params(design_params_pos);
            update(false);
            for (int j = 0;j < _joints.size();j++) {
                auto joint = _joints[j];
                dAjp0_dp1_fd[j](i) = (math::Ad(joint->_E_jp_0) - Ajp0[j]) / eps;
                dE0j_dp1_fd[j](i) = (joint->_E_0j - E0j[j]) / eps;
                dphi_dp1_fd[j].col(i) = (joint->_phi - phi[j]) / eps;
            }
        }
        for (int j = 0;j < _joints.size();j++) {
            print_error("Joint: dAjp0_dp1", dAjp0_dp1[j], dAjp0_dp1_fd[j]);
            print_error("Joint: dE0j_dp1", dE0j_dp1[j], dE0j_dp1_fd[j]);
            print_error("Joint: dphi_dp1", dphi_dp1[j], dphi_dp1_fd[j]);
        }

        set_design_params(design_params);
        update(true);
    }

    // test body design derivaitves
    // test for type 1
    if (_ndof_p1 > 0) {
        VectorX design_params = get_design_params();
        VectorX design_params_1 = design_params.head(_ndof_p1);
        update(true);
        
        std::vector<JacobianMatrixVector> dE0i_dp1;
        std::vector<MatrixX> dphi_dp1;
        std::vector<MatrixX> E0i;
        std::vector<VectorX> phi;
        for (auto body : _bodies) {
            dE0i_dp1.push_back(body->_dE0i_dp1);
            dphi_dp1.push_back(body->_dphi_dp1);
            E0i.push_back(body->_E_0i);
            phi.push_back(body->_phi);
        }

        std::vector<JacobianMatrixVector> dE0i_dp1_fd;
        std::vector<MatrixX> dphi_dp1_fd;
        for (auto body : _bodies) {
            dE0i_dp1_fd.push_back(JacobianMatrixVector(4, 4, _ndof_p1));
            dphi_dp1_fd.push_back(MatrixX::Zero(6, _ndof_p1));
        }
        for (int i = 0;i < design_params_1.size();i++) {
            VectorX design_params_pos = design_params;
            design_params_pos[i] += eps;
            set_design_params(design_params_pos);
            update(false);
            for (int j = 0;j < _bodies.size();j++) {
                auto body = _bodies[j];
                dE0i_dp1_fd[j](i) = (body->_E_0i - E0i[j]) / eps;
                dphi_dp1_fd[j].col(i) = (body->_phi - phi[j]) / eps;
            }
        }
        for (int j = 0;j < _bodies.size();j++) {
            print_error("Body: dE0i_dp1", dE0i_dp1[j], dE0i_dp1_fd[j]);
            print_error("Body: dphi_dp1", dphi_dp1[j], dphi_dp1_fd[j]);
        }

        set_design_params(design_params);
        update(true);
    }
    // test for type 2
    if (_ndof_p2 > 0) {
        VectorX design_params = get_design_params();

        update(true);

        std::vector<JacobianMatrixVector> dAij_dp2, dE0i_dp2;
        std::vector<MatrixX> dphi_dp2;
        std::vector<MatrixX> Aij, E0i;
        std::vector<VectorX> phi;
        for (auto body : _bodies) {
            if (body->_design_params_2._active) {
                dAij_dp2.push_back(body->_dAij_dp2);
                dE0i_dp2.push_back(body->_dE0i_dp2);
                dphi_dp2.push_back(body->_dphi_dp2);
                Aij.push_back(body->_A_ij);
                E0i.push_back(body->_E_0i);
                phi.push_back(body->_phi);
            }
        }

        std::vector<JacobianMatrixVector> dAij_dp2_fd, dE0i_dp2_fd;
        std::vector<MatrixX> dphi_dp2_fd;
        int i = 0;
        for (auto body : _bodies) {
            if (body->_design_params_2._active) {
                dAij_dp2_fd.push_back(JacobianMatrixVector(6, 6, body->_design_params_2._ndof));
                dE0i_dp2_fd.push_back(JacobianMatrixVector(4, 4, body->_design_params_2._ndof));
                dphi_dp2_fd.push_back(MatrixX::Zero(6, body->_design_params_2._ndof));

                for (int j = 0;j < body->_design_params_2._ndof;j++) {
                    int idx = body->_design_params_2._param_index[j];
                    VectorX design_params_pos = design_params;
                    design_params_pos[_ndof_p1 + idx] += eps;
                    set_design_params(design_params_pos);
                    update(false);
                    dAij_dp2_fd[i](j) = (body->_A_ij - Aij[i]) / eps;
                    dE0i_dp2_fd[i](j) = (body->_E_0i - E0i[i]) / eps;
                    dphi_dp2_fd[i].col(j) = (body->_phi - phi[i]) / eps;
                }

                i += 1;
            }
        }

        set_design_params(design_params);
        update(true);

        for (int j = 0;j < dAij_dp2_fd.size();j++) {
            print_error("Body: dAij_dp2", dAij_dp2[j], dAij_dp2_fd[j]);
            print_error("Body: dE0i_dp2", dE0i_dp2[j], dE0i_dp2_fd[j]);
            print_error("Body: dphi_dp2", dphi_dp2[j], dphi_dp2_fd[j]);
        }
    }

    // test joint jacobian design derivatives
    // test for type 1
    if (_ndof_p1 > 0) {
        VectorX design_params = get_design_params();
        VectorX design_params_1 = design_params.head(_ndof_p1);
        update(true);
        
        MatrixX J, Jdot;
        JacobianMatrixVector dJ_dq, dJdot_dq;
        JacobianMatrixVector dJ_dp1, dJdot_dp1, dJ_dp2, dJdot_dp2;
        computeJointJacobianWithDerivative(J, Jdot, dJ_dq, dJdot_dq, dJ_dp1, dJ_dp2, dJdot_dp1, dJdot_dp2);

        JacobianMatrixVector dJ_dp1_fd(_ndof_m, _ndof_r, _ndof_p1);
        JacobianMatrixVector dJdot_dp1_fd(_ndof_m, _ndof_r, _ndof_p1);
        for (int i = 0;i < design_params_1.size();i++) {
            VectorX design_params_pos = design_params;
            design_params_pos[i] += eps;
            set_design_params(design_params_pos);
            update(false);
            MatrixX J_pos, Jdot_pos;
            computeJointJacobian(J_pos, Jdot_pos);
            dJ_dp1_fd(i) = (J_pos - J) / eps;
            dJdot_dp1_fd(i) = (Jdot_pos - Jdot) / eps;
        }
        
        print_error("Jacobian: dJ_dp1", dJ_dp1, dJ_dp1_fd);
        print_error("Jacobian: dJdot_dp1", dJdot_dp1, dJdot_dp1_fd);

        set_design_params(design_params);
        update(true);
    }

    // test for type 2
    if (_ndof_p2 > 0) {
        VectorX design_params = get_design_params();
        VectorX design_params_2 = design_params.segment(_ndof_p1, _ndof_p2);
        update(true);
        
        MatrixX J, Jdot;
        JacobianMatrixVector dJ_dq, dJdot_dq;
        JacobianMatrixVector dJ_dp1, dJdot_dp1, dJ_dp2, dJdot_dp2;
        computeJointJacobianWithDerivative(J, Jdot, dJ_dq, dJdot_dq, dJ_dp1, dJ_dp2, dJdot_dp1, dJdot_dp2);

        JacobianMatrixVector dJ_dp2_fd(_ndof_m, _ndof_r, _ndof_p2);
        JacobianMatrixVector dJdot_dp2_fd(_ndof_m, _ndof_r, _ndof_p2);
        for (int i = 0;i < design_params_2.size();i++) {
            VectorX design_params_pos = design_params;
            design_params_pos[_ndof_p1 + i] += eps;
            set_design_params(design_params_pos);
            update(false);
            MatrixX J_pos, Jdot_pos;
            computeJointJacobian(J_pos, Jdot_pos);
            dJ_dp2_fd(i) = (J_pos - J) / eps;
            dJdot_dp2_fd(i) = (Jdot_pos - Jdot) / eps;
        }
        print_error("Jacobian: dJ_dp2", dJ_dp2, dJ_dp2_fd);
        print_error("Jacobian: dJdot_dp2", dJdot_dp2, dJdot_dp2_fd);

        set_design_params(design_params);
        update(true);
    }

    // test dfm_dp, dfr_dp
    if (_ndof_p > 0) {
        VectorX design_params = get_design_params();
        
        update(true);
        
        VectorX fm, fr;
        MatrixX Km, Dm, Kr, Dr, dfm_dp, dfr_dp;
        computeForceWithDerivative(fm, fr, Km, Dm, Kr, Dr, dfm_dp, dfr_dp);

        MatrixX dfm_dp_fd(_ndof_m, _ndof_p);
        MatrixX dfr_dp_fd(_ndof_r, _ndof_p);
        for (int i = 0;i < _ndof_p;i++) {
            VectorX design_params_pos = design_params;
            design_params_pos[i] += eps;
            set_design_params(design_params_pos);
            update(false);
            VectorX fm_pos, fr_pos;
            computeForce(fm_pos, fr_pos);
            dfm_dp_fd.col(i) = (fm_pos - fm) / eps;
            dfr_dp_fd.col(i) = (fr_pos - fr) / eps;
        }
        
        print_error("Robot: dfm_dp1", dfm_dp.middleCols(0, _ndof_p1), dfm_dp_fd.middleCols(0, _ndof_p1));
        print_error("Robot: dfm_dp2", dfm_dp.middleCols(_ndof_p1, _ndof_p2), dfm_dp_fd.middleCols(_ndof_p1, _ndof_p2));
        print_error("Robot: dfm_dp3", dfm_dp.middleCols(_ndof_p1 + _ndof_p2, _ndof_p3), dfm_dp_fd.middleCols(_ndof_p1 + _ndof_p2, _ndof_p3));
        print_error("Robot: dfm_dp4", dfm_dp.middleCols(_ndof_p1 + _ndof_p2 + _ndof_p3, _ndof_p4), dfm_dp_fd.middleCols(_ndof_p1 + _ndof_p2 + _ndof_p3, _ndof_p4));
        print_error("Robot: dfm_dp6", dfm_dp.middleCols(_ndof_p1 + _ndof_p2 + _ndof_p3 + _ndof_p4 + _ndof_p5, _ndof_p6), dfm_dp_fd.middleCols(_ndof_p1 + _ndof_p2 + _ndof_p3 + _ndof_p4 + _ndof_p5, _ndof_p6));
        print_error("Robot: dfr_dp", dfr_dp, dfr_dp_fd);

        set_design_params(design_params);
        update(true);
    }

    // test dvar_dp
    if (_ndof_p > 0) {
        VectorX design_params = get_design_params();
        
        update(true);

        VectorX variables;
        MatrixX dvar_dq, dvar_dp;
        computeVariablesWithDerivative(variables, dvar_dq, dvar_dp);

        MatrixX dvar_dp_fd(_ndof_var, _ndof_p);
        for (int i = 0;i < _ndof_p;i++) {
            VectorX design_params_pos = design_params;
            design_params_pos[i] += eps;
            set_design_params(design_params_pos);
            update(false);
            VectorX variables_pos;
            computeVariables(variables_pos);
            dvar_dp_fd.col(i) = (variables_pos - variables) / eps;
        }
        
        print_error("Robot: dvar_dp ", dvar_dp, dvar_dp_fd);
        set_design_params(design_params);
        update(true);
    }
}

// rendering related
void Robot::get_rendering_objects(
    std::vector<Matrix3Xf>& vertex_list, 
    std::vector<Matrix3Xi>& face_list,
    std::vector<opengl_viewer::Option>& option_list,
    std::vector<opengl_viewer::Animator*>& animator_list) {
    
    for (Body* body : _bodies) {
        body->get_rendering_objects(vertex_list, face_list, option_list, animator_list);
    }

    for (EndEffector* end_effector : _end_effectors) {
        end_effector->get_rendering_objects(vertex_list, face_list, option_list, animator_list);
    }

    for (VirtualObject* virtual_object : _virtual_objects) {
        virtual_object->get_rendering_objects(vertex_list, face_list, option_list, animator_list);
    }
}

void Robot::reset_time_report() {
    for (auto force : _forces) {
        force->reset_time_report();
    }
}

void Robot::print_time_report() {
    std::cerr << "|---------------------------------|" << std::endl;
    std::cerr << "|Force                            |" << std::endl;
    std::cerr << "|---------------------------------|" << std::endl;
    for (auto force : _forces) {
        force->print_time_report();
    }
}

}