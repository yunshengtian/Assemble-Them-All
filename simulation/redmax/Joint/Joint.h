#pragma once
#include "Common.h"
#include "Utils.h"
#include "Joint/JointDesignParameters.h"

namespace redmax {

class Body;
class Simulation;

class Joint {
public:
    enum Frame {
        LOCAL, // parent joint's frame
        WORLD
    };

    // simulation
    Simulation* _sim;

    // id for joint
    int _id;
    std::string _name = "";

    // structure
    Joint* _parent;
    vector<Joint*> _children;
    Body* _body;

    // indices
    vector<int> _index;         // index for dof in reduced coordinates

    // data
    // constants
    SE3 _E_pj_0, _E_jp_0;       // initial transformation from j to its parent joint
    SE3 _E_j0_0;
    dtype _Kr, _Dr;             // stiffness and damping coefficients of joint
    VectorX _q_rest;            // rest state of joints
    dtype _joint_limit_lower, _joint_limit_upper; // joint limits constraints
    dtype _joint_limit_stiffness; // stiffness coefficients for joint limits constraints

    // variables
    SE3 _Q, _Q_inv;             // transformation from joint j to its parent joint induced by q
    Matrix6 _A, _A_inv, _Adot;  // adjoint matrix of Q
    SE3 _E_pj, _E_jp;           // transformation between joint j and its parent joint
    SE3 _E_0j;                  // transformation between joint j and world frame
    Matrix6 _A_jp;              // adjoint matrix between joint j anod its parent joint
    MatrixX _S_j;               // joint jacobian, phi_j = _S_j * qdot
    MatrixX _S_j_dot;          
    se3 _phi;                   // joint twist

    // design parameters
    JointDesignParameters _design_params_1, _design_params_5;

    // q
    int _ndof;                  // number of dof
    VectorX _q, _qdot;  

    // derivatives
    JacobianMatrixVector _dSj_dq, _dSjdot_dq;     
    JacobianMatrixVector _dA_dq, _dAdot_dq; 
    JacobianMatrixVector _dQ_dq, _dEpj_dq; 

    // design derivatives
    JacobianMatrixVector _dAjp0_dp1;
    JacobianMatrixVector _dE0j_dp1;
    MatrixX _dphi_dp1;

    Joint(Simulation *sim, int id, Joint* parent, Matrix3 R, Vector3 p, int ndof, Frame frame = Frame::LOCAL);

    // init joint
    void init();

    // stiffness and damping setting
    void set_damping(dtype damping);
    void set_stiffness(dtype stiffness, VectorX q_rest);
    
    // joint limit setting
    void set_joint_limit(dtype lower, dtype upper, dtype joint_limit_stiffness);

    // update data
    virtual void update(bool design_gradient = false);

    // activate design parameters
    void activate_design_parameters_type_1(bool active = true);
    void activate_design_parameters_type_5(bool active = true);

    // compute joint force
    void computeJointForce(VectorX& fr);
    void computeJointForceWithDerivative(VectorX& fr, MatrixX& Kr, MatrixX& Dr);
    void computeJointForceWithDerivative(VectorX& fr, MatrixX& Kr, MatrixX& Dr, MatrixX& dfr_dp);

    // compute position/velocity in world frame from local position
    Vector3 position_in_world(Vector3& pos);
    Vector3 velocity_in_world(Vector3& pos);

    // test derivatives
    void test_derivatives();
    void test_derivatives_runtime();

    // check if q is valid
    virtual bool check_valid() { return true; }

    // reparameterize state
    virtual bool reparam() { return false; }

private:
    virtual void update_design_derivatives();
};

}