#pragma once
#include "Actuator/Actuator.h"

namespace redmax {

class Joint;

class ActuatorMotor : public Actuator {
public:
    enum ControlMode {
        FORCE, // input to FORCE mode is [-1, 1] mapping to [ctrl_min, ctrl_max]
        POS // input to POS mode is target_dof, the torque is P * (target_dof - current_dof) + D * (target_dof_vel - current_dof_vel)
    };

    Joint* _joint; // the joint to apply the motor force.
    ControlMode _control_mode; // the control mode
    VectorX _ctrl_P, _ctrl_D;
    // stored temporary variables
    VectorX _pos_error, _vel_error;
    
    ActuatorMotor(std::string name, Joint* joint, ControlMode control_mode = ControlMode::FORCE, 
                    VectorX ctrl_min = VectorX::Constant(1, (dtype)INT_MIN), VectorX ctrl_max = VectorX::Constant(1, (dtype)INT_MAX), 
                    VectorX ctrl_P = VectorX::Zero(1), VectorX ctrl_D = VectorX::Zero(1));
    ActuatorMotor(std::string name, Joint* joint, ControlMode control_mode = ControlMode::FORCE,
                    dtype ctrl_min = (dtype)INT_MIN, dtype ctrl_max = (dtype)INT_MAX, 
                    dtype ctrl_P = 0., dtype ctrl_D = 0.);

    void update_dofs(const VectorX& dofs, const VectorX& dofs_vel);

    void computeForce(VectorX& fm, VectorX& fr);
    void computeForceWithDerivative(VectorX& fm, VectorX& fr, MatrixX& Km, MatrixX& Dm, MatrixX& Kr, MatrixX& Dr);

    void compute_dfdu(MatrixX& dfm_du, MatrixX& dfr_du);
};

}