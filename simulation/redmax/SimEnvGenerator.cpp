#include "SimEnvGenerator.h"
#include "Robot.h"
#include "Body/BodyCuboid.h"
#include "Body/BodyMeshObj.h"
#include "Body/BodySphere.h"
#include "Joint/JointFixed.h"
#include "Joint/JointRevolute.h"
#include "Joint/JointPrismatic.h"
#include "Joint/JointFree2D.h"
#include "Joint/JointTranslational.h"
#include "Joint/JointSphericalEuler.h"
#include "Joint/JointFree3DEuler.h"
#include "Joint/JointFree3DExp.h"
#include "Force/ForceCable.h"
#include "Force/ForceGroundContact.h"
#include "Force/ForceCuboidCuboidContact.h"
#include "Force/ForceGeneralPrimitiveContact.h"
#include "Force/ForceGeneralSDFContact.h"
#include "Actuator/Actuator.h"
#include "Actuator/ActuatorMotor.h"
#include "Utils.h"
#include "Common.h"

namespace redmax {

Simulation* SimEnvGenerator::createSinglePendulumTest(std::string integrator) {
    // simulation options
    Simulation::Options *options = new Simulation::Options();
    options->_gravity = -980. * Vector3::UnitZ();
    options->_h = 0.01;
    options->_integrator = integrator;

    // construct simulation
    Simulation* sim = new Simulation(options, "single pendulum");

    // define robot
    Robot* robot = new Robot();
    Vector3 length;
    length(0) = 10; length(1) = 1; length(2) = 1;
    Vector3 pos;
    pos(0) = 5; pos(1) = 0; pos(2) = 0;
    Vector3 axis;
    axis(0) = 0; axis(1) = 1; axis(2) = 0;

    Joint* joint = new JointRevolute(sim, 0, axis, nullptr, Matrix3::Identity(), Vector3::Zero());
    joint->set_damping(1e3);

    BodyCuboid* body = new BodyCuboid(sim, joint, length, Matrix3::Identity(), pos, (dtype)1.0);

    robot->_root_joints.push_back(joint);

    robot->init();

    sim->addRobot(robot);
    sim->init();
    
    return sim;
}

Simulation* SimEnvGenerator::createSinglePendulumObjTest(std::string integrator) {
    // simulation options
    Simulation::Options *options = new Simulation::Options();
    options->_gravity = -980. * Vector3::UnitZ();
    options->_h = 0.01;
    options->_integrator = integrator;

    // construct simulation
    Simulation* sim = new Simulation(options, "single pendulum");

    // define robot
    Robot* robot = new Robot();
    Vector3 length;
    length(0) = 10; length(1) = 1; length(2) = 1;
    Vector3 pos;
    pos(0) = 5; pos(1) = 0; pos(2) = 0;
    Vector3 axis;
    axis(0) = 0; axis(1) = 1; axis(2) = 0;

    Joint* joint = new JointRevolute(sim, 0, axis, nullptr, Matrix3::Identity(), Vector3::Zero());
    joint->set_damping(1e3);

    // BodyMeshObj* body = new BodyMeshObj(sim, joint, 
    //     std::string(GRAPHICS_CODEBASE_SOURCE_DIR) + "/resources/meshes/cube.obj",
    //     Matrix3::Identity(), Vector3::Zero(), BodyMeshObj::OBJ_TO_JOINT,
    //     (dtype)1.0, (dtype)0.0);
    // BodyMeshObj* body = new BodyMeshObj(sim, joint, 
    //     std::string(GRAPHICS_CODEBASE_SOURCE_DIR) + "/resources/meshes/cube_1.obj",
    //     Matrix3::Identity(), Vector3::Zero(), BodyMeshObj::OBJ_TO_WOLRD,
    //     (dtype)1.0, (dtype)0.0);
    BodyMeshObj* body = new BodyMeshObj(sim, joint, 
        std::string(GRAPHICS_CODEBASE_SOURCE_DIR) + "/resources/meshes/cube.obj",
        Matrix3::Identity(), pos, BodyMeshObj::BODY_TO_JOINT,
        (dtype)1.0);

    robot->_root_joints.push_back(joint);

    robot->init();

    sim->addRobot(robot);
    sim->init();
    
    return sim;
}

Simulation* SimEnvGenerator::createMultiPendulumTest(int num_links, std::string integrator) {
    // simulation options
    Simulation::Options *options = new Simulation::Options();
    options->_gravity = -980. * Vector3::UnitZ();
    options->_h = 0.01;
    options->_integrator = integrator;
    
    // construct simulation
    Simulation* sim = new Simulation(options, "multi pendulum");

    Robot* robot = new Robot();
    Vector3 length; length(0) = 10; length(1) = 1; length(2) = 1;
    Vector3 body_pos; body_pos(0) = 5; body_pos(1) = 0; body_pos(2) = 0;
    Vector3 joint_pos; joint_pos(0) = 10; joint_pos(1) = 0; joint_pos(2) = 0;

    vector<Joint*> joints;
    joints.clear();
    for (int i = 0;i < num_links;++i) {
        Joint* joint;
        if (i == 0) {
            joint = new JointRevolute(sim, i, Vector3::UnitY(), nullptr, Matrix3::Identity(), Vector3::Zero());
        } else {
            joint = new JointRevolute(sim, i, Vector3::UnitY(), joints[i - 1], Matrix3::Identity(), joint_pos);
        }
        joint->set_damping((dtype)1e2);
        joints.push_back(joint);

        BodyCuboid* body = new BodyCuboid(sim, joint, length, Matrix3::Identity(), body_pos, (dtype)1.0);

    }
    
    robot->_root_joints.push_back(joints[0]);

    robot->init();
    
    sim->addRobot(robot);
    sim->init();
    
    return sim;
}

Simulation* SimEnvGenerator::createPrismaticTest(std::string integrator) {
    // simulation options
    Simulation::Options *options = new Simulation::Options();
    options->_gravity = -980. * Vector3::UnitZ();
    options->_h = 0.01;
    options->_integrator = integrator;

    // construct simulation
    Simulation* sim = new Simulation(options, "prismatic test");

    Robot* robot = new Robot();
    
    Joint* joint0 = new JointPrismatic(sim, 0, Vector3::UnitX(), nullptr, Matrix3::Identity(), Vector3::Zero());
    BodyCuboid* body0 = new BodyCuboid(sim, joint0, Vector3(20, 1, 1), Matrix3::Identity(), Vector3::Zero(), 1.);
    Joint* joint1 = new JointRevolute(sim, 1, Vector3::UnitY(), joint0, Matrix3::Identity(), Vector3(-10, 0, 0));
    joint1->_q(0) = constants::pi / 2.;
    BodyCuboid* body1 = new BodyCuboid(sim, joint1, Vector3(1, 1, 10), Matrix3::Identity(), Vector3(0, 0, -5.), 1.);

    robot->_root_joints.push_back(joint0);

    robot->init();

    sim->addRobot(robot);
    sim->init();

    return sim;
}

Simulation* SimEnvGenerator::createCableTest(std::string integrator) {
    // simulation options
    Simulation::Options *options = new Simulation::Options();
    options->_gravity = -980. * Vector3::UnitZ();
    options->_h = 0.01;
    options->_integrator = integrator;

    // construct simulation
    Simulation* sim = new Simulation(options, "cable test");

    Robot* robot = new Robot();

    Joint* joint0 = new JointRevolute(sim, 0, Vector3::UnitY(), nullptr, Matrix3::Identity(), Vector3::Zero());
    joint0->_q(0) = constants::pi / 2.;
    BodyCuboid* body0 = new BodyCuboid(sim, joint0, Vector3(10, 1, 1), Matrix3::Identity(), Vector3(5, 0, 0), 1);
    
    Joint* joint1 = new JointRevolute(sim, 1, Vector3::UnitY(), joint0, Matrix3::Identity(), Vector3(10, 0, 0));
    joint1->_q(0) = -constants::pi / 2.;
    BodyCuboid* body1 = new BodyCuboid(sim, joint1, Vector3(10, 1, 1), Matrix3::Identity(), Vector3(5, 0, 0), 1);
    
    Joint* joint2 = new JointPrismatic(sim, 2, Vector3::UnitX(), nullptr, Matrix3::Identity(), Vector3(10, 0, 0));
    joint2->set_stiffness(1e4, Vector1::Zero());
    joint2->set_damping(1e3);
    BodyCuboid* body2 = new BodyCuboid(sim, joint2, Vector3(1, 1, 1), Matrix3::Identity(), Vector3(0, 0, 0), 1);
    

    robot->_root_joints.push_back(joint0);
    robot->_root_joints.push_back(joint2);

    ForceCable* force = new ForceCable(sim);
    force->set_stiffness(1e6);
    force->set_damping(1e3);
    force->addBodyPoint(body2, Vector3::Zero());
    force->addBodyPoint(body0, Vector3(-4, 0, 1));
    force->addBodyPoint(body1, Vector3(-4, 0, 1));
    robot->add_force(force);

    robot->init();

    sim->addRobot(robot);
    sim->init();
    
    return sim;
}

Simulation* SimEnvGenerator::createFree2DTest(std::string integrator) {
    // simulation options
    Simulation::Options *options = new Simulation::Options();
    options->_gravity = -980. * Vector3::UnitY();
    options->_h = 5e-3;
    options->_integrator = integrator;

    // viewer options
    Simulation::ViewerOptions *viewer_options = new Simulation::ViewerOptions();
    viewer_options->_fps = 60;
    viewer_options->_speed = 0.5;
    viewer_options->_camera_up = Vector3::UnitY();
    viewer_options->_camera_pos = Vector3(0., 0., 2.);
    viewer_options->_ground = false;

    // construct simulation
    Simulation* sim = new Simulation(options, viewer_options, "Free2D joint test");

    Robot* robot = new Robot();

    Joint* joint0 = new JointFree2D(sim, 0, nullptr, Matrix3::Identity(), Vector3::Zero());
    joint0->_q << -10, -10, 0;
    joint0->_qdot << 50, 200, 20;
    BodyCuboid* body0 = new BodyCuboid(sim, joint0, Vector3(1, 1, 1), Matrix3::Identity(), Vector3::Zero(), 1.);
    
    robot->_root_joints.push_back(joint0);
    
    robot->init();

    sim->addRobot(robot);
    sim->init();

    return sim;
}

Simulation* SimEnvGenerator::createFree3DEulerTest(std::string integrator) {
    // simulation options
    Simulation::Options *options = new Simulation::Options();
    options->_h = 5e-3;
    options->_integrator = integrator;

    // viewer options
    Simulation::ViewerOptions *viewer_options = new Simulation::ViewerOptions();
    viewer_options->_fps = 60;
    viewer_options->_speed = 0.5;
    viewer_options->_ground = true;
    viewer_options->_E_g.setIdentity();

    // construct simulation
    Simulation* sim = new Simulation(options, viewer_options, "Free3D(euler) joint test");

    Robot* robot = new Robot();

    // Joint* joint0 = new JointSphericalEuler(0, nullptr, Matrix3::Identity(), Vector3(0, 0, 6.));
    // Joint* joint0 = new JointTranslational(0, nullptr, Matrix3::Identity(), Vector3(0, 0, 6.));
    Joint* joint0 = new JointFree3DEuler(sim, 0, nullptr, Matrix3::Identity(), Vector3(0, 0, 6.));
    // Joint* joint0 = new JointFixed(0, nullptr, Matrix3::Identity(), Vector3(0, 0, 6.));
    joint0->_qdot << 30., 0., 10., 5., 0., 0;
    BodyCuboid* body0 = new BodyCuboid(sim, joint0, Vector3(4, 4, 4), Matrix3::Identity(), Vector3::Zero(), 1.);
    
    Matrix4 E_g = Matrix4::Identity();
    ForceGroundContact* force0 = new ForceGroundContact(sim, body0, E_g, 1e6, 1e3, 0.2, 3e1);
    robot->add_force(force0);

    robot->_root_joints.push_back(joint0);
    
    robot->init();

    sim->addRobot(robot);
    sim->init();

    return sim;
}

Simulation* SimEnvGenerator::createFree3DExpTest(std::string integrator) {
    // simulation options
    Simulation::Options *options = new Simulation::Options();
    options->_h = 5e-3;
    options->_integrator = integrator;

    // viewer options
    Simulation::ViewerOptions *viewer_options = new Simulation::ViewerOptions();
    viewer_options->_fps = 60;
    viewer_options->_speed = 0.5;
    viewer_options->_ground = true;
    viewer_options->_E_g.setIdentity();

    // construct simulation
    Simulation* sim = new Simulation(options, viewer_options, "Free3D(euler) joint test");

    Robot* robot = new Robot();

    // Joint* joint0 = new JointSphericalEuler(0, nullptr, Matrix3::Identity(), Vector3(0, 0, 6.));
    // Joint* joint0 = new JointTranslational(0, nullptr, Matrix3::Identity(), Vector3(0, 0, 6.));
    Joint* joint0 = new JointFree3DExp(sim, 0, nullptr, Matrix3::Identity(), Vector3(0, 0, 6.));
    // Joint* joint0 = new JointFixed(0, nullptr, Matrix3::Identity(), Vector3(0, 0, 6.));
    joint0->_qdot << 30., 0., 10., 5., 0., 0;
    BodyCuboid* body0 = new BodyCuboid(sim, joint0, Vector3(4, 4, 4), Matrix3::Identity(), Vector3::Zero(), 1.);
    
    Matrix4 E_g = Matrix4::Identity();
    ForceGroundContact* force0 = new ForceGroundContact(sim, body0, E_g, 1e6, 1e3, 0.2, 3e1);
    robot->add_force(force0);

    robot->_root_joints.push_back(joint0);
    
    robot->init();

    sim->addRobot(robot);
    sim->init();

    return sim;
}

Simulation* SimEnvGenerator::createGroundContactTest(std::string integrator) {
    // simulation options
    Simulation::Options *options = new Simulation::Options();
    options->_gravity = -980. * Vector3::UnitY();
    options->_h = 5e-4;
    options->_integrator = integrator;

    // viewer options
    Simulation::ViewerOptions *viewer_options = new Simulation::ViewerOptions();
    viewer_options->_fps = 30;
    viewer_options->_speed = 0.5;
    viewer_options->_camera_up = Vector3::UnitY();
    viewer_options->_camera_pos = Vector3(0.75, 0.75, 1.25);
    viewer_options->_ground = true;
    viewer_options->_E_g.topLeftCorner(3, 3) = Eigen::AngleAxis<dtype>(-constants::pi / 2., Vector3::UnitX()).matrix();
    viewer_options->_E_g.topRightCorner(3, 1).setZero();
    // viewer_options->_record = true;
    // viewer_options->_record_folder = "";

    // construct simulation
    Simulation* sim = new Simulation(options, viewer_options, "Free2D joint with ground");

    Robot* robot = new Robot();

    Joint* joint0 = new JointFree2D(sim, 0, nullptr, Matrix3::Identity(), Vector3::Zero());
    joint0->_q << -2, 2, 0;
    joint0->_qdot << 5, 70, 5;
    
    BodyCuboid* body0 = new BodyCuboid(sim, joint0, Vector3(3, 1, 1), Matrix3::Identity(), Vector3::Zero(), 1.);
    // BodySphere* body0 = new BodySphere(sim, joint0, 1., Matrix3::Identity(), Vector3::Zero(), 1.);

    Matrix3 R = Eigen::AngleAxis<dtype>(-constants::pi / 2., Vector3::UnitX()).matrix();
    Matrix4 E_g = Matrix4::Identity();
    E_g.topLeftCorner(3, 3) = R;
    ForceGroundContact* force0 = new ForceGroundContact(sim, body0, E_g, 1e6, 1e4, 0.8, 3e1);
    robot->add_force(force0);

    robot->_root_joints.push_back(joint0);
    
    robot->init();

    sim->addRobot(robot);
    sim->init();

    return sim;
}

Simulation* SimEnvGenerator::createBoxContactTest(std::string integrator) {
    // simulation options
    Simulation::Options *options = new Simulation::Options();
    options->_gravity = -980. * Vector3::UnitY();
    options->_h = 5e-4;
    options->_integrator = integrator;

    // viewer options
    Simulation::ViewerOptions *viewer_options = new Simulation::ViewerOptions();
    viewer_options->_fps = 60;
    viewer_options->_speed = 0.5;
    viewer_options->_camera_up = Vector3::UnitY();
    viewer_options->_camera_pos = Vector3(0.8, 0.75, 1.25);
    // viewer_options->_camera_pos = Vector3(0.5, 0.45, 0.8);
    // viewer_options->_camera_lookat = Vector3(0.2, 0., 0.);
    viewer_options->_ground = true;
    viewer_options->_E_g.topLeftCorner(3, 3) = Eigen::AngleAxis<dtype>(-constants::pi / 2., Vector3::UnitX()).matrix();
    viewer_options->_E_g.topRightCorner(3, 1).setZero();
    // viewer_options->_record = true;
    // viewer_options->_record_folder = "";

    // construct simulation
    Simulation* sim = new Simulation(options, viewer_options, "Free2D Box-Box contact with ground");

    Robot* robot = new Robot();

    Joint* joint0 = new JointFree2D(sim, 0, nullptr, Matrix3::Identity(), Vector3(0.0, 0.5, 0.0));
    // joint0->_qdot << 35., 35., 0.5;
    joint0->_qdot << 45., 0., 0.;
    BodyCuboid* body0 = new BodyCuboid(sim, joint0, Vector3(1, 1, 1), Matrix3::Identity(), Vector3::Zero(), 1.);
    
    Joint* joint1 = new JointFree2D(sim, 1, nullptr, Matrix3::Identity(), Vector3(3.0, 0.5, 0.0));
    BodyCuboid* body1 = new BodyCuboid(sim, joint1, Vector3(1, 1, 1), Matrix3::Identity(), Vector3::Zero(), 1.);
    
    robot->_root_joints.push_back(joint0);
    robot->_root_joints.push_back(joint1);

    // add ground contact
    Matrix3 R = Eigen::AngleAxis<dtype>(-constants::pi / 2., Vector3::UnitX()).matrix();
    Matrix4 E_g = Matrix4::Identity();
    E_g.topLeftCorner(3, 3) = R;
    ForceGroundContact* force0 = new ForceGroundContact(sim, body0, E_g, 1e5, 1e2, 0.2, 3e1);
    robot->add_force(force0);
    ForceGroundContact* force1 = new ForceGroundContact(sim, body1, E_g, 1e5, 1e2, 0.2, 3e1);
    robot->add_force(force1);

    // add box-box contact
    ForceCuboidCuboidContact* force2 = new ForceCuboidCuboidContact(sim, body0, body1, 1e4, 0.0);
    robot->add_force(force2);

    robot->init();

    sim->addRobot(robot);
    sim->init();

    return sim;
}

// TorqueFingerDemo
// contains a torque-driven finger
// goal: the finger to achieve target pose
Simulation* SimEnvGenerator::createTorqueFingerDemo(std::string integrator) {
    // simulation options
    Simulation::Options *options = new Simulation::Options();
    options->_gravity = -980. * Vector3::UnitY();
    options->_h = 1e-2;
    options->_integrator = integrator;

    // viewer options
    Simulation::ViewerOptions *viewer_options = new Simulation::ViewerOptions();
    viewer_options->_fps = 30;
    viewer_options->_speed = 0.5;
    viewer_options->_camera_up = Vector3::UnitY();
    viewer_options->_camera_pos = Vector3(0.6, 1.3, 2.0);
    // viewer_options->_camera_pos = Vector3(0.5, 0.45, 0.8);
    viewer_options->_camera_lookat = Vector3(0.1, 0.8, 0.);
    viewer_options->_ground = true;
    viewer_options->_E_g.topLeftCorner(3, 3) = Eigen::AngleAxis<dtype>(-constants::pi / 2., Vector3::UnitX()).matrix();
    viewer_options->_E_g.topRightCorner(3, 1).setZero();
    // viewer_options->_record = true;
    // viewer_options->_record_folder = "";

    // construct simulation
    Simulation* sim = new Simulation(options, viewer_options, "Torque-driven finger demo");

    Robot* robot = new Robot();

    // MCP
    Joint* joint0 = new JointRevolute(sim, 0, Vector3(0, 0, -1), nullptr, Matrix3::Identity(), Vector3(0., 10., 0.));
    joint0->set_damping(1e4);
    BodyCuboid* body0 = new BodyCuboid(sim, joint0, Vector3(4, 1, 1), Matrix3::Identity(), Vector3(2., 0., 0.), 1.);
    Actuator* actuator0 = new ActuatorMotor("actuator", joint0, ActuatorMotor::ControlMode::FORCE, -1e5, 1e5);
    robot->add_actuator(actuator0);

    // PIP
    Joint* joint1 = new JointRevolute(sim, 1, Vector3(0, 0, -1), joint0, Matrix3::Identity(), Vector3(4., 0., 0.));
    joint1->set_damping(1e4);
    BodyCuboid* body1 = new BodyCuboid(sim, joint1, Vector3(2, 1, 1), Matrix3::Identity(), Vector3(1., 0., 0.), 1.);
    Actuator* actuator1 = new ActuatorMotor("actuator", joint1, ActuatorMotor::ControlMode::FORCE, -1e5, 1e5);
    robot->add_actuator(actuator1);

    // DIP
    Joint* joint2 = new JointRevolute(sim, 2, Vector3(0, 0, -1), joint1, Matrix3::Identity(), Vector3(2., 0., 0.));
    joint2->set_damping(1e4);
    BodyCuboid* body2 = new BodyCuboid(sim, joint2, Vector3(1, 1, 1), Matrix3::Identity(), Vector3(0.5, 0., 0.), 1.);
    Actuator* actuator2 = new ActuatorMotor("actuator", joint2, ActuatorMotor::ControlMode::FORCE, -1e5, 1e5);
    robot->add_actuator(actuator2);

    robot->_root_joints.push_back(joint0);

    robot->init();

    sim->addRobot(robot);
    sim->init();

    return sim;
}

// TorqueFingerFlickDemo
// contains a torque-driven finger and a box with 2D free joint
// goal: the finger to flick the box to target position
Simulation* SimEnvGenerator::createTorqueFingerFlickDemo(std::string integrator) {
    // simulation options
    Simulation::Options *options = new Simulation::Options();
    options->_gravity = -980. * Vector3::UnitY();
    options->_h = 5e-4;
    options->_integrator = integrator;

    // viewer options
    Simulation::ViewerOptions *viewer_options = new Simulation::ViewerOptions();
    viewer_options->_fps = 30;
    viewer_options->_speed = 0.5;
    viewer_options->_camera_up = Vector3::UnitY();
    viewer_options->_camera_pos = Vector3(0.8, 0.6, 2.0);
    // viewer_options->_camera_pos = Vector3(0.5, 0.45, 0.8);
    viewer_options->_camera_lookat = Vector3(0.3, 0.4, 0.);
    viewer_options->_ground = true;
    viewer_options->_E_g.topLeftCorner(3, 3) = Eigen::AngleAxis<dtype>(-constants::pi / 2., Vector3::UnitX()).matrix();
    viewer_options->_E_g.topRightCorner(3, 1).setZero();
    // viewer_options->_record = true;
    // viewer_options->_record_folder = "";

    // construct simulation
    Simulation* sim = new Simulation(options, viewer_options, "Torque-driven finger flick box demo");

    Robot* robot = new Robot();

    // MCP
    Joint* joint0 = new JointRevolute(sim, 0, Vector3(0, 0, -1), nullptr, Matrix3::Identity(), Vector3(0., 3.2, 0.));
    joint0->set_damping(1e4);
    BodyCuboid* body0 = new BodyCuboid(sim, joint0, Vector3(4, 1, 1), Matrix3::Identity(), Vector3(2., 0., 0.), 1.);
    Actuator* actuator0 = new ActuatorMotor("actuator", joint0, ActuatorMotor::ControlMode::FORCE, -1e5, 1e5);
    robot->add_actuator(actuator0);

    // PIP
    Joint* joint1 = new JointRevolute(sim, 1, Vector3(0, 0, -1), joint0, Matrix3::Identity(), Vector3(4., 0., 0.));
    joint1->set_damping(1e4);
    BodyCuboid* body1 = new BodyCuboid(sim, joint1, Vector3(2, 1, 1), Matrix3::Identity(), Vector3(1., 0., 0.), 1.);
    Actuator* actuator1 = new ActuatorMotor("actuator", joint1, ActuatorMotor::ControlMode::FORCE, -1e5, 1e5);
    robot->add_actuator(actuator1);

    // DIP
    Joint* joint2 = new JointRevolute(sim, 2, Vector3(0, 0, -1), joint1, Matrix3::Identity(), Vector3(2., 0., 0.));
    joint2->set_damping(1e4);
    BodyCuboid* body2 = new BodyCuboid(sim, joint2, Vector3(1, 1, 1), Matrix3::Identity(), Vector3(0.5, 0., 0.), 1.);
    Actuator* actuator2 = new ActuatorMotor("actuator", joint2, ActuatorMotor::ControlMode::FORCE, -1e5, 1e5);
    robot->add_actuator(actuator2);

    // Box
    Joint* box_joint = new JointFree2D(sim, 0, nullptr, Matrix3::Identity(), Vector3(5.2, 0.5, 0.0));
    BodyCuboid* box = new BodyCuboid(sim, box_joint, Vector3(1, 1, 1), Matrix3::Identity(), Vector3::Zero(), 1.);
    box->set_color(Vector3(0., 0.2, 0.4));

    // virtual goal
    JointFixed* goal_joint = new JointFixed(sim, 10, nullptr, Matrix3::Identity(), Vector3(5.2 + 10.5, 0.5, 0.));
    BodyCuboid* goal = new BodyCuboid(sim, goal_joint, Vector3(0.2, 0.2, 0.2), Matrix3::Identity(), Vector3::Zero(), 1.);
    goal->set_color(Vector3(0., 1., 0.));

    robot->_root_joints.push_back(joint0);
    robot->_root_joints.push_back(box_joint);
    robot->_root_joints.push_back(goal_joint);

    // add ground contact
    Matrix3 R = Eigen::AngleAxis<dtype>(-constants::pi / 2., Vector3::UnitX()).matrix();
    Matrix4 E_g = Matrix4::Identity();
    E_g.topLeftCorner(3, 3) = R;
    ForceGroundContact* force_ground_0 = new ForceGroundContact(sim, body0, E_g, 1e4, 1e2, 0.2, 3e1);
    robot->add_force(force_ground_0);
    ForceGroundContact* force_ground_1 = new ForceGroundContact(sim, body1, E_g, 1e4, 1e2, 0.2, 3e1);
    robot->add_force(force_ground_1);
    ForceGroundContact* force_ground_2 = new ForceGroundContact(sim, body2, E_g, 1e4, 1e2, 0.2, 3e1);
    robot->add_force(force_ground_2);
    ForceGroundContact* force_ground_3 = new ForceGroundContact(sim, box, E_g, 1e4, 1e2, 0.2, 3e1);
    robot->add_force(force_ground_3);

    // add box-box contact
    ForceCuboidCuboidContact* force0 = new ForceCuboidCuboidContact(sim, body0, box, 1e4, 0.0);
    robot->add_force(force0);

    ForceCuboidCuboidContact* force1 = new ForceCuboidCuboidContact(sim, body1, box, 1e4, 0.0);
    robot->add_force(force1);

    ForceCuboidCuboidContact* force2 = new ForceCuboidCuboidContact(sim, body2, box, 1e4, 0.0);
    robot->add_force(force2);

    robot->init();

    sim->addRobot(robot);
    sim->init();

    return sim;
}

Simulation* SimEnvGenerator::createSphereGroundContactTest(std::string integrator) {
    // simulation options
    Simulation::Options *options = new Simulation::Options();
    options->_h = 5e-4;
    options->_integrator = integrator;

    // viewer options
    Simulation::ViewerOptions *viewer_options = new Simulation::ViewerOptions();
    viewer_options->_ground = true;
    viewer_options->_E_g.topRightCorner(3, 1).setZero();

    // construct simulation
    Simulation* sim = new Simulation(options, viewer_options, "Free2D joint with ground");

    Robot* robot = new Robot();

    Matrix3 joint0_R;
    joint0_R.col(0) = Vector3::UnitX();
    joint0_R.col(1) = -Vector3::UnitZ();
    joint0_R.col(2) = Vector3::UnitY();
    Joint* joint0 = new JointFree2D(sim, 0, nullptr, joint0_R, Vector3(0, 0, 1.));
    joint0->_qdot << 5., 0., 0.;
    // Joint* joint0 = new JointFree2D(0, nullptr, joint0_R, Vector3(0, 0, 0.5));
    // joint0->_qdot << 50., 0., 0;
    
    // BodyCuboid* body0 = new BodyCuboid(sim, joint0, Vector3(1, 1, 1), Matrix3::Identity(), Vector3::Zero(), 1.);
    BodySphere* body0 = new BodySphere(sim, joint0, 1., Matrix3::Identity(), Vector3::Zero(), 1.);
    
    Matrix4 E_g = Matrix4::Identity();
    ForceGroundContact* force0 = new ForceGroundContact(sim, body0, E_g, 1e5, 1e3, 0.3, 3e1);
    robot->add_force(force0);

    robot->_root_joints.push_back(joint0);
    
    robot->init();

    sim->addRobot(robot);
    sim->init();

    return sim;
}

}