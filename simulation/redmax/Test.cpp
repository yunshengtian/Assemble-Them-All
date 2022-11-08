#include "Test.h"
#include <Eigen/Dense>
#include "Robot.h"
#include "Body/Body.h"
#include "Utils.h"
#include "SimEnvGenerator.h"
#include "Joint/Joint.h"
#include "Body/BodyCuboid.h"
#include "Body/BodySphere.h"
#include "Body/BodySDFObj.h"
#include "Body/BodyBVHObj.h"
#include "Body/BodyCylinder.h"
#include "Joint/JointFixed.h"
#include "Joint/JointSphericalExp.h"
#include "Force/ForceGeneralPrimitiveContact.h"
#include "Force/ForceGeneralSDFContact.h"
#include "Force/ForceGeneralSDFContactWithSuction.h"
#include "Force/ForceGeneralBVHContact.h"

using namespace redmax;

Eigen::VectorXd Test::random_generated_vector(int len, Eigen::VectorXd lower_bound, Eigen::VectorXd upper_bound) {
    Eigen::VectorXd res(len);
    for (int i = 0;i < len;++i) {
        double a = (double)rand() / RAND_MAX;
        res[i] = a * lower_bound[i] + (1.0 - a) * upper_bound[i];
    }
    return res;
}

Eigen::VectorXd Test::random_generated_vector(int len, double lower_bound, double upper_bound) {
    Eigen::VectorXd res(len);
    for (int i = 0;i < len;++i) {
        double a = (double)rand() / RAND_MAX;
        res[i] = a * lower_bound + (1.0 - a) * upper_bound;
    }
    return res;
}

Eigen::MatrixXd Test::random_generated_matrix(int n, int m, double lower_bound, double upper_bound) {
    Eigen::MatrixXd res(n, m);
    for (int i = 0;i < n;++i)
        for (int j = 0;j < m;++j) {
            double a = (double)rand() / RAND_MAX;
            res(i, j) = a * lower_bound + (1.0 - a) * upper_bound;
        }
    return res;
}

void Test::test_xml_reader() {
    Simulation* sim = new Simulation("../../assets/assemble_box_sdf.xml");

    sim->_viewer_options->_loop = true;
    sim->_viewer_options->_infinite = true;

    sim->init();

    sim->print_ctrl_info();

    sim->print_design_params_info();
    
    auto qdot_init = sim->get_qdot_init();
    qdot_init[sim->_ndof_r - 2] = 5.;
    sim->set_qdot_init(qdot_init);

    sim->reset();

    int t0 = clock();
    sim->forward(300, false, false);

    int t1 = clock();

    std::cerr << "time = " << (t1 - t0) / 1000.0 << "ms" << std::endl;

    sim->replay();
}

void Test::test_speed() {
    Simulation* sim = new Simulation("../../assets/assemble_box_sdf.xml");

    sim->init();

    sim->reset();

    int t0 = clock();
    for (int i = 0;i < 5000;i++)
        sim->forward(1, false, false);

    int t1 = clock();

    std::cerr << "time = " << (t1 - t0) / 1000.0 << "ms" << std::endl;
}

// for passive simulation without actuators
void Test::test_forward() {
    // Simulation* sim = SimEnvGenerator::createSinglePendulumTest();
    // Simulation* sim = SimEnvGenerator::createMultiPendulumTest(4);
    // Simulation* sim = SimEnvGenerator::createPrismaticTest();
    // Simulation* sim = SimEnvGenerator::createCableTest();
    // Simulation* sim = SimEnvGenerator::createFree2DTest();
    // Simulation* sim = SimEnvGenerator::createGroundContactTest("BDF1");
    // Simulation* sim = SimEnvGenerator::createBoxContactTest("BDF1");
    // Simulation* sim = SimEnvGenerator::createSinglePendulumObjTest();
    // Simulation* sim = SimEnvGenerator::createSphereGroundContactTest("BDF2");
    // Simulation* sim = SimEnvGenerator::createFree3DEulerTest();
    Simulation* sim = SimEnvGenerator::createFree3DExpTest();

    sim->_viewer_options->_loop = true;
    sim->_viewer_options->_infinite = true;

    sim->init();
    
    int num_steps = 1000;

    int t0 = clock();
    sim->reset();
    sim->forward(num_steps, false);
    int t1 = clock();

    std::cerr << "time = " << (t1 - t0) / 1000.0 << "ms" << std::endl;

    std::cerr << "q = " << sim->get_q().transpose() << ", qdot = " << sim->get_qdot().transpose() << std::endl;

    sim->replay();
}

void Test::test_derivative() {
    // Simulation* sim = SimEnvGenerator::createSinglePendulumTest();
    // Simulation* sim = SimEnvGenerator::createMultiPendulumTest(4);
    // Simulation* sim = SimEnvGenerator::createPrismaticTest();
    // Simulation* sim = SimEnvGenerator::createCableTest();
    // Simulation* sim = SimEnvGenerator::createFree2DTest();
    // Simulation* sim = SimEnvGenerator::createGroundContactTest();
    // Simulation* sim = SimEnvGenerator::createBoxContactTest();
    // Simulation* sim = new Simulation("../../assets/box_stack.xml");

    // sim->init();

    // srand(time(0));
    // sim->set_q(Test::random_generated_vector(sim->_ndof_r, -1., 1.));
    // sim->set_qdot(Test::random_generated_vector(sim->_ndof_r, -1., 1.));

    // sim->test_derivative();

    // // test JointSphericalExp Derivatives
    // {
    //     JointSphericalExp* joint = new JointSphericalExp(sim, 0, nullptr, Matrix3::Identity(), Vector3::Zero());
    //     joint->test_derivatives();
    //     delete joint;
    // }

    // // test BodySphere Derivatives
    // {
    //     JointFixed* joint = new JointFixed(sim, 0, nullptr, Matrix3::Identity(), Vector3::Zero());
    //     BodySphere* body_sphere = new BodySphere(nullptr, joint, 2., Matrix3::Identity(), Vector3::Zero(), 1.);
    //     body_sphere->test_collision_derivatives();
    //     delete joint;
    //     delete body_sphere;
    // }

    // // test BodyCuboid Derivatives
    // {
    //     JointFixed* joint = new JointFixed(sim, 0, nullptr, Matrix3::Identity(), Vector3::Zero());
    //     BodyCuboid* body_cuboid = new BodyCuboid(nullptr, joint, Vector3(1., 1., 1.), Matrix3::Identity(), Vector3::Zero(), 1.);
    //     body_cuboid->test_collision_derivatives();
    //     delete joint;
    //     delete body_cuboid;
    // }

    // test general-primitive contact Derivatives
    // {
    //     Simulation* sim = new Simulation("../../assets/finger_mesh.xml");

    //     dynamic_cast<ForceGeneralPrimitiveContact*>(sim->_robot->_forces[2])->test_derivatives();
        
    //     delete sim;
    // }

    // // test BodySDFObj Derivatives
    // {
    //     JointFixed* joint = new JointFixed(sim, 0, nullptr, Matrix3::Identity(), Vector3::Zero());
    //     BodySDFObj* body_sdf = new BodySDFObj(sim, joint, "../../assets/box/cube.obj", Matrix3::Identity(), Vector3::Zero(), 0.05, 0.01, BodyMeshObj::OBJ_TO_JOINT);
    //     body_sdf->test_collision_derivatives();
    //     delete joint;
    //     delete body_sdf;
    // }

    // // test BodyCylinder Derivatives
    // {
    //     JointFixed* joint = new JointFixed(sim, 0, nullptr, Matrix3::Identity(), Vector3::Zero());
    //     BodyCylinder* body_cylinder = new BodyCylinder(nullptr, joint, 1.0, 0.5, Matrix3::Identity(), Vector3::Zero(), 1.);
    //     body_cylinder->test_collision_derivatives();
    //     delete joint;
    //     delete body_cylinder;
    // }

    // test general-BVH contact Derivatives
    // {
    //     Simulation* sim = new Simulation("../../assets/box_stack_bvh.xml");
    //     dynamic_cast<ForceGeneralBVHContact*>(sim->_robot->_forces[2])->test_derivatives();

    //     JointFixed* joint = new JointFixed(sim, 0, nullptr, Matrix3::Identity(), Vector3::Zero());
    //     BodyBVHObj* body_bvh = new BodyBVHObj(sim, joint, "../../assets/box/cube.obj", "../../assets/box/cube.obj", Matrix3::Identity(), Vector3::Zero(), 0.05, 0.01, BodyMeshObj::OBJ_TO_JOINT);
    //     body_bvh->test_collision_derivatives();
    //     delete joint;
    //     delete body_bvh;

    //     delete sim;
    // }

    // // test general-SDF contact with suction Derivatives
    // {
    //     Simulation* sim = new Simulation("../../assets/box_lift.xml");
    //     dynamic_cast<ForceGeneralSDFContactWithSuction*>(sim->_robot->_forces[2])->test_derivatives();
    //     delete sim;
    // }
}

void Test::test_design_derivative() {
    // Simulation* sim = new Simulation("../../assets/double_pendulum_design_1.xml");
    // Simulation* sim = new Simulation("../../assets/double_pendulum_design.xml");
    // Simulation* sim = new Simulation("../../assets/double_pendulum_contact_design.xml");
    Simulation* sim = new Simulation("../../assets/abstract_finger.xml");
    // Simulation* sim = new Simulation("../../assets/new_finger.xml");

    sim->_viewer_options->_loop = true;
    sim->_viewer_options->_infinite = true;

    sim->init();

    sim->print_ctrl_info();

    sim->print_design_params_info();

    VectorX design_params = sim->get_design_params();

    // srand(time(0));

    // for (int i = 0;i < sim->_robot->_ndof_p1 / 12;i++) {
    //     int idx = i * 12;
    //     Matrix4 E = math::compose_E(design_params.segment(idx, 12));
    //     Matrix4 dE = Matrix4::Identity();
    //     Eigen::AngleAxis<dtype> angle_axis = Eigen::AngleAxis<dtype>((dtype)rand() / RAND_MAX - 1.0, random_generated_vector(3, -1., 1.).normalized());
    //     dE.topLeftCorner(3, 3) = angle_axis.toRotationMatrix();
    //     dE.topRightCorner(3, 1) = random_generated_vector(3, -1., 1.);
    //     E = E * dE;
    //     design_params.segment(idx, 12) = math::flatten_E(E);
    // }

    // for (int i = 0;i < sim->_robot->_ndof_p2 / 12;i++) {
    //     int idx = sim->_robot->_ndof_p1 + i * 12;
    //     Matrix4 E = math::compose_E(design_params.segment(idx, 12));
    //     Matrix4 dE = Matrix4::Identity();
    //     Eigen::AngleAxis<dtype> angle_axis = Eigen::AngleAxis<dtype>((dtype)rand() / RAND_MAX - 1.0, random_generated_vector(3, -1., 1.).normalized());
    //     dE.topLeftCorner(3, 3) = angle_axis.toRotationMatrix();
    //     dE.topRightCorner(3, 1) = random_generated_vector(3, -1., 1.);
    //     E = E * dE;
    //     design_params.segment(idx, 12) = math::flatten_E(E);
    // }

    // for (int i = 0;i < sim->_robot->_ndof_p3;i++) {
    //     int idx = sim->_robot->_ndof_p1 + sim->_robot->_ndof_p2 + i;
    //     design_params(idx) += (dtype)rand() / RAND_MAX - 0.5;
    // }

    // for (int i = 0;i < sim->_robot->_ndof_p4 / 4;i++) {
    //     int idx = sim->_robot->_ndof_p1 + sim->_robot->_ndof_p2 + sim->_robot->_ndof_p3 + i * 4;
    //     for (int j = 0;j < 3;j++) {
    //         dtype ratio = (dtype)rand() / RAND_MAX + 0.5;
    //         design_params(idx + j) *= ratio;
    //     }
    //     dtype l = fabs(design_params(idx + 1) - design_params(idx + 2));
    //     dtype r = design_params(idx + 1) + design_params(idx + 2);
    //     design_params(idx + 3) = (dtype)rand() / RAND_MAX * (r - l) + l;
    // }
    
    sim->set_design_params(design_params);
    
    // VectorX q_init = sim->get_q_init();
    // VectorX qdot_init = sim->get_qdot_init();
    // q_init(3) = 0.3;
    // qdot_init(0) = 10., qdot_init(1) = -50., qdot_init(2) = 5.;
    // sim->set_q_init(q_init);
    // sim->set_qdot_init(qdot_init);

    sim->reset(true, true);
    
    int t0 = clock();
    dtype f = 0.;
    int num_steps = 100;
    for (int ii = 0;ii < num_steps;ii ++) {
        sim->forward(1, false, true);
        VectorX var = sim->get_variables();
        f += var.sum();
    }
    sim->_backward_info.set_flags(false, false, true, false);
    sim->_backward_info._df_dq0 = VectorX::Zero(sim->_ndof_r);
    sim->_backward_info._df_dqdot0 = VectorX::Zero(sim->_ndof_r);
    sim->_backward_info._df_dq = VectorX::Zero(sim->_ndof_r * num_steps);
    sim->_backward_info._df_dvar = VectorX::Ones(sim->_ndof_var * num_steps);
    sim->_backward_info._df_dp = VectorX::Zero(sim->_ndof_p);
    sim->backward();
    int t1 = clock();

    std::cerr << "time = " << (t1 - t0) / 1000.0 << "ms" << std::endl;

    design_params = sim->get_design_params();

    std::cerr << "design params = " << design_params.transpose() << std::endl;
    std::cerr << "f = " << f << std::endl;
    std::cerr << "df_dp = " << sim->_backward_results._df_dp.transpose() << std::endl;
    VectorX df_dp = sim->_backward_results._df_dp;

    // dtype eps = 1e-2;
    // for (int trial = 0;trial < 10;trial++) {
    //     VectorX df_dp_fd = VectorX::Zero(sim->_ndof_p);
    //     for (int i = 0;i < sim->_ndof_p;i++) {
    //         // if (i % 12 < 9) {
    //         //     df_dp(i) = 0;
    //         //     continue;
    //         // }
    //         VectorX design_params_pos = design_params;
    //         design_params_pos[i] += eps;
    //         sim->set_design_params(design_params_pos);
    //         sim->reset(false);
    //         dtype f_pos = 0.;
    //         for (int ii = 0;ii < num_steps;ii++) {
    //             sim->forward(1, false, false);
    //             VectorX var = sim->get_variables();
    //             f_pos += var.sum();
    //         }
    //         df_dp_fd(i) = (f_pos - f) / eps;
    //     }

    //     // std::cerr << "df_dp_fd = " << df_dp_fd.transpose() << std::endl;

    //     print_error_full("df_dp", df_dp, df_dp_fd);
    //     eps /= 10.;
    // }

    sim->replay();
}

void Test::test_control_derivative() {
    // Simulation* sim = new Simulation("../../assets/double_pendulum_design_1.xml");
    // Simulation* sim = new Simulation("../../assets/double_pendulum_contact_design.xml");
    Simulation* sim = new Simulation("../../assets/abstract_finger.xml");

    sim->_viewer_options->_loop = true;
    sim->_viewer_options->_infinite = true;

    sim->init();

    sim->print_ctrl_info();

    sim->reset(true, false);
    
    int t0 = clock();
    dtype f = 0.;
    int num_steps = 1;
    VectorX action = VectorX::Zero(num_steps * sim->_ndof_u);
    action << -3.65349773,  2.90075299,  2.12093891,  2.06099197;
    VectorX u = action;
    for (int i = 0;i < u.size();i++)
        u(i) = tanh(u(i));

    for (int ii = 0;ii < num_steps;ii ++) {
        sim->set_u(u.segment(ii * sim->_ndof_u, (ii + 1) * sim->_ndof_u));
        sim->forward(1, false, false);
        VectorX var = sim->get_variables();
        f += var.sum();
    }
    sim->_backward_info.set_flags(false, false, false, true);
    sim->_backward_info._df_dq = VectorX::Zero(sim->_ndof_r * num_steps);
    sim->_backward_info._df_du = VectorX::Zero(sim->_ndof_u * num_steps);
    sim->_backward_info._df_dvar = VectorX::Ones(sim->_ndof_var * num_steps);
    // sim->_backward_info._df_du = VectorX::Ones(sim->_ndof_u * num_steps);
    // sim->_backward_info._df_dvar = VectorX::Zero(sim->_ndof_var * num_steps);
    sim->backward();
    int t1 = clock();

    std::cerr << "time = " << (t1 - t0) / 1000.0 << "ms" << std::endl;

    std::cerr << "f = " << f << std::endl;
    
    VectorX df_du = sim->_backward_results._df_du;
    VectorX df_daction = df_du;
    for (int i = 0;i < df_daction.size();i++) {
        df_daction[i] *= (1. - tanh(action(i)) * tanh(action(i)));
    }

    std::cerr << "df_du: " << df_du.transpose() << std::endl;
    {
        dtype eps = 1e-2;
        for (int trial = 0;trial < 10;trial++) {
            VectorX df_du_fd = VectorX::Zero(sim->_ndof_u * num_steps);
            for (int i = 0;i < num_steps * sim->_ndof_u;i++) {
                VectorX u_pos = u;
                u_pos(i) += eps;

                sim->set_u(u_pos);
                sim->reset(false, false);
                dtype f_pos = 0.;
                for (int ii = 0;ii < num_steps;ii++) {
                    sim->set_u(u_pos.segment(ii * sim->_ndof_u, (ii + 1) * sim->_ndof_u));
                    sim->forward(1, false, false);
                    VectorX var = sim->get_variables();
                    f_pos += var.sum();
                }
                // f_pos = u_pos.sum();
                df_du_fd(i) = (f_pos - f) / eps;
            }
            print_error_full("eps = " + to_string(eps) + ", df_du", df_du, df_du_fd);
            eps /= 10.;
        }
    }

    std::cerr << "df_daction: " << df_daction.transpose() << std::endl;

    {
        dtype eps = 1e-2;
        for (int trial = 0;trial < 10;trial++) {
            VectorX df_daction_fd = VectorX::Zero(sim->_ndof_u * num_steps);
            for (int i = 0;i < num_steps * sim->_ndof_u;i++) {
                VectorX action_pos = action;
                action_pos(i) += eps;
                VectorX u_pos = action_pos;
                for (int j = 0;j < u_pos.size();j++)
                    u_pos(j) = tanh(u_pos(j));

                sim->set_u(u_pos);
                sim->reset(false, false);
                dtype f_pos = 0.;
                for (int ii = 0;ii < num_steps;ii++) {
                    sim->set_u(u_pos.segment(ii * sim->_ndof_u, (ii + 1) * sim->_ndof_u));
                    sim->forward(1, false, false);
                    VectorX var = sim->get_variables();
                    f_pos += var.sum();
                }
                // f_pos = u_pos.sum();
                df_daction_fd(i) = (f_pos - f) / eps;
            }
            print_error_full("eps = " + to_string(eps) + ", df_daction", df_daction, df_daction_fd);
            eps /= 10.;
        }
    }

    sim->replay();
}

void Test::test_torque_finger() {
    Simulation* sim = SimEnvGenerator::createTorqueFingerDemo("BDF1");

    sim->init();

    int num_steps = 50;
    int ndof_r = sim->_ndof_r;
    int ndof_u = sim->_ndof_u;

    VectorX ctrl_min, ctrl_max;
    sim->get_ctrl_range(ctrl_min, ctrl_max);

    Vector3 q_goal(0., constants::pi / 2., constants::pi / 4.);
    Vector3 P_q = Vector3(10., 1.5, 1.5) * 1e5;

    int t0 = clock();
    sim->reset(false);

    for (int i = 0;i < num_steps;i++) {
        VectorX q = sim->get_q();
        VectorX error = q_goal - q;
        VectorX u = error.cwiseProduct(P_q);
        
        sim->set_u(u);

        sim->forward(1, false);
    }

    int t1 = clock();

    std::cerr << "time = " << (t1 - t0) / 1000.0 << "ms" << std::endl;

    std::cerr << "q = " << sim->get_q().transpose() << ", qdot = " << sim->get_qdot().transpose() << std::endl;

    sim->replay();
}

void Test::test_torque_finger_flick() {
    Simulation* sim = SimEnvGenerator::createTorqueFingerFlickDemo("BDF1");

    sim->init();

    int num_steps = 1000;
    int ndof_r = sim->_ndof_r;
    int ndof_u = sim->_ndof_u;

    VectorX ctrl_min, ctrl_max;
    sim->get_ctrl_range(ctrl_min, ctrl_max);

    Vector3 q_goal(0., 0., 0.);
    Vector3 P_q = Vector3(10., 1.5, 1.5) * 1e5;

    int t0 = clock();
    Vector6 q0 = (Vector6() << 0., constants::pi / 2., constants::pi / 4., 0., 0., 0.).finished();
    sim->set_q_init(q0);

    sim->reset(false);

    for (int i = 0;i < num_steps;i++) {
        VectorX q = sim->get_q();
        VectorX error = q_goal - q.head(3);
        VectorX u = error.cwiseProduct(P_q);
        
        sim->set_u(u);

        sim->forward(1, false);
    }

    int t1 = clock();

    std::cerr << "time = " << (t1 - t0) / 1000.0 << "ms" << std::endl;

    std::cerr << "q = " << sim->get_q().transpose() << ", qdot = " << sim->get_qdot().transpose() << std::endl;

    sim->replay();  
}

void Test::test_pendulum_optimization() {
    Simulation* sim = SimEnvGenerator::createMultiPendulumTest(2, "BDF2");

    sim->init();

    sim->_viewer_options->_loop = false;
    sim->_viewer_options->_infinite = true;
    bool test_gradient = false;

    int num_steps = 100;

    Vector2 q_goal(constants::pi / 2., constants::pi / 4.);

    Vector2 q0(0., 0.);
    Vector2 qdot0(0., 0.);
    for (int opt_iter = 0;opt_iter < 2000;opt_iter ++) {
        sim->set_state_init(q0, qdot0);

        sim->reset(true);
        sim->forward(num_steps, false);

        VectorX q = sim->get_q();
        VectorX qdot = sim->get_qdot();

        if (opt_iter == 0) { // render the initial guess
            sim->replay(); 
        }

        // objective f = |q - q_goal|^2
        dtype f = (q - q_goal).dot(q - q_goal);
        std::cerr << "iter = " << opt_iter << ", f = " << f << std::endl;

        // set backward info
        sim->_backward_info.set_flags(true, true, false, false);

        int ndof_r = sim->_ndof_r;

        // set terminal derivatives
        sim->_backward_info._df_dq0 = VectorX::Zero(ndof_r);
        sim->_backward_info._df_dqdot0 = VectorX::Zero(ndof_r);
        sim->_backward_info._df_dq = VectorX::Zero(ndof_r * num_steps);
        sim->_backward_info._df_dq.tail(ndof_r) = 2. * (q - q_goal);

        // backward
        sim->backward();

        // test gradient
        if (test_gradient) {
            dtype eps = 0.01;
            for (int ii = 0;ii < 10;ii ++) {
                VectorX df_dq0_fd(ndof_r);
                for (int i = 0;i < ndof_r;i++) {
                    VectorX q0_pos = q0;
                    q0_pos(i) += eps;
                    sim->set_state_init(q0_pos, qdot0);
                    sim->reset(false);
                    sim->forward(num_steps, false);
                    VectorX q_pos = sim->get_q();
                    dtype f_pos = (q_pos - q_goal).dot(q_pos - q_goal);

                    VectorX q0_neg = q0;
                    q0_neg(i) -= eps;
                    sim->set_state_init(q0_neg, qdot0);
                    sim->reset(false);
                    sim->forward(num_steps, false);
                    VectorX q_neg = sim->get_q();
                    dtype f_neg = (q_neg - q_goal).dot(q_neg - q_goal);

                    df_dq0_fd(i) = (f_pos - f_neg) / eps / 2.;
                }
                VectorX df_dqdot0_fd(ndof_r);
                for (int i = 0;i < ndof_r;i++) {
                    VectorX qdot0_pos = qdot0;
                    qdot0_pos(i) += eps;
                    sim->set_state_init(q0, qdot0_pos);
                    sim->reset(false);
                    sim->forward(num_steps, false);
                    VectorX q_pos = sim->get_q();
                    dtype f_pos = (q_pos - q_goal).dot(q_pos - q_goal);
                    
                    VectorX qdot0_neg = qdot0;
                    qdot0_neg(i) -= eps;
                    sim->set_state_init(q0, qdot0_neg);
                    sim->reset(false);
                    sim->forward(num_steps, false);
                    VectorX q_neg = sim->get_q();
                    dtype f_neg = (q_neg - q_goal).dot(q_neg - q_goal);

                    df_dqdot0_fd(i) = (f_pos - f_neg) / eps / 2.;
                }

                std::cerr << "eps = " << eps << std::endl;
                print_error("df_dq0", sim->_backward_results._df_dq0, df_dq0_fd);
                print_error("df_dqdot0", sim->_backward_results._df_dqdot0, df_dqdot0_fd);

                eps /= 2.;
            }
        }

        dtype step = 0.02;
        q0 -= step * sim->_backward_results._df_dq0;
        qdot0 -= step * sim->_backward_results._df_dqdot0;
    }

    sim->set_state_init(q0, qdot0);

    sim->reset(false);
    sim->forward(num_steps, false);

    std::cerr << "q0 = " << q0.transpose() << ", qdot0 = " << qdot0.transpose() << std::endl;
    std::cerr << "q = " << sim->get_q().transpose() << ", qdot = " << sim->get_qdot().transpose() << std::endl;

    sim->replay(); // render the optimized values
}

void Test::test_torque_finger_optimization() {
    Simulation* sim = SimEnvGenerator::createTorqueFingerDemo("BDF2");

    sim->init();

    sim->_viewer_options->_loop = true;
    sim->_viewer_options->_infinite = true;
    bool test_gradient = false;

    int num_steps = 50;
    int ndof_r = sim->_ndof_r;
    int ndof_u = sim->_ndof_u;

    VectorX ctrl_min, ctrl_max;
    sim->get_ctrl_range(ctrl_min, ctrl_max);

    Vector3 q_goal(0., constants::pi / 2., constants::pi / 4.);
    
    VectorX u = VectorX::Zero(ndof_u * num_steps);

    int t0 = clock();
    for (int opt_iter = 0;opt_iter < 1000;opt_iter ++) {
        sim->reset(true);
        std::vector<VectorX> qs;
        qs.clear();
        for (int i = 0;i < num_steps;i++) {
            sim->set_u(u.segment(ndof_u * i, ndof_u));
            sim->forward(1, false);
            qs.push_back(sim->get_q());
        }

        if (opt_iter == 0) { // render the initial guess
            sim->replay(); 
        }

        // objective f = sum(|q(i) - q_goal|^2)
        dtype f = 0.;
        for (int i = 0;i < num_steps;i++) {
            f += (qs[i] - q_goal).dot(qs[i] - q_goal);
        }
        std::cerr << "iter = " << opt_iter << ", f = " << f << std::endl;

        // set backward info
        sim->_backward_info.set_flags(false, false, false, true);

        // set terminal derivatives
        sim->_backward_info._df_du = VectorX::Zero(ndof_u * num_steps);
        sim->_backward_info._df_dq = VectorX::Zero(ndof_r * num_steps);
        for (int i = 0;i < num_steps;i++) {
            sim->_backward_info._df_dq.segment(ndof_r * i, ndof_r) = 2. * (qs[i] - q_goal);
        }

        // backward
        sim->backward();

        // test gradient
        if (test_gradient) {
            // dtype eps = 1e-7;
            dtype eps = 0.01;
            for (int ii = 0;ii < 10;ii++) {
                VectorX df_du_fd(ndof_u * num_steps);
                for (int k = 0;k < ndof_u * num_steps;k++) {
                    VectorX u_pos = u;
                    u_pos(k) += eps;
                    sim->reset(false);
                    dtype f_pos = 0.;
                    for (int i = 0;i < num_steps;i++) {
                        sim->set_u(u_pos.segment(ndof_u * i, ndof_u));
                        sim->forward(1, false);
                        VectorX q = sim->get_q();
                        f_pos += (q - q_goal).dot(q - q_goal);
                    }
                    df_du_fd(k) = (f_pos - f) / eps;
                }

                std::cerr << "eps = " << eps << std::endl;
                print_error("df_du", sim->_backward_results._df_du, df_du_fd);

                eps /= 2.;
            }
        }

        dtype step = 0.05;
        u -= step * sim->_backward_results._df_du;
        // for (int i = 0;i < num_steps;i++)
        //     for (int j = 0;j < ndof_u;j++)
        //         u(i * ndof_u + j) = min(ctrl_max(j), max(ctrl_min(j), u(i * ndof_u + j)));
    }

    int t1 = clock();
    std::cerr << "time = " << (t1 - t0) / 1000.0 << "ms" << std::endl;

    sim->reset(false);
    for (int i = 0;i < num_steps;i++) {
        sim->set_u(u.segment(ndof_u * i, ndof_u));
        sim->forward(1, false);
    }

    sim->replay();
}

void Test::test_torque_finger_flick_optimization() {
    Simulation* sim = SimEnvGenerator::createTorqueFingerFlickDemo("BDF2");

    sim->init();

    sim->_viewer_options->_loop = true;
    sim->_viewer_options->_infinite = true;
    bool test_gradient = false;

    int num_steps = 1000;
    int ndof_r = sim->_ndof_r;
    int ndof_u = sim->_ndof_u;

    VectorX ctrl_min, ctrl_max;
    sim->get_ctrl_range(ctrl_min, ctrl_max);

    Vector6 q0 = (Vector6() << 0., constants::pi / 2., constants::pi / 4., 0., 0., 0.).finished();
    sim->set_q_init(q0);

    dtype x_goal = 5.5;
    
    VectorX u = VectorX::Zero(ndof_u * num_steps);

    // initial guess
    Vector3 q_goal(0., 0., 0.);
    Vector3 P_q = Vector3(10., 2.0, 3.0);

    sim->reset(false);

    for (int i = 0;i < num_steps;i++) {
        VectorX q = sim->get_q();
        VectorX error = q_goal - q.head(3);
        VectorX ui = error.cwiseProduct(P_q);
        // ui.setZero();
        // for (int j = 0;j < ui.size();j++)
        //     ui[j] = min(1., max(ui[j], -1.));

        u.segment(ndof_u * i, ndof_u) = ui;

        sim->set_u(ui);

        sim->forward(1, false);
    }

    std::cerr << "u = " << u.head(10).transpose() << std::endl;
    std::cerr << u.maxCoeff() << ", " << u.minCoeff() << std::endl;
    std::cerr << "u.norm = " << u.norm() << std::endl;

    std::cerr << "q = " << sim->get_q().transpose() << ", qdot = " << sim->get_qdot().transpose() << std::endl;
    int t0 = clock();
    for (int opt_iter = 0;opt_iter < 1;opt_iter ++) {
        sim->reset(true);
        for (int i = 0;i < num_steps;i++) {
            sim->set_u(u.segment(ndof_u * i, ndof_u));
            sim->forward(1, false);
        }

        VectorX q = sim->get_q();

        if (opt_iter == 0) { // render the initial guess
            sim->replay(); 
        }

        // objective f = sum(|q(3) - x_goal|^2)
        dtype f = (q(3) - x_goal) * (q(3) - x_goal);
        
        std::cerr << "iter = " << opt_iter << ", f = " << f << std::endl;

        // set backward info
        sim->_backward_info.set_flags(false, false, false, true);

        // set terminal derivatives
        sim->_backward_info._df_du = VectorX::Zero(ndof_u * num_steps);
        sim->_backward_info._df_dq = VectorX::Zero(ndof_r * num_steps);
        sim->_backward_info._df_dq(ndof_r * (num_steps - 1) + 3) = 2. * (q(3) - x_goal);

        // backward
        sim->backward();

        // // test gradient
        // if (test_gradient) {
        //     dtype eps = 1e-7;
        //     VectorX df_du_fd(ndof_u * num_steps);
        //     for (int k = 0;k < ndof_u * num_steps;k++) {
        //         VectorX u_pos = u;
        //         u_pos(k) += eps;
        //         sim->reset(false);
        //         dtype f_pos = 0.;
        //         for (int i = 0;i < num_steps;i++) {
        //             sim->set_u(u_pos.segment(ndof_u * i, ndof_u));
        //             sim->forward(1, false);
        //             VectorX q_pos = sim->get_q();
        //             f_pos += (q_pos(3) - x_goal) * (q_pos(3) - x_goal);
        //         }
        //         df_du_fd(k) = (f_pos - f) / eps;
        //     }

        //     // print_error("df_du", sim->_backward_results._df_du, df_du_fd);
        //     std::cerr << "df_du_fd = " << df_du_fd.head(10).transpose() << std::endl;    
        // }

        std::cerr << "df_du = " << sim->_backward_results._df_du.head(10).transpose() << std::endl;
        std::cerr << "df_du = " << sim->_backward_results._df_du.norm() << std::endl;
        dtype step = 0.2;
        u -= step * sim->_backward_results._df_du;
        // for (int i = 0;i < num_steps;i++)
        //     for (int j = 0;j < ndof_u;j++)
        //         u(i * ndof_u + j) = min(1., max(-1., u(i * ndof_u + j)));
    }

    int t1 = clock();
    std::cerr << "time = " << (t1 - t0) / 1000.0 << "ms" << std::endl;

    sim->reset(false);
    for (int i = 0;i < num_steps;i++) {
        sim->set_u(u.segment(ndof_u * i, ndof_u));
        sim->forward(1, false);
    }

    std::cerr << "q = " << sim->get_q().transpose() << ", qdot = " << sim->get_qdot().transpose() << std::endl;

    sim->replay();
}

void Test::test_ground_contact_optimization() {
    Simulation* sim = SimEnvGenerator::createGroundContactTest("BDF2");

    // target x,y position
    Vector2 pos_goal(1.0, 0.5);
    
    // Add virtual goal
    JointFixed* joint = new JointFixed(sim, 10, nullptr, Matrix3::Identity(), Vector3(pos_goal(0), pos_goal(1), 0.));
    BodyCuboid* target = new BodyCuboid(sim, joint, Vector3(0.2, 0.2, 0.2), Matrix3::Identity(), Vector3::Zero(), 1.);
    target->set_color(Vector3(0., 1., 0.));

    sim->_robot->_root_joints.push_back(joint);
    sim->_robot->init();

    // start optimization
    sim->init();

    sim->_viewer_options->_loop = true;
    sim->_viewer_options->_infinite = true;
    bool test_gradient = false;

    int num_steps = 1000;
    int ndof_r = sim->_ndof_r;

    // optimize the initial state for cube 1 (i.e. qdot[0]~qdot[2])
    Vector3 q0 = (Vector3() << -2, 2, 0).finished();
    Vector3 qdot0 = (Vector3() << 5, 50, 50).finished();

    int t0 = clock();
    for (int opt_iter = 0;opt_iter < 200;opt_iter ++) {
        sim->set_state_init(q0, qdot0);

        sim->reset(true);
        sim->forward(num_steps, false);
        
        VectorX q = sim->get_q();

        if (opt_iter == 0) {// render the initial guess
            sim->replay();
        }

        // objective f = |q[0:2] - pos_goal|^2
        dtype f = (q.head(2) - pos_goal).dot(q.head(2) - pos_goal);
        std::cerr << "iter = " << opt_iter << ", f = " << f << std::endl;
        
        // set backward info
        sim->_backward_info.set_flags(false, true, false, false);

        int ndof_r = sim->_ndof_r;

        // set terminal derivatives
        sim->_backward_info._df_dqdot0 = VectorX::Zero(ndof_r);
        sim->_backward_info._df_dq = VectorX::Zero(ndof_r * num_steps);
        sim->_backward_info._df_dq.segment(ndof_r * (num_steps - 1), 2) = 2. * (q.head(2) - pos_goal);

        // backward
        sim->backward();

        // test gradient
        if (test_gradient) {
            dtype eps = 1e-7;
            VectorX df_dqdot0_fd(ndof_r);
            for (int i = 0;i < ndof_r;i++) {
                VectorX qdot0_pos = qdot0;
                qdot0_pos(i) += eps;
                sim->set_state_init(q0, qdot0_pos);
                sim->reset(false);
                sim->forward(num_steps, false);
                VectorX q_pos = sim->get_q();
                dtype f_pos = (q_pos.head(2) - pos_goal).dot(q_pos.head(2) - pos_goal);
                df_dqdot0_fd(i) = (f_pos - f) / eps;
            }
            
            print_error("df_dqdot0", sim->_backward_results._df_dqdot0, df_dqdot0_fd);
        }

        dtype step = 0.002;
        qdot0.head(3) -= step * sim->_backward_results._df_dqdot0.head(3);

        qdot0(0) = min((dtype)50., max((dtype)-50., qdot0(0)));
        qdot0(1) = min((dtype)50., max((dtype)-50., qdot0(1)));
        qdot0(2) = min((dtype)50., max((dtype)-50., qdot0(2)));
        
        std::cerr << "df_dqdot0 = " << sim->_backward_results._df_dqdot0.transpose() << std::endl;
        std::cerr << "qdot0 = " << qdot0.head(3).transpose() << std::endl;
    }

    int t1 = clock();
    std::cerr << "time = " << (t1 - t0) / 1000.0 << "ms" << std::endl;

    sim->set_state_init(q0, qdot0);

    sim->reset(false);
    sim->forward(num_steps, false);

    std::cerr << "qdot0 = " << qdot0.transpose() << std::endl;
    std::cerr << "q = " << sim->get_q().transpose() << std::endl;

    sim->replay();
}

void Test::test_box_contact_optimization() {
    Simulation* sim = SimEnvGenerator::createBoxContactTest("BDF2");

    // target x position for the cube 2 (i.e. target q[3])
    dtype x_goal = 3.5;

    // Add virtual goal
    JointFixed* joint = new JointFixed(sim, 10, nullptr, Matrix3::Identity(), Vector3(x_goal + 3., 0.5, 0.));
    BodyCuboid* target = new BodyCuboid(sim, joint, Vector3(0.2, 0.2, 0.2), Matrix3::Identity(), Vector3::Zero(), 1.);
    target->set_color(Vector3(0., 1., 0.));

    sim->_robot->_root_joints.push_back(joint);
    sim->_robot->init();

    sim->init();

    sim->_viewer_options->_loop = false;
    sim->_viewer_options->_infinite = true;
    bool test_gradient = false;

    int num_steps = 400;
    int ndof_r = sim->_ndof_r;

    // optimize the initial state for cube 1 (i.e. qdot[0]~qdot[2])
    Vector6 q0 = (Vector6() << 0., 0., 0., 0., 0., 0.).finished();
    Vector6 qdot0 = (Vector6() << 35., 0., 0., 0., 0., 0.).finished();
    for (int opt_iter = 0;opt_iter < 200;opt_iter ++) {
        sim->set_state_init(q0, qdot0);

        sim->reset(true);
        sim->forward(num_steps, false);
        
        VectorX q = sim->get_q();

        // if (opt_iter == 0) {// render the initial guess
        //     std:cerr << "q = " << q.transpose() << std::endl;
        //     sim->replay();
        // }

        // objective f = |q[3] - x_goal|^2
        dtype f = (q[3] - x_goal) * (q[3] - x_goal);
        std::cerr << "iter = " << opt_iter << ", f = " << f << std::endl;
        
        // set backward info
        sim->_backward_info.set_flags(false, true, false, false);

        int ndof_r = sim->_ndof_r;

        // set terminal derivatives
        sim->_backward_info._df_dq0 = VectorX::Zero(ndof_r);
        sim->_backward_info._df_dqdot0 = VectorX::Zero(ndof_r);
        sim->_backward_info._df_dq = VectorX::Zero(ndof_r * num_steps);
        sim->_backward_info._df_dq[ndof_r * (num_steps - 1) + 3] = 2. * (q[3] - x_goal);

        // backward
        sim->backward();

        // test gradient
        if (test_gradient) {
            dtype eps = 1e-7;
            VectorX df_dq0_fd(ndof_r);
            for (int i = 0;i < ndof_r;i++) {
                VectorX q0_pos = q0;
                q0_pos(i) += eps;
                sim->set_state_init(q0_pos, qdot0);
                sim->reset(false);
                sim->forward(num_steps, false);
                VectorX q_pos = sim->get_q();
                dtype f_pos = (q_pos[3] - x_goal) * (q_pos[3] - x_goal);
                df_dq0_fd(i) = (f_pos - f) / eps;
            }
            VectorX df_dqdot0_fd(ndof_r);
            for (int i = 0;i < ndof_r;i++) {
                VectorX qdot0_pos = qdot0;
                qdot0_pos(i) += eps;
                sim->set_state_init(q0, qdot0_pos);
                sim->reset(false);
                sim->forward(num_steps, false);
                VectorX q_pos = sim->get_q();
                dtype f_pos = (q_pos[3] - x_goal) * (q_pos[3] - x_goal);
                df_dqdot0_fd(i) = (f_pos - f) / eps;
            }
            
            // print_error("df_dq0", sim->_backward_results._df_dq0, df_dq0_fd);
            print_error("df_dqdot0", sim->_backward_results._df_dqdot0, df_dqdot0_fd);
        }

        dtype step = 0.5;
        qdot0(0) -= step * sim->_backward_results._df_dqdot0(0);
        qdot0(1) -= step * sim->_backward_results._df_dqdot0(1);

        // qdot0(0) = min(60., max(30., qdot0(0)));
        
        std::cerr << "df_dqdot0 = " << sim->_backward_results._df_dqdot0.head(3).transpose() << std::endl;
        std::cerr << "qdot0 = " << qdot0.head(3).transpose() << std::endl;
    }

    sim->set_state_init(q0, qdot0);

    sim->reset(false);
    sim->forward(num_steps, false);

    std::cerr << "q0 = " << q0.transpose() << ", qdot0 = " << qdot0.transpose() << std::endl;
    std::cerr << "q = " << sim->get_q().transpose() << ", qdot = " << sim->get_qdot().transpose() << std::endl;

    sim->replay();
}

void Test::test_tactile() {
    Simulation* sim = new Simulation("../../assets/tactile_pad.xml");

    sim->_viewer_options->_loop = true;
    sim->_viewer_options->_infinite = true;

    sim->init();

    sim->print_ctrl_info();

    // sim->print_design_params_info();
    
    // auto qdot_init = sim->get_qdot_init();
    // qdot_init[sim->_ndof_r - 2] = 5.;
    // sim->set_qdot_init(qdot_init);

    sim->reset();

    int t0 = clock();
    sim->forward(300, false, false);

    int t1 = clock();

    std::cerr << "time = " << (t1 - t0) / 1000.0 << "ms" << std::endl;

    sim->replay();
}
