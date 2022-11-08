#pragma once
#include "Common.h"
#include "Simulation.h"
                    
class Test {
public:
    static Eigen::VectorXd random_generated_vector(int len, Eigen::VectorXd lower_bound, Eigen::VectorXd upper_bound);

    static Eigen::VectorXd random_generated_vector(int len, double lower_bound, double upper_bound);

    static Eigen::MatrixXd random_generated_matrix(int n, int m, double lower_bound, double upper_bound);

    static void test_xml_reader();

    static void test_speed();
    
    static void test_forward();

    static void test_derivative();

    static void test_design_derivative();
    
    static void test_control_derivative();
    
    static void test_torque_finger();

    static void test_torque_finger_flick();

    static void test_pendulum_optimization();

    static void test_torque_finger_optimization();

    static void test_torque_finger_flick_optimization();

    static void test_ground_contact_optimization();

    static void test_box_contact_optimization();

    static void test_tactile();
};