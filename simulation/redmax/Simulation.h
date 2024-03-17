#pragma once
#include "Common.h"
#include "Utils.h"
#include "BackwardData.h"
#include "pugixml.hpp"

namespace redmax {

class Robot;
class SimViewer;
class Joint;
class Body;
class Tactile;

class Simulation {
private:
    // state history
    std::vector<VectorX> _q_his, _qdot_his;

    // utils
    VectorX str_to_eigen(std::string str);
    VectorXi str_to_eigen_int(std::string str);
    std::vector<Vector3> parse_contact_points(std::string str);

    // compute M and f matrices
    void computeMatrices(MatrixX& M, VectorX& f); // evaluate force and mass matrix
    void computeMatrices(MatrixX& M, VectorX& f, JacobianMatrixVector& dM_dq, MatrixX& K, MatrixX& D); // evaluate derivatives to q
    void computeMatrices(MatrixX& M, VectorX& f, JacobianMatrixVector& dM_dq, MatrixX& K, MatrixX& D, MatrixX& df_du); // evaluate derivatives to q and to optimization parameters
    void computeMatrices(
        MatrixX& M, VectorX& f,
        JacobianMatrixVector& dM_dq, MatrixX& K, MatrixX& D,
        MatrixX& df_du,
        JacobianMatrixVector& dM_dp, MatrixX& df_dp); // evaluate derivatives to q and to optimization parameters

    // compute variables
    void computeVariablesWithDerivative(VectorX& variables, MatrixX& dvar_dq);
    void computeVariablesWithDerivative(VectorX& variables, MatrixX& dvar_dq, MatrixX& dvar_dp);

    // temporary variables for evaluation functions
    VectorX _q1, _qdot1, _q0, _qdot0, _q_alpha, _qdot_alpha;
    dtype _h;

    // Viewer
    std::shared_ptr<SimViewer> _viewer;
    int _viewer_step;

    void evaluate_g_BDF1(const VectorX& q1, VectorX& g);
    void evaluate_g_with_derivatives_BDF1(const VectorX& q1, VectorX& g, MatrixX& H, bool save_backward_info = false);
    void evaluate_g_SDIRK2b(const VectorX& q1, VectorX& g);
    void evaluate_g_with_derivatives_SDIRK2b(const VectorX& q1, VectorX& g, MatrixX& H, bool save_backward_info = false);
    void evaluate_g_BDF2(const VectorX& q2, VectorX& g);
    void evaluate_g_with_derivatives_BDF2(const VectorX& q2, VectorX& g, MatrixX& H, bool save_backward_info = false);

    // newton's method
    typedef void (Simulation::*Func)(const VectorX&, VectorX&);
    typedef void (Simulation::*Func_With_Derivatives)(const VectorX&, VectorX&, MatrixX&, bool);
    void newton(VectorX& q, Func func, Func_With_Derivatives func_with_derivatives);
    bool _newton_converge = true;

    // integration methods
    void integration_Euler(const VectorX q0, const VectorX qdot0, const dtype h, VectorX& q1, VectorX& qdot1);
    void integration_Midpoint(const VectorX q0, const VectorX qdot0, const dtype h, VectorX& q1, VectorX& qdot1);
    void integration_RK4(const VectorX q0, const VectorX qdot0, const dtype h, VectorX& q1, VectorX& qdot1);
    void integration_BDF1(const VectorX q0, const VectorX qdot0, const dtype h, VectorX& q1, VectorX& qdot1);
    void integration_SDIRK2(const VectorX q0, const VectorX qdot0, const dtype h, VectorX& q1, VectorX& qdot1);
    void integration_BDF2(const VectorX q0, const VectorX qdot0, const VectorX q1, const VectorX qdot1, const dtype h, VectorX& q2, VectorX& qdot2);

    // backward methods
    void backward_BDF1();
    void backward_BDF2();

    // inverse dynamics
    VectorX inverse_dynamics_BDF1(const VectorX q0, const VectorX qdot0, const VectorX q1, bool reset = true);
    VectorX inverse_dynamics_BDF2(const VectorX q0, const VectorX qdot0, const VectorX q1, const VectorX qdot1, const VectorX q2, bool reset = true);
    std::vector<VectorX> inverse_dynamics_his_BDF1(bool forward, bool reset = true);
    std::vector<VectorX> inverse_dynamics_his_BDF2(bool forward, bool reset = true);

    // xml file related
    std::string _asset_folder;
    
public:
    // -------------------- Constants -----------------------
    class Options {
    public:
        Vector3 _gravity;
        dtype _h;
        string _integrator; // [euler, rk4, BDF1, BDF2]
        string _unit; // ["cm-g", "m-kg"]

        Options(Vector3 gravity = -980. * Vector3::UnitZ(), dtype h = 0.02, string integrator = "BDF1", string unit = "cm-g"): 
            _gravity(gravity), _h(h), _integrator(integrator), _unit(unit) {}
    };
    
    Options *_options;

    class ViewerOptions {
    public:
        int _fps;
        dtype _speed;
        Vector4 _bg_color;
        Vector3 _camera_pos;
        Vector3 _camera_up;
        Vector3 _camera_lookat;
        bool _ground;
        Matrix4 _E_g;
        bool _record;
        std::string _record_folder;
        bool _loop;     // whether loop the replay
        bool _infinite; // whether to replay until close the window

        ViewerOptions() {
            _fps = 30;
            _speed = 1.;
            _bg_color = Vector4(0., 0., 0., 1.);
            _camera_pos = Vector3(1., -1.25, 0.75);
            _camera_up = Vector3::UnitZ();
            _camera_lookat = Vector3::Zero();
            _ground = false;
            _E_g.topLeftCorner(3, 3).setIdentity();
            _E_g.topRightCorner(3, 1) = Vector3::UnitZ() * -2.;
            _record = false;
            _loop = true;
            _infinite = true;
        }
    };

    ViewerOptions * _viewer_options;

    class TimeReport {
    public:
        long long _time_solver, _time_save_backward, _time_backward;
        long long _time_compute_matrices, _time_compose_matrices;
        long long _time_compute_dJ, _time_compute_df;
        long long _time_dM_dp, _time_df_dp;
        long long _time_dM_dp1, _time_dM_dp2, _time_dM_dp4;

        void reset() {
            _time_solver = _time_save_backward = _time_backward = 0;
            _time_compute_matrices = _time_compose_matrices = 0;
            _time_compute_dJ = _time_compute_df = 0;
            _time_dM_dp = _time_df_dp = 0;
            _time_dM_dp1 = _time_dM_dp2 = _time_dM_dp4 = 0;
        }
    };

    TimeReport _time_report;

    // -------------------- forward dynamics related -------------------
    std::string _name;

    // robot related
    Robot* _robot;
    std::map<std::string, Joint*> _joint_map;
    std::map<std::string, Body*> _body_map;

    bool _ground;
    Matrix4 _E_g; // ground transform matrix

    int _ndof_r, _ndof_m, _ndof_u, _ndof_var;
    int _ndof_p, _ndof_p1, _ndof_p2, _ndof_p3, _ndof_p4, _ndof_p5, _ndof_p6;

    // controller parameters
    VectorX _phi;

    // states
    VectorX _q_init, _qdot_init;

    // reserved variables for matrices
    JacobianMatrixVector _dM_dp;
    MatrixX _df_dp;

    // backward related
    bool _backward_flag, _backward_design_params_flag; // whether do backward after forward
    BackwardInfo _backward_info;
    BackwardResults _backward_results;

    // verbose output
    bool _verbose;

    // constructors
    Simulation(Options *options, std::string name = "");
    Simulation(Options *options, ViewerOptions *viewer_options, std::string name = "");
    Simulation(std::string xml_file_path, bool verbose = false);
    Simulation(std::string xml_string, std::string asset_folder, bool verbose = false);
    
    // destructor
    ~Simulation();

    Joint* parse_from_xml_file(pugi::xml_node root, pugi::xml_node node, \
                                Joint* parent_joint, int &joint_cnt, bool verbose = false);

    void addRobot(Robot* robot) {
        _robot = robot;
    }

    // init simulation
    void init(bool verbose = false);    

    // init states
    void set_state_init(const VectorX q_init, const VectorX qdot_init);
    void set_q_init(const VectorX q_init);
    void set_qdot_init(const VectorX qdot_init);
    const VectorX get_q_init();
    const VectorX get_qdot_init();

    // states
    void set_state(const VectorX q, const VectorX qdot);
    void set_q(const VectorX q);
    void set_qdot(const VectorX qdot);
    const VectorX get_q();
    const VectorX get_qdot();
    void reparam();

    // history
    const std::vector<VectorX> get_q_his();
    const std::vector<VectorX> get_qdot_his();
    void set_state_his(std::vector<VectorX> q_his, std::vector<VectorX> qdot_his);

    // variables
    const VectorX get_variables();

    // control variables
    void set_u(const VectorX& u);
    void get_ctrl_range(VectorX& ctrl_min, VectorX& ctrl_max);
    void print_ctrl_info();

    // joint related
    Joint* get_joint(const string name);
    VectorX get_joint_q(const string name);
    VectorX get_joint_qdot(const string name);
    VectorX get_joint_qm(const string name);
    VectorX get_joint_qmdot(const string name);
    void set_joint_q(const string name, VectorX q);
    void set_joint_qm(const string name, VectorX qm);
    void zero_joint_q(const string name);
    void zero_joint_qdot(const string name);
    
    // body related
    Body* get_body(const string name);
    dtype get_body_mass(const string name);
    SE3 get_body_E0i(const string name);
    SE3 get_body_Ei0(const string name);
    Matrix3X get_body_vertices(const string name, bool world_frame = false);
    Matrix3Xi get_body_faces(const string name);
    void set_body_external_force(const string name, const VectorX& force);
    dtype get_body_distance(const string name_from, const string name_to);
    bool body_in_contact(const string name_a, const string name_b, dtype eps = 1e-5);
    std::vector<std::string> get_contact_bodies(const string name);
    void clear_contact_bodies(const string name);
    void clear_all_contact_bodies();
    void clear_all_saved_sdfs();

    // tactile data
    std::vector<dtype> get_tactile_depth(const string name);
    std::vector<Vector2i> get_tactile_image_pos(const string name);

    // design parameters
    void set_design_params(const VectorX &design_params);
    VectorX get_design_params();
    void print_design_params_info();

    // set contact coefficient scale
    void set_contact_scale(dtype scale);

    // rendering mesh for abstract bodies
    void set_rendering_mesh_vertices(const std::vector<Matrix3X> &Vs);
    void set_rendering_mesh(const std::vector<Matrix3X> &Vs, const std::vector<Matrix3Xi> &Fs);

    // virtual objects
    void update_virtual_object(std::string name, VectorX data);

    // updata robot to propagate the state
    void update_robot(bool design_gradient = false);

    void test_derivatives_runtime();
    void test_design_derivatives_runtime();

    // reset the simulation
    void reset(bool backward_flag = false, bool backward_design_params_flag = true);
    
    // verbose defines if log state history
    void forward(int num_steps, bool verbose = false, bool test_derivatives = false);

    bool is_converged() const;
    
    // backward to compute the derivatives
    void backward();

    // compute inverse dynamics of the forward trajectory
    std::vector<VectorX> inverse_dynamics_his(bool forward, bool reset = true);

    // rendering
    void get_rendering_objects(
            std::vector<Matrix3Xf>& vertex_list, 
            std::vector<Matrix3Xi>& face_list,
            std::vector<opengl_viewer::Option>& option_list,
            std::vector<opengl_viewer::Animator*>& animator_list);

    void init_viewer();
    bool advance_viewer_step(int num_steps);
    void replay();

    // export simulation replay to a folder
    void export_replay(std::string folder);

    void print_time_report();
};

}