#pragma once
#include "Common.h"
#include "Utils.h"

namespace redmax {

class Joint;
class EndEffector;

class EndEffectorAnimator : public opengl_viewer::Animator {
public:
    EndEffectorAnimator(EndEffector* end_effector) : _end_effector(end_effector) {}

    const Eigen::Matrix4f AnimatedModelMatrix(const float t);

    EndEffector* _end_effector;
};

// EndEffector is defined as local position in joint frame, and the observation value is the position in world frame.
class EndEffector {
public:
    // structure
    Joint* _joint;

    // local position
    Vector3 _pos;

    // index in variables
    int _ndof;                      // number of variable dofs
    std::vector<int> _index;        // the indices in the variable vector

    // rendering
    Vector3f _color;
    dtype _radius;

    EndEffector(Joint* joint, Vector3 pos, dtype radius);

    ~EndEffector();

    // set rendering color
    void set_color(Vector3 color);

    // get variables
    void computeVariables(VectorX& variables);
    void computeVariablesWithDerivative(VectorX& variables, MatrixX& dvar_dq);
    void computeVariablesWithDerivative(VectorX& variables, MatrixX& dvar_dq, MatrixX& dvar_dp);

    // rendering
    void get_rendering_objects(
        std::vector<Matrix3Xf>& vertex_list, 
        std::vector<Matrix3Xi>& face_list,
        std::vector<opengl_viewer::Option>& option_list,
        std::vector<opengl_viewer::Animator*>& animator_list);

protected:
    EndEffectorAnimator* _animator;
};

}