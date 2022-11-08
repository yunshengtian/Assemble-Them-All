#include "EndEffector/EndEffector.h"
#include "Utils.h"
#include "Joint/Joint.h"
#include "Simulation.h"

namespace redmax {

const Eigen::Matrix4f EndEffectorAnimator::AnimatedModelMatrix(const float t) {
    Matrix4f model_matrix = _end_effector->_joint->_E_0j.cast<float>();;
    Vector3f pos_world = (_end_effector->_joint->position_in_world(_end_effector->_pos)).cast<float>();;
    
    model_matrix.topRightCorner(3, 1) = pos_world / 10.; // scale for better visualization

    return model_matrix;
}

EndEffector::EndEffector(Joint *joint, Vector3 pos, dtype radius): 
    _joint(joint), _pos(pos), _radius(radius) {
    _ndof = 3;
    _color << 1.f, 0.f, 0.f;
}

EndEffector::~EndEffector() {}

void EndEffector::set_color(Vector3 color) {
    _color = color.cast<float>();
}

void EndEffector::computeVariables(VectorX& variables) {
    Vector3 pos = _joint->position_in_world(_pos);
    
    for (int i = 0;i < 3;i++)
        variables(_index[i]) = pos(i);
}

void EndEffector::computeVariablesWithDerivative(VectorX& variables, MatrixX& dvar_dq) {
    Vector3 pos = _joint->position_in_world(_pos);
    
    for (int i = 0;i < 3;i++)
        variables(_index[i]) = pos(i);
    
    // compute dpos_dq = dE0j_dq * (p, 1)' 
    //                 = (dE0p_dq * Epj + E0p * dEpj_dq) * (p, 1)'
    //                 = dE0p_dq * Epj * (p, 1)' + E0p * dEpj_dq * (p, 1)'
    // tmp = ... * E_(p-1)p * Epj * (p, 1)'
    Vector4 tmp = (Vector4() << _pos, 1).finished();
    for (auto now = _joint;now != nullptr;now = now->_parent) {
        for (int i = 0;i < now->_index.size();i++) {
            if (now->_parent != nullptr)
                dvar_dq.block(_index[0], now->_index[i], 3, 1) += (now->_parent->_E_0j * (now->_dEpj_dq(i) * tmp)).head(3);
            else {
                dvar_dq.block(_index[0], now->_index[i], 3, 1) += (now->_dEpj_dq(i) * tmp).head(3);
            }
        }
        tmp = now->_E_pj * tmp;
    }
}

void EndEffector::computeVariablesWithDerivative(VectorX& variables, MatrixX& dvar_dq, MatrixX& dvar_dp) {
    Vector3 pos = _joint->position_in_world(_pos);
    
    for (int i = 0;i < 3;i++)
        variables(_index[i]) = pos(i);
    
    // compute dpos_dq = dE0j_dq * (p, 1)' 
    //                 = (dE0p_dq * Epj + E0p * dEpj_dq) * (p, 1)'
    //                 = dE0p_dq * Epj * (p, 1)' + E0p * dEpj_dq * (p, 1)'
    // tmp = ... * E_(p-1)p * Epj * (p, 1)'
    Vector4 tmp = (Vector4() << _pos, 1).finished();
    for (auto now = _joint;now != nullptr;now = now->_parent) {
        for (int i = 0;i < now->_index.size();i++) {
            if (now->_parent != nullptr)
                dvar_dq.block(_index[0], now->_index[i], 3, 1) += (now->_parent->_E_0j * (now->_dEpj_dq(i) * tmp)).head(3);
            else {
                dvar_dq.block(_index[0], now->_index[i], 3, 1) += (now->_dEpj_dq(i) * tmp).head(3);
            }
        }
        tmp = now->_E_pj * tmp;
    }

    // design derivatives
    for (auto ancestor = _joint; ancestor != nullptr; ancestor = ancestor->_parent) {
        if (ancestor->_design_params_1._active) {
            for (int k = 0;k < ancestor->_design_params_1._ndof;k++) {
                int idx = ancestor->_design_params_1._param_index(k);
                dvar_dp.block(_index[0], idx, 3, 1) = _joint->_dE0j_dp1(idx).topLeftCorner(3, 3) * _pos + _joint->_dE0j_dp1(idx).topRightCorner(3, 1);
            }
            
        } 
    }
}

// rendering
void EndEffector::get_rendering_objects(
    std::vector<Matrix3Xf>& vertex_list, 
    std::vector<Matrix3Xi>& face_list,
    std::vector<opengl_viewer::Option>& option_list,
    std::vector<opengl_viewer::Animator*>& animator_list) {

    Matrix3Xf vertex;
    Matrix3Xi face;
    Matrix2Xf uv;
    opengl_viewer::Option object_option;

    opengl_viewer::ReadFromObjFile(
        std::string(GRAPHICS_CODEBASE_SOURCE_DIR) + "/resources/meshes/sphere.obj",
        vertex, face, uv);

    for (int i = 0;i < 3;i++)
        vertex.row(i) *= (float)_radius / 10.;

    object_option.SetBoolOption("smooth normal", false);
    object_option.SetVectorOption("ambient", _color(0), _color(1), _color(2));
    object_option.SetVectorOption("diffuse", 0.4f, 0.2368f, 0.1036f);
    object_option.SetVectorOption("specular", 0.774597f, 0.458561f, 0.200621f);
    object_option.SetFloatOption("shininess", 76.8f);

    _animator = new EndEffectorAnimator(this);

    vertex_list.push_back(vertex);
    face_list.push_back(face);
    option_list.push_back(object_option);
    animator_list.push_back(_animator);
}

}