#include "BodyMeshObj.h"
#include "Joint/Joint.h"
#include "tiny_obj_loader.h"
#include "Eigen/Eigenvalues"
#include "Simulation.h"

namespace redmax {

BodyMeshObj::BodyMeshObj(
    Simulation* sim, Joint* joint,
    std::string filename,
    Matrix3 R, Vector3 p,
    TransformType transform_type,
    dtype density,
    Vector3 scale,
    bool adaptive_sample) 
    : Body(sim, joint, density) {
    
    _filename = filename;

    load_mesh(_filename, scale);

    process_mesh();

    if (transform_type == BODY_TO_JOINT) {
        set_transform(R, p);
    } else if (transform_type == OBJ_TO_WOLRD) {
        // E_ji = E_j0 * E_0o * E_oi
        Matrix4 E_0o = math::SE(R, p);
        Matrix4 E_ji = _joint->_E_j0_0 * E_0o * _E_oi;
        set_transform(E_ji.topLeftCorner(3, 3), E_ji.topRightCorner(3, 1));
    } else if (transform_type == OBJ_TO_JOINT) {
        // E_ji = E_jo * E_oi
        Matrix3 R_ji = R * _E_oi.topLeftCorner(3, 3);
        Vector3 p_ji = R * _E_oi.topRightCorner(3, 1) + p;
        set_transform(R_ji, p_ji);
    } else {
        std::cerr << "[BodyMeshObj::BodyMeshObj] Undefined transform_type: " << transform_type << std::endl;
    }

    precompute_bounding_box();

    precompute_contact_points(adaptive_sample);
}

void BodyMeshObj::load_mesh(std::string filename, Vector3 scale) {
    std::vector<tinyobj::shape_t> obj_shape;
    std::vector<tinyobj::material_t> obj_material;
    tinyobj::attrib_t attrib;
    std::string err;
    tinyobj::LoadObj(&attrib, &obj_shape, &obj_material, &err, filename.c_str());

    int num_vertices = (int)attrib.vertices.size() / 3;
    _V.resize(3, num_vertices);
    for (int i = 0;i < num_vertices;i++) {
        _V.col(i) = Vector3(attrib.vertices[i * 3] * scale[0], 
            attrib.vertices[i * 3 + 1] * scale[1],
            attrib.vertices[i * 3 + 2] * scale[2]);
    }
    
    int num_elements = (int)obj_shape[0].mesh.indices.size() / 3;
    _F.resize(3, num_elements);
    for (int i = 0;i < num_elements;i++) {
        _F.col(i) = Vector3i(obj_shape[0].mesh.indices[i * 3].vertex_index,
            obj_shape[0].mesh.indices[i * 3 + 1].vertex_index,
            obj_shape[0].mesh.indices[i * 3 + 2].vertex_index);
    }
}

void BodyMeshObj::process_mesh() {
    // compute mass properties
    dtype volume;
    Vector3 COM;
    Matrix3 I;
    compute_mass_property(_V, _F, volume, COM, I);

    // compute mass
    _mass = volume * _density;

    // computed I assume mass = 1, to need to multiply by mass
    I *= _mass;

    // get the principal axes for inertia tensor by eigenvalue decomposition
    // https://en.wikipedia.org/wiki/Moment_of_inertia#Principal_axes
    _Inertia.setZero();
    Eigen::SelfAdjointEigenSolver<Matrix3> eigensolver(I);
    Vector3 eig_values = eigensolver.eigenvalues();
    Matrix3 eig_vectors = eigensolver.eigenvectors();
    _Inertia.head(3) = eig_values;
    _Inertia(3) = _Inertia(4) = _Inertia(5) = _mass;

    Matrix4 E = Matrix4::Identity();
    E.topLeftCorner(3, 3) = eig_vectors;
    E.topRightCorner(3, 1) = COM;

    // check for right-handedness
    Vector3 x = E.block(0, 0, 3, 1);
    Vector3 y = E.block(0, 1, 3, 1);
    Vector3 z = E.block(0, 2, 3, 1);
    if (x.cross(y).dot(z) < 0.0)
        E.block(0, 2, 3, 1) *= -1;
    
    // check if the rotation part is valid
    Matrix3 res = E.topLeftCorner(3, 3) * E.topLeftCorner(3, 3).transpose();
    if ((res - Matrix3::Identity()).norm() > 1e-6) {
        std::cerr << "invalid rotational part: " << std::endl << E.topLeftCorner(3, 3) << std::endl;
    }
    
    // transform the mesh into body frame
    _E_oi = E;
    _E_io = math::Einv(E);

    Matrix3 R_io = _E_io.topLeftCorner(3, 3);
    Vector3 p_io = _E_io.topRightCorner(3, 1);

    _V = (R_io * _V).colwise() + p_io;

    // std::cerr << "body " << _name << ", I = " << _Inertia.transpose() << std::endl;
}

/**
 * http://melax.github.io/volint.html
 **/
void BodyMeshObj::compute_mass_property(
    const Matrix3X &V, const Matrix3Xi &F, /*input*/
    dtype &volume, Vector3 &COM,          /*output*/
    Matrix3 &I) {
    
    int num_vertices = V.cols();
    int num_faces = F.cols();

    // compute COM and volume
    volume = 0.;
    COM = Vector3::Zero();
    for (int i = 0;i < num_faces;i++) {
        Matrix3 A;
        A.col(0) = V.col(F(0, i));
        A.col(1) = V.col(F(1, i));
        A.col(2) = V.col(F(2, i));
        dtype vol = A.determinant();

        volume += vol;
        COM += vol * (A.col(0) + A.col(1) + A.col(2));
    }
    
    COM /= volume * 4.;
    volume /= 6.;
    
    // compute inertia tensor
    // assume mass = 1.
    Vector3 diag = Vector3::Zero();
    Vector3 offd = Vector3::Zero();
    for (int i = 0;i < num_faces;i++) {
        Matrix3 A;
        A.col(0) = V.col(F(0, i)) - COM;
        A.col(1) = V.col(F(1, i)) - COM;
        A.col(2) = V.col(F(2, i)) - COM;
        A.transposeInPlace();
        dtype d = A.determinant();

        for (int j = 0;j < 3;j++) {
            int j1 = (j + 1) % 3;
            int j2 = (j + 2) % 3;
            diag[j] += (A(0, j) * A(1, j) + A(1, j) * A(2, j) + A(2, j) * A(0, j) + 
                        A(0, j) * A(0, j) + A(1, j) * A(1, j) + A(2, j) * A(2, j)) * d; // divide by 60.0f later;
            offd[j] += (A(0, j1) * A(1, j2) + A(1, j1) * A(2, j2) + A(2, j1) * A(0, j2)  +
                        A(0, j1) * A(2, j2) + A(1, j1) * A(0, j2) + A(2, j1) * A(1, j2)  +
                        A(0, j1) * A(0, j2) * 2 + A(1, j1) * A(1, j2) * 2 + A(2, j1) * A(2, j2) * 2) * d; // divide by 120.0f later
        }
    }

    diag /= volume * 60.;
    offd /= volume * 120.;
    I = (Matrix3() << diag(1) + diag(2), -offd(2), -offd(1),
                        -offd(2), diag(0) + diag(2), -offd(0),
                        -offd(1), -offd(0), diag(0) + diag(1)).finished();
}

void BodyMeshObj::precompute_bounding_box() {
    _bounding_box.first = _V.rowwise().minCoeff();
    _bounding_box.second = _V.rowwise().maxCoeff();
}

void BodyMeshObj::precompute_contact_points(bool adaptive_sample) {
    // simple random sample on the surface vertices
    // TODO: more sophisticated sampling
    srand(1000);
    _contact_points.clear();

    if (adaptive_sample)
        _sample_rate = std::max((int)std::floor(_V.cols() / 1000.0), 1); // adjust sampling rate based on number of vertices
    else
        _sample_rate = 1;

    for (int i = 0;i < _V.cols();i++) {
        int p = rand() % _sample_rate;
        if (p == 0) {
            _contact_points.push_back(_V.col(i));
        }
    }
    
    if (adaptive_sample) {
        // random sample on large faces according to http://extremelearning.com.au/evenly-distributing-points-in-a-triangle/
        // std::cerr << "before: " << _contact_points.size() << std::endl;
        float g = 1.32471795724474602596;
        float alpha1 = 1.0f / g;
        float alpha2 = alpha1 * alpha1;
        float seed = 0.5;
        float min_area = 0.5; // NOTE: need to adjusted case by case (this value is because in our case all objects are within 10x10x10)
        for (int i = 0;i < _F.cols();i++) {
            int p = rand() % _sample_rate;
            if (p == 0) {
                Vector3 A = _V.col(_F(0, i));
                Vector3 B = _V.col(_F(1, i));
                Vector3 C = _V.col(_F(2, i));
                dtype a = (B - C).norm();
                dtype b = (C - A).norm();
                dtype c = (A - B).norm();
                dtype s = 0.5 * (a + b + c);
                dtype area = std::sqrt(s * (s - a) * (s - b) * (s - c));
                if (area > min_area) {
                    int n_sample = (int)std::floor(area / min_area);
                    for (int j = 0;j < n_sample;j++) {
                        float z1 = std::fmod((seed + alpha1 * (j + 1)), 1.0f);
                        float z2 = std::fmod((seed + alpha2 * (j + 1)), 1.0f);
                        if (z1 + z2 > 1.0f) {
                            z1 = 1.0f - z1;
                            z2 = 1.0f - z2;
                        }
                        Vector3 P = (1 - z1 - z2) * A + z1 * B + z2 * C;
                        _contact_points.push_back(P);
                    }
                }
            } 
        }
        // std::cerr << "after: " << _contact_points.size() << std::endl;
    }
}

Matrix3X BodyMeshObj::get_vertices() const {
    return _V;
}

Matrix3Xi BodyMeshObj::get_faces() const {
    return _F;
}

void BodyMeshObj::get_rendering_objects(
    std::vector<Matrix3Xf>& vertex_list, 
    std::vector<Matrix3Xi>& face_list,
    std::vector<opengl_viewer::Option>& option_list,
    std::vector<opengl_viewer::Animator*>& animator_list) {
    
    Matrix3Xf vertex = _V.cast<float>();

    if (_sim->_options->_unit == "cm-g") 
        vertex /= 10.;
    else
        vertex *= 10.;
        
    opengl_viewer::Option object_option;

    object_option.SetBoolOption("smooth normal", false);
    // object_option.SetVectorOption("ambient", 0.25f, 0.148f, 0.06475f);
    object_option.SetVectorOption("ambient", _color(0), _color(1), _color(2));
    object_option.SetVectorOption("diffuse", 0.4f, 0.2368f, 0.1036f);
    object_option.SetVectorOption("specular", 0.774597f, 0.658561f, 0.400621f);
    object_option.SetFloatOption("shininess", 76.8f);

    _animator = new BodyAnimator(this);

    vertex_list.push_back(vertex);
    face_list.push_back(_F);
    option_list.push_back(object_option);
    animator_list.push_back(_animator);
}

bool BodyMeshObj::filter_single(Vector3 xi) {
    if ((xi.array() < _bounding_box.first.array()).any()) {
        return false;
    }
    if ((xi.array() > _bounding_box.second.array()).any()) {
        return false;
    }
    return true;
}

std::vector<int> BodyMeshObj::filter(Matrix3X xi) {
    std::vector<int> filter_indices;
    for (int i = 0;i < xi.cols();++i) {
        if (filter_single(xi.col(i))) {
            filter_indices.push_back(i);
        }
    }
    return filter_indices;
}

}