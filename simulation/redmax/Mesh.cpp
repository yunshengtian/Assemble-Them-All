#include "Mesh.h"

namespace redmax {

Mesh::Mesh(Matrix3X V, Matrix3Xi F) {
    _V = V;
    _F = F;

    precompute_bounding_box();
}

void Mesh::precompute_bounding_box() {
    _bounding_box.first = _V.rowwise().minCoeff();
    _bounding_box.second = _V.rowwise().maxCoeff();
}

dtype Mesh::min_distance(Matrix3X X) {
    return distance(X).minCoeff();
}

bool Mesh::filter_single(Vector3 xi) {
    if ((xi.array() < _bounding_box.first.array()).any()) {
        return false;
    }
    if ((xi.array() > _bounding_box.second.array()).any()) {
        return false;
    }
    return true;
}

std::vector<int> Mesh::filter(Matrix3X xi) {
    std::vector<int> filter_indices;
    for (int i = 0;i < xi.cols();++i) {
        if (filter_single(xi.col(i))) {
            filter_indices.push_back(i);
        }
    }
    return filter_indices;
}

SDFMesh::SDFMesh(Matrix3X V, Matrix3Xi F, dtype dx, std::string load_path, std::string save_path): Mesh(V, F) {
    bool loaded = false;
    if (load_path != "") {
        loaded = load_SDF(load_path);
    }
    if (!loaded) {
        precompute_SDF(dx);
        if (save_path != "") {
            save_SDF(save_path);
        }
    }
}

void SDFMesh::precompute_SDF(dtype dx) {

    //start with a massive inside out bound box.
    sdfgen::Vec3f min_box(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max()), 
        max_box(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());

    std::vector<sdfgen::Vec3f> vertList;
    std::vector<sdfgen::Vec3ui> faceList;

    int num_vertices = _V.cols();
    for (int i = 0;i < num_vertices;i++) {
        sdfgen::Vec3f point(_V(0, i), _V(1, i), _V(2, i));
        vertList.push_back(point);
        sdfgen::update_minmax(point, min_box, max_box);
    }

    int num_faces = _F.cols();
    for (int i = 0;i < num_faces;i++) {
        faceList.push_back(sdfgen::Vec3ui(_F(0, i), _F(1, i), _F(2, i)));
    }

    int min_size = 20;
    Vector3 bbox_size = _bounding_box.second - _bounding_box.first;
    _dx = std::min(dx, bbox_size(0) / min_size);
    _dy = std::min(dx, bbox_size(1) / min_size);
    _dz = std::min(dx, bbox_size(2) / min_size);

    // add padding around the box
    sdfgen::Vec3f unit((float)_dx, (float)_dy, (float)_dz);
    int padding = 1;
    min_box -= (float)padding * unit;
    max_box += (float)padding * unit;
    sdfgen::Vec3ui sizes = sdfgen::Vec3ui((max_box - min_box + sdfgen::Vec3f(0.5 * _dx, 0.5 * _dy, 0.5 * _dz)) / sdfgen::Vec3f(_dx, _dy, _dz)) + sdfgen::Vec3ui(1, 1, 1);

    // make SDF
    sdfgen::Array3f SDF;
    sdfgen::make_level_set3(faceList, vertList, min_box, (float)_dx, (float)_dy, (float)_dz, sizes[0], sizes[1], sizes[2], SDF);

    // type conversion
    _min_box = Vector3(min_box[0], min_box[1], min_box[2]);
    _max_box = Vector3(max_box[0], max_box[1], max_box[2]);
    int ni = SDF.ni, nj = SDF.nj, nk = SDF.nk;
    _SDF.resize(ni, nj, nk);
    for (int i = 0; i < ni; ++i) for (int j = 0; j < nj; ++j) for (int k = 0; k < nk; ++k) {
        _SDF(i, j, k) = (dtype)SDF(i, j, k);
    }
}

bool SDFMesh::load_SDF(std::string filename) {
    // std::string inname = filename.substr(0, filename.size()-4) + std::string(".sdf");
    // std::ifstream infile(inname.c_str());
    _SDF_path = filename;
    std::ifstream infile(filename.c_str());
    bool success = infile.good();
    if (success) {
        int ni, nj, nk;
        infile >> ni >> nj >> nk;
        infile >> _min_box[0] >> _min_box[1] >> _min_box[2];
        infile >> _max_box[0] >> _max_box[1] >> _max_box[2];
        infile >> _dx >> _dy >> _dz;
        _SDF.resize(ni, nj, nk);
        for(unsigned int i = 0; i < _SDF.a.size(); ++i) {
            infile >> _SDF.a[i];
        }
    }
    infile.close();
    return success;
}

void SDFMesh::save_SDF(std::string filename) {
    // std::string outname = _filename.substr(0, _filename.size()-4) + std::string(".sdf");    
    // std::ofstream outfile(outname.c_str());
    _SDF_path = filename;
    std::ofstream outfile(filename.c_str());
    outfile << _SDF.ni << " " << _SDF.nj << " " << _SDF.nk << std::endl;
    outfile << _min_box[0] << " " << _min_box[1] << " " << _min_box[2] << std::endl;
    outfile << _max_box[0] << " " << _max_box[1] << " " << _max_box[2] << std::endl;
    outfile << _dx << " " << _dy << " " << _dz << std::endl;
    for(unsigned int i = 0; i < _SDF.a.size(); ++i) {
      outfile << _SDF.a[i] << std::endl;
    }
    outfile.close();
}

void SDFMesh::clear_saved_SDF() {
    // std::string name = _filename.substr(0, _filename.size()-4) + std::string(".sdf");    
    remove(_SDF_path.c_str());
}

dtype SDFMesh::query_SDF(Vector3 x) {

    bool inside = true;
    for (int i = 0;i < 3;i++) {
        inside = inside && x(i) >= _min_box[i] && x(i) <= _max_box[i];
    }

    if (inside) {
        // compute grid location and interpolation coefficients
        dtype ci = (x(0) - _min_box[0]) / _dx, cj = (x(1) - _min_box[1]) / _dy, ck = (x(2) - _min_box[2]) / _dz;
        int i = std::floor(ci), j = std::floor(cj), k = std::floor(ck);
        ci = ci - i, cj = cj - j, ck = ck - k;

        // boundary handling
        int ni = _SDF.ni, nj = _SDF.nj, nk = _SDF.nk;
        int i_pos = (i == ni - 1) ? i : i + 1;
        int j_pos = (j == nj - 1) ? j : j + 1;
        int k_pos = (k == nk - 1) ? k : k + 1;

        // interpolate along i-axis
        dtype vi_j0k0 = (1 - ci) * _SDF(i, j, k) + ci * _SDF(i_pos, j, k);
        dtype vi_j1k0 = (1 - ci) * _SDF(i, j_pos, k) + ci * _SDF(i_pos, j_pos, k);
        dtype vi_j0k1 = (1 - ci) * _SDF(i, j, k_pos) + ci * _SDF(i_pos, j, k_pos);
        dtype vi_j1k1 = (1 - ci) * _SDF(i, j_pos, k_pos) + ci * _SDF(i_pos, j_pos, k_pos);

        // interpolate along j-axis
        dtype vij_k0 = (1 - cj) * vi_j0k0 + cj * vi_j1k0;
        dtype vij_k1 = (1 - cj) * vi_j0k1 + cj * vi_j1k1;

        // interpolate along k-axis
        dtype v = (1 - ck) * vij_k0 + ck * vij_k1;

        return v;

    } else {

        // // find closest boundary point
        // Vector3 x_b;
        // for (int i = 0; i < 3; ++i) x_b(i) = std::min(std::max(x(i), _min_box[i]), _max_box[i]);

        // // approximate by boundary distance + boundary SDF
        // dtype v = (x - x_b).norm() + query_SDF(x_b);
        // return v;

        return (dtype)std::numeric_limits<float>::max();
    }
}

VectorX SDFMesh::distance(Matrix3X X) {
    VectorX d;
    d.resize(X.cols(), 1);
    for (int i = 0; i < X.cols(); ++i) {
        d(i) = query_SDF(X.col(i));
    }
    return d;
}

BVHMesh::BVHMesh(Matrix3X V, Matrix3Xi F): Mesh(V, F) {
    _BVH = BVHEngine::create(3);
    _BVH->set_mesh(_V.transpose(), _F.transpose());
    _BVH->build();
}

VectorX BVHMesh::distance(Matrix3X X) {
    VectorX d = VectorX::Constant(X.cols(), 10.);;

    std::vector<int> filtered_indices = filter(X);

    if (filtered_indices.size() > 0) {
        Matrix3X x_filtered(3, filtered_indices.size());
        for (int i = 0;i < filtered_indices.size();++i)
            x_filtered.col(i) = X.col(filtered_indices[i]);
        
        VectorX d_filtered;
        VectorXi closest_faces;
        MatrixXr closest_points;
        MatrixXr closest_face_normals;
        _BVH->lookup_signed(x_filtered.transpose(), d_filtered, closest_faces, closest_points, closest_face_normals);

        for (int i = 0; i < d_filtered.size(); ++i) {
            d(filtered_indices[i]) = d_filtered(i);
        }
    }

    return d;
}

}