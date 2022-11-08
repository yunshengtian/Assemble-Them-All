#pragma once

#include "Common.h"
#include "array3.h"
#include "vec.h"
#include "makelevelset3.h"
#include "Body/BVH/BVHEngine.h"

namespace redmax {

class Mesh {
public:
    Mesh(Matrix3X V, Matrix3Xi F);

    Matrix3X _V;
    Matrix3Xi _F;
    std::pair<Vector3, Vector3> _bounding_box;

    void precompute_bounding_box();
    bool filter_single(Vector3 xi);
    std::vector<int> filter(Matrix3X xi);
    
    virtual VectorX distance(Matrix3X X) = 0;
    dtype min_distance(Matrix3X X);


};

class SDFMesh: public Mesh {
public:
    SDFMesh(Matrix3X V, Matrix3Xi F, dtype dx, std::string load_path = "", std::string save_path = "");

    std::string _SDF_path;
    dtype _dx, _dy, _dz;
    Vector3 _min_box, _max_box;
    sdfgen::Array3<dtype> _SDF;

    VectorX distance(Matrix3X X);

    void clear_saved_SDF();

private:
    void precompute_SDF(dtype dx);
    bool load_SDF(std::string path);
    void save_SDF(std::string path);
    dtype query_SDF(Vector3 x);
};

class BVHMesh: public Mesh {
public:
    BVHMesh(Matrix3X V, Matrix3Xi F);
    
    std::shared_ptr<BVHEngine> _BVH;

    VectorX distance(Matrix3X X);
};

}