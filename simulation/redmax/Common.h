#pragma once
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <vector>
#include <map>
#include <unordered_map>
#include <iostream>
#include <cmath>
#include <ctime>
#include <fstream>
#include <memory>
#include <chrono>
#include <sstream>
#include <iomanip>
#include <chrono>
#include "opengl_viewer.h"

using namespace std;
using namespace Eigen;

namespace redmax {

using dtype = double;

typedef Eigen::Matrix<dtype, 1, 1> Vector1;
typedef Eigen::Matrix<dtype, 2, 1> Vector2;
typedef Eigen::Matrix<dtype, 3, 1> Vector3;
typedef Eigen::Matrix<dtype, 4, 1> Vector4;
typedef Eigen::Matrix<dtype, 5, 1> Vector5;
typedef Eigen::Matrix<dtype, 6, 1> Vector6;
typedef Eigen::Matrix<dtype, Eigen::Dynamic, 1> VectorX;
typedef Eigen::Matrix<dtype, 1, 1> RowVector1;
typedef Eigen::Matrix<dtype, 1, 2> RowVector2;
typedef Eigen::Matrix<dtype, 1, 3> RowVector3;
typedef Eigen::Matrix<dtype, 1, 4> RowVector4;
typedef Eigen::Matrix<dtype, 1, 5> RowVector5;
typedef Eigen::Matrix<dtype, 1, 6> RowVector6;
typedef Eigen::Matrix<dtype, 1, Eigen::Dynamic> RowVectorX;
typedef Eigen::Matrix<dtype, 3, 3> Matrix3;
typedef Eigen::Matrix<dtype, 4, 4> Matrix4;
typedef Eigen::Matrix<dtype, 5, 5> Matrix5;
typedef Eigen::Matrix<dtype, 6, 6> Matrix6;
typedef Eigen::Matrix<dtype, 3, 6> Matrix36;
typedef Eigen::Matrix<dtype, Eigen::Dynamic, Eigen::Dynamic> MatrixX;
typedef Eigen::Matrix<int  , Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXir;
typedef Eigen::Matrix<dtype, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXr;
typedef Eigen::DiagonalMatrix<dtype, Eigen::Dynamic> DiagonalMatrixX;
typedef Eigen::Matrix<dtype, 3, Eigen::Dynamic> Matrix3X;
typedef Eigen::SparseMatrix<dtype> SparseMatrixX;
typedef Eigen::Triplet<dtype> TripletX;

typedef Matrix3 SO3;
typedef Vector3 so3;
typedef Matrix4 SE3;
typedef Vector6 se3;

inline void throw_error(std::string message) {
    std::cerr << "[Error] " << message << std::endl;
    throw message;
}

}