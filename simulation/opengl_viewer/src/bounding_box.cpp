// @ Copyright 2016 Massachusetts Institute of Technology.
// 
// This program is free software; you can redistribute it and / or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
// MA 02110-1301, USA.
#include "bounding_box.h"
#include "geometry.h"

namespace opengl_viewer {

BoundingBox::BoundingBox()
  : min_corner_(Eigen::Vector3f::Zero()),
  box_size_(Eigen::Vector3f::Zero()) {}

void BoundingBox::Initialize(const Eigen::Vector3f& min_corner,
  const Eigen::Vector3f& box_size) {
  min_corner_ = min_corner;
  box_size_ = box_size_;
}

void BoundingBox::Initialize(const Eigen::Matrix3Xf& vertex) {
  if (vertex.cols() == 0) return;
  min_corner_ = vertex.rowwise().minCoeff();
  box_size_ = vertex.rowwise().maxCoeff() - min_corner_;
}

void BoundingBox::Rotate(const float angle, const Eigen::Vector3f& axis) {
  *this = opengl_viewer::Rotate(angle, axis) * (*this);
}

void BoundingBox::Scale(const Eigen::Vector3f& scale) {
  min_corner_ = min_corner_.cwiseProduct(scale);
  box_size_ = box_size_.cwiseProduct(scale);
}

void BoundingBox::Scale(const float scale) {
  this->Scale(Eigen::Vector3f(scale, scale, scale));
}

void BoundingBox::Translate(const Eigen::Vector3f& translate) {
  min_corner_ += translate;
}

void BoundingBox::Translate(const float x, const float y, const float z) {
  this->Translate(Eigen::Vector3f(x, y, z));
}

const Eigen::Matrix3Xf BoundingBox::AllVertices() const {
  const Eigen::Matrix<float, 3, 8> scale = (Eigen::Matrix<float, 3, 8>()
    << 0, 0, 0, 0, 1, 1, 1, 1,
    0, 0, 1, 1, 0, 0, 1, 1,
    0, 1, 0, 1, 0, 1, 0, 1).finished();
  return (box_size_.asDiagonal() * scale).colwise() + min_corner_;
}

const BoundingBox operator*(const Eigen::Matrix4f& transform,
  const BoundingBox& box) {
  BoundingBox new_box;
  new_box.Initialize(HomogeneousTransform(transform, box.AllVertices()));
  return new_box;
}

}