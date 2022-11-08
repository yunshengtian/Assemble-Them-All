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
#ifndef _OPENGL_VIEWER_BOUNDING_BOX_H_
#define _OPENGL_VIEWER_BOUNDING_BOX_H_

#include "Eigen/Dense"

namespace opengl_viewer {

class BoundingBox {
public:
  BoundingBox();

  void Initialize(const Eigen::Vector3f& min_corner,
    const Eigen::Vector3f& box_size);
  void Initialize(const Eigen::Matrix3Xf& vertex);

  const Eigen::Vector3f min_corner() const { return min_corner_; }
  const Eigen::Vector3f box_size() const { return box_size_; }

  void Rotate(const float angle, const Eigen::Vector3f& axis);
  void Scale(const Eigen::Vector3f& scale);
  void Scale(const float scale);
  void Translate(const Eigen::Vector3f& translate);
  void Translate(const float x, const float y, const float z);

  const Eigen::Matrix3Xf AllVertices() const;

private:
  Eigen::Vector3f min_corner_;
  Eigen::Vector3f box_size_;
};

const BoundingBox operator*(const Eigen::Matrix4f& transform,
  const BoundingBox& box);

}
 
#endif