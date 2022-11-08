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
#include "static_opengl_shape.h"
#include "geometry.h"

namespace opengl_viewer {

StaticOpenglShape::StaticOpenglShape()
  : OpenglShape(),
  model_matrix_(Eigen::Matrix4f::Identity()),
  vertex_in_world_(Eigen::Matrix3Xf::Zero(3, 0)),
  bounding_box_in_world_(),
  normal_matrix_(Eigen::Matrix3f::Identity()) {}

void StaticOpenglShape::Initialize(const Eigen::Matrix3Xf& vertex_in_model,
  const Eigen::Matrix3Xi& face, const Option& options) {
  OpenglShape::Initialize(vertex_in_model, face, options);

  if (options.HasMatrixOption("model matrix")) {
    model_matrix_ = options.GetMatrixOption("model matrix");
  }
  normal_matrix_ = model_matrix_.topLeftCorner(3, 3).inverse().transpose();
  vertex_in_world_ = HomogeneousTransform(model_matrix_, vertex_in_model);
  bounding_box_in_world_ = model_matrix_ * bounding_box_in_model();
}

}