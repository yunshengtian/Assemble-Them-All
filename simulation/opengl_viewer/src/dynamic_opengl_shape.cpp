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
#include "dynamic_opengl_shape.h"
#include "animator.h"

namespace opengl_viewer {

DynamicOpenglShape::DynamicOpenglShape()
  : OpenglShape(), animator_(NULL) {}

DynamicOpenglShape:: ~DynamicOpenglShape() {
  animator_ = NULL;
}

void DynamicOpenglShape::Initialize(const Eigen::Matrix3Xf& vertex_in_model,
  const Eigen::Matrix3Xi& face, Animator* const animator,
  const Option& options) {
  OpenglShape::Initialize(vertex_in_model, face, options);
  animator_ = animator;
}

const Eigen::Matrix4f DynamicOpenglShape::ModelMatrix(const float t) const {
  return animator_->AnimatedModelMatrix(t);
}

}