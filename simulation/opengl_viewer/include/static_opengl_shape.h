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
#ifndef _OPENGL_VIEWER_STATIC_OPENGL_SHAPE_H_
#define _OPENGL_VIEWER_STATIC_OPENGL_SHAPE_H_

#include "opengl_shape.h"

namespace opengl_viewer {

class StaticOpenglShape : public OpenglShape {
public:
  StaticOpenglShape();

  void Initialize(const Eigen::Matrix3Xf& vertex_in_model,
    const Eigen::Matrix3Xi& face, const Option& options = Option());

  const Eigen::Matrix4f ModelMatrix(const float t) const {
    return model_matrix_;
  }
  const Eigen::Matrix3Xf VertexInWorld(const float t) const {
    return vertex_in_world_;
  }
  const Eigen::Matrix3f NormalMatrix(const float t) const {
    return normal_matrix_;
  }
  const BoundingBox BoundingBoxInWorld(const float t) const {
    return bounding_box_in_world_;
  }

private:
  Eigen::Matrix4f model_matrix_;
  // Some helper matrices.
  Eigen::Matrix3Xf vertex_in_world_;
  BoundingBox bounding_box_in_world_;
  Eigen::Matrix3f normal_matrix_;
};

}

#endif