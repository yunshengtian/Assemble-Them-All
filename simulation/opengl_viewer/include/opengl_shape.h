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
#ifndef _OPENGL_VIEWER_OPENGL_SHAPE_H_
#define _OPENGL_VIEWER_OPENGL_SHAPE_H_

#include "Eigen/Dense"
#include "bounding_box.h"
#include "option.h"

namespace opengl_viewer {

class OpenglShape {
public:
  OpenglShape();
  virtual ~OpenglShape();

  void Initialize(const Eigen::Matrix3Xf& vertex_in_model,
    const Eigen::Matrix3Xi& face, const Option& options = Option());
  void Update(const Option& options);

  const Eigen::Matrix3Xf vertex_in_model() const { return vertex_in_model_; }
  const BoundingBox bounding_box_in_model() const {
    return bounding_box_in_model_;
  }
  const Eigen::Matrix3Xi face() const { return face_; }

  const Eigen::Vector3f ambient() const { return ambient_; }
  const Eigen::Vector3f diffuse() const { return diffuse_; }
  const Eigen::Vector3f specular() const { return specular_; }
  const float shininess() const { return shininess_; }

  // Wraps up OpenGL commands.
  void BindVertexBuffer(const int attribute_id) const;
  void BindNormalBuffer(const int attribute_id) const;
  void BindTexture(const int attribute_id, const int texture_id) const;
  void BindElementBuffer() const;
  void Draw() const;

  // Time relevant functions. For static objects simply return constant values.
  virtual const Eigen::Matrix4f ModelMatrix(const float t) const = 0;
  // Declare them virtual so the static object can use the cached value.
  virtual const Eigen::Matrix3Xf VertexInWorld(const float t) const;
  virtual const Eigen::Matrix3f NormalMatrix(const float t) const;
  virtual const BoundingBox BoundingBoxInWorld(const float t) const;

private:
  // Geometry related members.
  unsigned int vertex_buffer_id_, normal_buffer_id_;
  unsigned int element_buffer_id_;

  Eigen::Matrix3Xf vertex_in_model_;
  Eigen::Matrix3Xi face_;
  BoundingBox bounding_box_in_model_;

  // Materials.
  Eigen::Vector3f ambient_, diffuse_, specular_;
  float shininess_;

  // Texture related members.
  unsigned int texture_id_;
  unsigned int uv_buffer_id_;
};

}

#endif
