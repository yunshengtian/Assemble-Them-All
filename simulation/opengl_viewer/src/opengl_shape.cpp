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
#include "opengl_shape.h"
#include "GL/glew.h"
#include "geometry.h"

namespace opengl_viewer {

OpenglShape::OpenglShape()
  : vertex_buffer_id_(0), normal_buffer_id_(0), element_buffer_id_(0),
  vertex_in_model_(Eigen::Matrix3Xf::Zero(3, 0)),
  face_(Eigen::Matrix3Xi::Zero(3, 0)),
  bounding_box_in_model_(),
  ambient_(Eigen::Vector3f::Zero()),
  diffuse_(Eigen::Vector3f::Zero()),
  specular_(Eigen::Vector3f::Zero()),
  shininess_(0.0f),
  texture_id_(0), uv_buffer_id_(0) {}

OpenglShape::~OpenglShape() {
  if (vertex_buffer_id_)
    glDeleteBuffers(1, &vertex_buffer_id_);
  if (element_buffer_id_)
    glDeleteBuffers(1, &element_buffer_id_);
  if (normal_buffer_id_)
    glDeleteBuffers(1, &normal_buffer_id_);
  // Determine if we need to delete texture.
  if (texture_id_) {
    glDeleteTextures(1, &texture_id_);
    glDeleteBuffers(1, &uv_buffer_id_);
  }
}

void OpenglShape::Initialize(const Eigen::Matrix3Xf& vertex_in_model,
  const Eigen::Matrix3Xi& face, const Option& options) {
  // Initialize the ambient buffer.
  ambient_ = options.HasVectorOption("ambient")
    ? options.GetVectorOption("ambient") : Eigen::Vector3f::Zero();
  // Initialize the diffuse buffer.
  diffuse_ = options.HasVectorOption("diffuse")
    ? options.GetVectorOption("diffuse") : Eigen::Vector3f::Zero();
  // Initialize the specular buffer.
  specular_ = options.HasVectorOption("specular")
    ? options.GetVectorOption("specular") : Eigen::Vector3f::Zero();
  // Initialize the shininess.
  shininess_ = options.HasFloatOption("shininess")
    ? options.GetFloatOption("shininess") : 1.0f;

  // Initialize geometry.
  bool smooth_normal = true;
  if (options.HasBoolOption("smooth normal")) {
    smooth_normal = options.GetBoolOption("smooth normal");
  }
  // Based on the value of smooth normal, we should determine if we need to
  // flatten vertices and faces or not.
  Eigen::Matrix2Xf uv = options.HasMatrixOption("uv") ?
    options.GetMatrixOption("uv") :
    Eigen::Matrix2Xf::Zero(2, vertex_in_model.cols());
  vertex_in_model_ = vertex_in_model; face_ = face;
  if (!smooth_normal) {
    FlattenMesh(vertex_in_model_, face_, uv);
  }
  bounding_box_in_model_.Initialize(vertex_in_model_);

  // Initialize the vertex buffer.
  glGenBuffers(1, &vertex_buffer_id_);
  glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer_id_);
  const int vertex_num = static_cast<int>(vertex_in_model_.cols());
  glBufferData(GL_ARRAY_BUFFER, 3 * vertex_num * sizeof(float),
    vertex_in_model_.data(), GL_STATIC_DRAW);

  // Initialize the element buffer.
  glGenBuffers(1, &element_buffer_id_);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, element_buffer_id_);
  const int element_num = static_cast<int>(face_.cols()) * 3;
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, element_num * sizeof(int),
    face_.data(), GL_STATIC_DRAW);

  // Initialize the normal buffer.
  glGenBuffers(1, &normal_buffer_id_);
  glBindBuffer(GL_ARRAY_BUFFER, normal_buffer_id_);
  Eigen::Matrix3Xf normal_in_model = Eigen::Matrix3Xf::Zero(3, vertex_num);
  if (smooth_normal) {
    normal_in_model = SmoothVertexNormal(vertex_in_model_, face_);
  } else {
    // Compute the normal.
    const int face_num = static_cast<int>(face_.cols());
    for (int i = 0; i < face_num; ++i) {
      const Eigen::Vector3f v0 = vertex_in_model_.col(3 * i),
        v1 = vertex_in_model_.col(3 * i + 1),
        v2 = vertex_in_model_.col(3 * i + 2);
      const Eigen::Vector3f normal = (v1 - v0).cross(v2 - v1).normalized();
      normal_in_model.middleCols(3 * i, 3) = normal.replicate(1, 3);
    }
  }
  glBufferData(GL_ARRAY_BUFFER, 3 * vertex_num * sizeof(float),
    normal_in_model.data(), GL_STATIC_DRAW);

  // Initialize textures if necessary.
  const bool has_texture_ = options.HasMatrixOption("texture") &&
    options.HasMatrixOption("uv") &&
    options.HasIntOption("texture row num") &&
    options.HasIntOption("texture col num");
  Eigen::Matrix3Xf texture = Eigen::Matrix3Xf::Ones(3, 1);
  int row_num = 1, col_num = 1;
  if (has_texture_) {
    texture = options.GetMatrixOption("texture");
    row_num = options.GetIntOption("texture row num");
    col_num = options.GetIntOption("texture col num");
  }
  // Initialize the texture.
  glGenTextures(1, &texture_id_);
  glBindTexture(GL_TEXTURE_2D, texture_id_);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, col_num, row_num, 0, GL_RGB,
    GL_FLOAT, texture.data());
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  int mag_filter = GL_LINEAR;
  if (options.HasStringOption("texture mag filter")) {
    const std::string mag_filter_str =
      options.GetStringOption("texture mag filter");
    if (mag_filter_str == "nearest") mag_filter = GL_NEAREST;
    // In all other cases we use linear.
  }
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, mag_filter);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
    GL_LINEAR_MIPMAP_LINEAR);
  glGenerateMipmap(GL_TEXTURE_2D);

  // Initialize the uv buffer.
  glGenBuffers(1, &uv_buffer_id_);
  glBindBuffer(GL_ARRAY_BUFFER, uv_buffer_id_);
  glBufferData(GL_ARRAY_BUFFER, 2 * vertex_num * sizeof(float),
    uv.data(), GL_STATIC_DRAW);
}

void OpenglShape::BindVertexBuffer(const int attribute_id) const {
  glEnableVertexAttribArray(attribute_id);
  glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer_id_);
  glVertexAttribPointer(attribute_id,
    3,                  // Size.
    GL_FLOAT,           // Type.
    GL_FALSE,           // Normalized.
    0,                  // Stride.
    NULL                // Array buffer offset.
  );
}

void OpenglShape::Update(const Option& options) {
  // So far we only allow for the color update.
  if (options.HasVectorOption("ambient"))
    ambient_ = options.GetVectorOption("ambient");
  if (options.HasVectorOption("diffuse"))
    diffuse_ = options.GetVectorOption("diffuse");
  if (options.HasVectorOption("specular"))
    specular_ = options.GetVectorOption("specular");
  if (options.HasFloatOption("shininess"))
    shininess_ = options.GetFloatOption("shininess");
  // TODO: allow for more options in the future.
}

void OpenglShape::BindNormalBuffer(const int attribute_id) const {
  glEnableVertexAttribArray(attribute_id);
  glBindBuffer(GL_ARRAY_BUFFER, normal_buffer_id_);
  glVertexAttribPointer(attribute_id, 3, GL_FLOAT, GL_FALSE, 0, NULL);
}

void OpenglShape::BindTexture(const int attribute_id,
  const int texture_id) const {
  glEnableVertexAttribArray(attribute_id);
  glBindBuffer(GL_ARRAY_BUFFER, uv_buffer_id_);
  glVertexAttribPointer(attribute_id, 2, GL_FLOAT, GL_FALSE, 0, NULL);
  // Prepare for the textures.
  glActiveTexture(texture_id);
  glBindTexture(GL_TEXTURE_2D, texture_id_);
}

void OpenglShape::BindElementBuffer() const {
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, element_buffer_id_);
}

void OpenglShape::Draw() const {
  glDrawElements(GL_TRIANGLES, static_cast<GLsizei>(face_.size()),
    GL_UNSIGNED_INT, NULL);
}

const Eigen::Matrix3Xf OpenglShape::VertexInWorld(const float t) const {
  return HomogeneousTransform(ModelMatrix(t), vertex_in_model_);
}

const Eigen::Matrix3f OpenglShape::NormalMatrix(const float t) const {
  const Eigen::Matrix3f model_matrix = ModelMatrix(t).topLeftCorner(3, 3);
  return model_matrix.inverse().transpose();
}

const BoundingBox OpenglShape::BoundingBoxInWorld(const float t) const {
  return ModelMatrix(t) * bounding_box_in_model_;
}

}
