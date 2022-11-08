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
#include "opengl_light.h"
#include "GL/glew.h"

namespace opengl_viewer {

OpenglLight::OpenglLight()
  : frame_buffer_id_(0), depth_map_id_(0),
  depth_map_size_(1024),
  max_dist_(0.0f),
  ambient_(Eigen::Vector3f::Zero()),
  diffuse_(Eigen::Vector3f::Zero()),
  specular_(Eigen::Vector3f::Zero()) {}

OpenglLight::~OpenglLight() {
  glDeleteFramebuffers(1, &frame_buffer_id_);
  glDeleteTextures(1, &depth_map_id_);
}

void OpenglLight::Initialize(const Option& options) {
  if (options.HasVectorOption("ambient")) {
    ambient_ = options.GetVectorOption("ambient");
  }
  if (options.HasVectorOption("diffuse")) {
    diffuse_ = options.GetVectorOption("diffuse");
  }
  if (options.HasVectorOption("specular")) {
    specular_ = options.GetVectorOption("specular");
  }

  // Initialize the texture.
  if (options.HasIntOption("depth map size")) {
    depth_map_size_ = options.GetIntOption("depth map size");
  }
  glGenFramebuffers(1, &frame_buffer_id_);
  glGenTextures(1, &depth_map_id_);
  glBindTexture(GL_TEXTURE_CUBE_MAP, depth_map_id_);
  // Depth texture.
  for (int j = 0; j < 6; ++j) {
    glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + j, 0, GL_DEPTH_COMPONENT,
      depth_map_size_, depth_map_size_, 0, GL_DEPTH_COMPONENT, GL_FLOAT, 0);
  }
  glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
  // Attach depth texture as FBO's depth buffer.
  glBindFramebuffer(GL_FRAMEBUFFER, frame_buffer_id_);
  glFramebufferTexture(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,
    depth_map_id_, 0);
  glDrawBuffer(GL_NONE);
  glReadBuffer(GL_NONE);
  glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void OpenglLight::BindTexture(const int texture_id) const {
  glActiveTexture(texture_id);
  glBindTexture(GL_TEXTURE_CUBE_MAP, depth_map_id_);
}

void OpenglLight::BindFrameBuffer() const {
  glBindFramebuffer(GL_FRAMEBUFFER, frame_buffer_id_);
}

}