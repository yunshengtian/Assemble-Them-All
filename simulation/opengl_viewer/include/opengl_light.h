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
#ifndef _OPENGL_VIEWER_OPENGL_LIGHT_H_
#define _OPENGL_VIEWER_OPENGL_LIGHT_H_

#include "Eigen/Dense"
#include "option.h"

namespace opengl_viewer {

class OpenglLight {
public:
  OpenglLight();
  virtual ~OpenglLight();

  void Initialize(const Option& options = Option());

  const int depth_map_size() const { return depth_map_size_; }
  const float max_dist() const { return max_dist_; }
  void set_max_dist(const float max_dist) { max_dist_ = max_dist; }
  const Eigen::Vector3f ambient() const { return ambient_; }
  const Eigen::Vector3f diffuse() const { return diffuse_; }
  const Eigen::Vector3f specular() const { return specular_; }

  void BindTexture(const int texture_id) const;
  void BindFrameBuffer() const;

  virtual const Eigen::Vector3f PositionInWorld(const float t) const = 0;

private:
  // For its depth map.
  unsigned int frame_buffer_id_, depth_map_id_;
  int depth_map_size_;
  float max_dist_;

  Eigen::Vector3f ambient_, diffuse_, specular_;
};


}

#endif