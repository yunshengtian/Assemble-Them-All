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
#ifndef _OPENGL_VIEWER_DYNAMIC_OPENGL_SHAPE_H_
#define _OPENGL_VIEWER_DYNAMIC_OPENGL_SHAPE_H_

#include "opengl_shape.h"

namespace opengl_viewer {

class Animator;

class DynamicOpenglShape : public OpenglShape {
public:
  DynamicOpenglShape();
  ~DynamicOpenglShape();

  void Initialize(const Eigen::Matrix3Xf& vertex_in_model,
    const Eigen::Matrix3Xi& face, Animator* const animator,
    const Option& options = Option());

  const Eigen::Matrix4f ModelMatrix(const float t) const;

private:
  Animator* animator_;
};

}

#endif