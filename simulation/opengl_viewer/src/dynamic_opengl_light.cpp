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
#include "dynamic_opengl_light.h"
#include "animator.h"

namespace opengl_viewer {

DynamicOpenglLight::DynamicOpenglLight()
  : OpenglLight(), animator_(NULL) {}

DynamicOpenglLight::~DynamicOpenglLight() {
  animator_ = NULL;
}

void DynamicOpenglLight::Initialize(Animator* const animator,
  const Option& options) {
  OpenglLight::Initialize(options);
  animator_ = animator;
}

const Eigen::Vector3f DynamicOpenglLight::PositionInWorld(
  const float t) const {
  return animator_->AnimatedModelMatrix(t).topRightCorner(3, 1);
}

}