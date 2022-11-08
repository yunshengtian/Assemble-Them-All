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
#include "static_opengl_light.h"

namespace opengl_viewer {

StaticOpenglLight::StaticOpenglLight()
  : OpenglLight(), position_in_world_(Eigen::Vector3f::Zero()) {}

void StaticOpenglLight::Initialize(const Eigen::Vector3f& position_in_world,
  const Option& options) {
  OpenglLight::Initialize(options);
  position_in_world_ = position_in_world;
}

}