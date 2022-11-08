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
#version 330 core
layout(triangles) in;
layout(triangle_strip, max_vertices = 18) out;

uniform mat4 shadow_matrix[6];

out vec4 fragment_position;

void main() {
  for (int i = 0; i < 6; ++i) {
    // Built-in variable that specifies to which face we render.
    gl_Layer = i;
    for (int j = 0; j < 3; ++j) {
      fragment_position = gl_in[j].gl_Position;
      gl_Position = shadow_matrix[i] * fragment_position;
      EmitVertex();
    }
    EndPrimitive();
  }
} 