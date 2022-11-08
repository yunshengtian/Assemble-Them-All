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
#version 330

uniform mat4 project_matrix;

in vec2 position;
in vec2 uv;

in vec4 color;
out vec2 fragment_uv;
out vec4 fragment_color;

void main() {
  fragment_uv = uv;
  fragment_color = color;
  gl_Position = project_matrix * vec4(position.xy, 0, 1);
}