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
#ifndef _OPENGL_VIEWER_SHADER_PROGRAM_H_
#define _OPENGL_VIEWER_SHADER_PROGRAM_H_

#include <string>
#include "Eigen/Dense"

namespace opengl_viewer {

class Viewer;

class ShaderProgram {
public:
  ShaderProgram() : id_(0) {}

  // Load a vertex/fragment/geometry shader and returns if it is successful.
  // The geometry shader is optional.
  const bool InitializeFromFile(const std::string& vertex_file_path,
    const std::string& fragment_file_path,
    const std::string& geometry_file_path = "");
  const bool InitializeFromSource(const std::string& vertex_source,
    const std::string& fragment_source,
    const std::string& geometry_source = "");

  // Wrap up the call to glUseProgram.
  void UseShaderProgram() const;

  // Set uniforms.
  int GetUniformLocation(const std::string& name) const;
  int GetAttribLocation(const std::string& name) const;
  void SetUniform1i(const std::string& name, const int value) const;
  void SetUniform1f(const std::string& name, const float value) const;
  void SetUniform3f(const std::string& name,
    const Eigen::Vector3f& value) const;
  void SetUniformMatrix3f(const std::string& name,
    const Eigen::Matrix3f& value) const;
  void SetUniformMatrix4f(const std::string& name,
    const Eigen::Matrix4f& value) const;

  void Cleanup();

private:
  unsigned int id_;
};

const std::string GeneratePhongModelVertexShader(const Viewer& viewer);
const std::string GeneratePhongModelFragShader(const Viewer& viewer);

}

#endif