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
#include "shader_program.h"
#include <fstream>
#include <iostream>
#include <vector>
#include "GL/glew.h"
#include "geometry.h"
#include "viewer.h"

namespace opengl_viewer {

static const std::string ReadSourceFromFile(const std::string& file_path) {
  // Read source code.
  std::ifstream shader_file(file_path);
  if (!shader_file.is_open()) return "";

  std::string source_code = "";
  std::string line = "";
  while (std::getline(shader_file, line)) {
    source_code += "\n" + line;
  }
  return source_code;
}

static const GLuint InitializeShaderFromSource(const std::string& source_code,
  const int shader_type) {
  // Consider invalid if source_code is empty.
  if (source_code.empty()) return 0;

  // Compile.
  GLuint shader_id = glCreateShader(shader_type);
  const char* source_code_ptr = source_code.c_str();
  glShaderSource(shader_id, 1, &source_code_ptr, NULL);
  glCompileShader(shader_id);

  // Check.
  GLint result = GL_FALSE;
  int info_log_length = 0;
  glGetShaderiv(shader_id, GL_COMPILE_STATUS, &result);
  glGetShaderiv(shader_id, GL_INFO_LOG_LENGTH, &info_log_length);
  if (!result) {
    std::vector<char> error_message(info_log_length + 1);
    glGetShaderInfoLog(shader_id, info_log_length, NULL,
      &error_message[0]);
    std::cout << "Error in " << source_code << ": "
      << std::string(error_message.data()) << std::endl;
    return 0;
  }

  return shader_id;
}

const bool ShaderProgram::InitializeFromFile(
  const std::string& vertex_file_path, const std::string& fragment_file_path,
  const std::string& geometry_file_path) {
  const std::string vertex_source = ReadSourceFromFile(vertex_file_path),
    fragment_source = ReadSourceFromFile(fragment_file_path),
    geometry_source = ReadSourceFromFile(geometry_file_path);
  if (vertex_source.empty() || fragment_source.empty()) return false;

  return InitializeFromSource(vertex_source, fragment_source, geometry_source);
}

const bool ShaderProgram::InitializeFromSource(
  const std::string& vertex_source, const std::string& fragment_source,
  const std::string& geometry_source) {
  const GLuint vertex_id = InitializeShaderFromSource(vertex_source,
    GL_VERTEX_SHADER);
  if (!vertex_id) return false;
  const GLuint fragment_id = InitializeShaderFromSource(fragment_source,
    GL_FRAGMENT_SHADER);
  if (!fragment_id) return false;
  // Optionally, load the geometry shader.
  const bool has_geometry_shader = !geometry_source.empty();
  const GLuint geometry_id = InitializeShaderFromSource(geometry_source,
    GL_GEOMETRY_SHADER);
  if (has_geometry_shader && !geometry_id) return false;

  // Linking.
  GLuint program_id = glCreateProgram();
  glAttachShader(program_id, vertex_id);
  glAttachShader(program_id, fragment_id);
  if (has_geometry_shader) glAttachShader(program_id, geometry_id);
  glLinkProgram(program_id);

  // Check the program.
  GLint result = GL_FALSE;
  glGetProgramiv(program_id, GL_LINK_STATUS, &result);
  int info_log_length = 0;
  glGetProgramiv(program_id, GL_INFO_LOG_LENGTH, &info_log_length);
  if (!result) {
    std::vector<char> error_message(info_log_length + 1);
    glGetProgramInfoLog(program_id, info_log_length, NULL, &error_message[0]);
    std::cout << "Error in linking the shader program:" << std::endl
      << "  " << vertex_source << std::endl
      << "  " << fragment_source << std::endl;
    if (has_geometry_shader)
      std::cout << "  " << geometry_source << std::endl;
    std::cout << "Error message: "
      << std::string(error_message.data()) << std::endl;
    return false;
  }

  glDetachShader(program_id, vertex_id);
  glDeleteShader(vertex_id);
  glDetachShader(program_id, fragment_id);
  glDeleteShader(fragment_id);
  if (has_geometry_shader) {
    glDetachShader(program_id, geometry_id);
    glDeleteShader(geometry_id);
  }

  // Now we are assured program_id is valid, assign it to id_.
  id_ = program_id;
  return true;
}

void ShaderProgram::UseShaderProgram() const {
  glUseProgram(id_);
}

int ShaderProgram::GetUniformLocation(const std::string& name) const {
  return glGetUniformLocation(id_, name.c_str());
}

int ShaderProgram::GetAttribLocation(const std::string& name) const {
  return glGetAttribLocation(id_, name.c_str());
}

void ShaderProgram::SetUniform1i(const std::string& name, const int value)
  const {
  glUniform1i(GetUniformLocation(name), value);
}

void ShaderProgram::SetUniform1f(const std::string& name, const float value)
  const {
  glUniform1f(GetUniformLocation(name), value);
}

void ShaderProgram::SetUniform3f(const std::string& name,
  const Eigen::Vector3f& value) const {
  glUniform3fv(GetUniformLocation(name), 1, value.data());
}

void ShaderProgram::SetUniformMatrix3f(const std::string& name,
  const Eigen::Matrix3f& value) const {
  glUniformMatrix3fv(GetUniformLocation(name), 1, GL_FALSE, value.data());
}

void ShaderProgram::SetUniformMatrix4f(const std::string& name,
  const Eigen::Matrix4f& value) const {
  glUniformMatrix4fv(GetUniformLocation(name), 1, GL_FALSE, value.data());
}

void ShaderProgram::Cleanup() {
  if (id_) {
    glDeleteProgram(id_);
    id_ = 0;
  }
}

const std::string GeneratePhongModelVertexShader(const Viewer& viewer) {
  std::string source = "";
  source += "#version 330 core\n"
    "\n"
    "layout(location = 0) in vec3 vertex;\n"
    "layout(location = 1) in vec3 normal;\n"
    "layout(location = 2) in vec2 object_uv;\n"
    "\n"
    "out vec3 vertex_world;\n"
    "out vec3 normal_world;\n"
    "out vec2 uv;\n"
    "\n"
    "uniform mat4 model_matrix;\n"
    "uniform mat4 view_matrix;\n"
    "uniform mat4 projection_matrix;\n"
    "uniform mat3 normal_matrix;\n"
    "\n"
    "void main() {\n"
    "  gl_Position = projection_matrix * view_matrix * model_matrix\n"
    "  * vec4(vertex, 1.0f);\n"
    "\n"
    "  // Textures.\n"
    "  uv = object_uv;\n"
    "\n"
    "  // Vertex position in world space.\n"
    "  vertex_world = (model_matrix * vec4(vertex, 1)).xyz;\n"
    "\n"
    "  // Normal of the vertex.\n"
    "  normal_world = normal_matrix * normal;\n"
    "}\n";
  return source;
}

const std::string GeneratePhongModelFragShader(const Viewer& viewer) {
  const bool enable_shadow = viewer.options().GetBoolOption("shadow");
  std::string source = "";
  // Version information.
  source += "#version 330 core\n";

  // Inputs, outputs and uniforms.
  source += "\n"
    "in vec3 vertex_world;\n"
    "in vec3 normal_world;\n"
    "in vec2 uv;\n"
    "out vec3 color;\n"
    "uniform mat4 view_matrix;\n";

  // Light information.
  const int light_num = viewer.NumOfPointLights();
  const Option& option = viewer.options();
  const float shadow_sampling_angle =
    DegreeToRadian(option.GetFloatOption("shadow sampling angle"));
  const float shadow_acne_bias = option.GetFloatOption("shadow acne bias");
  const int shadow_sampling_number =
    option.GetIntOption("shadow sampling number");

  const std::string light_num_str = std::to_string(light_num),
    shadow_sampling_angle_str = std::to_string(shadow_sampling_angle),
    shadow_acne_bias_str = std::to_string(shadow_acne_bias),
    shadow_sampling_number_str = std::to_string(shadow_sampling_number);

  source += "\n"
    "struct Light {\n"
    "  vec3 position;\n"
    "  vec3 ambient, diffuse, specular;\n"
    "  float max_dist;\n"
    "};\n"
    "uniform Light point_light[" + light_num_str + "];\n";
  if (enable_shadow)
    source +=
      "uniform samplerCube shadow_map_sampler[" + light_num_str + "];\n";

  // Material.
  source += "\n"
    "struct Material {\n"
    "  vec3 ambient, diffuse, specular;\n"
    "  float shininess;\n"
    "};\n"
    "uniform Material material;\n"
    "uniform sampler2D texture_sampler;\n";

  // Main function.
  source += "\n"
    "void main() {\n"
    "  vec3 n = normalize(normal_world);\n"
    "  mat3 R = transpose(mat3(view_matrix[0].xyz, view_matrix[1].xyz,\n"
    "    view_matrix[2].xyz));\n"
    "  vec3 eye_direction = -R * view_matrix[3].xyz - vertex_world;\n"
    "  vec3 normalized_eye_direction = normalize(eye_direction);\n"
    "\n";

  // Compute the visibility for each light source.
  source += "  float visible[" + light_num_str + "];\n";
  for (int i = 0; i < light_num; ++i) {
    const std::string i_str = std::to_string(i);
    source += "\n  // Light " + i_str + ".\n";
    source += "  {\n"
      "    Light light = point_light[" + i_str + "];\n"
      "    vec3 shadow_dir = vertex_world - light.position;\n"
      "    float depth = length(shadow_dir);\n"
      "    float bias = " + shadow_acne_bias_str + "\n"
      "      * max(1.0 - dot(n, normalize(shadow_dir)), 1.0f);\n"
      "\n";
    // No need to do soft shadow if sampler_num <= 1.
    if (enable_shadow) {
      if (shadow_sampling_number <= 1) {
        source +=
          "    float shadow_depth = texture(shadow_map_sampler["
            + i_str + "], shadow_dir).r * light.max_dist;\n"
          "    if (depth - bias <= shadow_depth) {\n"
          "      visible[" + i_str + "] = 1.0f;\n"
          "    } else {\n"
          "      visible[" + i_str + "] = 0.0f;\n"
          "    }\n";
      } else {
        source +=
          "    float visibility = 0.0f;\n"
          "    float disk_radius = " + shadow_sampling_angle_str + " * depth;\n"
          "    float step = 2.0f * disk_radius / "
            + std::to_string(shadow_sampling_number - 1) + ";\n"
          "    vec3 shadow_origin = shadow_dir - disk_radius * vec3(1, 1, 1);\n"
          "    for (int i = 0; i < " + shadow_sampling_number_str + "; ++i)\n"
          "      for (int j = 0; j < " + shadow_sampling_number_str + "; ++j)\n"
          "        for (int k = 0; k < " + shadow_sampling_number_str + "; ++k) {\n"
          "          vec3 sample_dir = shadow_origin + vec3(i, j, k) * step;\n"
          "          float shadow_depth = texture(shadow_map_sampler["
            + i_str + "], sample_dir).r * light.max_dist;\n"
          "          if (depth - bias <= shadow_depth) {\n"
          "            visibility += 1.0f;\n"
          "          }\n"
          "        }\n"
          "    visibility /= pow(" + shadow_sampling_number_str + ", 3.0f);\n"
          "    visible[" + i_str + "] = visibility;\n";
      }
    } else {
      source += "visible[" + i_str + "] = 1.0;\n";
    }
    source += "  }\n";
  }

  // Compute the color.
  source += "  color = vec3(0.0f, 0.0f, 0.0f);\n";
  source += "  for (int i = 0; i < " + light_num_str + "; ++i) {\n"
    "    float visibility = visible[i];\n"
    "    Light light = point_light[i];\n"
    "    vec3 light_direction = light.position - vertex_world;\n"
    "    float distance = length(light_direction);\n"
    "\n"
    "    // Ambient color.\n"
    "    color += material.ambient * light.ambient;\n"
    "\n"
    "    // Diffuse color.\n"
    "    vec3 normalized_light_direction = normalize(light_direction);\n"
    "    float cos_theta = clamp(dot(n, normalized_light_direction), 0.0f, 1.0f);\n"
    "    color += visibility * cos_theta * material.diffuse * light.diffuse\n"
    "      / (distance * distance);\n"
    "\n"
    "    // Specular color.\n"
    "    vec3 reflect_light_direction = reflect(-normalized_light_direction, n);\n"
    "    float cos_alpha = clamp(dot(n, reflect_light_direction), 0.0f, 1.0f);\n"
    "    color += visibility * pow(cos_alpha, material.shininess)\n"
    "      * material.specular * light.specular / (distance * distance);\n";
  source += "  }\n";  // End of the for loop.

  // Its own texture.
  source += "  color *= texture(texture_sampler, uv).rgb;\n";

  source += "}\n";  // End of main.

  return source;
}

}
