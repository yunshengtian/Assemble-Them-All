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
#ifndef _OPENGL_VIEWER_IMGUI_WRAPPER_H_
#define _OPENGL_VIEWER_IMGUI_WRAPPER_H_

#include "shader_program.h"

struct GLFWwindow;
struct ImDrawData;

namespace opengl_viewer {

class ImGuiWrapper {
public:
  ImGuiWrapper();
  virtual ~ImGuiWrapper() {}

  // If you want to use ImGui, this is the only function you need to override.
  virtual void SetupUi() {}

  friend class Viewer;

private:
  void Initialize(GLFWwindow* window);
  void NewFrame();
  void Render();
  void Cleanup();

  void CreateDeviceObjects();
  void CreateFontsTexture();
  void RenderDrawLists(ImDrawData* draw_data);

  void MouseWheelScrollCallback(const float y_offset);
  void MouseButtonCallback(const int button, const int action);
  void KeyCallback(const int key, const int action);
  void CharCallback(const unsigned int ch);

  GLFWwindow*  window_;
  ShaderProgram ui_shader_;

  // Used to compute dt.
  float time_stamp_;

  // Mouse states.
  bool mouse_pressed_[3];
  float mouse_wheel_;

  // Members relevant to OpenGL shaders.
  unsigned int vertex_array_id_, vertex_buffer_id_, element_id_;
  unsigned int font_texture_id_;
  int tex_location_, project_matirx_location_;
  int pos_location_, uv_location_, color_location_;
};

} // opengl_viewer

#endif  // _OPENGL_VIEWER_IMGUI_WRAPPER_H_