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
//
//
// A C++ OpenGL wrapper.
// Valid options (categorized by type). Default values are listed and other
// possible values are parenthesized if necessary. Float types are used instead
// of double, as I noticed some double OpenGL types (e.g. dvec3) are not well
// supported in my laptop.
// Int:
// - height: 768.
// - shadow sampling number: 2 by default. This is the number of samples along
//   each axis, i.e., by default it will generate 8 samples.
// - width: 1024.
//
// Float:
// - camera aspect ratio: 4.0f / 3.0f.
// - camera field of view (degree): 90.0f.
// - camera pan speed (unit / pixel): 0.01f.
// - camera range max: 100.0f.
// - camera range min: 0.01f.
// - camera rotate speed (degree / pixel): 0.2f.
// - camera zoom speed (unit per mouse wheel): 1.0f.
// - shadow acne bias: 0.005f by default. Larger bias fixes acne but results in
//   blurry shadows.
// - shadow sampling angle (degree): 1.0f by default. This angle defines a cone
//   where the soft shadow sampling will take place.
//
// Bool:
// - shadow: true.
//
// String:
// - window name: viewer.
//
// Vector:
// - background color: (0.0f, 0.0f, 0.0f, 1.0f).
// - camera look at: (0.0f, 0.0f, 0.0f). This point is also the anchor of the
//     default rotation interaction: when you press the middle button in your
//     mouse the camera will be rotated around this point.
// - camera pos: (3.0f, 4.0f, 5.0f).
// - camera up: (0.0f, 1.0f, 0.0f).
//
// Pointer (for your customized data types):
// - imgui wrapper: a pointer to a subclass of ImGuiWrapper.
// - keyboard handler: a point to a subclass of KeyboardHandler.
// - mouse handler: a pointer to a subclass of MouseHandler.
// - timer: a pointer to a subclass of Timer.
//
// Interacting with the camera (similar to Blender):
// - Scroll your mouse wheel to zoom in/out (towards the center of the image).
// - Shift + press the mouse wheel to translate.
// - Press the mouse wheel to rotate.
//
// Order of commands:
// - Viewer& v = Viewer::GetViewer().
// - Create an Option object and set your options.
// - Call v.Initialize(option).
// - Call v.AddStaticOject, AddDynamicObject, AddStaticPointLight, etc.
// - Call v.Run. It will start a loop and terminate when you press ESC or close
//   the window.
// - Call v.Cleanup() to free memory, etc.
#ifndef _OPENGL_VIEWER_VIEWER_H_
#define _OPENGL_VIEWER_VIEWER_H_

// OpenGL header.
#include <unordered_map>
#include <vector>
#include "GL/glew.h"
#include "glfw3.h"
#include "option.h"
#include "shader_program.h"

namespace opengl_viewer {

class Animator;
class ImGuiWrapper;
class KeyboardHandler;
class MouseHandler;
class OpenglLight;
class OpenglShape;
class SamplingAnimator;
class Timer;

class Viewer {
public:
  static Viewer& GetViewer();

  void Initialize(const Option& option);

  // For Python binding only.
  // This sets up default timer, keyboard handlerm, etc.
  void RegisterPyBindingDefaultComponents(const int fps);

  // Use AddStaticObject to add an object that does not change during
  // visualization. For example: the ground of your virtual world.
  // The object should be a triangle mesh defined in its body frame.
  // Valid options include:
  // model matrix (Matrix4f): an identity matrix by default.
  // smooth normal (bool): determine if vertex normal should be averaged from
  //   adjacent faces. True by default.
  // ambient (Vector3f): the ambient coefficient. Default value is (0, 0, 0).
  // diffuse (Vector3f): see above.
  // specular (Vector3f): see above.
  // shininess (float): default value is 1.
  // uv (Matrix2Xf): texture coordinates.
  // texture (Matrix3Xf): each column contains RGBA value for a pixel, which
  //   are organized in a column-major order.
  // texture row num (int):
  // texture col num (int): the row and column number of the texture image.
  // texture mag filter (string): "linear" or "nearest". By default we use
  //   linear.
  const int AddStaticObject(
    // Each column is a vertex in the world frame.
    const Eigen::Matrix3Xf& vertex,
    // Each column is a triangle, and must have consistent ordering.
    const Eigen::Matrix3Xi& face,
    const Option& options = Option()
  );
  // For Python Binding.
  const int AddStaticObject(
    const std::vector<std::vector<float>>& vertex,
    const std::vector<std::vector<int>>& face,
    const Option& options = Option()
  );
  void UpdateStaticObject(const int object_id, const Option& options);
  // Similar to AddStaticObject, but does not support model matrix. Instead the
  // model matrix is extracted from Animator::AnimatedModelMatrix.
  const int AddDynamicObject(
    const Eigen::Matrix3Xf& vertex,
    const Eigen::Matrix3Xi& face,
    Animator* const animator,
    const Option& options = Option()
  );
  // For Python Binding.
  const int AddDynamicObject(
    const std::vector<std::vector<float>>& vertex,
    const std::vector<std::vector<int>>& face,
    SamplingAnimator* const animator,
    const Option& options = Option()
  );
  // Use the id returned by AddStaticObject/AddDynamicObject to remove it.
  void RemoveObject(const int object_id);

  // Add a static point light source in the world.
  // Valid options include:
  // ambient (Vector3f): (0, 0, 0) by default.
  // diffuse (Vector3f): see above.
  // specular (Vector3f): see above.
  // depth map size (int): 1024 by default.
  void AddStaticPointLight(
    const Eigen::Vector3f& position,
    const Option& options = Option()
  );
  // For Python Binding.
  void AddStaticPointLight(
    const std::vector<float>& position,
    const Option& options = Option()
  );
  // Same as above. The position is extracted from the animator.
  void AddDynamicPointLight(
    Animator* const animator,
    const Option& options = Option()
  );
  const int NumOfPointLights() const {
    return static_cast<int>(point_lights_.size());
  }
  void Run();
  void Cleanup();

  // Const reference to data members.
  const Option& options() const { return options_; }

  // User interaction.
  void MouseWheelScrollCallback(const float y_offset);
  void MouseButtonCallback(const int button, const int action);
  void MouseCursorPosCallback(const float x_pos, const float y_pos);
  void KeyCallback(const int key, const int action);
  void CharCallback(const unsigned int ch);

private:
  // Since we use singleton, we disallow all constructors.
  Viewer();
  Viewer(const Viewer&);
  void operator=(const Viewer&);

  void InitializeOptions(const Option& option);
  const std::string ToLowerCase(const std::string& name) const;
  void RenderShadow(const float t);

  // Parameters.
  Option options_;

  GLFWwindow* window_;

  // Camera matrices.
  Eigen::Matrix4f view_matrix_;
  Eigen::Vector3f rotation_anchor_point_;
  Eigen::Matrix4f projection_matrix_;

  // Shader program.
  ShaderProgram phong_shader_, point_light_depth_shader_;
  // Vertex array id.
  GLuint vertex_array_id_;

  enum VertexAttribute { kVertex = 0, kNormal, kTexture };

  std::unordered_map<int, OpenglShape*> objects_;
  int next_object_id_;
  std::vector<OpenglLight*>point_lights_;

  // Mouse button states.
  bool mouse_wheel_pressed_;
  // Cursor states.
  Eigen::Vector2i mouse_cursor_last_position_;
  // Shift state.
  bool shift_pressed_;
  MouseHandler* mouse_handler_;
  KeyboardHandler* keyboard_handler_;

  // Timer.
  Timer* timer_;

  // ImGui member.
  ImGuiWrapper* imgui_wrapper_;

  // Static data member.
  static Viewer viewer_;
};

}

#endif
