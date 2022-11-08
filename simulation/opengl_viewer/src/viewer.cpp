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
#include "viewer.h"
#include "imgui.h"
#include "imgui_wrapper.h"
#include <iostream>
#include <fstream>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include "dynamic_opengl_light.h"
#include "dynamic_opengl_shape.h"
#include "geometry.h"
#include "image.h"
#include "keyboard_handler.h"
#include "mouse_handler.h"
#include "pybinding_default_keyboard_handler.h"
#include "pybinding_default_timer.h"
#include "sampling_animator.h"
#include "static_opengl_light.h"
#include "static_opengl_shape.h"
#include "timer.h"

namespace opengl_viewer {

// Max number of point light.
static const int kMaxPointLightNum = 32;
Viewer Viewer::viewer_ = Viewer();

Viewer& Viewer::GetViewer() {
  return viewer_;
}

static void MouseWheelScrollCallback(GLFWwindow* window, double x_offset,
  double y_offset) {
  Viewer::GetViewer().MouseWheelScrollCallback(static_cast<float>(y_offset));
}

static void MouseButtonCallback(GLFWwindow* window, int button, int action,
  int mods) {
  Viewer::GetViewer().MouseButtonCallback(button, action);
}

static void MouseCursorPosCallback(GLFWwindow* window, double x_pos,
  double y_pos) {
  Viewer::GetViewer().MouseCursorPosCallback(static_cast<float>(x_pos),
    static_cast<float>(y_pos));
}

static void KeyCallback(GLFWwindow* window, int key, int scancode, int action,
  int mods) {
  Viewer::GetViewer().KeyCallback(key, action);
}

static void CharCallback(GLFWwindow* window, unsigned int ch) {
  Viewer::GetViewer().CharCallback(ch);
}

Viewer::Viewer()
  : window_(NULL),
  view_matrix_(Eigen::Matrix4f::Zero()),
  rotation_anchor_point_(Eigen::Vector3f::Zero()),
  projection_matrix_(Eigen::Matrix4f::Zero()),
  vertex_array_id_(0),
  objects_(0),
  next_object_id_(0),
  point_lights_(0),
  mouse_wheel_pressed_(false),
  mouse_cursor_last_position_(0, 0),
  shift_pressed_(false),
  mouse_handler_(NULL),
  keyboard_handler_(NULL),
  timer_(NULL),
  imgui_wrapper_(NULL) {
  // Initialize all parameters with default values.
  // Int parameters.
  options_.SetIntOption("height", 768);
  options_.SetIntOption("shadow sampling number", 2);
  options_.SetIntOption("width", 1024);

  // Float parameters.
  options_.SetFloatOption("camera aspect ratio", 4.0f / 3.0f);
  options_.SetFloatOption("camera field of view", 90.0f);
  options_.SetFloatOption("camera pan speed", 0.01f);
  options_.SetFloatOption("camera range max", 100.0f);
  options_.SetFloatOption("camera range min", 0.01f);
  options_.SetFloatOption("camera rotate speed", 0.2f);
  options_.SetFloatOption("camera zoom speed", 1.0f);
  options_.SetFloatOption("shadow acne bias", 0.005f);
  options_.SetFloatOption("shadow sampling angle", 1.0f);

  // Boolean parameters.
  options_.SetBoolOption("shadow", true);
  options_.SetBoolOption("record", false);

  // String parameters.
  options_.SetStringOption("record folder", "");
  options_.SetStringOption("window name", "viewer");

  // Vector parameters.
  options_.SetVectorOption("background color", 0.0f, 0.0f, 0.0f, 1.0f);
  options_.SetVectorOption("camera look at", 0.0f, 0.0f, 0.0f);
  options_.SetVectorOption("camera pos", 3.0f, 4.0f, 5.0f);
  options_.SetVectorOption("camera up", 0.0f, 1.0f, 0.0f);

  // Pointer parameters.
  options_.SetPointerOption("imgui wrapper", NULL);
  options_.SetPointerOption("keyboard handler", NULL);
  options_.SetPointerOption("mouse handler", NULL);
  options_.SetPointerOption("timer", NULL);
}

void Viewer::RegisterPyBindingDefaultComponents(const int fps) {
  // TODO: fix this memory leak.
  PyBindingDefaultKeyboardHandler* keyboard = new PyBindingDefaultKeyboardHandler();
  PyBindingDefaultTimer* timer = new PyBindingDefaultTimer();
  timer->Initialize(fps, keyboard);

  keyboard_handler_ = dynamic_cast<KeyboardHandler*>(keyboard);
  timer_ = dynamic_cast<Timer*>(timer);
}

void Viewer::Initialize(const Option& option) {
  // Initialize options.
  InitializeOptions(option);

  // Initialize GLFW.
  if (!glfwInit()) {
    std::cout << "Error: Failed to initialize GLFW. Consider updating your "
      "graphics card driver." << std::endl;
    exit(0);
  }
  glfwWindowHint(GLFW_SAMPLES, 4);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

  // Open a window and create its OpenGL context.
  const int height = options_.GetIntOption("height");
  const int width = options_.GetIntOption("width");
  const std::string& name = options_.GetStringOption("window name");
  const bool record = options_.GetBoolOption("record");
  if (record) glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);
  window_ = glfwCreateWindow(width, height, name.c_str(), NULL, NULL);
  if (!window_) {
    std::cout << "Error: Failed to create GLFW window. Consider updating your "
      "graphics card driver." << std::endl;
    glfwTerminate();
    exit(0);
  }
  glfwMakeContextCurrent(window_);

  // Initialize GLEW.
  // Enabling experimental features to resolve seg faults.
  // https://stackoverflow.com/questions/8302625/segmentation-fault-at-glgenvertexarrays-1-vao
  glewExperimental = GL_TRUE;
  if (glewInit() != GLEW_OK) {
    glfwTerminate();
    exit(0);
  }

  // Initialize ImGuiWrapper.
  imgui_wrapper_ = reinterpret_cast<ImGuiWrapper*>(
    options_.GetPointerOption("imgui wrapper"));
  if (imgui_wrapper_) imgui_wrapper_->Initialize(window_);

  // User interaction.
  glfwSetInputMode(window_, GLFW_STICKY_KEYS, GL_FALSE);
  glfwSetInputMode(window_, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
  glfwSetInputMode(window_, GLFW_STICKY_MOUSE_BUTTONS, GL_FALSE);
  glfwSetScrollCallback(window_, opengl_viewer::MouseWheelScrollCallback);
  glfwSetMouseButtonCallback(window_, opengl_viewer::MouseButtonCallback);
  glfwSetCursorPosCallback(window_, opengl_viewer::MouseCursorPosCallback);
  glfwSetKeyCallback(window_, opengl_viewer::KeyCallback);
  glfwSetCharCallback(window_, opengl_viewer::CharCallback);
  // Set the mouse at the center of the screen.
  glfwPollEvents();
  // glfwSetCursorPos(window_, width / 2, height / 2);
  mouse_cursor_last_position_ = Eigen::Vector2i(width / 2, height / 2);

  // Background color.
  const Eigen::Vector4f background_color =
    options_.GetVectorOption("background color");
  glClearColor(background_color(0), background_color(1),
    background_color(2), background_color(3));

  // Enable depth test.
  glEnable(GL_DEPTH_TEST);
  // Accept fragment if it closer to the camera than the former one.
  glDepthFunc(GL_LESS);

  // Initialize VAO.
  glGenVertexArrays(1, &vertex_array_id_);
  glBindVertexArray(vertex_array_id_);

  // Initialize the view matrix.
  const Eigen::Vector3f& camera_pos = options_.GetVectorOption("camera pos");
  rotation_anchor_point_ = options_.GetVectorOption("camera look at");
  const Eigen::Vector3f& camera_up = options_.GetVectorOption("camera up");
  view_matrix_ = LookAt(camera_pos, rotation_anchor_point_, camera_up);

  // Initialize the projection matrix.
  const float aspect_ratio = options_.GetFloatOption("camera aspect ratio");
  const float field_of_view = options_.GetFloatOption("camera field of view");
  const float z_max = options_.GetFloatOption("camera range max");
  const float z_min = options_.GetFloatOption("camera range min");
  projection_matrix_ = Perspective(aspect_ratio, field_of_view, z_min, z_max);

  // Initialize the timer.
  timer_ = reinterpret_cast<Timer*>(options_.GetPointerOption("timer"));

  // Initialize the user interaction object.
  mouse_handler_ = reinterpret_cast<MouseHandler*>(
    options_.GetPointerOption("mouse handler"));
  keyboard_handler_ = reinterpret_cast<KeyboardHandler*>(
    options_.GetPointerOption("keyboard handler"));
}

const int Viewer::AddStaticObject(const Eigen::Matrix3Xf& vertex,
  const Eigen::Matrix3Xi& face, const Option& options) {
  StaticOpenglShape* new_object = new StaticOpenglShape();
  new_object->Initialize(vertex, face, options);
  objects_[next_object_id_++] = new_object;
  return next_object_id_ - 1;
}

const int Viewer::AddStaticObject(
  const std::vector<std::vector<float>>& vertex,
  const std::vector<std::vector<int>>& face,
  const Option& options) {
  const int v_num = static_cast<int>(vertex[0].size());
  const int f_num = static_cast<int>(face[0].size());
  Eigen::Matrix3Xf vertices(3, v_num);
  Eigen::Matrix3Xi faces(3, f_num);
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < v_num; ++j)
      vertices(i, j) = vertex[i][j];
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < f_num; ++j)
      faces(i, j) = face[i][j];
  return AddStaticObject(vertices, faces, options);
}

void Viewer::UpdateStaticObject(const int object_id, const Option& options) {
  objects_[object_id]->Update(options);
}

const int Viewer::AddDynamicObject(const Eigen::Matrix3Xf& vertex,
  const Eigen::Matrix3Xi& face, Animator* const animator,
  const Option& options) {
  DynamicOpenglShape* new_object = new DynamicOpenglShape();
  new_object->Initialize(vertex, face, animator, options);
  objects_[next_object_id_++] = new_object;
  return next_object_id_ - 1;
}

const int Viewer::AddDynamicObject(
  const std::vector<std::vector<float>>& vertex,
  const std::vector<std::vector<int>>& face,
  SamplingAnimator* const animator,
  const Option& options) {
  const int v_num = static_cast<int>(vertex[0].size());
  const int f_num = static_cast<int>(face[0].size());
  Eigen::Matrix3Xf vertices(3, v_num);
  Eigen::Matrix3Xi faces(3, f_num);
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < v_num; ++j)
      vertices(i, j) = vertex[i][j];
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < f_num; ++j)
      faces(i, j) = face[i][j];
  return AddDynamicObject(vertices, faces, dynamic_cast<Animator*>(animator), options);
}

void Viewer::RemoveObject(const int object_id) {
  assert(object_id >= 0 && object_id < next_object_id_);
  if (objects_.find(object_id) != objects_.end()) {
    delete objects_[object_id];
    objects_.erase(object_id);
  }
}

void Viewer::AddStaticPointLight(const Eigen::Vector3f& position,
  const Option& options) {
  if (static_cast<int>(point_lights_.size()) < kMaxPointLightNum) {
    StaticOpenglLight* light = new StaticOpenglLight();
    light->Initialize(position, options);
    point_lights_.push_back(light);
  }
}

// For Python Binding.
void Viewer::AddStaticPointLight(const std::vector<float>& position,
  const Option& options) {
  Eigen::Vector3f pos;
  for (int i = 0; i < 3; ++i) pos(i) = position[i];
  AddStaticPointLight(pos, options);
}

void Viewer::AddDynamicPointLight(Animator* const animator,
  const Option& options) {
  if (static_cast<int>(point_lights_.size()) < kMaxPointLightNum) {
    DynamicOpenglLight* light = new DynamicOpenglLight();
    light->Initialize(animator, options);
    point_lights_.push_back(light);
  }
}

void Viewer::Run() {
  // Load shader programs.
  const std::string shader_folder = std::string(GRAPHICS_CODEBASE_SOURCE_DIR)
    + "/opengl_viewer/shader/";
  phong_shader_.InitializeFromSource(GeneratePhongModelVertexShader(*this),
    GeneratePhongModelFragShader(*this));
  const std::string geometry_shader = options_.GetBoolOption("shadow") ?
    shader_folder + "geometry_point_light_depth.glsl" : "";
  point_light_depth_shader_.InitializeFromFile(
    shader_folder + "vertex_point_light_depth.glsl",
    shader_folder + "fragment_point_light_depth.glsl",
    geometry_shader);

  // Use our shader.
  phong_shader_.UseShaderProgram();
  // Render on the whole framebuffer, complete from the lower left corner to
  // the upper right.
  int width = 0, height = 0;
  glfwGetFramebufferSize(window_, &width, &height);

  // Now set the light information.
  const int active_point_light_num =
    static_cast<int>(point_lights_.size());

  // Loop.
  int frame_idx = -1;
  float last_time = -1.0f;
  do {
    // Setup ImGui.
    if (imgui_wrapper_) {
      imgui_wrapper_->NewFrame();
      imgui_wrapper_->SetupUi();
    }

    // Bind the view and projection matrix.
    phong_shader_.SetUniformMatrix4f("view_matrix", view_matrix_);
    phong_shader_.SetUniformMatrix4f("projection_matrix",
      projection_matrix_);

    // Render the shadow.
    const float current_time = timer_ ? timer_->CurrentTime() : 0.0f;
    // Update frame count.
    if (current_time > last_time) {
      ++frame_idx;
    }

    if (current_time < -10.) {
        break;
    }
    
    last_time = current_time;

    if (options_.GetBoolOption("shadow")) {
        RenderShadow(current_time);
    }
    // Render the scene.
    glViewport(0, 0, width, height);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    phong_shader_.UseShaderProgram();
    // According to OpenGL samplers have to be uniforms.
    for (int i = 0; i < active_point_light_num; ++i) {
      const OpenglLight* light = point_lights_[i];
      const std::string name = "point_light[" + std::to_string(i) + "]";
      phong_shader_.SetUniform3f(name + ".position",
        light->PositionInWorld(current_time));
      phong_shader_.SetUniform3f(name + ".ambient", light->ambient());
      phong_shader_.SetUniform3f(name + ".diffuse", light->diffuse());
      phong_shader_.SetUniform3f(name + ".specular", light->specular());
      phong_shader_.SetUniform1f(name + ".max_dist", light->max_dist());
      // Bind the depth texture starting from texture 1.
      light->BindTexture(GL_TEXTURE0 + i);
      phong_shader_.SetUniform1i("shadow_map_sampler[" + std::to_string(i)
        + "]", i);
    }
    for (const auto& pair : objects_) {
      const OpenglShape* object = pair.second;
      phong_shader_.SetUniformMatrix4f("model_matrix",
        object->ModelMatrix(current_time));
      phong_shader_.SetUniformMatrix3f("normal_matrix",
        object->NormalMatrix(current_time));

      // Materials.
      phong_shader_.SetUniform3f("material.ambient", object->ambient());
      phong_shader_.SetUniform3f("material.diffuse", object->diffuse());
      phong_shader_.SetUniform3f("material.specular", object->specular());
      phong_shader_.SetUniform1f("material.shininess", object->shininess());

      // 1st attribute buffer: vertex positions.
      object->BindVertexBuffer(VertexAttribute::kVertex);
      // 2nd attribute buffer: normals.
      object->BindNormalBuffer(VertexAttribute::kNormal);
      // 3th attribute buffer: texture coordinates.
      object->BindTexture(VertexAttribute::kTexture,
        GL_TEXTURE0 + active_point_light_num);
      phong_shader_.SetUniform1i("texture_sampler", active_point_light_num);
      // Index buffer.
      object->BindElementBuffer();
      // Draw the triangles.
      object->Draw();

      // Cleanup texture binds.
      glActiveTexture(GL_TEXTURE0 + active_point_light_num);
      glBindTexture(GL_TEXTURE_2D, 0);
    }

    // Cleanup texture binds.
    for (int i = 0; i < active_point_light_num; ++i) {
      glActiveTexture(GL_TEXTURE0 + i);
      glBindTexture(GL_TEXTURE_CUBE_MAP, 0);
    }

    // Render GUI.
    if (imgui_wrapper_) imgui_wrapper_->Render();

    // Swap buffers.
    glfwSwapBuffers(window_);
    const std::string record_folder = options_.GetStringOption("record folder");
    if (!record_folder.empty()) {
      // Save screenshots.
      const GLenum format = GL_RGB;
      const int pixel_nbytes = 3;
      std::vector<GLubyte> pixels(pixel_nbytes * width * height);
      glReadPixels(0, 0, width, height, format, GL_UNSIGNED_BYTE, pixels.data());
      // Flip pixels in place.
      for (int i = 0; i < height / 2; ++i)
        for (int j = 0; j < width; ++j)
          for (int k = 0; k < 3; ++k) {
            const int old_pixel = pixel_nbytes * (i * width + j) + k;
            const int new_pixel = pixel_nbytes * ((height - 1 - i) * width + j) + k;
            GLubyte color = pixels[old_pixel];
            pixels[old_pixel] = pixels[new_pixel];
            pixels[new_pixel] = color;
          }
      const std::string png_name = record_folder + "/" + std::to_string(frame_idx) + ".png";
      stbi_write_png(png_name.c_str(), width, height, pixel_nbytes, pixels.data(), pixel_nbytes * width);
    }
    glfwPollEvents();
  } while (glfwGetKey(window_, GLFW_KEY_ESCAPE) != GLFW_PRESS
    && glfwWindowShouldClose(window_) == 0);
  // Check if the ESC key was pressed or the window was closed.
}

void Viewer::Cleanup() {
  // Shaders.
  phong_shader_.Cleanup();
  point_light_depth_shader_.Cleanup();

  // ImGui.
  if (imgui_wrapper_) {
    imgui_wrapper_->Cleanup();
    imgui_wrapper_ = NULL;
  }

  // Clean up VBO.
  for (const auto& pair : objects_) {
    delete pair.second;
  }
  objects_.clear();

  // Clean up VAO and textures.
  if (vertex_array_id_)
    glDeleteVertexArrays(1, &vertex_array_id_);

  // Clean up depth textures for the static point light.
  for (const auto& light : point_lights_) {
    delete light;
  }
  point_lights_.clear();

  // Close OpenGL window and terminate GLFW.
  if (window_) {
    glfwTerminate();
    window_ = NULL;
  }
}

void Viewer::InitializeOptions(const Option& option) {
  // Int.
  std::unordered_set<std::string> names = option.GetAllIntOptionNames();
  for (const std::string& name : names) {
    if (!options_.HasIntOption(name)) {
      std::cout << "Warning: Unsupported name: " << name << std::endl;
    } else {
      options_.SetIntOption(name, option.GetIntOption(name));
    }
  }
  // Float.
  names = option.GetAllFloatOptionNames();
  for (const std::string& name : names) {
    if (!options_.HasFloatOption(name)) {
      std::cout << "Warning: Unsupported name: " << name << std::endl;
    } else {
      options_.SetFloatOption(name, option.GetFloatOption(name));
    }
  }
  // Bool.
  names = option.GetAllBoolOptionNames();
  for (const std::string& name : names) {
    if (!options_.HasBoolOption(name)) {
      std::cout << "Warning: Unsupported name: " << name << std::endl;
    } else {
      options_.SetBoolOption(name, option.GetBoolOption(name));
    }
  }
  // String.
  names = option.GetAllStringOptionNames();
  for (const std::string& name : names) {
    if (!options_.HasStringOption(name)) {
      std::cout << "Warning: Unsupported name: " << name << std::endl;
    } else {
      options_.SetStringOption(name, option.GetStringOption(name));
    }
  }
  // Vector.
  names = option.GetAllVectorOptionNames();
  for (const std::string& name : names) {
    if (!options_.HasVectorOption(name)) {
      std::cout << "Warning: Unsupported name: " << name << std::endl;
    } else {
      options_.SetVectorOption(name, option.GetVectorOption(name));
    }
  }
  // Matrix.
  names = option.GetAllMatrixOptionNames();
  for (const std::string& name : names) {
    if (!options_.HasMatrixOption(name)) {
      std::cout << "Warning: Unsupported name: " << name << std::endl;
    } else {
      options_.SetMatrixOption(name, option.GetMatrixOption(name));
    }
  }
  // Pointer.
  names = option.GetAllPointerOptionNames();
  for (const std::string& name : names) {
    if (!options_.HasPointerOption(name)) {
      std::cout << "Warning: Unsupported name: " << name << std::endl;
    } else {
      options_.SetPointerOption(name, option.GetPointerOption(name));
    }
  }
}

const std::string Viewer::ToLowerCase(const std::string& name) const {
  std::string lower_name = name;
  const std::size_t length = lower_name.size();
  for (size_t i = 0; i < length; ++i) {
    if (lower_name[i] >= 'A' && lower_name[i] <= 'Z') {
      lower_name[i] -= 'A' - 'a';
    }
  }
  return lower_name;
}

void Viewer::RenderShadow(const float t) {
  // View matrix used in building the cube map.
  const Eigen::Vector3f unit_x = Eigen::Vector3f::UnitX(),
    unit_y = Eigen::Vector3f::UnitY(), unit_z = Eigen::Vector3f::UnitZ();
  const Eigen::Matrix<float, 3, 6> look_at = (Eigen::Matrix<float, 3, 6>()
    << unit_x, -unit_x, unit_y, -unit_y, unit_z, -unit_z).finished();
  const Eigen::Matrix<float, 3, 6> up = (Eigen::Matrix<float, 3, 6>()
    << -unit_y, -unit_y, unit_z, -unit_z, -unit_y, -unit_y).finished();

  // For each point light, we build a cube map to render shadows.
  const int point_light_num = static_cast<int>(point_lights_.size());
  for (int i = 0; i < point_light_num; ++i) {
    OpenglLight* light = point_lights_[i];
    const Eigen::Vector3f light_pos = light->PositionInWorld(t);
    // Analyze the distance between each static point light and objects.
    float near_plane = -1.0f, far_plane = 0.0f;
    for (const auto& pair : objects_) {
      const OpenglShape* object = pair.second;
      PointToBoundingBoxFrustum(light_pos, object->BoundingBoxInWorld(t),
        near_plane, far_plane);
    }
    // Perturb them by a little bit.
    near_plane *= 0.99f; far_plane *= 1.01f;
    light->set_max_dist(std::sqrt(3.0f) * far_plane);

    const int depth_map_size = light->depth_map_size();
    glViewport(0, 0, depth_map_size, depth_map_size);
    light->BindFrameBuffer();
    // Clear the depth buffer.
    glClear(GL_DEPTH_BUFFER_BIT);

    point_light_depth_shader_.UseShaderProgram();
    point_light_depth_shader_.SetUniform3f("light_position", light_pos);
    point_light_depth_shader_.SetUniform1f("max_dist", light->max_dist());
    for (int j = 0; j < 6; ++j) {
      const Eigen::Matrix4f view_matrix = LookAt(light_pos,
        light_pos + look_at.col(j), up.col(j));
      const Eigen::Matrix4f projection_matrix = Perspective(1.0f, 90.0f,
        near_plane, far_plane);
      const Eigen::Matrix4f shadow_matrix = projection_matrix * view_matrix;
      point_light_depth_shader_.SetUniformMatrix4f("shadow_matrix["
        + std::to_string(j) + "]", shadow_matrix);
    }

    // Now loop over all objects.
    for (const auto& pair : objects_) {
      const OpenglShape* object = pair.second;
      point_light_depth_shader_.SetUniformMatrix4f("model_matrix",
        object->ModelMatrix(t));
      object->BindVertexBuffer(VertexAttribute::kVertex);
      object->BindElementBuffer();
      // Draw the triangles.
      object->Draw();
    }
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
  }
}

void Viewer::MouseWheelScrollCallback(const float y_offset) {
  // Disable it when we are panning/rotation.
  if (mouse_wheel_pressed_) return;
  // +1 means zooming in, -1 zooming out.
  const float zoom_speed = options_.GetFloatOption("camera zoom speed");
  view_matrix_.topRightCorner(3, 1) +=
    Eigen::Vector3f::UnitZ() * zoom_speed * y_offset;

  // Forward it to ImGuiWrapper.
  if (imgui_wrapper_) imgui_wrapper_->MouseWheelScrollCallback(y_offset);

  // Custom code.
  if (mouse_handler_) {
    mouse_handler_->MouseWheelScrollCallback(y_offset);
  }
}

void Viewer::MouseButtonCallback(const int button, const int action) {
  if (button == GLFW_MOUSE_BUTTON_LEFT) {
    mouse_wheel_pressed_ = (action == GLFW_PRESS || action == GLFW_REPEAT);
  }

  // Forward it to ImGuiWrapper.
  if (imgui_wrapper_) imgui_wrapper_->MouseButtonCallback(button, action);

  // Custom code.
  if (mouse_handler_) {
    mouse_handler_->MouseButtonCallback(button, action);
  }
}

void Viewer::MouseCursorPosCallback(const float x_pos, const float y_pos) {
  const Eigen::Vector2i mouse_cursor_current_input(x_pos, y_pos);
  const Eigen::Vector2i mouse_diff =
    mouse_cursor_current_input - mouse_cursor_last_position_;

  if (mouse_wheel_pressed_) {
    // Extract R and p from the view matrix.
    // [R', -R'p; 0, 0, 0, 1].
    Eigen::Matrix3f R = view_matrix_.topLeftCorner(3, 3).transpose();
    Eigen::Vector3f p = -R * view_matrix_.topRightCorner(3, 1);
    if (shift_pressed_) {
      // Panning.
      const float panning_speed = options_.GetFloatOption("camera pan speed");
      const Eigen::Vector3f translate = R.leftCols(2) * panning_speed *
        Eigen::Vector2f(mouse_diff.x(), -mouse_diff.y());
      p -= translate;
      rotation_anchor_point_ -= translate;
    } else {
      // Rotation.
      const float rotate_speed =
        DegreeToRadian(options_.GetFloatOption("camera rotate speed"));
      const float angle_around_camera_x = -rotate_speed * mouse_diff.y();
      const float angle_around_camera_y = -rotate_speed * mouse_diff.x();
      // Build the local coordinates.
      Eigen::Vector3f rel_pos = p - rotation_anchor_point_;
      // Rotate around y axis first.
      const Eigen::Matrix3f rotate_y = Eigen::AngleAxisf(angle_around_camera_y,
        options_.GetVectorOption("camera up")).matrix();
      R = rotate_y * R;
      rel_pos = rotate_y * rel_pos;
      // Rotate around x axis.
      const Eigen::Matrix3f rotate_x =
        Eigen::AngleAxisf(angle_around_camera_x, R.col(0)).matrix();
      R = rotate_x * R;
      rel_pos = rotate_x * rel_pos;
      // Update p.
      p = rel_pos + rotation_anchor_point_;
    }
    view_matrix_.topLeftCorner(3, 3) = R.transpose();
    view_matrix_.topRightCorner(3, 1) = -R.transpose() * p;
  }
  mouse_cursor_last_position_ = mouse_cursor_current_input;

  // Custom code.
  if (mouse_handler_) {
    mouse_handler_->MouseCursorPosCallback(x_pos, y_pos);
  }
}

void Viewer::KeyCallback(const int key, const int action) {
  if (key == GLFW_KEY_LEFT_SHIFT) {
    shift_pressed_ = (action == GLFW_PRESS || action == GLFW_REPEAT);
  }

  // Forward it to ImGuiWrapper.
  if (imgui_wrapper_) imgui_wrapper_->KeyCallback(key, action);

  // Custom code.
  if (keyboard_handler_) {
    keyboard_handler_->KeyCallback(key, action);
  }
}

void Viewer::CharCallback(const unsigned int ch) {
  if (imgui_wrapper_) {
    imgui_wrapper_->CharCallback(ch);
  }
}

}
