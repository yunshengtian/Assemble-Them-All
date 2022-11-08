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
#include "imgui_wrapper.h"
#include "GL/glew.h"
#include "glfw3.h"
#ifdef _WIN32
#undef APIENTRY
#define GLFW_EXPOSE_NATIVE_WIN32
#define GLFW_EXPOSE_NATIVE_WGL
#include "glfw3native.h"
#endif
#include "imgui.h"

namespace opengl_viewer {

static const char* GetClipboardText(void* user_data) {
  return glfwGetClipboardString(static_cast<GLFWwindow*>(user_data));
}

static void SetClipboardText(void* user_data, const char* text) {
  glfwSetClipboardString(static_cast<GLFWwindow*>(user_data), text);
}

ImGuiWrapper::ImGuiWrapper()
  : window_(NULL), ui_shader_(),
  time_stamp_(0.0f),
  mouse_pressed_{false, false, false},
  mouse_wheel_(0.0f),
  vertex_array_id_(0), vertex_buffer_id_(0), element_id_(0),
  font_texture_id_(0),
  tex_location_(0), project_matirx_location_(0),
  pos_location_(0), uv_location_(0), color_location_(0) {}

void ImGuiWrapper::Initialize(GLFWwindow* window) {
  window_ = window;

  // Setup key mapping.
  ImGuiIO& io = ImGui::GetIO();
  io.KeyMap[ImGuiKey_Tab] = GLFW_KEY_TAB;
  io.KeyMap[ImGuiKey_LeftArrow] = GLFW_KEY_LEFT;
  io.KeyMap[ImGuiKey_RightArrow] = GLFW_KEY_RIGHT;
  io.KeyMap[ImGuiKey_UpArrow] = GLFW_KEY_UP;
  io.KeyMap[ImGuiKey_DownArrow] = GLFW_KEY_DOWN;
  io.KeyMap[ImGuiKey_PageUp] = GLFW_KEY_PAGE_UP;
  io.KeyMap[ImGuiKey_PageDown] = GLFW_KEY_PAGE_DOWN;
  io.KeyMap[ImGuiKey_Home] = GLFW_KEY_HOME;
  io.KeyMap[ImGuiKey_End] = GLFW_KEY_END;
  io.KeyMap[ImGuiKey_Delete] = GLFW_KEY_DELETE;
  io.KeyMap[ImGuiKey_Backspace] = GLFW_KEY_BACKSPACE;
  io.KeyMap[ImGuiKey_Enter] = GLFW_KEY_ENTER;
  io.KeyMap[ImGuiKey_Escape] = GLFW_KEY_ESCAPE;
  // Hot keys used in a text editor (Ctrl + A, + C, + V, etc).
  io.KeyMap[ImGuiKey_A] = GLFW_KEY_A;
  io.KeyMap[ImGuiKey_C] = GLFW_KEY_C;
  io.KeyMap[ImGuiKey_V] = GLFW_KEY_V;
  io.KeyMap[ImGuiKey_X] = GLFW_KEY_X;
  io.KeyMap[ImGuiKey_Y] = GLFW_KEY_Y;
  io.KeyMap[ImGuiKey_Z] = GLFW_KEY_Z;

  io.RenderDrawListsFn = NULL;
  io.SetClipboardTextFn = SetClipboardText;
  io.GetClipboardTextFn = GetClipboardText;
  io.ClipboardUserData = window_;
#ifdef _WIN32
  io.ImeWindowHandle = glfwGetWin32Window(window_);
#endif

  if (!font_texture_id_) {
    CreateDeviceObjects();
  }
}

void ImGuiWrapper::NewFrame() {
  ImGuiIO& io = ImGui::GetIO();

  // Setup display size (every frame to accommodate for window resizing).
  int width, height;
  int display_width, display_height;
  glfwGetWindowSize(window_, &width, &height);
  glfwGetFramebufferSize(window_, &display_width, &display_height);
  io.DisplaySize = ImVec2(static_cast<float>(width),
    static_cast<float>(height));
  io.DisplayFramebufferScale = ImVec2(
    width > 0 ? (display_width * 1.0f / width) : 0,
    height > 0 ? (display_height * 1.0f / height) : 0);

  // Setup time step.
  const float current_time = static_cast<float>(glfwGetTime());
  io.DeltaTime = time_stamp_ > 0.0f ? (current_time - time_stamp_)
    : (1.0f / 60.0f);
  time_stamp_ = current_time;

  // Setup inputs.
  if (glfwGetWindowAttrib(window_, GLFW_FOCUSED)) {
    // Mouse position in screen coordinates
    // (set to -1,-1 if no mouse / on another screen, etc.)
    double mouse_x, mouse_y;
    glfwGetCursorPos(window_, &mouse_x, &mouse_y);
    io.MousePos = ImVec2(static_cast<float>(mouse_x),
      static_cast<float>(mouse_y));
  } else {
    io.MousePos = ImVec2(-1.0f, -1.0f);
  }

  for (int i = 0; i < 3; i++) {
    // If a mouse press event came, always pass it as "mouse held this frame",
    // so we don't miss click-release events that are shorter than 1 frame.
    io.MouseDown[i] = mouse_pressed_[i]
      || glfwGetMouseButton(window_, i) != 0;
    mouse_pressed_[i] = false;
  }

  io.MouseWheel = mouse_wheel_;
  mouse_wheel_ = 0.0f;

  // Hide OS mouse cursor if ImGui is drawing it.
  glfwSetInputMode(window_, GLFW_CURSOR, io.MouseDrawCursor
    ? GLFW_CURSOR_HIDDEN : GLFW_CURSOR_NORMAL);

  // Start the frame.
  ImGui::NewFrame();
}

void ImGuiWrapper::Render() {
  ImGui::Render();
  RenderDrawLists(ImGui::GetDrawData());
}

void ImGuiWrapper::CreateDeviceObjects() {
  // Backup GL states.
  GLint last_texture, last_array_buffer, last_vertex_array;
  glGetIntegerv(GL_TEXTURE_BINDING_2D, &last_texture);
  glGetIntegerv(GL_ARRAY_BUFFER_BINDING, &last_array_buffer);
  glGetIntegerv(GL_VERTEX_ARRAY_BINDING, &last_vertex_array);

  // Load shaders.
  const std::string shader_folder = std::string(GRAPHICS_CODEBASE_SOURCE_DIR)
    + "/opengl_viewer/shader/";
  ui_shader_.InitializeFromFile(shader_folder + "vertex_imgui.glsl",
    shader_folder + "fragment_imgui.glsl");

  tex_location_ = ui_shader_.GetUniformLocation("texture");
  project_matirx_location_ = ui_shader_.GetUniformLocation("project_matrix");
  pos_location_ = ui_shader_.GetAttribLocation("position");
  uv_location_ = ui_shader_.GetAttribLocation("uv");
  color_location_ = ui_shader_.GetAttribLocation("color");

  glGenBuffers(1, &vertex_buffer_id_);
  glGenBuffers(1, &element_id_);

  glGenVertexArrays(1, &vertex_array_id_);
  glBindVertexArray(vertex_array_id_);
  glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer_id_);
  glEnableVertexAttribArray(pos_location_);
  glEnableVertexAttribArray(uv_location_);
  glEnableVertexAttribArray(color_location_);

#define OFFSETOF(TYPE, ELEMENT) \
  reinterpret_cast<GLvoid*>(&(static_cast<TYPE *>(0)->ELEMENT))
  glVertexAttribPointer(pos_location_, 2, GL_FLOAT, GL_FALSE,
    sizeof(ImDrawVert), OFFSETOF(ImDrawVert, pos));
  glVertexAttribPointer(uv_location_, 2, GL_FLOAT, GL_FALSE,
    sizeof(ImDrawVert), OFFSETOF(ImDrawVert, uv));
  glVertexAttribPointer(color_location_, 4, GL_UNSIGNED_BYTE, GL_TRUE,
    sizeof(ImDrawVert), OFFSETOF(ImDrawVert, col));
#undef OFFSETOF

  // Load font texture.
  CreateFontsTexture();

  // Restore modified GL state
  glBindTexture(GL_TEXTURE_2D, last_texture);
  glBindBuffer(GL_ARRAY_BUFFER, last_array_buffer);
  glBindVertexArray(last_vertex_array);
}

void ImGuiWrapper::CreateFontsTexture() {
  // Build texture atlas.
  ImGuiIO& io = ImGui::GetIO();
  unsigned char* pixels;
  int width, height;
  // Load as RGBA 32-bits (75% of the memory is wasted, but default font is so
  // small) because it is more likely to be compatible with user's existing
  // shaders. If your ImTextureId represent a higher-level concept than just a
  // GL texture id, consider calling GetTexDataAsAlpha8() instead to save on
  // GPU memory.
  io.Fonts->GetTexDataAsRGBA32(&pixels, &width, &height);

  // Upload texture to graphics system.
  GLint last_texture;
  glGetIntegerv(GL_TEXTURE_BINDING_2D, &last_texture);
  glGenTextures(1, &font_texture_id_);
  glBindTexture(GL_TEXTURE_2D, font_texture_id_);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA,
    GL_UNSIGNED_BYTE, pixels);

  // Store our identifier.
  io.Fonts->TexID = (void *)(intptr_t)font_texture_id_;

  // Restore state.
  glBindTexture(GL_TEXTURE_2D, last_texture);
}

void ImGuiWrapper::Cleanup() {
  if (vertex_array_id_) glDeleteVertexArrays(1, &vertex_array_id_);
  if (vertex_buffer_id_) glDeleteBuffers(1, &vertex_buffer_id_);
  if (element_id_) glDeleteBuffers(1, &element_id_);
  vertex_array_id_ = vertex_buffer_id_ = element_id_ = 0;

  if (font_texture_id_) {
    glDeleteTextures(1, &font_texture_id_);
    ImGui::GetIO().Fonts->TexID = 0;
    font_texture_id_ = 0;
  }

  ui_shader_.Cleanup();
  ImGui::Shutdown();
}

void ImGuiWrapper::RenderDrawLists(ImDrawData* draw_data) {
  // Avoid rendering when minimized, scale coordinates for retina displays
  // (screen coordinates != framebuffer coordinates)
  ImGuiIO& io = ImGui::GetIO();
  int framebuffer_width =
    static_cast<int>(io.DisplaySize.x * io.DisplayFramebufferScale.x);
  int framebuffer_height =
    static_cast<int>(io.DisplaySize.y * io.DisplayFramebufferScale.y);
  if (framebuffer_width == 0 || framebuffer_height == 0) return;
  draw_data->ScaleClipRects(io.DisplayFramebufferScale);

  // Backup GL states.
  GLint last_active_texture;
  glGetIntegerv(GL_ACTIVE_TEXTURE, &last_active_texture);
  glActiveTexture(GL_TEXTURE0);
  GLint last_program;
  glGetIntegerv(GL_CURRENT_PROGRAM, &last_program);
  GLint last_texture;
  glGetIntegerv(GL_TEXTURE_BINDING_2D, &last_texture);
  GLint last_array_buffer;
  glGetIntegerv(GL_ARRAY_BUFFER_BINDING, &last_array_buffer);
  GLint last_element_array_buffer;
  glGetIntegerv(GL_ELEMENT_ARRAY_BUFFER_BINDING, &last_element_array_buffer);
  GLint last_vertex_array;
  glGetIntegerv(GL_VERTEX_ARRAY_BINDING, &last_vertex_array);
  GLint last_blend_src_rgb;
  glGetIntegerv(GL_BLEND_SRC_RGB, &last_blend_src_rgb);
  GLint last_blend_dst_rgb;
  glGetIntegerv(GL_BLEND_DST_RGB, &last_blend_dst_rgb);
  GLint last_blend_src_alpha;
  glGetIntegerv(GL_BLEND_SRC_ALPHA, &last_blend_src_alpha);
  GLint last_blend_dst_alpha;
  glGetIntegerv(GL_BLEND_DST_ALPHA, &last_blend_dst_alpha);
  GLint last_blend_equation_rgb;
  glGetIntegerv(GL_BLEND_EQUATION_RGB, &last_blend_equation_rgb);
  GLint last_blend_equation_alpha;
  glGetIntegerv(GL_BLEND_EQUATION_ALPHA, &last_blend_equation_alpha);
  GLint last_viewport[4];
  glGetIntegerv(GL_VIEWPORT, last_viewport);
  GLint last_scissor_box[4];
  glGetIntegerv(GL_SCISSOR_BOX, last_scissor_box);
  GLboolean last_enable_blend = glIsEnabled(GL_BLEND);
  GLboolean last_enable_cull_face = glIsEnabled(GL_CULL_FACE);
  GLboolean last_enable_depth_test = glIsEnabled(GL_DEPTH_TEST);
  GLboolean last_enable_scissor_test = glIsEnabled(GL_SCISSOR_TEST);

  // Setup render state: Alpha-blending enabled, no face culling, no depth
  // testing, scissor enabled.
  glEnable(GL_BLEND);
  glBlendEquation(GL_FUNC_ADD);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glDisable(GL_CULL_FACE);
  glDisable(GL_DEPTH_TEST);
  glEnable(GL_SCISSOR_TEST);

  // Setup viewport, orthographic projection matrix.
  glViewport(0, 0, static_cast<GLsizei>(framebuffer_width),
    static_cast<GLsizei>(framebuffer_height));
  const float ortho_projection[4][4] = {
    { 2.0f / io.DisplaySize.x,  0.0f, 0.0f, 0.0f },
    { 0.0f, 2.0f / -io.DisplaySize.y, 0.0f, 0.0f },
    { 0.0f, 0.0f, -1.0f, 0.0f },
    { -1.0f, 1.0f, 0.0f, 1.0f },
  };
  ui_shader_.UseShaderProgram();
  glUniform1i(tex_location_, 0);
  glUniformMatrix4fv(project_matirx_location_, 1,
    GL_FALSE, &ortho_projection[0][0]);
  glBindVertexArray(vertex_array_id_);

  for (int n = 0; n < draw_data->CmdListsCount; n++) {
    const ImDrawList* cmd_list = draw_data->CmdLists[n];
    const ImDrawIdx* idx_buffer_offset = 0;

    glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer_id_);
    glBufferData(GL_ARRAY_BUFFER,
      static_cast<GLsizeiptr>(cmd_list->VtxBuffer.Size) * sizeof(ImDrawVert),
      static_cast<const GLvoid*>(cmd_list->VtxBuffer.Data), GL_STREAM_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, element_id_);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,
      static_cast<GLsizeiptr>(cmd_list->IdxBuffer.Size) * sizeof(ImDrawIdx),
      static_cast<const GLvoid*>(cmd_list->IdxBuffer.Data), GL_STREAM_DRAW);

    for (int cmd_i = 0; cmd_i < cmd_list->CmdBuffer.Size; cmd_i++) {
      const ImDrawCmd* pcmd = &cmd_list->CmdBuffer[cmd_i];
      if (pcmd->UserCallback) {
        pcmd->UserCallback(cmd_list, pcmd);
      } else {
        glBindTexture(GL_TEXTURE_2D,
          static_cast<GLuint>(reinterpret_cast<intptr_t>(pcmd->TextureId)));
        glScissor(static_cast<int>(pcmd->ClipRect.x),
          static_cast<int>(framebuffer_height - pcmd->ClipRect.w),
          static_cast<int>(pcmd->ClipRect.z - pcmd->ClipRect.x),
          static_cast<int>(pcmd->ClipRect.w - pcmd->ClipRect.y));
        glDrawElements(GL_TRIANGLES, static_cast<GLsizei>(pcmd->ElemCount),
          sizeof(ImDrawIdx) == 2 ? GL_UNSIGNED_SHORT : GL_UNSIGNED_INT,
          idx_buffer_offset);
      }
      idx_buffer_offset += pcmd->ElemCount;
    }
  }

  // Restore modified GL state
  glUseProgram(last_program);
  glBindTexture(GL_TEXTURE_2D, last_texture);
  glActiveTexture(last_active_texture);
  glBindVertexArray(last_vertex_array);
  glBindBuffer(GL_ARRAY_BUFFER, last_array_buffer);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, last_element_array_buffer);
  glBlendEquationSeparate(last_blend_equation_rgb, last_blend_equation_alpha);
  glBlendFuncSeparate(last_blend_src_rgb, last_blend_dst_rgb,
    last_blend_src_alpha, last_blend_dst_alpha);
  if (last_enable_blend) glEnable(GL_BLEND);
  else glDisable(GL_BLEND);
  if (last_enable_cull_face) glEnable(GL_CULL_FACE);
  else glDisable(GL_CULL_FACE);
  if (last_enable_depth_test) glEnable(GL_DEPTH_TEST);
  else glDisable(GL_DEPTH_TEST);
  if (last_enable_scissor_test) glEnable(GL_SCISSOR_TEST);
  else glDisable(GL_SCISSOR_TEST);
  glViewport(last_viewport[0], last_viewport[1],
    static_cast<GLsizei>(last_viewport[2]),
    static_cast<GLsizei>(last_viewport[3]));
  glScissor(last_scissor_box[0], last_scissor_box[1],
    static_cast<GLsizei>(last_scissor_box[2]),
    static_cast<GLsizei>(last_scissor_box[3]));
}

void ImGuiWrapper::MouseWheelScrollCallback(const float y_offset) {
  // Use fractional mouse wheel, 1.0 unit 5 lines.
  mouse_wheel_ += static_cast<float>(y_offset);
}

void ImGuiWrapper::MouseButtonCallback(const int button, const int action) {
  if (action == GLFW_PRESS && button >= 0 && button < 3)
    mouse_pressed_[button] = true;
}

void ImGuiWrapper::KeyCallback(const int key, const int action) {
  ImGuiIO& io = ImGui::GetIO();
  if (action == GLFW_PRESS)
    io.KeysDown[key] = true;
  if (action == GLFW_RELEASE)
    io.KeysDown[key] = false;

  io.KeyCtrl = io.KeysDown[GLFW_KEY_LEFT_CONTROL] || io.KeysDown[GLFW_KEY_RIGHT_CONTROL];
  io.KeyShift = io.KeysDown[GLFW_KEY_LEFT_SHIFT] || io.KeysDown[GLFW_KEY_RIGHT_SHIFT];
  io.KeyAlt = io.KeysDown[GLFW_KEY_LEFT_ALT] || io.KeysDown[GLFW_KEY_RIGHT_ALT];
  io.KeySuper = io.KeysDown[GLFW_KEY_LEFT_SUPER] || io.KeysDown[GLFW_KEY_RIGHT_SUPER];
}

void ImGuiWrapper::CharCallback(unsigned int ch) {
  ImGuiIO& io = ImGui::GetIO();
  if (ch > 0u && ch < 0x10000u)
    io.AddInputCharacter(static_cast<unsigned short>(ch));
}

} // opengl_viewer