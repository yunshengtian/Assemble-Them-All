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
#ifndef _OPENGL_VIEWER_PYBINDING_DEFAULT_KEYBOARD_HANDLER_H_
#define _OPENGL_VIEWER_PYBINDING_DEFAULT_KEYBOARD_HANDLER_H_

#include "keyboard_handler.h"

namespace opengl_viewer {

class PyBindingDefaultKeyboardHandler : public KeyboardHandler {
public:
  PyBindingDefaultKeyboardHandler()
  : KeyboardHandler(), paused_(false), reset_(false) {}

  void Initialize(const bool paused) { paused_ = paused; }

  void KeyCallback(const int key, const int action) {
    if (key == GLFW_KEY_P && action == GLFW_PRESS) paused_ = !paused_;

    if (key == GLFW_KEY_R && action == GLFW_PRESS) reset_ = true;
  }

  const bool paused() const { return paused_; }
  const bool reset() const { return reset_; }
  void set_reset(const bool reset) { reset_ = reset; }

private:
  bool paused_;
  bool reset_;
};

}

#endif