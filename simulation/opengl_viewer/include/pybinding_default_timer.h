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
#ifndef _OPENGL_VIEWER_PYBINDING_DEFAULT_TIMER_H_
#define _OPENGL_VIEWER_PYBINDING_DEFAULT_TIMER_H_

#include "pybinding_default_keyboard_handler.h"
#include "timer.h"

namespace opengl_viewer {

class PyBindingDefaultTimer : public Timer {
public:
  PyBindingDefaultTimer() : fps_(25), dt_(0.04f), current_time_(0.0f),
      keyboard_handler_(nullptr) {}

  void Initialize(const int fps, PyBindingDefaultKeyboardHandler* keyboard_handler) {
    fps_ = fps;
    dt_ = 1.0f / fps;
    current_time_ = 0.0f;
    keyboard_handler_ = keyboard_handler;
  }

  const float CurrentTime() {
    if (keyboard_handler_->reset()) {
      current_time_ = 0;
      keyboard_handler_->set_reset(false);
    }
    const float current_time = current_time_;
    if (!keyboard_handler_->paused()) current_time_ += dt_;
    return current_time;
  }

private:
  int fps_;
  float dt_;
  float current_time_;
  PyBindingDefaultKeyboardHandler* keyboard_handler_;
};

}

#endif