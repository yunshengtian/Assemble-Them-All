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
#ifndef _OPENGL_VIEWER_LINEAR_TIMER_H_
#define _OPENGL_VIEWER_LINEAR_TIMER_H_

#include "timer.h"

namespace opengl_viewer {

class LinearTimer : public Timer {
public:
  LinearTimer() : fps_(25), dt_(0.04f), current_time_(0.0f) {}

  void Initialize(const int fps) {
    fps_ = fps;
    dt_ = 1.0f / fps;
    current_time_ = 0.0f;
  }

  const float CurrentTime() {
    const float current_time = current_time_; 
    current_time_ += dt_;
    return current_time;
  }

private:
  int fps_;
  float dt_;
  float current_time_;
};

}

#endif