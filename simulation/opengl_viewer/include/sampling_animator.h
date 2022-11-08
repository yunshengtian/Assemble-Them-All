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
#ifndef _OPENGL_VIEWER_SAMPLING_ANIMATOR_H_
#define _OPENGL_VIEWER_SAMPLING_ANIMATOR_H_

#include <map>
#include "Eigen/Dense"
#include "animator.h"

namespace opengl_viewer {

class SamplingAnimator : public Animator {
public:
  explicit SamplingAnimator(const int sample_buffer)
    : Animator(), sample_buffer_(sample_buffer) {}

  void AddSample(const float t, const Eigen::Matrix4f& transform) {
    samples_[t] = transform;
    while (sample_buffer_ >= 0 && static_cast<int>(samples_.size()) > sample_buffer_)
      samples_.erase(samples_.begin());
  }

  // For Python binding.
  void AddSample(const float t, const std::vector<std::vector<float>>& transform) {
    Eigen::Matrix4f mat;
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j)
        mat(i, j) = transform[i][j];
    AddSample(t, mat);
  }

  const Eigen::Matrix4f AnimatedModelMatrix(const float t) {
    auto ptr = samples_.lower_bound(t);
    if (ptr == samples_.end()) {
      return samples_.rbegin()->second;
    } else if (ptr->first == t || ptr == samples_.begin()) {
      return ptr->second;
    } else {
      std::advance(ptr, -1);
      return ptr->second;
    }
  }

private:
  int sample_buffer_;
  std::map<float, Eigen::Matrix4f> samples_;
};

}

#endif
