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
#ifndef _OPENGL_VIEWER_IMAGE_H_
#define _OPENGL_VIEWER_IMAGE_H_

#include <vector>
#include "Eigen/Dense"

namespace opengl_viewer {

class Image {
public:
  Image() : rgb_data_(3, 0), row_num_(0), col_num_(0) {}

  // I assume rgb_data is in the following order:
  // r00, g00, b00, r01, g01, b01, r02, g02, b02, ...
  // where rij means the red channel for pixel in i-th row and j-th column.
  void Initialize(const int row_num, const int col_num,
    const float rgb_data[]);
  // rgb_data[i][j][0..2] is the red, green and blue channel for pixel in i-th
  // row and j-th column.
  void Initialize(const std::vector<std::vector<Eigen::Vector3f>>& rgb_data);
  // Load an image file.
  void Initialize(const std::string& image_path);

  const Eigen::Matrix3Xf rgb_data() const { return rgb_data_; }
  const int row_num() const { return row_num_; }
  const int col_num() const { return col_num_; }

private:
  Eigen::Matrix3Xf rgb_data_;
  int row_num_, col_num_;
};

}

#endif