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
#include "image.h"
// As requested in stb_image.h.
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

namespace opengl_viewer {

void Image::Initialize(const int row_num, const int col_num,
  const float rgb_data[]) {
  row_num_ = row_num;
  col_num_ = col_num;
  rgb_data_ = Eigen::Map<Eigen::Matrix3Xf>(
    const_cast<float *>(rgb_data), 3, row_num * col_num);
}

void Image::Initialize(
  const std::vector<std::vector<Eigen::Vector3f>>& rgb_data) {
  row_num_ = static_cast<int>(rgb_data.size());
  col_num_ = row_num_ == 0 ? 0 : static_cast<int>(rgb_data.at(0).size());
  rgb_data_ = Eigen::Matrix3Xf::Zero(3, row_num_ * col_num_);
  for (int i = 0; i < row_num_; ++i)
    for (int j = 0; j < col_num_; ++j)
      rgb_data_.col(i * col_num_ + j) = rgb_data[i][j];
}

void Image::Initialize(const std::string& image_path) {
  int channels = 0;
  // Enforce the number of output channels to be RGB.
  unsigned char* rgb_data = stbi_load(image_path.c_str(),
    &col_num_, &row_num_, &channels, 3);
  if (rgb_data) {
    rgb_data_ = Eigen::Matrix3Xf::Zero(3, row_num_ * col_num_);
    for (int i = 0; i < row_num_; ++i)
      for (int j = 0; j < col_num_; ++j)
        for (int k = 0; k < 3; ++k) {
          rgb_data_(k, i * col_num_ + j) =
            static_cast<float>(rgb_data[3 * (i * col_num_ + j) + k]) / 255.0f;
        }
  }
  stbi_image_free(rgb_data);
}

}