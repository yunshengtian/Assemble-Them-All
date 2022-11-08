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
#ifndef _OPENGL_VIEWER_OPTION_H_
#define _OPENGL_VIEWER_OPTION_H_

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "Eigen/Dense"

namespace opengl_viewer {

class Option {
public:
  // Clear all options.
  void Clear() {
    ClearIntOption();
    ClearFloatOption();
    ClearBoolOption();
    ClearStringOption();
    ClearVectorOption();
    ClearMatrixOption();
    ClearPointerOption();
  }

  // Int options.
  void SetIntOption(const std::string& name, const int value) {
    int_options_[name] = value;
  }
  void ClearIntOption(const std::string& name) {
    int_options_.erase(name);
  }
  void ClearIntOption() { int_options_.clear(); }
  const bool HasIntOption(const std::string& name) const {
    return int_options_.find(name) != int_options_.end();
  }
  const int GetIntOption(const std::string& name) const {
    return int_options_.at(name);
  }
  const std::unordered_set<std::string> GetAllIntOptionNames() const {
    std::unordered_set<std::string> name;
    for (const auto& pair : int_options_) {
      name.insert(pair.first);
    }
    return name;
  }

  // Float options.
  void SetFloatOption(const std::string& name, const float value) {
    float_options_[name] = value;
  }
  void ClearFloatOption(const std::string& name) {
    float_options_.erase(name);
  }
  void ClearFloatOption() { float_options_.clear(); }
  const bool HasFloatOption(const std::string& name) const {
    return float_options_.find(name) != float_options_.end();
  }
  const float GetFloatOption(const std::string& name) const {
    return float_options_.at(name);
  }
  const std::unordered_set<std::string> GetAllFloatOptionNames() const {
    std::unordered_set<std::string> name;
    for (const auto& pair : float_options_) {
      name.insert(pair.first);
    }
    return name;
  }

  // Bool options.
  void SetBoolOption(const std::string& name, const bool value) {
    bool_options_[name] = value;
  }
  void ClearBoolOption(const std::string& name) {
    bool_options_.erase(name);
  }
  void ClearBoolOption() { bool_options_.clear(); }
  const bool HasBoolOption(const std::string& name) const {
    return bool_options_.find(name) != bool_options_.end();
  }
  const bool GetBoolOption(const std::string& name) const {
    return bool_options_.at(name);
  }
  const std::unordered_set<std::string> GetAllBoolOptionNames() const {
    std::unordered_set<std::string> name;
    for (const auto& pair : bool_options_) {
      name.insert(pair.first);
    }
    return name;
  }

  // String options.
  void SetStringOption(const std::string& name, const std::string& value) {
    string_options_[name] = value;
  }
  void ClearStringOption(const std::string& name) {
    string_options_.erase(name);
  }
  void ClearStringOption() { string_options_.clear(); }
  const bool HasStringOption(const std::string& name) const {
    return string_options_.find(name) != string_options_.end();
  }
  const std::string GetStringOption(const std::string& name) const {
    return string_options_.at(name);
  }
  const std::unordered_set<std::string> GetAllStringOptionNames() const {
    std::unordered_set<std::string> name;
    for (const auto& pair : string_options_) {
      name.insert(pair.first);
    }
    return name;
  }

  // Vector options.
  void SetVectorOption(const std::string& name, const Eigen::VectorXf& value) {
    vector_options_[name] = value;
  }
  void SetVectorOption(const std::string& name, const std::vector<float>& value) {
    Eigen::VectorXf val(value.size());
    for (auto i = 0; i < val.size(); ++i) val(i) = value[i];
    SetVectorOption(name, val);
  }
  // More convenient way for initialization.
  void SetVectorOption(const std::string& name,
    const float x, const float y) {
    vector_options_[name] = Eigen::Vector2f(x, y);
  }
  void SetVectorOption(const std::string& name,
    const float x, const float y, const float z) {
    vector_options_[name] = Eigen::Vector3f(x, y, z);
  }
  void SetVectorOption(const std::string& name,
    const float x, const float y, const float z, const float w) {
    vector_options_[name] = Eigen::Vector4f(x, y, z, w);
  }
  void ClearVectorOption(const std::string& name) {
    vector_options_.erase(name);
  }
  void ClearVectorOption() { vector_options_.clear(); }
  const bool HasVectorOption(const std::string& name) const {
    return vector_options_.find(name) != vector_options_.end();
  }
  const Eigen::VectorXf GetVectorOption(const std::string& name) const {
    return vector_options_.at(name);
  }
  // For Python binding only.
  const std::vector<float> GetVectorOptionPyBinding(const std::string& name) const {
    const Eigen::VectorXf val = GetVectorOption(name);
    return std::vector<float>(val.data(), val.data() + val.size());
  }
  const std::unordered_set<std::string> GetAllVectorOptionNames() const {
    std::unordered_set<std::string> name;
    for (const auto& pair : vector_options_) {
      name.insert(pair.first);
    }
    return name;
  }

  // Matrix options.
  void SetMatrixOption(const std::string& name, const Eigen::MatrixXf& value) {
    matrix_options_[name] = value;
  }
  void SetMatrixOption(const std::string& name, const std::vector<std::vector<float>>& value) {
    Eigen::MatrixXf val(value.size(), value[0].size());
    for (int i = 0; i < static_cast<int>(value.size()); ++i)
      for (int j = 0; j < static_cast<int>(value[0].size()); ++j)
        val(i, j) = value[i][j];
    matrix_options_[name] = val;
  }
  void ClearMatrixOption(const std::string& name) {
    matrix_options_.erase(name);
  }
  void ClearMatrixOption() { matrix_options_.clear(); }
  const bool HasMatrixOption(const std::string& name) const {
    return matrix_options_.find(name) != matrix_options_.end();
  }
  const Eigen::MatrixXf GetMatrixOption(const std::string& name) const {
    return matrix_options_.at(name);
  }
  const std::vector<std::vector<float>> GetMatrixOptionPyBinding(const std::string& name) const {
    const Eigen::MatrixXf val = matrix_options_.at(name);
    std::vector<std::vector<float>> value(val.rows(), std::vector<float>(val.cols(), 0));
    for (int i = 0; i < static_cast<int>(val.rows()); ++i)
      for (int j = 0; j < static_cast<int>(val.cols()); ++j)
        value[i][j] = val(i, j);
    return value;
  }
  const std::unordered_set<std::string> GetAllMatrixOptionNames() const {
    std::unordered_set<std::string> name;
    for (const auto& pair : matrix_options_) {
      name.insert(pair.first);
    }
    return name;
  }

  // Pointer options.
  void SetPointerOption(const std::string& name, void* const value) {
    pointer_options_[name] = value;
  }
  void ClearPointerOption(const std::string& name) {
    pointer_options_.erase(name);
  }
  void ClearPointerOption() { pointer_options_.clear(); }
  const bool HasPointerOption(const std::string& name) const {
    return pointer_options_.find(name) != pointer_options_.end();
  }
  void* const GetPointerOption(const std::string& name) const {
    return pointer_options_.at(name);
  }
  const std::unordered_set<std::string> GetAllPointerOptionNames() const {
    std::unordered_set<std::string> name;
    for (const auto& pair : pointer_options_) {
      name.insert(pair.first);
    }
    return name;
  }

private:
  std::unordered_map<std::string, int> int_options_;
  std::unordered_map<std::string, float> float_options_;
  std::unordered_map<std::string, bool> bool_options_;
  std::unordered_map<std::string, std::string> string_options_;
  std::unordered_map<std::string, Eigen::VectorXf> vector_options_;
  std::unordered_map<std::string, Eigen::MatrixXf> matrix_options_;
  std::unordered_map<std::string, void*> pointer_options_;
};

}

#endif