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
#include "geometry.h"
#include <fstream>
#include <vector>
#include "bounding_box.h"

namespace opengl_viewer {

const Eigen::Matrix4f LookAt(const Eigen::Vector3f& position,
  const Eigen::Vector3f& lookat, const Eigen::Vector3f& up) {
  // x, y and z axes.
  const Eigen::Vector3f view_z = (position - lookat).normalized();
  const Eigen::Vector3f view_x = up.cross(view_z).normalized();
  const Eigen::Vector3f view_y = view_z.cross(view_x).normalized();
  const Eigen::Matrix3f R = (Eigen::Matrix3f() << view_x, view_y, view_z).finished();
  return (Eigen::Matrix4f() << R.transpose(), -R.transpose() * position,
    Eigen::RowVector3f::Zero(), 1.0f).finished();
}

const Eigen::Matrix4f Perspective(const float aspect_ratio,
  const float field_of_view, const float z_min, const float z_max) {
  const float tan_half_fov = std::tan(DegreeToRadian(field_of_view) / 2.0f);
  Eigen::Matrix4f matrix = Eigen::Matrix4f::Zero();
  matrix(0, 0) = 1.0f / (aspect_ratio * tan_half_fov);
  matrix(1, 1) = 1.0f / tan_half_fov;
  matrix(2, 2) = -(z_max + z_min) / (z_max - z_min);
  matrix(3, 2) = -1.0f;
  matrix(2, 3) = -(2.0f * z_min * z_max) / (z_max - z_min);
  return matrix;
}

void Sphere(const float radius, const int slices, const int stacks,
  Eigen::Matrix3Xf& vertex, Eigen::Matrix3Xi& face, Eigen::Matrix2Xf& uv) {
  // Check the arguments.
  assert(radius >= 0.0f);
  assert(slices >= 3);
  assert(stacks >= 2);
  // We first draw a unit sphere then scale it.
  // theta \in [0, 2pi), phi \in [-\pi / 2, pi / 2].
  const float dtheta = 2.0f * kPi / slices;
  const float dphi = kPi / stacks;
  vertex = Eigen::Matrix3Xf::Zero(3, slices * 2 + (stacks - 1) * (slices + 1));
  face = Eigen::Matrix3Xi::Zero(3, (stacks - 1) * slices * 2);
  uv = vertex.topRows(2);
  vertex.leftCols(slices) =
    Eigen::Vector3f(0.0f, 0.0f, 1.0f).replicate(1, slices);
  vertex.rightCols(slices) =
    Eigen::Vector3f(0.0f, 0.0f, -1.0f).replicate(1, slices);
  // Initialize u and v.
  const float u_step = 1.0f / slices;
  const float v_step = 1.0f / stacks;
  for (int i = 0; i < slices; ++i) {
    const float u_coord = u_step * (i + 0.5f);
    uv.col(i) = Eigen::Vector2f(u_coord, 1.0f);
    uv.col(slices + (stacks - 1) * (slices + 1) + i)
      = Eigen::Vector2f(u_coord, 0.0f);
  }

  // The top-z triangles.
  int cur_vertex_id = slices, cur_face_id = 0;
  for (int i = 0; i < slices; ++i) {
    vertex.col(cur_vertex_id++) = Eigen::Vector3f(
      cos(kPi / 2.0f - dphi) * cos(i * dtheta),
      cos(kPi / 2.0f - dphi) * sin(i * dtheta),
      sin(kPi / 2.0f - dphi)
    );
    uv.col(cur_vertex_id - 1) = Eigen::Vector2f(u_step * i, 1.0f - v_step);
    face.col(cur_face_id++) =
      Eigen::Vector3i(i, cur_vertex_id - 1, cur_vertex_id);
  }
  vertex.col(cur_vertex_id++) = vertex.col(slices);
  uv.col(cur_vertex_id - 1) = Eigen::Vector2f(1.0f, 1.0f - v_step);

  // The middle layers.
  for (int i = 0; i < stacks - 2; ++i) {
    const float phi = kPi / 2.0f - dphi * (i + 2);
    const int first_vertex_id = cur_vertex_id;
    const float v_coord = (stacks - i - 2) * v_step;
    for (int j = 0; j < slices; ++j) {
      vertex.col(cur_vertex_id++) = Eigen::Vector3f(
        cos(phi) * cos(j * dtheta),
        cos(phi) * sin(j * dtheta),
        sin(phi)
      );
      uv.col(cur_vertex_id - 1) = Eigen::Vector2f(u_step * j, v_coord);
      face.col(cur_face_id++) = Eigen::Vector3i(
        cur_vertex_id - 1 - slices - 1,
        cur_vertex_id - 1,
        cur_vertex_id
      );
      face.col(cur_face_id++) = Eigen::Vector3i(
        cur_vertex_id - 1 - slices - 1,
        cur_vertex_id,
        cur_vertex_id - slices - 1
      );
    }
    vertex.col(cur_vertex_id++) = vertex.col(first_vertex_id);
    uv.col(cur_vertex_id - 1) = Eigen::Vector2f(1.0f, v_coord);
  }

  // The bottom-z layer.
  for (int i = 0; i < slices; ++i) {
    face.col(cur_face_id++) = Eigen::Vector3i(
      cur_vertex_id + i - slices - 1,
      cur_vertex_id + i,
      cur_vertex_id + i - slices
    );
  }
  assert(cur_vertex_id + slices == static_cast<int>(vertex.cols()));
  assert(cur_face_id == static_cast<int>(face.cols()));
  // Now scale the sphere.
  vertex *= radius;
}

void Cylinder(const float base, const float top, const float height,
  const int slices, const int stacks, Eigen::Matrix3Xf& vertex,
  Eigen::Matrix3Xi& face) {
  const float dtheta = 2 * kPi / slices;
  vertex = Eigen::Matrix3Xf::Zero(3, slices * (stacks + 1) + 2);
  vertex.col(0) = Eigen::Vector3f(0.0f, 0.0f, 0.0f);
  int cur_vertex_id = 1;
  // Add all vertices first.
  for (int i = 0; i < stacks + 1; ++i) {
    const float s = i * 1.0f / stacks;
    const float radius = (1.0f - s) * base + s * top;
    for (int j = 0; j < slices; ++j) {
      vertex.col(cur_vertex_id) = Eigen::Vector3f(
        cos(j * dtheta) * radius,
        sin(j * dtheta) * radius,
        height * s
      );
      ++cur_vertex_id;
    }
  }
  vertex.col(cur_vertex_id) = Eigen::Vector3f(0.0f, 0.0f, height);
  // Next add all faces.
  std::vector<Eigen::Vector3i> face_list(0);
  // Base.
  if (base > 0.0f) {
    for (int i = 0; i < slices; ++i) {
      face_list.push_back(Eigen::Vector3i(0,
        (i == slices - 1) ? 1 : (i + 2), i + 1));
    }
  }
  // Body.
  for (int i = 0; i < stacks; ++i) {
    for (int j = 0; j < slices; ++j) {
      const int i00 = 1 + i * slices + j;
      const int i10 = (j == slices - 1) ? (1 + i * slices) : i00 + 1;
      const int i01 = i00 + slices;
      const int i11 = i10 + slices;
      face_list.push_back(Eigen::Vector3i(i00, i10, i11));
      face_list.push_back(Eigen::Vector3i(i00, i11, i01));
    }
  }
  // Top.
  if (top > 0.0f) {
    for (int i = 0; i < slices; ++i) {
      face_list.push_back(Eigen::Vector3i(cur_vertex_id,
        1 + stacks * slices + i,
        1 + stacks * slices + ((i == slices - 1) ? 0 : (i + 1))
      ));
    }
  }
  const int face_num = static_cast<int>(face_list.size());
  face = Eigen::Matrix3Xi::Zero(3, face_num);
  for (int i = 0; i < face_num; ++i) {
    face.col(i) = face_list[i];
  }
}

void Cone(const float base, const float height,
  const int slices, const int stacks, Eigen::Matrix3Xf& vertex,
  Eigen::Matrix3Xi& face) {
  Cylinder(base, 0.0f, height, slices, stacks, vertex, face);
}

void Cube(const float size, Eigen::Matrix3Xf& vertex, Eigen::Matrix3Xi& face) {
  vertex = size / 2.0f * (Eigen::Matrix<float, 3, 8>()
    << -1.0f, 1.0f, -1.0f, 1.0f, -1.0f, 1.0f, -1.0f, 1.0f,
    -1.0f, -1.0f, 1.0f, 1.0f, -1.0f, -1.0f, 1.0f, 1.0f,
    -1.0f, -1.0f, -1.0f, -1.0f, 1.0f, 1.0f, 1.0f, 1.0f).finished();
  face = (Eigen::Matrix<int, 3, 12>()
    << 2, 1, 4, 2, 1, 4, 6, 5, 3, 6, 5, 3,
    1, 2, 2, 4, 4, 1, 5, 6, 6, 3, 3, 5,
    0, 3, 0, 6, 0, 5, 7, 4, 7, 2, 7, 1).finished();
}

const Eigen::Matrix4f Rotate(const Eigen::Quaternionf& q) {
  return (Eigen::Matrix4f() << q.matrix(), Eigen::Vector3f::Zero(),
    0.0f, 0.0f, 0.0f, 1.0f).finished();
}

const Eigen::Matrix4f Rotate(const float angle, const Eigen::Vector3f& axis) {
  return (Eigen::Matrix4f() <<
    Eigen::AngleAxisf(DegreeToRadian(angle), axis.normalized()).matrix(),
    Eigen::Vector3f::Zero(), 0.0f, 0.0f, 0.0f, 1.0f).finished();
}

const Eigen::Matrix4f Rotate(const float angle,
  const float x, const float y, const float z) {
  return Rotate(angle, Eigen::Vector3f(x, y, z));
}

const Eigen::Matrix4f Scale(const Eigen::Vector3f& scale) {
  return (Eigen::Vector4f() << scale, 1.0f).finished().asDiagonal();
}

const Eigen::Matrix4f Scale(const float scale) {
  return Scale(Eigen::Vector3f(scale, scale, scale));
}

const Eigen::Matrix4f Scale(const float x, const float y, const float z) {
  return Scale(Eigen::Vector3f(x, y, z));
}

const Eigen::Matrix4f Translate(const Eigen::Vector3f& translate) {
  return (Eigen::Matrix4f() << Eigen::Matrix3f::Identity(), translate,
    0.0f, 0.0f, 0.0f, 1.0f).finished();
}

const Eigen::Matrix4f Translate(const float x, const float y, const float z) {
  return Translate(Eigen::Vector3f(x, y, z));
}

void Rotate(const Eigen::Quaternionf& q, Eigen::Matrix3Xf& vertex) {
  vertex = q.matrix() * vertex;
}

void Rotate(const float angle, const Eigen::Vector3f& axis,
  Eigen::Matrix3Xf& vertex) {
  vertex = Eigen::AngleAxisf(DegreeToRadian(angle), axis.normalized()).matrix()
    * vertex;
}

void Rotate(const float angle, const float x, const float y, const float z,
  Eigen::Matrix3Xf& vertex) {
  Rotate(angle, Eigen::Vector3f(x, y, z), vertex);
}

void Scale(const Eigen::Vector3f& scale, Eigen::Matrix3Xf& vertex) {
  vertex = scale.asDiagonal() * vertex;
}

void Scale(const float scale, Eigen::Matrix3Xf& vertex) {
  Scale(Eigen::Vector3f(scale, scale, scale), vertex);
}

void Scale(const float x, const float y, const float z,
  Eigen::Matrix3Xf& vertex) {
  Scale(Eigen::Vector3f(x, y, z), vertex);
}

void Translate(const Eigen::Vector3f& translate, Eigen::Matrix3Xf& vertex) {
  vertex = vertex.colwise() + translate;
}

void Translate(const float x, const float y, const float z,
  Eigen::Matrix3Xf& vertex) {
  Translate(Eigen::Vector3f(x, y, z), vertex);
}

void DecomposeTransform(const Eigen::Matrix4f& transform,
  Eigen::Vector3f& translation, Eigen::Matrix3f& rotation_u,
  Eigen::Vector3f& scale, Eigen::Matrix3f& rotation_v) {
  assert(transform(3, 0) == 0.0);
  assert(transform(3, 1) == 0.0);
  assert(transform(3, 2) == 0.0);
  assert(transform(3, 3) == 1.0);
  // Determine translation.
  translation = transform.col(3).head(3);
  // SVD.
  const Eigen::JacobiSVD<Eigen::Matrix3f> svd(transform.topLeftCorner(3, 3),
    Eigen::ComputeFullU | Eigen::ComputeFullV);
  rotation_u = svd.matrixU();
  scale = svd.singularValues();
  rotation_v = svd.matrixV().transpose();
  if (rotation_u.determinant() < 0.0) {
    rotation_u.col(2) *= -1;
    rotation_v.col(2) *= -1;
  }
  assert(rotation_u.determinant() > 0.0 && rotation_v.determinant() > 0.0);
}

const Eigen::Matrix3Xf HomogeneousTransform(const Eigen::Matrix4f& transform,
  const Eigen::Matrix3Xf& vertex) {
  const Eigen::Vector4f translation = transform.col(3);
  const Eigen::Matrix4Xf new_homo_vertex =
    (transform.leftCols(3) * vertex).colwise() + translation;
  return new_homo_vertex.topRows(3).cwiseQuotient(
    new_homo_vertex.row(3).replicate(3, 1));
}

void MergeMeshes(const std::vector<Eigen::Matrix3Xf>& vertices,
  const std::vector<Eigen::Matrix3Xi>& faces,
  Eigen::Matrix3Xf& vertex, Eigen::Matrix3Xi& face) {
  assert(vertices.size() == faces.size());
  const int mesh_num = static_cast<int>(vertices.size());
  assert(mesh_num > 0);

  // First, collect the number of vertices of each mesh.
  Eigen::VectorXi cum_vertex_num = Eigen::VectorXi::Zero(mesh_num),
    vertex_num = Eigen::VectorXi::Zero(mesh_num),
    cum_face_num = Eigen::VectorXi::Zero(mesh_num),
    face_num = Eigen::VectorXi::Zero(mesh_num);
  for (int i = 0; i < mesh_num; ++i) {
    vertex_num(i) = static_cast<int>(vertices[i].cols());
    face_num(i) = static_cast<int>(faces[i].cols());
    cum_vertex_num(i) = vertex_num(i) + (i == 0 ? 0 : cum_vertex_num(i - 1));
    cum_face_num(i) = face_num(i) + (i == 0 ? 0 : cum_face_num(i - 1));
  }
  // Now do the merge.
  vertex = Eigen::Matrix3Xf::Zero(3, cum_vertex_num(mesh_num - 1));
  face = Eigen::Matrix3Xi::Zero(3, cum_face_num(mesh_num - 1));
  for (int i = 0; i < mesh_num; ++i) {
    vertex.middleCols(i == 0 ? 0 : cum_vertex_num(i - 1), vertex_num(i))
      = vertices[i];
    face.middleCols(i == 0 ? 0 : cum_face_num(i - 1), face_num(i))
      = faces[i].array() + (i == 0 ? 0 : cum_vertex_num(i - 1));
  }
}

void FlattenMesh(Eigen::Matrix3Xf& vertex, Eigen::Matrix3Xi& face,
  Eigen::Matrix2Xf& uv) {
  assert(uv.size() == 0 || uv.cols() == vertex.cols());
  const int face_num = static_cast<int>(face.cols());
  Eigen::Matrix3Xf flatten_vertex = Eigen::Matrix3Xf::Zero(3, 3 * face_num);
  Eigen::Matrix3Xi flatten_face = Eigen::Matrix3Xi::Zero(3, face_num);

  // Loop over all triangles.
  for (int i = 0; i < face_num; ++i) {
    // Update face.
    flatten_face(0, i) = 3 * i;
    flatten_face(1, i) = 3 * i + 1;
    flatten_face(2, i) = 3 * i + 2;
    // Update vertex.
    const int i0 = face(0, i), i1 = face(1, i), i2 = face(2, i);
    flatten_vertex.col(3 * i) = vertex.col(i0);
    flatten_vertex.col(3 * i + 1) = vertex.col(i1);
    flatten_vertex.col(3 * i + 2) = vertex.col(i2);
  }

  // Update uv if necessary.
  if (uv.size()) {
    Eigen::Matrix2Xf flatten_uv = Eigen::Matrix2Xf::Zero(2, 3 * face_num);
    for (int i = 0; i < face_num; ++i) {
      const int i0 = face(0, i), i1 = face(1, i), i2 = face(2, i);
      flatten_uv.col(3 * i) = uv.col(i0);
      flatten_uv.col(3 * i + 1) = uv.col(i1);
      flatten_uv.col(3 * i + 2) = uv.col(i2);
    }
    uv = flatten_uv;
  }

  // Write back.
  vertex = flatten_vertex;
  face = flatten_face;
}

const Eigen::Matrix3Xf SmoothVertexNormal(const Eigen::Matrix3Xf& vertex,
  const Eigen::Matrix3Xi& face) {
  const int vertex_num = static_cast<int>(vertex.cols());
  const int face_num = static_cast<int>(face.cols());
  Eigen::Matrix3Xf normal = Eigen::Matrix3Xf::Zero(3, vertex_num);
  // Loop over all triangles.
  for (int i = 0; i < face_num; ++i) {
    const int i0 = face(0, i), i1 = face(1, i), i2 = face(2, i);
    const Eigen::Vector3f v0 = vertex.col(i0), v1 = vertex.col(i1),
      v2 = vertex.col(i2);
    const Eigen::Vector3f weighted_normal = (v1 - v0).cross(v2 - v1);
    normal.col(i0) += weighted_normal;
    normal.col(i1) += weighted_normal;
    normal.col(i2) += weighted_normal;
  }
  // Normalization.
  for (int i = 0; i < vertex_num; ++i) {
    Eigen::Vector3f weighted_normal = normal.col(i);
    const float length = weighted_normal.norm();
    if (length < 1e-10f) weighted_normal.setZero();
    else weighted_normal /= length;
    normal.col(i) = weighted_normal;
  }
  return normal;
}

void PointToBoundingBoxFrustum(const Eigen::Vector3f& position,
  const BoundingBox& box, float& near_plane, float& far_plane) {
  // This turns out to be rather tricky. Updating far_plane is easier so let's
  // do it first:
  // Note that cubes are convex, so a frustum defined by far_plane bounds the
  // bounding box as long as all its 8 corners are inside the frustum.
  BoundingBox new_box = box;
  new_box.Translate(-position);
  const Eigen::Vector3f min_corner = new_box.min_corner();
  const Eigen::Vector3f box_size = new_box.box_size();
  const Eigen::Vector3f max_corner = min_corner + box_size;
  far_plane = std::max(far_plane, std::max(min_corner.cwiseAbs().maxCoeff(),
    max_corner.cwiseAbs().maxCoeff()));
  // If near_plane is already zero, then there is no need to update.
  if (near_plane == 0.0f) return;
  // Updating the near_plane is harder. First let's rule out the case where the
  // origin is inside the bounding box.
  if (min_corner.maxCoeff() <= 0.0f && max_corner.minCoeff() >= 0.0f) {
    near_plane = 0.0f;
    return;
  }
  // Now we are certain the origin is strictly outside the bounding box.
  const Eigen::Matrix<float, 6, 1> candidate = (Eigen::Matrix<float, 6, 1>()
    << min_corner.cwiseAbs(), max_corner.cwiseAbs()).finished();
  float candidate_near_plane = candidate.minCoeff();
  for (int i = 0; i < 6; ++i) {
    // Check if two cuboid overlap.
    // Cuboid A: new_box.
    // Cube B: a bounding box centered at (0, 0, 0) with size = 2 * s.
    const float s = candidate(i);
    const bool overlap = min_corner.maxCoeff() < s
      && max_corner.minCoeff() > -s;
    if (!overlap && s > candidate_near_plane) {
      candidate_near_plane = s;
    }
  }
  near_plane = (near_plane < 0.0f ? candidate_near_plane
    : std::min(near_plane, candidate_near_plane));
}

void ReadFromObjFile(const std::string& file_name, Eigen::Matrix3Xf& vertex,
  Eigen::Matrix3Xi& face, Eigen::Matrix2Xf& uv) {
  std::ifstream obj_file(file_name);
  std::string line = "";
  std::vector<float> all_vertices(0);
  std::vector<int> all_faces(0);
  std::vector<float> all_textures(0);
  int vertex_num = 0, face_num = 0, texture_num = 0;
  const std::string whitespace = "\r\n\t ";
  while (std::getline(obj_file, line)) {
    // Trim whitespace.
    size_t pos = line.find_first_not_of(whitespace);
    if (pos != std::string::npos)
      line = line.substr(pos);

    // Extract vertex, face and texture information.
    if (line.substr(0, 2) == "v ") {
      std::stringstream vertex_str(line.substr(2));
      float x, y, z;
      vertex_str >> x >> y >> z;
      all_vertices.push_back(x);
      all_vertices.push_back(y);
      all_vertices.push_back(z);
      ++vertex_num;
    } else if (line.substr(0, 3) == "vt ") {
      std::stringstream vt_str(line.substr(3));
      float u, v;
      vt_str >> u >> v;
      all_textures.push_back(u); all_textures.push_back(v);
      ++texture_num;
    } else if (line.substr(0, 2) == "f ") {
      // All possible strings:
      // - f v1 v2 v3
      // - f v1/vt1 v2/vt2 ...
      // - f v1/vt1/vn1 v2/vt2/vn2 ...
      // - f v1//vn1 v2//vn2 ...
      line = line.substr(2);
      std::vector<int> id_in_face(0);
      while (!line.empty()) {
        // Trim whitespace.
        pos = line.find_first_not_of(whitespace);
        if (pos == std::string::npos) break;
        line = line.substr(pos);
        // Read vertex index.
        std::stringstream f_str(line);
        int v_id; f_str >> v_id;
        id_in_face.push_back(v_id);
        // Proceed.
        pos = line.find_first_of(whitespace);
        if (pos == std::string::npos) break;
        line = line.substr(pos);
      }
      assert(id_in_face.size() >= 3u);
      const int id_num_in_face = static_cast<int>(id_in_face.size());
      for (int i = 1; i < id_num_in_face - 1; ++i) {
        all_faces.push_back(id_in_face[0]);
        all_faces.push_back(id_in_face[i]);
        all_faces.push_back(id_in_face[i + 1]);
        ++face_num;
      }
    }
  }
  assert(!texture_num || vertex_num == texture_num);
  vertex = Eigen::Map<Eigen::Matrix3Xf>(all_vertices.data(), 3, vertex_num);
  face = Eigen::Map<Eigen::Matrix3Xi>(all_faces.data(), 3, face_num);
  assert(face.maxCoeff() <= vertex_num);
  face = face.array() - 1;
  uv = Eigen::Map<Eigen::Matrix2Xf>(all_textures.data(), 2, texture_num);
}

void WriteToObjFile(const Eigen::Matrix3Xf& vertex,
  const Eigen::Matrix3Xi& face, const std::string& file_name) {
  std::ofstream obj_file;
  obj_file.open(file_name);
  const int vertex_num = static_cast<int>(vertex.cols());
  for (int i = 0; i < vertex_num; ++i) {
    obj_file << "v " << vertex(0, i) << " "
      << vertex(1, i) << " " << vertex(2, i) << std::endl;
  }
  const int face_num = static_cast<int>(face.cols());
  for (int i = 0; i < face_num; ++i) {
    obj_file << "f " << face(0, i) + 1 << " "
      << face(1, i) + 1 << " " << face(2, i) + 1 << std::endl;
  }
  obj_file.close();
}

}