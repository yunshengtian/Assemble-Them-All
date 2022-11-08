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
#ifndef _OPENGL_VIEWER_GEOMETRY_H_
#define _OPENGL_VIEWER_GEOMETRY_H_

#include <string>
#include <vector>
#include "Eigen/Dense"

namespace opengl_viewer {

const float kPi = 3.14159265358979323846264335f;

class BoundingBox;

inline const float DegreeToRadian(const float degree) {
  return degree / 180.0f * kPi;
}

inline const float RadianToDegree(const float radian) {
  return radian / kPi * 180.0f;
}

const Eigen::Matrix4f LookAt(const Eigen::Vector3f& position,
  const Eigen::Vector3f& lookat, const Eigen::Vector3f& up);

// Field of view is in degrees.
const Eigen::Matrix4f Perspective(const float aspect_ratio,
  const float field_of_view, const float z_min, const float z_max);

// Draw a sphere centered at (0, 0, 0) with given radius. slices is the number
// of subdivisions around the z axis, and stacks is the number of subdivisions
// along the z axis.
void Sphere(const float radius, const int slices, const int stacks,
  Eigen::Matrix3Xf& vertex, Eigen::Matrix3Xi& face, Eigen::Matrix2Xf& uv);

// https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/gluCylinder.xml
void Cylinder(const float base, const float top, const float height,
  const int slices, const int stacks, Eigen::Matrix3Xf& vertex,
  Eigen::Matrix3Xi& face);

void Cone(const float base, const float height, const int slices,
  const int stacks, Eigen::Matrix3Xf& vertex, Eigen::Matrix3Xi& face);

void Cube(const float size, Eigen::Matrix3Xf& vertex, Eigen::Matrix3Xi& face);

const Eigen::Matrix4f Rotate(const Eigen::Quaternionf& q);
// angle is in degrees.
const Eigen::Matrix4f Rotate(const float angle, const Eigen::Vector3f& axis);
const Eigen::Matrix4f Rotate(const float angle,
  const float x, const float y, const float z);
const Eigen::Matrix4f Scale(const Eigen::Vector3f& scale);
const Eigen::Matrix4f Scale(const float scale);
const Eigen::Matrix4f Scale(const float x, const float y, const float z);
const Eigen::Matrix4f Translate(const Eigen::Vector3f& translate);
const Eigen::Matrix4f Translate(const float x, const float y, const float z);

void Rotate(const Eigen::Quaternionf& q, Eigen::Matrix3Xf& vertex);
void Rotate(const float angle, const Eigen::Vector3f& axis,
  Eigen::Matrix3Xf& vertex);
void Rotate(const float angle, const float x, const float y, const float z,
  Eigen::Matrix3Xf& vertex);
void Scale(const Eigen::Vector3f& scale, Eigen::Matrix3Xf& vertex);
void Scale(const float scale, Eigen::Matrix3Xf& vertex);
void Scale(const float x, const float y, const float z,
  Eigen::Matrix3Xf& vertex);
void Translate(const Eigen::Vector3f& translate, Eigen::Matrix3Xf& vertex);
void Translate(const float x, const float y, const float z,
  Eigen::Matrix3Xf& vertex);
// Each valid transform (a combination of Translation, Rotate and Scale) can be
// decomposed into the following form:
//   [U\Sigma V, t; 0 0 0 1].
// where U\Sigma V is the SVD of the top left 3 x 3 block.
void DecomposeTransform(const Eigen::Matrix4f& transform,
  Eigen::Vector3f& translation, Eigen::Matrix3f& rotation_u,
  Eigen::Vector3f& scale, Eigen::Matrix3f& rotation_v);

const Eigen::Matrix3Xf HomogeneousTransform(const Eigen::Matrix4f& transform,
  const Eigen::Matrix3Xf& vertex);

// Mesh-related functions.
// Provide an empty uv if you want to ignore it.
void MergeMeshes(const std::vector<Eigen::Matrix3Xf>& vertices,
  const std::vector<Eigen::Matrix3Xi>& faces,
  Eigen::Matrix3Xf& vertex, Eigen::Matrix3Xi& face);
void FlattenMesh(Eigen::Matrix3Xf& vertex, Eigen::Matrix3Xi& face,
  Eigen::Matrix2Xf& uv);
// Given vertices and faces, compute a normal vector for each vertex. The
// normal is the area-weighted average of all normals from adjacent triangles.
const Eigen::Matrix3Xf SmoothVertexNormal(const Eigen::Matrix3Xf& vertex,
  const Eigen::Matrix3Xi& face);

// Given a point and a bounding box, modify the near_plane and far_plane so
// that the bounding box *AND* the old frustum is fully visible in the new
// frustum. For the first call you should pass a *strictly* negative near_plane
// and far_plane = 0.0f, so that I know I should use the bounding box to
// initialize the frustum.
void PointToBoundingBoxFrustum(const Eigen::Vector3f& position,
  const BoundingBox& box, float& near_plane, float& far_plane);

// Read *one single object* from the obj file.
void ReadFromObjFile(const std::string& file_name, Eigen::Matrix3Xf& vertex,
  Eigen::Matrix3Xi& face, Eigen::Matrix2Xf& uv);

void WriteToObjFile(const Eigen::Matrix3Xf& vertex,
  const Eigen::Matrix3Xi& face, const std::string& file_name);

}

#endif