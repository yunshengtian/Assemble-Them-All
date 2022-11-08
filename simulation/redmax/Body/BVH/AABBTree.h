/* This file is part of PyMesh. Copyright (c) 2018 by Qingnan Zhou */
#pragma once

#include <igl/AABB.h>
#include <igl/signed_distance.h>

#include "Common.h"
#include "Body/BVH/BVHEngine.h"

namespace redmax {
namespace IGL {

template<int DIM>
class AABBTree : public BVHEngine {};

template<>
class AABBTree<2> : public BVHEngine {
    public:
        using Ptr = std::shared_ptr<AABBTree>;
        using Tree = igl::AABB<MatrixXr, 2>;

    public:
        virtual ~AABBTree() = default;

        virtual void build() override {
            m_tree.init(m_vertices, m_faces);
        }

        virtual void lookup(const MatrixXr& points,
                VectorX& squared_distances,
                VectorXi& closest_faces,
                MatrixXr& closest_points) const override {
            m_tree.squared_distance(m_vertices, m_faces, points,
                    squared_distances, closest_faces, closest_points);
        }

        virtual void lookup_signed(const MatrixXr& points,
                VectorX& signed_distances,
                VectorXi& closest_faces,
                MatrixXr& closest_points,
                MatrixXr& closest_face_normals) const override {
            throw NotImplementedError("Signed distance in 2D is not yet supported.");
        }

    private:
        Tree m_tree;
};

template<>
class AABBTree<3> : public BVHEngine {
    public:
        using Ptr = std::shared_ptr<AABBTree>;
        using Tree = igl::AABB<MatrixXr, 3>;

    public:
        virtual ~AABBTree() = default;

        virtual void build() override {
            m_tree.init(m_vertices, m_faces);
        }

        virtual void lookup(const MatrixXr& points,
                VectorX& squared_distances,
                VectorXi& closest_faces,
                MatrixXr& closest_points) const override {
            m_tree.squared_distance(m_vertices, m_faces, points,
                    squared_distances, closest_faces, closest_points);
        }

        virtual void lookup_signed(const MatrixXr& points,
                VectorX& signed_distances,
                VectorXi& closest_faces,
                MatrixXr& closest_points,
                MatrixXr& closest_face_normals) const override {
          igl::signed_distance_pseudonormal(points, m_vertices, m_faces, m_tree,
              m_face_normals, m_vertex_normals, m_edge_normals, m_edge_map,
              signed_distances, closest_faces, closest_points, closest_face_normals);
        }

    private:
        Tree m_tree;
};

}
}
