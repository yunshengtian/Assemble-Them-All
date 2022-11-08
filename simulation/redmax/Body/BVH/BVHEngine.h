/* This file is part of PyMesh. Copyright (c) 2018 by Qingnan Zhou */
#pragma once

#include <memory>
#include <string>
#include <vector>

#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_edge_normals.h>

#include "Common.h"
#include "Body/BVH/Exception.h"

namespace redmax {

/**
 * Boundary Volume Hierarchy is commonly used for accelerate spatial queries and
 * intersections tests.  This class unify multiple BVH implementations under the
 * same interface.
 */
class BVHEngine {
    public:
        typedef std::shared_ptr<BVHEngine> Ptr;
        static Ptr create(size_t dim);

    public:
        virtual ~BVHEngine() = default;

    public:
        void set_mesh(const MatrixXr& vertices, const MatrixXir& faces) {
            if (faces.cols() != 3) {
                throw NotImplementedError(
                        "Only Triangle mesh is supported by BVH engine");
            }
            m_vertices = vertices;
            m_faces = faces;
            igl::per_face_normals(m_vertices, m_faces, m_face_normals);
            igl::per_vertex_normals(m_vertices, m_faces, igl::PER_VERTEX_NORMALS_WEIGHTING_TYPE_ANGLE, m_face_normals, m_vertex_normals);
            igl::per_edge_normals(m_vertices, m_faces, igl::PER_EDGE_NORMALS_WEIGHTING_TYPE_UNIFORM, m_face_normals, m_edge_normals, m_edges, m_edge_map);
            // std::cerr << "debug BVHEngine" << std::endl;
            // std::cerr << m_vertices.rows() << " " << m_vertices.cols() << std::endl;
            // std::cerr << m_faces.rows() << " " << m_faces.cols() << std::endl;
            // std::cerr << m_face_normals.rows() << " " << m_face_normals.cols() << std::endl;
            // std::cerr << m_vertex_normals.rows() << " " << m_vertex_normals.cols() << std::endl;
            // std::cerr << m_edge_normals.rows() << " " << m_edge_normals.cols() << std::endl;
            // std::cerr << m_edges.rows() << " " << m_edges.cols() << std::endl;
            // std::cerr << m_edge_map.size() << std::endl;
        }

        /**
         * Build BVH for the given mesh.
         */
        virtual void build() {
            throw NotImplementedError("BVH algorithm is not implemented");
        }

        /**
         * For each point in points, lookup the closest points on mesh and the
         * corresponding distances and faces.
         */
        virtual void lookup(const MatrixXr& points,
                VectorX& squared_distances,
                VectorXi& closest_faces,
                MatrixXr& closest_points) const {
            throw NotImplementedError("BVH algorithm is not implemented");
        }

        /**
         * For each point in points, lookup the closest points on mesh and the
         * corresponding signed (un-squared) distances and faces.
         * Warning: only work with IGL engine
         */
        virtual void lookup_signed(const MatrixXr& points,
                VectorX& signed_distances,
                VectorXi& closest_faces,
                MatrixXr& closest_points,
                MatrixXr& closest_face_normals) const {
            throw NotImplementedError("BVH algorithm is not implemented");
        }

    protected:
        MatrixXr m_vertices;
        MatrixXir m_faces;
        MatrixXr m_face_normals;
        MatrixXr m_vertex_normals;
        MatrixXr m_edge_normals;
        MatrixXir m_edges;
        VectorXi m_edge_map;
};

}