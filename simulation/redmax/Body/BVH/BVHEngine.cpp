/* This file is part of PyMesh. Copyright (c) 2018 by Qingnan Zhou */

#include "Body/BVH/BVHEngine.h"
#include "Body/BVH/AABBTree.h"

using namespace redmax;

BVHEngine::Ptr BVHEngine::create(size_t dim) {
    switch (dim) {
        case 2:
            return std::make_shared<IGL::AABBTree<2>>();
        case 3:
            return std::make_shared<IGL::AABBTree<3>>();
        default:
            throw NotImplementedError("Only 2D and 3D meshes are supported");
    }
}
