#pragma once
#include "Common.h"
#include "Utils.h"
#include "CollisionDetection/Contact.h"

namespace redmax {

class Body;
class BodyCuboid;
class BodyPrimitiveShape;
class BodySDFObj;
class BodyBVHObj;

// detect the collision between ground and a cuboid body
// return the list of the contact points represented in body frame.
void collision_detection_ground_object(const Matrix4 E_g, const Body* body, std::vector<Contact>& contacts);

void collision_detection_cuboid_cuboid(const BodyCuboid* cuboid1, const BodyCuboid* cuboid2, std::vector<Contact>& contacts);

// detect the collision between a general body and a primitive body
// the general body should be able to give a list of contact points on the surface.
// the primitive body is supposed to have an anlytical distance field
void collision_detection_general_primitive(const Body* contact_body, const BodyPrimitiveShape* primitive_body, std::vector<Contact>& contacts);

// detect the collision between a general body and a SDF body
void collision_detection_general_SDF(const Body* contact_body, const BodySDFObj* SDF_body, std::vector<Contact>& contacts);

// detect the collision between a general body and a BVH body
void collision_detection_general_BVH(const Body* contact_body, const BodyBVHObj* BVH_body, std::vector<Contact>& contacts);

}