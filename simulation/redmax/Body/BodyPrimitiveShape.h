#pragma once
#include "Body/Body.h"

namespace redmax {

class BodyPrimitiveShape : public Body {
public:
    BodyPrimitiveShape(Simulation* sim, Joint* joint, dtype density)
        :Body(sim, joint, density) {}

    BodyPrimitiveShape(Simulation* sim, Joint* joint, Matrix3 R_ji, Vector3 p_ji, dtype density)
        :Body(sim, joint, R_ji, p_ji, density) {}

    /**
     * analytical distance function
     * @param:
     * xw: the location of the query position in world frame
     * @return: contact distance of xw
    **/
    virtual dtype distance(Vector3 xw) = 0;

    /**
     * analytical distance function return dist, normal, x2i, ddot, tdot
     * @return
     * dist: the distance of xw
     * normal: the normalized contact normal
     * xi2: the contact position on this body
     * ddot: the magnitude of the velocity on normal direction
     * tdot: the tangential velocity
     **/
    virtual void collision(Vector3 xw, Vector3 xw_dot, /* input */
                            dtype &d, Vector3 &n,  /* output */
                            dtype &ddot, Vector3 &tdot,
                            Vector3 &xi2) = 0;

    /**
     * analytical distance function return dist, normal, x2i, ddot, tdot and derivatives
     * @return
     * dist: the distance of xw
     * normal: the normalized contact normal
     * xi2: the contact position on this body
     * ddot: the magnitude of the velocity on normal direction
     * tdot: the tangential velocity
     * derivatives
     **/
    virtual void collision(Vector3 xw, Vector3 xw_dot, /* input */
                            dtype &d, Vector3 &n,  /* output */
                            dtype &ddot, Vector3 &tdot,
                            Vector3 &xi2,
                            RowVector3 &dd_dxw, RowVector6 &dd_dq2, /* derivatives for d */ 
                            Matrix3 &dn_dxw, Matrix36 &dn_dq2, /* derivatives for n */
                            RowVector3 &dddot_dxw, RowVector3 &dddot_dxwdot, /* derivatives for ddot */
                            RowVector6 &dddot_dq2, RowVector6 &dddot_dphi2,
                            Matrix3 &dtdot_dxw, Matrix3 &dtdot_dxwdot, /* derivatives for tdot */
                            Matrix36 &dtdot_dq2, Matrix36 &dtdot_dphi2,
                            Matrix3 &dxi2_dxw, Matrix36 &dxi2_dq2 /* derivatives for xi2 */) = 0;

    virtual void test_collision_derivatives() {};
    virtual void test_collision_derivatives_runtime(Vector3 xw, Vector3 xw_dot) {};
};

}