#include "Joint/JointFree2D.h"
#include "Utils.h"

namespace redmax {

JointFree2D::JointFree2D(Simulation *sim, int id, Joint* parent, Matrix3 R_pj_0, Vector3 p_pj_0, Joint::Frame frame) 
    : Joint(sim, id, parent, R_pj_0, p_pj_0, 3, frame) {}

/*  notes section 3.2.7
    Q(q) = [R  p]
           [0  1]
    A(q) = [R     0]
           [[p]R  R]
    dR_dq3 = [-sin(q3) -cos(q3) 0]
             [cos(q3)  -sin(q3) 0]
             [0        0        0]
    dA_dq1 = [0      0]
             [[e1]R  0]
    dA_dq2 = [0      0]
             [[e2]R  0]
    dA_dq3 = [dR_dq3          0]
             [[p]dR_dq3  dR_dq3]
    Rdot = [-sin(q3)*qdot3  -cos(q3)*qdot3  0]
           [cos(q3)*qdot3   -sin(q3)*qdot3  0]
           [0               0               0]
    Adot = [Rdot                0]
           [[pdot]R+[p]Rdot  Rdot]
    dRdot_dq3 = [-cos(q3)*qdot3   sin(q3)*qdot3  0]
                [-sin(q3)*qdot3  -cos(q3)*qdot3  0]
                [0                0              0]
    dAdot_dq1 = [0        0]
                [[e1]Rdot 0]
    dAdot_dq2 = [0        0]
                [[e2]Rdot 0]
    dAdot_dq3 = [dRdot_dq3                          0]
                [[pdot]dR_dq3+[p]dRdot_dq3  dRdot_dq3]
    S(Q) = [0        0       0]
           [0        0       0]
           [0        0       1]
           [cos(q3)  sin(q3) 0]
           [-sin(q3) cos(q3) 0]
           [0        0       0]
    dS_dq1 = dS_dq2 = 0
    dS_dq3 = [0         0        0]
             [0         0        0]
             [0         0        0]
             [-sin(q3)  cos(q3)  0]
             [-cos(q3) -sin(q3)  0]
             [0         0        0]
    Sdot = [0               0             0]
           [0               0             0]
           [0               0             0]
           [-qdot3*sin(q3)  qdot3*cos(q3) 0]
           [-qdot3*cos(q3) -qdot3*sin(q3) 0]
           [0               0             0]
    dSdot_dq1 = dSdot_dq2 = 0
    dSdot_dq3 = [0               0             0]
                [0               0             0]
                [0               0             0]
                [-qdot3*cos(q3) -qdot3*sin(q3) 0]
                [qdot3*sin(q3)  -qdot3*cos(q3) 0]
                [0               0             0]
    dQ_dq1 = [0 ex]
             [0  0]
    dQ_dq2 = [0 ey]
             [0  0]
    dQ_dq3 = [dR_dq3 0]
             [0      0]
*/
void JointFree2D::update(bool design_gradient) {
    Vector3 p(_q(0), _q(1), 0);
    dtype r = _q(2);
    Vector3 pdot(_qdot(0), _qdot(1), 0);
    dtype rdot = _qdot(2);
    dtype s = sin(r), c = cos(r);
    Matrix3 pbrac = math::skew(p);
    Matrix3 pdotbrac = math::skew(pdot);

    Matrix3 R = Matrix3::Identity();
    R.topLeftCorner(2, 2) << c, -s, 
                            s, c;
    
    _Q.setIdentity();
    _Q.topLeftCorner(3, 3) = R;
    _Q.topRightCorner(3, 1) = p;

    _A = math::Ad(_Q);

    _S_j.setZero();
    _S_j(2, 2) = 1;
    _S_j.block(3, 0, 2, 2) << c, s,
                             -s, c;
    
    Matrix3 dR_dq3 = Matrix3::Zero();
    dR_dq3.topLeftCorner(2, 2) << -s, -c,
                                   c, -s;
    
    Matrix3 Rdot = dR_dq3 * rdot;

    _Adot << Rdot, Matrix3::Zero(), pdotbrac * R + pbrac * Rdot, Rdot;
    
    _S_j_dot.setZero();
    _S_j_dot.block(3, 0, 2, 2) << -rdot * s, rdot * c, -rdot * c, -rdot * s;

    // derivatives
    _dSjdot_dq.setZero();
    _dSjdot_dq(2).block(3, 0, 2, 2) << -rdot * c, -rdot * s, rdot * s, -rdot * c;

    Matrix3 e1brac = math::skew(Vector3::UnitX());
    Matrix3 e2brac = math::skew(Vector3::UnitY());
    
    _dA_dq.setZero();
    _dA_dq(0).bottomLeftCorner(3, 3) = e1brac * R;
    _dA_dq(1).bottomLeftCorner(3, 3) = e2brac * R;
    _dA_dq(2) << dR_dq3, Matrix3::Zero(), 
                 pbrac * dR_dq3, dR_dq3;
    
    _dSj_dq.setZero();
    _dSj_dq(2).block(3, 0, 2, 2) << -s, c,
                                    -c, -s;
    
    Matrix3 dRdot_dq3 = Matrix3::Zero();
    dRdot_dq3.topLeftCorner(2, 2) << -rdot * c, rdot * s,
                                     -rdot * s, -rdot * c;
    
    _dAdot_dq.setZero();
    _dAdot_dq(0).bottomLeftCorner(3, 3) = e1brac * Rdot;
    _dAdot_dq(1).bottomLeftCorner(3, 3) = e2brac * Rdot;
    _dAdot_dq(2) << dRdot_dq3, Matrix3::Zero(),
                    pdotbrac * dR_dq3 + pbrac * dRdot_dq3, dRdot_dq3;
    
    _dQ_dq.setZero();
    _dQ_dq(0)(0, 3) = 1.;
    _dQ_dq(1)(1, 3) = 1.;
    _dQ_dq(2).topLeftCorner(3, 3) = dR_dq3;
    
    Joint::update(design_gradient);
}

}