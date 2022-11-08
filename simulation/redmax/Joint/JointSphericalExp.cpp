#include "Joint/JointSphericalExp.h"
#include "Utils.h"

namespace redmax {

JointSphericalExp::JointSphericalExp(
    Simulation *sim, int id, Joint * parent, Matrix3 R_pj_0, Vector3 p_pj_0, Joint::Frame frame) 
    : Joint(sim, id, parent, R_pj_0, p_pj_0, 3, frame) {
}

/* follow overleaf notes
*/
void JointSphericalExp::update(bool design_gradient) {
    Vector3 q = _q;
    Matrix3 R = math::exp(q);
    Matrix3 I = Matrix3::Identity();
    Vector3 e[3];
    for (int j = 0;j < 3;j++)
        e[j] = Vector3::Unit(j);

    JacobianMatrixVector dR_dq(3, 3, 3); 
    Matrix3 Rdot;

    if (q.norm() < math::eps) {     // case 1: |q| < eps
        R = I;

        _Q.setIdentity();
        _Q.topLeftCorner(3, 3) = R;

        dR_dq.setZero();
        Rdot.setZero();
        for (int j = 0;j < 3;j++) {
            dR_dq(j) = math::skew(e[j]);
            Rdot += math::skew(e[j]) * _qdot(j);
        }

        _S_j.setZero(); _S_j.topRows(3).setIdentity();
        // _S_j_dot.setZero();
        for (int j = 0;j < 3;j++) {
            _S_j_dot.col(j).head(3) = math::unskew(Rdot.transpose() * math::skew(e[j]));
        }

        _A.setZero();
        _A.topLeftCorner(3, 3) = R; _A.bottomRightCorner(3, 3) = R;

        _Adot.setZero();
        _Adot.topLeftCorner(3, 3) = Rdot; _Adot.bottomRightCorner(3, 3) = Rdot;

        _dSj_dq.setZero();
        for (int j = 0;j < 3;j++)
            for (int k = 0;k < 3;k++)
                _dSj_dq(k).col(j).head(3) = math::unskew(math::skew(e[k]).transpose() * math::skew(e[j]));

        _dSjdot_dq.setZero();

        _dA_dq.setZero();
        for (int k = 0;k < 3;k++) {
            _dA_dq(k).topLeftCorner(3, 3) = dR_dq(k);
            _dA_dq(k).bottomRightCorner(3, 3) = dR_dq(k);
        }

        _dQ_dq.setZero();
        for (int k = 0;k < 3;k++)
            _dQ_dq(k).topLeftCorner(3, 3) = dR_dq(k);

    } else {                        // case 2: |q| >= eps
        _Q.setIdentity();
        _Q.topLeftCorner(3, 3) = R;

        dtype q2norm = q.dot(q);

        // compute components
        dtype d, ddot;
        Matrix3 A[3], B[3], C[3], Adot[3], Bdot[3], Cdot[3];
        RowVector3 dd_dq, dddot_dq;
        JacobianMatrixVector dA_dq[3], dB_dq[3], dC_dq[3], dAdot_dq[3], dBdot_dq[3], dCdot_dq[3];
        for (int i = 0;i < 3;i++) {
            dA_dq[i] = JacobianMatrixVector(3, 3, 3);
            dB_dq[i] = JacobianMatrixVector(3, 3, 3);
            dC_dq[i] = JacobianMatrixVector(3, 3, 3);
            dAdot_dq[i] = JacobianMatrixVector(3, 3, 3);
            dBdot_dq[i] = JacobianMatrixVector(3, 3, 3);
            dCdot_dq[i] = JacobianMatrixVector(3, 3, 3);
        }
        
        // compute d, A, B, C
        d = 1. / q2norm;
        for (int j = 0;j < 3;j++) {
            C[j] = math::skew(math::skew(q) * ((I - R) * e[j]));
            B[j] = q(j) * math::skew(q);
            A[j] = (B[j] + C[j]) * d;
        }

        // compute dR_dq and Rdot
        Rdot.setZero();
        for (int j = 0;j < 3;j++) {
            dR_dq(j) = A[j] * R;
            Rdot += dR_dq(j) * _qdot(j);
        }

        // compute ddot, Adot, Bdot, Cdot
        ddot = -2 * q.dot(_qdot) / (q2norm * q2norm);
        for (int j = 0;j < 3;j++) {
            Cdot[j] = math::skew(math::skew(_qdot) * ((I - R) * e[j]) - math::skew(q) * (Rdot * e[j]));
            Bdot[j] = _qdot(j) * math::skew(q) + q(j) * math::skew(_qdot);
            Adot[j] = (Bdot[j] + Cdot[j]) * d + (B[j] + C[j]) * ddot;
        }

        // compute components' derivatives
        for (int k = 0;k < 3;k++) {
            dd_dq(k) = -2 * q(k) / (q2norm * q2norm);

            for (int j = 0;j < 3;j++) {
                dB_dq[j](k) = q(j) * math::skew(e[k]);
                if (j == k) {
                    dB_dq[j](k) += math::skew(q);
                }

                dC_dq[j](k) = math::skew(math::skew(e[k]) * ((I - R) * e[j])) + math::skew(-math::skew(q) * (dR_dq(k) * e[j]));

                dA_dq[j](k) = (dB_dq[j](k) + dC_dq[j](k)) * d + (B[j] + C[j]) * dd_dq(k);
            }
        }

        // compute dRdot_dq
        JacobianMatrixVector dRdot_dq(3, 3, 3); dRdot_dq.setZero();
        for (int k = 0;k < 3;k++) {
            for (int j = 0;j < 3;j++) {
                Matrix3 d2R_dqjdqk = dA_dq[j](k) * R + A[j] * dR_dq(k);
                dRdot_dq(k) += d2R_dqjdqk * _qdot[j];
            }
        }

        // compute components' dot's derivatives
        for (int k = 0;k < 3;k++) {
            dddot_dq(k) = (-2. * _qdot(k) * q2norm + 8. * q.dot(_qdot) * q(k)) / pow(q2norm, 3);

            for (int j = 0;j < 3;j++) {
                dBdot_dq[j](k) = _qdot(j) * math::skew(e[k]);
                if (j == k) 
                    dBdot_dq[j](k) += math::skew(_qdot);
                
                dCdot_dq[j](k) = math::skew(-math::skew(_qdot) * (dR_dq(k) * e[j]) - math::skew(e[k]) * (Rdot * e[j]) - math::skew(q) * (dRdot_dq(k) * e[j]));

                dAdot_dq[j](k) = (dBdot_dq[j](k) + dCdot_dq[j](k)) * d + (Bdot[j] + Cdot[j]) * dd_dq(k) + (dB_dq[j](k) + dC_dq[j](k)) * ddot + (B[j] + C[j]) * dddot_dq(k);
            }
        }

        // compute our required values and their derivatives
        _S_j.setZero();
        for (int j = 0;j < 3;j++) {
            _S_j.col(j).head(3) = math::unskew(R.transpose() * dR_dq(j));
            _S_j_dot.col(j).head(3) = math::unskew(Rdot.transpose() * A[j] * R + R.transpose() * Adot[j] * R + R.transpose() * A[j] * Rdot);

            _A.setZero(); 
            _A.topLeftCorner(3, 3) = R; _A.bottomRightCorner(3, 3) = R;

            _Adot.setZero();
            _Adot.topLeftCorner(3, 3) = Rdot; _Adot.bottomRightCorner(3, 3) = Rdot;

            // derivatives
            for (int k = 0;k < 3;k++) {
                Matrix3 d2R_dqjdqk = dA_dq[j](k) * R + A[j] * dR_dq(k);

                _dSj_dq(k).col(j).head(3) = math::unskew(dR_dq(k).transpose() * dR_dq(j) + R.transpose() * d2R_dqjdqk);

                _dSjdot_dq(k).col(j).head(3) = math::unskew(
                    dRdot_dq(k).transpose() * A[j] * R + Rdot.transpose() * dA_dq[j](k) * R + Rdot.transpose() * A[j] * dR_dq(k) 
                    + dR_dq(k).transpose() * Adot[j] * R + R.transpose() * dAdot_dq[j](k) * R + R.transpose() * Adot[j] * dR_dq(k)
                    + dR_dq(k).transpose() * A[j] * Rdot + R.transpose() * dA_dq[j](k) * Rdot + R.transpose() * A[j] * dRdot_dq(k));
            }
        }

        _dQ_dq.setZero();
        _dA_dq.setZero();
        _dAdot_dq.setZero();
        for (int k = 0;k < 3;k++) {
            _dA_dq(k).topLeftCorner(3, 3) = dR_dq(k); _dA_dq(k).bottomRightCorner(3, 3) = dR_dq(k);

            _dAdot_dq(k).topLeftCorner(3, 3) = dRdot_dq(k); _dAdot_dq(k).bottomRightCorner(3, 3) = dRdot_dq(k);

            _dQ_dq(k).topLeftCorner(3, 3) = dR_dq(k);
        }
    }

    Joint::update(design_gradient);
}

void JointSphericalExp::inner_update() {
    Vector3 q = _q;
    Matrix3 R = math::exp(q);
    Matrix3 I = Matrix3::Identity();
    Vector3 e[3];
    for (int j = 0;j < 3;j++)
        e[j] = Vector3::Unit(j);

    JacobianMatrixVector dR_dq(3, 3, 3); 

    if (q.norm() < math::eps) {     // case 1: |q| < eps
        R = I;

        _Q.setIdentity();
        _Q.topLeftCorner(3, 3) = R;

        _S_j.setZero(); _S_j.topRows(3).setIdentity();
    } else {                        // case 2: |q| >= eps
        _Q.setIdentity();
        _Q.topLeftCorner(3, 3) = R;

        dtype q2norm = q.dot(q);

        // compute components
        dtype d;
        Matrix3 A[3], B[3], C[3];
        
        // compute d, A, B, C
        d = 1. / q2norm;
        for (int j = 0;j < 3;j++) {
            C[j] = math::skew(math::skew(q) * ((I - R) * e[j]));
            B[j] = q(j) * math::skew(q);
            A[j] = (B[j] + C[j]) * d;
        }

        // compute dR_dq
        for (int j = 0;j < 3;j++) {
            dR_dq(j) = A[j] * R;
        }

        // compute our required values and their derivatives
        _S_j.setZero();
        for (int j = 0;j < 3;j++) {
            _S_j.col(j).head(3) = math::unskew(R.transpose() * dR_dq(j));
        }
    }
}

bool JointSphericalExp::reparam() {
    // if q.norm is close to 2pi, then reparameterize as q = (1 - 2pi / q.norm) * q
    dtype qnorm = _q.norm();
    if (qnorm > (dtype)1.5 * constants::pi) {
        Matrix3 S0 = _S_j.topLeftCorner(3, 3);
        // std::cerr << "q0 = " << _q.transpose() << std::endl;
        // std::cerr << "phi0 = " << (S0 * _qdot).transpose() << std::endl;
        // std::cerr << "S = " << _S_j << std::endl;
        // std::cerr << "Q0 = " << _Q << std::endl;

        while (qnorm > (dtype)1.5 * constants::pi) {
            _q = ((dtype)1.0 - 2.0 * constants::pi / qnorm) * _q;
            qnorm = _q.norm();
        }
        
        inner_update();
        Matrix3 S1 = _S_j.topLeftCorner(3, 3);
        _qdot = S1.inverse() * S0 * _qdot;
        // std::cerr << "q1 = " << _q.transpose() << std::endl;
        // std::cerr << "phi1 = " << (S1 * _qdot).transpose() << std::endl;
        // std::cerr << "S = " << _S_j << std::endl;
        // std::cerr << "Q1 = " << _Q << std::endl;
        // inner_update();
        return true;
    } else {
        return false;
    }
}

}