#include "Force/ForceSpringMultiPointGeneric.h"
#include "Body/Body.h"
#include "Simulation.h"
#include "Joint/Joint.h"

namespace redmax {

ForceSpringMultiPointGeneric::ForceSpringMultiPointGeneric(Simulation* sim) : Force(sim) {
    _bodies.clear();
    _xls.clear();
}

void ForceSpringMultiPointGeneric::addBodyPoint(Body* body, Vector3 xl) {
    _bodies.push_back(body);
    _xls.push_back(xl);
}

void ForceSpringMultiPointGeneric::computeForce(VectorX& fm, VectorX& fr, bool verbose) {
    int n = _bodies.size();
    std::vector<Matrix4> E;
    std::vector<MatrixX> gamma, J;
    std::vector<Vector6> phi;
    std::vector<Matrix3> R;
    std::vector<Vector3> p;
    std::vector<Vector3> xw, vw;
    for (int k = 0;k < n;k++) {
        E.push_back(_bodies[k]->_E_0i);
        phi.push_back(_bodies[k]->_phi);
        gamma.push_back(math::gamma(_xls[k]));
        R.push_back(E[k].topLeftCorner(3, 3));
        p.push_back(E[k].topRightCorner(3, 1));
        J.push_back(R[k] * gamma[k]);
        xw.push_back(R[k] * _xls[k] + p[k]);
        vw.push_back(J[k] * phi[k]);
    }

    // compute force vector
    // compute normalized force vector fn, l, and l_dot
    // l = |dx| = |x(k+1) - x(k)|
    // l_dot = dx' * dx_dot / l
    VectorX fn = VectorX::Zero(6 * n);
    dtype l = 0., l_dot = 0.;
    for (int k = 0;k < n - 1;k++) {
        Vector3 dx = xw[k + 1] - xw[k];
        Vector3 dx_dot = vw[k + 1] - vw[k];
        dtype lk = dx.norm();
        l += lk;
        l_dot += dx.dot(dx_dot) / lk;
        VectorX fk(12); fk << J[k].transpose() * dx, -J[k + 1].transpose() * dx;
        fn.segment(6 * k, 12) += fk / lk;
    }

    // compute spring force
    dtype fs;
    computeSpringForce(l, l_dot, fs);
    
    // compute f = fs * fn
    VectorX f = fs * fn;

    // fill fm
    for (int i = 0;i < n;i++)
        fm.segment(_bodies[i]->_index[0], 6) += f.segment(i * 6, 6);
}

void ForceSpringMultiPointGeneric::computeForceWithDerivative(VectorX& fm, VectorX& fr, MatrixX& Km, MatrixX& Dm, MatrixX& Kr, MatrixX& Dr, bool verbose) {
    int n = _bodies.size();
    std::vector<Matrix4> E;
    std::vector<MatrixX> gamma, J;
    std::vector<Vector6> phi;
    std::vector<Matrix3> R;
    std::vector<Vector3> p;
    std::vector<Vector3> xw, vw;
    for (int k = 0;k < n;k++) {
        E.push_back(_bodies[k]->_E_0i);
        phi.push_back(_bodies[k]->_phi);
        gamma.push_back(math::gamma(_xls[k]));
        R.push_back(E[k].topLeftCorner(3, 3));
        p.push_back(E[k].topRightCorner(3, 1));
        J.push_back(R[k] * gamma[k]);
        xw.push_back(R[k] * _xls[k] + p[k]);
        vw.push_back(J[k] * phi[k]);
    }

    // compute force vector
    // compute normalized force vector fn, l, and l_dot
    // l = |dx| = |x(k+1) - x(k)|
    // l_dot = dx' * dx_dot / l
    VectorX fn = VectorX::Zero(6 * n);
    dtype l = 0., l_dot = 0.;
    for (int k = 0;k < n - 1;k++) {
        Vector3 dx = xw[k + 1] - xw[k];
        Vector3 dx_dot = vw[k + 1] - vw[k];
        dtype lk = dx.norm();
        l += lk;
        l_dot += dx.dot(dx_dot) / lk;
        VectorX fk(12); fk << J[k].transpose() * dx, -J[k + 1].transpose() * dx;
        fn.segment(6 * k, 12) += fk / lk;
    }

    // compute spring force
    dtype fs, dfs_dl, dfs_dldot;
    computeSpringForceWithDerivative(l, l_dot, fs, dfs_dl, dfs_dldot, verbose);
    
    // compute f = fs * fn
    VectorX f = fs * fn;

    // fill fm
    for (int i = 0;i < n;i++)
        fm.segment(_bodies[i]->_index[0], 6) += f.segment(i * 6, 6);

    // compute K and D
    // K = dfs_dq * fn + fs * Kn
    Matrix3 I = Matrix3::Identity();

    // Kn = sum(dfk_dq) = K1 + K2
    MatrixX Kn = MatrixX::Zero(6 * n, 6 * n); // normalized vector part
    MatrixX dfs_dq = MatrixX::Zero(1, 6 * n);
    MatrixX dfs_dqdot = MatrixX::Zero(1, 6 * n);
    Matrix3 ex_skew = math::skew(Vector3::UnitX());
    Matrix3 ey_skew = math::skew(Vector3::UnitY());
    Matrix3 ez_skew = math::skew(Vector3::UnitZ());

    for (int k = 0;k < n - 1;k++) {
        Vector3 dx = xw[k + 1] - xw[k];
        Vector3 dx_dot = vw[k + 1] - vw[k];
        dtype lk = dx.norm();
        Vector3 dxnor = dx / lk;
        MatrixX dxnorT = dxnor.transpose();
        Vector3 Gphi0 = gamma[k] * phi[k];
        Vector3 Gphi1 = gamma[k + 1] * phi[k + 1];

        // K scalar term
        MatrixX Jx(3, 12); Jx << -J[k], J[k + 1];
        MatrixX dl_dq = dxnorT * Jx; // (1, 12)
        MatrixX dldot_dq = ((I - dxnor * dxnorT) / lk * dx_dot).transpose() * Jx; // (1, 12)
        dldot_dq(0, 0) += dxnor.dot(-R[k] * (ex_skew * Gphi0));
        dldot_dq(0, 6) += dxnor.dot(-R[k + 1] * (ex_skew * Gphi1));
        dldot_dq(0, 1) += dxnor.dot(-R[k] * (ey_skew * Gphi0));
        dldot_dq(0, 7) += dxnor.dot(-R[k + 1] * (ey_skew * Gphi1));
        dldot_dq(0, 2) += dxnor.dot(-R[k] * (ex_skew * Gphi0));
        dldot_dq(0, 8) += dxnor.dot(-R[k + 1] * (ex_skew * Gphi1));

        // contribute to dfs_dq
        dfs_dq.block(0, k * 6, 1, 12).noalias() += dfs_dl * dl_dq + dfs_dldot * dldot_dq;

        // K normalized vector term: K1
        VectorX fx(12); fx << J[k].transpose() * dx, -J[k + 1].transpose() * dx;
        Vector3 dk = -dx / (lk * lk * lk);
        MatrixX dlinv_dq(1, 12); dlinv_dq << dk.transpose() * J[k], -dk.transpose() * J[k + 1];
        MatrixX K1 = fx * dlinv_dq;

        // K normalized vector term: K2
        MatrixX K2(12, 12);
        Matrix3 x0_skew = math::skew(_xls[k]);
        Matrix3 x1_skew = math::skew(_xls[k + 1]);
        Matrix3 R1R0 = R[k + 1].transpose() * R[k];
        Matrix3 R0R1 = R1R0.transpose();
        K2.block(3, 0, 3, 3) = math::skew(R[k].transpose() * (p[k] - xw[k + 1]));
        K2.block(0, 0, 3, 3) = x0_skew * K2.block(3, 0, 3, 3);
        K2.block(9, 0, 3, 3) = R1R0 * x0_skew;
        K2.block(6, 0, 3, 3) = x1_skew * K2.block(9, 0, 3, 3);
        K2.block(0, 3, 3, 3) = x0_skew;
        K2.block(3, 3, 3, 3) = I;
        K2.block(9, 3, 3, 3) = -R1R0;
        K2.block(6, 3, 3, 3) = x1_skew * K2.block(9, 3, 3, 3);
        K2.block(3, 6, 3, 3) = R0R1 * x1_skew;
        K2.block(0, 6, 3, 3) = x0_skew * K2.block(3, 6, 3, 3);
        K2.block(9, 6, 3, 3) = math::skew(R[k + 1].transpose() * (p[k + 1] - xw[k]));
        K2.block(6, 6, 3, 3) = x1_skew * K2.block(9, 6, 3, 3);
        K2.block(3, 9, 3, 3) = -R0R1;
        K2.block(0, 9, 3, 3) = x0_skew * K2.block(3, 9, 3, 3);
        K2.block(6, 9, 3, 3) = x1_skew;
        K2.block(9, 9, 3, 3) = I;
        K2 /= lk;
        
        // copy to Kn
        Kn.block(k * 6, k * 6, 12, 12) += K1 + K2;

        // D: dfs_dqdot
        MatrixX d = dfs_dldot * dxnorT;
        dfs_dqdot.block(0, k * 6, 1, 6).noalias() += -d * J[k];
        dfs_dqdot.block(0, (k + 1) * 6, 1, 6).noalias() += d * J[k + 1];
    }

    MatrixX Ks = fn * dfs_dq;
    MatrixX K = Ks - fs * Kn; // check the sign
    MatrixX D = fn * dfs_dqdot;

    // copy to global
    for (int i = 0;i < n;i++)
        for (int j = 0;j < n;j++) {
            Km.block(_bodies[i]->_index[0], _bodies[j]->_index[0], 6, 6).noalias() += K.block(i * 6, j * 6, 6, 6);
            Dm.block(_bodies[i]->_index[0], _bodies[j]->_index[0], 6, 6).noalias() += D.block(i * 6, j * 6, 6, 6);
        }
}

}