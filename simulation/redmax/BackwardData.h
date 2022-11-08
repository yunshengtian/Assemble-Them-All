#include "Common.h"

namespace redmax {

class BackwardInfo {
public:
    // computed values in forward
    // For BDF2 integrator, the item at index 0 is for time step alpha
    std::vector<MatrixX> _H_inv; // inv(dg(k)/dq(k)), (ndof_r x ndof_r)
    std::vector<Eigen::PartialPivLU<MatrixX> > _H_lu;
    std::vector<MatrixX> _M; // mass matrix, (ndof_r x ndof_r)
    std::vector<MatrixX> _D; // damping matrix, (ndof_r x ndof_r)
    std::vector<MatrixX> _dg_dp; // dg(k)/dp
    std::vector<MatrixX> _dg_du; // dg(k)/du(k).
    std::vector<MatrixX> _dvar_dq; // dvar(t)/dq(t), to be noticed: var only save for full step, not for alpha step
    std::vector<MatrixX> _dvar_dp; // dvar(t)/dp
    
    // terminal derivatives set by task before backward
    // p = (q0, qdot0, u)
    VectorX _df_dq0;         // ndof_r    
    VectorX _df_dqdot0;      // ndof_r
    VectorX _df_dp;          // ndof_p (design parameters)
    VectorX _df_du;          // ndof_u * T
    // q = q(1), ..., q(T)
    VectorX _df_dq;          // ndof_r * T
    // var = var(1), ..., var(T)
    VectorX _df_dvar;        // ndof_v * T

    bool _flag_q0, _flag_qdot0, _flag_p, _flag_u;   // flag to indicate the active optimization variables

    BackwardInfo() {
        clear();
    }

    void clear() {
        _H_inv.clear();
        _H_lu.clear();
        _M.clear();
        _D.clear();
        _dg_dp.clear();
        _dg_du.clear();
        _dvar_dq.clear();
        _dvar_dp.clear();
    }

    void set_flags(bool flag_q0, bool flag_qdot0, bool flag_p, bool flag_u) {
        _flag_q0 = flag_q0;
        _flag_qdot0 = flag_qdot0;
        _flag_p = flag_p;
        _flag_u = flag_u;
    }
};

class BackwardResults {
public:
    VectorX _df_dq0;
    VectorX _df_dqdot0;
    VectorX _df_dp;
    VectorX _df_du;
};

}