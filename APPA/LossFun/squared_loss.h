#ifndef LOSSFUN_SQUARED_LOSS_H
#define LOSSFUN_SQUARED_LOSS_H

#include "AbsLoss.h"

/*
 * 0.5||Ax-b||^2  = sum_{i=1}^m {0.5 ||a_ix -b_i||^2}
 * w = Ax -b 
 */

template<typename L, typename D>
class squared_loss:public AbsLoss<L, D>{

public:

    squared_loss(const std::vector<D>& b_value):
            AbsLoss<L, D>(b_value){}
    ~squared_loss(){}
    //grad(x) = AT(Ax-b)
    void cord_grad_update(std::vector<D>& cord_grad, const ProbData<L, D>* const data, const std::vector<D>& w,
                          const std::vector<D>& x, const L& blocksize, const std::vector<L>& cord);
    /* w_y = A*y, w_z = A*z-b;
     * grad(theta_y*y + z) = A(theta_y *w_y + w_z)
     */
    void cord_grad_sum_update(std::vector<D>& cord_grad_sum, const ProbData<L, D>* const data, const std::vector<D>& w_y,
                              const std::vector<D>& y, const std::vector<D>& w_z, const std::vector<D>& z,
                              const D& theta_y, const L& blocksize, const std::vector<L>& cord);
    void fun_value_compute(D& fun_value, const ProbData<L, D>* const data, const std::vector<D>& w, const std::vector<D>& x);
    void v_set(std::vector<D>& v, const ProbData<L, D>* const data, const L& blocksize);
    void Lip_set(std::vector<D>& Lip, const ProbData<L, D>* const data);
};

template<typename L, typename D>
void squared_loss<L, D>::cord_grad_update(std::vector<D>& cord_grad, const ProbData<L, D>* const data, const std::vector<D>& w,
                                               const std::vector<D>& x, const L& blocksize, const std::vector<L>& cord){

    L tmp_cord;
    for (L row = 0; row < blocksize; ++row) {
        tmp_cord = cord[row];
        cord_grad[row] = 0;
        for (L col = data->AT_row_ptr[tmp_cord]; col < data->AT_row_ptr[tmp_cord + 1]; ++col)
            cord_grad[row] += data->AT_value[col] * w[data->AT_col_idx[col]];
    }

}

template<typename L, typename D>
void squared_loss<L, D>::cord_grad_sum_update(std::vector<D>& cord_grad_sum, const ProbData<L, D>* const data,
                                                   const std::vector<D>& w_y, const std::vector<D>& y,
                                                   const std::vector<D>& w_z, const std::vector<D>& z,
                                                   const D& theta_y, const L& blocksize, const std::vector<L>& cord) {
    L tmp_cord,tmp_idx;
    for (L row = 0; row < blocksize; ++row) {
        tmp_cord = cord[row];
        cord_grad_sum[row] = 0;
        for (L col = data->AT_row_ptr[tmp_cord]; col < data->AT_row_ptr[tmp_cord + 1]; ++col) {
            tmp_idx = data->AT_col_idx[col];
            cord_grad_sum[row] += data->AT_value[col] * (theta_y * w_y[tmp_idx] + w_z[tmp_idx]);

        }
    }

}

template<typename L, typename D>
void squared_loss<L, D>::fun_value_compute(D& fun_value,const ProbData<L, D> *const data, const std::vector<D>& w, const std::vector<D>& x){
    fun_value = 0;
    for (L row = 0; row < data->m; ++row){
            fun_value += w[row]* w[row];
    }
    fun_value *= 0.5;
};

template<typename L, typename D>
void squared_loss<L, D>::v_set(std::vector<D> &v, const ProbData<L, D> *const data, const L &blocksize) {
    AbsLoss<L, D>::v_initial_set(v, data, blocksize);
}

template<typename L, typename D>
void squared_loss<L, D>::Lip_set(std::vector<D> &Lip, const ProbData<L, D> *const data) {
    AbsLoss<L, D>::Lip_initial_set(Lip,data);
}


#endif //LOSSFUN_ABSLOSS_H

