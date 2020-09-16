#ifndef ALGORITHM_CG_H
#define ALGORITHM_CG_H

#include "AbsAlgorithm.h"
#include "cmd_line.h"
#include "ParamPNCG.h"

/** Solving Linear System:
 *  H*delta_x = -g,
 *  where H = ( Ai_J'*D(wi)*Ai_J + Ae_J'*Ae_J ) +sub_alpha *I_J
 *  	  g = ( Ai_J'*[wi]_{+} + Ae_J'*we ) + sub_alpha * (x_J -sub_d_J)
 *  	  w = A*x -sub_b;
 *  	  w += A*(delta_x)
 *        H * x = -g;
 */

template<typename L, typename D>
class CG:public AbsAlgorithm<L, D>{
public:
     L J,cg_iter;
protected:
    L m,n,mplus,nplus;
    std::vector<D> r,p;
    std::vector<D> Hp;
    D r_square,r_square_before;
    D pHp,alpha,beta, sub_alpha;
    D cg_tol;
    std::vector<L> JATD_row_ptr;
    std::vector<L> JATD_col_idx;
    std::vector<D> JATD_value;

public:
    CG(ProbData<L, D>* const data_inst, const std::string& data_file_path="");
    virtual ~CG() {}

    void step();
    void Hp_compute();
    void LS_active_solver(const D& sub_alpha, const std::vector<L> &active_index, const std::vector<D> &w,
                          const std::vector<D> &grad_J, bool result_iflog = false);

    void information_print();
};

template<typename L, typename D>
CG<L, D>::CG(ProbData<L, D> *const data_inst,const std::string& data_file_path):
        AbsAlgorithm<L, D>(data_inst,data_file_path) {
    m = data_inst->m;
    n = data_inst->n;
    nplus = data_inst->nplus;
    mplus = data_inst->mplus;
}

template<typename L, typename D>
void CG<L, D>::step() {

    beta = r_square/r_square_before;
    r_square_before = r_square;
    for (L row = 0; row < J; ++row)
        p[row] = r[row] + beta*p[row];
    Hp_compute();
    pHp = 0;
    for (L row = 0; row < J; ++row)
        pHp += p[row] * Hp[row];
    alpha = r_square_before/pHp;
    for (L row = 0; row < J; ++row) {
        AbsAlgorithm<L, D>::x[row] += alpha * p[row];
        r[row] -= alpha * Hp[row];
    }
    r_square = l2norm_square(J,r);

    ++cg_iter;
    ++ALGparam.total_cg_iterations;

    if(AbsAlgorithm<L, D>::result_iflog == true) {
        std::cout << "#iter=" << cg_iter << ", pHp = " << pHp << ", alpha=" << alpha;
        std::cout << ", beta =" << beta << ", r_square=" << r_square << std::endl;
    }
}

template<typename L, typename D>
void CG<L, D>::Hp_compute() {
    std::vector<D> vec(m,0);
    D tmp;
    for (L row = 0; row < J; ++row) {
        tmp = p[row];
        for (L col = JATD_row_ptr[row]; col < JATD_row_ptr[row+1]; ++col) {
            vec[JATD_col_idx[col]] +=  JATD_value[col] * tmp;
        }
    }
    D sum;
    for (L row = 0; row < J; ++row) {
        sum = 0;
        for (L col = JATD_row_ptr[row]; col < JATD_row_ptr[row+1]; ++col) {
            sum += JATD_value[col] * vec[JATD_col_idx[col]];
        }
        Hp[row] = sub_alpha * p[row] + sum;
    }
}

template<typename L, typename D>
void CG<L, D>::LS_active_solver(const D& sub_alpha, const std::vector<L> &active_index, const std::vector<D> &w, const std::vector<D> &grad_J, bool result_iflog) {
    AbsAlgorithm<L, D>::result_iflog = result_iflog;
    this->sub_alpha = sub_alpha;
    J = active_index.size();
    JATD_row_ptr = std::vector<L>(J+1);
    JATD_col_idx.clear();
    JATD_value.clear();
    //std::vector<L> test_A_col_idx(m,0);
    L row_index, col_index, length;
    for (L row = 0; row < J; ++row) {
        row_index = active_index[row];
        length = 0;
        for (L col = AbsAlgorithm<L,D>::data->AT_row_ptr[row_index]; col < AbsAlgorithm<L,D>::data->AT_row_ptr[row_index+1]; ++col) {
            col_index = AbsAlgorithm<L,D>::data->AT_col_idx[col];
            ALGparam.test_A_col_idx[col_index] = 1;
            if(col_index >= mplus || w[col_index] > 0) {
                JATD_col_idx.push_back(col_index);
                JATD_value.push_back(AbsAlgorithm<L,D>::data->AT_value[col]);
                ++length;
            }
        }
        JATD_row_ptr[row+1] = JATD_row_ptr[row] + length;
    }
    AbsAlgorithm<L, D>::A_effect_nnz = JATD_row_ptr.back();
    AbsAlgorithm<L,D>::real_A_row_idx.clear();
    for (L row = 0; row < m; ++row) {
        if(ALGparam.test_A_col_idx[row] ==1)
            AbsAlgorithm<L,D>::real_A_row_idx.push_back(row);
    }
   // Print(J,m,JATD_row_ptr,JATD_col_idx,JATD_value);
    r = std::vector<D>(J);
    p = std::vector<D>(J);
    AbsAlgorithm<L, D>::x = std::vector<D>(J);
    Hp = std::vector<D>(J);
    D tmp;
    for (L row = 0; row < J; ++row) {
        tmp = -grad_J[row];
        r[row] = tmp;
        p[row] = tmp;
    }
    r_square_before = l2norm_square(J,r);
    Hp_compute();
    pHp = 0;
    for (L row = 0; row < J; ++row)
        pHp += p[row] * Hp[row];
    alpha = r_square_before/pHp;
    for (L row = 0; row < J; ++row) {
        AbsAlgorithm<L, D>::x[row] += alpha * p[row];
        r[row] -= alpha * Hp[row];
    }
    r_square = l2norm_square(J,r);

    //cg_tol = std::min<D>(0.1 * r_square_before, r_square_before*r_square_before);
    cg_tol = std::min<D>(pow(PNCGparam.eta,2) * r_square_before, pow(r_square_before,1+PNCGparam.tau));

    cg_iter = 0;
    if(AbsAlgorithm<L, D>::result_iflog == true)
        information_print();
    while((r_square > cg_tol) && (cg_iter < J))
        step();
    //Print("cg_iter", cg_iter);
}

template<typename L, typename D>
void CG<L, D>::information_print() {
    printf("\033[31mMethod Information:\033[0m \n");
    std::cout << "method type=Conjugate Gradient Method for Linear System(Active Index)" << std::endl;
    std::cout << "sub_alpha=" << sub_alpha << ", active_index_size=" << J  << std::endl;
    std::cout << "real_A_row_idx_size = " << AbsAlgorithm<L, D>::real_A_row_idx.size();
    std::cout << ", JATD_norm" << sqrt(l2norm_square<L, D>(JATD_value)) << std::endl;
    std::cout << "cg_tol=" << cg_tol << ", r_square_initial=" << r_square_before << std::endl;
}


#endif //ALGORITHM_RESCG_H
