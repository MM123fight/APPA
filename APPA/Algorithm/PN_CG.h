#ifndef ALGORITHM_PN_CG_H
#define ALGORITHM_PN_CG_H
/* Silence_index = { i| x_i = 0, grad_i >0 }
 * Active_index = { i| x_i > 0 or grad_i <= 0 }
 */
#include "AbsAlgorithm.h"
//#include "CG.h"
#include "CG.h"

template<typename L, typename D>
class PN_CG:public AbsAlgorithm<L, D>{
private:
    L m,n,mi,nb;
    D sub_alpha;
    std::vector<D> grad_J;
    std::vector<L> active_index;
    std::vector<D> delta_x_J;
    std::vector<L> real_A_row_idx;
    L J;
    CG<L, D>* cg;
    L w_effect_size;
    unsigned int right_per_check = 1;
    unsigned int per_check = 1;
    unsigned int t;
    D precision_output_bound = 0.1;
public:
    PN_CG(ProbData<L, D>* const data_inst, const unsigned int& loss_type, const unsigned int& reg_type,
          const std::string& data_file_path="");
    ~PN_CG() {
        delete cg;
    }
    void step();
    void PPA_solver(const D& sigma, const std::vector<D> &x_initial, const std::vector<D> &lambda_over_sigma_initial, const L &blocksize,
                    const D &varepsilon, const D &delta, bool result_iflog = false);

    void PPA_pass_check();
    void information_print();
    void result_print();
};

template<typename L, typename D>
PN_CG<L, D>::PN_CG(ProbData<L, D> *const data_inst, const unsigned int &loss_type,const unsigned int &reg_type,
                   const std::string& data_file_path):
        AbsAlgorithm<L, D>(data_inst, loss_type,reg_type,data_file_path) {
    m = data_inst->m;
    n = data_inst->n;
    mi = data_inst->mi;
    nb = data_inst->nb;
    AbsAlgorithm<L,D>::grad = std::vector<D>(n);
    AbsAlgorithm<L,D>::w = std::vector<D>(m);
    AbsAlgorithm<L,D>::x = std::vector<D>(n);

    cg = new CG<L,D>(data_inst);
}

template<typename L, typename D>
void PN_CG<L, D>::step() {
    for (int i = 0; i < per_check; ++i) {
        
        active_index.clear();
        grad_J.clear();

        D grad_tmp;
        for (L row = 0; row < nb; ++row) {
            grad_tmp = AbsAlgorithm<L, D>::grad[row];
            if (grad_tmp <= 0 || AbsAlgorithm<L, D>::x[row] > 0) {
                active_index.push_back(row);
                grad_J.push_back(grad_tmp);
            }
        }
        for (L row = nb; row < n; ++row) {
            active_index.push_back(row);
            grad_J.push_back(AbsAlgorithm<L, D>::grad[row]);
        }
        J = active_index.size();

        //grad_J = AbsAlgorithm<L, D>::grad;
        ALGparam.test_A_col_idx = std::vector<L>(m,0);
        cg->LS_active_solver(sub_alpha, active_index, AbsAlgorithm<L, D>::w, grad_J, false);
        delta_x_J = cg->x;

        real_A_row_idx = cg->real_A_row_idx;
        L real_A_row_idx_size = real_A_row_idx.size();
        D tmp, temp1, temp2;
        L tmp_idx;

        D zeta = 1.;
        D cond,tmp1, tmp2, tmp3;
        D x_old, x_new, f_diff;

        //tmp1_vec: x_J - d_J
        std::vector<D> tmp1_vec(J);
        for (L row = 0; row < J; ++row) {
            tmp1_vec[row] = AbsAlgorithm<L, D>::x[active_index[row]]- AbsAlgorithm<L, D>::loss->d[active_index[row]];
        }

        //tmp3: -mu*<grad_J, delta_x_J>
        tmp3 = 0;
        for (L row = 0; row < J; ++row) {
            tmp3 += grad_J[row] * delta_x_J[row];
        }
        tmp3 *= -PNCGparam.mu;
        std::vector<D> delta_w_tmp(m);
        std::vector<D> delta_x_true(J);
        for (t = 0; t < PNCGparam.max_num_linesearch; ++t) {
            //delta_x_true = [x + zeta * delta_x]_{+} - x
            //tmp1: sub_alpha * <x -d, delta_x_true>
            //tmp2: 0.5 * sub_alpha * || delta_x_true||^2
            tmp1 = 0.;
            tmp2 = 0.;
            //Print("delta_x_J", l2norm_square<D>(delta_x_J));
            
            for (L row = 0; row < J; ++row) {
                tmp_idx = active_index[row];
                x_old = AbsAlgorithm<L, D>::x[tmp_idx];
                x_new = x_old + zeta * delta_x_J[row];
                if (tmp_idx < nb && x_new <= 0.){
                    tmp = -x_old;
                } else{
                    tmp = zeta * delta_x_J[row];
                }
                delta_x_true[row] = tmp;
                tmp1 += tmp1_vec[row] * tmp;
                tmp2 += tmp*tmp;
            }
            tmp1 *= sub_alpha;
            tmp2 *= (0.5 * sub_alpha);

            delta_w_tmp = std::vector<D>(m,0.);
            AbsAlgorithm<L, D>::loss->w_update(delta_w_tmp, AbsAlgorithm<L, D>::data, delta_x_true, J, active_index);

            //f_diff: 0.5*||[A(x+delta_x_true) -b]_{+}||^2 - 0.5*||[Ax-b]_{+}||^2
            f_diff = 0.;
            for (L row = 0; row < real_A_row_idx_size; ++row) {
                tmp_idx = real_A_row_idx[row];
                tmp = delta_w_tmp[tmp_idx];
                temp1 = AbsAlgorithm<L, D>::w[tmp_idx];
                temp2 = tmp + temp1; 
                if (tmp_idx < mi) {
                    if(temp1 > 0. && temp2 > 0. ){
                        f_diff += tmp*tmp + 2*tmp*temp1;
                    } else if(temp1 > 0. && temp2 <= 0.){
                        f_diff -= temp1*temp1;
                    } else if (temp1 <= 0. && temp2 > 0.){
                        f_diff += temp2*temp2;
                    }
                        
                }else{
                    f_diff += tmp*tmp + 2*tmp*temp1;
                }
            }
            f_diff *= 0.5;
            cond =  f_diff + tmp2 + tmp1;
            if (cond <= 0)
                break;
            else
                zeta *= PNCGparam.beta;
        }
        if( t == PNCGparam.max_num_linesearch ){
            std::cerr << "reach max_num_linesearch" << std::endl;
            PNCGparam.reach_linesearch = true;
            return;
        }

        for (L row = 0; row < J; ++row) {
            AbsAlgorithm<L, D>::x[active_index[row]] += delta_x_true[row];
        }
        

        for (L row = 0; row < m; ++row) {
            AbsAlgorithm<L, D>::w[row] += delta_w_tmp[row];
        }

        AbsAlgorithm<L, D>::loss->grad_compute(AbsAlgorithm<L, D>::grad, AbsAlgorithm<L, D>::data,
                                               AbsAlgorithm<L, D>::w,
                                               AbsAlgorithm<L, D>::x);
        ++AbsAlgorithm<L, D>::iter;
        ++ALGparam.total_pn_iterations;

        std::cout << "pn_iter=" << AbsAlgorithm<L, D>::iter  << ", cg_iter=" << cg->cg_iter << ", J=" << J;
        std::cout <<  ", zeta=" << zeta << ", A_effect_nnz_percent=" << std::setprecision(3);
        std::cout << cg->A_effect_nnz/(D)AbsAlgorithm<L, D>::data->A_row_ptr.back() * 100 << "%";      

    }
    AbsAlgorithm<L, D>::PPA_left_value_compute();
    std::cout << ", grad_norm=" << AbsAlgorithm<L, D>::grad_norm;

}

template<typename L, typename D>
void PN_CG<L, D>::PPA_solver(const D& sigma, const std::vector<D> &x_initial, const std::vector<D> &lambda_over_sigma_initial, const L &blocksize,
                              const D &varepsilon, const D &delta, bool result_iflog) {
    AbsAlgorithm<L, D>::stop_type = 1;
    AbsAlgorithm<L, D>::result_iflog = result_iflog;
    AbsAlgorithm<L, D>::iter = 0;
    AbsAlgorithm<L, D>::time = 0;
    D start = clock();
    AbsAlgorithm<L, D>::PPA_sub_solver_initial(sub_alpha,sigma,x_initial,lambda_over_sigma_initial);

    D precision = varepsilon * sub_alpha;
    AbsAlgorithm<L, D>::time += (clock() - start)/CLOCKS_PER_SEC;
    if(AbsAlgorithm<L, D>::result_iflog == true) {
        information_print();
        printf("\033[31mPPA Information:\033[0m \n");
        std::cout << "sigma=" << sigma << ", varepsilon=" << varepsilon << ", delta=" << delta << std::endl;
        printf("\033[31mmin(varepsilon/(sigma*sigma),right) is the stopping criterion for PPA\033[0m \n");
        std::cout << "#pass= " << AbsAlgorithm<L, D>::iter << ", obj=";
        std::cout << AbsAlgorithm<L, D>::fun_value << ", grad_norm=" << AbsAlgorithm<L, D>::left << std::endl;
    }
    AbsAlgorithm<L, D>::right = 0;
    start = clock();
    std::cout << "initial_grad_norm=" << AbsAlgorithm<L, D>::grad_norm << std::endl;
    while(AbsAlgorithm<L, D>::left > AbsAlgorithm<L, D>::right) {
        for (int i = 0; i < right_per_check; ++i) {
            step();
            if(PNCGparam.reach_linesearch == true)
                return;
            /*
            if(t == PNCGparam.max_num_linesearch){
                return;
            }*/
        }
        while ( AbsAlgorithm<L, D>::left > precision ) {
            std::cout << std::endl;
            step();
            if(PNCGparam.reach_linesearch == true)
                return;
            /*
            if(t == PNCGparam.max_num_linesearch){
                return;
            }
             */
            if(AbsAlgorithm<L, D>::iter > PNCGparam.pn_max_iter)
                break;
        }

        AbsAlgorithm<L,D>::PPA_right_update(AbsAlgorithm<L, D>::right,sigma,delta,x_initial,lambda_over_sigma_initial);
        std::cout <<  ", right=" << AbsAlgorithm<L, D>::right << std::endl;
    }
    AbsAlgorithm<L, D>::time += (clock() - start) / CLOCKS_PER_SEC;

    if(result_iflog == true) {
        result_print();
    }


}

template<typename L, typename D>
void PN_CG<L, D>::PPA_pass_check() {
    std::cout << "#pass= " << AbsAlgorithm<L, D>::iter << ", obj=";
    std::cout << AbsAlgorithm<L, D>::fun_value << ", grad_norm=" << AbsAlgorithm<L, D>::grad_norm << ", grad_norm_infty=";
    std::cout << AbsAlgorithm<L, D>::grad_norm_infty << ", right=" << AbsAlgorithm<L, D>::right << ", active_idx_size=" << J;
    std::cout << ", w_effect_size=" << w_effect_size << ", A_effect_nnz_percent=" ;
    std::cout << std::setprecision(3) << cg->A_effect_nnz/(D)AbsAlgorithm<L, D>::data->A_row_ptr.back() * 100 << "%";
    std::cout << std::setprecision(ALGparam.prec) << ", time=" << AbsAlgorithm<L, D>::time << std::endl;
}

template<typename L, typename D>
void PN_CG<L, D>::information_print() {
    printf("\033[31mMethod Information:\033[0m \n");
    std::cout << "method type=Projection Newton-Conjugate Gradient Method" << std::endl;
    printf("\033[31mLoss Information:\033[0m \n");
    std::cout << "alpha=" << AbsAlgorithm<L, D>::loss->alpha << ", b_norm=";
    std::cout << sqrt(l2norm_square<L, D>(AbsAlgorithm<L, D>::loss->b)) << ", d_norm=";
    std::cout << sqrt(l2norm_square<L, D>(AbsAlgorithm<L, D>::loss->d)) << ", A_total_nnz=";
    std::cout << AbsAlgorithm<L, D>::data->A_row_ptr.back() << std::endl;
}

template<typename L, typename D>
void PN_CG<L, D>::result_print() {
    printf("\033[31mResult:\033[0m \n");
    std::cout << "initial_obj=" << AbsAlgorithm<L, D>::initial_fun_value;
    std::cout << ", initial_grad_norm=" << AbsAlgorithm<L, D>::initial_grad_norm << std::endl;
    std::cout << "final_obj=" << AbsAlgorithm<L, D>::fun_value;
    std::cout << ", final_grad_norm=" <<  AbsAlgorithm<L, D>::grad_norm << ", right=" << AbsAlgorithm<L, D>::right << std::endl;
    std::cout << "total_pass=" << AbsAlgorithm<L, D>::iter << ", time=" << AbsAlgorithm<L, D>::time << std::endl;
}

#endif //ALGORITHM_PN_CG_H
