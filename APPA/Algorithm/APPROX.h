#ifndef ALGORITHM_APPROX_H
#define ALGORITHM_APPROX_H

#include "AbsAlgorithm.h"
#include "cmd_line.h"

template<typename L, typename D>
class APPROX:public AbsAlgorithm<L, D>{
protected:
    L m,n,nb,tau,s;
    L num_of_block_per_bacth;
    L per_check = 1;
    L right_per_check  = 1;
    L J,jj;
    std::vector<D> y;
    std::vector<D> z;
    std::vector<D> w_y;
    std::vector<D> w_z;
    std::vector<D> cord_grad_sum;
    std::vector<D> delta_y;
    std::vector<D> delta_z;
    std::vector<L> cord;
    std::vector<D> zeros_m;
    std::vector<D> zeros_n;
    std::vector<D> x_before;
    std::vector<L> active_set;
    D right;
    D theta, theta_square, theta_square_before;
    D x_diff_norm_square;
    D Lip_sum,C, mu_guess, res_bound;
    L K_iter;
    L variant_K, t, K0;
    std::vector<L> full_block;
    D precision_output_bound = 1e-1;
public:
    APPROX(ProbData<L, D>* const data_inst, const unsigned int& loss_type, const unsigned int& reg_type,
           const std::string& data_file_path="");
    virtual ~APPROX() {
    }

    void sub_initial();
    void sub_step();
    void step();
    void restart_step();

    void initialization();
    void PPA_solver(const D& sigma,const std::vector<D> &x_initial, const std::vector<D> &lambda_over_sigma_initial, const L &blocksize,
                    const D &varepsilon, const D &delta, bool result_iflog = false);

    void x_w_update();
    void x_diff_norm_square_update(D& x_diff_norm_square);

    void information_print();
    void PPA_pass_check();
    void result_print();

};

template<typename L, typename D>
APPROX<L, D>::APPROX(ProbData<L, D> *const data_inst, const unsigned int &loss_type,const unsigned int &reg_type,
                     const std::string& data_file_path):
        AbsAlgorithm<L, D>(data_inst, loss_type,reg_type,data_file_path) {
    tau = ALGparam.blocksize;
    AbsAlgorithm<L, D>::check_blocksize(tau);
    srand((unsigned)time(NULL));
    m = data_inst->m;
    n = data_inst->n;
    nb = data_inst->nb;
    zeros_m = std::vector<D>(m,0);
    zeros_n = std::vector<D>(n,0);
    y = std::vector<D>(n);
    w_y = std::vector<D>(m);
    z = std::vector<D>(n);
    w_z = std::vector<D>(m);
    AbsAlgorithm<L,D>::grad = std::vector<D>(n);
    AbsAlgorithm<L,D>::w = std::vector<D>(m);
    AbsAlgorithm<L,D>::x = std::vector<D>(n);
    AbsAlgorithm<L,D>::Lip = std::vector<D>(n);
    x_before = std::vector<D>(n);
    full_block = std::vector<L>(n);
    for (L row = 0; row < n; ++row)
        full_block[row] = row;
    num_of_block_per_bacth = n/tau;
    cord_grad_sum = std::vector<D>(tau);
    delta_y = std::vector<D>(tau);
    delta_z = std::vector<D>(tau);
    cord = std::vector<L>(tau);
    K0 = ceil(20*exp(1)*n/tau);

}

template<typename L, typename D>
void APPROX<L, D>::sub_initial(){
    theta = tau/(D)n;
    theta_square = theta * theta;
    y = zeros_n;
    w_y = zeros_m;
    z = AbsAlgorithm<L,D>::x;
    w_z = AbsAlgorithm<L,D>::w;
}

template<typename L, typename D>
void APPROX<L, D>::sub_step() {
    //y: u in paper; z: \tilde{z} in paper
    D tmp = n* theta/(D)tau;
    D tmp_u = -(1- tmp)/theta_square;

    //calculate coordinate grad of f at point theta^2*y +z
    AbsAlgorithm<L,D>::loss->cord_grad_sum_update(cord_grad_sum,AbsAlgorithm<L,D>::data,w_y,y,w_z,z,theta_square,tau,cord);
    //calculate z and delta_z: t in paper;
    AbsAlgorithm<L, D>::reg_cord_cal(z,delta_z,tmp,cord_grad_sum,tau,cord);
    //update y
    for (L row = 0; row < tau ; ++row){
        delta_y[row] = tmp_u * delta_z[row];
        y[cord[row]] += delta_y[row];
    }
    //update w_z and w_y
    AbsAlgorithm<L,D>::loss->w_update(w_z,AbsAlgorithm<L,D>::data,delta_z,tau,cord);
    AbsAlgorithm<L,D>::loss->w_update(w_y,AbsAlgorithm<L,D>::data,delta_y,tau,cord);
    //update theta
    theta = 0.5*(sqrt(theta_square*theta_square + 4*theta_square) - theta_square);
    theta_square = theta * theta;
    ++s;
}

template<typename L, typename D>
void APPROX<L, D>::step() {
    D start = clock();
    for (L i = 0; i < per_check; ++i) {
        random_shuffle(full_block.begin(), full_block.end());
        for (L i = 0; i < num_of_block_per_bacth ; ++i) {
            for (L row = 0; row < tau; ++row)
                cord[row] = full_block[row + i * tau];
            sub_step();
        }
        ++AbsAlgorithm<L, D>::iter;
    }
    x_w_update();
    AbsAlgorithm<L, D>::left_value_compute();
    AbsAlgorithm<L, D>::time += (clock() - start) / CLOCKS_PER_SEC;
    if(AbsAlgorithm<L, D>::result_iflog == true)
        AbsAlgorithm<L, D>::result_log();
}

template<typename L, typename D>
void APPROX<L, D>::restart_step(){
    D start = clock();
    for (L i = 0; i < per_check; ++i) {
        random_shuffle(full_block.begin(), full_block.end());
        for (L i = 0; i < num_of_block_per_bacth ; ++i) {
            for (L row = 0; row < tau; ++row)
                cord[row] = full_block[row + i * tau];
            if (s == 0) {
                sub_initial();
            }
            sub_step();
            if (s == AbsAlgorithm<L, D>::K) {
                x_w_update();
                s =0;
            }
        }
        ++AbsAlgorithm<L, D>::iter;
    }
    x_w_update();
    AbsAlgorithm<L, D>::time += (clock() - start) / CLOCKS_PER_SEC;
    if(AbsAlgorithm<L, D>::result_iflog == true)
        AbsAlgorithm<L, D>::result_log();
}

template<typename L, typename D>
void APPROX<L, D>::initialization() {
    AbsAlgorithm<L, D>::set_v_mu(tau);
    AbsAlgorithm<L,D>::loss->Lip_set(AbsAlgorithm<L,D>::Lip,AbsAlgorithm<L,D>::data);
    Lip_sum = 0.;
    for (L row = 0; row < n; ++row) {
        Lip_sum += AbsAlgorithm<L, D>::Lip[row];
    }
    AbsAlgorithm<L, D>::Krcd_compute(AbsAlgorithm<L, D>::K,tau);
    AbsAlgorithm<L, D>::iter = 0;
    AbsAlgorithm<L, D>::time = 0.;
    variant_K  = K0;
    //Restart iterations.
    K_iter = 0;
    //When s = K, go to next ResAPPROX iter.
    s = 0;
    //When t = variant_K, go to next variant ResAPPROX iter.
    t = 0;
    typename std::vector<D>::iterator maxPosition = std::max_element(AbsAlgorithm<L,D>::v.begin(),
                                                                     AbsAlgorithm<L,D>::v.end());
    random_shuffle(full_block.begin(), full_block.end());
    active_set = full_block;
    J = n;
    jj=0;
}

template<typename L, typename D>
void APPROX<L, D>::PPA_solver(const D& sigma, const std::vector<D> &x_initial, const std::vector<D> &lambda_over_sigma_initial, const L &blocksize,
                              const D &varepsilon, const D &delta, bool result_iflog) {
    AbsAlgorithm<L, D>::stop_type = 1;
    AbsAlgorithm<L, D>::method_type = "restart";
    AbsAlgorithm<L, D>::result_iflog = result_iflog;
    D sub_alpha;
    D start = clock();
    AbsAlgorithm<L, D>::PPA_sub_solver_initial(sub_alpha,sigma,x_initial,lambda_over_sigma_initial);
    D precision = varepsilon * sub_alpha;
    initialization();
    AbsAlgorithm<L,D>::time += (clock() - start)/CLOCKS_PER_SEC;
    if(AbsAlgorithm<L, D>::result_iflog == true) {
        information_print();
        printf("\033[31mPPA Information:\033[0m \n");
        std::cout << "sigma=" << sigma << ", varepsilon=" << varepsilon << ", delta=" << delta << std::endl;
        printf("\033[31mmin(varepsilon/(sigma*sigma),right) is the stopping criterion for PPA\033[0m \n");
    }
    D right_update_time, total_right_update_time = 0;
    right = 0;
    PPA_pass_check();
    while(AbsAlgorithm<L, D>::left > right) {
        restart_step();

        AbsAlgorithm<L,D>::loss->grad_compute(AbsAlgorithm<L,D>::grad,AbsAlgorithm<L,D>::data, AbsAlgorithm<L,D>::w,
                                              AbsAlgorithm<L,D>::x);
        D tmp_x;
        for (L row = 0; row < nb; ++row){
            tmp_x = AbsAlgorithm<L,D>::x[row]- 1./Lip_sum * AbsAlgorithm<L,D>::grad[row];
            if(tmp_x <= 0)
                AbsAlgorithm<L,D>::x[row] = 0;
            else
                AbsAlgorithm<L,D>::x[row] = tmp_x;
        }
        for (L row = nb; row < n; ++row){
            AbsAlgorithm<L,D>::x[row] -= 1./Lip_sum * AbsAlgorithm<L,D>::grad[row];
        }

        AbsAlgorithm<L,D>::loss->grad_compute(AbsAlgorithm<L,D>::grad,AbsAlgorithm<L,D>::data, AbsAlgorithm<L,D>::w,
                                              AbsAlgorithm<L,D>::x);
        AbsAlgorithm<L, D>::PPA_left_value_compute();
        PPA_pass_check();

        while (AbsAlgorithm<L, D>::left > precision) {
            restart_step();
            AbsAlgorithm<L,D>::loss->grad_compute(AbsAlgorithm<L,D>::grad,AbsAlgorithm<L,D>::data, AbsAlgorithm<L,D>::w,
                                                  AbsAlgorithm<L,D>::x);
            for (L row = 0; row < nb; ++row){
                tmp_x = AbsAlgorithm<L,D>::x[row]- 1./Lip_sum * AbsAlgorithm<L,D>::grad[row];
                if(tmp_x <= 0)
                    AbsAlgorithm<L,D>::x[row] = 0;
                else
                    AbsAlgorithm<L,D>::x[row] = tmp_x;
            }
            for (L row = nb; row < n; ++row){
                  AbsAlgorithm<L,D>::x[row] -= 1./Lip_sum * AbsAlgorithm<L,D>::grad[row];
            }

            AbsAlgorithm<L,D>::loss->grad_compute(AbsAlgorithm<L,D>::grad,AbsAlgorithm<L,D>::data, AbsAlgorithm<L,D>::w,
                                              AbsAlgorithm<L,D>::x);
            AbsAlgorithm<L, D>::PPA_left_value_compute();
            PPA_pass_check();
        }
        start = clock();
        AbsAlgorithm<L,D>::PPA_right_update(right,sigma,delta,x_initial,lambda_over_sigma_initial);
        right_update_time = (clock() - start) / CLOCKS_PER_SEC;
        if (s == 0) {
            sub_initial();
        }

        D tmp = n* theta/(D)tau;
        D tmp_u = -(1- tmp)/theta_square;

        //calculate coordinate grad of f at point theta^2*y +z
        AbsAlgorithm<L,D>::loss->cord_grad_sum_update(cord_grad_sum,AbsAlgorithm<L,D>::data,w_y,y,w_z,z,theta_square,tau,cord);

        /*
        //active strategy
        //----------------------------------------------------------------------
        L tmp_idx = cord[0];
        if(cord_grad_sum[0] > 0 && -1e-3 <= z[tmp_idx] <= 0) {
            std::swap(active_set[jj],active_set[J-1]);
            --J;
            --jj;
            continue;
        }
        ++jj;
        //----------------------------------------------------------------------
        */

        //calculate z and delta_z: t in paper;
        AbsAlgorithm<L, D>::reg_cord_cal(z,delta_z,tmp,cord_grad_sum,tau,cord);
        //update y
        for (L row = 0; row < tau ; ++row){
            delta_y[row] = tmp_u * delta_z[row];
            y[cord[row]] += delta_y[row];
        }
        //update w_z and w_y
        AbsAlgorithm<L,D>::loss->w_update(w_z,AbsAlgorithm<L,D>::data,delta_z,tau,cord);
        AbsAlgorithm<L,D>::loss->w_update(w_y,AbsAlgorithm<L,D>::data,delta_y,tau,cord);
        //update theta
        theta = 0.5*(sqrt(theta_square*theta_square + 4*theta_square) - theta_square);
        theta_square = theta * theta;
        ++s;

        if (s == AbsAlgorithm<L, D>::K) {
            x_w_update();
            s =0;
        }  //AbsAlgorithm<L, D>::time += (clock() - start) / CLOCKS_PER_SEC;
        AbsAlgorithm<L, D>::time += right_update_time;
        total_right_update_time += right_update_time;
    }

    if(result_iflog == true) {
        result_print();
        std::cout << "total_right_update_time=" << total_right_update_time << std::endl;
    }
}

template<typename L, typename D>
void APPROX<L, D>::x_w_update() {
    theta_square_before = theta_square/(1-theta_square);
    for (L row = 0; row < n; ++row) {
        AbsAlgorithm<L, D>::x[row] = theta_square_before * y[row] + z[row];
    }
    for (L row = 0; row < m; ++row) {
        AbsAlgorithm<L, D>::w[row] = theta_square_before * w_y[row] + w_z[row];
    }
}

template<typename L, typename D>
void APPROX<L, D>::x_diff_norm_square_update(D& x_diff_norm_square) {
    x_diff_norm_square = 0.;
    D tmp;
    for (L row = 0; row < n; ++row) {
        tmp = AbsAlgorithm<L, D>::x[row] - x_before[row];
        x_diff_norm_square += tmp*tmp;
    }
}

template<typename L, typename D>
void APPROX<L, D>::information_print() {
    printf("\033[31mMethod Information:\033[0m \n");
    std::cout << "method type=Fixed Period Restart APPROX" << std::endl;
    std::cout << "restart period=" << AbsAlgorithm<L, D>::K << ", block size=" << tau << std::endl;

    printf("\033[31mLoss Information:\033[0m \n");
    std::cout << "alpha=" << AbsAlgorithm<L, D>::loss->alpha << ", b_norm=";
    std::cout << sqrt(l2norm_square<L, D>(AbsAlgorithm<L, D>::loss->b)) << ", d_norm=";
    std::cout << sqrt(l2norm_square<L, D>(AbsAlgorithm<L, D>::loss->d)) <<std::endl;
}

template<typename L, typename D>
void APPROX<L, D>::PPA_pass_check() {
    std::cout << std::setprecision(10) << "#pass= " << AbsAlgorithm<L, D>::iter << ", obj=";
    std::cout << AbsAlgorithm<L, D>::fun_value << ", grad_norm=" << AbsAlgorithm<L, D>::left;
    std::cout << ", right=" << right << ", time=" << AbsAlgorithm<L, D>::time << std::endl;
    if( AbsAlgorithm<L, D>::left <= precision_output_bound){
        printf("\033[31m---------------------\033[0m");
        std::cout << "reach precision " << to_string(precision_output_bound) << ": total_pass= " << AbsAlgorithm<L, D>::iter;
        printf("\033[31m---------------------\033[0m \n");
        precision_output_bound /=10;
    }
}

template<typename L, typename D>
void APPROX<L, D>::result_print() {
    printf("\033[31mResult:\033[0m \n");
    std::cout << "initial_obj=" << AbsAlgorithm<L, D>::initial_fun_value;
    std::cout << ", initial_grad_norm=" << AbsAlgorithm<L, D>::initial_grad_norm << std::endl;
    std::cout << "final_obj=" << AbsAlgorithm<L, D>::fun_value;
    std::cout << ", final_grad_norm=" <<  AbsAlgorithm<L, D>::grad_norm << ", right=" << right << std::endl;
    std::cout << "restart_iter=" << AbsAlgorithm<L, D>::iter/AbsAlgorithm<L, D>::K;
    std::cout << ", total_pass=" << AbsAlgorithm<L, D>::iter << ", time=" << AbsAlgorithm<L, D>::time << std::endl;
}

#endif //ALGORITHM_RESAPPROX_H
