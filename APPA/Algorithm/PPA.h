#ifndef ALGORITHM_PPA_H
#define ALGORITHM_PPA_H
#include "AbsAlgorithmSub.h"
#include "AbsAlgorithm.h"
#include "ParamPPA.h"

//min ||x_{s,t+1} - x_{s,t}||
//initial point for next iteration: x_{s,t}      
template<typename L, typename D>
class PPA:public AbsAlgorithmSub<L, D>{
private:
    L m,mi,n,nb,tau;
    L t;

    //parameters for PPA
    D sigma;
    D varepsilon, varepsilon_stage1;
    D gap;
    L PPA_iter = 0;
    D PPA_time = 0.;
    D Lip, max_Ai;
    D prob = 0.05;
    ParamPPA<D>* param;
    std::vector<D> zeros_n;

    bool precision_output_aver = false;
    D precision_output_aver_bound = 1e-1;
    bool change_method = false;
    bool sigma_first = false;
    bool sigma_second = false;
    std::vector<D> x_diff_over_sigma;
    std::vector<D> lambda_over_sigma;
    D c_norm, b_norm;
    D mpdg,mpdg_bound;
    D p_aver,d_aver,g_aver, g_aver_over, mpdg_aver;

    D mpdg_aver_min;
    D mpdg_min, gap_min;
    std::vector<D> x_min;
    std::vector<D> lambda_over_sigma_min;


public:
    PPA(ProbData<L, D>* const data_inst, const unsigned int& sub_loss_type,
        const unsigned int& sub_reg_type,const unsigned int& sub_method_type,
        const std::string& data_file_path="");

    virtual ~PPA(){}

    void sub_step();
    void restart_step();

    void solver(const std::vector<D>& x_initial,const std::vector<D>& lambda_initial,const L& sub_method_blocksize,
                const L& max_iter, const D& precision, ParamPPA<D>* const param,
                const bool& result_iflog = false, const std::string& method_type = "");

    void left_res_initial();
    void left_res_compute();

    void PPA_intial_result_log();
    void PPA_result_log();
    void result_print();
};

template<typename L, typename D>
PPA<L, D>::PPA(ProbData<L, D>* const data_inst, const unsigned int& sub_loss_type,
               const unsigned int& sub_reg_type,const unsigned int& sub_method_type,
               const std::string& data_file_path):
        AbsAlgorithmSub<L,D>(data_inst,sub_loss_type, sub_reg_type, sub_method_type,data_file_path){
    m = data_inst->m;
    mi = data_inst->mi;
    n = data_inst->n;
    nb = data_inst->nb;
    zeros_n = std::vector<D>(n);
    x_diff_over_sigma = std::vector<D>(n);
    lambda_over_sigma = std::vector<D>(m);

    AbsAlgorithmSub<L, D>::x = std::vector<D>(n);
    AbsAlgorithmSub<L, D>::w = std::vector<D>(m);
    AbsAlgorithmSub<L, D>::lambda = std::vector<D>(m);
    AbsAlgorithmSub<L, D>::w_dual = std::vector<D>(n);
    b_norm = data_inst->b_norm + 1;
    c_norm = data_inst->c_norm + 1;
    Lip = pow(data_inst->A_norm_F,2);
    this->max_Ai = data_inst->max_Ai;

    x_min = std::vector<D>(n);
    lambda_over_sigma_min = std::vector<D>(m);

}

//PPA_step
template<typename L, typename D>
void PPA<L, D>::sub_step() {

    PNCGparam.pn_max_iter = ceil(log(2*pow((1+param->delta)*Lip*sigma/param->delta,2)*(PPA_iter+1)/prob))*sigma*max_Ai;
    std::cout << "pn_max_iter=" << PNCGparam.pn_max_iter << std::endl;
    if(PNCGparam.reach_linesearch==true){
        std::cout << "transfer to PN" << std::endl;
        PNCGparam.reach_linesearch=false;
        AbsAlgorithmSub<L, D>::method_transfer();
    }
    AbsAlgorithmSub<L, D>::sub_method->PPA_solver(sigma,AbsAlgorithmSub<L, D>::x,lambda_over_sigma,
                                                  ALGparam.blocksize,varepsilon,param->delta,false);
    if(PNCGparam.reach_linesearch==true){
        std::cout << "transfer back to approx" << std::endl;
        AbsAlgorithmSub<L, D>::method_transfer_back();
        AbsAlgorithmSub<L, D>::sub_method->PPA_solver(sigma,AbsAlgorithmSub<L, D>::x,lambda_over_sigma,
                                                  ALGparam.blocksize,varepsilon,param->delta,false);
    }
    

    for (L row = 0; row < n; ++row)
        x_diff_over_sigma[row] = (AbsAlgorithmSub<L, D>::sub_method->x[row] - AbsAlgorithmSub<L, D>::x[row])/sigma;
    AbsAlgorithmSub<L, D>::x = AbsAlgorithmSub<L, D>::sub_method->x;
    /*
    D tmp;
    for (L row = 0; row < m; ++row){
        tmp = AbsAlgorithmSub<L, D>::sub_method->w[row];
        AbsAlgorithmSub<L, D>::w[row] = tmp - lambda_over_sigma[row];
        if((row < mi)&&(tmp < 0))
            tmp = 0;
        lambda_over_sigma[row] = tmp;
    }
     */
    D tmp;
    for (L row = 0; row < mi; ++row){
        tmp = AbsAlgorithmSub<L, D>::sub_method->w[row];
        AbsAlgorithmSub<L, D>::w[row] = tmp - lambda_over_sigma[row];
        if(tmp < 0) {
            lambda_over_sigma[row] = 0;
        }else{
            lambda_over_sigma[row] = tmp;
        }
    }

    for (L row = mi; row < m; ++row){
        tmp = AbsAlgorithmSub<L, D>::sub_method->w[row];
        AbsAlgorithmSub<L, D>::w[row] = tmp - lambda_over_sigma[row];
        lambda_over_sigma[row] = tmp;
    }

    ++t;
    ++PPA_iter;
}

template<typename L, typename D>
void PPA<L, D>::restart_step() {
    D start, minus_start, minus_time;
    start = clock();
    inv_update<L,D>(varepsilon, varepsilon_stage1,param->inv_idx_varepsilon,t);
    sub_step();
    mpdg = AbsAlgorithmSub<L, D>::sub_method->right;
    mpdg_bound = param->C * mpdg;
    mpdg_min = mpdg;
    left_res_compute();
    PPA_time += (clock() - start)/(D)CLOCKS_PER_SEC;
    PPA_result_log();
    while( (mpdg <=  mpdg_bound) && (mpdg_aver > AbsAlgorithmSub<L, D>::precision) ){
        start = clock();
        //do PPA(K): with param->delta and sigma to be constant
        /********************************************************************************************/
        inv_update<L,D>(varepsilon, varepsilon_stage1,param->inv_idx_varepsilon,t);
        sub_step();
        mpdg = AbsAlgorithmSub<L, D>::sub_method->right;
        left_res_compute();
        /********************************************************************************************/
        mpdg_bound = std::min<D>(param->C * mpdg, param->rho * mpdg_bound);
        if(ALGparam.type==1){
        if (mpdg_aver < mpdg_aver_min){
            mpdg_aver_min = mpdg_aver;
            x_min = AbsAlgorithmSub<L, D>::x;
            lambda_over_sigma_min = lambda_over_sigma;
        }  
        }

        PPA_time += (clock() - start)/(D)CLOCKS_PER_SEC;
        PPA_result_log();
        if(change_method == false && p_aver <= ALGparam.tol_trans){
            printf("\033[31m------------------------------------------------------\033[0m");
            if(ALGparam.sub_method == 1) {
                std::cout << "Change method to PN_CG";
            }
            else if(ALGparam.sub_method == 2){
                std::cout << "Change method to PN_CG_T";
            }
            printf("\033[31m------------------------------------------------------\033[0m \n");
            AbsAlgorithmSub<L, D>::method_transfer();
            change_method = true;
        }
    }
    printf("\033[31m-----------------------------------------------------------------\033[0m \n");
    printf("\033[31m-----------------------------------------------------------------\033[0m \n");
    start = clock();
    ++AbsAlgorithmSub<L, D>::iter;
    if(mpdg_aver <= AbsAlgorithmSub<L, D>::precision){
        result_print();
        exit(0);
    }
    //if(sigma > 2){
       // param->param_set(0.9,1);
    //}

    /*
    if(mpdg_aver <= 0.1 && sigma_first == false){
        param->param_set(0.7,2);
        sigma_first = true;
        printf("\033[31mParameters of PPA change as :\033[0m \n");
        std::cout << "rho="<< param->rho << ", delta="<< param->delta << ", alpha=" << param->alpha;
        std::cout << ", C=" << param->C << ", exp_idx_kappa=" << param->exp_idx_kappa << std::endl;
    }

    if(mpdg_aver <= 1e-3 && sigma_first == false){
        param->param_set(0.5,5);
        sigma_first = true;
        printf("\033[31mParameters of PPA change as :\033[0m \n");
        std::cout << "rho="<< param->rho << ", delta="<< param->delta << ", alpha=" << param->alpha;
        std::cout << ", C=" << param->C << ", exp_idx_kappa=" << param->exp_idx_kappa << std::endl;
    }

    
    if(mpdg_aver <= 1e-6 && sigma_second == false){
        param->param_set(0.9,5);
        sigma_second = true;
        printf("\033[31mParameters of PPA change as :\033[0m \n");
        std::cout << "rho="<< param->rho << ", delta="<< param->delta << ", alpha=" << param->alpha;
        std::cout << ", C=" << param->C << ", exp_idx_kappa=" << param->exp_idx_kappa << std::endl;
    }
    */


    if (ALGparam.type == 1){
        mpdg_aver = mpdg_aver_min;
        AbsAlgorithmSub<L, D>::x = x_min;
    }

    if (ALGparam.type == 0){
        for (L row = 0; row < m; ++row) {
          lambda_over_sigma[row] /= param->exp_idx_kappa;
        }
    } else{
    for (L row = 0; row < m; ++row) {
        lambda_over_sigma_min[row] /= param->exp_idx_kappa;
        lambda_over_sigma[row] = lambda_over_sigma_min[row];
    }
    }    
    sigma *= param->exp_idx_kappa;
    varepsilon_stage1 *= param->exp_idx_varepsilon;

    t = 0;
    PPA_time += (clock() - start)/(D)CLOCKS_PER_SEC;
}

template<typename L, typename D>
void PPA<L, D>::solver(const std::vector<D>& x_initial, const std::vector<D>& lambda_initial, const L& sub_method_blocksize,
                       const L& max_iter, const D& precision,ParamPPA<D>* const param,
                       const bool& result_iflog, const std::string& method_type){
    tau = sub_method_blocksize;

    printf("\033[31mMethod type:\033[0m \n");
    std::cout <<"APPA:";
    if(ALGparam.sub_method == 1)
        std::cout << "PN_CG";
    else
        std::cout << "PN_CG_T";
    if(ALGparam.solve_with_scale == true)
        std::cout << ", scaling";
    else
        std::cout << ", non-scaling";
    if(ALGparam.type == 1){
        std::cout << ", min aver";
    }else if(ALGparam.type == 0){
        std::cout << ", last point";
    }else{
        std::cout << ", there is no such type";
    }
    std::cout << std::endl;


    if(result_iflog == true) {
        //AbsAlgorithmSub<L, D>::full_result_path = AbsAlgorithmSub<L, D>::data->data_file_path + "/result/APPA_sub_bs" + to_string(tau);
        AbsAlgorithmSub<L, D>::full_result_path = AbsAlgorithmSub<L, D>::data->data_file_path + "/result/APPA_h" + to_string(ALGparam.rho) + "_e" + to_string(ALGparam.exp_idx_kappa) + "_p" + to_string(ALGparam.type);
        if(ALGparam.solve_from_dual == true){
            AbsAlgorithmSub<L, D>::full_result_path += "_d";
        }
        //if(ALGparam.varepsilon < 1e16){
            AbsAlgorithmSub<L, D>::full_result_path += "_v" + to_string(ALGparam.varepsilon);
        //}

	AbsAlgorithmSub<L, D>::logfile_open();
    }

    AbsAlgorithmSub<L, D>::x = x_initial;
    AbsAlgorithmSub<L, D>::max_iter = max_iter;
    AbsAlgorithmSub<L, D>::precision = precision;
    this->param = param;
    AbsAlgorithmSub<L, D>::result_iflog = result_iflog;
    // If we have different version PPA, we may need method_type to discuss.
    // So far, we only focus on the APPA and method_type is useless.
    AbsAlgorithmSub<L, D>::method_type = method_type;

    varepsilon_stage1 = param->varepsilon;
    varepsilon = param->varepsilon;
    sigma = param->alpha * param->kappa;

    // w = Ax-b w_dual = AT*lambda+c
    AbsAlgorithmSub<L, D>::loss->w_initial(AbsAlgorithmSub<L, D>::w,AbsAlgorithmSub<L, D>::data, x_initial,
                                           AbsAlgorithmSub<L, D>::data->b);
    if(l2norm_square<L,D>(lambda_initial)==0 && l2norm_square<L,D>(x_initial) == 0) {
        // Do one dual update instead of do the primal update first.
        AbsAlgorithmSub<L, D>::lambda = AbsAlgorithmSub<L, D>::w;
        for (L row = 0; row < m; ++row) {
            AbsAlgorithmSub<L, D>::lambda[row] *= sigma;
        }
    }
    else {
        AbsAlgorithmSub<L, D>::lambda = lambda_initial;
    }
    AbsAlgorithmSub<L, D>::loss->w_dual_initial(AbsAlgorithmSub<L, D>::w_dual,AbsAlgorithmSub<L, D>::data,
                                                AbsAlgorithmSub<L, D>::lambda,AbsAlgorithmSub<L, D>::data->c);
    // In our code, we use lambda_over_sigma = lambda/sigma instead of lambda for convenience.
    for (L row = 0; row < m; ++row) {
        lambda_over_sigma[row] = AbsAlgorithmSub<L, D>::lambda[row]/sigma;
    }

    left_res_initial();

    if(ALGparam.type == 1){
        mpdg_aver_min = mpdg_aver;
        x_min = AbsAlgorithmSub<L, D>::x;
        lambda_over_sigma_min = lambda_over_sigma;
    }

    t = 0;
    //initial_check();
    printf("\033[31miter = PPA_iter, #inner = pass for every PPA_iter, p_aver = relative primal feasibility, d_aver = relative dual feasibility, \033[0m \n");
    printf("\033[31mgap = primal dual gap :\033[0m \n");
    PPA_intial_result_log();
    AbsAlgorithmSub<L, D>::iter = 0;
    if(change_method == false && p_aver <= ALGparam.tol_trans){

        printf("\033[31m------------------------------------------------------\033[0m");
        if(ALGparam.sub_method == 1) {
            std::cout << "Change method to PN_CG";
        }
        else if(ALGparam.sub_method == 2){
            std::cout << "Change method to PN_CG_T";
        }
        printf("\033[31m------------------------------------------------------\033[0m \n");
        AbsAlgorithmSub<L, D>::method_transfer();
        change_method = true;
    }
    while (mpdg_aver > precision  && AbsAlgorithmSub<L, D>::iter < max_iter){
        restart_step();
    }

}

template<typename L, typename D>
void PPA<L, D>::left_res_initial(){
    AbsAlgorithmSub<L, D>::sub_method->grad = zeros_n;
    x_diff_over_sigma = zeros_n;
    left_res_compute();
    D tmp;
    d_aver = 0;
    for (L row = 0; row < nb; ++row) {
        tmp = -AbsAlgorithmSub<L, D>::w_dual[row];
        if(tmp > 0)
            d_aver += tmp * tmp;
    }
    for (L row = nb; row < n; ++row) {
        tmp = -AbsAlgorithmSub<L, D>::w_dual[row];
        d_aver += tmp * tmp;
    }


    d_aver = sqrt(d_aver);
    d_aver /= c_norm;

    mpdg_aver = std::max(g_aver,std::max(p_aver,d_aver));
}

template<typename L, typename D>
void PPA<L, D>::left_res_compute() {
    //gap
    D tmp;
    //<c, x>
    AbsAlgorithmSub<L, D>::fun_value = 0;
    g_aver_over = 0;
    for (L row = 0; row < n; ++row)
        g_aver_over += AbsAlgorithmSub<L, D>::data->c[row] * AbsAlgorithmSub<L, D>::x[row];
    AbsAlgorithmSub<L, D>::fun_value = g_aver_over;
    g_aver_over = abs(g_aver_over);
    // <b, lambda>
    gap = 0;
    for (L row = 0; row < m; ++row)
        gap += AbsAlgorithmSub<L, D>::data->b[row]*lambda_over_sigma[row];
    gap *= sigma;
    g_aver_over += abs(gap) + 1;
    gap += AbsAlgorithmSub<L, D>::fun_value;
    tmp = abs(gap);
    g_aver = tmp/g_aver_over;
    //primal_res
    p_aver = 0;
    for (L row = 0; row < mi; ++row){
        tmp = AbsAlgorithmSub<L, D>::w[row];
        if(tmp > 0)
            p_aver += tmp*tmp;
    }
    for (L row = mi; row < m; ++row){
        tmp = abs(AbsAlgorithmSub<L, D>::w[row]);
        p_aver += tmp*tmp;
    }
    //dual_res
    d_aver = 0;
    for (L row = 0; row < nb; ++row) {
        tmp = -(sigma * AbsAlgorithmSub<L, D>::sub_method->grad[row] - x_diff_over_sigma[row]);
        if(tmp > 0)
            d_aver += tmp * tmp;
    }
    for (L row = nb; row < n; ++row) {
        tmp = -(sigma * AbsAlgorithmSub<L, D>::sub_method->grad[row] - x_diff_over_sigma[row]);
        d_aver += tmp * tmp;
    }

    p_aver = sqrt(p_aver);
    p_aver /= b_norm;
    d_aver = sqrt(d_aver);
    d_aver /= c_norm;

    mpdg_aver = std::max(g_aver,std::max(p_aver,d_aver));

}

template<typename L, typename D>
void PPA<L, D>::PPA_intial_result_log() {
    AbsAlgorithmSub<L, D>::logFile << std::setprecision(10);
    AbsAlgorithmSub<L, D>::logFile << PPA_iter << "\t"  << g_aver << "\t" << mpdg_aver << "\t" << 0. << std::endl;
    
    std::cout << std::setprecision(ALGparam.prec);
    std::cout << "iter= " << PPA_iter << ", obj= " << AbsAlgorithmSub<L, D>::fun_value;
    std::cout << ", p_aver=" << p_aver << ", d_aver=" << d_aver << ", g_aver=" << g_aver;
    std::cout << ", gap=" << gap;
    std::cout << ", sigma=" << sigma  << ", param->delta=" << param->delta<< ", time = " << PPA_time << std::endl;
}

template<typename L, typename D>
void PPA<L, D>::PPA_result_log() {
    AbsAlgorithmSub<L, D>::logFile << std::setprecision(10);
    AbsAlgorithmSub<L, D>::logFile << PPA_iter << "\t"  << g_aver << "\t" << mpdg_aver << "\t" << PPA_time << std::endl;
    
    std::cout << std::setprecision(ALGparam.prec);
    std::cout << "iter= " << PPA_iter << ", #inner=" << AbsAlgorithmSub<L, D>::sub_method->iter;
    std::cout << ", obj= " << AbsAlgorithmSub<L, D>::fun_value;
    std::cout << ", p_aver=" << p_aver << ", d_aver=" << d_aver << ", g_aver=" << g_aver;
    std::cout << ", gap=" << gap;
    std::cout << ", sigma=" << sigma  << ", param->delta=" << param->delta<< ", time = " << PPA_time << std::endl;

    //std::cout << "x_norm" << sqrt(l2norm_square<L, D>(AbsAlgorithmSub<L, D>::x)) << std::endl;

    if(precision_output_aver == false && mpdg_aver <= precision_output_aver_bound){
        printf("\033[31m--------------------------\033[0m");
        std::cout << "mpdg_aver reach precision " << to_string(precision_output_aver_bound);
        std::cout << ": total_pn_iterations= " << ALGparam.total_pn_iterations << ": total_cg_iterations= " << ALGparam.total_cg_iterations;
        printf("\033[31m--------------------------\033[0m \n");
        precision_output_aver_bound /=10;
        /*
        if(precision_output_aver_bound == 1e-4){
            for (L row = 0; row < m; ++row) {
                lambda_over_sigma[row] *= sigma;
            }
            sigma = param->alpha * param->kappa;
            for (L row = 0; row < m; ++row) {
                lambda_over_sigma[row] /= sigma;
            }

        }

        if(precision_output_aver_bound == 1e-7){
            for (L row = 0; row < m; ++row) {
                lambda_over_sigma[row] *= sigma;
            }
            sigma = param->alpha * param->kappa;
            for (L row = 0; row < m; ++row) {
                lambda_over_sigma[row] /= sigma;
            }

        }
        */
        
    }
}

template<typename L, typename D>
void PPA<L, D>::result_print() {
    Print();
    Print("AdapPPA result:");
    std::cout << std::setprecision(ALGparam.prec);
    std::cout << "sigma=" << sigma << std::endl;
    std::cout << "Outer_iter="<< AbsAlgorithmSub<L, D>::iter << ", PPA_iter=" << PPA_iter ;
    std::cout << ", time=" << PPA_time << std::endl;
    std::cout << "obj=" << AbsAlgorithmSub<L, D>::fun_value << ", gap=" << gap << std::endl;
}

#endif //ALGORITHM_PPA_H
