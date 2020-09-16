#ifndef ALGORITHM_ABSALGORITHM_H
#define ALGORITHM_ABSALGORITHM_H

#include "../Problem/ProblemHeader.h"
#include "../LossFun/LossFunHeader.h"
#include "../use/useheader.h"
#include "cmd_line.h"

template<typename L, typename D>
class AbsAlgorithm{
protected:
    ProbData<L, D>* data = NULL;
    AbsLoss<L, D>* loss = NULL;
    //random algorithm
    D mu;
    std::vector<D> v;
    std::vector<D> Lip;
    unsigned int loss_type;
    unsigned int stop_type;
    unsigned int reg_type;
    std::ofstream logFile;
    std::string data_file_path;
    std::string full_result_path;
    bool result_iflog;
    D reg_alpha = 1.;
    L K,blocksize;
    L max_iter;
    D precision;
    D grad_norm_infty;
    std::string method_type;
    D Lip_initial_value;
    D Lip_value;
public:
    L num_of_batch_per_check_for_gradnorm = 1;
    D initial_grad_norm,grad_norm;
    D initial_fun_value, fun_value;
    D primal_dual_gap;
    std::vector<D> grad;
    std::vector<D> w;//res_primal;
    std::vector<D> x;
    std::vector<D> lambda;
    std::vector<D> w_dual;
    std::vector<L> real_A_row_idx;
    std::vector<D> res_grad;
    L A_effect_nnz;
    L iter;
    D time;
    D fun_value_opt = 0;
    D left,right;
    AbsAlgorithm(ProbData<L, D>* const data_inst, const unsigned int& loss_type, const unsigned int& reg_type,
                 const std::string& data_file_path=""):
            data(data_inst),loss_type(loss_type), reg_type(reg_type), data_file_path(data_file_path){
        std::cout << std::setprecision(ALGparam.prec);
        Lip_initial_value = l2norm_square<L, D>(data->A_value);
        res_grad = std::vector<D>(data->n);

    }
    AbsAlgorithm(ProbData<L, D>* const data_inst, const std::string& data_file_path=""):
            data(data_inst), data_file_path(data_file_path){
    }

    virtual ~AbsAlgorithm() {
        delete loss;
        if(logFile.is_open()) {
            logFile.close();
        }
    }

    //Set loss functions and  regularization functions
    void set_loss();
    void set_loss(const std::vector<D>&b);
    void set_loss(const D& alpha, const std::vector<D>&b, const std::vector<D>& d);
    void set_reg(const D& alpha);

    // f(x) + g(x) ~ ESO(v); K: restart period for APPROX
    void check_blocksize(const L& blocksize);
    void Krcd_compute(L& K, const L& blocksize);
    void set_v_mu(const L& blocksize);


    //Computation of grad and fun value
    void grad_norm_update();
    virtual void left_value_compute();
    void reg_cord_cal(std::vector<D> &x,std::vector<D>& delta_x, const D& coefficient ,
                      const std::vector<D> &cord_grad, const L& blocksize, const std::vector<L>& cord);
    void res_grad_compute();
    void res_fun_value_compute(D& fun_value);

    //Linear System Solver
    virtual void LS_active_solver(const D& sub_alpha, const std::vector<L> &active_index, const std::vector<D> &w,
                                                   const std::vector<D> &grad_J, bool result_iflog = false){}

    //PPA functions
    void PPA_sub_solver_initial(D& sub_alpha, const D& sigma, const std::vector<D>& x_initial,
                                const std::vector<D> &lambda_over_sigma_initial);
    void PPA_left_value_compute();
    void PPA_right_update(D& right, const D& sigma, const D& delta, const std::vector<D>& x_initial,
                          const std::vector<D>& lambda_over_sigma_initial);
    virtual void PPA_solver(const D& sigma,const std::vector<D> &x_initial, const std::vector<D> &lambda_over_sigma_initial, const L &blocksize,
                            const D &varepsilon, const D &delta, bool result_iflog){}
    virtual void lp_solver(D beta, const std::vector<D> &x_initial, bool result_iflog = false){}

    //log functions
    void logfile_open();
    void result_log();
};

template<typename L, typename D>
void AbsAlgorithm<L, D>::set_loss() {
    loss = new AbsLoss<L, D>();
}

template<typename L, typename D>
void AbsAlgorithm<L, D>::set_loss(const std::vector<D>&b) {
    if(loss == NULL)
        loss = new squared_loss<L, D>(b);
    else{
        loss->set_b(b);
    }
}

template<typename L, typename D>
void AbsAlgorithm<L, D>::set_loss(const D& alpha, const std::vector<D>&b, const std::vector<D>& d) {
    if(loss == NULL) {
        switch (loss_type) {
            case 1:
                loss = new plus_squared_loss<L, D>(alpha, b, d);
                break;
            default:
                Print("There is no such stop type!");
                break;
        }
    } else{
        loss->set_b(b);
        loss->set_d(d);
        loss->alpha = alpha;
    }
}

template<typename L, typename D>
void AbsAlgorithm<L, D>::set_reg(const D& alpha) {
    reg_alpha = alpha;
}

template<typename L, typename D>
void AbsAlgorithm<L, D>::check_blocksize(const L& blocksize) {
    if (blocksize > data->n) {
        Print("The size of block is larger than n! Please choose a new block size");
        exit(0);
    }
}

template<typename L, typename D>
void AbsAlgorithm<L, D>::Krcd_compute(L& K, const L& blocksize) {
    if(mu > 0)
        K = ceil(2 * exp(1) * data->n * (sqrt(1 + 1. / mu) - 1) / (D) blocksize + 1);
    else {
        Print("Please use restart APPROX for non-strongly convex problem");
    }
}

template<typename L, typename D>
void AbsAlgorithm<L, D>::set_v_mu(const L& blocksize){
    loss->v_set(v,data,blocksize);
    //nested dependent names need "typename" to verify.
    typename std::vector<D>::iterator maxPosition = std::max_element(v.begin(),v.end());
    mu = loss->alpha / (*maxPosition);
}

template<typename L, typename D>
void AbsAlgorithm<L, D>::grad_norm_update() {
    D tmp;
    grad_norm = 0;
    grad_norm_infty = 0;
    for (L row = 0; row < data->n; ++row) {
        tmp = abs(res_grad[row]);
        if( tmp > grad_norm_infty )
            grad_norm_infty = tmp;
        grad_norm +=  tmp*tmp;
    }
    grad_norm = sqrt(grad_norm);
}

template<typename L, typename D>
void AbsAlgorithm<L, D>::left_value_compute() {
    switch (stop_type) {
        case 1:
            loss->grad_compute(grad,data, w, x);
            res_grad_compute();
            grad_norm_update();
            left = grad_norm;
            //calculate funtion value
            loss->fun_value_compute(fun_value, data, w, x);
            res_fun_value_compute(fun_value);
            break;
        case 2:
            loss->fun_value_compute(fun_value, data, w, x);
            res_fun_value_compute(fun_value);
            left = fun_value - fun_value_opt;
            break;
        default:
            Print("There is no such stop type!");
            break;
    }
}

template<typename L, typename D>
void AbsAlgorithm<L, D>::reg_cord_cal(std::vector<D> &x,std::vector<D>& delta_x, const D& coefficient ,
                                      const std::vector<D> &cord_grad, const L& blocksize, const std::vector<L>& cord) {
    L tmp_cord;
    switch (reg_type) {
        case 1://indicator fun reg
            for (L row = 0; row < blocksize ; ++row) {
                tmp_cord = cord[row];
                delta_x[row] = -cord_grad[row] / (coefficient * v[tmp_cord]);
                if ((tmp_cord < data->nplus) && (x[tmp_cord] + delta_x[row] < 0)) {
                    delta_x[row] = -x[tmp_cord];
                    x[tmp_cord] = 0;
                } else {
                    x[tmp_cord] += delta_x[row];
                }
            }
            break;
        case 2://L1 reg
            D tmp_coef, tmp;
            //Print("check 2");
            for (L row = 0; row < blocksize ; ++row) {
                tmp_cord = cord[row];
                tmp_coef = coefficient * v[tmp_cord];
                tmp = tmp_coef * x[tmp_cord] - cord_grad[row];
                if (tmp < -1) {
                    delta_x[row] = (-cord_grad[row] + 1) / tmp_coef;
                    x[tmp_cord] = (tmp + 1) / tmp_coef;
                } else if (tmp > 1) {
                    delta_x[row] = (-cord_grad[row] - 1) / tmp_coef;
                    x[tmp_cord] = (tmp - 1) / tmp_coef;
                } else {
                    delta_x[row] = -x[tmp_cord];
                    x[tmp_cord] = 0;
                }
            }
            break;
        case 3://no regularization
            for (L row = 0; row < blocksize ; ++row) {
                tmp_cord = cord[row];
                delta_x[row] = -cord_grad[row] / (coefficient * v[tmp_cord]);
                x[tmp_cord] += delta_x[row];
            }
            break;
        default:
            Print("There is no such reg type!");
            break;
    }
}

template<typename L, typename D>
void AbsAlgorithm<L, D>::res_grad_compute() {
    res_grad = grad;
    switch (reg_type) {
        case 1:
            for (L row = 0; row < data->nplus; ++row) {
                if ((x[row] == 0) && (res_grad[row] > 0))
                    res_grad[row] = 0;
            }
            break;
        case 2:
            for (L row = 0; row < data->n; ++row) {
                if (x[row] > 0)
                    res_grad[row] += reg_alpha;
                else if(x[row] < 0)
                    res_grad[row] -= reg_alpha;
                else{
                    if(res_grad[row] > 1)
                        res_grad[row] -= 1;
                    else if(grad[row] < -1)
                        res_grad[row] += 1;
                    else
                        res_grad[row] = 0;
                }
            }
            break;
        case 3:
            break;
        default:
            Print("There is no such reg type!");
            break;
    }
}

template<typename L, typename D>
void AbsAlgorithm<L, D>::res_fun_value_compute(D& fun_value){
    switch (reg_type) {
        case 1:
            break;
        case 2:
            break;
        case 3:
            for (L row = 0; row < data->n; ++row)
                fun_value += reg_alpha * abs(x[row]);
            break;
        default:
            Print("There is no such reg type!");
            break;
    }

}

template<typename L, typename D>
void AbsAlgorithm<L, D>::PPA_sub_solver_initial(D& sub_alpha, const D& sigma, const std::vector<D>& x_initial,
                                                const std::vector<D> &lambda_over_sigma_initial) {
    L m = data->m;
    L n = data->n;
    sub_alpha = 1/(sigma * sigma);
    //Lip_value = Lip_initial_value + sub_alpha * n; //It seems that this one is for APPROX
    std::vector<D> sub_b(m);
    std::vector<D> sub_d(n);
    for (L row = 0; row < m; ++row) {
        sub_b[row] = data->b[row] - lambda_over_sigma_initial[row];
    }
    for (L row = 0; row < n; ++row) {
        sub_d[row] = x_initial[row] - sigma*data->c[row];
    }
    set_loss(sub_alpha,sub_b,sub_d);
    /*

    Print("--------------------------------------------");
    Print("--------------------------------------------");
    Print("sigma", sigma);
    Print("sub_alpha", sub_alpha);
    Print("sub_b", sub_b);
    Print("sub_d", sub_d);
     */
    x = x_initial;
    loss->w_initial(w,data,x,loss->b);
    /*
    Print("initial_x",x);
    Print("lambda_over_sigma_initial", lambda_over_sigma_initial);
    Print("intial_w",w);
     */
    loss->grad_compute(grad,data,w,x);
    //Print("initial_grad", grad);
    res_grad_compute();
    //Print("intial_res_grad",grad);
    loss->fun_value_compute(fun_value,data,w,x);
    grad_norm_update();
    initial_grad_norm = grad_norm;
    initial_fun_value = fun_value ;
    left = grad_norm;

    //data->data_print();
    //Print("--------------------------------------------");
    //Print("--------------------------------------------");

    //PPA_left_value_compute();
}

template<typename L, typename D>
void AbsAlgorithm<L, D>::PPA_left_value_compute() {
    res_grad_compute();
    grad_norm_update();
    left = grad_norm;
    loss->fun_value_compute(fun_value, data, w, x);
}

template<typename L, typename D>
void AbsAlgorithm<L, D>::PPA_right_update(D& right, const D& sigma, const D& delta, const std::vector<D>& x_initial,
                                          const std::vector<D>& lambda_over_sigma_initial) {
    right = 0;
    D tmp;
    for (L row = 0; row < data->n; ++row) {
        tmp = (x[row] - x_initial[row])/sigma;
        right += tmp*tmp;
    }

    for (L row = 0; row < data->mi; ++row) {
        tmp = w[row];
        if(tmp < 0) {
            tmp = 0;
        }
        tmp -= lambda_over_sigma_initial[row];
        right += tmp * tmp;
    }

    for (L row = data->mi; row < data->m; ++row) {
        tmp = w[row]- lambda_over_sigma_initial[row];
        right += tmp * tmp;
    }
    right = sqrt(right);
    right *= (delta/sigma);
}

template<typename L, typename D>
void AbsAlgorithm<L, D>::logfile_open() {
    logFile.open(full_result_path.c_str());
    if (logFile.fail()) {
        Print("!!! Cannot open experiment result file: ", full_result_path);
        exit(0);
    }
}

template<typename L, typename D>
void AbsAlgorithm<L, D>::result_log(){
    logFile << std::setprecision(ALGparam.prec) << iter << "\t" << fun_value << "\t" << left << "\t" << time << std::endl;
    std::cout << "iter= " << iter << ", obj=" << fun_value << ", grad_norm=" << left;
    std::cout << ", time = " << time << std::endl;
};

#endif //ALGORITHM_ABSALGORITHM_H
