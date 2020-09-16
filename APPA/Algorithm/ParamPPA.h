#ifndef ALGORITHM_PARAMPPA_H
#define ALGORITHM_PARAMPPA_H

#include "../use/useheader.h"
#include "../Problem/ProblemHeader.h"

template<typename L, typename D>
inline void exp_update(D& varepsilon, const D& varepsilon_initial, const D& exp_idx, const L& iter);

template<typename L, typename D>
inline void inv_update(D& varepsilon, const D& varepsilon_initial, const D& inv_idx, const L& iter);

template<typename D>
class ParamPPA {
public:
    //exp_idx < 1; cur_idx >1;
    //initial varepsilon and the parameters for updating varepsilon
    //varepsilon *= exp_idx_varepsilon;
    //D varepsilon = 1e10;
    D varepsilon = 1.;
    D exp_idx_varepsilon = 0.9;
    D inv_idx_varepsilon = 1.1;

    //initial delta and the parameters for updating delta
    D delta;

    /* kappa: initial guess for kappa of max(Hoffman constant, growth condition parameter)
     * sigma = alpha * kappa: the parameter for PPA, i.e, P = (I+sigma*T)^{-1}
     * exp_idx_sigma: we also can say exp_idx_kappa, every time we double kappa, then we also double sigma
     * sigma = alpha * kappa;
     */
    D kappa;
    D alpha;
    D sigma;
    D exp_idx_kappa;

    D C, rho;

    ParamPPA(const D& rho, const D& delta, const D& sigma, const D& exp_idx_kappa, const D& varepsilon){
        this->rho = rho;
        this->delta = delta;
        //alpha = sqrt(pow(1./(rho - (2+rho)*delta),2)-1);
        alpha = sqrt(pow((1+delta)/(rho*(1-delta) - delta),2)-1);
        this->sigma = sigma;
        this->exp_idx_kappa = exp_idx_kappa;
        this->varepsilon = varepsilon;
        //sigma = alpha * kappa;
        kappa = sigma/alpha;
        C = (1+delta)/(1-delta)/(1-1./sqrt(1+alpha*alpha));
    }
    ~ParamPPA(){}

    void param_set(const D& rho, const D& exp_idx_kappa){
        this->rho = rho;
        delta = 0.9*rho/(rho+1);
        alpha = sqrt(pow((1+delta)/(rho*(1-delta) - delta),2)-1);
        C = (1+delta)/(1-delta)/(1-1./sqrt(1+alpha*alpha));
        this->exp_idx_kappa = exp_idx_kappa;
    }

    void printparam(){
        printf("\033[31mParameters of PPA:\033[0m \n");
        std::cout << std::setprecision(4);
        std::cout << "varepsilon_initial="<< varepsilon << ", exp_idx_varepsilon=" << exp_idx_varepsilon;
        std::cout << ", inv_idx_varepsilon=" << inv_idx_varepsilon << std::endl;
        std::cout << "sigma_initial=" << sigma << ", kappa_initial=" << kappa << std::endl;
        std::cout << "rho="<< rho << ", delta="<< delta << ", alpha=" << alpha << ", C=" << C;
        std::cout << ", exp_idx_kappa=" << exp_idx_kappa << std::endl;
    }
};

template<typename L, typename D>
inline void exp_update(D& varepsilon, const D& varepsilon_initial, const D& exp_idx, const L& iter){
    varepsilon = varepsilon_initial * pow(exp_idx,iter);
}

template<typename L, typename D>
inline void inv_update(D& varepsilon, const D& varepsilon_initial, const D& inv_idx, const L& iter){
    varepsilon = varepsilon_initial/pow(iter+1,inv_idx);
}


#endif //ALGORITHM_PARAMPPA_H
