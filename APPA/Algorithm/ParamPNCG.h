#ifndef ALGORITHM_PARAMPNCG_H
#define ALGORITHM_PARAMPNCG_H

#include "../use/useheader.h"
#include "../Problem/ProblemHeader.h"


template<typename L, typename D>
class ParamPNCG {
public:
    //exp_idx < 1; cur_idx >1;
    //initial varepsilon and the parameters for updating varepsilon
    //varepsilon *= exp_idx_varepsilon;
    D mu = 0.01;
    D eta = 0.1;
    D tau = 0.2;
    D beta = 0.5;
    D enlargement = 100;
    L max_num_linesearch = 50;

    L pn_max_iter;
    bool reach_linesearch = false;

    ParamPNCG(){
    }
    ~ParamPNCG(){}

    void printparam(){
        printf("\033[31mParameters of PN_CG:\033[0m \n");
        std::cout << "mu="<< mu << ", eta=" << eta << ", tau=" << tau << ", beta=" << beta;
        std::cout << ", enlargement=" << enlargement;
        std::cout << ", max_num_linesearch=" << max_num_linesearch << std::endl;
    }
};
ParamPNCG<unsigned int, double> PNCGparam;



#endif //ALGORITHM_PARAMPNCG_H
