#ifndef PROBLEM_PROBPARAM_H
#define PROBLEM_PROBPARAM_H

template <typename L, typename D>
class ProblemParam{

public:
    std::string data_dir = "../data";
    std::string data_name;
    L LP_type;

    L prec = 15;
    L k;
    L idx = 0;
    D lambda = 1.;
    D delta = 0.01;
    D epsilon_multiple = 0.001;
    bool svm_scale = true;
    bool scale = false;


};

ProblemParam<int,double> Probparam;


#endif //ALGORITHM_CMD_LINE_H
