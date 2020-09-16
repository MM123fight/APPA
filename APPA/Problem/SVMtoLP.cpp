#include "ProbData.h"
#include "SVMtoLP_cmd_line.h"

int main(int argc, char** argv) {
    if( argc < 1+3 ){
        SVM_exit_with_help();
    }
    SVM_parse_cmd_line( argc, argv );

    std::string root_directory = Probparam.data_dir;
    std::string experiment = Probparam.data_name;
    ProbData<int, double> *inst = new ProbData<int, double>(root_directory, experiment);

    double lambda = Probparam.lambda;
    int k = Probparam.k;
    int SVMtoLP_type = Probparam.LP_type;
    if(SVMtoLP_type == 1) {
        inst->SVMtoLP(lambda, k);
        Print("Type: SVMtoLP");
    }else if(SVMtoLP_type == 2){
        inst->SVMtoLP_inequ(lambda,k);
        Print("Type: SVMtoLP_inequ");
    }else if(SVMtoLP_type == 3){
        inst->SVMtoLPd_inequ(lambda,k);
        Print("Type: SVMtoLPd_inequ");
    }else{
        Print("No such Type");
    }
    inst->dimension_print();

    delete inst;
    return 0;
}