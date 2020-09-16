#include "ProbData.h"
#include "NMFtoLP_cmd_line.h"

int main(int argc, char** argv) {
    if( argc < 1+2 ){
        NMF_exit_with_help();
    }
    NMF_parse_cmd_line( argc, argv );

    std::string root_directory = Probparam.data_dir;
    std::string experiment = Probparam.data_name;
    ProbData<int, double> *inst = new ProbData<int, double>(root_directory, experiment);

    int NMFtoLP_type = Probparam.LP_type;
    double epsilon = Probparam.epsilon_multiple * inst->m;
    if(NMFtoLP_type == 1) {
        inst->NMFtoLP(epsilon);
        Print("Type: NMFtoLP");
    }else if(NMFtoLP_type == 2){
        inst->NMFtoLP_inequ(epsilon);
        Print("Type: NMFtoLP_inequ");
    }else if(NMFtoLP_type == 3){
        inst->NMFtoLPd_inequ(epsilon);
        Print("Type: NMFtoLPd_inequ");
    }else{
        Print("No such Type");
    }
    inst->dimension_print();

    delete inst;
    return 0;
}