#include "ProbData.h"
#include "SCIMEtoLP_cmd_line.h"

int main(int argc, char** argv) {
    if( argc < 1+2 ){
        SCIME_exit_with_help();
    }
    SCIME_parse_cmd_line( argc, argv );

    std::string root_directory = Probparam.data_dir;
    std::string experiment = Probparam.data_name;

    ProbData<int, double> *inst = new ProbData<int, double>(root_directory, experiment);

    int idx = Probparam.idx;
    double delta = Probparam.delta;
    int SCIMEtoLP_type = Probparam.LP_type;
    if(SCIMEtoLP_type == 1) {
        inst->SCIMEtoLP(idx, delta);
        Print("Type: SCIMEtoLP");
    }else if(SCIMEtoLP_type == 2){
        inst->SCIMEtoLP_inequ(idx, delta);
        Print("Type: SCIMEtoLP_inequ");
    }else if(SCIMEtoLP_type == 3){
        inst->SCIMEtoLPd_inequ(idx, delta);
        Print("Type: SCIMEtoLPd_inequ");
    }else{
        Print("No such Type");
    }

    inst->dimension_print();

    delete inst;
    return 0;
}