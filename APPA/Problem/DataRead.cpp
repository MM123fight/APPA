#include "ProbData.h"
#include "SVMtoLP_cmd_line.h"

int main(int argc, char** argv) {
    if( argc < 1+1 ){
        SVM_exit_with_help();
    }
    SVM_parse_cmd_line( argc, argv );

    std::string root_directory = Probparam.data_dir;
    std::string experiment = Probparam.data_name;
    ProbData<int, double> *inst = new ProbData<int, double>(root_directory, experiment, "LP");
    double A_norm_F;
    //inst->data_print();
    inst->A_norm_F_compute(A_norm_F);
    Print("A_norm_F_square", A_norm_F*A_norm_F);
	inst->LP_scale();
    inst->A_norm_F_compute(A_norm_F);
    Print("A_norm_F_square", A_norm_F*A_norm_F);
	//Print(inst->scale_matrix_R);
	//inst->data_print();
    delete inst;
    return 0;
}
