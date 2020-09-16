#include "ProbData.h"

int main(int argc, char** argv) {
    if( argc < 1+1 ){
	std::cout << "LPtoLP_inequ data_name" << std::endl;
        exit(0);
    }
    Probparam.data_name = argv[1];
    std::string root_directory = Probparam.data_dir;
    std::string experiment = Probparam.data_name;
    ProbData<int, double> *inst = new ProbData<int, double>(root_directory, experiment, "LP");
   
    inst->LPtoLP_inequ();
    inst->split_A();
    inst->result_path = root_directory + "/" + experiment + "_inequ";
    //inst->scs_data_log();
    inst->data_log();
    
    delete inst;
    return 0;
}
