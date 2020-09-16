#include "ProbData.h"
#include "information_cmd_line.h"

int main(int argc, char** argv) {
    if( argc < 1+1 ){
        information_exit_with_help();
    }
    information_parse_cmd_line( argc, argv );

    std::string root_directory = Probparam.data_dir;
    std::string experiment = Probparam.data_name;

    ProbData<int, double> *inst = new ProbData<int, double>(root_directory, experiment);
    inst->information_log();

    delete inst;
    return 0;
}