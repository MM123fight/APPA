#ifndef PROBLEM_SVMTOLP_CMD_LINE_H
#define PROBLEM_SVMTOLP_CMD_LINE_H

void SVM_exit_with_help(){

    std::cout << "Usage: ./YourAlgorithm (options) -k [#classes] -t [SVMtoLP_type] [data_name]" << std::endl;
    std::cout << "-k #classes: The number of classes for SVM" <<std::endl;
    std::cout << "-t SVMtoLP_type: 1 for SVMtoLP, 2 for SVMtoLP_inequ, 3 for SVMtoLPd_inequ" << std::endl;

    std::cout << "options:" << std::endl;
    std::cout << "-l lambda: The parameter for SVM penalty(default:1.)" << std::endl;
    std::cout << "-v svm-scale:(default = true)" << std::endl;
    std::cout << "-s scale: (default = nonscale)" << std::endl;

    exit(0);
}

void SVM_parse_cmd_line(int argc, char** argv){

    int i;
    for(i=1;i<argc;i++){
        if( argv[i][0] != '-' )
            break;
        if( ++i >= argc )
            SVM_exit_with_help();

        switch(argv[i-1][1]){
            case 'v': Probparam.svm_scale = false;
                --i;
                break;
            case 's': Probparam.scale = true;
                --i;
                break;
            case 'k': Probparam.k = atoi(argv[i]);
                break;
            case 'l': Probparam.lambda = atof(argv[i]);
                break;
            case 't': Probparam.LP_type = atoi(argv[i]);
                break;
            default:
                std::cout << "unknown option: -" << argv[i-1][1] << std::endl;
                exit(0);
        }
    }

    if(i>=argc)
        SVM_exit_with_help();
    Probparam.data_name = argv[i];
}

#endif //ALGORITHM_CMD_LINE_H
