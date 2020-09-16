#ifndef PROBLEM_NMFTOLP_CMD_LINE_H
#define PROBLEM_NMFTOLP_CMD_LINE_H

void NMF_exit_with_help(){

    std::cout << "Usage: ./YourAlgorithm (options) -t[type] [data_name]" << std::endl;
    std::cout << "-t NMFtoLP_type: 1 for NMFtoLP, 2 for NMFtoLP_inequ, 3 for NMFtoLPd_inequ" << std::endl;
    std::cout << "options:" << std::endl;
    std::cout << "-l epsilon_multiple: epsilon = epsilon_multiple * p_n" <<std::endl;
    exit(0);
}

void NMF_parse_cmd_line(int argc, char** argv){

    int i;
    for(i=1;i<argc;i++){
        if( argv[i][0] != '-' )
            break;
        if( ++i >= argc )
            NMF_exit_with_help();

        switch(argv[i-1][1]){
            case 't': Probparam.LP_type = atoi(argv[i]);
                break;
            case 'l': Probparam.epsilon_multiple = atof(argv[i]);
                break;
            default:
                std::cout << "unknown option: -" << argv[i-1][1] << std::endl;
                exit(0);
        }
    }

    if(i>=argc)
        NMF_exit_with_help();
    Probparam.data_name = argv[i];
}

#endif //ALGORITHM_CMD_LINE_H
