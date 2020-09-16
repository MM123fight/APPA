#ifndef PROBLEM_SCIMETOLP_CMD_LINE_H
#define PROBLEM_SCIMETOLP_CMD_LINE_H

void SCIME_exit_with_help(){

    std::cout << "Usage: ./YourAlgorithm (options) -t [type] [data_name]" << std::endl;
    std::cout << "-t SCIMEtoLP_type: 1 for SCIMEtoLP, 2 for SCIMEtoLP_inequ, 3 for SCIMEtoLPd_inequ" << std::endl;

    std::cout << "options:" << std::endl;
    std::cout << "-i The deleted index for SCIME:(default 0)" <<std::endl;
    std::cout << "-s scale: (default = nonscale)" << std::endl;
    std::cout << "-l delta: the tuning parameter(default 0.01)" << std::endl;

    exit(0);
}

void SCIME_parse_cmd_line(int argc, char** argv){

    int i;
    for(i=1;i<argc;i++){
        if( argv[i][0] != '-' )
            break;
        if( ++i >= argc )
            SCIME_exit_with_help();

        switch(argv[i-1][1]){
            case 'i': Probparam.idx = atoi(argv[i]);
                break;
            case 't': Probparam.LP_type = atoi(argv[i]);
                break;
            case 'l': Probparam.delta = atof(argv[i]);
                break;
            case 's': Probparam.scale = true;
                --i;
                break;
            default:
                std::cout << "unknown option: -" << argv[i-1][1] << std::endl;
                exit(0);
        }
    }

    if(i>=argc)
        SCIME_exit_with_help();
    Probparam.data_name = argv[i];
}

#endif //ALGORITHM_CMD_LINE_H
