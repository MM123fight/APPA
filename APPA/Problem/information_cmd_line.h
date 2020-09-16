#ifndef PROBLEM_INFORMATION_CMD_LINE_H
#define PROBLEM_INFORMATION_CMD_LINE_H

void information_exit_with_help(){
    std::cout << "Usage: ./YourAlgorithm (options)  [data_name]" << std::endl;
    exit(0);
}

void information_parse_cmd_line(int argc, char** argv){

    int i;
    for(i=1;i<argc;i++){
        if( argv[i][0] != '-' )
            break;
        if( ++i >= argc )
            information_exit_with_help();
    }

    if(i>=argc)
        information_exit_with_help();
    Probparam.data_name = argv[i];
}

#endif //ALGORITHM_CMD_LINE_H
