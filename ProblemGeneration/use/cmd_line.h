#ifndef ALGORITHM_CMD_LINE_H
#define ALGORITHM_CMD_LINE_H

template <typename L, typename D>
class AlgorithmParam{

public:

    std::string data_name; 
    std::string data_dir = "/data/zhengqu";
    bool solve_from_dual;
    bool solve_with_scale;
    bool res_with_scale_norm;
    bool PCG = false;
    D tol;
    D tol_trans;
    L max_iter,blocksize;
    D beta;
    L sub_method = 1;
    //The following vector is only used for the convenience when I give the code for PN_CG_T;
    std::vector<L> test_A_col_idx;
    unsigned int prec = 6;
    unsigned int print_prec = 6;
    unsigned int total_cg_iterations = 0;
    unsigned int total_pn_iterations = 0;
    D scale_column_norm;
    //type = 0:
    //type = 1:
    //type = 2:
    unsigned int type = 1;
    D rho;
    D exp_idx_kappa;
    D varepsilon;
    D deltac;
    AlgorithmParam(){
        rho = 0.7;
        exp_idx_kappa = 5;
	deltac = 0.9;
        varepsilon = 1e16;
        solve_from_dual = false;
        solve_with_scale = false;
        res_with_scale_norm = false;
		tol = 1e-6;
        max_iter = 1000000;
        blocksize = 1;
        tol_trans = 1e16;
        beta = 1;
        scale_column_norm = 1.;
    }
};

AlgorithmParam<unsigned int, double> ALGparam;

void exit_with_help(){

    std::cout << "Usage: ./YourAlgorithm (options) [data_name]" << std::endl;
    std::cout << "options:" << std::endl;
    std::cout << "-d solve_dual: solve dual to obtain solution of primal problem (default: No)" <<std::endl;
    std::cout << "-t tol: tolerance of primal and dual infeasibility for terminiation (default: 1e-3)" << std::endl;
    std::cout << "-a tol_trans: tolerance of approx to PN(default: 1e16)" << std::endl;
    std::cout << "-z sub_size: block size of the randomized method for subproblem (default: 1)" << std::endl;
    std::cout << "-m max_iter: maximum number of outer iterations (default 1000)" << std::endl;
    std::cout << "-b beta: The penalty parameter for LADMM" << std::endl;
    std::cout << "-s scale: solve with scaling" << std::endl;
    std::cout << "-w way: the method type to solve the subproblem" << std::endl;
    std::cout << "-r res: the residual with corresponding b_norm and c_norm before scaling or after scaling;" << std::endl;
    std::cout << "-c column_norm: the column_norm for scaling(default 1)" << std::endl;
    std::cout << "-p type: type for min(default aver_min)" << std::endl;
    std::cout << "-h rho: default 0.7" << std::endl;
    std::cout << "-e exp_idx_kappa: default 5" << std::endl;
    std::cout << "-f PCG: default false" << std::endl;
    std::cout << "-v varepsilon: default 1" << std::endl;
    std::cout << "-g deltac: default 0.7" << std::endl;

    exit(0);
}

void parse_cmd_line(int argc, char** argv){

    int i;
    for(i=1;i<argc;i++){
        if( argv[i][0] != '-' )
            break;
        if( ++i >= argc )
            exit_with_help();

        switch(argv[i-1][1]){

            case 'd': ALGparam.solve_from_dual = true;
                --i;
                break;
		    case 's': ALGparam.solve_with_scale = true;
                --i;
                break;
            case 'r': ALGparam.res_with_scale_norm = true;
                --i;
                break;
            case 'f': ALGparam.PCG = false;
                --i;
                break;
            case 't': ALGparam.tol = atof(argv[i]);
                break;
            case 'a': ALGparam.tol_trans = atof(argv[i]);
                break;
            case 'h': ALGparam.rho = atof(argv[i]);
                break;
            case 'e': ALGparam.exp_idx_kappa = atof(argv[i]);
                break;
            case 'z': ALGparam.blocksize = atoi(argv[i]);
                break;
            case 'p': ALGparam.type = atoi(argv[i]);
                break;
            case 'm': ALGparam.max_iter = atoi(argv[i]);
                break;
            case 'b': ALGparam.beta = atof(argv[i]);
                break;
            case 'w': ALGparam.sub_method = atoi(argv[i]);
                break;
            case 'c': ALGparam.scale_column_norm= atof(argv[i]);
                break;
            case 'v': ALGparam.varepsilon= atof(argv[i]);
                break;
	    case 'g': ALGparam.deltac = atof(argv[i]);
                break;
            default:
                std::cout << "unknown option: -" << argv[i-1][1] << std::endl;
                exit(0);
        }
    }

    if(i>=argc)
        exit_with_help();
    ALGparam.data_name = argv[i];
}

#endif //ALGORITHM_CMD_LINE_H
