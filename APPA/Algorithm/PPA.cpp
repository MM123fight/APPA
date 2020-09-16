#include "PPA.h"

int main(int argc, char** argv){

    if( argc < 1+1 ){
        exit_with_help();
    }
    parse_cmd_line( argc, argv );

    std::string root_directory = ALGparam.data_dir;
    std::string experiment = ALGparam.data_name;
    ProbData<unsigned int, double> *inst = new ProbData<unsigned int, double>(root_directory, experiment, "LP");

    if(inst->mi < inst->nb){
        ALGparam.solve_from_dual == true;
        inst->solve_from_dual();
    }
    

    if(ALGparam.solve_with_scale == true) {
        //Print("The data before scaling:");
        //inst->data_print();
        inst->LP_scale();
        //Print("The data after scaling:");
    }else{
        inst->constant_scale();
    }
    if(ALGparam.res_with_scale_norm == true){
        inst->b_norm = sqrt(l2norm_square<unsigned int, double>(inst->b));
        inst->c_norm = sqrt(l2norm_square<unsigned int, double>(inst->c));
    }
    //inst->data_print();

    //Set the parameters for PPA!
    double rho = ALGparam.rho;
    double delta = 0.9*rho/(rho+1);
    double alpha = sqrt(pow((1+delta)/(rho*(1-delta) - delta),2)-1);
    double C = (1+delta)/(1-delta)/(1-1./sqrt(1+alpha*alpha));
    double exp_idx_kappa = ALGparam.exp_idx_kappa;
    double kappa = 1./(inst->A_norm_F);
    double sigma = alpha * kappa;
    double varepsilon = ALGparam.varepsilon;

    //------------------------------------
    //double sigma = 0.1;

    ParamPPA<double>* param = new ParamPPA<double>(rho,delta,sigma,exp_idx_kappa,varepsilon);
    AbsAlgorithmSub<unsigned int, double> *method = new PPA<unsigned int, double>(inst,1,1,1);


    std::vector<double> x(inst->n, 0);
    std::vector<double> lambda(inst->m, 0);

    //Information and parameters print.
    inst->dimension_print();
    PNCGparam.printparam();
    param->printparam();

    method->solver(x,lambda,ALGparam.blocksize, ALGparam.max_iter, ALGparam.tol,param,true,"APPA");

    delete method;
    delete param;
    delete inst;
    return 0;
}
