#ifndef PROBDATA_H
#define PROBDATA_H

#include "LPuse.h"
#include "../Algorithm/cmd_line.h"

/* (P)LP problem:
 * min <c,x>   s.t.
 * Ae*x = be
 * Ai*x <= bi
 * x_[nb] >= 0
 * where A:m*n, b:R^m, x:R^n, c:R^n
 */
/* we compare APPA with scs: https://github.com/cvxgrp/scs
 * and LPsparse:http://ianyen.site/LPsparse/
 */

template <typename L, typename D>
class ProbData{
public:

    L mi,me,m;
    L nb,nf,n;
    L nplus,mplus;

    //Ai: mi*n
    std::vector<L> Ai_row_ptr = std::vector<L>();
    std::vector<L> Ai_col_idx = std::vector<L>();
    std::vector<D> Ai_value = std::vector<D>();

    // Ae: me*n
    std::vector<L> Ae_row_ptr = std::vector<L>();
    std::vector<L> Ae_col_idx = std::vector<L>();
    std::vector<D> Ae_value = std::vector<D>();

    // A = [Ai;Ae], m*n
    std::vector<L> A_row_ptr = std::vector<L>();
    std::vector<L> A_col_idx = std::vector<L>();
    std::vector<D> A_value = std::vector<D>();

    // AT: n*m
    std::vector<L> AT_row_ptr = std::vector<L>();
    std::vector<L> AT_col_idx = std::vector<L>();
    std::vector<D> AT_value = std::vector<D>();

    // bi: R^mi, be:R^me
    // b = [bi;be], R^m
    std::vector<D> bi = std::vector<D>();
    std::vector<D> be = std::vector<D>();
    std::vector<D> b = std::vector<D>();

    // c: R^n
    std::vector<D> c = std::vector<D>();
    std::vector<D> Ai_norm_square = std::vector<D>();
    D A_norm_F,b_norm,c_norm;
    D max_Ai;

    std::string data_type;
    std::string data_file_path;
    std::string result_path;
    L N = 999;

    ProbData(const std::string &root_path_name, const std::string &data_file_name, std::string data_type = "");
    virtual ~ProbData(){}

    void Ai_clear();
    void Ae_clear();

    // Split A into Ai and Ae; b into bi and be.
    // Combine Ai and Ae into A; bi and be into b;
    void split_A();
    void combine_A();

    // ai: the i-th row of matrix A
    // Ai: the i-th column of matrix A
    void ai_norm_square_compute(std::vector<D>& ai_norm);
    void Ai_norm_square_compute(std::vector<D>& Ai_norm);
    void A_norm_F_compute(D& A_norm_F);
    void max_Ai_compute(D& max_Ai);

    void scale();
    void scale_sum_one();
    void svm_scale();
    void LP_scale();
    void constant_scale();

    // LPtoLPd: obtains the dual LP to (1)
    // LPtoLP_inequ: Ae*x = be, Ai*x <= bi, x[nb] >=0 -> Ae*x = be, Ai*x <= bi
    // LP_iequ2LP_only_inequ: Ae*x = be, Ai*x <= bi -> Ai*x <= bi
    void LPtoLPd();
    void LPtoLP_inequ();
    void LP_iequ2LP_only_inequ();
    //XXXtoLP: for general LP as (P)
    //XXXtoLP_inequ: for standard LP as: min <c,x> s.t. Ax <= b;
    //XXtoLPd_inequ: the standard LP of the dual LP to (P).
    void SVMtoLP_prepare(const D& lambda, const L& k);
    void SVMtoLP(const D& lambda, const L& k);
    void SVMtoLP_inequ(const D& lambda, const L& k);
    void SVMtoLPd_inequ(const D& lambda, const L& k);

    void SCIMEtoLP_prepare(const L& idx, const D& delta);
    void SCIMEtoLP(const L &idx, const D& delta);
    void SCIMEtoLP_inequ(const L &idx, const D& delta);
    void SCIMEtoLPd_inequ(const L &idx, const D& delta);

    void NMFtoLP_prepare(const D& epsilon);
    void NMFtoLP(const D& epsilon);
    void NMFtoLP_inequ(const D& epsilon);
    void NMFtoLPd_inequ(const D& epsilon);

    //information_log: check the information for SVM data before transforming into LP.
    //scs_data_log: log the data form for scs;
    //data_log:log the data form for APPA and LPsparse;
    void information_log();
    void scs_data_log();
    void data_log();

    //If APPA want to solve the dual LP to (P), then transform the data firstly with solve_from_dual();
    void solve_from_dual();
    void dimension_print();
    void data_print();
};

template <typename L, typename D>
ProbData<L, D>::ProbData(const std::string &root_path_name, const std::string &data_file_name, std::string data_type) {
    std::cout << std::setprecision(Probparam.prec);
    this->data_type = data_type;
    data_file_path = root_path_name + "/"+ data_file_name;
    if(data_type == "") {
        //Read from SVM data file, and then you can transform the data to the corresponding LP.
        //For example: data/arcene2/arcene2
        MatrixRead(data_file_path+"/" + data_file_name, m, n, A_row_ptr, A_col_idx, A_value, b);
    } else if(data_type == "LP"){
        //Read data from LP form, for example as data set: data/arcene2/SVMtoLP_lambda1
        //Read mi,me,nb,nf;
        readMeta(data_file_path+"/meta", mi, me, nb, nf);
        n = nb+nf;
        //Read Ai,Ae,bi,be,c;
        readMat(data_file_path+"/A", mi, n,Ai_row_ptr, Ai_col_idx, Ai_value);
        readMat(data_file_path+"/Aeq", me, n,Ae_row_ptr, Ae_col_idx, Ae_value);
        VectorRead(data_file_path+"/b",mi,bi);
        VectorRead(data_file_path+"/beq",me,be);
        VectorRead(data_file_path+"/c",n,c);
        //A=[Ai;Ae]; b = [bi;be];
        combine_A();
        Mtranspose(AT_row_ptr,AT_col_idx,AT_value,m,n,A_row_ptr,A_col_idx,A_value);
        Ai_norm_square_compute(Ai_norm_square);
        A_norm_F_compute(A_norm_F);
        max_Ai_compute(max_Ai);
        b_norm = sqrt(l2norm_square<L, D>(b));
        c_norm = sqrt(l2norm_square<L, D>(c));
    }else{
        Print("There is no such data type");
        exit(0);
    }
}

template<typename L, typename D>
void ProbData<L, D>::Ai_clear() {
    Ai_row_ptr.clear();
    Ai_col_idx.clear();
    Ai_value.clear();
    bi.clear();
}

template<typename L, typename D>
void ProbData<L, D>::Ae_clear() {
    Ae_row_ptr.clear();
    Ae_col_idx.clear();
    Ae_value.clear();
    be.clear();
}

template<typename L, typename D>
void ProbData<L, D>::split_A() {
    L nnz_Ai = A_row_ptr[mi];
    Ai_row_ptr = std::vector<L>(A_row_ptr.begin(),A_row_ptr.begin()+mi+1);
    Ai_col_idx = std::vector<L>(A_col_idx.begin(),A_col_idx.begin()+nnz_Ai);
    Ai_value = std::vector<D>(A_value.begin(),A_value.begin()+nnz_Ai);

    Ae_row_ptr = std::vector<L>(A_row_ptr.begin()+mi+1,A_row_ptr.end());
    for (L row = 0; row < Ae_row_ptr.size(); ++row) {
        Ae_row_ptr[row] -= nnz_Ai;
    }
    Ae_row_ptr.insert(Ae_row_ptr.begin(),0);
    Ae_col_idx = std::vector<L>(A_col_idx.begin()+nnz_Ai,A_col_idx.end());
    Ae_value = std::vector<D>(A_value.begin()+nnz_Ai,A_value.end());

    bi = std::vector<D>(b.begin(),b.begin()+mi);
    be = std::vector<D>(b.begin()+mi,b.end());
}

template<typename L, typename D>
void ProbData<L, D>::combine_A() {
    m = mi+me;
    nplus = nb;
    mplus = mi;

    A_row_ptr = Ai_row_ptr;
    A_col_idx = Ai_col_idx;
    A_value = Ai_value;
    b = bi;
    A_col_idx.insert(A_col_idx.end(), Ae_col_idx.begin(),Ae_col_idx.end());
    A_value.insert(A_value.end(), Ae_value.begin(), Ae_value.end());
    L tmp = A_row_ptr.back();
    for (L row = 1; row < me+1; ++row)
        A_row_ptr.push_back(Ae_row_ptr[row]+tmp);
    b.insert(b.end(),be.begin(),be.end());
}

template<typename L, typename D>
void ProbData<L, D>::ai_norm_square_compute(std::vector<D>& ai_norm) {
    ai_norm = std::vector<D>(m,0);
    D tmp;
    D tmp_value;
    for (L row = 0; row < m; ++row) {
        tmp = 0;
        for (L col = A_row_ptr[row]; col < A_row_ptr[row+1]; ++col) {
            tmp_value = A_value[col];
            tmp += tmp_value * tmp_value;
        }
        ai_norm[row] = tmp;
    }

}

template<typename L, typename D>
void ProbData<L, D>::Ai_norm_square_compute(std::vector<D>& Ai_norm) {
    Ai_norm = std::vector<D>(n,0);
    D tmp;
    D tmp_value;
    for (L row = 0; row < n; ++row) {
        tmp = 0;
        for (L col = AT_row_ptr[row]; col < AT_row_ptr[row+1]; ++col) {
            tmp_value = AT_value[col];
            tmp += tmp_value * tmp_value;
        }
        Ai_norm[row] = tmp;
    }
}

template<typename L, typename D>
void ProbData<L, D>::A_norm_F_compute(D &A_norm_F) {
    A_norm_F = 0;
    for (L row = 0; row < n; ++row) {
        A_norm_F += Ai_norm_square[row];
    }
    A_norm_F = sqrt(A_norm_F);
}

template<typename L, typename D>
void ProbData<L, D>::max_Ai_compute(D &max_Ai) {
    max_Ai = 0;
    for (L row = 0; row < n; ++row) {
        if(max_Ai < Ai_norm_square[row])
            max_Ai = Ai_norm_square[row];
    }
    max_Ai = sqrt(max_Ai);
}



template<typename L, typename D>
void ProbData<L, D>::scale() {
    Mtranspose(AT_row_ptr,AT_col_idx,AT_value,m,n,A_row_ptr,A_col_idx,A_value);
    D norm,tmp;
    std::vector<D> std_var(n);
    for (L row = 0; row < n; ++row) {
        norm = 0;
        for (L col = AT_row_ptr[row]; col < AT_row_ptr[row + 1]; ++col) {
            tmp = AT_value[col];
            norm += tmp*tmp;
        }
        norm = sqrt(norm);
        std_var[row] = norm;
        for (L col = AT_row_ptr[row]; col < AT_row_ptr[row + 1]; ++col) {
            AT_value[col] /= norm;
        }
    }
    Mtranspose(A_row_ptr,A_col_idx,A_value,n,m,AT_row_ptr,AT_col_idx,AT_value);
    /*
    D b_norm = sqrt(l2norm_square<L, D>(b));
    for (L row = 0; row < b.size(); ++row) {
        b[row] /= b_norm;
    }
    D c_norm = sqrt(l2norm_square<L, D>(c));
    for (L row = 0; row < c.size(); ++row) {
        c[row] /= c_norm;
    }
     */
}

template<typename L, typename D>
void ProbData<L, D>::scale_sum_one() {
    for (L row = 0; row < A_value.size(); ++row) {
        D noise = rand()%(N+1)/(D)(N+1);
        noise -= 0.5;
        noise *= 1e-3;
        A_value[row] += noise;
    }
    bool neg = false;
    Mtranspose(AT_row_ptr,AT_col_idx,AT_value,m,n,A_row_ptr,A_col_idx,A_value);
    D sum,tmp;
    for (L row = 0; row < n; ++row) {
        sum = 0;
        for (L col = AT_row_ptr[row]; col < AT_row_ptr[row + 1]; ++col) {
            tmp = AT_value[col];
            if(tmp < 0){
                AT_value[col] = - tmp;
                neg = true;
            }
            sum += abs(tmp);
        }
        for (L col = AT_row_ptr[row]; col < AT_row_ptr[row + 1]; ++col) {
            AT_value[col] /= sum;
        }
    }
    Mtranspose(A_row_ptr,A_col_idx,A_value,n,m,AT_row_ptr,AT_col_idx,AT_value);
    if(neg == true){
        Print("You have neg entries");
    }


}

template<typename L, typename D>
void ProbData<L, D>::svm_scale() {
    std::vector<L> AT_row_ptr_scale;
    std::vector<L> AT_col_idx_scale;
    std::vector<D> AT_value_scale;
    std::vector<L> AT_row_ptr_scale_tmp(n+1,0);
    std::vector<L> AT_col_idx_scale_tmp;
    std::vector<D> AT_value_scale_tmp;
    Mtranspose(AT_row_ptr_scale,AT_col_idx_scale,AT_value_scale,m,n,A_row_ptr,A_col_idx,A_value);
    std::vector<D> mean(n);
    std::vector<D> std_var(n);
    D mean_val,norm,tmp, row_ptr_tmp;
    for (L row = 0; row < n; ++row) {
        mean_val = 0;
        for (L col = AT_row_ptr_scale[row]; col < AT_row_ptr_scale[row + 1]; ++col) {
            mean_val += AT_value_scale[col];
        }
        mean_val /= m;
        norm = 0;
        row_ptr_tmp = 0;
        for (L col = AT_row_ptr_scale[row]; col < AT_row_ptr_scale[row + 1]; ++col) {
            tmp = AT_value_scale[col] - mean_val;
            if(tmp != 0) {
                ++row_ptr_tmp;
                norm += tmp * tmp;
                AT_value_scale_tmp.push_back(tmp);
                AT_col_idx_scale_tmp.push_back(AT_col_idx_scale[col]);
            }
        }
        AT_row_ptr_scale_tmp[row+1] = AT_row_ptr_scale_tmp[row] + row_ptr_tmp;
        std_var[row] = sqrt(norm);
        mean[row] = mean_val;
    }
    for (L row = 0; row < n; ++row) {
        tmp = std_var[row];
        for (L col = AT_row_ptr_scale_tmp[row]; col < AT_row_ptr_scale_tmp[row + 1]; ++col) {
            AT_value_scale_tmp[col] /= tmp;
        }
    }
    Mtranspose(A_row_ptr,A_col_idx,A_value,n,m,AT_row_ptr_scale_tmp,AT_col_idx_scale_tmp,AT_value_scale_tmp);
    VectorWrite(data_file_path + "/mean", mean);
    VectorWrite(data_file_path + "/std_var", std_var);
}

template<typename L, typename D>
void ProbData<L, D>::LPtoLPd() {
    AT_row_ptr = A_row_ptr;
    AT_col_idx = A_col_idx;
    AT_value = A_value;
    for (L row = 0; row < AT_value.size(); ++row) {
        AT_value[row] = -AT_value[row];
    }
    std::vector<D> tmp_vec;
    tmp_vec = c;
    c = b;
    b = tmp_vec;
    Mtranspose(A_row_ptr, A_col_idx, A_value, m, n, AT_row_ptr, AT_col_idx, AT_value);
    L tmp;
    tmp = m, m = n, n = tmp;
    tmp = mi, mi = nb, nb = tmp;
    tmp = me, me = nf, nf = tmp;
    mplus = mi, nplus = nb;
}

template<typename L, typename D>
void ProbData<L, D>::LPtoLP_inequ(){
    std::vector<D> zeros(nb,0);
    std::vector<D> ones(nb,-1);
    std::vector<L> ones_idx(nb);
    for (L row = 0; row < nb; ++row)
        ones_idx[row] = row;
    A_col_idx.insert(A_col_idx.begin(),ones_idx.begin(),ones_idx.end());
    A_value.insert(A_value.begin(),ones.begin(),ones.end());
    b.insert(b.begin(),zeros.begin(),zeros.end());
    for (L row = 0; row < m+1; ++row)
        A_row_ptr[row] += nb;
    A_row_ptr.insert(A_row_ptr.begin(),ones_idx.begin(),ones_idx.end());
    mi += nb, nb = 0, nf = n, m = mi+me, nplus = nb, mplus = mi;
    Mtranspose(AT_row_ptr,AT_col_idx,AT_value,m,n,A_row_ptr,A_col_idx,A_value);
}

template<typename L, typename D>
void ProbData<L, D>::LP_iequ2LP_only_inequ() {
    L nnz_Ai = A_row_ptr[mi];
    L nnz = A_row_ptr[m];
    L increase = nnz - nnz_Ai;
    Ae_row_ptr = std::vector<L>(A_row_ptr.begin()+mi+1,A_row_ptr.end());
    for (L row = 0; row < Ae_row_ptr.size(); ++row) {
        Ae_row_ptr[row] += increase;
    }
    Ae_col_idx = std::vector<L>(A_col_idx.begin()+nnz_Ai,A_col_idx.end());
    std::vector<D> neg_Ae_value = std::vector<D>(A_value.begin()+nnz_Ai,A_value.end());
    for (L row = 0; row < neg_Ae_value.size(); ++row) {
        neg_Ae_value[row] = -neg_Ae_value[row];
    }

    A_row_ptr.insert(A_row_ptr.end(),Ae_row_ptr.begin(),Ae_row_ptr.end());
    A_col_idx.insert(A_col_idx.end(),Ae_col_idx.begin(),Ae_col_idx.end());
    A_value.insert(A_value.end(),neg_Ae_value.begin(),neg_Ae_value.end());

    std::vector<D> neg_be = std::vector<D>(b.begin()+mi,b.end());
    for (L row = 0; row < neg_be.size(); ++row) {
        neg_be[row] = - neg_be[row];
    }
    b.insert(b.end(),neg_be.begin(),neg_be.end());

    m += me;
    mi = m;
    me = 0;
    mplus = m;
}

template<typename L, typename D>
void ProbData<L, D>::SVMtoLP_prepare(const D& lambda, const L& k) {
    if(Probparam.svm_scale == true) {
        if(A_row_ptr.back()/(D)(m*n) == 1){
            svm_scale();
        }else {
            scale();
        }
        information_log();
        Print("The data is svm_scaled");
    }
    std::vector<D> sign = class_sign<L, D>(b);
    if(sign.size()!=k){
        Print("The number of classes is not accordant with k");
        exit(0);
    }
    Ai_clear();
    me = 0;
    mi = (k-1) * m;
    Ai_row_ptr = std::vector<L>(mi + 1);
    std::vector<D> neg_A_value(A_value.size());
    for (L row = 0; row < A_value.size(); ++row) {
        neg_A_value[row] = - A_value[row];
    }
   // Print("Ai_row_ptr_size", Ai_row_ptr.size());
    Ai_row_ptr[0] = 0;
    for (L row = 0; row < m; ++row) {
        for (L col = 0; col < k-1; ++col) {
            Ai_row_ptr[row*(k-1) + col+1] =  Ai_row_ptr[row*(k-1) + col] + 4*(A_row_ptr[row+1] - A_row_ptr[row])+1;
        }
    }
    L start, end;
    L size;
    std::vector<L> fix_idx, fix_idx_later;
    std::vector<L> variant_idx;
    std::vector<L> idx;
    L class_cord_b;
    for (L row = 0; row < m; ++row) {
        start = A_row_ptr[row];
        end = A_row_ptr[row+1];
        idx.assign(A_col_idx.begin()+start, A_col_idx.begin()+end);
        fix_idx = idx;
        size = end - start;
        class_cord_b = check_class(b[row],k,sign);
        for (L col = 0; col < size; ++col) {
            fix_idx[col] += (class_cord_b-1)*n;
        }
        fix_idx_later = fix_idx;
        for (L col = 0; col < size; ++col) {
            fix_idx_later[col] += k*n;
        }
        for (L i = 0;  i < class_cord_b-1; ++i) {
            Ai_value.insert(Ai_value.end(), A_value.begin() + start, A_value.begin() + end);
            Ai_value.insert(Ai_value.end(), neg_A_value.begin() + start, neg_A_value.begin() + end);
            Ai_value.insert(Ai_value.end(), neg_A_value.begin() + start, neg_A_value.begin() + end);
            Ai_value.insert(Ai_value.end(), A_value.begin() + start, A_value.begin() + end);
            Ai_value.push_back(-1.);
            variant_idx = idx;
            for (L col = 0; col < size; ++col) {
                variant_idx[col] += i*n;
            }
            Ai_col_idx.insert(Ai_col_idx.end(),variant_idx.begin(), variant_idx.end());
            Ai_col_idx.insert(Ai_col_idx.end(),fix_idx.begin(),fix_idx.end());
            for (L col = 0; col < size; ++col) {
                variant_idx[col] += k*n;
            }
            Ai_col_idx.insert(Ai_col_idx.end(),variant_idx.begin(),variant_idx.end());
            Ai_col_idx.insert(Ai_col_idx.end(),fix_idx_later.begin(),fix_idx_later.end());
            Ai_col_idx.push_back(2*k*n+row);
        }
        for (L i = class_cord_b; i < k; ++i) {
            Ai_value.insert(Ai_value.end(), neg_A_value.begin() + start, neg_A_value.begin() + end);
            Ai_value.insert(Ai_value.end(), A_value.begin() + start, A_value.begin() + end);
            Ai_value.insert(Ai_value.end(), A_value.begin() + start, A_value.begin() + end);
            Ai_value.insert(Ai_value.end(), neg_A_value.begin() + start, neg_A_value.begin() + end);
            Ai_value.push_back(-1.);

            Ai_col_idx.insert(Ai_col_idx.end(),fix_idx.begin(),fix_idx.end());
            variant_idx = idx;
            for (L col = 0; col < size; ++col) {
                variant_idx[col] += i*n;
            }
            Ai_col_idx.insert(Ai_col_idx.end(),variant_idx.begin(),variant_idx.end());
            Ai_col_idx.insert(Ai_col_idx.end(),fix_idx_later.begin(),fix_idx_later.end());
            for (L col = 0; col < size; ++col) {
                variant_idx[col] += k*n;
            }
            Ai_col_idx.insert(Ai_col_idx.end(),variant_idx.begin(), variant_idx.end());
            Ai_col_idx.push_back(2*k*n+row);
        }

    }
    c = std::vector<D>(2*k*n,lambda);
    std::vector<D> ones_vec(m,1);
    c.insert(c.end(),ones_vec.begin(),ones_vec.end());
    bi = std::vector<D>(mi,-1.);
    n = 2*k*n+m,nb = n,nplus = n, nf = 0;
    m = mi,mplus = m, me = 0;
    A_row_ptr = Ai_row_ptr;
    A_col_idx = Ai_col_idx;
    A_value = Ai_value;
    b = bi;
    //Print(Ai_row_ptr);
    //Print(mi,n,Ai_row_ptr,Ai_col_idx,Ai_value);
    if(Probparam.scale == true){
        scale();
        Print("The data is std_scaled");
    }
}

template<typename L, typename D>
void ProbData<L, D>::SVMtoLP(const D& lambda, const L& k) {
    SVMtoLP_prepare(lambda,k);
    split_A();
    //Mtranspose(AT_row_ptr,AT_col_idx,AT_value,m,n,A_row_ptr,A_col_idx,A_value);
    result_path =  data_file_path + "/SVMtoLP";
    data_log();
}

template<typename L, typename D>
void ProbData<L, D>::SVMtoLP_inequ(const D& lambda, const L& k) {
    SVMtoLP_prepare(lambda,k);
    LPtoLP_inequ();
    split_A();
    result_path =  data_file_path + "/SVMtoLP_inequ";
    //data_log();
    scs_data_log();
}

template<typename L, typename D>
void ProbData<L, D>::SVMtoLPd_inequ(const D& lambda, const L& k) {
    SVMtoLP_prepare(lambda,k);
    LPtoLPd();
    LPtoLP_inequ();
    split_A();
    result_path =  data_file_path + "/SVMtoLPd_inequ";
    //data_log();
    scs_data_log();
}

template<typename L, typename D>
void ProbData<L, D>::SCIMEtoLP_prepare(const L &idx, const D& delta) {
    L p_n = m;
    L p_d = n;
    if(idx >= p_d){
        Print("Your index exceed over p_d");
        exit(0);
    }
    scale();
    std::vector<D> mean(p_d);
    D tmp_mean;
    for (L row = 0; row < p_d; ++row) {
        tmp_mean = 0;
        for (L col = AT_row_ptr[row]; col < AT_row_ptr[row+1]; ++col) {
            tmp_mean += AT_value[col];
        }
        mean[row] = tmp_mean / (D)sqrt((D)p_n);
    }
    D mean_i = mean[idx];
    mean.erase(mean.begin()+idx);
    L start = AT_row_ptr[idx];
    L end = AT_row_ptr[idx+1];
    std::vector<D> Ai(p_n,0);
    for (L row = start; row < end; ++row) {
        Ai[AT_col_idx[row]] =  AT_value[row];
    }
    AT_col_idx.erase(AT_col_idx.begin()+start,AT_col_idx.begin()+end);
    AT_value.erase(AT_value.begin()+start,AT_value.begin()+end);
    L size = AT_row_ptr[idx+1] - AT_row_ptr[idx];
    for (L row = idx+1; row < p_d+1; ++row) {
        AT_row_ptr[row] -= size;
    }
    AT_row_ptr.erase(AT_row_ptr.begin() + idx);

    p_d -= 1;
    std::vector<D> bi(p_d);
    D bi_tmp;
    for (L row = 0; row < p_d; ++row) {
        bi_tmp = 0;
        for (L col = AT_row_ptr[row]; col < AT_row_ptr[row+1]; ++col) {
            bi_tmp += AT_value[col] * Ai[AT_col_idx[col]];
        }
        bi[row] = bi_tmp-mean_i * mean[row];
    }

    mi = 2*p_d;
    me = p_n+1;
    nb = 2*p_d;
    nf = p_n+1;
    m = mi + me;
    n = nb + nf;
    //obtain b and c;
    c = std::vector<D>(n,0);
    for (L row = 0; row < nb; ++row) {
        c[row] = 1;
    }
    b = std::vector<D>(m,0);
    for (L row = 0; row < p_d; ++row) {
	 b[row] = bi[row] + delta;
    }
    for (L row = 0; row < p_d; ++row) {
        b[row+p_d] = -bi[row] + delta;
    }
    //obtain matrix
    Mtranspose(A_row_ptr,A_col_idx,A_value,p_d,p_n,AT_row_ptr,AT_col_idx,AT_value);
    std::vector<L> A1T_row_ptr = A_row_ptr;
    std::vector<L> A1T_col_idx = A_col_idx;
    std::vector<D> A1T_value = A_value;
    A1T_row_ptr.push_back(A1T_row_ptr.back()+p_d);
    for (L row = 0; row < p_d; ++row) {
        A1T_col_idx.push_back(row);
        A1T_value.push_back(-mean[row]);
    }
    std::vector<L> A2_row_ptr = A1T_row_ptr;
    std::vector<L> A2_col_idx = A1T_col_idx;
    std::vector<D> A2_value = A_value;
    A2_value.insert(A2_value.end(),mean.begin(),mean.end());

    std::vector<L> A1_row_ptr;
    std::vector<L> A1_col_idx;
    std::vector<D> A1_value;
    Mtranspose(A1_row_ptr,A1_col_idx,A1_value,p_n+1,p_d,A1T_row_ptr,A1T_col_idx,A1T_value);


    std::vector<D> neg_A1_value(A1_value.size());
    std::vector<D> neg_A2_value(A2_value.size());
    for (L row = 0; row < A1_value.size(); ++row) {
        neg_A1_value[row] = - A1_value[row];
    }
    for (L row = 0; row < A2_value.size(); ++row) {
        neg_A2_value[row] = - A2_value[row];
    }

    A1_value.insert(A1_value.end(),neg_A1_value.begin(), neg_A1_value.end());
    for (L row = 0; row < p_n+1; ++row) {
        start = A2_row_ptr[row];
        end = A2_row_ptr[row+1];
        A1_value.insert(A1_value.end(),A2_value.begin()+start, A2_value.begin()+end);
        A1_value.insert(A1_value.end(),neg_A2_value.begin()+start, neg_A2_value.begin()+end);
        A1_value.push_back(-1);
    }

    for (L row = 0; row < A1_col_idx.size(); ++row) {
        A1_col_idx[row] += mi;
    }
    std::vector<L> A2_col_idx_tmp = A2_col_idx;
    A1_col_idx.insert(A1_col_idx.end(),A1_col_idx.begin(),A1_col_idx.end());
    for (L row = 0; row < A2_col_idx_tmp.size(); ++row) {
        A2_col_idx_tmp[row] += p_d;
    }
    for (L row = 0; row < p_n+1; ++row) {
        start = A2_row_ptr[row];
        end = A2_row_ptr[row+1];
        A1_col_idx.insert(A1_col_idx.end(),A2_col_idx.begin()+start, A2_col_idx.begin()+end);
        A1_col_idx.insert(A1_col_idx.end(),A2_col_idx_tmp.begin()+start, A2_col_idx_tmp.begin()+end);
        A1_col_idx.push_back(mi+row);
    }

    std::vector<L> A1_row_ptr_tmp = A1_row_ptr;
    L nnz = A1_row_ptr.back();
    for (L row = 0; row < A1_row_ptr_tmp.size(); ++row) {
        A1_row_ptr_tmp[row] += nnz;
    }
    A1_row_ptr.insert(A1_row_ptr.end(),A1_row_ptr_tmp.begin()+1,A1_row_ptr_tmp.end());
    nnz *= 2;
    for (L row = 1; row < p_n+2; ++row) {
       A2_row_ptr[row] = 2 *  A2_row_ptr[row] + row + nnz;
    }
    A1_row_ptr.insert(A1_row_ptr.end(), A2_row_ptr.begin()+1, A2_row_ptr.end());
    m = me + mi;
    n = nb + nf;
    mplus = mi, nplus = nb;

    A_row_ptr = A1_row_ptr;
    A_col_idx = A1_col_idx;
    A_value = A1_value;
     
}


template<typename L, typename D>
void ProbData<L, D>::SCIMEtoLP(const L &idx, const D& delta) {
    SCIMEtoLP_prepare(idx, delta);
    //Mtranspose(AT_row_ptr, AT_col_idx, AT_value, m, n, A_row_ptr, A_col_idx, A_value);
    split_A();
    result_path = data_file_path + "/SCIMEtoLP";
    data_log();
}

template<typename L, typename D>
void ProbData<L, D>::SCIMEtoLP_inequ(const L &idx, const D& delta) {
    SCIMEtoLP_prepare(idx, delta);
    LPtoLP_inequ();
    LP_iequ2LP_only_inequ();
    Mtranspose(AT_row_ptr,AT_col_idx,AT_value,m,n,A_row_ptr,A_col_idx,A_value);
    result_path = data_file_path + "/SCIMEtoLP_inequ";
    //data_log();
    scs_data_log();
}

template<typename L, typename D>
void ProbData<L, D>::SCIMEtoLPd_inequ(const L &idx, const D& delta) {
    SCIMEtoLP_prepare(idx, delta);
    LPtoLPd();
    LPtoLP_inequ();
    LP_iequ2LP_only_inequ();
    Mtranspose(AT_row_ptr,AT_col_idx,AT_value,m,n,A_row_ptr,A_col_idx,A_value);
    result_path = data_file_path + "/SCIMEtoLPd_inequ";
    //data_log();
    scs_data_log();
}

template<typename L, typename D>
void ProbData<L, D>::NMFtoLP_prepare(const D& epsilon) {
    L p_n = m;
    L p_d = n;
    scale_sum_one();

    // Ae * x = be;
    std::vector<L> M_row_ptr = AT_row_ptr;
    std::vector<L> M_col_idx = AT_col_idx;
    std::vector<D> M_value = AT_value;
    std::vector<D> ones_2(2);
    std::vector<L> ones_2_idx(2);
    ones_2[0] = 1; ones_2[1] =-1;
    ones_2_idx[0] = p_n * p_n + p_d;
    ones_2_idx[1] = p_n * p_n + 2*p_d;
    L ptr_tmp;
    for (L row = p_d; row > 0; --row) {
        ones_2_idx[0] -= 1;
        ones_2_idx[1] -= 1;
        ptr_tmp = AT_row_ptr[row];
        M_value.insert(M_value.begin()+ptr_tmp,ones_2.begin(),ones_2.end());
        M_col_idx.insert(M_col_idx.begin()+ptr_tmp,ones_2_idx.begin(),ones_2_idx.end());
    }
    for (L row = 0; row < p_d; ++row) {
        M_row_ptr[row+1] += 2*(row+1);
    }

    Ae_row_ptr = M_row_ptr;
    Ae_col_idx = M_col_idx;
    Ae_value = M_value;
    L M_length = M_col_idx.size();
    std::vector<L> M_row_ptr_tmp = M_row_ptr;
    for (L i = 0; i < p_n-1; ++i) {
        Ae_value.insert(Ae_value.end(),M_value.begin(),M_value.end());

        for (L row = 0; row < p_d; ++row) {
            for (L col = M_row_ptr[row]; col < M_row_ptr[row+1] - 2; ++col) {
                M_col_idx[col] += p_n;
            }
            for (L col = M_row_ptr[row+1]-2; col < M_row_ptr[row+1] ; ++col) {
                M_col_idx[col] += 2*p_d;
            }
        }
        Ae_col_idx.insert(Ae_col_idx.end(),M_col_idx.begin(),M_col_idx.end());

        for (L row = 0; row < M_row_ptr_tmp.size(); ++row) {
            M_row_ptr_tmp[row] += M_length;
        }
        Ae_row_ptr.insert(Ae_row_ptr.end(),M_row_ptr_tmp.begin()+1,M_row_ptr_tmp.end());

    }
    n = p_n * p_n + 2 * p_d * p_n;
    me = p_d * p_n;
    be = std::vector<D>(me);
    L idx;
    for (L row = 0; row < p_n; ++row) {
        idx = row * p_d;
        for (L col = A_row_ptr[row]; col < A_row_ptr[row+1]; ++col) {
            be[idx + A_col_idx[col]] = A_value[col];
        }
    }
    std::vector<L> AiT_row_ptr(n+1,0);
    std::vector<L> AiT_col_idx;
    std::vector<D> AiT_value;
    std::vector<D> neg_ones_pn(p_n,-1);
    std::vector<D> neg_ones_pn_tmp;
    std::vector<L> idx_tmp(p_n);
    for (L row = 0; row < p_n; ++row) {
        idx_tmp[row] = row;
    }
    for (L i = 0; i < p_n; ++i) {
        neg_ones_pn_tmp = neg_ones_pn;
        neg_ones_pn_tmp[i] = 1;
        idx = i*p_n;
        for (L j = 0; j < p_n; ++j) {
            if(j == i){
                AiT_row_ptr[idx+j+1] =  AiT_row_ptr[idx+j] + p_n;
                AiT_value.insert(AiT_value.end(),neg_ones_pn_tmp.begin(),neg_ones_pn_tmp.end());
                AiT_col_idx.insert(AiT_col_idx.end(),idx_tmp.begin(), idx_tmp.end());
            }
            else{
                AiT_row_ptr[idx+j+1] = AiT_row_ptr[idx+j] + 1;
                AiT_value.push_back(1.);
                AiT_col_idx.insert(AiT_col_idx.end(),j*p_n+i);
            }

        }
        for (L row = 0; row < p_n; ++row) {
            idx_tmp[row] += p_n;
        }
    }
    for (L row = p_n*p_n; row < n; ++row) {
        AiT_row_ptr[row+1] = AiT_row_ptr[p_n*p_n];
    }
    mi = p_n * p_n;
    bi = std::vector<D>(mi, 0);
    for (L row = 0; row < p_n; ++row) {
        bi[row * p_n + row] = 1;
    }
    Mtranspose(Ai_row_ptr,Ai_col_idx,Ai_value,n,mi,AiT_row_ptr,AiT_col_idx,AiT_value);

    std::vector<D> ones_y(2*p_n*p_d,1);
    Ai_value.insert(Ai_value.end(),ones_y.begin(),ones_y.end());
    std::vector<D> idx_y(2*p_n*p_d);
    for (L row = 0; row < 2 * p_n * p_d; ++row) {
        idx_y[row] = p_n * p_n + row;
    }
    Ai_col_idx.insert(Ai_col_idx.end(),idx_y.begin(),idx_y.end());
    for (L row = 0; row < p_n; ++row) {
        Ai_row_ptr.push_back(Ai_row_ptr.back()+2*p_d);
        bi.push_back(epsilon);
    }
    mi += p_n;
    m = mi + me; mplus = mi;
    nb = n; nf = 0; nplus = 0;
    A_value =Ai_value;
    A_value.insert(A_value.end(),Ae_value.begin(),Ae_value.end());
    A_col_idx = Ai_col_idx;
    A_col_idx.insert(A_col_idx.end(),Ae_col_idx.begin(),Ae_col_idx.end());
    A_row_ptr = Ai_row_ptr;
    std::vector<L> Ae_row_ptr_tmp = Ae_row_ptr;
    L size = Ai_row_ptr.back();
    for (L row = 0; row < Ae_row_ptr.size(); ++row) {
        Ae_row_ptr_tmp[row] += size;
    }
    A_row_ptr.insert(A_row_ptr.end(),Ae_row_ptr_tmp.begin()+1,Ae_row_ptr_tmp.end());
    b = bi;
    b.insert(b.end(),be.begin(),be.end());
    c = std::vector<D>(n,0);
    for (L row = 0; row < p_n; ++row) {
            D noise = rand() % (N + 1) / (D) (N + 1);
            noise -= 0.5;
            noise *= 1e-3;
            c[row*p_n+row] = 1+noise;
    }
}

template<typename L, typename D>
void ProbData<L, D>::NMFtoLP(const D& epsilon) {
    NMFtoLP_prepare(epsilon);
    result_path = data_file_path + "/NMFtoLP";
    data_log();
}

template<typename L, typename D>
void ProbData<L, D>::NMFtoLP_inequ(const D& epsilon) {
    NMFtoLP_prepare(epsilon);
    LPtoLP_inequ();
    LP_iequ2LP_only_inequ();
    Mtranspose(AT_row_ptr,AT_col_idx,AT_value,m,n,A_row_ptr,A_col_idx,A_value);
    result_path = data_file_path + "/NMFtoLP_inequ";
    //data_log();
    scs_data_log();
}

template<typename L, typename D>
void ProbData<L, D>::NMFtoLPd_inequ(const D& epsilon) {
    NMFtoLP_prepare(epsilon);
    LPtoLPd();
    LPtoLP_inequ();
    LP_iequ2LP_only_inequ();
    Mtranspose(AT_row_ptr,AT_col_idx,AT_value,m,n,A_row_ptr,A_col_idx,A_value);
    result_path = data_file_path + "/NMFtoLPd_inequ";
    //data_log();
    scs_data_log();

}

template<typename L, typename D>
void ProbData<L, D>::information_log() {
    std::string information_path = data_file_path+"/information";
    std::ofstream fout;
    fout.open(information_path.c_str());
    if (fout.fail()) {
        std::cerr << "Cannot open data file: " << information_path << std::endl;
        exit(0);
    }
    printf("\033[31mThe SVM data:\033[0m \n");
    D percent = A_row_ptr.back()/(D)(m*n)*100;
    std::cout << "m=" << m << ", n=" << n << std::endl;
    std::cout << "nnz=" << A_row_ptr.back() << ", nnz_percent=" << std::setprecision(3) << percent << "%";
    std::cout << std::setprecision(Probparam.prec) << std::endl;

    fout << "m=" << m << ", n=" << n << std::endl;
    fout << "nnz=" << A_row_ptr.back() << ", nnz_percent=" << std::setprecision(3) << percent << "%";
    fout << std::setprecision(Probparam.prec) << std::endl;
    fout.close();
}

template<typename L, typename D>
void ProbData<L, D>::data_log() {
    writeMeta(result_path +"/meta",mi,me,nb,nf);
    VectorWrite(result_path + "/c",c);
    writeMat(result_path +"/A",mi,n,Ai_row_ptr,Ai_col_idx,Ai_value);
    writeMat(result_path +"/Aeq",me,n,Ae_row_ptr,Ae_col_idx,Ae_value);
    VectorWrite(result_path + "/b",bi);
    VectorWrite(result_path + "/beq",be);
}

template<typename L, typename D>
void ProbData<L, D>::scs_data_log() {
    std::string file_meta = result_path + "/scs_meta";
    char *path_meta = (char*)file_meta.c_str();
    scs_writeMeta(path_meta,mi,me,nb,nf);

    std::string file_b = result_path + "/scs_b";
    char *path_b  =  (char*)file_b.c_str();
    scs_writeVec(path_b,m,b);

    std::string file_c = result_path + "/scs_c";
    char *path_c  =  (char*)file_c.c_str();
    scs_writeVec(path_c,n,c);

    std::string file_AT_row_ptr = result_path + "/scs_AT_row_ptr";
    char *path_AT_row_ptr  =  (char*)file_AT_row_ptr.c_str();
    scs_writeVec(path_AT_row_ptr, n+1, AT_row_ptr);

    int nnz = AT_row_ptr.back();
    std::string file_AT_col_idx = result_path + "/scs_AT_col_idx";
    char *path_AT_col_idx  =  (char*)file_AT_col_idx.c_str();
    scs_writeVec(path_AT_col_idx, nnz, AT_col_idx);

    std::string file_AT_value = result_path + "/scs_AT_value";
    char *path_AT_value  =  (char*)file_AT_value.c_str();
    scs_writeVec(path_AT_value, nnz, AT_value);
}

template<typename L, typename D>
void ProbData<L, D>::solve_from_dual() {
    LPtoLPd();
}

template<typename L, typename D>
void ProbData<L, D>::LP_scale() {
    D norm,tmp;
    for (L row = 0; row < n; ++row) {
        norm = 0;
        for (L col = AT_row_ptr[row]; col < AT_row_ptr[row + 1]; ++col) {
            tmp = AT_value[col];
            norm += tmp*tmp;
        }
        norm = sqrt(norm);
        norm /= ALGparam.scale_column_norm;
        if(norm !=0) {
            for (L col = AT_row_ptr[row]; col < AT_row_ptr[row + 1]; ++col) {
                AT_value[col] /= norm;
            }
            c[row] /= norm;
        }
    }
    Mtranspose(A_row_ptr,A_col_idx,A_value,n,m,AT_row_ptr,AT_col_idx,AT_value);
    Ai_norm_square_compute(Ai_norm_square);
    A_norm_F_compute(A_norm_F);
    max_Ai_compute(max_Ai);
    split_A();
    result_path = data_file_path +"_scale";
    data_log();
}

template<typename L, typename D>
void ProbData<L, D>::constant_scale() {
    for (L row = 0; row < n; ++row) {
        for (L col = AT_row_ptr[row]; col < AT_row_ptr[row + 1]; ++col) {
            AT_value[col] *= ALGparam.scale_column_norm;
        }
        c[row] *= ALGparam.scale_column_norm;
    }
    Mtranspose(A_row_ptr,A_col_idx,A_value,n,m,AT_row_ptr,AT_col_idx,AT_value);
    Ai_norm_square_compute(Ai_norm_square);
    A_norm_F_compute(A_norm_F);
    max_Ai_compute(max_Ai);
}

template<typename L, typename D>
void ProbData<L, D>::dimension_print() {
    printf("\033[31mDimension of matrix A:\033[0m \n");
    std::cout << "m=" << m << ", mi=" << mi << ", me=" << me << ", mplus=" << mplus << std::endl;
    std::cout << "n=" << n << ", nb=" << nb << ", nf=" << nf << ", nplus=" << nplus << std::endl;
    std::cout << "nnz=" << A_row_ptr.back() << std::endl;
}



template <typename L, typename D>
void ProbData<L, D>::data_print(){
    Print("c",c);
    Print("b",b);
    Print("Matrix A:");
    Print(m,n,A_row_ptr,A_col_idx,A_value);
    Print("Matrix AT:");
    Print(n,m,AT_row_ptr,AT_col_idx,AT_value);
}

#endif //PROBDATA_H
