#ifndef USE_MATRIX_H
#define USE_MATRIX_H
#include "stdafx.h"

template <typename L, typename D>
void Mtranspose(std::vector<L> &MT_row_ptr, std::vector<L> &MT_col_idx,std::vector<D> &MT_value, const L &n, const L &d, const std::vector<L> &M_row_ptr, const std::vector<L> &M_col_idx,
                const std::vector<D> &M_value);
template <typename L,typename D>
inline D l2norm_square(const std::vector<D> &x);
template <typename L, typename D>
inline D l2norm_square(const L &n, const std::vector<D> &x);

//Know the information of A and calculate the informmation of A^T
template <typename L, typename D>
void Mtranspose(std::vector<L> &MT_row_ptr, std::vector<L> &MT_col_idx,std::vector<D> &MT_value, const L &n, const L &d, const std::vector<L> &M_row_ptr, const std::vector<L> &M_col_idx,
                const std::vector<D> &M_value){

    MT_value.clear();
    MT_col_idx.clear();
    MT_row_ptr.clear();


    std::vector<std::vector<L>> tmp_idx(d);
    std::vector<std::vector<D>> tmp_value(d);
    for(L i = 0; i < n; ++i){
        for(L j = M_row_ptr[i]; j < M_row_ptr[i+1];++j){
            tmp_idx[M_col_idx[j]].push_back(i);
            tmp_value[M_col_idx[j]].push_back(M_value[j]);
        }
    }

    MT_row_ptr.resize(d+1,0);
    for(L i = 0; i < d; ++i){
        MT_value.insert(MT_value.end(), tmp_value[i].begin(), tmp_value[i].end());
        MT_col_idx.insert(MT_col_idx.end(), tmp_idx[i].begin(), tmp_idx[i].end());
        MT_row_ptr[i+1] = MT_col_idx.size();
    }

}

template <typename L,typename D>
inline D l2norm_square(const std::vector<D> &x){
    D tmp = 0;
    for (L row = 0;  row< x.size(); ++row)
        tmp += x[row]*x[row];
    return tmp;
};

template <typename L, typename D>
inline D l2norm_square(const L &n, const std::vector<D> &x){
    D tmp = 0;
    for (L row = 0;  row< n; ++row){
        tmp += x[row]*x[row];
    }
    return tmp;
};

#endif //USE_MATRIX_H
