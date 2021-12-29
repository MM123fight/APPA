#ifndef USE_DATARW_H
#define USE_DATARW_H

#include "stdafx.h"
#include "Print.h"
#include "cmd_line.h"

template <typename L, typename D>
void MatrixRead(const std::string &full_data_path, L &n, L &d,
                std::vector<L> &A_row_ptr, std::vector<L> &A_col_idx, std::vector<D> &A_value, std::vector<D> &b);
template <typename L, typename D>
void VectorRead(const std::string &full_data_path, L &d, std::vector<D> &c);

template <typename D>
void VectorWrite(const std::string &full_data_path,const std::vector<D> &c, const unsigned int& precision = 15);
template <typename D>
void RowVecWrite(const std::string &full_data_path,const std::vector<D> &vec);

template <typename L, typename D>
void SvmWrite(const std::string& full_data_path, const std::vector<L> &row_ptr, 
    const std::vector<L> &col_idx, const std::vector<D> &value, const std::vector<D> & b);

void Makefile(const std::string &full_data_path);

template <typename L, typename D>
void MatrixRead(const std::string &full_data_path, L &n, L &d,
                std::vector<L> &A_row_ptr, std::vector<L> &A_col_idx, std::vector<D> &A_value, std::vector<D> &b){
    A_row_ptr.clear();
    A_col_idx.clear();
    A_value.clear();
    b.clear();
    n = 0;
    d = 0;

    std::ifstream fin(full_data_path.c_str());
    if (fin.fail()) {
        std::cout << "Cannot open data file: " << full_data_path << std::endl;
        exit(0);
    }
    Print("Loading Data...");
    std::string line;
    char x;
    L idx;
    D val;
    A_row_ptr.push_back((L)0);
    while (getline(fin, line)) {
        std::istringstream iss(line);
        iss >> val;
        b.push_back(val);
        while (iss >> idx >> x >> val) {
            A_col_idx.push_back(idx - 1);
            A_value.push_back(val);
        }
        d = std::max<L>(d, A_col_idx.back());
        A_row_ptr.push_back(A_col_idx.size());
        ++n;
    }
    ++d;
    fin.close();
    Print("Finish loading data from:", full_data_path);

};

template <typename L, typename D>
void VectorRead(const std::string &full_data_path, L &d, std::vector<D> &c){
    c.clear();
    std::ifstream fin(full_data_path.c_str());
    if (fin.fail()) {
        std::cout << "Cannot open data file of vector c for LP: " << full_data_path << std::endl;
        exit(0);
    }

    std::string line;
    char x;
    D val;
    while (getline(fin, line)) {
        std::istringstream iss(line);
        iss >> val;
        c.push_back(val);
    }
    fin.close();
    //Print("size", c.size());
    if(c.size() != d) {
        std::cout << "The size of vector from " << full_data_path << " is not accordant with m" << std::endl;
        exit(0);
    }
};

template <typename D>
void VectorWrite(const std::string &full_data_path,const std::vector<D> &c, const unsigned int& precision) {

    std::ofstream fout;
    fout.open(full_data_path.c_str());
    fout << std::setprecision(ALGparam.prec);
    for (int i = 0; i < c.size(); ++i) {
        fout << std::setprecision(precision) << c[i] << "\n";
    }
    fout.close();
}

template <typename D>
void RowVecWrite(const std::string &full_data_path,const std::vector<D> &vec) {

    std::ofstream fout;
    fout.open(full_data_path.c_str());
    fout << std::setprecision(ALGparam.prec);
    for (int i = 0; i < vec.size(); ++i) {
        fout << vec[i] << "\t";
    }
    fout.close();
}

template <typename L, typename D>
void SvmWrite(const std::string& full_data_path, const std::vector<L> &row_ptr, 
    const std::vector<L> &col_idx, const std::vector<D> &value, const std::vector<D> & b){
    std::ofstream fout;
    fout.open(full_data_path.c_str());
    fout << std::setprecision(ALGparam.prec);
    for(L row = 0; row < b.size(); ++row){
        fout << b[row] << " ";
        for(L col = row_ptr[row]; col < row_ptr[row+1]; ++col)
            fout << col_idx[col]+1 << ":" << value[col] << " ";
        fout << std::endl;
    }
    fout.close();
}

void Makefile(const std::string &full_data_path) {

    std::ofstream fout;
    fout.open(full_data_path.c_str());
    fout.close();
}

#endif //USE_DATARW_H
