#ifndef SDNA_LPUSE_H
#define SDNA_LPUSE_H

#include "../use/useheader.h"
#include "ProbParam.h"

//Read from SVM form:
template <typename L, typename D>
void readMatrix(const std::string &full_data_path, const L &m, const L &n, std::vector<L> &A_row_ptr,
                std::vector<L> &A_col_idx, std::vector<D> &A_value );

//For LPsparse form:
template <typename L>
void readMeta(const std::string &full_data_path ,L &mi, L &me, L &nb, L &nf);
template <typename L>
void writeMeta(const std::string &full_data_path ,const L &mi, const L &me, const L &nb, const L &nf);
template <typename L, typename D>
void readMat(const std::string &full_data_path, const L &m, const L &n, std::vector<L> &A_row_ptr,
             std::vector<L> &A_col_idx, std::vector<D> &A_value );
template <typename L, typename D>
void writeMat(const std::string &full_data_path, const L &m, const L &n, const std::vector<L> &A_row_ptr,
              const std::vector<L> &A_col_idx, const std::vector<D> &A_value );
template <typename L>
void writeMat(const std::string &full_data_path, const L &n);

//for scs form:
void scs_writeMeta(char * filename , const int &mi, const  int &me, const int &nb, const int &nf);
void scs_writeVec(char * filename , const int &size, const std::vector<int> &vec);
void scs_writeVec(char * filename , const int &size, const std::vector<double> &vec);
double * scs_readVec_double(char * filename , const int &size);
int * scs_readVec_int(char * filename , const int &size);


//Check the number of class of SVM and change to the corresponding sign.
template <typename L, typename D>
L check_class(const D &y, const L &k, const std::vector<D> &sign);
template <typename L, typename D>
std::vector<D> class_sign( const std::vector<D> &y);


template <typename L, typename D>
void readMatrix(const std::string &full_data_path, const L &m, const L &n, std::vector<L> &A_row_ptr,
                std::vector<L> &A_col_idx, std::vector<D> &A_value ){

    L i,j;
    D val;
    std::ifstream fin(full_data_path.c_str());
    if( fin.fail() ){
        std::cerr << "fail to open " << full_data_path << std::endl;
        exit(0);
    }

    L row_idx = 1;
    L tmp = 0;
    A_row_ptr.push_back(0);
    while( !fin.eof() ){
        fin >> val;

        if( fin.eof() )
            break;
        if( i-1 >= m || j-1 >= n ){
            std::cerr << "index:" << "(" << i-1 << ", " << j-1 << ") out of bound when reading " << full_data_path << std::endl;
            exit(0);
        }
        A_col_idx.push_back(j-1);
        A_value.push_back(val);
        while(i != row_idx){
            A_row_ptr.push_back(tmp);
            ++row_idx;
        }
        ++tmp;
    }
    if(tmp != 0)
        A_row_ptr.push_back(tmp);
    fin.close();
}

template <typename L>
void readMeta(const std::string &full_data_path ,L &mi, L &me, L &nb, L &nf){
    std::ifstream fin(full_data_path.c_str());
    if (fin.fail()) {
        std::cerr << "Cannot open data file: " << full_data_path << std::endl;
        exit(0);
    }
    std::string tmp;
    fin >> tmp >> nb;
    fin >> tmp >> nf;
    fin >> tmp >> mi;
    fin >> tmp >> me;
    fin.close();
}

template <typename L>
void writeMeta(const std::string &full_data_path ,const L &mi, const L &me, const L &nb, const L &nf){
    std::ofstream fout;
    fout.open(full_data_path.c_str());
    if (fout.fail()) {
        std::cerr << "Cannot open data file: " << full_data_path << std::endl;
        exit(0);
    }
    fout << "nb" << " "<< nb << "\n";
    fout << "nf" << " "<< nf << "\n";
    fout << "mi" << " "<< mi << "\n";
    fout << "me" << " "<< me;
    fout.close();
}

template <typename L, typename D>
void readMat(const std::string &full_data_path, const L &m, const L &n, std::vector<L> &A_row_ptr,
             std::vector<L> &A_col_idx, std::vector<D> &A_value ){

    L i,j;
    D val;
    std::ifstream fin(full_data_path.c_str());
    if( fin.fail() ){
        std::cerr << "fail to open " << full_data_path << std::endl;
        exit(0);
    }
    fin >> i >> j >> val; //filter one line
    if( i != m || j != n ){
        std::cerr << "dimension in " << full_data_path << " does not match that in meta file" << std::endl;
        exit(0);
    }
    L row_idx = 1;
    L tmp = 0;
    A_row_ptr.push_back(0);
    while( !fin.eof() ){
        fin >> i >> j >> val;

        if( fin.eof() )
            break;
        if( i-1 >= m || j-1 >= n ){
            std::cerr << "index:" << "(" << i-1 << ", " << j-1 << ") out of bound when reading " << full_data_path << std::endl;
            exit(0);
        }
        A_col_idx.push_back(j-1);
        A_value.push_back(val);
        while(i != row_idx){
            A_row_ptr.push_back(tmp);
            ++row_idx;
        }
        ++tmp;
        fin.get(); // 读取最后的回车符
        if(fin.peek() == '\n') break;
    }
    if(tmp != 0)
        A_row_ptr.push_back(tmp);
    fin.close();
}

template <typename L, typename D>
void writeMat(const std::string &full_data_path, const L &m, const L &n, const std::vector<L> &A_row_ptr,
              const std::vector<L> &A_col_idx, const std::vector<D> &A_value ){

    L i,j;
    D val;
    std::ofstream fout;
    fout.open(full_data_path.c_str());
    fout << std::setprecision(Probparam.prec);
    if( fout.fail() ){
        std::cerr << "fail to open " << full_data_path << std::endl;
        exit(0);
    }
    fout << m << " " << n << " " << "0.0" << "\n"; //filter one line
    for (L row = 0; row < m; ++row) {
        for (L col = A_row_ptr[row]; col < A_row_ptr[row+1]; ++col)
            fout << row+1 << " " << A_col_idx[col] + 1 << " " << A_value[col] << "\n";
    }
    fout.close();
}

template <typename L>
void writeMat(const std::string &full_data_path, const L &n){

    std::ofstream fout;
    fout.open(full_data_path.c_str());
    if( fout.fail() ){
        std::cerr << "fail to open " << full_data_path << std::endl;
        exit(0);
    }
    fout  << "0 " << n << " " << "0.0" << "\n"; //filter one line

    fout.close();
}



void scs_writeMeta(char * filename , const int &mi, const  int &me, const int &nb, const int &nf){
    FILE* fout = fopen(filename, "wb");
    if (!fout) {
        printf("Error openning file %s\n", filename);
        exit(0);
    }
    printf("Reading data from %s\n", filename);
    fwrite(&nb, sizeof(int), 1, fout);
    fwrite(&nf, sizeof(int), 1, fout);
    fwrite(&mi, sizeof(int), 1, fout);
    fwrite(&me, sizeof(int), 1, fout);
    fclose(fout);
}

void scs_writeVec(char * filename , const int &size, const std::vector<int> &vec){
    FILE* fout = fopen(filename, "wb");
    if (!fout) {
        printf("Error openning file %s\n", filename);
        exit(0);
    }
    int * vec_tmp;
    vec_tmp = (int *)calloc(size, sizeof(int));
    for (int i = 0; i < size; ++i) {
        vec_tmp[i] = vec[i];
    }
    fwrite(vec_tmp, sizeof(int), size , fout);
    fclose(fout);
}

void scs_writeVec(char * filename , const int &size, const std::vector<double> &vec){
    FILE* fout = fopen(filename, "wb");
    if (!fout) {
        printf("Error openning file %s\n", filename);
        exit(0);
    }
    double * vec_tmp;
    vec_tmp = (double *)calloc(size, sizeof(double));
    for (int i = 0; i < size; ++i) {
        vec_tmp[i] = vec[i];
    }
    fwrite(vec_tmp, sizeof(double), size , fout);
    fclose(fout);
}

int * scs_readVec_int(char * filename , const int &size){
    FILE* fin = fopen(filename, "rb");
    if (!fin) {
        printf("Error openning file %s\n", filename);
        exit(0);
    }
    int *vec = (int *)calloc(size, sizeof(int));
    rewind(fin);
    Print("finalfinal");
    fread(vec, sizeof(int),size, fin);
    for (int i = 0; i < size; ++i) {
        printf("%d",vec[i]);
    }
    printf("\n");
    printf("finalfinal");
    fclose(fin);
    return vec;
}

double * scs_readVec_double(char * filename , const int &size){
    FILE* fin = fopen(filename, "rb");
    if (!fin) {
        printf("Error openning file %s\n", filename);
        exit(0);
    }
    double *vec = (double *)calloc(size, sizeof(double));
    rewind(fin);
    Print("finalfinal");
    fread(vec, sizeof(double),size, fin);
    for (int i = 0; i < size; ++i) {
        printf("%f",vec[i]);
    }
    printf("\n");
    printf("finalfinal");
    fclose(fin);
    return vec;
}

template <typename L, typename D>
L check_class(const D &y, const L &k, const std::vector<D> &sign){
    if(sign.size()!=k){
        Print("The size of class is not accordant");
        exit(0);
    }
    L begin = 0;
    L end = k -1;
    L check_point = (begin+end)/2;
    L class_type;
    while(check_point!=begin){
        if(y == sign[check_point]) {
            class_type = check_point + 1;
            break;
        }
        else if(y < sign[check_point])
            end = check_point;
        else
            begin = check_point;
        check_point = (begin+end)/2;
    }
    if(y==sign[check_point]){
        class_type = check_point+1;
    }else if (y==sign[end]){
        class_type = end+1;
    }else{
        Print("It does not belong to any class");
        exit(0);
    }

    return class_type;
};

template <typename L, typename D>
std::vector<D> class_sign( const std::vector<D> &y){
    std::vector<D> sign;
    L iter = 0;
    sign.push_back(y[0]);
    L sign_size = 1;
    for (L i = 1; i < y.size(); ++i) {
        if(y[i] > sign[iter]){
            ++iter;
            while((y[i] > sign[iter])&&(iter<sign_size))
                ++iter;
            if(((iter<sign_size)&&(y[i]!= sign[iter]))||(iter == sign_size)) {
                sign.insert(sign.begin() + iter, y[i]);
                ++sign_size;
            }
        }else if(y[i] < sign[iter]){
            --iter;
            while((y[i] < sign[iter])&&(iter>=0))
                --iter;
            if(((iter>=0)&&(y[i]!= sign[iter]))||(iter == -1)) {
                sign.insert(sign.begin() + iter + 1, y[i]);
                ++sign_size;
            }
            if(iter == -1)
                iter = 0;
        }
    }
    return sign;
};

#endif //SDNA_LPUSE_H
