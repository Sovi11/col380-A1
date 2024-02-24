#include <vector>
#include <iostream>
#include <cmath>


#define debug(x) std::cout<<#x<<" is :\n" ;
#define check_here(x) std::cout<<"Here"<<x<<"\n" ;
#define what_is(x) std::cout<<#x<<" is "<<x<<"\n" ;
#define printing_matrix(x) debug(x) ;print_matrix(x) ;

template<typename T>
void swap(std::vector<T>& a, std::vector<T>& b) {
    // Ensure the vectors are of the same size
    if (a.size() != b.size())
        return; // or throw an exception

    for (size_t i = 0; i < a.size(); ++i)
        std::swap(a[i], b[i]);
}

void initialize_pi(std::vector<int> &pi){
    int n= pi.size() ;
    for (int i = 0; i < n; ++i)
        pi[i] = i;
}


std::vector<std::vector<int>> expand_pi(const std::vector<int> &pi){
    std::vector<std::vector<int>> ans ; 
    for (int i = 0 ; i < pi.size() ; i++){
        std::vector<int> v_temp ; 
        for (int j = 0 ; j < pi.size() ; j++){
            if (j == pi[i]){
                v_temp.push_back(1) ;
            }
            else{
                v_temp.push_back(0) ;
            }
        }
        ans.push_back(v_temp) ;
    }
    return ans ;

}

template<typename T1, typename T2>
std::vector<std::vector<double>> matrix_mult(T1** a,T2** b, int n){
    int m= n; 
    int k = n ; 
    std::vector<std::vector<double>> ans_matrix ; 
    for (int i = 0 ; i < n ; i++){
        std::vector<double> v_temp ; 
        for (int j=0 ; j <  k ; j++){
            double ans = 0 ; 
            for (int rr =0 ; rr < m ;  rr ++){
                ans += a[i][rr]*b[rr][j] ;
            }
            v_temp.push_back(ans) ;
        }
        ans_matrix.push_back(v_temp) ;
    }
    return ans_matrix ;
}

void initialize_lu(double** l, double** u, int n){
    for (int i = 0; i < n; ++i) {
        u[i] = new double[n];
        l[i] = new double[n];
        for (int j = 0; j < n; ++j) {
            if (i == j) {
                l[i][j] = 1; // Diagonal elements of L are 1
                u[i][j] = 0; // Diagonal elements of U are 0
            }
            else {
                l[i][j] = 0;
                u[i][j] = 0;
            }
        }
    }
}

template<typename T5, typename T6>
std::vector<std::vector<double>> matrix_mult_vec(const std::vector<std::vector<T5>> &a,const std::vector<std::vector<T6>> &b){
    int n= a.size() ; 
    int m= a[0].size() ; 
    int k = b[0].size() ; 
    std::vector<std::vector<double>> ans_matrix ; 
    for (int i = 0 ; i < n ; i++){
        std::vector<double> v_temp ; 
        for (int j=0 ; j <  k ; j++){
            double ans = 0 ; 
            for (int rr =0 ; rr < m ;  rr ++){
                ans += a[i][rr]*b[rr][j] ;
            }
            v_temp.push_back(ans) ;
        }
        ans_matrix.push_back(v_temp) ;
    }
    return ans_matrix ;
}

template<typename T3>
void print_matrix(const std::vector<std::vector<T3>> &a){
    for (int i = 0 ; i < a.size() ; i++){
        for (int j = 0 ; j < a[0].size() ; j++){
            std::cout<<a[i][j]<<" " ;
        }
        std::cout<<"\n" ;
    }
}

std::vector<std::vector<double>> subtract_matrices(const std::vector<std::vector<double>> &a , const std::vector<std::vector<double>> &b){
    int n = a.size() ; 
    int m = a[0].size() ; 
    std::vector<std::vector<double>> ans ; 
    for (int i = 0 ; i < n ; i++){
        std::vector<double> v_temp ; 
        for (int j = 0 ; j < m ; j++){
            v_temp.push_back(a[i][j]-b[i][j]) ;
        }
        ans.push_back(v_temp) ;
    }
    return ans ;
}


double find_norm(const std::vector<std::vector<double>> &a){
    // we want to find L1,2 norm here 
    double ans = 0 ; 
    for (int j = 0 ; j < a[0].size() ; j++){
        double temp_ans = 0 ;
        for (int i = 0 ; i < a.size() ; i++){
            temp_ans += a[i][j]*a[i][j] ;
        }
        ans+= sqrt(temp_ans) ;
    }

    return ans ;
}

void read_file(std::vector<std::vector<double>> &a, std::vector<std::vector<double>> &a1, int n, std::ifstream &file){
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
        {
            int temp;
            file>>temp;
            a[i][j] = (temp);
            a1[i][j] = (temp);
        }
    }
}
