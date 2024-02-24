#include <iostream>
#include <vector>
#include <cmath>
#include <pthread.h>
#include <chrono>
#include <fstream>
#include "functions.cpp"

bool debug_mode = false ;



void luDecomposition(std::vector<std::vector<double>>& a, std::vector<int>& pi, double** l, double** u ) {

    int n = a.size();
    
    initialize_pi(pi) ;
    initialize_lu(l,u,n) ;

    for (int k = 0; k < n; ++k) {
        // Find pivot element
        double max = 0;
        int k_prime = k;
        for (int i = k; i < n; ++i) {
            if (std::abs(a[i][k]) > max) {
                max = std::abs(a[i][k]);
                k_prime = i;
            }
        }

        if (max == 0) {
            std::cerr << "Error: Singular matrix\n";
            return;
        }

        std::swap(pi[k], pi[k_prime]);

        swap(a[k], a[k_prime]);

        for (int i = 0; i < k; ++i) {
            std::swap(l[k][i], l[k_prime][i]);
        }
        u[k][k] = a[k][k];

        for (int i = k + 1; i < n; ++i) {
            l[i][k] = a[i][k] / u[k][k];
            u[k][i] = a[k][i];
        }

        for (int i = k+1; i < n; ++i) {
            for (int j = k+1; j < n; ++j) {
                a[i][j] -= l[i][k] * u[k][j]; // parallelizable
            }
        }
    }
}


int main(int argc, char* argv[]) {
    int matrix_dimension = std::stoi(argv[1]) ;
    int n = matrix_dimension;
    what_is(matrix_dimension)
    std::ifstream file("input.txt");
    std::vector<std::vector<double>> a(n, std::vector<double>(n));
    std::vector<std::vector<double>> a1(n, std::vector<double>(n));
    read_file(a,a1,n,file) ;
    std::vector<int> pi(n);
    double** l = new double*[n];
    double** u = new double*[n];
    auto start_time = std::chrono::high_resolution_clock::now();
    luDecomposition(a, pi, l, u );
    auto end_time = std::chrono::high_resolution_clock::now();
    int checking = std::stoi(argv[2]) ;
    if (checking == 1){
    std::vector<std::vector<double>> lu = matrix_mult(l,u,n) ;
    std::vector<std::vector<int>> pi_p =  expand_pi(pi) ;
    std::vector<std::vector<double>> pi_PA = matrix_mult_vec( pi_p,a1);
    std::vector<std::vector<double>> residue = subtract_matrices(pi_PA,lu) ;

    std::cout<<"Norm is "<<find_norm(residue)<<std::endl  ;
    }


    std::chrono::duration<double> time_taken = end_time - start_time ;
    std::cout<<"Time taken is "<<time_taken.count()<<std::endl ;


    return 0;
}

