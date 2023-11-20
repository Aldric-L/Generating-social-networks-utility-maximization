//
//  MatrixOperations.hpp
//  AKML Project
//
//  Created by Aldric Labarthe on 20/11/2023.
//

#ifndef MatrixOperations_hpp
#define MatrixOperations_hpp

#include "Matrix.hpp"

namespace akml {
    template <typename element_type, std::size_t ROWS, std::size_t COLUMNS>
    inline Matrix<element_type, COLUMNS, ROWS> transpose(const Matrix<element_type, ROWS, COLUMNS>& old_matrix){
        Matrix<element_type, COLUMNS, ROWS> newself;
        for (std::size_t i=1; i <= ROWS; i++){
            for (std::size_t j=1; j <= COLUMNS; j++){
                newself(j,i) = old_matrix.read(j,i);
            }
        }
        return newself;
    }

    template <typename MATRIX_INNER_TYPE, std::size_t DIM>
    inline MATRIX_INNER_TYPE inner_product(akml::Matrix<MATRIX_INNER_TYPE, DIM, 1>& a, akml::Matrix<MATRIX_INNER_TYPE, 1, DIM>& b){
        MATRIX_INNER_TYPE total(0);
        for (std::size_t i=0; i <= DIM; i++){
            total += a.read(i, 1) * b.read(1, i);
        }
        return total;
    }

    template <typename MATRIX_INNER_TYPE, std::size_t DIM>
    inline MATRIX_INNER_TYPE inner_product(akml::Matrix<MATRIX_INNER_TYPE, 1, DIM>& a, akml::Matrix<MATRIX_INNER_TYPE, DIM, 1>& b){
        return inner_product(b, a);
    }

    template <typename MATRIX_INNER_TYPE, std::size_t DIM>
    inline MATRIX_INNER_TYPE inner_product(akml::Matrix<MATRIX_INNER_TYPE, DIM, 1>& a, akml::Matrix<MATRIX_INNER_TYPE, DIM, 1>& b){
        akml::Matrix<MATRIX_INNER_TYPE, 1, DIM> b_transpose = transpose<MATRIX_INNER_TYPE, DIM, 1>(b);
        return inner_product(a, b);
    }

    template <typename element_type, std::size_t ROWS, std::size_t COLUMNS>
    inline Matrix<element_type, ROWS, COLUMNS> hadamard_product(Matrix<element_type, ROWS, COLUMNS> A, Matrix<element_type, ROWS, COLUMNS> B){
        Matrix<element_type, ROWS, COLUMNS> product;
        for (std::size_t i=1; i <= A.rows; i++){
            for (std::size_t j=1; j <= B.columns; j++){
                product(i, j) = A(i, j) * B(i, j);
            }
        }
        return product;
    };

    template <typename element_type, const std::size_t R1, const std::size_t C1, const std::size_t R2, const std::size_t C2>
    inline Matrix<element_type, R1, C2> matrix_product(Matrix<element_type, R1, C1> A, Matrix<element_type, R2, C2> B){
        if (A.columns != B.rows)
            throw std::invalid_argument("Attempting to perform a non-defined matrix product.");
       
        return Matrix<element_type, R1, C2>::product(A, B);
        
    };

    template <typename element_type, std::size_t ROWS, std::size_t COLUMNS>
    inline void cout_matrix(Matrix<element_type, ROWS, COLUMNS> A){
        return Matrix<element_type, ROWS, COLUMNS>::cout(A);
        
    };
}

#endif /* MatrixOperations_h */
