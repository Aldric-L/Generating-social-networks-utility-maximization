//
//  MatrixOperations.hpp
//  AKML Project
//
//  Created by Aldric Labarthe on 20/11/2023.
//

#include "MatrixInterface.hpp"
#include "StaticMatrix.hpp"
#include "DynamicMatrix.hpp"
#include "Matrix.hpp"

#ifndef MatrixOperations_hpp
#define MatrixOperations_hpp

namespace akml {
    template <typename MATRIX_TYPE, typename element_type>
    inline MATRIX_TYPE transform(MATRIX_TYPE matrix, std::function<element_type(element_type, std::size_t, std::size_t)> transfunc){
        matrix.transform(transfunc);
        return matrix;
    }

    template <typename MATRIX_TYPE, typename element_type>
    inline MATRIX_TYPE transform(MATRIX_TYPE matrix, std::function<element_type(element_type)> transfunc){
        matrix.transform(transfunc);
        return matrix;
    }

    
    template <typename element_type>
    inline DynamicMatrix<element_type> transpose(DynamicMatrix<element_type> old_matrix){
        if (!old_matrix.isInitialized())
            throw std::invalid_argument("Matrix provided is not initialized.");
        
        element_type* localdata;
        localdata = new element_type[(old_matrix.getNRows())*(old_matrix.getNColumns())];
        for (std::size_t i=0; i < old_matrix.getNColumns(); i++){
            for (std::size_t j=0; j < old_matrix.getNRows(); j++){
                *(localdata+i*(old_matrix.getNRows())+j) = old_matrix.read(j+1, i+1);

            }
        }
        old_matrix.deleteInternStorage();
        old_matrix.getStorage() = localdata;
        old_matrix.getStorageEnd() = localdata+(old_matrix.getNRows())*(old_matrix.getNColumns());
        std::size_t old_rows = old_matrix.getNRows();
        old_matrix.setNRows(old_matrix.getNColumns());
        old_matrix.setNColumns(old_rows);
        return old_matrix;
    }

    template <typename element_type, std::size_t ROWS, std::size_t COLUMNS>
    inline Matrix<element_type, COLUMNS, ROWS> transpose(Matrix<element_type, ROWS, COLUMNS>& old_matrix){
        if (!old_matrix.isInitialized())
            throw std::invalid_argument("Matrix provided is not initialized.");
        
        Matrix<element_type, COLUMNS, ROWS> newmatrix;
        for (std::size_t i=0; i < old_matrix.getNColumns(); i++){
            for (std::size_t j=0; j < old_matrix.getNRows(); j++){
                *(newmatrix.getStorage()+i*(old_matrix.getNRows())+j) = old_matrix.read(j+1, i+1);
            }
        }
        return newmatrix;
    }

    template <typename MATRIX_INNER_TYPE>
    inline MATRIX_INNER_TYPE inner_product(const akml::MatrixInterface<MATRIX_INNER_TYPE>& a, const akml::MatrixInterface<MATRIX_INNER_TYPE>& b){
        if (!a.isInitialized() || !b.isInitialized())
            throw std::invalid_argument("Matrix provided is not initialized.");
        
        bool A_iscol=false;
        bool B_iscol=true;
        std::size_t dim(0);
        if (a.getNRows() == 1 && b.getNColumns()==1 && a.getNColumns() == b.getNRows()){
        }else if (a.getNColumns() == 1 && b.getNRows()==1 && a.getNRows() == b.getNColumns()){
            A_iscol = true;B_iscol=false;dim=a.getNRows();
        }else if (a.getNRows() == 1 && b.getNRows()==1 && a.getNColumns() == b.getNColumns()){
            A_iscol = false;B_iscol=false;dim=a.getNColumns();
        }else if (a.getNColumns() == 1 && b.getNColumns()==1 && a.getNRows() == b.getNRows()){
            A_iscol = true;B_iscol=true;dim=a.getNRows();
        }else{
            throw std::invalid_argument("Matrices are not columns or lines");
        }
        
        MATRIX_INNER_TYPE total(0);
        for (std::size_t i=0; i <= dim; i++){
            total += (A_iscol ? a.read(i, 1) : a.read(1, i)) * (B_iscol ? b.read(i, 1) : b.read(1, i));
        }
        return total;
    }

    template <typename MATRIX_TYPE>
    inline MATRIX_TYPE hadamard_product(MATRIX_TYPE A, MATRIX_TYPE B){
        if (!A.isInitialized() || !B.isInitialized())
            throw std::invalid_argument("Matrix provided is not initialized.");
        if (A.getNColumns() != B.getNColumns() || A.getNRows() != B.getNRows())
            throw std::invalid_argument("Attempting to perform a product on non-equally sized matrix.");
        MATRIX_TYPE product;
        for (std::size_t i=1; i <= A.getNRows(); i++){
            for (std::size_t j=1; j <= B.getNColumns(); j++){
                product(i, j) = A(i, j) * B(i, j);
            }
        }
        return product;
    };

    template <typename element_type>
    inline DynamicMatrix<element_type> matrix_product(const DynamicMatrix<element_type>& A, const DynamicMatrix<element_type>& B){
        return DynamicMatrix<element_type>::product(A, B);
    };

    template <typename element_type, std::size_t ROWS, std::size_t TEMPDIM, std::size_t COLUMNS>
    inline StaticMatrix<element_type,ROWS, COLUMNS> matrix_product(const StaticMatrix<element_type, ROWS, TEMPDIM>& A, const StaticMatrix<element_type, TEMPDIM, COLUMNS>& B){
        return StaticMatrix<element_type, ROWS, COLUMNS>::product(A, B);
    };

    template <typename element_type, std::size_t ROWS, std::size_t TEMPDIM, std::size_t COLUMNS>
    inline Matrix<element_type,ROWS, COLUMNS> matrix_product(const Matrix<element_type, ROWS, TEMPDIM>& A, const Matrix<element_type, TEMPDIM, COLUMNS>& B){
        return Matrix<element_type, ROWS, COLUMNS>::product(A, B);
    };


    template <typename element_type>
    inline void cout_matrix(const MatrixInterface<element_type>& A){
        return MatrixInterface<element_type>::cout(A);
        
    };

    template <typename element_type>
    inline std::size_t arg_max(const element_type* begin, const element_type* end) {
        return static_cast<std::size_t>(std::distance(begin, std::max_element(begin, end)));
    }

    template <typename element_type>
    inline std::size_t arg_min(const element_type* begin, const element_type* end) {
        return static_cast<std::size_t>(std::distance(begin, std::min_element(begin, end)));
    }

    template <typename element_type>
    inline std::size_t arg_max(const MatrixInterface<element_type>& matrix, const bool absval_mode=false) {
        if (!matrix.isInitialized())
            throw std::invalid_argument("Matrix provided is not initialized.");
        
        if (matrix.getNRows() != 1 && matrix.getNColumns() != 1)
            throw std::invalid_argument("Arg_max only applies on column vectors, not on every matrix.");
        
        if (absval_mode){
            element_type* newstorage = new element_type[matrix.getNRows() * matrix.getNColumns()];
            for (std::size_t i(0); i < matrix.getNRows() * matrix.getNColumns(); i++){
                *(newstorage+i) = std::abs(*(matrix.getStorage()+i));
            }
            std::size_t result = arg_max(newstorage, newstorage+matrix.getNRows() * matrix.getNColumns());
            delete[] newstorage;
            return result;
        }
        return arg_max(matrix.getStorage(), matrix.getStorageEnd());
    }

    template <typename element_type>
    inline std::size_t arg_min(const MatrixInterface<element_type>& matrix, const bool absval_mode=false) {
        if (!matrix.isInitialized())
            throw std::invalid_argument("Matrix provided is not initialized.");
        
        if (matrix.getNRows() != 1 && matrix.getNColumns() != 1)
            throw std::invalid_argument("Arg_min only applies on column vectors, not on every matrix.");
        
        if (absval_mode){
            element_type* newstorage = new element_type[matrix.getNRows() * matrix.getNColumns()];
            for (std::size_t i(0); i < matrix.getNRows() * matrix.getNColumns(); i++){
                *(newstorage+i) = std::abs(*(matrix.getStorage()+i));
            }
            std::size_t result = arg_min(newstorage, newstorage+matrix.getNRows() * matrix.getNColumns());
            delete[] newstorage;
            return result;
        }
        
        return arg_min(matrix.getStorage(), matrix.getStorageEnd());
    }


    /*
     * __localComputeDijkstra - THIS METHOD SHOULD NEVER BE USED
     * It has been created to avoid code duplication due to matrix types.
     * Please use the dijkstra_distance_algorithm who does all the verifications for you
     * Freely adjusted from the article of geeksforgeeks.org
     */
    template <typename element_type>
    inline void __localComputeDijkstra(const MatrixInterface<element_type>& matrix, const std::size_t& from, MatrixInterface<element_type>* dist){
        akml::DynamicMatrix<element_type> shortestPathTreeSet(matrix.getNColumns(), 1);
     
        for (std::size_t i = 0; i < matrix.getNColumns(); i++){
            (*dist)[{i, 0}] = ULONG_MAX;
            shortestPathTreeSet[{i,0}] = false;
        }
             
        (*dist)[{from, 0}] = 0;
     
        for (std::size_t count = 0; count < matrix.getNColumns() - 1; count++) {
            std::size_t u = 0, min = ULONG_MAX;
         
            for (std::size_t v = 0; v < matrix.getNColumns(); v++){
                if (shortestPathTreeSet[{v, 0}] == false && (*dist)[{v, 0}] <= min){
                    min = (*dist)[{v, 0}];
                    u = v;
                }
            }
            shortestPathTreeSet[{u, 0}] = true;
     
            for (std::size_t v = 0; v < matrix.getNColumns(); v++)
                if (!shortestPathTreeSet[{v, 0}] && matrix[{u, v}] > 0
                    && (*dist)[{u, 0}] != INT_MAX
                    && (*dist)[{u, 0}] + matrix[{u, v}] < (*dist)[{v, 0}])
                    (*dist)[{v, 0}] = (*dist)[{u, 0}] + matrix[{u, v}];
        }
    }

    /*
     * dijkstra_distance_algorithm - Method to find the length of the shortest path between two vertices
     */
    template <typename element_type>
    inline akml::DynamicMatrix<element_type> dijkstra_distance_algorithm(const MatrixInterface<element_type>& matrix, const std::size_t from) {
        if (!matrix.isInitialized())
            throw std::invalid_argument("Matrix provided is not initialized.");
        
        if (matrix.getNRows() <= 1 || matrix.getNColumns() != matrix.getNRows() || matrix.getNColumns() < from )
            throw std::invalid_argument("Error with matrix dimension.");
        
        akml::DynamicMatrix<element_type> dist(matrix.getNColumns(), 1);
        __localComputeDijkstra(matrix, from, &dist);
        return dist;
    }

    template <typename element_type, std::size_t MATRIX_DIM>
    inline akml::Matrix<element_type, MATRIX_DIM, 1> dijkstra_distance_algorithm(const Matrix<element_type, MATRIX_DIM, 1>& matrix, const std::size_t from) {
        if (!matrix.isInitialized())
            throw std::invalid_argument("Matrix provided is not initialized.");
        
        if (matrix.getNRows() <= 1 || matrix.getNColumns() != matrix.getNRows() || matrix.getNColumns() < from )
            throw std::invalid_argument("Error with matrix dimension.");
        
        akml::Matrix<element_type, MATRIX_DIM, 1> dist;
        __localComputeDijkstra(matrix, from, &dist);
        return dist;
    }

    template <typename element_type, std::size_t MATRIX_DIM>
    inline akml::StaticMatrix<element_type, MATRIX_DIM, 1> dijkstra_distance_algorithm(const StaticMatrix<element_type, MATRIX_DIM, 1>& matrix, const std::size_t from) {
        if (!matrix.isInitialized())
            throw std::invalid_argument("Matrix provided is not initialized.");
        
        if (matrix.getNRows() <= 1 || matrix.getNColumns() != matrix.getNRows() || matrix.getNColumns() < from )
            throw std::invalid_argument("Error with matrix dimension.");
        
        akml::StaticMatrix<element_type, MATRIX_DIM, 1> dist;
        __localComputeDijkstra(matrix, from, &dist);
        return dist;
    }

}
#endif /* MatrixOperations_h */
