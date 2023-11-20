//
//  Matrix.hpp
//  AKML Project
//
//  Created by Aldric Labarthe on 08/09/2023.
//

#ifndef Matrix_hpp
#define Matrix_hpp

namespace akml {

template <typename element_type>
class MatrixPrototype {
public:
    const std::size_t rows, columns;
    
    inline MatrixPrototype(const std::size_t rows, const std::size_t columns) : columns(columns), rows(rows) {};
    
    inline element_type read(const size_t row, const size_t column) const { return 0; }
    
    inline friend void operator<<(std::ostream& flux, MatrixPrototype<element_type> matrix){
        std::cout << "MatrixPrototype M(" << matrix.rows << ";" << matrix.columns << ")" << std::endl;
    }
    
};

template <typename element_type, std::size_t ROWS, std::size_t COLUMNS>
class Matrix : public MatrixPrototype<element_type> {
public:
    static inline std::function<element_type(element_type, std::size_t, std::size_t)> NO_ACTION_TRANSFORM = [](element_type x, std::size_t row, std::size_t column) {return x;};

    static inline std::function<element_type(element_type, std::size_t, std::size_t)> IDENTITY_TRANSFORM = [](element_type x, std::size_t row, std::size_t column) { return (column == row) ? 1 : 0;};

    static inline std::function<element_type(element_type, std::size_t, std::size_t)> RANDOM_TRANSFORM = [](element_type x, std::size_t row, std::size_t column) { std::random_device rd;  std::mt19937 gen(rd()); std::normal_distribution<double> distribution(0.0,3); return distribution(gen); };
    
    static inline Matrix<element_type, ROWS, COLUMNS> IDENTITY = Matrix<element_type, ROWS, COLUMNS>::transform(Matrix<element_type, ROWS, COLUMNS>::EMPTY, Matrix<element_type, ROWS, COLUMNS>::IDENTITY_TRANSFORM);

    static inline Matrix<element_type, ROWS, COLUMNS> EMPTY = Matrix<element_type, ROWS, COLUMNS> (true);
    
    std::array <std::array <element_type, COLUMNS>, ROWS> m_data;
    
    inline Matrix(const bool fromscratch=false) : MatrixPrototype<element_type>(ROWS, COLUMNS) {
        (!fromscratch) ? m_data = Matrix<element_type, ROWS, COLUMNS>::EMPTY.m_data : Matrix<element_type, ROWS, COLUMNS>::create(this); }
    
    inline element_type& operator()(size_t row, size_t column)
    {
        return m_data[row-1][column-1];
    }
    inline element_type read(const size_t row, const size_t column) const
    {
        return m_data[row-1][column-1];
    }
    inline element_type read(const std::array<std::size_t, 2> row_and_col)
    {
        return m_data[row_and_col[0]-1][row_and_col[1]-1];
    }
    
    inline element_type& operator[](std::array<std::size_t, 2> row_and_col)
    {
        return m_data[row_and_col[0]][row_and_col[1]];
    }
    
    inline void operator=(std::array <std::array <element_type, COLUMNS>, ROWS> data)
    {
        m_data = data;
    }
    
    inline Matrix<element_type, ROWS, COLUMNS>& operator=(Matrix<element_type, ROWS, COLUMNS> other)
    {
        m_data = other.m_data;
        return *this;
    }
    
    inline Matrix<element_type, ROWS, COLUMNS>& operator=(MatrixPrototype<element_type> other)
    {
        for (std::size_t i=0; i < other.rows; i++){
            for (std::size_t j=0; j < other.columns; j++){
                m_data[i][j] = other.read(i, j);
            }
        }
        return *this;
    }
    
    inline friend void operator<<(std::ostream& flux, Matrix<element_type, ROWS, COLUMNS> matrix){
        return Matrix<element_type, ROWS, COLUMNS>::cout(matrix);
    }
    
    inline bool operator==(const Matrix<element_type, ROWS, COLUMNS>& rhs) const {
        return m_data == rhs.m_data;
    }
    
    /*inline static Matrix<element_type, ROWS, COLUMNS>& create(){
     Matrix<element_type, ROWS, COLUMNS> matrix;
     for (std::size_t i=0; i <= ROWS; i++){
     for (std::size_t j=0; j <= COLUMNS; j++){
     matrix.m_data[0][0] = 0;
     }
     }
     return matrix;
     }*/
    
    inline static Matrix<element_type, ROWS, COLUMNS>& create(Matrix<element_type, ROWS, COLUMNS>* matrix){
        for (std::size_t i=0; i <= ROWS; i++){
            for (std::size_t j=0; j <= COLUMNS; j++){
                matrix->m_data[i][j] = 0;
            }
        }
        return *matrix;
    }
    
    inline Matrix(std::array <std::array <element_type, COLUMNS>, ROWS> data) : MatrixPrototype<element_type>(ROWS, COLUMNS) { m_data = data; }
    
    inline Matrix<element_type, ROWS, COLUMNS>& operator+=(const Matrix<element_type, ROWS, COLUMNS>& mat){
        for (std::size_t i=1; i <= ROWS; i++){
            for (std::size_t j=1; j <= COLUMNS; j++){
                m_data[i-1][j-1] = m_data[i-1][j-1] + mat.read(i, j);
            }
        }
        return *this;
    }
    
    inline friend Matrix<element_type, ROWS, COLUMNS> operator+(Matrix<element_type, ROWS, COLUMNS> selfmat, const Matrix<element_type, ROWS, COLUMNS>& mat){
        selfmat += mat;
        return selfmat;
    }
    
    inline Matrix<element_type, ROWS, COLUMNS>& operator-=(const Matrix<element_type, ROWS, COLUMNS>& mat){
        for (std::size_t i=1; i <= ROWS; i++){
            for (std::size_t j=1; j <= COLUMNS; j++){
                m_data[i-1][j-1] = m_data[i-1][j-1] - mat.read(i, j);
            }
        }
        return *this;
    }
    
    inline friend Matrix<element_type, ROWS, COLUMNS> operator-(Matrix<element_type, ROWS, COLUMNS> selfmat, const Matrix<element_type, ROWS, COLUMNS>& mat){
        selfmat -= mat;
        return selfmat;
    }
    
    template <const std::size_t R1, const std::size_t C1, const std::size_t R2, const std::size_t C2>
    //inline static Matrix<element_type, ROWS, COLUMNS> product(Matrix<element_type, ROWS, R2> A, Matrix<element_type, R2, COLUMNS> B){
    inline static Matrix<element_type, ROWS, COLUMNS> product(Matrix<element_type, R1, C1> A, Matrix<element_type, R2, C2> B){
        if (A.columns != B.rows)
            throw std::invalid_argument("Attempting to perform a non-defined matrix product.");
        
        Matrix<element_type, ROWS, COLUMNS> product;
        element_type *temp = nullptr;
        for (std::size_t i=1; i <= A.rows; i++){
            for (std::size_t j=1; j <= B.columns; j++){
                temp = new element_type;
                for (size_t k=1; k <= A.columns; k++){
                    *temp += A(i, k) * B(k, j);
                }
                product(i, j) = *temp;
                delete temp;
            }
        }
        
        return product;
        
    };
    
    inline Matrix<element_type, ROWS, COLUMNS>& transform(std::function<element_type(element_type, std::size_t, std::size_t)> transfunc){
        Matrix<element_type, ROWS, COLUMNS>::transform(*this, transfunc); return *this;
    }
    
    inline static Matrix<element_type, ROWS, COLUMNS> transform(Matrix<element_type, ROWS, COLUMNS>& matrix, std::function<element_type(element_type, std::size_t, std::size_t)> transfunc){
        for (std::size_t i=0; i < matrix.rows; i++){
            for (std::size_t j=0; j < matrix.columns; j++){
                matrix.m_data[i][j] = transfunc(matrix.m_data[i][j], i, j);
            }
        }
        return matrix;
    }
    
    inline Matrix<element_type, ROWS, COLUMNS>& transform(std::function<element_type(element_type)> transfunc){
        Matrix<element_type, ROWS, COLUMNS>::transform(*this, transfunc); return *this;
    }
    
    inline static Matrix<element_type, ROWS, COLUMNS> transform(Matrix<element_type, ROWS, COLUMNS>& matrix, std::function<element_type(element_type)> transfunc){
        for (std::size_t i=0; i < matrix.rows; i++){
            for (std::size_t j=0; j < matrix.columns; j++){
                matrix.m_data[i][j] = transfunc(matrix.m_data[i][j]);
            }
        }
        return matrix;
    }
    
    inline void transpose(){
        if (COLUMNS != ROWS)
            throw std::invalid_argument("Attempting to perform a direct transposition on a non-squared matrix.");
        
        std::array <std::array <element_type, COLUMNS>, ROWS> new_m_data;
        for (std::size_t i=0; i <= ROWS; i++){
            for (std::size_t j=0; j <= COLUMNS; j++){
                new_m_data[j][i] = m_data[i][j];
            }
        }
        m_data = new_m_data;
    }
    
    inline static void cout(Matrix<element_type, ROWS, COLUMNS> matrix){
        //std::cout << "Printing Matrix <" << matrix.rows << " " << matrix.columns << ">" << std::endl;
        std::cout << "[";
        for (int i(0); i < matrix.rows; i++){
            if (i != 0)
                std::cout << " ";
            std::cout << "[ ";
            for (int j(0); j<matrix.columns; j++){
                std::cout << matrix(i+1,j+1) << " ";
                if (j != matrix.columns-1)
                    std::cout << "; ";
            }
            if (i != matrix.rows-1)
                std::cout << "]," << std::endl;
            else
                std::cout << "]";
        }
        std::cout << "]" << std::endl;
    }
};

}

#endif /* Matrix_hpp */
