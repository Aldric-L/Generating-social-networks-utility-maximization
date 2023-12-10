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
    
    inline element_type read(const size_t row, const size_t column) const { return static_cast<element_type>(0); }
    
    inline friend void operator<<(std::ostream& flux, MatrixPrototype<element_type> matrix){
        std::cout << "MatrixPrototype M(" << matrix.rows << ";" << matrix.columns << ")" << std::endl;
    }
    
};

template <typename element_type, std::size_t ROWS, std::size_t COLUMNS>
class Matrix : public MatrixPrototype<element_type> {
protected:
    element_type* m_data = nullptr;
    
    inline void create(){
        m_data = new element_type[ROWS*COLUMNS];
        for (std::size_t i(0); i < ROWS*COLUMNS; i++){
            *(m_data + i) = static_cast<element_type>(0);
        }
     }
    
public:
    static inline std::function<element_type(element_type, std::size_t, std::size_t)> NO_ACTION_TRANSFORM = [](element_type x, std::size_t row, std::size_t column) {return x;};

    static inline std::function<element_type(element_type, std::size_t, std::size_t)> IDENTITY_TRANSFORM = [](element_type x, std::size_t row, std::size_t column) { return (column == row) ? static_cast<element_type>(1) : static_cast<element_type>(0);};

    static inline std::function<element_type(element_type, std::size_t, std::size_t)> RANDOM_TRANSFORM = [](element_type x, std::size_t row, std::size_t column) { std::random_device rd;  std::mt19937 gen(rd()); std::normal_distribution<double> distribution(0.0,3); return distribution(gen); };
    
    /*static std::function<element_type(element_type, std::size_t, std::size_t)> NO_ACTION_TRANSFORM;

    static std::function<element_type(element_type, std::size_t, std::size_t)> IDENTITY_TRANSFORM;

    static std::function<element_type(element_type, std::size_t, std::size_t)> RANDOM_TRANSFORM;*/
    
    //static Matrix<element_type, ROWS, COLUMNS> EMPTY;
    
    //static Matrix<element_type, ROWS, COLUMNS> IDENTITY;
    
    //static inline akml::Matrix<element_type, ROWS, COLUMNS> EMPTY = akml::Matrix<element_type, ROWS, COLUMNS> (true);

    //static inline akml::Matrix<element_type, ROWS, COLUMNS> IDENTITY = akml::Matrix<element_type, ROWS, COLUMNS> ([](element_type x, std::size_t row, std::size_t column) { return (column == row) ? static_cast<element_type>(1) : static_cast<element_type>(0);});
    
    inline Matrix(const bool fromscratch=false) : MatrixPrototype<element_type>(ROWS, COLUMNS) {
        Matrix<element_type, ROWS, COLUMNS>::create();
        //std::cout << "A matrix is born " << this << " - Mdata ref " << m_data << std::endl;
        return;
        /*if(fromscratch || !Matrix<element_type, ROWS, COLUMNS>::EMPTY.isInitialized()){
            Matrix<element_type, ROWS, COLUMNS>::create();
        }else {
            m_data = new element_type[ROWS*COLUMNS];
            std::copy(Matrix<element_type, ROWS, COLUMNS>::EMPTY.m_data, Matrix<element_type, ROWS, COLUMNS>::EMPTY.m_data + ROWS*COLUMNS, m_data);
        }*/
    }

    //Column-based constructor
    inline Matrix(const std::array<akml::Matrix<element_type, ROWS, 1>, COLUMNS>& cols) : MatrixPrototype<element_type>(ROWS, COLUMNS) {
        m_data = new element_type[ROWS*COLUMNS];
        for (std::size_t line(0); line < ROWS; line++){
            for (std::size_t col(0); col < COLUMNS; col++){
                *(m_data + line*COLUMNS+col) = cols[col].read(line+1, 1);
            }
        }
        //std::cout << "A matrix is born " << this << " - Mdata ref " << m_data << std::endl;
    }
    
    inline Matrix(const std::array <std::array <element_type, COLUMNS>, ROWS>& data) : MatrixPrototype<element_type>(ROWS, COLUMNS) {
        m_data = new element_type[ROWS*COLUMNS];
        for (std::size_t line(0); line < ROWS; line++){
            for (std::size_t col(0); col < COLUMNS; col++){
                *(m_data + line*COLUMNS+col) = data.read(line+1, col+1);
            }
        }
        //std::cout << "A matrix is born " << this << " - Mdata ref " << m_data << std::endl;
    }
    
    //copy constructor
    inline Matrix(const Matrix<element_type, ROWS, COLUMNS>& other) : MatrixPrototype<element_type>(ROWS, COLUMNS) {
        m_data = new element_type[ROWS*COLUMNS];
        std::copy(other.getStorage(), other.getStorage() + ROWS * COLUMNS, m_data);
        //std::cout << "A matrix is born " << this << " - Mdata ref " << m_data << std::endl;
    }
    
    inline Matrix(std::function<element_type(element_type, std::size_t, std::size_t)>& transfunc) : MatrixPrototype<element_type>(ROWS, COLUMNS) {
        Matrix<element_type, ROWS, COLUMNS>::create();
        /*if(!Matrix<element_type, ROWS, COLUMNS>::EMPTY.isInitialized()){
            Matrix<element_type, ROWS, COLUMNS>::create();
        }else {
            m_data = new element_type[ROWS*COLUMNS];
            std::copy(Matrix<element_type, ROWS, COLUMNS>::EMPTY.m_data, Matrix<element_type, ROWS, COLUMNS>::EMPTY.m_data + ROWS*COLUMNS, m_data);
        }*/
        transform(transfunc);
        //std::cout << "A matrix is born " << this << " - Mdata ref " << m_data << std::endl;
    }
    
    inline Matrix(std::function<element_type(element_type)>& transfunc) : MatrixPrototype<element_type>(ROWS, COLUMNS) {
        Matrix<element_type, ROWS, COLUMNS>::create();
        /*if(!Matrix<element_type, ROWS, COLUMNS>::EMPTY.isInitialized()){
            Matrix<element_type, ROWS, COLUMNS>::create();
        }else {
            m_data = new element_type[ROWS*COLUMNS];
            std::copy(Matrix<element_type, ROWS, COLUMNS>::EMPTY.m_data, Matrix<element_type, ROWS, COLUMNS>::EMPTY.m_data + ROWS*COLUMNS, m_data);
        }*/
        transform(transfunc);
        //std::cout << "A matrix is born " << this << " - Mdata ref " << m_data << std::endl;
    }
    
    
    inline ~Matrix(){
        //std::cout << "Deleter called on " << this << " - Mdata ref " << m_data << std::endl;
        if (m_data != nullptr)
            delete[] m_data;
        m_data = nullptr;
    }
    
    inline Matrix<element_type, ROWS, COLUMNS>& operator=(Matrix<element_type, ROWS, COLUMNS>&& other){
        if (this != &other) {
            delete[] m_data;  // Delete old data before copying
            m_data = new element_type[ROWS * COLUMNS];
            std::copy(other.getStorage(), other.getStorage() + ROWS * COLUMNS, m_data);
        }
        return *this;
    }
    
    inline Matrix(Matrix<element_type, ROWS, COLUMNS>&& other) : MatrixPrototype<element_type>(ROWS, COLUMNS) {
        m_data = new element_type[ROWS*COLUMNS];
        m_data = std::move(other.getStorage());
    }
    
    inline bool isInitialized() const { return !(m_data==nullptr); }
    inline element_type* getStorage() { return m_data; }
    inline element_type* getStorage() const { return m_data; }

    
    inline element_type& operator[](const std::array<std::size_t, 2> row_and_col) {
        return *(m_data + row_and_col[0]*COLUMNS+row_and_col[1]);
    }
    
    inline element_type& operator[](const std::array<std::size_t, 2> row_and_col) const {
        return *(m_data + row_and_col[0]*COLUMNS+row_and_col[1]);
    }
    
    inline element_type read(const size_t row, const size_t column) const {
        return *(m_data + (row-1)*COLUMNS+(column-1));
    }
    
    inline element_type& operator()(size_t row, size_t column) { return operator[]({row-1, column-1}); }
    inline element_type read(const std::array<std::size_t, 2> row_and_col) { return operator[]({row_and_col[0]-1, row_and_col[1]-1}); }
    
    
    inline void operator=(const std::array <std::array <element_type, COLUMNS>, ROWS>& data)
    {
        if (!isInitialized())
            m_data = new element_type[ROWS*COLUMNS];
        
        for (std::size_t line(0); line < ROWS; line++){
            for (std::size_t col(0); col < COLUMNS; col++){
                *(m_data + line*COLUMNS+col) = data[line][col];
            }
        }
    }
    
    inline Matrix<element_type, ROWS, COLUMNS>& operator=(const Matrix<element_type, ROWS, COLUMNS>& other)
    {
        if (this != &other){
            if (!other.isInitialized())
                throw std::invalid_argument("Matrix provided is not initialized.");
            
            if (!isInitialized())
                m_data = new element_type[ROWS*COLUMNS];
            
            //std::memcpy(m_data, other.m_data, ROWS * COLUMNS * sizeof(element_type));
            std::copy(other.m_data, other.m_data + ROWS * COLUMNS, m_data);
            
            /*for (std::size_t line(0); line < ROWS; line++){
                for (std::size_t col(0); col < COLUMNS; col++){
                    *(m_data + line*COLUMNS+col) = other.read(line+1, col+1);
                }
            }*/
        }
        
        return *this;
    }
    
    /*inline Matrix<element_type, ROWS, COLUMNS>& operator=(const MatrixPrototype<element_type>& other)
    {
        if (!other.isInitialized())
            throw std::invalid_argument("Matrix provided is not initialized.");
        
        for (std::size_t line(0); line < ROWS; line++){
            for (std::size_t col(0); col < COLUMNS; col++){
                *(m_data + line*COLUMNS+col) = other.read(line+1, col+1);
            }
        }
        return *this;
    }*/
    
    Matrix<element_type, ROWS, COLUMNS>& operator=(const MatrixPrototype<element_type>& other);
    
    inline friend void operator<<(std::ostream& flux, Matrix<element_type, ROWS, COLUMNS> matrix){
        if (!matrix.isInitialized())
            throw std::invalid_argument("Matrix provided is not initialized.");
        
        return Matrix<element_type, ROWS, COLUMNS>::cout(matrix);
    }
    
    inline bool operator==(const Matrix<element_type, ROWS, COLUMNS>& rhs) const {
        if (!rhs.isInitialized())
            throw std::invalid_argument("Matrix provided is not initialized.");
        
        return m_data == rhs.m_data;
    }
    
    /* DEPRECATED
    inline static Matrix<element_type, ROWS, COLUMNS>& create(Matrix<element_type, ROWS, COLUMNS>* matrix){
        for (std::size_t i=0; i <= ROWS; i++){
            for (std::size_t j=0; j <= COLUMNS; j++){
                matrix->m_data[i][j] = 0;
            }
        }
        return *matrix;
    }*/
    
    inline Matrix<element_type, ROWS, COLUMNS>& operator+=(const Matrix<element_type, ROWS, COLUMNS>& mat){
        if (!mat.isInitialized())
            throw std::invalid_argument("Matrix provided is not initialized.");
        
        for (std::size_t line(0); line < ROWS; line++){
            for (std::size_t col(0); col < COLUMNS; col++){
                *(m_data + line*COLUMNS+col) += mat.read(line+1, col+1);
            }
        }
        return *this;
    }
    
    inline friend Matrix<element_type, ROWS, COLUMNS> operator+(Matrix<element_type, ROWS, COLUMNS>& selfmat, const Matrix<element_type, ROWS, COLUMNS>& mat){
        if (!mat.isInitialized() || !selfmat.isInitialized())
            throw std::invalid_argument("Matrix provided is not initialized.");
        selfmat += mat;
        return selfmat;
    }
    
    inline friend Matrix<element_type, ROWS, COLUMNS> operator*(Matrix<element_type, ROWS, COLUMNS>& selfmat, const Matrix<element_type, ROWS, COLUMNS>& mat){
        if (selfmat.columns != mat.rows)
            throw std::invalid_argument("Attempting to perform a non-defined matrix product.");
        if (!mat.isInitialized() || !selfmat.isInitialized())
            throw std::invalid_argument("Matrix provided is not initialized.");
        selfmat *= mat;
        return selfmat;
    }
    
    inline Matrix<element_type, ROWS, COLUMNS>& operator*=(const Matrix<element_type, ROWS, COLUMNS>& mat){
        if (this->columns != mat.rows)
            throw std::invalid_argument("Attempting to perform a non-defined matrix product.");
        if (!mat.isInitialized())
            throw std::invalid_argument("Matrix provided is not initialized.");
        
        *this = Matrix<element_type, ROWS, COLUMNS>::product(*this, mat);
        return *this;
    }
    
    inline Matrix<element_type, ROWS, COLUMNS>& operator-=(const Matrix<element_type, ROWS, COLUMNS>& mat){
        if (!mat.isInitialized())
            throw std::invalid_argument("Matrix provided is not initialized.");
        
        for (std::size_t line(0); line < ROWS; line++){
            for (std::size_t col(0); col < COLUMNS; col++){
                *(m_data + line*COLUMNS+col) -= mat.read(line+1, col+1);
            }
        }
        return *this;
    }
    
    inline friend Matrix<element_type, ROWS, COLUMNS> operator-(Matrix<element_type, ROWS, COLUMNS>& selfmat, const Matrix<element_type, ROWS, COLUMNS>& mat){
        if (!mat.isInitialized() || !selfmat.isInitialized())
            throw std::invalid_argument("Matrix provided is not initialized.");
        
        selfmat -= mat;
        return selfmat;
    }
    
    template <const std::size_t R1, const std::size_t C1, const std::size_t R2, const std::size_t C2>
    //inline static Matrix<element_type, ROWS, COLUMNS> product(Matrix<element_type, ROWS, R2> A, Matrix<element_type, R2, COLUMNS> B){
    inline static Matrix<element_type, ROWS, COLUMNS> product(const Matrix<element_type, R1, C1>& A, const Matrix<element_type, R2, C2>& B){
        if (A.columns != B.rows)
            throw std::invalid_argument("Attempting to perform a non-defined matrix product.");
        if (!A.isInitialized() || !B.isInitialized())
            throw std::invalid_argument("Matrix provided is not initialized.");
        
        Matrix<element_type, ROWS, COLUMNS> product;
        element_type *temp = nullptr;
        for (std::size_t i=1; i <= A.rows; i++){
            for (std::size_t j=1; j <= B.columns; j++){
                temp = new element_type;
                for (size_t k=1; k <= A.columns; k++){
                    *temp += A.read(i, k) * B.read(k, j);
                }
                product(i, j) = *temp;
                delete temp;
            }
        }
        
        return product;
        
    };
    
    inline Matrix<element_type, ROWS, COLUMNS>& transform(std::function<element_type(element_type, std::size_t, std::size_t)> transfunc){
        for (std::size_t i(0); i < ROWS; i++){
            for (std::size_t j(0); j < COLUMNS; j++){
                operator[]({i, j}) = transfunc(operator[]({i, j}), i, j);
            }
        }
        return *this;
    }
    
    inline static Matrix<element_type, ROWS, COLUMNS> transform(const Matrix<element_type, ROWS, COLUMNS>& matrix, std::function<element_type(element_type, std::size_t, std::size_t)> transfunc){
        if (!matrix.isInitialized())
            throw std::invalid_argument("Matrix provided is not initialized.");
        
        Matrix<element_type, ROWS, COLUMNS> new_matrix;
        for (std::size_t i(0); i < matrix.rows; i++){
            for (std::size_t j(0); j < matrix.columns; j++){
                new_matrix[{i, j}] = transfunc(matrix.operator[]({i, j}), i, j);
            }
        }
        return new_matrix;
    }
    
    inline Matrix<element_type, ROWS, COLUMNS>& transform(std::function<element_type(element_type)> transfunc){
        for (std::size_t i(0); i < ROWS; i++){
            for (std::size_t j(0); j < COLUMNS; j++){
                operator[]({i, j}) = transfunc(operator[]({i, j}));
            }
        }
        return *this;
    }
    
    inline static Matrix<element_type, ROWS, COLUMNS> transform(Matrix<element_type, ROWS, COLUMNS>& matrix, std::function<element_type(element_type)> transfunc){
        if (!matrix.isInitialized())
            throw std::invalid_argument("Matrix provided is not initialized.");
        Matrix<element_type, ROWS, COLUMNS> new_matrix;
        for (std::size_t i(0); i < matrix.rows; i++){
            for (std::size_t j(0); j < matrix.columns; j++){
                new_matrix[{i, j}] = transfunc(matrix.operator[]({i, j}));
            }
        }
        return new_matrix;
    }
    
    inline Matrix<element_type, ROWS, COLUMNS> transpose(){
        if (COLUMNS != ROWS)
            throw std::invalid_argument("Attempting to perform a direct transposition on a non-squared matrix.");
        
        Matrix<element_type, ROWS, COLUMNS> new_matrix;
        for (std::size_t i(0); i < ROWS; i++){
            for (std::size_t j(0); j < COLUMNS; j++){
                new_matrix[{i, j}] = read(j+1, i+1);
            }
        }
        for (std::size_t i(0); i < ROWS; i++){
            for (std::size_t j(0); j < COLUMNS; j++){
                operator[]({i, j}) = new_matrix.read(i+1, j+1);
            }
        }
        return new_matrix;
    }
    
    inline static void cout(const Matrix<element_type, ROWS, COLUMNS>& matrix){
        //std::cout << "Printing Matrix <" << matrix.rows << " " << matrix.columns << ">" << std::endl;
        std::cout << "[";
        for (std::size_t i(0); i < matrix.rows; i++){
            if (i != 0)
                std::cout << " ";
            std::cout << "[ ";
            for (std::size_t j(0); j<matrix.columns; j++){
                std::cout << matrix.read(i+1,j+1) << " ";
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
