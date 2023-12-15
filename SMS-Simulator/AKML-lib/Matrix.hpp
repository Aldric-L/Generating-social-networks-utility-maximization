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
class MatrixInterface {
protected:
    element_type* m_data = nullptr;
    element_type* m_data_end = nullptr;
    std::size_t rows, columns;
    
public:
    static inline std::function<element_type(element_type, std::size_t, std::size_t)> NO_ACTION_TRANSFORM = [](element_type x, std::size_t row, std::size_t column) {return x;};
    
    static inline std::function<element_type(element_type, std::size_t, std::size_t)> IDENTITY_TRANSFORM = [](element_type x, std::size_t row, std::size_t column) { return (column == row) ? static_cast<element_type>(1) : static_cast<element_type>(0);};
    
    static inline std::function<element_type(element_type, std::size_t, std::size_t)> RANDOM_TRANSFORM = [](element_type x, std::size_t row, std::size_t column) { std::random_device rd;  std::mt19937 gen(rd()); std::normal_distribution<double> distribution(0.0,3); return distribution(gen); };
    
    MatrixInterface(const std::size_t rows, const std::size_t columns) : rows(rows), columns(columns) {};
    
    virtual void create() = 0;
    virtual void createInternStorage() = 0;
    virtual void deleteInternStorage() = 0;
    
    inline std::size_t getNColumns() const { return this->columns; }
    inline std::size_t getNRows() const { return this->rows; }
    
    inline bool isInitialized() const { return !(this->m_data==nullptr); }
    inline element_type*& getStorage() { return this->m_data; }
    inline element_type*& getStorageEnd() { return this->m_data_end; }
    inline element_type* getStorage() const { return this->m_data; }
    inline element_type* getStorageEnd() const { return this->m_data_end; }
    inline element_type* getInternElement(const std::size_t pos) const {
        if (pos >= (this->rows) * (this->columns))
            throw std::invalid_argument("Attempt to access to an out of reach element of a matrix");
        
        return (this->m_data + pos);
    }
    
    inline void setMDataPointer(element_type* array_start){
        this->m_data = array_start;
        this->m_data_end = this->m_data + (this->rows)*(this->columns);
    }
    
    inline element_type& operator[](const std::array<std::size_t, 2> row_and_col) const {
        if (this->rows <= row_and_col[0] || this->columns <= row_and_col[1])
            throw std::invalid_argument("Attempt to access to an out of reach element of a matrix");
        return *(this->m_data + row_and_col[0]*(this->columns)+row_and_col[1]);
    }
    
    inline element_type read(const size_t row, const size_t column) const {
        if (this->rows < row || this->columns < column)
            throw std::invalid_argument("Attempt to access to an out of reach element of a matrix");
        return *(this->m_data + (row-1)*(this->columns)+(column-1));
    }
    
    inline element_type& operator[](const std::array<std::size_t, 2> row_and_col) {
        return *(this->getInternElement(row_and_col[0]*(this->columns)+row_and_col[1]));
    }
    
    inline element_type read(const std::array<std::size_t, 2> row_and_col) const { return operator[]({row_and_col[0]-1, row_and_col[1]-1}); }
    
    inline element_type& operator()(const size_t row, const size_t column) { return operator[]({row-1, column-1}); }
    
    inline void transform(std::function<element_type(element_type, std::size_t, std::size_t)> transfunc){
        for (std::size_t i(0); i < (this->rows); i++){
            for (std::size_t j(0); j < (this->columns); j++){
                operator[]({i, j}) = transfunc(operator[]({i, j}), i, j);
            }
        }
    }
    
    inline void transform(std::function<element_type(element_type)> transfunc){
        for (std::size_t i(0); i < (this->rows)*(this->columns); i++){
            *(getInternElement(i)) = transfunc(*getInternElement(i));
        }
    }
    
    inline static void cout(const MatrixInterface<element_type>& matrix){
        //std::cout << "Printing Matrix <" << matrix.getNRows() << " " << matrix.getNColumns() << ">" << std::endl;
        std::cout << "[";
        for (std::size_t i(0); i < matrix.getNRows(); i++){
            if (i != 0)
                std::cout << " ";
            std::cout << "[ ";
            for (std::size_t j(0); j<matrix.getNColumns(); j++){
                std::cout << matrix.read(i+1,j+1) << " ";
                if (j != matrix.getNColumns()-1)
                    std::cout << "; ";
            }
            if (i != matrix.getNRows()-1)
                std::cout << "]," << std::endl;
            else
                std::cout << "]";
        }
        std::cout << "]" << std::endl;
    }
    
    inline bool operator==(const MatrixInterface<element_type>& rhs) const {
        if (!rhs.isInitialized())
            throw std::invalid_argument("Matrix provided is not initialized.");
        
        element_type* tar_other(rhs.getStorage());
        for (element_type* tar(getStorage()); tar < getStorageEnd(); tar++){
            if (*tar != *tar_other)
                return false;
            tar_other++;
        }
        
        return true;
    }
};



template <typename element_type>
class DynamicMatrix : public MatrixInterface<element_type>{
public:
    inline void setNColumns(std::size_t& c) { resize(this->rows, c); }
    inline void setNRows(std::size_t& r) { resize(r, this->columns); }
    
    inline void resize(std::size_t r, std::size_t c){
        if (this->isInitialized()){
            element_type* newstorage = new element_type[r*c];
            for (std::size_t line(0); line < r; line++){
                for (std::size_t col(0); col < c; col++){
                    if (col >= this->columns || line >= this->rows)
                        *(newstorage + line*(c)+col) = static_cast<element_type>(0);
                    else
                        *(newstorage + line*(c)+col) = this->operator[]({line, col});
                }
            }
            delete[] this->m_data;
            this->setMDataPointer(newstorage);
            //this->m_data = newstorage;
        }
        this->columns = c;
        this->rows = r;
    }
    
    inline void create(){
        this->createInternStorage();
        for (std::size_t i(0); i < (this->rows)*(this->columns); i++){
            *(this->m_data + i) = static_cast<element_type>(0);
        }
    }
    
    inline void createInternStorage(){
        if (this->m_data != nullptr || this->m_data_end != nullptr)
            this->deleteInternStorage();
        
        this->m_data = new element_type[(this->rows)*(this->columns)];
        this->m_data_end = this->m_data+(this->rows)*(this->columns);
    }
    
    inline void deleteInternStorage(){
        if (this->m_data != nullptr)
            delete[] this->m_data;
        this->m_data = nullptr;
        this->m_data_end = nullptr;
    }
    
    inline DynamicMatrix(const std::size_t rows, const std::size_t columns, const bool fromscratch=false) : MatrixInterface<element_type>(rows, columns) {
        (fromscratch) ? this->createInternStorage() : this->create();
        //std::cout << "A matrix is born " << this << " - Mdata ref " << m_data << std::endl;
    }
    
    //Column-based constructor
    template<std::size_t COLUMNS>
    inline DynamicMatrix(const std::array<akml::DynamicMatrix<element_type>, COLUMNS>& cols) : MatrixInterface<element_type>(cols[0].getNRows(), COLUMNS) {
        if (cols.size() == 0)
            throw std::invalid_argument("Empty initialize list.");
        
        this->createInternStorage();
        for (std::size_t line(0); line < this->rows; line++){
            for (std::size_t col(0); col < this->columns; col++){
                *(this->m_data + line*(this->columns)+col) = cols[col].read(line+1, 1);
            }
        }
        
        //std::cout << "A matrix is born " << this << " - Mdata ref " << m_data << std::endl;
    }
    
    template<std::size_t ROWS, std::size_t COLUMNS>
    inline DynamicMatrix(const std::array <std::array <element_type, COLUMNS>, ROWS>& data) : MatrixInterface<element_type>(ROWS, COLUMNS) {
        this->createInternStorage();
        for (std::size_t line(0); line < ROWS; line++){
            for (std::size_t col(0); col < COLUMNS; col++){
                *(this->m_data + line*COLUMNS+col) = data[line][col];
            }
        }
        //std::cout << "A matrix is born " << this << " - Mdata ref " << m_data << std::endl;
    }
    
    //copy constructor
    inline DynamicMatrix(const DynamicMatrix<element_type>& other) : MatrixInterface<element_type>(other.getNRows(), other.getNColumns()) {
        //deleteInternStorage();
        this->createInternStorage();
        std::copy(other.getStorage(), other.getStorageEnd(), this->m_data);
        //std::cout << "A matrix is born " << this << " - Mdata ref " << m_data << std::endl;
    }
    
    inline DynamicMatrix(const std::size_t rows, const std::size_t columns, std::function<element_type(element_type, std::size_t, std::size_t)>& transfunc) : MatrixInterface<element_type>(rows, columns) {
        this->create();
        this->transform(transfunc);
        //std::cout << "A matrix is born " << this << " - Mdata ref " << m_data << std::endl;
    }
    
    inline DynamicMatrix(const std::size_t rows, const std::size_t columns, std::function<element_type(element_type)>& transfunc) : MatrixInterface<element_type>(rows, columns) {
        this->create();
        this->transform(transfunc);
        //std::cout << "A matrix is born " << this << " - Mdata ref " << m_data << std::endl;
    }
    
    
    inline ~DynamicMatrix(){
        //std::cout << "Deleter called on " << this << " - Mdata ref " << m_data << std::endl;
        this->deleteInternStorage();
    }
    
    inline DynamicMatrix<element_type>& operator=(DynamicMatrix<element_type>&& other){
        if (this != &other) {
            if (other.getNColumns() != this->columns || other.getNRows() != this->rows)
                throw std::invalid_argument("Matrix should be equally sized to be assignable.");
            //this->deleteInternStorage(); this->createInternStorage();
            //std::copy(other.getStorage(), other.getStorageEnd(), this->m_data);
            this->setMDataPointer(other.getStorage());
            //this->m_data = other.getStorage();
            other.getStorage() = nullptr;
        }
        return *this;
    }
    
    inline DynamicMatrix(DynamicMatrix<element_type>&& other) : MatrixInterface<element_type>(other.getNRows(), other.getNColumns()) {
        //createInternStorage();
        //m_data = std::move(other.getStorage());
        //deleteInternStorage();
        //this->m_data = other.getStorage();
        this->setMDataPointer(other.getStorage());
        other.getStorage() = nullptr;
        // IS IT RELEVANT ???
        //delete[] other.getStorage();
        //other.getStorage() = nullptr;
        
    }
    
    template<std::size_t ROWS, std::size_t COLUMNS>
    inline void operator=(const std::array <std::array <element_type, COLUMNS>, ROWS>& data){
        if (COLUMNS != (this->columns) || ROWS != (this->rows))
            throw std::invalid_argument("Array and Matrix should be equally sized to be assignable.");
        
        if (!this->isInitialized())
            this->createInternStorage();
        
        for (std::size_t line(0); line < ROWS; line++){
            for (std::size_t col(0); col < COLUMNS; col++){
                *(this->getInternElement(line*COLUMNS+col)) = data[line][col];
            }
        }
    }
    
    inline DynamicMatrix<element_type>& operator=(const DynamicMatrix<element_type>& other){
        if (this != &other){
            if (other.getNColumns() != (this->columns) || other.getNRows() != (this->rows))
                throw std::invalid_argument("Matrices should be equally sized to be assignable.");
            
            if (!other.isInitialized())
                throw std::invalid_argument("DynamicMatrix provided is not initialized.");
            
            if (!this->isInitialized())
                this->createInternStorage();
            
            std::copy(other.getStorage(), other.getStorageEnd(), this->m_data);
        }
        
        return *this;
    }
    
    inline friend void operator<<(std::ostream& flux, DynamicMatrix<element_type> matrix){
        if (!matrix.isInitialized())
            throw std::invalid_argument("Matrix provided is not initialized.");
        
        return DynamicMatrix<element_type>::cout(matrix);
    }
    
    inline DynamicMatrix<element_type>& operator+=(const DynamicMatrix<element_type>& mat){
        if (mat.getNColumns() != (this->columns) || mat.getNRows() != (this->rows))
            throw std::invalid_argument("Matrices should be equally sized to be assignable.");
        
        if (!mat.isInitialized())
            throw std::invalid_argument("Matrix provided is not initialized.");
        
        for (std::size_t i(0); i < (this->rows)*(this->columns); i++){
            *(this->getInternElement(i)) += *(mat.getInternElement(i));
        }
        return *this;
    }
    
    inline friend DynamicMatrix<element_type> operator+(DynamicMatrix<element_type> selfmat, const DynamicMatrix<element_type>& mat){
        selfmat += mat;
        return selfmat;
    }
    
    inline friend DynamicMatrix<element_type> operator*(DynamicMatrix<element_type> selfmat, const DynamicMatrix<element_type>& mat){
        selfmat *= mat;
        return selfmat;
    }
    
    template <typename NumericType>
    inline friend DynamicMatrix<element_type> operator*(DynamicMatrix<element_type> selfmat, const NumericType& num){
        static_assert(std::is_arithmetic<NumericType>::value, "NumericType must be numeric");
        for (element_type* tar(selfmat.getStorage()); tar < selfmat.getStorageEnd(); tar++){
            *tar = *tar * num;
        }
        return selfmat;
    }
    
    template <typename NumericType>
    inline friend DynamicMatrix<element_type> operator*(const NumericType& num, DynamicMatrix<element_type> selfmat){
        return operator*(selfmat, num);
    }
    
    inline DynamicMatrix<element_type>& operator*=(const DynamicMatrix<element_type>& mat){
        if (!mat.isInitialized())
            throw std::invalid_argument("Matrix provided is not initialized.");
        
        *this = product(*this, mat);
        return *this;
    }
    
    inline DynamicMatrix<element_type>& operator-=(const DynamicMatrix<element_type>& mat){
        if (!mat.isInitialized())
            throw std::invalid_argument("Matrix provided is not initialized.");
        
        if (mat.getNColumns() != (this->columns) || mat.getNRows() != (this->rows))
            throw std::invalid_argument("Matrices should be equally sized to be assignable.");
        
        for (std::size_t i(0); i < (this->rows)*(this->columns); i++){
            *(this->getInternElement(i)) -= *(mat.getInternElement(i));
        }

        return *this;
    }
    
    inline friend DynamicMatrix<element_type> operator-(DynamicMatrix<element_type> selfmat, const DynamicMatrix<element_type>& mat){
        selfmat -= mat;
        return selfmat;
    }
    
    
    inline void transpose(){
        if ((this->columns) != (this->rows))
            throw std::invalid_argument("Attempting to perform a direct transposition on a non-squared matrix.");
        
        element_type* localdata;
        localdata = new element_type[(this->rows)*(this->columns)];
        
        for (std::size_t i(0); i < (this->rows); i++){
            for (std::size_t j(0); j < (this->columns); j++){
                *(localdata+i*(this->rows)+j) = this->read(j+1, i+1);
            }
        }
        delete[] this->m_data;
        this->setMDataPointer(localdata);
        //this->m_data = localdata;
    }
    
    inline static DynamicMatrix<element_type> product(const DynamicMatrix<element_type>& A, const DynamicMatrix<element_type>& B){
        if (A.getNColumns() != B.getNRows())
            throw std::invalid_argument("Attempting to perform a non-defined matrix product.");
        if (!A.isInitialized() || !B.isInitialized())
            throw std::invalid_argument("Matrix provided is not initialized.");
        
        DynamicMatrix<element_type> product(A.getNRows(), B.getNColumns());
        for (std::size_t i=1; i <= A.getNRows(); i++){
            for (std::size_t j=1; j <= B.getNColumns(); j++){
                for (size_t k=1; k <= A.getNColumns(); k++){
                    product(i, j) += A.read(i, k) * B.read(k, j);
                }
            }
        }
        
        return product;
        
    };
};

template <typename element_type, std::size_t ROWS, std::size_t COLUMNS>
class Matrix : public DynamicMatrix<element_type>{
public:
    inline void setNColumns(std::size_t& c) = delete;
    inline void setNRows(std::size_t& r) = delete;
    inline void resize(std::size_t r, std::size_t c) = delete;
    
    inline Matrix(const bool fromscratch=false) : DynamicMatrix<element_type>(ROWS, COLUMNS, fromscratch) {};
    
    //Column-based constructor
    inline Matrix(const std::array<akml::Matrix<element_type, ROWS, 1>, COLUMNS>& cols) : DynamicMatrix<element_type>(ROWS, COLUMNS, false) {
        if (cols.size() == 0)
            throw std::invalid_argument("Empty initialize list.");
        
        for (std::size_t line(0); line < this->rows; line++){
            for (std::size_t col(0); col < this->columns; col++){
                *(this->m_data + line*(this->columns)+col) = cols[col].read(line+1, 1);
            }
        }
        
        //std::cout << "A matrix is born " << this << " - Mdata ref " << m_data << std::endl;
    }
    
    inline Matrix(const std::array <std::array <element_type, COLUMNS>, ROWS>& data) : DynamicMatrix<element_type>(data) {}
    
    //copy constructor
    inline Matrix(const Matrix<element_type, ROWS, COLUMNS>& other) : DynamicMatrix<element_type>(ROWS, COLUMNS, true) {
        if (other.getNColumns() != COLUMNS || other.getNRows() != ROWS)
            throw std::invalid_argument("Matrix should be equally sized to be assignable.");
        
        std::copy(other.getStorage(), other.getStorageEnd(), this->m_data);
    }
    
    inline Matrix(std::function<element_type(element_type, std::size_t, std::size_t)>& transfunc) : DynamicMatrix<element_type>(ROWS, COLUMNS, transfunc) {}
    
    inline Matrix(std::function<element_type(element_type)>& transfunc) : DynamicMatrix<element_type>(ROWS, COLUMNS, transfunc) {}
 
    inline Matrix(Matrix<element_type, ROWS, COLUMNS>&& other) : DynamicMatrix<element_type>(ROWS, COLUMNS) {
        if (other.getNColumns() != COLUMNS || other.getNRows() != ROWS)
            throw std::invalid_argument("Matrix should be equally sized to be assignable.");
        
        //this->m_data = other.getStorage();
        this->setMDataPointer(other.getStorage());
        other.getStorage() = nullptr;
    }
    
    inline Matrix<element_type, ROWS, COLUMNS>& operator=(const Matrix<element_type, ROWS, COLUMNS>& other){
        if (!this->isInitialized())
            this->createInternStorage();
        
        if (other.isInitialized())
            std::copy(other.getStorage(), other.getStorageEnd(), this->m_data);
        return *this;
    }
    
    inline Matrix<element_type, ROWS, COLUMNS>& operator=(const DynamicMatrix<element_type>& other){
        DynamicMatrix<element_type>::operator=(other);
        return *this;
    }
     
    inline friend void operator<<(std::ostream& flux, Matrix<element_type, ROWS, COLUMNS> matrix){
        if (!matrix.isInitialized())
            throw std::invalid_argument("Matrix provided is not initialized.");
        
        return DynamicMatrix<element_type>::cout(matrix);
    }
    
    inline Matrix<element_type, ROWS, COLUMNS>& operator+=(const Matrix<element_type, ROWS, COLUMNS>& mat){
        if (!mat.isInitialized())
            throw std::invalid_argument("Matrix provided is not initialized.");
        
        if (mat.getNColumns() != (this->columns) || mat.getNRows() != (this->rows))
            throw std::invalid_argument("Matrices should be equally sized to be assignable.");
        
        for (std::size_t i(0); i < (this->rows)*(this->columns); i++){
            *(this->getInternElement(i)) += *(mat.getInternElement(i));
        }

        return *this;
    }
    
    inline friend Matrix<element_type, ROWS, COLUMNS> operator+(Matrix<element_type, ROWS, COLUMNS> selfmat, const Matrix<element_type, ROWS, COLUMNS>& mat){
        selfmat += mat;
        return selfmat;
    }
    
    inline friend Matrix<element_type, ROWS, COLUMNS> operator*(Matrix<element_type, ROWS, COLUMNS> selfmat, const Matrix<element_type, ROWS, COLUMNS>& mat){
        selfmat *= mat;
        return selfmat;
    }
    
    template <typename NumericType>
    inline friend Matrix<element_type, ROWS, COLUMNS> operator*(Matrix<element_type, ROWS, COLUMNS> selfmat, const NumericType& num){
        static_assert(std::is_arithmetic<NumericType>::value, "NumericType must be numeric");
        if (selfmat.getNColumns() != COLUMNS || selfmat.getNRows() != ROWS)
            throw std::invalid_argument("Matrix should be equally sized to be assignable.");
        for (element_type* tar(selfmat.getStorage()); tar < selfmat.getStorageEnd(); tar++){
            *tar = *tar * num;
        }
        return selfmat;
    }
    
    template <typename NumericType>
    inline friend Matrix<element_type, ROWS, COLUMNS> operator*(const NumericType& num, Matrix<element_type, ROWS, COLUMNS> selfmat){
        return operator*(selfmat, num);
    }
    
    inline Matrix<element_type, ROWS, COLUMNS>& operator*=(const Matrix<element_type, ROWS, COLUMNS>& mat){
        if (!mat.isInitialized())
            throw std::invalid_argument("Matrix provided is not initialized.");
        
        if (mat.getNColumns() != COLUMNS || mat.getNRows() != ROWS)
            throw std::invalid_argument("Matrix should be equally sized to be assignable.");
        
        *this = product(*this, mat);
        return *this;
    }
    
    inline Matrix<element_type, ROWS, COLUMNS>& operator-=(const Matrix<element_type, ROWS, COLUMNS>& mat){
        if (!mat.isInitialized())
            throw std::invalid_argument("Matrix provided is not initialized.");

        for (std::size_t i(0); i < (this->rows)*(this->columns); i++){
            *(this->getInternElement(i)) -= *(mat.getInternElement(i));
        }

        return *this;
    }
    
    inline friend Matrix<element_type, ROWS, COLUMNS> operator-(Matrix<element_type, ROWS, COLUMNS> selfmat, const Matrix<element_type, ROWS, COLUMNS>& mat){
        selfmat -= mat;
        return selfmat;
    }
    
    template <const std::size_t R, const std::size_t DIMTEMP, const std::size_t C>
    inline static Matrix<element_type, R, C> product(const Matrix<element_type, R, DIMTEMP>& A, const Matrix<element_type, DIMTEMP, C>& B){
        if (A.getNColumns() != B.getNRows())
            throw std::invalid_argument("Attempting to perform a non-defined matrix product.");
        if (!A.isInitialized() || !B.isInitialized())
            throw std::invalid_argument("Matrix provided is not initialized.");
        
        Matrix<element_type, R, C> product;
        for (std::size_t i=1; i <= A.getNRows(); i++){
            for (std::size_t j=1; j <= B.getNColumns(); j++){
                for (size_t k=1; k <= A.getNColumns(); k++){
                    product(i, j) += A.read(i, k) * B.read(k, j);
                }
            }
        }
        
        return product;
        
    };
};


template <typename element_type, std::size_t ROWS, std::size_t COLUMNS>
class StaticMatrix : public MatrixInterface<element_type>{
protected:
    element_type m_data_static_array[ROWS*COLUMNS];
    
public:
    inline void create(){
        this->createInternStorage();
        for (std::size_t i(0); i < (this->rows)*(this->columns); i++){
            *(this->m_data + i) = static_cast<element_type>(0);
        }
    }
    
    inline void createInternStorage(){
        if (this->m_data == nullptr)
            this->m_data = &m_data_static_array[0];
        if (this->m_data_end == nullptr)
            this->m_data_end = this->m_data+(this->rows)*(this->columns);
    }
    
    inline void deleteInternStorage(){}
    
    inline StaticMatrix(const bool fromscratch=false) : MatrixInterface<element_type>(ROWS, COLUMNS) {
        (fromscratch) ? this->createInternStorage() : this->create();
        //std::cout << "A matrix is born " << this << " - Mdata ref " << m_data << std::endl;
    }
    
    //Column-based constructor
    inline StaticMatrix(const std::array<akml::StaticMatrix<element_type, ROWS, 1>, COLUMNS>& cols) : MatrixInterface<element_type>(ROWS, COLUMNS) {
        if (cols.size() == 0)
            throw std::invalid_argument("Empty initialize list.");
        
        this->createInternStorage();
        for (std::size_t line(0); line < this->rows; line++){
            for (std::size_t col(0); col < this->columns; col++){
                *(this->m_data + line*(this->columns)+col) = cols[col].read(line+1, 1);
            }
        }
        
        //std::cout << "A matrix is born " << this << " - Mdata ref " << m_data << std::endl;
    }
    
    inline StaticMatrix(const std::array <std::array <element_type, COLUMNS>, ROWS>& data) : MatrixInterface<element_type>(ROWS, COLUMNS) {
        this->createInternStorage();
        for (std::size_t line(0); line < ROWS; line++){
            for (std::size_t col(0); col < COLUMNS; col++){
                *(this->m_data + line*COLUMNS+col) = data[line][col];
            }
        }
        //std::cout << "A matrix is born " << this << " - Mdata ref " << m_data << std::endl;
    }
    
    //copy constructor
    inline StaticMatrix(const StaticMatrix<element_type, ROWS, COLUMNS>& other) : MatrixInterface<element_type>(ROWS, COLUMNS) {
        //deleteInternStorage();
        this->createInternStorage();
        std::copy(other.getStorage(), other.getStorageEnd(), this->m_data);
        //std::cout << "A matrix is born " << this << " - Mdata ref " << m_data << std::endl;
    }
    
    inline StaticMatrix(std::function<element_type(element_type, std::size_t, std::size_t)>& transfunc) : MatrixInterface<element_type>(ROWS, COLUMNS) {
        this->create();
        this->transform(transfunc);
        //std::cout << "A matrix is born " << this << " - Mdata ref " << m_data << std::endl;
    }
    
    inline StaticMatrix(std::function<element_type(element_type)>& transfunc) : MatrixInterface<element_type>(ROWS, COLUMNS) {
        this->create();
        this->transform(transfunc);
        //std::cout << "A matrix is born " << this << " - Mdata ref " << m_data << std::endl;
    }
    
    inline StaticMatrix<element_type, ROWS, COLUMNS>& operator=(StaticMatrix<element_type,ROWS, COLUMNS>&& other){
        if (this != &other)
            std::copy(other.getStorage(), other.getStorageEnd(), this->m_data);
        return *this;
    }
    
    inline StaticMatrix(StaticMatrix<element_type, ROWS, COLUMNS>&& other) : MatrixInterface<element_type>(ROWS, COLUMNS) {
        createInternStorage();
        std::copy(other.getStorage(), other.getStorageEnd(), this->m_data);
    }
    
    inline void operator=(const std::array <std::array <element_type, COLUMNS>, ROWS>& data){
        for (std::size_t line(0); line < ROWS; line++){
            for (std::size_t col(0); col < COLUMNS; col++){
                *(this->getInternElement(line*COLUMNS+col)) = data[line][col];
            }
        }
    }
    
    inline StaticMatrix<element_type, ROWS, COLUMNS>& operator=(const StaticMatrix<element_type, ROWS, COLUMNS>& other){
        std::copy(other.getStorage(), other.getStorageEnd(), this->m_data);
        return *this;
    }
    
    inline friend void operator<<(std::ostream& flux, StaticMatrix<element_type, ROWS, COLUMNS> matrix){
        return StaticMatrix<element_type, ROWS, COLUMNS>::cout(matrix);
    }
    
    
    inline StaticMatrix<element_type, ROWS, COLUMNS>& operator+=(const StaticMatrix<element_type, ROWS, COLUMNS>& mat){
        for (std::size_t i(0); i < (this->rows)*(this->columns); i++){
            *(this->getInternElement(i)) += *(mat.getInternElement(i));
        }
        return *this;
    }
    
    inline friend StaticMatrix<element_type, ROWS, COLUMNS> operator+(StaticMatrix<element_type, ROWS, COLUMNS> selfmat, const StaticMatrix<element_type, ROWS, COLUMNS>& mat){
        selfmat += mat;
        return selfmat;
    }
    
    inline friend StaticMatrix<element_type, ROWS, COLUMNS> operator*(StaticMatrix<element_type, ROWS, COLUMNS> selfmat, const StaticMatrix<element_type, ROWS, COLUMNS>& mat){
        selfmat *= mat;
        return selfmat;
    }
    
    template <typename NumericType>
    inline friend StaticMatrix<element_type, ROWS, COLUMNS> operator*(StaticMatrix<element_type, ROWS, COLUMNS> selfmat, const NumericType& num){
        static_assert(std::is_arithmetic<NumericType>::value, "NumericType must be numeric");
        for (element_type* tar(selfmat.getStorage()); tar < selfmat.getStorageEnd(); tar++){
            *tar = *tar * num;
        }
        return selfmat;
    }
    
    template <typename NumericType>
    inline friend StaticMatrix<element_type, ROWS, COLUMNS> operator*(const NumericType& num, StaticMatrix<element_type, ROWS, COLUMNS> selfmat){
        return operator*(selfmat, num);

    }
    
    inline StaticMatrix<element_type, ROWS, COLUMNS>& operator*=(const StaticMatrix<element_type, ROWS, COLUMNS>& mat){
        *this = product(*this, mat);
        return *this;
    }
    
    inline StaticMatrix<element_type, ROWS, COLUMNS>& operator-=(const StaticMatrix<element_type, ROWS, COLUMNS>& mat){
        if (!mat.isInitialized())
            throw std::invalid_argument("Matrix provided is not initialized.");
        
        for (std::size_t i(0); i < (this->rows)*(this->columns); i++){
            *(this->getInternElement(i)) -= *(mat.getInternElement(i));
        }

        return *this;
    }
    
    inline friend StaticMatrix<element_type, ROWS, COLUMNS> operator-(StaticMatrix<element_type, ROWS, COLUMNS> selfmat, const StaticMatrix<element_type, ROWS, COLUMNS>& mat){
        selfmat -= mat;
        return selfmat;
    }
    
    
    inline void transpose(){
        if ((this->columns) != (this->rows))
            throw std::invalid_argument("Attempting to perform a direct transposition on a non-squared matrix.");
        
        element_type* localdata;
        localdata = new element_type[(this->rows)*(this->columns)];
        
        for (std::size_t i(0); i < (this->rows); i++){
            for (std::size_t j(0); j < (this->columns); j++){
                *(localdata+i*(this->rows)+j) = this->read(j+1, i+1);
            }
        }
        
        for (std::size_t i(0); i < (this->rows)*(this->columns); i++){
            m_data_static_array[i] = *(localdata+i);
        }
        
        delete[] localdata;
    }
    
    template <const std::size_t R, const std::size_t DIMTEMP, const std::size_t C>
    inline static StaticMatrix<element_type, R, C> product(const StaticMatrix<element_type, R, DIMTEMP>& A, const StaticMatrix<element_type, DIMTEMP, C>& B){
        StaticMatrix<element_type, R, C> product;
        for (std::size_t i=1; i <= A.getNRows(); i++){
            for (std::size_t j=1; j <= B.getNColumns(); j++){
                for (size_t k=1; k <= A.getNColumns(); k++){
                    product(i, j) += A.read(i, k) * B.read(k, j);
                }
            }
        }
        return product;
        
    };
};




};


#endif /* Matrix_hpp */
