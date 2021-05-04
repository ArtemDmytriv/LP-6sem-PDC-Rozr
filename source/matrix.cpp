#include "matrix.h"
#include <iomanip>
#include <algorithm>
#include <fstream>
#include <chrono>
#include <random>

int Matrix::out_width = 9;
int Matrix::out_precision = 1;

Matrix::Matrix(int n) : col(n), row(n) {
    if(n < 0) {
        std::cerr << "Matrix ctor shape error: r=" << row << ", c=" << col << std::endl;
        exit(-1);
    }
    if(n > 0) {
        arr = new double[n*n];
    }
    else {
        arr = nullptr;
    }
}

Matrix::Matrix(int row, int col) : row(row), col(col) {
    if (row < 0 || col < 0) {
        std::cerr << "Matrix ctor shape error: r=" << row << ", c=" << col << std::endl;
        exit(-1);
    }
    if(col > 0 && row > 0) {
        arr = new double[row * col];
    }
    else {
        arr = nullptr;
    }
}

Matrix::Matrix(const Matrix &op) : row(op.row), col(op.col) {
    //std::cout << "#Debug# Copy ctor called\n";
    if (row < 0 || col < 0) {
        std::cerr << "Matrix copy ctor shape error: r=" << row << ", c=" << col << std::endl;
        exit(-1);
    }
    if (col * row > 0) {
        this->arr = new double[row*col];
        std::copy(op.arr, op.arr + row * col, this->arr);
    }
    else {
        arr = nullptr;
    }
}

Matrix::Matrix(Matrix &&other) : arr(other.arr), row(other.row), col(other.col) {
    //std::cout << "#Debug# Move ctor called\n";
    other.arr = nullptr;
    other.row = 0;
    other.col = 0;
}

Matrix::~Matrix() {
    if (arr)
        delete [] arr;
    arr = nullptr;
}

void Matrix::init_matrix(double range_start, double range_end) {
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> dist(range_start, range_end);

    for (int i = 0; i < row*col; ++i) {
        arr[i] = dist(generator);
    }
}

void Matrix::init_matrix(int range_start, int range_end) {
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_int_distribution<int> dist(range_start, range_end);

    for (int i = 0; i < row*col; ++i) {
        arr[i] = dist(generator);
    }
}

void Matrix::init_matrix(const std::vector<std::vector<double>> &vec) {
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            at(i, j) = vec[i][j];
        }
    }
}

void Matrix::init_matrix(std::istream &is) {
    double num = 0.0;
    for (int i = 0; i < col*row; ++i) {
        is >> num;
        arr[i] = num;
    }
}

Matrix& Matrix::operator= (const Matrix &op) {
    //std::cout << "#Debug# Copy assign called\n";
    if (this->row * this->col != op.row * op.col) {
        if (this->arr) {
            delete [] arr;
            arr = nullptr;
        }
        this->arr = new double[op.row * op.col];
    }
    std::copy(op.arr, op.arr + op.row * op.col, this->arr);
    this->row = op.row;
    this->col = op.col;
    return *this;
}

Matrix& Matrix::operator= (Matrix &&other) {
    //std::cout << "#Debug# Move assign called\n";
    if (this != &other) {
        delete [] this->arr;
        this->arr = other.arr;
        this->row = other.row;
        this->col = other.col;
        other.arr = nullptr;
        other.row = 0;
        other.col = 0;
    }
    return *this;
}

Matrix operator+(const Matrix &op1, const Matrix &op2) {
    if (op1.col != op2.col || op1.row != op2.row) {
        std::cerr << "Cannot add: different dimmensions of matrix\n";
        return Matrix();
    }
    Matrix res(op1.row, op1.col);
    for (int i = 0; i < res.row*res.col; ++i) {
        res.at(i) = op1.at(i) + op2.at(i);
    }
    return res;
}

Matrix operator-(const Matrix &op1, const Matrix &op2) {
    if (op1.col != op2.col || op1.row != op2.row) {
        std::cerr << "Cannot add: different dimmensions of matrix\n";
        return Matrix();
    }
    Matrix res(op1.row, op1.col);
    for (int i = 0; i < res.row*res.col; ++i) {
        res.at(i) = op1.at(i) - op2.at(i);
    }
    return res;
}

Matrix operator*(const Matrix &op1, const Matrix &op2) {
    if (op1.col != op2.row) {
        std::cerr << "Cannot mult: incorrect dimensions of matrix\n";
        return Matrix();
    }
    Matrix res(op1.row, op2.col);
    for (int i = 0; i < res.row; ++i) {
        for (int j = 0; j < res.col; ++j) {
            res.at(i, j) = 0;
            for (int k = 0; k < op1.col; ++k)
                res.at(i, j) += op1.at(i, k) * op2.at(k, j);
        }
    }
    return res;
}

Matrix operator*(const Matrix &op1, double num) {
    Matrix res(op1.row, op1.col);
    std::transform(op1.arr, op1.arr + op1.row*op1.col, res.arr, [=](double x) {return x * num;});
    return res;
}

Matrix operator*(double num, const Matrix &op2) {
    return op2 * num;
}

std::ostream& operator<<(std::ostream &os, const Matrix &op) {
    for (int i = 0; i < op.row; ++i) {
        for (int j = 0; j < op.col; ++j) {
            os <<  std::setw(Matrix::out_width) << std::setprecision(Matrix::out_precision) << std::fixed << op.at(i, j);
        }
        os << '\n';
    }
    return os;
}

std::istream& operator>>(std::istream &is, Matrix &op) {
    for (int i = 0; i < op.get_row(); ++i) {
        for (int j = 0; j < op.get_col(); ++j) {
            is >> op.at(i, j);
        }
    }
    return is;
}

Matrix M_x_Mt(const Matrix &mat, const Matrix &transp) {
    if (mat.col != transp.col) {
        std::cerr << "Cannot mult: incorrect dimensions of matrix\n";
        return Matrix();
    }
    Matrix res(mat.row, transp.row);
    for (int i = 0; i < res.row; ++i) {
        for (int j = 0; j < res.col; ++j) {
            res.at(i, j) = 0;
            for (int k = 0; k < mat.col; ++k)
                res.at(i, j) += mat.at(i, k) * transp.at(j, k);
        }
    }
    return res;
}

Matrix Matrix::get_transpose() const {
    Matrix res(col, row);
    for (int i = 0; i < res.row; ++i) {
        for (int j = 0; j < res.col; ++j) {
            res.arr[i*res.col + j] = this->arr[j*res.row + i];
        }
    }
    return res;

}

void Matrix::set_shape(std::pair<int, int> shape) {
    if (shape.first < 0 || shape.second < 0) {
        std::cerr << "set_shape error: r=" << shape.first << ", c=" << shape.second << std::endl;
        return;
    }
    if (shape.first * shape.second > this->row * this->col) {
        if(arr)
            delete [] arr;
        arr = new double[shape.first * shape.second];
    }

    this->row = shape.first;
    this->col = shape.second;
}
