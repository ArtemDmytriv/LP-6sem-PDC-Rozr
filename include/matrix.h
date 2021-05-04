#pragma once
#include <iostream>
#include <vector>
#include <utility>

class Matrix {
private:
    int row,
        col;
    double *arr;

    void print(std::ostream &os = std::cout) const;
    double at(int n) const { return arr[n]; };
    double& at(int n) { return arr[n]; };
public: // Static vals
    static int out_width;
    static int out_precision;
public:
    explicit Matrix(int n = 0);
    Matrix(int row, int col);
    Matrix(const Matrix &op);
    Matrix(Matrix &&other);
    ~Matrix();
    void init_matrix(double range_start = 0.0, double range_end = 100.0);
    void init_matrix(int range_start, int range_end);
    void init_matrix(const std::vector<std::vector<double>> &vec);
    void init_matrix(std::istream &is);

    Matrix& operator= (const Matrix &op);
    Matrix& operator= (Matrix &&other);

    friend Matrix operator+(const Matrix &op1, const Matrix &op2);
    friend Matrix operator-(const Matrix &op1, const Matrix &op2);
    friend Matrix operator*(const Matrix &op1, const Matrix &op2);
    friend Matrix operator*(const Matrix &op1, double num);
    friend Matrix operator*(double num, const Matrix &op2);

    friend std::ostream& operator<<(std::ostream &os, const Matrix &op);
    friend std::istream& operator>>(std::istream &is, Matrix &op);

    friend Matrix M_x_Mt(const Matrix &mat, const Matrix &transp);

    void set_shape(std::pair<int, int> shape);

    Matrix get_transpose() const;
    // inline getters
    double at(int r, int c) const { return arr[r*col + c]; }
    double& at(int r, int c) { return arr[r*col + c]; }
    int get_row() const { return row; }
    int get_col() const { return col; }
    std::pair<int, int> get_shape() const { return {row, col}; }
    double* data() { return arr; }
    double* data() const { return arr; }
};

