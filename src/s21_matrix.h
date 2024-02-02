#ifndef SRC_S21_MATRIX_H_
#define SRC_S21_MATRIX_H_

#include <cmath>
#include <iostream>

#define DEFAULT_SIZE 1

class S21Matrix {
 private:
  int rows_, cols_;
  double** matrix_;

 public:
  S21Matrix();
  S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix& other);
  S21Matrix(S21Matrix&& other) noexcept;
  ~S21Matrix();

  int GetRows() const;
  void SetRows(int rows);
  int GetCols() const;
  void SetCols(int cols);

  bool EqMatrix(const S21Matrix& other) const;
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix& other);
  S21Matrix Transpose();
  S21Matrix CalcComplements();
  double Determinant();
  S21Matrix InverseMatrix();

  double& operator()(int i, int j) const;
  friend S21Matrix operator+(const S21Matrix& matrix1,
                             const S21Matrix& matrix2);
  friend S21Matrix operator-(const S21Matrix& matrix1,
                             const S21Matrix& matrix2);
  friend S21Matrix operator*(const S21Matrix& matrix1,
                             const S21Matrix& matrix2);
  friend S21Matrix operator*(const S21Matrix& matrix1, const double& num);
  friend S21Matrix operator*(const double& num, const S21Matrix matrix1);
  bool operator==(const S21Matrix& matrix1) const;
  S21Matrix& operator=(const S21Matrix& matrix1);
  S21Matrix& operator=(S21Matrix&& matrix1) noexcept;
  S21Matrix& operator+=(const S21Matrix& matrix1);
  S21Matrix& operator-=(const S21Matrix& matrix1);
  S21Matrix& operator*=(const S21Matrix& matrix1);
  S21Matrix& operator*=(const double num);

  void FreeMatrix(double*** matrix, int rows);
  void InizializeZeroMatrix(double*** matrix, int rows, int cols);
  void MallocMatrix(double*** matrix, int rows, int cols);
  double Minor(const S21Matrix& other);
  void RefreshMatrix(S21Matrix& new_matrix, const S21Matrix& matrix, int n,
                     int ki, int kj);
};

#endif  // SRC_S21_MATRIX_H_