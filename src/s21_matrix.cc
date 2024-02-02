#include "s21_matrix.h"

void S21Matrix::MallocMatrix(double*** matrix, int rows, int cols) {
  *matrix = new double*[rows];
  for (int i = 0; i < rows; i++) {
    (*matrix)[i] = new double[cols];
  }
}
void S21Matrix::InizializeZeroMatrix(double*** matrix, int rows, int cols) {
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      (*matrix)[i][j] = 0;
    }
  }
}
void S21Matrix::FreeMatrix(double*** matrix, int rows) {
  for (int i = 0; i < rows; i++) {
    delete[](*matrix)[i];
  }
  delete[](*matrix);
}

S21Matrix::S21Matrix() {
  rows_ = DEFAULT_SIZE;
  cols_ = DEFAULT_SIZE;
  MallocMatrix(&matrix_, rows_, cols_);
  InizializeZeroMatrix(&matrix_, rows_, cols_);
}

S21Matrix::S21Matrix(int rows, int cols) {
  if (rows > 0 && cols > 0) {
    rows_ = rows;
    cols_ = cols;
    MallocMatrix(&matrix_, rows_, cols_);
    InizializeZeroMatrix(&matrix_, rows_, cols_);
  } else {
    throw ::std::invalid_argument("Arguments aren`t correct!");
  }
}

S21Matrix::S21Matrix(const S21Matrix& other)
    : S21Matrix(other.rows_, other.cols_) {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }
}
S21Matrix::S21Matrix(S21Matrix&& other) noexcept {
  rows_ = other.rows_;
  cols_ = other.cols_;
  matrix_ = other.matrix_;
  other.matrix_ = nullptr;
  other.cols_ = 0;
  other.rows_ = 0;
}

S21Matrix::~S21Matrix() { FreeMatrix(&matrix_, rows_); }

// main function

bool S21Matrix::EqMatrix(const S21Matrix& other) const {
  bool equal = true;
  if (rows_ == other.cols_ && cols_ == other.cols_) {
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        if (fabs(matrix_[i][j] - other.matrix_[i][j]) > 1e-7) equal = false;
      }
    }
  } else {
    equal = false;
  }
  return equal;
}

void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (rows_ == other.rows_ && cols_ == other.cols_) {
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        matrix_[i][j] = matrix_[i][j] + other.matrix_[i][j];
      }
    }
  } else {
    throw std::invalid_argument("Arguments aren`t correct!");
  }
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (rows_ == other.rows_ && cols_ == other.cols_) {
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        matrix_[i][j] = matrix_[i][j] - other.matrix_[i][j];
      }
    }
  } else {
    throw std::invalid_argument("Arguments aren`t correct!");
  }
}

void S21Matrix::MulNumber(const double num) {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = matrix_[i][j] * num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (cols_ == other.rows_ && rows_ == other.cols_) {
    double** new_matrix;
    MallocMatrix(&new_matrix, cols_, other.cols_);
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < other.cols_; j++) {
        int counter = 0;
        double sum = 0;
        while (counter < cols_) {
          sum += matrix_[i][counter] * other.matrix_[counter][j];
          counter++;
        }
        new_matrix[i][j] = sum;
      }
    }
    FreeMatrix(&this->matrix_, cols_);
    matrix_ = new_matrix;
  } else {
    throw std::invalid_argument("Arguments aren`t correct!");
  }
}

S21Matrix S21Matrix::Transpose() {
  S21Matrix matrix = S21Matrix(cols_, rows_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix.matrix_[j][i] = matrix_[i][j];
    }
  }
  return matrix;
}

void S21Matrix::RefreshMatrix(S21Matrix& new_matrix, const S21Matrix& matrix,
                              int n, int ki, int kj) {
  int d1 = 0;
  int d2 = 0;
  for (int i = 0; i < n - 1; i++) {
    if (i == ki) d1 = 1;
    for (int j = 0; j < n - 1; j++) {
      if (j == kj) d2 = 1;
      new_matrix(i, j) = matrix(i + d1, j + d2);
    }
    d2 = 0;
  }
}

double S21Matrix::Minor(const S21Matrix& other) {
  double det = 0;
  if (other.GetRows() == 1)
    det = other(0, 0);
  else if (other.GetRows() == 2)
    det = other(0, 0) * other(1, 1) - other(1, 0) * other(0, 1);
  else {
    S21Matrix new_matrix(other.GetRows() - 1, other.GetRows() - 1);
    double k = 1;
    for (int i = 0; i < other.GetRows(); i++) {
      RefreshMatrix(new_matrix, other, other.GetRows(), 0, i);
      det += k * other(0, i) * Minor(new_matrix);
      k *= -1;
    }
  }
  return det;
}

double S21Matrix::Determinant() {
  if (rows_ == cols_) {
    return Minor(*this);
  } else {
    throw std::invalid_argument("Arguments aren`t correct!");
  }
}

S21Matrix S21Matrix::CalcComplements() {
  S21Matrix result(rows_, cols_);
  if (rows_ == cols_) {
    if (rows_ == 1) {
      result.matrix_[0][0] = matrix_[0][0];
    } else {
      double k = -1;
      for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
          S21Matrix matrix(rows_ - 1, cols_ - 1);
          RefreshMatrix(matrix, *this, rows_, i, j);
          double det = matrix.Determinant();
          k = pow(-1, i + j);
          result.matrix_[i][j] = k * det;
        }
      }
    }
  } else {
    throw std::invalid_argument("Arguments aren`t correct!");
  }
  return result;
}

S21Matrix S21Matrix::InverseMatrix() {
  double det = Determinant();
  if (det != 0) {
    S21Matrix result(rows_, cols_);
    if (cols_ != 1) {
      S21Matrix matrix = CalcComplements();
      result = matrix.Transpose();
      for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
          result.matrix_[i][j] *= (1.0 / det);
        }
      }
    } else {
      result.matrix_[0][0] = 1.0 / det;
    }
    return result;

  } else {
    throw std::invalid_argument("Determinant equal zero!");
  }
}

// get and set

int S21Matrix::GetRows() const { return rows_; }
int S21Matrix::GetCols() const { return cols_; }

void S21Matrix::SetRows(int rows) {
  if (rows > 0) {
    double** new_matrix;
    MallocMatrix(&new_matrix, rows, cols_);
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols_; j++) {
        if (i >= rows_) {
          new_matrix[i][j] = 0;
        } else {
          new_matrix[i][j] = matrix_[i][j];
        }
      }
    }
    FreeMatrix(&matrix_, rows_);
    rows_ = rows;
    matrix_ = new_matrix;
  } else {
    throw std::invalid_argument("Arguments aren`t correct!");
  }
}

void S21Matrix::SetCols(int cols) {
  if (cols > 0) {
    double** new_matrix;
    MallocMatrix(&new_matrix, rows_, cols);
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols; j++) {
        if (j >= cols_) {
          new_matrix[i][j] = 0;
        } else {
          new_matrix[i][j] = matrix_[i][j];
        }
      }
    }
    FreeMatrix(&matrix_, rows_);
    cols_ = cols;
    matrix_ = new_matrix;
  } else {
    throw std::invalid_argument("Arguments aren`t correct!");
  }
}

// operator

double& S21Matrix::operator()(int i, int j) const {
  if (i >= 0 && j >= 0 && i < rows_ && j < cols_) {
    return matrix_[i][j];
  } else {
    throw std::out_of_range("Arguments aren`t correct!");
  }
}

S21Matrix operator+(const S21Matrix& matrix1, const S21Matrix& matrix2) {
  S21Matrix result = S21Matrix(matrix1.GetRows(), matrix1.GetCols());
  result = matrix1;
  result.SumMatrix(matrix2);
  return result;
}
S21Matrix operator-(const S21Matrix& matrix1, const S21Matrix& matrix2) {
  S21Matrix result = S21Matrix(matrix1.GetRows(), matrix1.GetCols());
  result = matrix1;
  result.SubMatrix(matrix2);
  return result;
}
S21Matrix operator*(const S21Matrix& matrix1, const S21Matrix& matrix2) {
  S21Matrix result = S21Matrix(matrix1.GetRows(), matrix1.GetCols());
  result = matrix1;
  result.MulMatrix(matrix2);
  return result;
}
S21Matrix operator*(const S21Matrix& matrix1, const double& num) {
  S21Matrix result = S21Matrix(matrix1.GetRows(), matrix1.GetCols());
  result = matrix1;
  result.MulNumber(num);
  return result;
}
S21Matrix operator*(const double& num, const S21Matrix matrix1) {
  return matrix1 * num;
}
bool S21Matrix::operator==(const S21Matrix& matrix1) const {
  return EqMatrix(matrix1);
}

S21Matrix& S21Matrix::operator=(const S21Matrix& matrix1) {
  if (this != &matrix1) {
    FreeMatrix(&this->matrix_, rows_);
    rows_ = matrix1.rows_;
    cols_ = matrix1.cols_;
    MallocMatrix(&this->matrix_, rows_, cols_);
    for (int i = 0; i < matrix1.rows_; i++) {
      for (int j = 0; j < matrix1.cols_; j++) {
        matrix_[i][j] = matrix1.matrix_[i][j];
      }
    }
  }
  return *this;
}

S21Matrix& S21Matrix::operator=(S21Matrix&& matrix1) noexcept {
  if (this == &matrix1) return *this;
  FreeMatrix(&this->matrix_, rows_);
  this->matrix_ = matrix1.matrix_;
  this->rows_ = rows_;
  this->cols_ = cols_;
  matrix1.rows_ = 0;
  matrix1.cols_ = 0;
  matrix1.matrix_ = nullptr;
  return *this;
}

S21Matrix& S21Matrix::operator+=(const S21Matrix& matrix1) {
  SumMatrix(matrix1);
  return *this;
}
S21Matrix& S21Matrix::operator-=(const S21Matrix& matrix1) {
  SubMatrix(matrix1);
  return *this;
}
S21Matrix& S21Matrix::operator*=(const S21Matrix& matrix1) {
  MulMatrix(matrix1);
  return *this;
}
S21Matrix& S21Matrix::operator*=(const double num) {
  MulNumber(num);
  return *this;
}