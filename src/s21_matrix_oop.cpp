#include "s21_matrix_oop.h"

#include <cmath>
#include <iostream>
const double S21Matrix::kEps = 1.0e-6;

S21Matrix::S21Matrix() {
  _rows = 3;
  _cols = 3;
  matrix_ = new double*[_rows];
  for (int i = 0; i < _rows; i++) {
    matrix_[i] = new double[_cols];
    for (int j = 0; j < _cols; j++) {
      matrix_[i][j] = 0.0;
    }
  }
}

S21Matrix::S21Matrix(int rows, int cols) : _rows(rows), _cols(cols) {
  matrix_ = new double*[_rows];
  for (int i = 0; i < _rows; i++) {
    matrix_[i] = new double[_cols];
  }
  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < _cols; j++) {
      matrix_[i][j] = 0.0;
    }
  }
}

S21Matrix::S21Matrix(const S21Matrix& other) {
  this->_rows = other._rows;
  this->_cols = other._cols;

  this->matrix_ = new double*[_rows];
  for (int i = 0; i < _rows; i++) {
    this->matrix_[i] = new double[_cols];
  }
  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < _cols; j++) {
      this->matrix_[i][j] = other.matrix_[i][j];
    }
  }
}

S21Matrix::S21Matrix(S21Matrix&& other) {
  _rows = other._rows;
  _cols = other._cols;
  matrix_ = other.matrix_;

  other._rows = 0;
  other._cols = 0;
  other.matrix_ = nullptr;
}

int S21Matrix::get_rows() const { return _rows; }

int S21Matrix::get_cols() const { return _cols; }
void S21Matrix::set_rows(int rows) {
  if (rows < 1) {
    throw std::invalid_argument("The number of rows is less than one.");
  }
  // this->_rows = rows;
  if (rows != _rows) {
    S21Matrix tmp(rows, _cols);
    tmp.SetSwap(*this);
  }
}
void S21Matrix::set_cols(int cols) {
  this->_cols = cols;
  if (cols < 1) {
    throw std::invalid_argument("The number of rows is less than one.");
  }
  if (cols != _cols) {
    S21Matrix tmp(_rows, cols);
    tmp.SetSwap(*this);
  }
}

double** S21Matrix::getmatrix() const {
  double** res = new double*[_rows];
  for (int i = 0; i < _rows; i++) {
    res[i] = new double[_cols];
    for (int j = 0; j < _cols; j++) {
      res[i][j] = matrix_[i][j];
    }
  }
  return res;
  // Освобождение памяти
  for (int i = 0; i < _rows; ++i) {
    delete[] res[i];
  }
  delete[] res;
}

void S21Matrix::set_element(int rows, int cols, double value) {
  if (rows >= 0 && rows < _rows && cols >= 0 && cols < _cols) {
    matrix_[rows][cols] = value;
  } else
  // {(rows <= 0 || rows != _rows || cols <= 0 || cols != _cols)
  {
    throw "throwОшибка: недопустимые индексы!";
  }
}

double S21Matrix::get_element(int rows, int cols) const {
  if (rows >= 0 && rows < _rows && cols >= 0 && cols < _cols) {
    return matrix_[rows][cols];
  } else {
    std::cout << "Ошибка: недопустимые индексы!";
    return 0;
  }
}
double S21Matrix::getValue(int row, int col) const {
  if (row < 0 || row >= _rows || col < 0 || col >= _cols) {
    // Обработка ошибки: выбросить исключение или вернуть значение по умолчанию
    throw std::out_of_range("Invalid row or column index");
  }
  return matrix_[row][col];
}

S21Matrix& S21Matrix ::operator=(const S21Matrix& other) {
  if (this == &other) {
    return *this;
  }
  for (int i = 0; i < _rows; i++) {
    delete[] matrix_[i];
  }

  delete[] matrix_;

  _rows = other._rows;
  _cols = other._cols;

  matrix_ = new double*[_rows];
  for (int i = 0; i < _rows; i++) {
    matrix_[i] = new double[_cols];
  }

  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < _cols; j++) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }

  return *this;
}

double& S21Matrix ::operator()(int rows, int cols) const {
  return matrix_[rows][cols];
}

S21Matrix S21Matrix::operator+(const S21Matrix& other) {
  // S21Matrix res(_rows, _cols);
  // if(DiffSizeMatrix(other)){
  //     for (int i = 0; i<_rows; i++){
  //         for(int j = 0; j <_cols; j++){
  //             res.matrix_[i][j] = other.matrix_[i][j] + matrix_[i][j];

  //         }
  //     }

  // }
  // else {
  //     throw std::runtime_error("Матрицы разных размеров!");
  // }
  // return res;

  S21Matrix res(_rows, _cols);
  res.SumMatrix(other);
  return res;
}

S21Matrix S21Matrix::operator-(const S21Matrix& other) {
  // S21Matrix res(_rows, _cols);
  // if(DiffSizeMatrix(other)){
  //     for (int i = 0; i<_rows; i++){
  //         for(int j = 0; j <_cols; j++){
  //             res.matrix_[i][j] = other.matrix_[i][j] - matrix_[i][j];

  //         }
  //     }

  // }
  // else {
  //     throw std::runtime_error("Матрицы разных размеров!");
  // }
  // return res;

  S21Matrix res(_rows, _cols);
  res.SubMatrix(other);
  return res;
}

S21Matrix S21Matrix::operator*(double x) {
  // S21Matrix res(_rows, _cols);
  //     for (int i = 0; i<_rows; i++){
  //         for(int j = 0; j <_cols; j++){
  //         res.matrix_[i][j] = matrix_[i][j] * x;
  //     }
  // }
  // return res;

  this->MulNumber(x);
  return *this;
}

S21Matrix S21Matrix::operator*(const S21Matrix& other) {
  this->MulMatrix(other);
  return *this;
}
bool S21Matrix::operator==(const S21Matrix& other) const {
  if (EqMatrix(other)) {
    return true;

  } else {
    return false;
  }
}

S21Matrix& S21Matrix::operator+=(const S21Matrix& other) {
  SumMatrix(other);
  return (*this);
}

S21Matrix& S21Matrix::operator-=(const S21Matrix& other) {
  SubMatrix(other);
  return (*this);
}

S21Matrix& S21Matrix::operator*=(double x) noexcept {
  MulNumber(x);
  return (*this);
}

S21Matrix& S21Matrix::operator*=(const S21Matrix& other) {
  MulMatrix(other);
  return (*this);
}

bool S21Matrix::EqMatrix(const S21Matrix& other) const {
  if (_rows != other._rows || _cols != other._cols) {
    return false;  // Если размеры матриц не совпадают, они точно не равны
  }

  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < _cols; j++) {
      if (matrix_[i][j] != other.matrix_[i][j]) {
        return false;  // Если хотя бы один элемент не совпадает, матрицы не
                       // равны
      }
    }
  }

  return true;
}

void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (_rows != other._rows || _cols != other._cols) {
    throw std::invalid_argument("Pазмеры матриц не совпадают");
    // std::cout << "1"; // Если размеры матриц не совпадают, они точно не равны
  }
  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < _cols; j++) {
      matrix_[i][j] += other.matrix_[i][j];
    }
  }
}
void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (_rows != other._rows || _cols != other._cols) {
    // std::cout << "1"; // Если размеры матриц не совпадают, они точно не равны
    throw std::invalid_argument("Pазмеры матриц не совпадают");
  }
  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < _cols; j++) {
      matrix_[i][j] -= other.matrix_[i][j];
    }
  }
}

void S21Matrix::MulNumber(const double num) {
  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < _cols; j++) {
      matrix_[i][j] *= num;
    }
  }
}

void S21Matrix::SetSwap(S21Matrix& other) {
  std::swap(_rows, other._rows);
  std::swap(_cols, other._cols);
  std::swap(matrix_, other.matrix_);
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  S21Matrix res(_rows, other._cols);

  if (_cols != other._rows) {
    throw std::invalid_argument("The matrices are incompatible.");
  }
  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < other._cols; j++) {
      for (int k = 0; k < _cols; ++k) {
        res.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];
      }
    }
  }

  SetSwap(res);
}

S21Matrix S21Matrix::Transpose() {
  S21Matrix res(_rows, _cols);
  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < _cols; j++) {
      res.set_element(j, i, matrix_[i][j]);
    }
  }

  return res;
}

S21Matrix S21Matrix::Minor(int row, int col) const {
  S21Matrix minor(_cols - 1, _rows - 1);
  for (int i = 0, k = 0; i < _rows - 1; ++i, ++k) {
    for (int j = 0, l = 0; j < _cols - 1; ++j, ++l) {
      if (k == row) {
        ++k;
      }
      if (l == col) {
        ++l;
      }
      minor.matrix_[i][j] = matrix_[k][l];
    }
  }

  return (minor);
}

double S21Matrix::Determinant(void) const {
  if (_rows != _cols) {
    throw std::invalid_argument("The matrix is not square.");
  }

  double det = 0.0;
  if (_rows == 1) {
    det = matrix_[0][0];
  } else if (_rows == 2) {
    det = matrix_[0][0] * matrix_[1][1] - matrix_[0][1] * matrix_[1][0];
  } else if (_rows == 3) {
    det = matrix_[0][0] * matrix_[1][1] * matrix_[2][2] -
          matrix_[0][0] * matrix_[1][2] * matrix_[2][1] -
          matrix_[0][1] * matrix_[1][0] * matrix_[2][2] +
          matrix_[0][1] * matrix_[1][2] * matrix_[2][0] +
          matrix_[0][2] * matrix_[1][0] * matrix_[2][1] -
          matrix_[0][2] * matrix_[1][1] * matrix_[2][0];
  } else {
    for (int j = 0; j < _cols; ++j) {
      double sign = j % 2 ? -1.0 : 1.0;
      det += sign * matrix_[0][j] * Minor(0, j).Determinant();
    }
  }

  return (det);
}

S21Matrix S21Matrix::CalcComplements(void) const {
  if (_rows != _cols) {
    throw std::invalid_argument("The matrix is not square.");
  }

  S21Matrix complements(_rows, _cols);
  if (_rows == 1) {
    complements.matrix_[0][0] = 1.0;
  } else {
    for (int i = 0; i < _rows; ++i) {
      for (int j = 0; j < _cols; ++j) {
        double sign = (i + j) % 2 ? -1.0 : 1.0;
        complements(i, j) = sign * Minor(i, j).Determinant();
      }
    }
  }

  return (complements);
}

S21Matrix S21Matrix::InverseMatrix(void) const {
  if (_rows != _cols) {
    throw std::invalid_argument("The matrix is not square.");
  }

  double det = Determinant();
  if (fabs(det) < kEps) {
    throw std::invalid_argument(
        "The determinant is zero and there is no inverse matrix.");
  }
  det = 1.0 / det;

  return (CalcComplements().Transpose() * det);
}

S21Matrix::~S21Matrix() {
  for (int i = 0; i < _rows; i++) {
    delete[] matrix_[i];
  }
  delete[] matrix_;
}
