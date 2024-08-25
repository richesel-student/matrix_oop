#ifndef S21_MATRIX_OOP_H
#define S21_MATRIX_OOP_H
#include <exception>
#include <iostream>

class S21Matrix {
 private:
  int _rows, _cols;  
  double** matrix_;  

 public:
  static const double kEps;
  S21Matrix();                        // default constructor
  S21Matrix(int rows, int cols);      // parameterized constructor
  S21Matrix(const S21Matrix& other);  // copy constructor
  S21Matrix(S21Matrix&& other);       // move constructor
  S21Matrix& operator=(const S21Matrix& other);
  ~S21Matrix();

  bool EqMatrix(const S21Matrix& other) const;
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix& other);

  int get_rows() const;
  int get_cols() const;
  void set_rows(int rows);
  void set_cols(int cols);
  double** getmatrix() const;  // установка матрицы
  void set_element(int row, int col, double value);  // установка элемента
  double get_element(int rows, int cols) const;  // вывод элемента
  double getValue(int row, int col) const;
  void SetSwap(S21Matrix& other);
  S21Matrix Transpose();
  S21Matrix CalcComplements(void) const;
  S21Matrix Minor(int row, int col) const;
  double Determinant(void) const;
  S21Matrix InverseMatrix(void) const;

  S21Matrix operator+(const S21Matrix& other);
  S21Matrix operator-(const S21Matrix& other);
  S21Matrix operator*(double x);
  S21Matrix operator*(const S21Matrix& other);
  bool operator==(const S21Matrix& other) const;
  S21Matrix& operator+=(const S21Matrix& other);
  S21Matrix& operator-=(const S21Matrix& other);
  S21Matrix& operator*=(double x) noexcept;
  S21Matrix& operator*=(const S21Matrix& other);
  double& operator()(int rows, int cols) const;
};
#endif
