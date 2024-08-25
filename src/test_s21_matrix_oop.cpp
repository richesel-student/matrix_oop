#include <gtest/gtest.h>

#include <utility>

#include "s21_matrix_oop.h"

// Тестирование конструктора по умолчанию
TEST(S21MatrixTest, DefaultConstructor) {
  S21Matrix matrix;
  EXPECT_EQ(matrix.get_rows(), 3);
  EXPECT_EQ(matrix.get_cols(), 3);
}

// Тестирование конструктора  c параметрами
TEST(S21MatrixTest, ParenthesesTest) {
  S21Matrix matrix_1(1, 1), matrix_2(2, 2), matrix_3(3, 3);
  EXPECT_EQ(matrix_1.get_rows(), 1);
  EXPECT_EQ(matrix_1.get_cols(), 1);
  EXPECT_EQ(matrix_2.get_rows(), 2);
  EXPECT_EQ(matrix_2.get_cols(), 2);
  EXPECT_EQ(matrix_3.get_rows(), 3);
  EXPECT_EQ(matrix_3.get_cols(), 3);
}
TEST(S21MatrixTest, OperatorParenthesesTest) {
  S21Matrix matrix;
  int rows = 3;
  int cols = 4;

  // Создание динамического массива
  double** matrix_ = new double*[rows];
  for (int i = 0; i < rows; i++) {
    matrix_[i] = new double[cols];
  }
  // Заполнение матрицы значениями
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      matrix_[i][j] = 0.0;  // Примерное значение элемента матрицы
    }
  }

  // Создание объекта класса S21Matrix
  S21Matrix matrixObj(rows, cols);

  // Проверка значения элемента матрицы
  EXPECT_EQ(matrixObj(0, 0), 0.0);
  EXPECT_EQ(matrixObj(1, 2), 0.0);
  EXPECT_EQ(matrixObj(2, 3), 0.0);

  // Освобождение памяти
  for (int i = 0; i < rows; ++i) {
    delete[] matrix_[i];
  }
  delete[] matrix_;
}

TEST(S21MatrixTest, GetMatrixTest) {
  // Создаем матрицу 2x2
  S21Matrix matrix(2, 2);
  matrix.set_element(0, 0, 1.0);
  matrix.set_element(0, 1, 2.0);
  matrix.set_element(1, 0, 3.0);
  matrix.set_element(1, 1, 4.0);

  // Получаем указатель на двумерный массив с данными матрицы
  double** result = matrix.getmatrix();

  // Проверяем, что данные скопированы правильно
  EXPECT_DOUBLE_EQ(result[0][0], 1.0);
  EXPECT_DOUBLE_EQ(result[0][1], 2.0);
  EXPECT_DOUBLE_EQ(result[1][0], 3.0);
  EXPECT_DOUBLE_EQ(result[1][1], 4.0);

  // Освобождаем память, выделенную для результата
  for (int i = 0; i < 2; ++i) {
    delete[] result[i];
  }
  delete[] result;
}


TEST(SET_ROWS_COLS, ParenthesesTest) {
  S21Matrix matrix1;
  matrix1.set_rows(1);
  matrix1.set_cols(1);
  EXPECT_EQ(1, matrix1.get_rows());
  EXPECT_EQ(1, matrix1.get_cols());
  matrix1.set_rows(2);
  matrix1.set_cols(2);
  EXPECT_EQ(2, matrix1.get_rows());
  EXPECT_EQ(2, matrix1.get_rows());
  matrix1.set_rows(3);
  matrix1.set_cols(3);
  EXPECT_EQ(3, matrix1.get_rows());
  EXPECT_EQ(3, matrix1.get_rows());
}

TEST(EqMatrixTest, MatrixTRUE) {
  S21Matrix matrix1(2, 2);
  matrix1.set_element(0, 0, 1.0);
  matrix1.set_element(0, 1, 2.0);
  matrix1.set_element(1, 0, 3.0);
  matrix1.set_element(1, 1, 4.0);

  S21Matrix matrix2(2, 2);
  matrix2.set_element(0, 0, 1.0);
  matrix2.set_element(0, 1, 2.0);
  matrix2.set_element(1, 0, 3.0);
  matrix2.set_element(1, 1, 4.0);

  EXPECT_TRUE(matrix1.EqMatrix(matrix2));
}

TEST(EqMatrix, FALSE_Matrix) {
  S21Matrix matrix1(2, 2);
  matrix1.set_element(0, 0, 1.0);
  matrix1.set_element(0, 1, 2.0);
  matrix1.set_element(1, 0, 3.0);
  matrix1.set_element(1, 1, 4.0);

  S21Matrix matrix2(2, 2);
  matrix2.set_element(0, 0, 1.0);
  matrix2.set_element(0, 1, 2.0);
  matrix2.set_element(1, 0, 5.0);  // Различие в одной ячейке

  try {
    matrix2.set_element(
        1, 1, 6.0);  // Попытка установить значение в недопустимой ячейке
  } catch (const char* msg) {
    // Обработка исключения
    std::cout << "Поймано исключение: " << msg << std::endl;
  }

  EXPECT_FALSE(matrix1.EqMatrix(matrix2));
}
TEST(S21MatrixTest, SetElementValidIndex) {
  S21Matrix myMatrix(3, 3);
  myMatrix.set_element(1, 1, 4.5);
  EXPECT_EQ(4.5, myMatrix.get_element(1, 1));
}

TEST(CopyConstructorTest, CorrectlyCopiesMatrix) {
  S21Matrix original(2, 2);  // Создаем и заполняем исходную матрицу
  original.set_element(0, 0, 1);
  original.set_element(0, 1, 2);
  original.set_element(1, 0, 3);
  original.set_element(1, 1, 4);

  S21Matrix copy(original);  // Создаем копию исходной матрицы

  // Проверяем, что размерности копии совпадают с исходной матрицей
  EXPECT_EQ(original.get_rows(), copy.get_rows());
  EXPECT_EQ(original.get_cols(), copy.get_cols());

  // Проверяем, что значения в копии совпадают с исходной матрицей
  for (int i = 0; i < original.get_rows(); i++) {
    for (int j = 0; j < original.get_cols(); j++) {
      EXPECT_EQ(original.getValue(i, j), copy.getValue(i, j));
    }
  }
}

TEST(CopyConstructorTest, CorrectlyHandlesEmptyMatrix) {
  S21Matrix original(0, 0);  // Создаем пустую исходную матрицу
  S21Matrix copy(original);  // Создаем копию пустой матрицы

  // Проверяем, что размерности копии совпадают с исходной матрицей
  EXPECT_EQ(original.get_rows(), copy.get_rows());
  EXPECT_EQ(original.get_cols(), copy.get_cols());
}

TEST(CopyConstructorTest, CorrectlyHandlesLargeMatrix) {
  S21Matrix original(100, 100);  // Создаем большую исходную матрицу
  // Заполняем исходную матрицу какими-то значениями

  S21Matrix copy(original);  // Создаем копию большой матрицы

  // Проверяем, что размерности копии совпадают с исходной матрицей
  EXPECT_EQ(original.get_rows(), copy.get_rows());
  EXPECT_EQ(original.get_cols(), copy.get_cols());

  // Проверяем, что значения в копии совпадают с исходной матрицей
  for (int i = 0; i < original.get_rows(); i++) {
    for (int j = 0; j < original.get_cols(); j++) {
      EXPECT_EQ(original.getValue(i, j), copy.getValue(i, j));
    }
  }
}

TEST(MoveConstructorTest, CorrectlyMovesMatrix) {
  S21Matrix original(2, 2);  // Создаем исходную матрицу
  // Заполняем исходную матрицу значениями

  S21Matrix copy(original);  // Создаем копию исходной матрицы
  S21Matrix moved(std::move(copy));  // Перемещаем матрицу из копии

  // Проверяем, что размерности и значения перемещенной матрицы совпадают с
  // исходной
  ASSERT_EQ(original.get_rows(), moved.get_rows());
  ASSERT_EQ(original.get_cols(), moved.get_cols());
  for (int i = 0; i < original.get_rows(); i++) {
    for (int j = 0; j < original.get_cols(); j++) {
      ASSERT_EQ(original(i, j), moved(i, j));
    }
  }
  // Проверяем, что у копии после перемещения размерности равны 0
  ASSERT_EQ(0, copy.get_rows());
  ASSERT_EQ(0, copy.get_cols());
}

TEST(MatrixOperatorPlusAssign, SumMatrix) {
  int rows = 9;
  int cols = 10;
  S21Matrix m1(rows, cols);
  S21Matrix m2(rows, cols);
  S21Matrix m3(rows, cols);

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      double a = 1.23 + i * 16.123 + j * 4.123;
      double b = 2.1 + i * 19.123 + j * 14.123;
      m1(i, j) = a;
      m2(i, j) = b;
      m3(i, j) = a + b;
    }
  }

  m1 += m2;
  EXPECT_EQ(m1.get_rows(), m3.get_rows());
  EXPECT_EQ(m1.get_cols(), m3.get_cols());
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      EXPECT_EQ(m1(i, j), m3(i, j));
    }
  }
}

TEST(SumMatrixTest, Testsum) {
  S21Matrix matrixsum_1(2, 2), matrixsum_2(2, 2);
  matrixsum_1.set_element(0, 0, 1.0);
  matrixsum_1.set_element(0, 1, 2.0);
  matrixsum_1.set_element(1, 0, 3.0);
  matrixsum_1.set_element(1, 1, 4.0);
  matrixsum_2.set_element(0, 0, 1.0);
  matrixsum_2.set_element(0, 1, 2.0);
  matrixsum_2.set_element(1, 0, 3.0);
  matrixsum_2.set_element(1, 1, 4.0);

  S21Matrix expected(2, 2);
  expected.set_element(0, 0, 2.0);
  expected.set_element(0, 1, 4.0);
  expected.set_element(1, 0, 6.0);
  expected.set_element(1, 1, 8.0);

  matrixsum_1.SumMatrix(matrixsum_2);

  EXPECT_TRUE(expected.EqMatrix(matrixsum_1));
}

TEST(SubMatrixTest, Testsum) {
  S21Matrix matrixsub_1(2, 2), matrixsub_2(2, 2);
  matrixsub_1.set_element(0, 0, 1.0);
  matrixsub_1.set_element(0, 1, 2.0);
  matrixsub_1.set_element(1, 0, 3.0);
  matrixsub_1.set_element(1, 1, 4.0);
  matrixsub_2.set_element(0, 0, 1.0);
  matrixsub_2.set_element(0, 1, 2.0);
  matrixsub_2.set_element(1, 0, 3.0);
  matrixsub_2.set_element(1, 1, 4.0);
  S21Matrix expected_sub(2, 2);
  matrixsub_1.SubMatrix(matrixsub_2);
  expected_sub.set_element(0, 0, 0.0);
  expected_sub.set_element(0, 1, 0.0);
  expected_sub.set_element(1, 0, 0.0);
  expected_sub.set_element(1, 1, 0.0);
  EXPECT_TRUE(expected_sub.EqMatrix(matrixsub_1));
}

TEST(S21MatrixTest, SizemarixSumSub) {
  S21Matrix matrix1(2, 2);
  S21Matrix matrix2(3, 3);

  EXPECT_THROW(matrix1.SumMatrix(matrix2), std::invalid_argument);
  EXPECT_THROW(matrix1.SubMatrix(matrix2), std::invalid_argument);
}

TEST(MulMatrixTest, Testsum) {
  S21Matrix matrixmul_1(3, 3), matrixmul_resul(3, 3);
  matrixmul_1.set_element(0, 0, 4.0);
  matrixmul_1.set_element(0, 1, 5.0);
  matrixmul_1.set_element(0, 2, -7.0);
  matrixmul_1.set_element(1, 0, 7.3);
  matrixmul_1.set_element(1, 1, 14.6);
  matrixmul_1.set_element(1, 2, 3.1);
  matrixmul_1.set_element(2, 0, 1.0);
  matrixmul_1.set_element(2, 1, 32.0);
  matrixmul_1.set_element(2, 2, 32.0);
  matrixmul_1.MulNumber(0.0);
  EXPECT_TRUE(matrixmul_resul.EqMatrix(matrixmul_resul));
}


TEST(MatrixOperatorMinus, OtherRows) {
  S21Matrix m1(11, 12);
  S21Matrix m2(10, 12);

  EXPECT_THROW(m1 - m2, std::invalid_argument);
  EXPECT_THROW(m2 - m1, std::invalid_argument);
}

TEST(MatrixOperatorMinus, OtherCols) {
  S21Matrix m1(10, 11);
  S21Matrix m2(10, 12);

  EXPECT_THROW(m1 - m2, std::invalid_argument);
  EXPECT_THROW(m2 - m1, std::invalid_argument);
}

TEST(MatrixOperatorMinus, OtherRowsCols) {
  S21Matrix m1(11, 12);
  S21Matrix m2(10, 1);

  EXPECT_THROW(m1 - m2, std::invalid_argument);
  EXPECT_THROW(m2 - m1, std::invalid_argument);
}

TEST(MatrixMulNum, TestValue) {
  S21Matrix mat(2, 2);
  mat.set_element(0, 0, 1.0);
  mat.set_element(0, 1, 2.0);
  mat.set_element(1, 0, 3.0);
  mat.set_element(1, 1, 4.0);
  S21Matrix resmulnum = mat * 2.0;

  EXPECT_EQ(resmulnum.get_element(0, 0), 2.0);
  EXPECT_EQ(resmulnum.get_element(0, 1), 4.0);
  EXPECT_EQ(resmulnum.get_element(1, 0), 6.0);
  EXPECT_EQ(resmulnum.get_element(1, 1), 8.0);
}
TEST(MatrixMulNum, TestZero) {
  S21Matrix mat(2, 2);
  mat.set_element(0, 0, 1.0);
  mat.set_element(0, 1, 2.0);
  mat.set_element(1, 0, 3.0);
  mat.set_element(1, 1, 4.0);
  S21Matrix resmulnum = mat * 0.0;

  EXPECT_EQ(resmulnum.get_element(0, 0), 0.0);
  EXPECT_EQ(resmulnum.get_element(0, 1), 0.0);
  EXPECT_EQ(resmulnum.get_element(1, 0), 0.0);
  EXPECT_EQ(resmulnum.get_element(1, 1), 0.0);
}

TEST(MatrixMulNum, TestNegative) {
  S21Matrix mat(2, 2);
  mat.set_element(0, 0, 1.0);
  mat.set_element(0, 1, 2.0);
  mat.set_element(1, 0, 3.0);
  mat.set_element(1, 1, 4.0);
  S21Matrix resmulnum = mat * (-2.0);

  EXPECT_EQ(resmulnum.get_element(0, 0), -2.0);
  EXPECT_EQ(resmulnum.get_element(0, 1), -4.0);
  EXPECT_EQ(resmulnum.get_element(1, 0), -6.0);
  EXPECT_EQ(resmulnum.get_element(1, 1), -8.0);
}

TEST(Matrixoperatoreq, TestEqTrue) {
  S21Matrix mat_operatoreq(2, 2), mat_operatoreq2(2, 2);
  mat_operatoreq.set_element(0, 0, 1.0);
  mat_operatoreq.set_element(0, 1, 2.0);
  mat_operatoreq.set_element(1, 0, 3.0);
  mat_operatoreq.set_element(1, 1, 4.0);

  mat_operatoreq2.set_element(0, 0, 1.0);
  mat_operatoreq2.set_element(0, 1, 2.0);
  mat_operatoreq2.set_element(1, 0, 3.0);
  mat_operatoreq2.set_element(1, 1, 4.0);

  EXPECT_TRUE(mat_operatoreq == mat_operatoreq);
}

TEST(Matrixoperatoreq, TestEqFalse) {
  S21Matrix mat_operatoreq(2, 2), mat_operatoreq2(2, 2);
  mat_operatoreq.set_element(0, 0, 1.0);
  mat_operatoreq.set_element(0, 1, 2.0);
  mat_operatoreq.set_element(1, 0, 3.0);
  mat_operatoreq.set_element(1, 1, 4.0);

  mat_operatoreq2.set_element(0, 0, 2.0);
  mat_operatoreq2.set_element(0, 1, 7.0);
  mat_operatoreq2.set_element(1, 0, 15.0);
  mat_operatoreq2.set_element(1, 1, 21.0);

  EXPECT_TRUE(mat_operatoreq == mat_operatoreq);
}

TEST(SetMatrixT, Testset) {
  S21Matrix mat_val1(2, 2), mat_val2(2, 2), res_set(2, 2);
  mat_val1.set_element(0, 0, 12.0);
  mat_val1.set_element(0, 1, 2.0);
  mat_val1.set_element(1, 0, 3.0);
  mat_val1.set_element(1, 1, 4.0);
  mat_val1.SetSwap(mat_val2);
  EXPECT_TRUE(mat_val1 == res_set);
}

TEST(MatrixMulMatrix, IncompatibleMatrices) {
  S21Matrix m1(19, 16);
  S21Matrix m2(18, 17);

  EXPECT_THROW(m1.MulMatrix(m2), std::invalid_argument);

  S21Matrix a1(14, 15);
  S21Matrix a2(14, 15);
  EXPECT_THROW(m1.MulMatrix(m2), std::invalid_argument);
}

TEST(MatrixMulMatrix, SquareMatrices) {
  double a[3][3] = {{-1.0, 2.0, 5.0}, {3.0, 4.0, 6.0}, {-8.0, 2.0, 12.0}};
  double b[3][3] = {{-2.0, 2.0, 19.1}, {5.0, 7.0, 17.7}, {-1.0, 4.0, -13.56}};
  double c[3][3] = {
      {7.0, 32.0, -51.5}, {8.0, 58.0, 46.74}, {14.0, 46.0, -280.12}};

  S21Matrix m1(3, 3);
  S21Matrix m2(3, 3);
  S21Matrix m3(3, 3);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      m1(i, j) = a[i][j];
      m2(i, j) = b[i][j];
      m3(i, j) = c[i][j];
    }
  }

  m1.MulMatrix(m2);
  EXPECT_EQ(m1.get_rows(), m3.get_rows());
  EXPECT_EQ(m1.get_cols(), m3.get_cols());
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      EXPECT_NEAR(m1(i, j), m3(i, j), S21Matrix::kEps);
    }
  }
}

TEST(MatrixMulMatrix, RectangleMatrices) {
  double a[3][4] = {
      {-1.0, 2.0, 5.0, 78.45}, {3.0, 4.0, 6.0, 19.01}, {-8.0, 2.0, 12.0, 0.43}};
  double b[4][5] = {{-2.0, 2.0, 19.1, 0.5, 0.001},
                    {5.0, 7.0, 17.7, -0.9, -18.78},
                    {-1.0, 4.0, -13.56, 189.1, 19.43},
                    {18.1, 0.3, -17.1, 1983.14, 0.93}};
  double c[3][5] = {{1426.945, 55.535, -1392.995, 156520.533, 132.5475},
                    {352.081, 63.703, -278.331, 38831.9914, 59.1423},
                    {21.783, 46.129, -287.473, 3116.1502, 195.9919}};

  S21Matrix m1(3, 4);
  S21Matrix m2(4, 5);
  S21Matrix m3(3, 5);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 4; ++j) {
      m1(i, j) = a[i][j];
    }
  }
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 5; ++j) {
      m2(i, j) = b[i][j];
    }
  }
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 5; ++j) {
      m3(i, j) = c[i][j];
    }
  }

  m1.MulMatrix(m2);
  EXPECT_EQ(m1.get_rows(), m3.get_rows());
  EXPECT_EQ(m1.get_cols(), m3.get_cols());
  for (int i = 0; i < m3.get_rows(); ++i) {
    for (int j = 0; j < m3.get_cols(); ++j) {
      EXPECT_NEAR(m1(i, j), m3(i, j), S21Matrix::kEps);
    }
  }
}

TEST(TestTranspoe, test1) {
  S21Matrix mat_Transpoe(2, 2), res_Transpoe(2, 2);
  mat_Transpoe.set_element(0, 0, 1.0);
  mat_Transpoe.set_element(0, 1, 2.0);
  mat_Transpoe.set_element(1, 0, 3.0);
  mat_Transpoe.set_element(1, 1, 4.0);
  mat_Transpoe = mat_Transpoe.Transpose();
  // EXPECT_EQ(resmulnum.get_element(0, 0), -2.0);

  EXPECT_EQ(mat_Transpoe.get_element(0, 0), 1.0);
  EXPECT_EQ(mat_Transpoe.get_element(0, 1), 3.0);
  EXPECT_EQ(mat_Transpoe.get_element(1, 0), 2.0);
  EXPECT_EQ(mat_Transpoe.get_element(1, 1), 4.0);

  res_Transpoe.set_element(0, 0, 1.0);
  res_Transpoe.set_element(0, 1, 3.0);
  res_Transpoe.set_element(1, 0, 2.0);
  res_Transpoe.set_element(1, 1, 4.0);

  EXPECT_TRUE(mat_Transpoe == res_Transpoe);
}


TEST(TestTranspoe, test2) {
  S21Matrix mat_Transpoe(2, 3), res_Transpoe(3, 2);
  mat_Transpoe.set_element(0, 0, 1.0);
  mat_Transpoe.set_element(0, 1, 2.0);
  mat_Transpoe.set_element(0, 2, 3.0);
  mat_Transpoe.set_element(1, 0, 4.0);
  mat_Transpoe.set_element(1, 1, 5.0);
  mat_Transpoe.set_element(1, 2, 6.0);

  
  EXPECT_THROW(mat_Transpoe.Transpose(), const char*);
}
TEST(MatrixDeterminant, Matrix1x1) {
  double dt = 1183.2019381738;
  S21Matrix m(1, 1);

  m(0, 0) = dt;
  EXPECT_NEAR(dt, m.Determinant(), S21Matrix::kEps);
}

TEST(MatrixDeterminant, Matrix2x2) {
  double a[2][2] = {{179.38, 18.91}, {2.18821, 472.9428}};
  double dt = 84795.1004129;
  S21Matrix m(2, 2);

  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      m(i, j) = a[i][j];
    }
  }
  EXPECT_NEAR(dt, m.Determinant(), S21Matrix::kEps);
}

TEST(MatrixDeterminant, Matrix3x3) {
  double a[3][3] = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {5.0, 7.0, 9.0}};
  double dt = 0.0;
  S21Matrix m(3, 3);

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      m(i, j) = a[i][j];
    }
  }
  EXPECT_NEAR(dt, m.Determinant(), S21Matrix::kEps);
}

TEST(MatrixDeterminant, Matrix5x5) {
  double a[5][5] = {{3.0, 2.0, -6.0, 2.0, -6.0},
                    {-4.0, 17.0, 7.0, 17.0, 7.0},
                    {1.0, 2.0, 9.0, -3.0, 4.0},
                    {12.0, 3.0, 3.0, 2.0, 9.0},
                    {-1.0, -2.0, 4.0, 8.0, -1.0}};
  double dt = -158255.0;

  S21Matrix m(5, 5);
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
      m(i, j) = a[i][j];
    }
  }
  EXPECT_NEAR(dt, m.Determinant(), S21Matrix::kEps);
}

TEST(MatrixDeterminant, Matrix6x6) {
  double a[6][6] = {{1.1, 1.2, 1.3, 1.4, 1.5, 1.6},
                    {2.8, -2.9, -2.3, -2.4, 2.5, 2.7},
                    {3.33, 3.2, -3.87, 3.99, 3.47, -3.02},
                    {4.85, 4.23, 4.32, -4.18, 4.89, 4.23},
                    {5.12, 5.32, 5.28, 5.67, -5.73, 5.91},
                    {6.15, -6.53, 6.44, 6.32, 6.78, 6.98}};
  double dt = -77591.0 - (269266237810933.0 / 3733527061589101.0);

  S21Matrix m(6, 6);
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      m(i, j) = a[i][j];
    }
  }
  EXPECT_NEAR(dt, m.Determinant(), S21Matrix::kEps);
}

TEST(MatrixDeterminant, Matrix6x6_2) {
  double a[6][6] = {{0.0, 1.2, 1.3, 1.4, 1.5, 1.6},
                    {0.0, -2.9, -2.3, -2.4, 2.5, 2.7},
                    {0.0, 3.2, -3.87, 3.99, 3.47, -3.02},
                    {0.0, 4.23, 4.32, -4.18, 4.89, 4.23},
                    {0.0, 5.32, 5.28, 5.67, -5.73, 5.91},
                    {0.0, -6.53, 6.44, 6.32, 6.78, 6.98}};
  double dt = 0.0;

  S21Matrix m(6, 6);
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      m(i, j) = a[i][j];
    }
  }
  EXPECT_NEAR(dt, m.Determinant(), S21Matrix::kEps);
}

TEST(MatrixDeterminant, RectangleMatrix) {
  S21Matrix m1(1, 2);
  S21Matrix m2(115, 23);

  EXPECT_THROW(m1.Determinant(), std::invalid_argument);
  EXPECT_THROW(m2.Determinant(), std::invalid_argument);
}

// InverseMatrix

TEST(MatrixInverse, Matrix1x1) {
  double a = 125480.4;
  double b = 5.0 / 627402.0;

  S21Matrix m(1, 1);
  m(0, 0) = a;
  S21Matrix inverse = m.InverseMatrix();
  EXPECT_EQ(inverse.get_rows(), m.get_rows());
  EXPECT_EQ(inverse.get_cols(), m.get_cols());
  EXPECT_NEAR(inverse(0, 0), b, S21Matrix::kEps);
}

TEST(MatrixInverse, Matrix3x3) {
  double a[3][3] = {{2., 5., 7.}, {6., 3., 4.}, {5., -2., -3.}};
  double b[3][3] = {{1., -1., 1.}, {-38., 41., -34.}, {27., -29., 24.}};

  S21Matrix ma(3, 3);
  S21Matrix mb(3, 3);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      ma(i, j) = a[i][j];
      mb(i, j) = b[i][j];
    }
  }

  S21Matrix inverse = ma.InverseMatrix();
  EXPECT_TRUE(inverse == mb);
}

TEST(MatrixInverse, Matrix6x6) {
  double a[6][6] = {
      {12., 47., 58., 47., 42., 14.},   {47., 59., 63., 54., 89., 12.},
      {15., 56., -65., -97., 32., 16.}, {58., 75., -24., 45., -16., 8.},
      {78., 93., 17., 13., -11., 7.},   {84., 65., 19., -35., 54., 18.}};
  double b[6][6] = {{-0.01, -0.001, -0.009, 0.006, -0.005, 0.016},
                    {0.0, 0.008, 0.01, -0.007, 0.018, -0.018},
                    {0.009, -0.005, -0.004, -0.014, 0.011, 0.002},
                    {0.001, 0.003, -0.004, 0.012, -0.009, -0.001},
                    {-0.011, 0.017, 0.005, -0.001, -0.004, -0.005},
                    {0.075, -0.067, -0.012, 0.04, -0.061, 0.058}};

  S21Matrix ma(6, 6);
  S21Matrix mb(6, 6);
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      ma(i, j) = a[i][j];
      mb(i, j) = b[i][j];
    }
  }

  S21Matrix inverse = ma.InverseMatrix();
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      EXPECT_NEAR(inverse(i, j), b[i][j], 1.0e-3);
    }
  }
}

TEST(MatrixInverse, NoExist) {
  double a[3][3] = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {5.0, 7.0, 9.0}};

  S21Matrix m(3, 3);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      m(i, j) = a[i][j];
    }
  }

  EXPECT_THROW(m.InverseMatrix(), std::invalid_argument);
}

TEST(MatrixInverse, NoExist2) {
  S21Matrix m1(1, 1);
  m1(0, 0) = 0.0;
  EXPECT_THROW(m1.InverseMatrix(), std::invalid_argument);

  S21Matrix m2(13, 10);
  EXPECT_THROW(m2.InverseMatrix(), std::invalid_argument);

  S21Matrix m3(9, 10);
  EXPECT_THROW(m3.InverseMatrix(), std::invalid_argument);
}

// Tests for CalcComplemnts

TEST(MatrixCalcComplements, Matrix1x1) {
  S21Matrix m1(1, 1);
  S21Matrix m2(1, 1);

  m1(0, 0) = 0.0;
  S21Matrix comp1 = m1.CalcComplements();
  EXPECT_EQ(comp1.get_rows(), 1);
  EXPECT_EQ(comp1.get_cols(), 1);
  EXPECT_EQ(comp1(0, 0), 1.0);

  m2(0, 0) = 100.0;
  S21Matrix comp2 = m2.CalcComplements();
  EXPECT_EQ(comp2.get_rows(), 1);
  EXPECT_EQ(comp2.get_cols(), 1);
  EXPECT_EQ(comp2(0, 0), 1.0);
}

TEST(MatrixCalcComplements, Matrix2x2) {
  double a[2][2] = {{15.87, 78.98}, {47.25, -45.478}};
  double b[2][2] = {{-45.478, -47.25}, {-78.98, 15.87}};

  S21Matrix ma(2, 2);
  S21Matrix mb(2, 2);
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      ma(i, j) = a[i][j];
      mb(i, j) = b[i][j];
    }
  }

  S21Matrix comp = ma.CalcComplements();
  EXPECT_TRUE(comp == mb);
}

TEST(MatrixCalcComplements, Matrix3x3) {
  double a[3][3] = {{1., 2., 3.}, {0., 4., 2.}, {5., 2., 1.}};
  double b[3][3] = {{0., 10., -20.}, {4., -14., 8.}, {-8., -2., 4.}};

  S21Matrix ma(3, 3);
  S21Matrix mb(3, 3);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      ma(i, j) = a[i][j];
      mb(i, j) = b[i][j];
    }
  }

  S21Matrix comp = ma.CalcComplements();
  EXPECT_TRUE(comp == mb);
}

TEST(MatrixCalcComplements, Matrix4x4) {
  double a[4][4] = {
      {4., 5., 9., 8.}, {4., 1., 2., 3.}, {8., 7., 15., 4.}, {7., 6., 4., 9}};
  double b[4][4] = {{-145., -169., 109., 177.},
                    {252., -504., 72., 108.},
                    {47., 95., 25., -111.},
                    {24., 276., -132., -36.}};

  S21Matrix ma(4, 4);
  S21Matrix mb(4, 4);
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      ma(i, j) = a[i][j];
      mb(i, j) = b[i][j];
    }
  }

  S21Matrix comp = ma.CalcComplements();
  EXPECT_TRUE(comp == mb);
}

TEST(MatrixCalcComplements, Matrix5x5) {
  double a[5][5] = {{78., 951., 147., 47., 52.},
                    {76., 98., 78., 753., -89.},
                    {87., 457., 253., 984., -71.},
                    {47., 453., 786., 123., 357.},
                    {765., -896., 783., 478., 456}};
  double b[5][5] = {
      {892211883., -9088259207., 44376427597., -13166556043., -81751647719.},
      {97617917421., -13672761316., 251606522691., -104032036661.,
       -513616435766.},
      {-71997449493., 10510919457., -193843105045., 72992451018.,
       397773228858.},
      {25486500814., -1504267981., 29580687324., -14989913303., -80792756249.},
      {-12212500158., 1182045334., -9293332343., 4297527901., 19088191207}};

  S21Matrix ma(5, 5);
  S21Matrix mb(5, 5);
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
      ma(i, j) = a[i][j];
      mb(i, j) = b[i][j];
    }
  }

  S21Matrix comp = ma.CalcComplements();
  EXPECT_TRUE(comp == mb);
}

TEST(MatrixCalcComplements, RectangleMatrix) {
  S21Matrix m1(19, 18);
  S21Matrix m2(21, 43);

  EXPECT_THROW(m1.CalcComplements(), std::invalid_argument);
  EXPECT_THROW(m2.CalcComplements(), std::invalid_argument);
}

// Tests for operator-=

TEST(MatrixOperatorMinusAssign, OtherRows) {
  S21Matrix m1(11, 12);
  S21Matrix m2(10, 12);

  EXPECT_THROW(m1 -= m2, std::invalid_argument);
  EXPECT_THROW(m2 -= m1, std::invalid_argument);
}

TEST(MatrixOperatorMinusAssign, OtherCols) {
  S21Matrix m1(10, 11);
  S21Matrix m2(10, 12);

  EXPECT_THROW(m1 -= m2, std::invalid_argument);
  EXPECT_THROW(m2 -= m1, std::invalid_argument);
}

TEST(MatrixOperatorMinusAssign, OtherRowsCols) {
  S21Matrix m1(11, 12);
  S21Matrix m2(10, 1);

  EXPECT_THROW(m1 -= m2, std::invalid_argument);
  EXPECT_THROW(m2 -= m1, std::invalid_argument);
}

TEST(MatrixOperatorMinusAssign, SumMatrix) {
  int rows = 9;
  int cols = 10;
  S21Matrix m1(rows, cols);
  S21Matrix m2(rows, cols);
  S21Matrix m3(rows, cols);

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      double a = 1.23 + i * 16.123 + j * 4.123;
      double b = 2.1 + i * 19.123 + j * 14.123;
      m1(i, j) = a;
      m2(i, j) = b;
      m3(i, j) = a - b;
    }
  }

  m1 -= m2;
  EXPECT_EQ(m1.get_rows(), m3.get_rows());
  EXPECT_EQ(m1.get_cols(), m3.get_cols());
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      EXPECT_EQ(m1(i, j), m3(i, j));
    }
  }
}

// Tests for operator*= num

TEST(MatrixOperatorMulAssignNum, MulMatrixNum) {
  int rows = 9;
  int cols = 10;
  double ratio = -17.123;
  S21Matrix m1(rows, cols);
  S21Matrix m2(rows, cols);

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      m1(i, j) = m2(i, j) = 1.23 + i * 16.123 + j * 4.123;
      m2(i, j) *= ratio;
    }
  }

  m1 *= ratio;
  EXPECT_EQ(m1.get_rows(), m2.get_rows());
  EXPECT_EQ(m1.get_cols(), m2.get_cols());
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      EXPECT_EQ(m1(i, j), m2(i, j));
    }
  }
}

// Tests for operator*= matrix

TEST(MatrixOperatorMulAssignMatrix, IncompatibleMatrices) {
  S21Matrix m1(19, 16);
  S21Matrix m2(18, 17);

  EXPECT_THROW(m1 *= m2, std::invalid_argument);

  S21Matrix a1(14, 15);
  S21Matrix a2(14, 15);
  EXPECT_THROW(m1 *= m2, std::invalid_argument);
}

TEST(MatrixOperatorMulAssignMatrix, SquareMatrices) {
  double a[3][3] = {{-1.0, 2.0, 5.0}, {3.0, 4.0, 6.0}, {-8.0, 2.0, 12.0}};
  double b[3][3] = {{-2.0, 2.0, 19.1}, {5.0, 7.0, 17.7}, {-1.0, 4.0, -13.56}};
  double c[3][3] = {
      {7.0, 32.0, -51.5}, {8.0, 58.0, 46.74}, {14.0, 46.0, -280.12}};

  S21Matrix m1(3, 3);
  S21Matrix m2(3, 3);
  S21Matrix m3(3, 3);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      m1(i, j) = a[i][j];
      m2(i, j) = b[i][j];
      m3(i, j) = c[i][j];
    }
  }

  m1 *= m2;
  EXPECT_EQ(m1.get_rows(), m3.get_rows());
  EXPECT_EQ(m1.get_cols(), m3.get_cols());
  for (int i = 0; i < m3.get_cols(); ++i) {
    for (int j = 0; j < m3.get_cols(); ++j) {
      EXPECT_NEAR(m1(i, j), m3(i, j), S21Matrix::kEps);
    }
  }
}

TEST(MatrixOperatorMulAssignMatrix, RectangleMatrices) {
  double a[3][4] = {
      {-1.0, 2.0, 5.0, 78.45}, {3.0, 4.0, 6.0, 19.01}, {-8.0, 2.0, 12.0, 0.43}};
  double b[4][5] = {{-2.0, 2.0, 19.1, 0.5, 0.001},
                    {5.0, 7.0, 17.7, -0.9, -18.78},
                    {-1.0, 4.0, -13.56, 189.1, 19.43},
                    {18.1, 0.3, -17.1, 1983.14, 0.93}};
  double c[3][5] = {{1426.945, 55.535, -1392.995, 156520.533, 132.5475},
                    {352.081, 63.703, -278.331, 38831.9914, 59.1423},
                    {21.783, 46.129, -287.473, 3116.1502, 195.9919}};

  S21Matrix m1(3, 4);
  S21Matrix m2(4, 5);
  S21Matrix m3(3, 5);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 4; ++j) {
      m1(i, j) = a[i][j];
    }
  }
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 5; ++j) {
      m2(i, j) = b[i][j];
    }
  }
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 5; ++j) {
      m3(i, j) = c[i][j];
    }
  }

  m1 *= m2;
  EXPECT_EQ(m1.get_rows(), m3.get_rows());
  EXPECT_EQ(m1.get_cols(), m3.get_cols());
  for (int i = 0; i < m3.get_rows(); ++i) {
    for (int j = 0; j < m3.get_cols(); ++j) {
      EXPECT_NEAR(m1(i, j), m3(i, j), S21Matrix::kEps);
    }
  }
}

// Tests for operator+

TEST(MatrixOperatorPlus, OtherRows) {
  S21Matrix m1(11, 12);
  S21Matrix m2(10, 12);

  EXPECT_THROW(m1 + m2, std::invalid_argument);
  EXPECT_THROW(m2 + m1, std::invalid_argument);
}

TEST(MatrixOperatorPlus, OtherCols) {
  S21Matrix m1(10, 11);
  S21Matrix m2(10, 12);

  EXPECT_THROW(m1 + m2, std::invalid_argument);
  EXPECT_THROW(m2 + m1, std::invalid_argument);
}

TEST(MatrixOperatorPlus, OtherRowsCols) {
  S21Matrix m1(11, 12);
  S21Matrix m2(10, 1);

  EXPECT_THROW(m1 + m2, std::invalid_argument);
  EXPECT_THROW(m2 + m1, std::invalid_argument);
}


int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
