// Copyright 2022 Kitaev Pavel

#include <gtest/gtest.h>
#include <tbb/tbb.h>
#include "../../../modules/ver_1/dirichle/dirichle.h"

TEST(Sequential, Test_Sequential_Alg_10x10) {
  int size = 10;
  double eps = 0.01;
  double* matrix = new double[size * size];

  FillingTheMatrix(matrix, size);
  ASSERT_NO_THROW(SequentialAlg(matrix, size, eps));

  delete[] matrix;
}

TEST(Sequential, Test_Sequential_Alg_100x100) {
  int size = 100;
  double eps = 0.01;
  double* matrix = new double[size * size];

  FillingTheMatrix(matrix, size);
  ASSERT_NO_THROW(SequentialAlg(matrix, size, eps));

  delete[] matrix;
}

TEST(TBB, Test_TBB_Alg_10x10) {
  int size = 10;
  double eps = 0.01;
  double* matrix = new double[size * size];

  FillingTheMatrix(matrix, size);
  ASSERT_NO_THROW(ParallelAlgTBB(matrix, size, eps));

  delete[] matrix;
}

TEST(TBB, Test_TBB_Alg_100x100) {
  int size = 100;
  double eps = 0.01;
  double* matrix = new double[size * size];

  FillingTheMatrix(matrix, size);
  ASSERT_NO_THROW(ParallelAlgTBB(matrix, size, eps));

  delete[] matrix;
}

TEST(STD, Test_STD_Alg_10x10) {
  int size = 10;
  double eps = 0.01;
  double* matrix = new double[size * size];

  FillingTheMatrix(matrix, size);
  ASSERT_NO_THROW(ParallelAlgSTD(matrix, size, eps));

  delete[] matrix;
}

TEST(STD, Test_STD_Alg_100x100) {
  int size = 100;
  double eps = 0.01;
  double* matrix = new double[size * size];

  FillingTheMatrix(matrix, size);
  ASSERT_NO_THROW(ParallelAlgSTD(matrix, size, eps));

  delete[] matrix;
}

/*
TEST(Time_test_1500_1500, Test_SEQ_TBB_STD_Alg) {
  int size = 1000;
  double eps = 0.01;

  double* matrix_seq = new double[size * size];
  FillingTheMatrix(matrix_seq, size);
  double start_seq = clock();
  SequentialAlg(matrix_seq, size, eps);
  double end_seq = clock();
  double time_seq =
  static_cast<float>(end_seq - start_seq) / CLOCKS_PER_SEC;
  delete[] matrix_seq;

  double* matrix_tbb = new double[size * size];
  FillingTheMatrix(matrix_tbb, size);
  double start_tbb = clock();
  ParallelAlgTBB(matrix_tbb, size, eps);
  double end_tbb = clock();
  double time_tbb =
  static_cast<float>(end_tbb - start_tbb) / CLOCKS_PER_SEC;
  delete[] matrix_tbb;

  double* matrix_std = new double[size * size];
  FillingTheMatrix(matrix_std, size);
  double start_std = clock();
  ParallelAlgSTD(matrix_std, size, eps);
  double end_std = clock();
  double time_std =
  static_cast<float>(end_std - start_std) / CLOCKS_PER_SEC;
  delete[] matrix_std;

  std::cout << "SEQ: " << time_seq << std::endl;
  std::cout << "TBB: " << time_tbb << std::endl;
  std::cout << "STD: " << time_std << std::endl;

  SUCCEED();
}
*/
