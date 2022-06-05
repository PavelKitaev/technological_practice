// Copyright 2022 Kitaev Pavel

#include <tbb/tbb.h>
#include <iostream>
#include <cmath>
#include "../../../3rdparty/unapproved/unapproved.h"
#include "../../../modules/task_1/dirichle/dirichle.h"

void PrintMatrix(double* matrix, int size) {
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      std::cout << matrix[i * size + j] << " ";
    }

    std::cout << std::endl;
  }
}

void FillingTheMatrix(double* matrix, int size) {
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      if ((i == 0) || (i == size - 1) || (j == 0) || (j == size - 1)) {
        matrix[i * size + j] = 100;
      } else {
        matrix[i * size + j] = 0;
      }
    }
  }
}

void SequentialAlg(double* matrix, int size, double eps) {
  double* dm = new double[size]{ 0 };
  double dmax, temp, d;
  const int N = size - 2;

  do {
    dmax = 0;

    for (int wave = 1; wave < N + 1; wave++) {
      dm[wave-1] = 0;
      for (int i = 1; i < wave + 1; i++) {
        int j = wave + 1 - i;
        temp = matrix[size * i + j];

        matrix[size * i + j] = 0.25 * (matrix[size * i + j + 1] +
          matrix[size * i + j - 1] + matrix[size * (i + 1) + j] +
          matrix[size * (i - 1) + j]);

        d = fabs(temp - matrix[size * i + j]);
        if (dm[i - 1] < d) {
          dm[i - 1] = d;
        }
      }
    }

    for (int wave = N - 1; wave > 0; wave--) {
      for (int i = N - wave + 1; i < N + 1; i++) {
        int j = 2 * N - wave - i + 1;
        temp = matrix[size * i + j];

        matrix[size * i + j] = 0.25 * (matrix[size * i + j + 1] +
          matrix[size * i + j - 1] + matrix[size * (i + 1) + j] +
          matrix[size * (i - 1) + j]);

        d = fabs(temp - matrix[size * i + j]);
        if (dm[i - 1] < d) {
          dm[i - 1] = d;
        }
      }
    }

    for (int i = 0; i < size; i++) {
      if (dmax < dm[i]) {
        dmax = dm[i];
      }
    }
  } while (dmax > eps);

  delete[] dm;
}

void ParallelAlgTBB(double* matrix, int size, double eps) {
  double* dm = new double[size] { 0 };
  double dmax;
  const int N = size - 2;

  do {
    dmax = 0;

    tbb::parallel_for(tbb::blocked_range<int>(1, N + 1, 1),
      [&](const tbb::blocked_range<int>& range) {
        for (int wave = range.begin(); wave < range.end(); wave++) {
          dm[wave - 1] = 0;
          for (int i = 1; i < wave + 1; i++) {
            int j = wave + 1 - i;
            double temp = matrix[size * i + j];

            matrix[size * i + j] = 0.25 * (matrix[size * i + j + 1] +
              matrix[size * i + j - 1] + matrix[size * (i + 1) + j] +
              matrix[size * (i - 1) + j]);

            double d = fabs(temp - matrix[size * i + j]);

            if (dm[i - 1] < d) {
              dm[i - 1] = d;
            }
          }
        }
      });

    tbb::parallel_for(tbb::blocked_range<int>(0, N - 1, 1),
      [&](const tbb::blocked_range<int>& range) {
        for (int wave = range.end(); wave > range.begin(); wave--) {
          for (int i = N - wave + 1; i < N + 1; i++) {
            int j = 2 * N - wave - i + 1;
            double temp = matrix[size * i + j];

            matrix[size * i + j] = 0.25 * (matrix[size * i + j + 1] +
              matrix[size * i + j - 1] + matrix[size * (i + 1) + j] +
              matrix[size * (i - 1) + j]);

            double d = fabs(temp - matrix[size * i + j]);

            if (dm[i - 1] < d) {
              dm[i - 1] = d;
            }
          }
        }
      });

    for (int i = 0; i < size; i++) {
      if (dmax < dm[i]) {
        dmax = dm[i];
      }
    }
  } while (dmax > eps);

  delete[] dm;
}

void ParallelAlgSTD(double* matrix, int size, double eps) {
  const int th_num = std::thread::hardware_concurrency();
  std::thread* threads = new std::thread[th_num];

  double* dm = new double[size] { 0 };
  double dmax;

  const int delta = size / th_num;
  int residue = size % th_num;

  int i_start_inc, i_end_inc;
  int i_start_dec, i_end_dec;

  do {
    dmax = 0;
    i_start_inc = 1, i_end_inc = 0;
    i_start_dec = (size - 2) - 1, i_end_dec = (size - 2) - 1;

    for (int k = 0; k < th_num; k++) {
      i_end_inc += (delta);
      i_end_dec -= (delta);

      if (residue > 0) {
        i_end_dec--;

        i_end_inc++;
        residue--;
      }

      if (k == th_num - 1) {
        i_end_inc = size - 1;
        i_end_dec = size - 1;
      }

      threads[k] = std::thread([&](int i_start_inc, int i_end_inc,
        int i_start_dec, int i_end_dec) {
          for (int wave = i_start_inc; wave < i_end_inc; wave++) {
            dm[wave - 1] = 0;
            for (int i = 1; i < wave + 1; i++) {
              int j = wave + 1 - i;
              double temp = matrix[size * i + j];

              matrix[size * i + j] = 0.25 * (matrix[size * i + j + 1] +
                matrix[size * i + j - 1] + matrix[size * (i + 1) + j] +
                matrix[size * (i - 1) + j]);

              double d = fabs(temp - matrix[size * i + j]);

              if (dm[i - 1] < d) {
                dm[i - 1] = d;
              }
            }
          }

          for (int wave = i_start_dec; wave > i_end_dec; wave--) {
            for (int i = (size - 2) - wave + 1; i < (size - 2) + 1; i++) {
              int j = 2 * (size - 2) - wave - i + 1;

              double temp = matrix[size * i + j];

              matrix[size * i + j] = 0.25 * (matrix[size * i + j + 1] +
                matrix[size * i + j - 1] + matrix[size * (i + 1) + j] +
                matrix[size * (i - 1) + j]);

              double d = fabs(temp - matrix[size * i + j]);
              if (dm[i - 1] < d) {
                dm[i - 1] = d;
              }
            }
          }
        }, i_start_inc, i_end_inc, i_start_dec, i_end_dec);
      i_start_inc = i_end_inc;
      i_start_dec = i_end_dec;
    }

    for (int i = 0; i < size; i++) {
      if (dmax < dm[i]) {
        dmax = dm[i];
      }
    }

    for (int i = 0; i < th_num; i++) {
      threads[i].join();
    }
  } while (dmax > eps);

  delete[] threads;
  delete[] dm;
}
