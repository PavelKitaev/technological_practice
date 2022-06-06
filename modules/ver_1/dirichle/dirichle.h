// Copyright 2022 Kitaev Pavel

#ifndef MODULES_VER_1_DIRICHLE_DIRICHLE_H_
#define MODULES_VER_1_DIRICHLE_DIRICHLE_H_

void PrintMatrix(double* matrix, int size);
void FillingTheMatrix(double* matrix, int size);
void SequentialAlg(double* matrix, int size, double eps);
void ParallelAlgTBB(double* matrix, int size, double eps);
void ParallelAlgSTD(double* matrix, int size, double eps);

#endif  // MODULES_VER_1_DIRICHLE_DIRICHLE_H_
