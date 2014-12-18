#include "robustleastsquares.h"
#include "SuiteSparseQR.hpp"
#include <Eigen/SPQRSupport>
#include <iostream>

using namespace Eigen;

void RobustLeastSquares::solve(const SparseMatrix<double> &M, const Eigen::VectorXd &rhs, Eigen::VectorXd &result, bool min2norm)
{
    cholmod_common Common;
    cholmod_l_start(&Common);

    int rows = M.rows();
    int cols = M.cols();

    SparseMatrix<double, ColMajor, UF_long> mat(M);

    cholmod_sparse A = viewAsCholmod(mat);
    cholmod_dense *b = cholmod_l_zeros(rows, 1, CHOLMOD_REAL, &Common);
    cholmod_dense *x;

    memcpy(b->x, &rhs[0], rows*sizeof(double));

    if(min2norm)
        x = SuiteSparseQR_min2norm<double>(SPQR_ORDERING_DEFAULT, SPQR_DEFAULT_TOL, &A, b, &Common);
    else
        x = SuiteSparseQR<double>(&A, b, &Common);

    result.resize(cols);
    memcpy(&result[0], x->x, cols*sizeof(double));

    cholmod_free_dense(&b, &Common);
    cholmod_free_dense(&x, &Common);
    cholmod_l_finish(&Common);
}
