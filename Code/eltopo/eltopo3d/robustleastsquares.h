#ifndef ROBUSTLEASTSQUARES_H
#define ROBUSTLEASTSQUARES_H

#include <Eigen/Core>
#include <Eigen/Sparse>

class RobustLeastSquares
{
public:
    static void solve(const Eigen::SparseMatrix<double> &M, const Eigen::VectorXd &rhs, Eigen::VectorXd &result, bool min2norm=false);
};

#endif // ROBUSTLEASTSQUARES_H
