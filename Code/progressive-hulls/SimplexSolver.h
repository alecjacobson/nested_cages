#ifndef SIMPLEX_SOLVER_H
#define SIMPLEX_SOLVER_H
/*
    Simple Simplex Solver Class
    Copyright (C) 2012  Tamas Bolner
	For more information, visit: http://blog.bolner.hu/2012/08/22/simplex/
	
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

// SimplexSolver(SIMPLEX_MINIMIZE,f,C) solves:
//
//     min f' * x
//     subject A x >= b
//
// where C = [A b] is m by n+1 and f is n by 1. This corresponds to:
//
//     x = linprog(f,-A,-b);
//
// Alec Jacobson

#pragma once
#include <Eigen/Dense>
#include <cinttypes>


#define SIMPLEX_MINIMIZE 1
#define SIMPLEX_MAXIMIZE 2

class SimplexSolver {
private:
  Eigen::MatrixXd tableau;
	bool foundSolution;
	double optimum;
  Eigen::VectorXd solution;
	int64_t numberOfVariables;

	int64_t findPivot_min(int64_t column);
	bool simplexAlgorithm(int64_t variableNum);
	int64_t getPivotRow(int64_t column);

protected:

public:
	SimplexSolver(
    int mode, 
    const Eigen::VectorXd &objectiveFunction, 
    const Eigen::MatrixXd &constraints);
	bool hasSolution();
	double getOptimum();
  Eigen::VectorXd getSolution();
};
#endif
