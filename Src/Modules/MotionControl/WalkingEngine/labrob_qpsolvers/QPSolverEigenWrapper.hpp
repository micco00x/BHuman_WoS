#pragma once

#include "QPSolver.hpp"

// STL
#include <memory>
#include <iostream>

// Eigen
#include <Eigen/Core>

namespace labrob {
namespace qpsolvers {

template <typename Scalar, int numVariables, int numEqualityConstraints, int numInequalityConstraints>
class QPSolverEigenWrapper {
 public:
  QPSolverEigenWrapper(std::shared_ptr<QPSolver<Scalar, numVariables, numEqualityConstraints, numInequalityConstraints>> qp_solver_ptr) :
      qp_solver_ptr_(qp_solver_ptr) {

  }

  void
  solve(
      Eigen::Matrix<Scalar, numVariables, numVariables>& H,
      Eigen::Matrix<Scalar, numVariables, 1>& g,
      Eigen::Matrix<Scalar, numEqualityConstraints, numVariables>& A,
      Eigen::Matrix<Scalar, numEqualityConstraints, 1>& b,
      Eigen::Matrix<Scalar, numInequalityConstraints, numVariables>& C,
      Eigen::Matrix<Scalar, numInequalityConstraints, 1>& lg,
      Eigen::Matrix<Scalar, numInequalityConstraints, 1>& ug) {
    qp_solver_ptr_->solve(
        H.data(),
        g.data(),
        A.data(),
        b.data(),
        C.data(),
        lg.data(),
        ug.data()
    );
  }

  Eigen::Matrix<Scalar, numVariables, 1>
  get_solution() {
    Scalar* solution_ptr = qp_solver_ptr_->get_solution();
    Eigen::Matrix<Scalar, numVariables, 1> solution;
    for (int idx = 0; idx < numVariables; ++idx) {
      solution(idx) = *(solution_ptr + idx);
    }
    return solution;
  }

 private:
  std::shared_ptr<QPSolver<Scalar, numVariables, numEqualityConstraints, numInequalityConstraints>> qp_solver_ptr_;
}; // end class QPSolverEigenWrapper

} // end namespace qpsolvers
} // end namespace labrob