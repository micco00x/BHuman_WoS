#pragma once

namespace labrob {
namespace qpsolvers {

template <typename Scalar, int numVariables, int numEqualityConstraints, int numInequalityConstraints>
class QPSolver {
 public:
  virtual void solve(
      Scalar* H,
      Scalar* f,
      Scalar* A,
      Scalar* b,
      Scalar* C,
      Scalar* d_min,
      Scalar* d_max) = 0;
  virtual Scalar* get_solution() const = 0;
}; // end class QPSolver

} // end namespace labrob::qpsolvers
} // end namespace labrob