#pragma once

#include <memory>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>
#include <sys/time.h>

#include "labrob_qpsolvers/qpsolvers.hpp"

namespace mpcSolver{

template <int numVariables_, int numEqualityConstraints_, int numInequalityConstraints_>
class MPCSolver{
 public:
  MPCSolver(
      double mpc_timestep,
      double control_timestep,
      double single_support_duration,
      double double_support_duration,
      double com_target_height,
      double foot_constraint_square_width,
      const Eigen::Vector3d& com_position,
      const Eigen::Vector3d& zmp_position) :
      mpcTimeStep(mpc_timestep),
      controlTimeStep(control_timestep),
      singleSupportDuration(single_support_duration),
      doubleSupportDuration(double_support_duration),
      comTargetHeight(com_target_height),
      footContraintSquareWidth(foot_constraint_square_width),
      omega(std::sqrt(9.81 / comTargetHeight)),
      comPos(com_position),
      zmpPos(zmp_position),
      qp_solver_(std::make_shared<labrob::qpsolvers::QPOASESQPSolver<numVariables_, numEqualityConstraints_, numInequalityConstraints_>>()) {
    constexpr int N_ = numVariables_ / 3;

    S = std::round(singleSupportDuration / mpcTimeStep);
    D = std::round(doubleSupportDuration / mpcTimeStep);

    // Matrices for ZMP prediction
    p = Eigen::Matrix<double, N_, 1>::Ones();
    P = Eigen::Matrix<double, N_, N_>::Constant(mpcTimeStep);

    for (int i = 0; i < N_; ++i) {
      for (int j = 0; j < N_; ++j){
        if (j > i) P(i, j) = 0;
      }
    }
  }

  // Main method
  void solve(const std::vector<Eigen::VectorXd>& plan){

    constexpr int N_ = numVariables_ / 3;

    // Mapping
    Eigen::MatrixXd mapping(N_,plan.size());
    for (int j = 0; j < static_cast<int>(plan.size()); j++) {
      for (int i = 0; i < N_; i++) {
        mapping(i,j) = clamp(- (double)i/D + 1 + (double)(S-mpcIter+j*(S+D))/D, 0, 1) - clamp(- (double)i/D + 1 + (double)(S-mpcIter+(j-1)*(S+D))/D, 0, 1);
      }
    }
    for (int i = 1; i < N_; i++) {
      mapping(i,plan.size()-1) = std::max(mapping(i-1,plan.size()-1), mapping(i,plan.size()-1)); //last column must never decrease
    }

    // Generate moving constraint
    Eigen::VectorXd plan_x(plan.size());
    Eigen::VectorXd plan_y(plan.size());
    Eigen::VectorXd plan_z(plan.size());
    Eigen::VectorXd plan_th(plan.size());

    for (int j = 0; j < static_cast<int>(plan.size()); j++) {
      plan_x(j) = plan.at(j)(0);
      plan_y(j) = plan.at(j)(1);
      plan_z(j) = plan.at(j)(2);
      plan_th(j) = plan.at(j)(3);
    } 

    Eigen::VectorXd mc_x = mapping * plan_x;
    Eigen::VectorXd mc_y = mapping * plan_y;
    Eigen::VectorXd mc_z = mapping * plan_z;
    Eigen::VectorXd mc_th = mapping * plan_th;

    Eigen::Matrix<double, N_, N_> rotSin = Eigen::Matrix<double, N_, N_>::Zero();
    Eigen::Matrix<double, N_, N_> rotCos = Eigen::Matrix<double, N_, N_>::Zero();
    for (int i = 0; i < N_; i++) {
      rotSin(i,i) = sin(mc_th(i));
      rotCos(i,i) = cos(mc_th(i));
    }
    Eigen::Matrix<double, numVariables_, numVariables_> zmpRotMatrix = Eigen::Matrix<double, numVariables_, numVariables_>::Identity();
    //zmpRotMatrix <<  rotCos, rotSin, Eigen::MatrixXd::Zero(N_,N_),
    //                -rotSin, rotCos, Eigen::MatrixXd::Zero(N_,N_),
    //                 Eigen::MatrixXd::Zero(N_,2*N_), Eigen::MatrixXd::Identity(N_,N_);

    // Matrices for ZMP prediction
    P = Eigen::Matrix<double, N_, N_>::Ones() * mpcTimeStep;

    for(int i=0; i<N_;++i){
      for(int j=0;j<N_;++j){
        if (j>i) P(i,j)=0;
      }
    }

    AZmp.setZero();
    AZmp.block(     0,      0, N_, N_) = P;
    AZmp.block(    N_,     N_, N_, N_) = P;
    AZmp.block(2 * N_, 2 * N_, N_, N_) = P;
    AZmp = zmpRotMatrix * AZmp;

    Eigen::Matrix<double, numVariables_, 1> bZmpMinFixed;
    Eigen::Matrix<double, numVariables_, 1> bZmpMaxFixed;
    Eigen::Matrix<double, numVariables_, 1> bZmpMinState;
    Eigen::Matrix<double, numVariables_, 1> bZmpMaxState;

    bZmpMinFixed << mc_x - Eigen::Matrix<double, N_, 1>::Constant(footContraintSquareWidth / 2.0), 
                    mc_y - Eigen::Matrix<double, N_, 1>::Constant(footContraintSquareWidth / 2.0), 
                    mc_z - Eigen::Matrix<double, N_, 1>::Constant(footContraintSquareWidth / 2.0);

    bZmpMaxFixed << mc_x + Eigen::Matrix<double, N_, 1>::Constant(footContraintSquareWidth / 2.0), 
                    mc_y + Eigen::Matrix<double, N_, 1>::Constant(footContraintSquareWidth / 2.0), 
                    mc_z + Eigen::Matrix<double, N_, 1>::Constant(footContraintSquareWidth / 2.0);

    bZmpMinState << Eigen::Matrix<double, N_, 1>::Constant(zmpPos(0)), 
                    Eigen::Matrix<double, N_, 1>::Constant(zmpPos(1)), 
                    Eigen::Matrix<double, N_, 1>::Constant(zmpPos(2));

    bZmpMaxState << Eigen::Matrix<double, N_, 1>::Constant(zmpPos(0)), 
                    Eigen::Matrix<double, N_, 1>::Constant(zmpPos(1)), 
                    Eigen::Matrix<double, N_, 1>::Constant(zmpPos(2));

    bZmpMinState = zmpRotMatrix * bZmpMinState;
    bZmpMaxState = zmpRotMatrix * bZmpMaxState;

    bZmpMin = bZmpMinFixed - bZmpMinState;
    bZmpMax = bZmpMaxFixed - bZmpMaxState;

    Eigen::Matrix<double, N_, 1> b;
    Aeq.setZero();

    for(int i=0;i<N_;++i){
      b(i) = pow(exp(-omega*mpcTimeStep),i);
    }

    Aeq.block(0,      0, 1, N_) = (1.0 / omega) * (1.0 - exp(-omega * mpcTimeStep))*b.transpose();
    Aeq.block(1,     N_, 1, N_) = (1.0 / omega) * (1.0 - exp(-omega * mpcTimeStep))*b.transpose();
    Aeq.block(2, 2 * N_, 1, N_) = (1.0 / omega) * (1.0 - exp(-omega * mpcTimeStep))*b.transpose();

    beq << comPos(0) + comVel(0)/omega - zmpPos(0),
          comPos(1) + comVel(1)/omega - zmpPos(1),
          comPos(2) + comVel(2)/omega - (zmpPos(2) + comTargetHeight);

    costFunctionH.setZero();
    costFunctionH.block(     0,      0, N_, N_) = Eigen::Matrix<double, N_, N_>::Identity() + beta_ * P.transpose() * P;
    costFunctionH.block(    N_,     N_, N_, N_) = Eigen::Matrix<double, N_, N_>::Identity() + beta_ * P.transpose() * P;
    costFunctionH.block(2 * N_, 2 * N_, N_, N_) = Eigen::Matrix<double, N_, N_>::Identity() + beta_ * P.transpose() * P;

    costFunctionF.block(     0, 0, N_, 1) = beta_ * P.transpose() * (p * zmpPos.x() - mc_x);
    costFunctionF.block(    N_, 0, N_, 1) = beta_ * P.transpose() * (p * zmpPos.y() - mc_y);
    costFunctionF.block(2 * N_, 0, N_, 1) = beta_ * P.transpose() * (p * zmpPos.z() - mc_z);

    // Solve QP
    qp_solver_.solve(
        costFunctionH,
        costFunctionF,
        Aeq,
        beq,
        AZmp,
        bZmpMin,
        bZmpMax
    );
    Eigen::Matrix<double, numVariables_, 1> decisionVariables = qp_solver_.get_solution();

    // Split the QP solution in ZMP dot and footsteps
    Eigen::Matrix<double, N_, 1> zDotOptimalX;
    Eigen::Matrix<double, N_, 1> zDotOptimalY(N_);
    Eigen::Matrix<double, N_, 1> zDotOptimalZ(N_);

    zDotOptimalX = (decisionVariables.head(N_));
    zDotOptimalY = (decisionVariables.segment(N_,N_));
    zDotOptimalZ = (decisionVariables.segment(2*N_,N_));

    // Update the state based on the result of the QP
    Eigen::Vector3d nextStateX = updateState(zDotOptimalX(0),0,controlTimeStep);
    Eigen::Vector3d nextStateY = updateState(zDotOptimalY(0),1,controlTimeStep);
    Eigen::Vector3d nextStateZ = updateState(zDotOptimalZ(0),2,controlTimeStep);

    comPos << nextStateX(0),nextStateY(0),nextStateZ(0);
    comVel << nextStateX(1),nextStateY(1),nextStateZ(1);
    zmpPos << nextStateX(2),nextStateY(2),nextStateZ(2);

    controlIter = (controlIter + 1) % ((int) std::round(mpcTimeStep / controlTimeStep));
    if (controlIter == 0) {
      mpcIter = (mpcIter + 1) % (S + D);
      footstepCounter++;
    }
  }

  // Get stuff
  const Eigen::Vector3d& getOptimalCoMPosition() { return comPos; }
  const Eigen::Vector3d& getOptimalCoMVelocity() { return comVel; }
  const Eigen::Vector3d& getOptimalZMPPosition() { return zmpPos; }
  
  bool supportFootHasChanged() {
    if (controlIter==0 && footstepCounter>0) return true;
    else return false;
  }

  int getS() { return S; }
  int getD() { return D; }
  int getControlIteration() { return controlIter; }
  unsigned int getFootstepCounter() { return footstepCounter; }
  int getMPCIteration() const { return mpcIter; }

  double getControlTimestep() const { return controlTimeStep;}
  double getMPCTimestep() const { return mpcTimeStep; }
  double getSingleSupportDuration() const { return singleSupportDuration; }

  // Update the state
  Eigen::Vector3d updateState(double zmpDot, int dim, double timeStep) {
    // Update the state along the dim-th direction (0,1,2) = (x,y,z)

    double ch = cosh(omega*timeStep);
    double sh = sinh(omega*timeStep);

    Eigen::Matrix3d A_upd = Eigen::Matrix3d::Zero();
    Eigen::Vector3d B_upd = Eigen::Vector3d::Zero();
    A_upd<<ch,sh/omega,1-ch,omega*sh,ch,-omega*sh,0,0,1;
    B_upd<<timeStep-sh/omega,1-ch,timeStep;

    Eigen::Vector3d currentState = Eigen::Vector3d(comPos(dim),comVel(dim),zmpPos(dim));

    if (dim == 2) return A_upd*(currentState + Eigen::Vector3d(0.0,0.0,comTargetHeight)) + B_upd*zmpDot - Eigen::Vector3d(0.0,0.0,comTargetHeight);
    
    return A_upd*currentState + B_upd*zmpDot;
  }

 private:
  double clamp(double N_, double n_min, double n_max) {return std::max(n_min,std::min(n_max,N_));}

  // Constant parameters
  int S,D;
  unsigned int footstepCounter = 0;
  
  double mpcTimeStep;
  double controlTimeStep;
  double singleSupportDuration, doubleSupportDuration;
  double comTargetHeight;
  double footContraintSquareWidth;
  double omega;

  // Parameters for the current iteration
  bool supportFoot;
  double simulationTime;
  int mpcIter = 0, controlIter = 0;

  // Matrices for prediction
  Eigen::Matrix<double, numVariables_ / 3, 1> p;
  Eigen::Matrix<double, numVariables_ / 3, numVariables_ / 3> P;

  // Matrices for cost function
  Eigen::Matrix<double, numVariables_, numVariables_> costFunctionH;
  Eigen::Matrix<double, numVariables_, 1> costFunctionF;
  double beta_ = 1000.0;

  // Matrices for stability constraint
  Eigen::Matrix<double, numEqualityConstraints_, numVariables_> Aeq;
  Eigen::Matrix<double, numEqualityConstraints_, 1> beq;

  //Matrices for balance constraint
  Eigen::Matrix<double, numInequalityConstraints_, numVariables_> AZmp;
  Eigen::Matrix<double, numInequalityConstraints_, 1> bZmpMax;
  Eigen::Matrix<double, numInequalityConstraints_, 1> bZmpMin;

  // State
  Eigen::Vector3d comPos;
  Eigen::Vector3d comVel = Eigen::Vector3d::Zero();
  Eigen::Vector3d zmpPos;

  labrob::qpsolvers::QPSolverEigenWrapper<double, numVariables_, numEqualityConstraints_, numInequalityConstraints_> qp_solver_;

};

}

