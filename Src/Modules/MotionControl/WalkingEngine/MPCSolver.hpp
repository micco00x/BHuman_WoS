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
      double prediction_time,
      double single_support_duration,
      double double_support_duration,
      double com_target_height,
      double foot_constraint_square_width,
      const Eigen::Vector3d& com_position,
      const Eigen::Vector3d& zmp_position) :
      mpcTimeStep(mpc_timestep),
      controlTimeStep(control_timestep),
      predictionTime(prediction_time),
      singleSupportDuration(single_support_duration),
      doubleSupportDuration(double_support_duration),
      comTargetHeight(com_target_height),
      footContraintSquareWidth(foot_constraint_square_width),
      omega(std::sqrt(9.81 / comTargetHeight)),
      comPos(com_position),
      zmpPos(zmp_position),
      qp_solver_(std::make_shared<labrob::qpsolvers::QPOASESQPSolver<numVariables_, numEqualityConstraints_, numInequalityConstraints_>>()) {
    generatePredictionMatrices();
  }

  // Main method
  void solve(const std::vector<Eigen::VectorXd>& plan){

    // Save iteration parameters
    this->plan = plan;

    // Mapping
    Eigen::MatrixXd mapping(N,plan.size());
    for (int j = 0; j < static_cast<int>(plan.size()); j++) {
      for (int i = 0; i < N; i++) {
        mapping(i,j) = clamp(- (double)i/D + 1 + (double)(S-mpcIter+j*(S+D))/D, 0, 1) - clamp(- (double)i/D + 1 + (double)(S-mpcIter+(j-1)*(S+D))/D, 0, 1);
      }
    }
    for (int i = 1; i < N; i++) {
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

    Eigen::MatrixXd rotSin = Eigen::MatrixXd::Zero(N,N);
    Eigen::MatrixXd rotCos = Eigen::MatrixXd::Zero(N,N);
    for (int i = 0; i < N; i++) {
      rotSin(i,i) = sin(mc_th(i));
      rotCos(i,i) = cos(mc_th(i));
    }
    Eigen::MatrixXd zmpRotMatrix(3*N,3*N);
    //zmpRotMatrix <<  rotCos, rotSin, Eigen::MatrixXd::Zero(N,N),
    //                -rotSin, rotCos, Eigen::MatrixXd::Zero(N,N),
    //                 Eigen::MatrixXd::Zero(N,2*N), Eigen::MatrixXd::Identity(N,N);
    zmpRotMatrix = Eigen::MatrixXd::Identity(3*N,3*N);

    // Matrices for ZMP prediction
    Eigen::MatrixXd P = Eigen::MatrixXd::Ones(N,N)*mpcTimeStep;

    for(int i=0; i<N;++i){
      for(int j=0;j<N;++j){
        if (j>i) P(i,j)=0;
      }
    }

    AZmp.setZero();
    AZmp.block(0,0,N,N) = P;
    AZmp.block(N,N,N,N) = P;
    AZmp.block(2*N,2*N,N,N) = P;
    AZmp = zmpRotMatrix * AZmp;

    Eigen::VectorXd bZmpMinFixed(3*N);
    Eigen::VectorXd bZmpMaxFixed(3*N);
    Eigen::VectorXd bZmpMinState(3*N);
    Eigen::VectorXd bZmpMaxState(3*N);

    bZmpMinFixed << mc_x - Eigen::VectorXd::Ones(N)*footContraintSquareWidth/2, 
                    mc_y - Eigen::VectorXd::Ones(N)*footContraintSquareWidth/2, 
                    mc_z - Eigen::VectorXd::Ones(N)*footContraintSquareWidth/2;

    bZmpMaxFixed << mc_x + Eigen::VectorXd::Ones(N)*footContraintSquareWidth/2, 
                    mc_y + Eigen::VectorXd::Ones(N)*footContraintSquareWidth/2, 
                    mc_z + Eigen::VectorXd::Ones(N)*footContraintSquareWidth/2;

    bZmpMinState << Eigen::VectorXd::Ones(N) * zmpPos(0), 
                    Eigen::VectorXd::Ones(N) * zmpPos(1), 
                    Eigen::VectorXd::Ones(N) * zmpPos(2);

    bZmpMaxState << Eigen::VectorXd::Ones(N) * zmpPos(0), 
                    Eigen::VectorXd::Ones(N) * zmpPos(1), 
                    Eigen::VectorXd::Ones(N) * zmpPos(2);

    bZmpMinState = zmpRotMatrix * bZmpMinState;
    bZmpMaxState = zmpRotMatrix * bZmpMaxState;

    bZmpMin = bZmpMinFixed - bZmpMinState;
    bZmpMax = bZmpMaxFixed - bZmpMaxState;

    Eigen::VectorXd b(N);
    Aeq = Eigen::MatrixXd::Zero(3,3*N);

    for(int i=0;i<N;++i){
      b(i) = pow(exp(-omega*mpcTimeStep),i);
    }

    Aeq.block(0,0,1,N)       = (1/omega)*(1-exp(-omega*mpcTimeStep))*b.transpose();
    Aeq.block(1,N,1,N)     = (1/omega)*(1-exp(-omega*mpcTimeStep))*b.transpose();
    Aeq.block(2,2*N,1,N) = (1/omega)*(1-exp(-omega*mpcTimeStep))*b.transpose();

    beq = Eigen::VectorXd(3);
    beq << comPos(0) + comVel(0)/omega - zmpPos(0),
          comPos(1) + comVel(1)/omega - zmpPos(1),
          comPos(2) + comVel(2)/omega - (zmpPos(2) + comTargetHeight);

    costFunctionH = Eigen::MatrixXd::Zero(3*N,3*N);
    costFunctionF = Eigen::VectorXd::Zero(3*N);
    costFunctionH.block(0,0,N,N) = Eigen::MatrixXd::Identity(N,N);
    costFunctionH.block(N,N,N,N) = Eigen::MatrixXd::Identity(N,N);
    costFunctionH.block(2*N,2*N,N,N) = Eigen::MatrixXd::Identity(N,N);

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
    Eigen::VectorXd decisionVariables = qp_solver_.get_solution();

    // Split the QP solution in ZMP dot and footsteps
    Eigen::VectorXd zDotOptimalX(N);
    Eigen::VectorXd zDotOptimalY(N);
    Eigen::VectorXd zDotOptimalZ(N);

    zDotOptimalX = (decisionVariables.head(N));
    zDotOptimalY = (decisionVariables.segment(N,N));
    zDotOptimalZ = (decisionVariables.segment(2*N,N));

    // Update the state based on the result of the QP
    Eigen::Vector3d nextStateX = updateState(zDotOptimalX(0),0,controlTimeStep);
    Eigen::Vector3d nextStateY = updateState(zDotOptimalY(0),1,controlTimeStep);
    Eigen::Vector3d nextStateZ = updateState(zDotOptimalZ(0),2,controlTimeStep);

    comPos << nextStateX(0),nextStateY(0),nextStateZ(0);
    comVel << nextStateX(1),nextStateY(1),nextStateZ(1);
    zmpPos << nextStateX(2),nextStateY(2),nextStateZ(2);
    optimalFootsteps = plan.at(1);
    //optimalFootsteps << plan_x(1),plan_y(1),plan_z(1),plan_th(1);

    ++controlIter;
    mpcIter = floor(controlIter*controlTimeStep/mpcTimeStep);
    if(mpcIter>=S+D){
      controlIter = 0;
      mpcIter = 0;
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

  // Variable footstep timing
  void generatePredictionMatrices() {
    N = round(predictionTime/mpcTimeStep);
    S = round(singleSupportDuration/mpcTimeStep);
    D = round(doubleSupportDuration/mpcTimeStep);
    M = ceil(N/(S+D));
  
    // Matrices for feet rotations
    predictedRotations = Eigen::VectorXd::Zero(M+1);

    // Matrices for ZMP prediction
    p =  Eigen::VectorXd::Ones(N);
    P =  Eigen::MatrixXd::Ones(N,N)*mpcTimeStep;

    for(int i=0; i<N;++i){
      for(int j=0;j<N;++j){
        if (j>i) P(i,j)=0;
      }
    }
  }

  // Update the state
  Eigen::Vector3d updateState(double zmpDot, int dim, double timeStep) {
    // Update the state along the dim-th direction (0,1,2) = (x,y,z)

    double ch = cosh(omega*timeStep);
    double sh = sinh(omega*timeStep);

    Eigen::Matrix3d A_upd = Eigen::MatrixXd::Zero(3,3);
    Eigen::Vector3d B_upd = Eigen::VectorXd::Zero(3);
    A_upd<<ch,sh/omega,1-ch,omega*sh,ch,-omega*sh,0,0,1;
    B_upd<<timeStep-sh/omega,1-ch,timeStep;

    Eigen::Vector3d currentState = Eigen::Vector3d(comPos(dim),comVel(dim),zmpPos(dim));

    if (dim == 2) return A_upd*(currentState + Eigen::Vector3d(0.0,0.0,comTargetHeight)) + B_upd*zmpDot - Eigen::Vector3d(0.0,0.0,comTargetHeight);
    
    return A_upd*currentState + B_upd*zmpDot;
  }

 private:
  double clamp(double n, double n_min, double n_max) {return std::max(n_min,std::min(n_max,n));}

  // Constant parameters
  int N,S,D,M;
  unsigned int footstepCounter = 0;
  
  double mpcTimeStep;
  double controlTimeStep;
  double predictionTime;
  double singleSupportDuration, doubleSupportDuration;
  double comTargetHeight;
  double footContraintSquareWidth;
  double omega;

  // Parameters for the current iteration
  bool supportFoot;
  double simulationTime;
  int mpcIter = 0, controlIter = 0;

  // Matrices for prediction
  Eigen::VectorXd p;
  Eigen::MatrixXd P;

  // Matrices for cost function
  Eigen::Matrix<double, numVariables_, numVariables_> costFunctionH;
  Eigen::Matrix<double, numVariables_, 1> costFunctionF;

  // Matrices for stability constraint
  Eigen::Matrix<double, numEqualityConstraints_, numVariables_> Aeq;
  Eigen::Matrix<double, numEqualityConstraints_, 1> beq;

  //Matrices for balance constraint
  Eigen::Matrix<double, numInequalityConstraints_, numVariables_> AZmp;
  Eigen::Matrix<double, numInequalityConstraints_, 1> bZmpMax;
  Eigen::Matrix<double, numInequalityConstraints_, 1> bZmpMin;

  // Solution of the QP for determining orientations
  Eigen::VectorXd predictedRotations;

  // State
  Eigen::Vector3d comPos;
  Eigen::Vector3d comVel = Eigen::Vector3d::Zero();
  Eigen::Vector3d zmpPos;
  Eigen::Vector4d optimalFootsteps = Eigen::Vector4d::Zero();

  // Support foot
  Eigen::Vector4d currentSupportFoot = Eigen::Vector4d::Zero();
  Eigen::Vector4d nextSupportFoot;

  std::vector<Eigen::VectorXd> plan;

  labrob::qpsolvers::QPSolverEigenWrapper<double, numVariables_, numEqualityConstraints_, numInequalityConstraints_> qp_solver_;

};

}

