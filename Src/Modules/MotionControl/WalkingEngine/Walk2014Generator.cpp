/*
 * Copyright 2010 The University of New South Wales (UNSW).
 *
 * This file is part of the 2010 team rUNSWift RoboCup entry. You may
 * redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version as
 * modified below. As the original licensors, we add the following
 * conditions to that license:
 *
 * In paragraph 2.b), the phrase "distribute or publish" should be
 * interpreted to include entry into a competition, and hence the source
 * of any derived work entered into a competition must be made available
 * to all parties involved in that competition under the terms of this
 * license.
 *
 * In addition, if the authors of a derived work publish any conference
 * proceedings, journal articles or other academic papers describing that
 * derived work, then appropriate academic citations to the original work
 * must be included in that publication.
 *
 * This rUNSWift source is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this source code; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

/**
 * @file Walk2014Generator.cpp
 *
 * This file implements the UNSW 2014 walk generator. It was refactored to fit into the B-Human framework.
 *
 * The period of each foot-step is set by T. T generates the forcing function by alternatively lifting each foot.
 * Control:
 * The change of support foot is driven by the ZMP switching sign, but must be > T/2. If > 3*T we switch to try
 * to revive.
 * The front-back sway is controlled by ankle tilts proportional to the y gyro.
 * The user specifies forward(m), left(m), turn(radians) to activate the walk.
 * If these values are all zero the robot stands.
 *
 * @author Bernhard Hengst
 * @date 18 Jan 2014
 *
 * @author Thomas Röfer
 * @date 2 Oct 2018
 */

#include "Walk2014Generator.h"
#include "Platform/BHAssert.h"
#include "Platform/SystemCall.h"
#include "Representations/Sensing/RobotModel.h"
#include "Tools/Debugging/Annotation.h"
#include "Tools/Math/Constants.h"
#include "Tools/Math/Rotation.h"
#include "Tools/Motion/InverseKinematic.h"
#include "Tools/Range.h"
#include <algorithm>

#include <chrono>
#include <iostream>
#include <fstream>

MAKE_MODULE(Walk2014Generator, motionControl);

static const double mmPerM = 1000.0;

std::ofstream com_file("/tmp/com.txt");
std::ofstream zmp_file("/tmp/zmp.txt");
std::ofstream lsole_file("/tmp/lsole.txt");
std::ofstream rsole_file("/tmp/rsole.txt");
std::ofstream supp_file("/tmp/supp.txt");

template <class T>
Eigen::Matrix<T, 2, 2>
Rz_planar(T theta) {
  T c = std::cos(theta);
  T s = std::sin(theta);
  Eigen::Matrix<T, 2, 2> R;
  R << c, -s,
       s,  c;
  return R;
}

Eigen::Matrix3d Rz(double theta) {
  Eigen::Matrix3d R;
  double ctheta = std::cos(theta);
  double stheta = std::sin(theta);
  R << ctheta, -stheta, 0.0,
       stheta,  ctheta, 0.0,
          0.0,     0.0, 1.0;
  return R;
}

Eigen::Vector4d T_inv(const Eigen::Vector4d& T) {
  Eigen::Vector3d p = T.head<3>();
  Eigen::Matrix3d RzT = Rz(-T.w());
  Eigen::Vector4d T_r;
  T_r << -RzT * p, -T.w();
  return T_r;
}

Eigen::Vector4d T_mul(const Eigen::Vector4d& T1, const Eigen::Vector4d& T2) {
  Eigen::Matrix3d R1 = Rz(T1.w());
  Eigen::Vector3d p1 = T1.head<3>();
  Eigen::Vector3d p2 = T2.head<3>();
  Eigen::Vector4d T_r;
  T_r << R1 * p2 + p1, T1.w() + T2.w();
  return T_r;
}

template <class T>
T
angle_difference(T alpha, T beta) {
  Eigen::Matrix<T, 2, 2> R_diff = Rz_planar<T>(alpha - beta);
  return std::atan2(R_diff(1, 0), R_diff(0, 0));
}

Walk2014Generator::Walk2014Generator() {
  double z_torso = 0.24809; // taken from http://doc.aldebaran.com/2-1/family/robots/links_robot.html
  Pose T_torso_w = Pose(
      Eigen::Vector3d(0.0, 0.0, z_torso),
      Eigen::Matrix3d::Identity()
  );
  Eigen::Vector3d p_lsole_torso = Eigen::Vector3d(0.0, 0.05, -z_torso);
  Eigen::Matrix3d R_lsole_torso = Eigen::Matrix3d::Identity();
  Eigen::Vector3d p_rsole_torso = Eigen::Vector3d(0.0, -0.05, -z_torso);
  Eigen::Matrix3d R_rsole_torso = Eigen::Matrix3d::Identity();
  Pose T_lsole_torso(p_lsole_torso, R_lsole_torso);
  Pose T_rsole_torso(p_rsole_torso, R_rsole_torso);
  Pose T_lsole_w = T_torso_w * T_lsole_torso;
  Pose T_rsole_w = T_torso_w * T_rsole_torso;

  Eigen::Vector4d qL_init, qR_init;
  qL_init << T_lsole_w.position, T_lsole_w.orientation.eulerAngles(2, 1, 0).x();
  qR_init << T_rsole_w.position, T_rsole_w.orientation.eulerAngles(2, 1, 0).x();
  starting_configuration_ = Configuration(
    qL_init,
    qR_init,
    Foot::RIGHT,
    0.0
  );

  // Setup MPC solver:
  const Eigen::Vector3d& p_lsole_w = T_lsole_w.position;
  const Eigen::Vector3d& p_rsole_w = T_rsole_w.position;
  Eigen::Vector3d p_zmp_w = (p_lsole_w + p_rsole_w) / 2.0;
  Eigen::Vector3d p_com_w =
      p_zmp_w + Eigen::Vector3d(0.0, 0.0, com_target_height_);
  mpc_solver_ptr_ = std::make_shared<mpcSolver::MPCSolver<
      numVariables_, numEqualityConstraints_, numInequalityConstraints_>>(
    mpc_timestep_,
    controller_timestep_,
    single_support_duration_,
    double_support_duration_,
    com_target_height_,
    mpc_foot_constraint_size_,
    p_com_w,
    p_zmp_w
  );

  std::cerr << "Connecting to server..." << std::endl;
  if (tcp_client_.connectToServer(ip_addr_.c_str(), port_)) {
    tcp_client_.subscribeToFootstepPlan(&Walk2014Generator::footstepPlanCallback, this);
  } else {
    std::cerr << "Cannot connect to server." << std::endl;
  }
  std::cerr << "Done." << std::endl;
}

void
Walk2014Generator::footstepPlanCallback(const FootstepPlan& footstep_plan) {
  const std::lock_guard<std::mutex> lock(footstepPlanMutex_);
  footstep_plan_ = footstep_plan;

  std::cerr << "Footstep plan received." << std::endl;
  for (const auto& configuration : footstep_plan_) {
    std::cerr << configuration.to_string() << std::endl;
  }
}

void Walk2014Generator::update(WalkGenerator& generator)
{
  // Use other parameters in demo games to take care of the robots
  if(theGlobalOptions.slowWalk)
  {
    maxSpeed.translation.x() = std::min(maxSpeed.translation.x(), slowMaxSpeed.translation.x());
    maxSpeed.translation.y() = std::min(maxSpeed.translation.y(), slowMaxSpeed.translation.y());
    maxSpeed.rotation = std::min(maxSpeed.rotation, slowMaxSpeed.rotation);
    maxSpeedBackwards = std::min(maxSpeedBackwards, slowMaxSpeedBackwards);
    maxAcceleration.x() = std::min(maxAcceleration.x(), slowMaxAcceleration.x());
    maxAcceleration.y() = std::min(maxAcceleration.y(), slowMaxAcceleration.y());
  }

  generator.reset = [this, &generator]() { reset(generator); };
  generator.calcJoints =
    [this, &generator](const Pose2f& speed, const Pose2f& target, WalkGenerator::WalkMode walkMode,
                       const std::function<Pose3f(float phase)>& getKickFootOffset)
  {
    calcJoints(generator, speed, target, walkMode, getKickFootOffset);
  };
  generator.maxSpeed = Pose2f(maxSpeed.rotation * odometryScale.rotation,
                              maxSpeed.translation.x() * odometryScale.translation.x(),
                              maxSpeed.translation.y() * odometryScale.translation.y());

  filteredGyroX = gyroLowPassRatio * filteredGyroX + (1.f - gyroLowPassRatio) * theInertialSensorData.gyro.x();
  filteredGyroY = gyroLowPassRatio * filteredGyroY + (1.f - gyroLowPassRatio) * theInertialSensorData.gyro.y();
}

void Walk2014Generator::reset(WalkGenerator& generator)
{
  generator.stepDuration = 0.f;
  generator.t = 0.f;
  generator.speed = Pose2f();
  generator.upcomingOdometryOffset = Pose2f();
  walkState = standing;
  timeWhenStandBegan = theFrameInfo.time;
  if(theRobotInfo.penalty != PENALTY_NONE && theRobotInfo.penalty != PENALTY_SPL_ILLEGAL_MOTION_IN_SET)
    filteredGyroX = filteredGyroY = 0;
  forward = lastForward = 0.f;
  forwardL = forwardL0 = 0.f;
  forwardR = forwardR0 = 0.f;
  left = lastLeft = 0.f;
  leftL = leftR = 0_deg;
  turnRL = turnRL0 = 0_deg;
  swingAngle = 0_deg;
  switchPhase = 0.f;
  maxFootHeight = maxFootHeight0 = 0.f;
  weightShiftStatus = weightDidNotShift;
  prevForwardL = prevForwardR = 0.f;
  prevLeftL = prevLeftR = 0_deg;
  prevTurn = 0_deg;
  weightShiftMisses = 0;
  slowWeightShifts = 0;
  torsoTilt = 0;
}

void Walk2014Generator::calcJoints(WalkGenerator& generator,
                                   const Pose2f& speed,
                                   const Pose2f& target,
                                   WalkGenerator::WalkMode walkMode,
                                   const std::function<Pose3f(float phase)>& getKickFootOffset)
{
  const std::lock_guard<std::mutex> lock(footstepPlanMutex_);

  std::string walking_state_str;
  if (walking_state_ == WalkingState::Standing) {
    walking_state_str = "Standing";
  } else if (walking_state_ == WalkingState::Starting) {
    walking_state_str = "Starting";
  } else if (walking_state_ == WalkingState::Walking) {
    walking_state_str = "Walking";
  } else {
    walking_state_str = "Stopping";
  }

  // Update walking data:
  if (mpc_iter_ == 0 && control_iter_ == 0) {
    // Retrieve data from footstep plan:
    if (walking_state_ == WalkingState::Walking && footstep_plan_.size() >= 2) {
      target_configuration_ = footstep_plan_[1];
    } else {
      target_configuration_ = starting_configuration_;
      target_configuration_.setSwingFootTrajectoryHeight(0.0);
    }

    // Send target configuration to footstep planner:
    if (tcp_client_.sendConfiguration(target_configuration_)) {
      std::cerr << "Sending: " << target_configuration_.to_string() << std::endl;
    } else {
      std::cerr << "Cannot send configuration to footstep planner.";
    }

    // Setup swing foot trajectory:
    swing_foot_trajectory_ = 
        [&](double s) {
          Eigen::Vector4d T0, Tf;
          T0 = starting_configuration_.getSwingFootConfiguration();
          if (walking_state_ == WalkingState::Standing ||
              walking_state_ == WalkingState::Starting) {
            Tf = starting_configuration_.getSwingFootConfiguration();
          } else {
            Tf = target_configuration_.getSupportFootConfiguration();
          }
          double h_z = target_configuration_.h_z_;
          double s_h = 0.5;
          double z0 = T0.z();
          double zf = Tf.z();
          double a = (h_z - z0 + s_h * (z0 - zf)) / (s_h * (s_h - 1.0));
          double b = zf - z0 - a;
          double c = z0;
          double s_0 = 0.2;
          double s_f = 0.8;
          double s_bar = (s - s_0) / (s_f - s_0);
          double swing_x = T0.x() + (Tf.x() - T0.x()) * s_bar;
          double swing_y = T0.y() + (Tf.y() - T0.y()) * s_bar;
          double swing_z = a * s * s + b * s + c;
          double swing_theta = T0.w() + angle_difference(Tf.w(), T0.w()) * s_bar;
          if (s_0 <= s && s <= s_f) {
            return Pose(
                Eigen::Vector3d(swing_x, swing_y, swing_z),
                Rz(swing_theta)
            );
          } else if (s < s_0) {
            return Pose(
                Eigen::Vector3d(T0.x(), T0.y(), swing_z),
                Rz(T0.w())
            );
          } else {
            return Pose(
                Eigen::Vector3d(Tf.x(), Tf.y(), swing_z),
                Rz(Tf.w())
            );
          }
        };

    const auto& qSupport = starting_configuration_.getSupportFootConfiguration();
    const auto& qSwing = starting_configuration_.getSwingFootConfiguration();
    const auto& qTarget = target_configuration_.getSupportFootConfiguration();
    // Update MPC params:
    Eigen::Vector4d qMiddle = (qSupport + qSwing) / 2.0;
    if (walking_state_ == WalkingState::Standing) {
      mpc_plan_.clear();
      mpc_plan_.push_back(qMiddle);
      mpc_plan_.push_back(qMiddle);
      mpc_plan_.push_back(qMiddle);
    } else if (walking_state_ == WalkingState::Starting) {
      mpc_plan_.clear();
      if (footstep_plan_.size() >= 2) {
        mpc_plan_.push_back(qMiddle);
        mpc_plan_.push_back(footstep_plan_[0].getSupportFootConfiguration());
        mpc_plan_.push_back(footstep_plan_[1].getSupportFootConfiguration());
      } else if (footstep_plan_.size() >= 1) {
        mpc_plan_.push_back(qMiddle);
        mpc_plan_.push_back(footstep_plan_[0].getSupportFootConfiguration());
        mpc_plan_.push_back(footstep_plan_[0].getSupportFootConfiguration());
      } else {
        mpc_plan_.push_back(qMiddle);
        mpc_plan_.push_back(qMiddle);
        mpc_plan_.push_back(qMiddle);
      }
    } else if (walking_state_ == WalkingState::Walking) {
      mpc_plan_.clear();
      if (footstep_plan_.size() >= 3) {
        mpc_plan_.push_back(footstep_plan_[0].getSupportFootConfiguration());
        mpc_plan_.push_back(footstep_plan_[1].getSupportFootConfiguration());
        mpc_plan_.push_back(footstep_plan_[2].getSupportFootConfiguration());
      } else if (footstep_plan_.size() >= 2) {
        const auto& qSupportFinal = footstep_plan_[1].getSupportFootConfiguration();
        const auto& qSwingFinal   = footstep_plan_[1].getSwingFootConfiguration();
        Eigen::Vector4d qMiddleFinal = (qSupportFinal + qSwingFinal) / 2.0;
        mpc_plan_.push_back(footstep_plan_[0].getSupportFootConfiguration());
        mpc_plan_.push_back(qMiddleFinal);
        mpc_plan_.push_back(qMiddleFinal);
      }
    }
  }
      
  //std::cerr << "WalkingState::" << walking_state_str << std::endl;
  //std::cerr << "\tcontrol iter=" << control_iter_ << std::endl;
  //std::cerr << "\tmpc iter=" << mpc_iter_ << std::endl;
  //std::cerr << "\tstarting configuration: " << starting_configuration_.to_string() << std::endl;
  //std::cerr << "\ttarget configuration: " << target_configuration_.to_string() << std::endl;

  std::chrono::time_point<std::chrono::system_clock> t0 = std::chrono::system_clock::now();
  mpc_solver_ptr_->solve(mpc_plan_);
  std::chrono::time_point<std::chrono::system_clock> tf = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = tf - t0;
  //std::cerr << "MPC solved in " << elapsed_seconds.count() << std::endl;

  // Compute pose of support and swing foot:
  const auto& p_com_w_desired = mpc_solver_ptr_->getOptimalCoMPosition();
  const auto& p_zmp_w_desired = mpc_solver_ptr_->getOptimalZMPPosition();
  Pose T_supp_w_t0(
      starting_configuration_.getSupportFootConfiguration().head<3>(),
      Rz(starting_configuration_.getSupportFootConfiguration().w())
  );
  Eigen::Vector3d p_com_supp_desired = T_supp_w_t0.inv() * p_com_w_desired;
  Pose T_torso_supp_desired(p_com_supp_desired, Eigen::Matrix3d::Identity());

  double t = controller_timestep_ * control_iter_ + mpc_timestep_ * mpc_iter_;
  double s = 1.0;
  if (t < single_support_duration_) s = t / single_support_duration_;
  auto T_swing_w_desired = swing_foot_trajectory_(s);
  auto T_swing_supp_desired = T_supp_w_t0.inv() * T_swing_w_desired;

  Pose T_left_torso_desired;
  Pose T_right_torso_desired;

  if (starting_configuration_.getSupportFoot() == Foot::LEFT) {
    T_left_torso_desired = T_torso_supp_desired.inv();
    T_right_torso_desired = T_left_torso_desired * T_swing_supp_desired;
  } else {
    T_right_torso_desired = T_torso_supp_desired.inv();
    T_left_torso_desired  = T_right_torso_desired * T_swing_supp_desired;
  }

  // Setup data for IK:
  Eigen::Matrix3f R_left_torso = T_left_torso_desired.orientation.cast<float>();
  Eigen::Vector3f p_left_torso = T_left_torso_desired.position.cast<float>() * mmPerM;
  Eigen::Matrix3f R_right_torso = T_right_torso_desired.orientation.cast<float>();
  Eigen::Vector3f p_right_torso = T_right_torso_desired.position.cast<float>() * mmPerM;
  Pose3f leftFoot  = Pose3f(R_left_torso, p_left_torso);
  Pose3f rightFoot = Pose3f(R_right_torso, p_right_torso);

  // Log data:
  com_file << p_com_w_desired.transpose() << std::endl;
  zmp_file << p_zmp_w_desired.transpose() << std::endl;
  lsole_file << T_left_torso_desired.position.transpose() << std::endl;
  rsole_file << T_right_torso_desired.position.transpose() << std::endl;
  if (control_iter_ == 0 && mpc_iter_ == 0) supp_file << T_supp_w_t0.position.transpose() << " " << T_supp_w_t0.orientation.eulerAngles(2, 1, 0).x() << std::endl;

  // Update state of iters:
  control_iter_ = (control_iter_ + 1) %
      ((int) std::round(mpc_timestep_ / controller_timestep_));
  if (delay_) --delay_;
  if (control_iter_ == 0) {
    mpc_iter_ = (mpc_iter_ + 1) % (S_ + D_);
    // Update walking state:
    if (mpc_iter_ == 0 && walking_state_ == WalkingState::Standing &&
        !footstep_plan_.empty() && !delay_) {
      walking_state_ = WalkingState::Starting;
      starting_configuration_ = footstep_plan_.front();
    } else if (mpc_iter_ == 0 && walking_state_ == WalkingState::Starting) {
      walking_state_ = WalkingState::Walking;
    } else if (mpc_iter_ == 0 && walking_state_ == WalkingState::Walking) {
      footstep_plan_.pop_front();
      starting_configuration_ = footstep_plan_.front();
    } else if (mpc_iter_ == S_ && walking_state_ == WalkingState::Walking) {
      if (footstep_plan_.size() <= 2) {
        walking_state_ = WalkingState::Stopping;
      }
    } else if (mpc_iter_ == 0 && walking_state_ == WalkingState::Stopping) {
      footstep_plan_.pop_front();
      starting_configuration_ = footstep_plan_.front();
      walking_state_ = WalkingState::Standing;
      footstep_plan_.clear();
    }
  }

  /****************************************************************************/

  // 9.3 Inverse kinematics
  VERIFY(InverseKinematic::calcLegJoints(leftFoot, rightFoot, Vector2f::Zero(), generator.jointRequest, theRobotDimensions) || SystemCall::getMode() == SystemCall::logfileReplay);

  // 10. Set joint values and stiffness
  int stiffness = walking_state_ == WalkingState::Standing ? StiffnessData::useDefault : walkStiffness;
  for(uint8_t i = Joints::firstLegJoint; i < Joints::numOfJoints; ++i)
    generator.jointRequest.stiffnessData.stiffnesses[i] = stiffness;

  // Arm Swing
  generator.jointRequest.angles[Joints::lShoulderRoll] =  30_deg/2.0f;
  generator.jointRequest.angles[Joints::rShoulderRoll] = -30_deg/2.0f;
  generator.jointRequest.angles[Joints::lShoulderPitch] = 90_deg + 0.5f * generator.jointRequest.angles[Joints::rHipPitch];
  generator.jointRequest.angles[Joints::rShoulderPitch] = 90_deg + 0.5f * generator.jointRequest.angles[Joints::lHipPitch];
  generator.jointRequest.angles[Joints::lElbowRoll] =  0.5f*generator.jointRequest.angles[Joints::rHipPitch];
  generator.jointRequest.angles[Joints::rElbowRoll] = -0.5f*generator.jointRequest.angles[Joints::lHipPitch];
  generator.jointRequest.angles[Joints::lElbowYaw] = -90_deg;
  generator.jointRequest.angles[Joints::rElbowYaw] =  90_deg;
  generator.jointRequest.angles[Joints::lWristYaw] =  0_deg;
  generator.jointRequest.angles[Joints::rWristYaw] = -0_deg;

  // 7. Sagittal balance
  filteredGyroY = gyroLowPassRatio * filteredGyroY + (1.f - gyroLowPassRatio) * theInertialSensorData.gyro.y();
  //Angle balanceAdjustment = walkState == standing ? 0.f : filteredGyroY *
  //                          (filteredGyroY > 0 ? gyroForwardBalanceFactor : gyroBackwardBalanceFactor); // adjust ankle tilt in proportion to filtered gryoY
  //generator.jointRequest.angles[generator.isLeftPhase ? Joints::rAnklePitch : Joints::lAnklePitch] += balanceAdjustment;
  Angle balanceAdjustment = filteredGyroY * 0.05f;
  generator.jointRequest.angles[Joints::lAnklePitch] += balanceAdjustment;
  generator.jointRequest.angles[Joints::rAnklePitch] += balanceAdjustment;

  // Lateral balance
  filteredGyroX = gyroLowPassRatio * filteredGyroX + (1.f - gyroLowPassRatio) * theInertialSensorData.gyro.x();
  //balanceAdjustment = filteredGyroX * gyroSidewaysBalanceFactor;
  balanceAdjustment = filteredGyroX * 0.1f;
  generator.jointRequest.angles[Joints::lAnkleRoll] += balanceAdjustment;
  generator.jointRequest.angles[Joints::rAnkleRoll] += balanceAdjustment;

  // Head can move freely
  generator.jointRequest.angles[Joints::headPitch] = generator.jointRequest.angles[Joints::headYaw] = JointAngles::ignore;
}

