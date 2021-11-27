#pragma once

// STL
#include <string>

// Eigen
#include <Eigen/Core>

#include "Foot.hpp"

class Configuration {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  Configuration() = default;

  Configuration(
      const Eigen::Vector4d& qL,
      const Eigen::Vector4d& qR,
      Foot support_foot,
      double h_z) :
      qL_(qL),
      qR_(qR),
      support_foot_(support_foot),
      h_z_(h_z) {

  }

  bool isApprox(const Configuration& rhs, double precision=0.001) const {
    if (support_foot_ != rhs.support_foot_) {
      return false;
    }

    return qL_.isApprox(rhs.qL_, precision) &&
           qR_.isApprox(rhs.qR_, precision);
  }

  Foot getSupportFoot() const {
    return support_foot_;
  }

  const Eigen::Vector4d& getSupportFootConfiguration() const {
    if (support_foot_ == Foot::LEFT) {
      return qL_;
    } else {
      return qR_;
    }
  }

  const Eigen::Vector4d& getSwingFootConfiguration() const {
    if (support_foot_ == Foot::LEFT) {
      return qR_;
    } else {
      return qL_;
    }
  }

  void setSupportFoot(Foot support_foot) {
    support_foot_ = support_foot;
  }

  void setSwingFootTrajectoryHeight(double h_z) {
    h_z_ = h_z;
  }

  std::string to_string() const {
    std::string support_foot_str = "LEFT";
    if (support_foot_ == Foot::RIGHT) {
      support_foot_str = "RIGHT";
    }
    return std::string("<(") +
        std::to_string(qL_.x()) + ", " +
        std::to_string(qL_.y()) + ", " +
        std::to_string(qL_.z()) + ", " + 
        std::to_string(qL_.w()) + "), " +
        "(" +
        std::to_string(qR_.x()) + ", " +
        std::to_string(qR_.y()) + ", " +
        std::to_string(qR_.z()) + ", " +
        std::to_string(qR_.w()) + "), " +
        std::to_string(h_z_) + ", " +
        support_foot_str + ">";
  }

  Eigen::Vector4d qL_;
  Eigen::Vector4d qR_;
  Foot support_foot_;
  double h_z_;
};