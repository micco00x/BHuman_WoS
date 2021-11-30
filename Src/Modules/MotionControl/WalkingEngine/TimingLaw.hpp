#pragma once

#include <cmath>

namespace labrob {

class TimingLaw {
 public:
  virtual double eval(double t) const = 0;
  virtual double eval_dt(double t) const = 0;
}; // end class TimingLaw

class LinearTimingLaw : public TimingLaw {
 public:
  LinearTimingLaw(double duration) : duration_(duration) {

  }

  double eval(double t) const override {
    return t / duration_;
  }

  double eval_dt(double t) const override {
    return 1.0 / duration_;
  }

 private:
  double duration_;
}; // end class LinearTimingLaw

/**
 * \brief Class for timing law with trapezoidal acceleration.
 */
class TrapezoidalAccelerationTimingLaw : public TimingLaw {
 public:
  /**
   * \brief Constructor of a timing law with trapezoidal acceleration.
   * \param duration Duration of the timing law.
   * \param alpha Define where to have maximum acceleration. 
   */
  TrapezoidalAccelerationTimingLaw(double duration, double alpha)
  : duration_(duration), alpha_(alpha) {
    t0_ = 0.0;
    ta_ = t0_ + alpha_ * duration_;
    tb_ = t0_ + (1.0 - alpha_) * duration_;
    tf_ = t0_ + duration_;

    a_max_ = std::pow(
        (std::pow(ta_, 2.0) / 2.0 - t0_ * ta_ + std::pow(t0_, 2.0) / 2.0) + 
        (ta_ - t0_) * (tf_ - ta_) -
        (std::pow(tf_, 2.0) / 2.0 - tb_ * tf_ + std::pow(tb_, 2.0) / 2.0), -1.0);
  }

  /**
   * \brief Evaluate the timing law at a given time
   * \param t Time used to evaluate the timing law.
   */
  double eval(double t) const override {
    if (t <= ta_) {
      return s_0(t);
    } else if(t <= tb_) {
      return s_1(t);
    }
    return s_2(t);
  }

  /**
   * \brief Evaluate the timing law derivative at a given time
   * \param t Time used to evaluate the timing law derivative.
   */
  double eval_dt(double t) const override {
    if (t <= ta_) {
      return sdot_0(t);
    } else if (t <= tb_) {
      return sdot_0(ta_);
    }
    return sdot_2(t);
  }

 private:
  /**
   * Evaluate the timing law at \param t with \param t less than
   * or equal to \param ta_.
   */
  double s_0(double t) const {
    return a_max_ * 
        (std::pow(t, 2.0) / 2.0 - t0_ * t + std::pow(t0_, 2.0) / 2.0);
  }

  /**
   * Evaluate the timing law at \param t with \param t between
   * \param ta_ and \param tb_.
   */
  double s_1(double t) const {
    return s_0(ta_) + a_max_ * (ta_ - t0_) * (t - ta_);
  }

  /**
   * Evaluate the timing law at \param t with \param t greater than
   * or equal to \param tb_.
   */
  double s_2(double t) const {
    return s_1(t) - a_max_ *
        (std::pow(t, 2.0) / 2.0 - tb_ * t + std::pow(tb_, 2.0) / 2.0);
  }

  /**
   * Evaluate the timing law derivative at \param t with \param t less than
   * or equal to \param ta_.
   */
  double sdot_0(double t) const {
    return a_max_ * (t - t0_);
  }

  /**
   * Evaluate the timing law derivative at \param t with \param t greater than
   * or equal to \param tb_.
   */
  double sdot_2(double t) const {
    return sdot_0(ta_) - a_max_ * (t - tb_);
  }

  double duration_;
  double alpha_;
  double t0_, ta_, tb_, tf_;
  double a_max_;

}; // end class TrapezoidalAccelerationTimingLaw

} // end namespace labrob