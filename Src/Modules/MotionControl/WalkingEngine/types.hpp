#pragma once

// STL
#include <deque>

// Eigen
#include <Eigen/Core>

#include "Configuration.hpp"

using FootstepPlan = std::deque<Configuration, Eigen::aligned_allocator<Configuration>>;