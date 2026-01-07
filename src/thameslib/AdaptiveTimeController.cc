/**
@file AdaptiveTimeController.cc
@brief Implementation of adaptive time stepping controller for THAMES.
*/

#include "AdaptiveTimeController.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>

AdaptiveTimeController::AdaptiveTimeController(const AdaptiveTimeConfig &config)
    : config_(config), dt_current_(config.dt_initial), consecutive_successes_(0),
      consecutive_failures_(0), last_pci_(0.0), last_dxm_(1.0e-5),
      total_successes_(0), total_failures_(0), total_steps_(0) {}

double AdaptiveTimeController::getNextTimestep() const { return dt_current_; }

void AdaptiveTimeController::recordSuccess(long int iterDone, double pci,
                                           double dxm) {
  total_successes_++;
  total_steps_++;
  consecutive_failures_ = 0;
  consecutive_successes_++;
  last_pci_ = pci;
  last_dxm_ = dxm;

  // Update history
  updateHistory(iterDone, dt_current_, true);

  // Classify result and adjust timestep
  GEMSResultType result = classifyResult(iterDone, 0);
  adjustTimestepOnSuccess(result);

  if (config_.verbose) {
    std::clog << "AdaptiveTime: SUCCESS"
              << " iter=" << iterDone << " pci=" << std::scientific
              << std::setprecision(2) << pci << " ratio=" << std::fixed
              << std::setprecision(1) << (dxm > 0 ? pci / dxm : -1)
              << " dt=" << std::scientific << std::setprecision(3) << dt_current_
              << "h" << std::endl;
  }
}

double AdaptiveTimeController::recordFailure(int gemsErrorCode) {
  total_failures_++;
  total_steps_++;
  consecutive_successes_ = 0;
  consecutive_failures_++;

  // Update history
  updateHistory(0, dt_current_, false);

  // Classify and adjust
  GEMSResultType result = classifyResult(0, gemsErrorCode);
  adjustTimestepOnFailure(result);

  if (config_.verbose) {
    std::clog << "AdaptiveTime: FAILURE"
              << " code=" << gemsErrorCode
              << " consec=" << consecutive_failures_
              << " new_dt=" << std::scientific << std::setprecision(3)
              << dt_current_ << "h" << std::endl;
  }

  return dt_current_;
}

void AdaptiveTimeController::reset() {
  dt_current_ = config_.dt_initial;
  consecutive_successes_ = 0;
  consecutive_failures_ = 0;
  last_pci_ = 0.0;
  last_dxm_ = 1.0e-5;
  iteration_history_.clear();
  timestep_history_.clear();
  success_history_.clear();
  // Note: Don't reset total_successes_, total_failures_, total_steps_
  // Those track lifetime statistics across resets
}

bool AdaptiveTimeController::hasExceededMaxFailures() const {
  return consecutive_failures_ >= config_.max_consecutive_failures;
}

GEMSResultType AdaptiveTimeController::classifyResult(long int iterDone,
                                                      int errorCode) const {
  if (errorCode == 0) {
    // Success - classify by difficulty
    if (iterDone < config_.target_iterations) {
      return GEMSResultType::SUCCESS_EASY;
    } else if (iterDone < config_.warning_iterations) {
      return GEMSResultType::SUCCESS_NORMAL;
    } else {
      return GEMSResultType::SUCCESS_HARD;
    }
  } else {
    // Failure - classify by type based on GEMS NodeStatusCH codes:
    // 2 = OK_GEM_AIA (success, but we shouldn't be here)
    // 3 = BAD_GEM_AIA (warning)
    // 4 = ERR_GEM_AIA (failure with automatic IA)
    // 6 = OK_GEM_SIA (success)
    // 7 = BAD_GEM_SIA (warning)
    // 8 = ERR_GEM_SIA (failure with smart IA)
    // 9 = T_ERROR_GEM (terminal error)

    if (errorCode == 9) {
      return GEMSResultType::FAILURE_TERMINAL;
    } else if (errorCode == 4 || errorCode == 8) {
      // Standard GEMS failures - typically stiffness/convergence issues
      return GEMSResultType::FAILURE_STIFFNESS;
    } else {
      // Other failures (3, 7 are warnings but treated as structural)
      return GEMSResultType::FAILURE_STRUCTURAL;
    }
  }
}

double AdaptiveTimeController::computeMovingAverage() const {
  if (iteration_history_.empty())
    return 0.0;

  double sum = std::accumulate(iteration_history_.begin(),
                               iteration_history_.end(), 0.0);
  return sum / static_cast<double>(iteration_history_.size());
}

void AdaptiveTimeController::updateHistory(long int iterations, double timestep,
                                           bool success) {
  // Only add to iteration history on success
  if (success) {
    iteration_history_.push_back(iterations);
    while (static_cast<int>(iteration_history_.size()) > config_.history_size) {
      iteration_history_.pop_front();
    }
  }

  timestep_history_.push_back(timestep);
  while (static_cast<int>(timestep_history_.size()) > config_.history_size) {
    timestep_history_.pop_front();
  }

  success_history_.push_back(success);
  // Keep longer success history for rate calculation
  while (static_cast<int>(success_history_.size()) >
         config_.history_size * 4) {
    success_history_.pop_front();
  }
}

void AdaptiveTimeController::adjustTimestepOnSuccess(GEMSResultType resultType) {
  switch (resultType) {
  case GEMSResultType::SUCCESS_EASY:
    // GEMS found it easy - consider growing timestep
    if (consecutive_successes_ >= config_.successes_for_growth) {
      dt_current_ *= config_.growth_factor;
      if (config_.verbose) {
        std::clog << "  -> Growing timestep (easy success streak)"
                  << std::endl;
      }
    }
    break;

  case GEMSResultType::SUCCESS_NORMAL:
    // Normal convergence - maintain current timestep
    break;

  case GEMSResultType::SUCCESS_HARD:
    // GEMS struggled - preemptively reduce timestep
    dt_current_ *= config_.hard_shrink_factor;
    consecutive_successes_ = 0; // Reset growth counter
    if (config_.verbose) {
      std::clog << "  -> Shrinking timestep (hard success)" << std::endl;
    }
    break;

  default:
    break;
  }

  dt_current_ = clampTimestep(dt_current_);
}

void AdaptiveTimeController::adjustTimestepOnFailure(GEMSResultType resultType) {
  double old_dt = dt_current_;

  switch (resultType) {
  case GEMSResultType::FAILURE_STIFFNESS:
    // Stiffness-related failure - standard reduction
    dt_current_ *= config_.shrink_factor;
    break;

  case GEMSResultType::FAILURE_STRUCTURAL:
    // Structural failure - more aggressive reduction
    dt_current_ *= config_.shrink_factor * config_.shrink_factor;
    break;

  case GEMSResultType::FAILURE_TERMINAL:
    // Terminal error - very aggressive reduction
    dt_current_ *= 0.1;
    break;

  default:
    dt_current_ *= config_.shrink_factor;
    break;
  }

  dt_current_ = clampTimestep(dt_current_);

  // Clear iteration history on failure - past performance is no longer relevant
  // because we're now in a different regime
  iteration_history_.clear();

  if (config_.verbose && dt_current_ == config_.dt_min && old_dt > config_.dt_min) {
    std::clog << "  -> WARNING: Hit minimum timestep!" << std::endl;
  }
}

double AdaptiveTimeController::clampTimestep(double dt) const {
  return std::max(config_.dt_min, std::min(dt, config_.dt_max));
}

std::string AdaptiveTimeController::getStatusString() const {
  std::ostringstream oss;
  oss << std::scientific << std::setprecision(3);
  oss << "AdaptiveTimeController: "
      << "dt=" << dt_current_ << "h"
      << ", consec_success=" << consecutive_successes_
      << ", consec_fail=" << consecutive_failures_ << ", total=" << total_steps_
      << " (" << total_successes_ << " ok, " << total_failures_ << " fail)";

  if (!iteration_history_.empty()) {
    oss << ", avg_iter=" << std::fixed << std::setprecision(0)
        << computeMovingAverage();
  }

  double successRate = 0.0;
  if (total_steps_ > 0) {
    successRate = 100.0 * static_cast<double>(total_successes_) /
                  static_cast<double>(total_steps_);
    oss << ", success_rate=" << std::fixed << std::setprecision(1)
        << successRate << "%";
  }

  return oss.str();
}

void AdaptiveTimeController::getStatistics(double &avgIterations,
                                           double &successRate,
                                           double &avgTimestep) const {
  avgIterations = computeMovingAverage();

  if (total_steps_ > 0) {
    successRate =
        static_cast<double>(total_successes_) / static_cast<double>(total_steps_);
  } else {
    successRate = 0.0;
  }

  if (!timestep_history_.empty()) {
    avgTimestep = std::accumulate(timestep_history_.begin(),
                                  timestep_history_.end(), 0.0) /
                  static_cast<double>(timestep_history_.size());
  } else {
    avgTimestep = dt_current_;
  }
}

void AdaptiveTimeController::setConfig(const AdaptiveTimeConfig &config) {
  config_ = config;
  // Ensure current timestep is within new bounds
  dt_current_ = clampTimestep(dt_current_);
}
