/**
@file AdaptiveTimeController.h
@brief Adaptive time stepping controller for THAMES hydration simulations.

This controller adjusts simulation timestep based on GEMS solver feedback
to minimize failures while maximizing efficiency.

Key features:
- PI-like control based on iteration counts
- Smoothing via moving average of recent history
- Distinguishes stiffness failures from structural failures
- Configurable bounds and parameters

@author Jeffrey W. Bullard
@author Claude (Anthropic)
@date January 2026
*/

#ifndef ADAPTIVETIMECONTROLLER_H
#define ADAPTIVETIMECONTROLLER_H

#include <deque>
#include <string>
#include <vector>

/**
@brief GEMS solver result classification for adaptive stepping decisions.

These categories help the controller decide how aggressively to adjust
the timestep based on how the GEMS solver performed.
*/
enum class GEMSResultType {
  SUCCESS_EASY,       ///< Converged quickly (iterations < target)
  SUCCESS_NORMAL,     ///< Converged normally
  SUCCESS_HARD,       ///< Converged but struggled (iterations > warning threshold)
  FAILURE_STIFFNESS,  ///< Failed due to stiffness (GEMS codes 4, 8)
  FAILURE_STRUCTURAL, ///< Failed due to structural issues (other codes)
  FAILURE_TERMINAL    ///< Terminal error (GEMS code 9)
};

/**
@brief Configuration parameters for adaptive time stepping.

These parameters control how aggressively the controller adjusts timesteps
and what thresholds trigger changes. Default values are suitable for most
cement hydration simulations.
*/
struct AdaptiveTimeConfig {

  // Timestep bounds (hours)
  double dt_min = 1.0e-6;    ///< Minimum timestep (~3.6 ms)
  double dt_max = 1.0;       ///< Maximum timestep (1 hour)
  double dt_initial = 0.001; ///< Initial timestep (~3.6 seconds)

  // Growth/shrink factors
  double growth_factor = 1.2;      ///< Multiply dt by this on easy success
  double shrink_factor = 0.5;      ///< Multiply dt by this on stiffness failure
  double hard_shrink_factor = 0.7; ///< Multiply dt on hard success (preemptive)

  // Iteration thresholds. GEMS3K IPM solver reports the number of
  // iterations it needed to converge each equilibration call; these two
  // thresholds partition that count into a three-state signal.
  // Empirically calibrated against GEMS3K's default IPM tuning in
  // Sessions 23-38. If a different IPM solver is ever swapped in, or if
  // GEMS3K's solver parameters change substantially, these will need to
  // be retuned.
  long int target_iterations = 500;   ///< Below this = "easy" convergence
  long int warning_iterations = 5000; ///< Above this = "struggling"

  // Consecutive success requirement for growth
  int successes_for_growth = 3; ///< Must succeed this many times before growing

  // History size for moving average calculations
  int history_size = 5;

  // Maximum consecutive failures before giving up
  int max_consecutive_failures = 50;

  // Logging control
  bool verbose = false;
};

/**
@class AdaptiveTimeController
@brief Controls timestep size based on GEMS solver feedback.

This class implements a feedback-based adaptive time stepping algorithm
for THAMES hydration simulations. It monitors GEMS solver performance
(iteration counts, convergence metrics, error codes) and adjusts the
timestep to balance efficiency with reliability.

The algorithm uses a PI-like control strategy:
- On easy success (low iterations): grow timestep if consistent
- On hard success (high iterations): shrink timestep preemptively
- On failure: shrink timestep based on failure type

Usage:
@code
AdaptiveTimeController controller;

while (simulating) {
    double dt = controller.getNextTimestep();

    // ... run kinetics and GEMS ...

    if (gemsSucceeded) {
        controller.recordSuccess(iterDone, pci, dxm);
        // advance simulation
    } else {
        double newDt = controller.recordFailure(errorCode);
        // restore state and retry
    }
}
@endcode

@see AdaptiveTimeConfig for configuration options
@see GEMSResultType for result classification
*/
class AdaptiveTimeController {

public:
  /**
  @brief Construct controller with configuration.
  @param config Configuration parameters (uses defaults if not specified)
  */
  explicit AdaptiveTimeController(
      const AdaptiveTimeConfig &config = AdaptiveTimeConfig());

  /**
  @brief Destructor.
  */
  ~AdaptiveTimeController() = default;

  // =========================================================================
  // Main Interface
  // =========================================================================

  /**
  @brief Get the next proposed timestep.

  Returns the current adaptive timestep, which has been adjusted
  based on recent solver performance.

  @return Proposed timestep in hours
  */
  double getNextTimestep() const;

  /**
  @brief Record a successful GEMS calculation.

  Updates internal state based on how hard GEMS worked to converge.
  May adjust future timesteps proactively.

  @param iterDone Total iterations performed by GEMS
  @param pci Dikin's criterion value at convergence
  @param dxm Convergence threshold used
  */
  void recordSuccess(long int iterDone, double pci, double dxm);

  /**
  @brief Record a failed GEMS calculation.

  Shrinks timestep based on failure type:
  - Stiffness failures (codes 4, 8): moderate reduction (shrink_factor)
  - Structural failures: aggressive reduction (shrink_factor^2)
  - Terminal errors (code 9): very aggressive reduction

  @param gemsErrorCode GEMS NodeStatusCH error code
  @return New timestep to use for retry
  */
  double recordFailure(int gemsErrorCode);

  /**
  @brief Reset controller state.

  Clears history and resets timestep to initial value.
  Call this when starting a new simulation or after major state changes.
  */
  void reset();

  /**
  @brief Check if controller has hit maximum consecutive failures.
  @return true if max failures reached, simulation should abort
  */
  bool hasExceededMaxFailures() const;

  // =========================================================================
  // Diagnostics
  // =========================================================================

  /**
  @brief Get current controller state as string for logging.
  @return Diagnostic string with current state information
  */
  std::string getStatusString() const;

  /**
  @brief Get statistics about recent performance.
  @param[out] avgIterations Average iterations over recent history
  @param[out] successRate Success rate over all attempts
  @param[out] avgTimestep Average timestep used recently
  */
  void getStatistics(double &avgIterations, double &successRate,
                     double &avgTimestep) const;

  /**
  @brief Get current timestep value.
  @return Current timestep in hours
  */
  double getCurrentTimestep() const { return dt_current_; }

  /**
  @brief Get number of consecutive successes.
  @return Count of consecutive successful GEMS calculations
  */
  int getConsecutiveSuccesses() const { return consecutive_successes_; }

  /**
  @brief Get number of consecutive failures.
  @return Count of consecutive failed GEMS calculations
  */
  int getConsecutiveFailures() const { return consecutive_failures_; }

  /**
  @brief Get total number of successful calculations.
  @return Lifetime count of successes
  */
  long int getTotalSuccesses() const { return total_successes_; }

  /**
  @brief Get total number of failed calculations.
  @return Lifetime count of failures
  */
  long int getTotalFailures() const { return total_failures_; }

  // =========================================================================
  // Configuration
  // =========================================================================

  /**
  @brief Update configuration parameters.
  @param config New configuration to use
  */
  void setConfig(const AdaptiveTimeConfig &config);

  /**
  @brief Get current configuration.
  @return Reference to current configuration
  */
  const AdaptiveTimeConfig &getConfig() const { return config_; }

  /**
  @brief Set minimum timestep bound.
  @param dt_min Minimum timestep in hours
  */
  void setMinTimestep(double dt_min) { config_.dt_min = dt_min; }

  /**
  @brief Set maximum timestep bound.
  @param dt_max Maximum timestep in hours
  */
  void setMaxTimestep(double dt_max) { config_.dt_max = dt_max; }

  /**
  @brief Enable or disable verbose logging.
  @param verbose true to enable logging
  */
  void setVerbose(bool verbose) { config_.verbose = verbose; }

  /**
  @brief Set initial timestep based on kinetic rates.

  Calculates an appropriate initial timestep to limit the relative change
  in dissolving phases to a specified maximum. This provides a physics-based
  initial timestep rather than using a hard-coded default.

  @param maxRate Maximum dissolution rate across all phases [1/hour]
  @param maxRelativeChange Maximum allowed relative change per timestep (default 5%)
  */
  void setInitialTimestepFromKinetics(double maxRate,
                                      double maxRelativeChange = 0.05);

private:
  /// Configuration parameters
  AdaptiveTimeConfig config_;

  /// Current timestep (hours)
  double dt_current_;

  /// Count of consecutive successful calculations
  int consecutive_successes_;

  /// Count of consecutive failed calculations
  int consecutive_failures_;

  /// Last PCI value from GEMS
  double last_pci_;

  /// Last DXM value from GEMS
  double last_dxm_;

  /// Recent iteration counts for moving average
  std::deque<long int> iteration_history_;

  /// Recent timesteps used
  std::deque<double> timestep_history_;

  /// Recent success/failure record
  std::deque<bool> success_history_;

  /// Lifetime count of successful calculations
  long int total_successes_;

  /// Lifetime count of failed calculations
  long int total_failures_;

  /// Total calculation attempts
  long int total_steps_;

  // =========================================================================
  // Helper Methods
  // =========================================================================

  /**
  @brief Classify a GEMS result for decision-making.
  @param iterDone Iteration count (0 for failures)
  @param errorCode GEMS error code (0 for success)
  @return Classification of the result
  */
  GEMSResultType classifyResult(long int iterDone, int errorCode) const;

  /**
  @brief Compute moving average of recent iteration counts.
  @return Average iterations, or 0 if no history
  */
  double computeMovingAverage() const;

  /**
  @brief Update history buffers with new data point.
  @param iterations Iteration count from this calculation
  @param timestep Timestep used for this calculation
  @param success Whether calculation succeeded
  */
  void updateHistory(long int iterations, double timestep, bool success);

  /**
  @brief Adjust timestep after a successful calculation.
  @param resultType Classification of the success
  */
  void adjustTimestepOnSuccess(GEMSResultType resultType);

  /**
  @brief Adjust timestep after a failed calculation.
  @param resultType Classification of the failure
  */
  void adjustTimestepOnFailure(GEMSResultType resultType);

  /**
  @brief Clamp timestep to configured bounds.
  @param dt Proposed timestep
  @return Timestep clamped to [dt_min, dt_max]
  */
  double clampTimestep(double dt) const;
};

#endif // ADAPTIVETIMECONTROLLER_H
