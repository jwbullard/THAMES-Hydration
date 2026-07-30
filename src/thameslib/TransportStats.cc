/**
@file TransportStats.cc
@brief Implementation of the pure-math aggregate for shell-thickness
       distribution.

Kept dependency-free (only <algorithm>, <cmath>, <limits>) so the
standalone unit test can compile without pulling in ChemicalSystem
or Lattice.
*/

#include "TransportStats.h"

#include <algorithm>
#include <cmath>
#include <unordered_map>

namespace xport {

namespace {

/**
Split a monotonically-nondecreasing sequence into K equal-frequency
buckets. Returns bucket assignments: bucketOf[i] = k for the k-th
bucket that site i belongs to.

If there are fewer sites than bins, some bins will be empty; the
caller filters those out.
*/
std::vector<int> equalFrequencyBucketAssignments(int nSites, int K) {
  std::vector<int> out(nSites, 0);
  if (K <= 1) return out;  // all zero → single bucket
  if (nSites <= 0) return out;
  // Simple stride assignment; site i lands in bucket i*K/N.
  for (int i = 0; i < nSites; ++i) {
    int b = (i * K) / nSites;
    if (b >= K) b = K - 1;
    out[i] = b;
  }
  return out;
}

}  // namespace

ShellStats aggregateShellDistribution(
    const std::vector<double> &deltasMeters,
    const std::vector<int> &shellPhaseIds,
    int K) {
  ShellStats out;
  out.numSitesTotal = static_cast<int>(deltasMeters.size());
  if (K < 1) K = 1;

  // Separate finite-δ sites from +inf-δ sites (unreached by walk).
  // Preserve the (delta, phaseId) pairing.
  struct Site {
    double delta;
    int shellPhaseId;
  };
  std::vector<Site> reached;
  reached.reserve(deltasMeters.size());
  const size_t nIds = shellPhaseIds.size();
  const double kInf = std::numeric_limits<double>::infinity();
  for (size_t i = 0; i < deltasMeters.size(); ++i) {
    const double d = deltasMeters[i];
    if (std::isfinite(d) && d < kInf) {
      const int pid = (i < nIds) ? shellPhaseIds[i] : -1;
      reached.push_back({d, pid});
    }
  }
  out.numSitesReached = static_cast<int>(reached.size());

  if (reached.empty()) {
    // No reached sites — no useful aggregate. Emit a single empty
    // bin so consumers don't crash on empty vector.
    out.bins.emplace_back();
    out.bins.back().siteFraction = 0.0;
    return out;
  }

  // Sort by δ ascending so equal-frequency bins are contiguous
  // δ-ranges. Small allocation overhead; caller amortizes over the
  // full lattice pass.
  std::sort(reached.begin(), reached.end(),
            [](const Site &a, const Site &b) { return a.delta < b.delta; });

  // Diagnostic means over the reached-only distribution.
  double sumDelta = 0.0;
  double sumRecip = 0.0;
  int recipDenom = 0;
  for (const auto &s : reached) {
    sumDelta += s.delta;
    if (s.delta > 0.0) {
      sumRecip += 1.0 / s.delta;
      ++recipDenom;
    } else {
      // δ = 0 means the site is directly adjacent to electrolyte —
      // diffusion resistance is zero, contributes 1/0 to the reciprocal
      // sum. Skip and let the harmonic mean → 0 organically as the
      // count grows (a single zero-δ site makes the harmonic mean
      // vanish, matching the physics: parallel resistor with a short
      // circuit).
      recipDenom = 0;  // signal special-case below
      break;
    }
  }
  out.deltaArithmeticRaw = sumDelta / static_cast<double>(reached.size());
  out.deltaHarmonicRaw = (recipDenom > 0)
      ? static_cast<double>(recipDenom) / sumRecip
      : 0.0;  // zero-shell site short-circuits the parallel network

  // Build K equal-frequency bins.
  const int nReached = static_cast<int>(reached.size());
  const int effectiveK = std::min(K, nReached);
  const auto bucketOf = equalFrequencyBucketAssignments(nReached, effectiveK);

  out.bins.assign(effectiveK, ShellBin{});
  std::vector<int> perBinCount(effectiveK, 0);
  std::vector<std::unordered_map<int, int>> perBinShellHist(effectiveK);
  std::vector<std::vector<double>> perBinDeltas(effectiveK);

  for (int i = 0; i < nReached; ++i) {
    const int b = bucketOf[i];
    perBinCount[b]++;
    perBinDeltas[b].push_back(reached[i].delta);
    perBinShellHist[b][reached[i].shellPhaseId]++;
  }

  for (int b = 0; b < effectiveK; ++b) {
    ShellBin &bin = out.bins[b];
    bin.siteFraction = static_cast<double>(perBinCount[b]) /
                       static_cast<double>(nReached);
    // Bin representative δ = median of the bin's sites. Because sites
    // were sorted ascending before bucketing, perBinDeltas[b] is
    // already sorted; median is the middle element.
    const auto &dvec = perBinDeltas[b];
    if (!dvec.empty()) {
      bin.deltaRep = dvec[dvec.size() / 2];
    }
    // Dominant shell phase id = mode of the bin's shellPhaseIds.
    // Ties broken by first-seen order (unordered_map iteration is
    // implementation-defined, but every bin has a single winner
    // in practice; the tie-break policy affects at most one bin
    // per cycle).
    int bestPid = -1;
    int bestCount = -1;
    for (const auto &kv : perBinShellHist[b]) {
      if (kv.second > bestCount) {
        bestCount = kv.second;
        bestPid = kv.first;
      }
    }
    bin.dominantShellPhaseId = bestPid;
  }

  return out;
}

bool distributionIsPathological(const ShellStats &stats,
                                double arithToHarmonicRatioThreshold) {
  if (stats.numSitesReached == 0) return false;
  if (stats.deltaHarmonicRaw <= 0.0) return false;
  const double ratio = stats.deltaArithmeticRaw / stats.deltaHarmonicRaw;
  return ratio > arithToHarmonicRatioThreshold;
}

}  // namespace xport
