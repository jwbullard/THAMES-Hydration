/**
@file SumTree.h
@brief Fenwick (binary indexed) tree over integer weights, used for the
weighted random site picks in Lattice::growPhase, dissolvePhase and
transformPhase.

Each pick used to be a linear cumulative scan over the whole candidate vector
(growthVector[].affinityInt or dissolutionVector[].wmcInt). Both the vector
length and the number of picks per cycle grow as N^3 for an N^3 lattice, so
the per-cycle cost grew as N^6. The tree answers a pick and absorbs a weight
change in O(log n), which brings the per-cycle cost to about N^3.

The tree mirrors the weights stored in the local candidate vectors. Every
append (push), swap-with-last-and-pop removal (swapRemove), weight change
(set) and full rebuild of those vectors must be mirrored by the matching call
here; THAMES_SUMTREE_CHECK verifies that at run time.

Selection rule. With integer weights w_i and total S, a pick with a uniform
deviate rng in [0, 1) returns the smallest index i whose prefix sum
w_0 + ... + w_i is >= rng*S, which is the rule the old scan implemented in
floating point. The old scan accumulated w_i/S in double precision, so a draw
within rounding error of a boundary could land on the neighbouring index;
the tree uses exact integer prefix sums. Both sample the same distribution.
The two were shown to give byte-identical output at N = 32, 50 and 100 when
the tree deferred to the old scan on such near-ties (2026-10-12, before the
old scan was removed); since then runs agree with the old code statistically,
as an independent replicate would, rather than draw for draw.

Ported from the MSMSE-109162 revision (2026-10), where it was developed on
commit c6f84fb.

Environment variable THAMES_SUMTREE_CHECK (bit flags, default 0):
  1  per pick, verify that the tree's weights, size and total equal the
     vector's (O(n); aborts on a mismatch)
  2  per pick, recompute the pick by an exact integer prefix scan over the
     vector and abort if it differs from the tree's
When set, a one-line summary per function is written to stderr at exit.
*/

#ifndef THAMES_SUMTREE_H
#define THAMES_SUMTREE_H

#include <cstdio>
#include <cstdlib>
#include <vector>

namespace sumtree {

inline int checkMode() {
  static int mode = [] {
    const char *e = std::getenv("THAMES_SUMTREE_CHECK");
    return e ? std::atoi(e) : 0;
  }();
  return mode;
}

/// Draw statistics, reported at exit only when THAMES_SUMTREE_CHECK is set.
struct Stats {
  const char *name;
  long long draws = 0;
  long long maxN = 0;
  explicit Stats(const char *n) : name(n) {}
  ~Stats() {
    if (checkMode() != 0 && draws > 0)
      std::fprintf(stderr, "SumTree %s: checked draws=%lld maxn=%lld\n", name,
                   draws, maxN);
  }
};

inline Stats &growStats() {
  static Stats s("grow");
  return s;
}
inline Stats &dissStats() {
  static Stats s("diss");
  return s;
}
inline Stats &transStats() {
  static Stats s("transform");
  return s;
}

class SumTree {
public:
  /// rng*S exceeds the total, which floating-point rounding of rng*S can
  /// cause; callers clamp it as the old scan did.
  static const int PAST_END = -2;

  int size() const { return n_; }
  long long total() const { return total_; }
  long long weight(int i) const { return w_[i]; }

  /// Rebuild from a vector of structs in O(n); f(elem) returns its weight.
  template <class V, class F> void rebuild(const V &vec, F f) {
    n_ = static_cast<int>(vec.size());
    cap_ = 1;
    while (cap_ < n_)
      cap_ <<= 1;
    w_.assign(cap_, 0);
    total_ = 0;
    for (int i = 0; i < n_; i++) {
      w_[i] = f(vec[i]);
      total_ += w_[i];
    }
    build();
  }

  /// Set the weight of entry i (0-based, i < size()).
  void set(int i, long long val) {
    long long d = val - w_[i];
    if (d == 0)
      return;
    w_[i] = val;
    total_ += d;
    for (int k = i + 1; k <= cap_; k += k & (-k))
      t_[k] += d;
  }

  /// Append an entry with weight val.
  void push(long long val) {
    if (n_ == cap_) {
      cap_ = (cap_ < 1024) ? 1024 : 2 * cap_;
      w_.resize(cap_, 0);
      build();
    }
    n_++;
    set(n_ - 1, val);
  }

  /// Mirror of "v[i] = v.back(); v.pop_back();"
  void swapRemove(int i) {
    int last = n_ - 1;
    if (i != last)
      set(i, w_[last]);
    set(last, 0);
    n_--;
  }

  /// Smallest index whose prefix sum is >= rng*S, or PAST_END. Weights must
  /// be non-negative (growth affinities lie in [0, 1] and wmc counts are
  /// non-negative by construction); a negative weight means the mirror is
  /// broken, and the run is stopped rather than sampled wrongly.
  int pick(double rng, double S) const {
    double target = rng * S;
    int pos = 0;
    long long acc = 0;
    for (int step = cap_; step > 0; step >>= 1) {
      int nx = pos + step;
      if (nx <= cap_ && static_cast<double>(acc + t_[nx]) < target) {
        pos = nx;
        acc += t_[nx];
      }
    }
    if (pos >= n_)
      return PAST_END;
    if (w_[pos] < 0) {
      std::fprintf(stderr,
                   "SumTree: negative weight %lld at index %d; candidate "
                   "weights must be non-negative\n",
                   w_[pos], pos);
      std::abort();
    }
    return pos;
  }

private:
  void build() {
    t_.assign(cap_ + 1, 0);
    for (int i = 1; i <= cap_; i++) {
      t_[i] += w_[i - 1];
      int j = i + (i & (-i));
      if (j <= cap_)
        t_[j] += t_[i];
    }
  }

  int n_ = 0;
  int cap_ = 1;
  long long total_ = 0;
  std::vector<long long> w_;
  std::vector<long long> t_;
};

/// THAMES_SUMTREE_CHECK != 0: verify the tree mirrors the vector (bit 1) and
/// that the pick equals an exact integer prefix scan (bit 2). Aborts on any
/// discrepancy, naming the calling function.
template <class V, class F>
inline void check(const SumTree &tree, const V &vec, F f, long long runSum,
                  double rng, double S, int pickedPos, Stats &st,
                  const char *where) {
  st.draws++;
  if (static_cast<long long>(vec.size()) > st.maxN)
    st.maxN = static_cast<long long>(vec.size());
  const int mode = checkMode();
  if (mode & 1) {
    if (static_cast<int>(vec.size()) != tree.size()) {
      std::fprintf(stderr, "SumTree %s: size mismatch tree=%d vec=%zu\n",
                   where, tree.size(), vec.size());
      std::abort();
    }
    for (int i = 0; i < tree.size(); i++) {
      if (tree.weight(i) != f(vec[i])) {
        std::fprintf(stderr,
                     "SumTree %s: weight mismatch at %d tree=%lld vec=%lld\n",
                     where, i, tree.weight(i), (long long)f(vec[i]));
        std::abort();
      }
    }
    if (tree.total() != runSum) {
      std::fprintf(stderr, "SumTree %s: tree total %lld != running sum %lld\n",
                   where, tree.total(), runSum);
      std::abort();
    }
  }
  if (mode & 2) {
    const double target = rng * S;
    long long acc = 0;
    int exact = SumTree::PAST_END;
    for (int i = 0; i < static_cast<int>(vec.size()); i++) {
      acc += f(vec[i]);
      if (static_cast<double>(acc) >= target) {
        exact = i;
        break;
      }
    }
    if (exact != pickedPos) {
      std::fprintf(stderr, "SumTree %s: pick %d != exact prefix scan %d\n",
                   where, pickedPos, exact);
      std::abort();
    }
  }
}

} // namespace sumtree

#endif
