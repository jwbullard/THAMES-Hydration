// Standalone unit test for SumTree.h, the Fenwick tree behind the weighted
// site picks in Lattice::growPhase, dissolvePhase and transformPhase.
// Header only; no THAMES dependencies. Compile with build_and_run.sh.
//
// Checks, under random sequences of the operations Lattice performs (append,
// swap-with-last removal, weight update, rebuild):
//   1. the tree mirrors the vector exactly (every weight, the size, the total);
//   2. every pick satisfies the exact integer rule
//      prefix(i-1) < rng*S <= prefix(i);
//   3. where a pick differs from the floating-point cumulative scan Lattice
//      used before the tree (kept here only as a reference), rng*S lies
//      within rounding error of a boundary: the two sample the same
//      distribution and differ only on near-ties;
//   4. draws placed exactly on boundaries obey the exact rule;
//   5. zero-weight entries are never picked.

#include "SumTree.h"

#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

namespace {

int gFailed = 0;

void expect(bool cond, const char *tag) {
  std::printf("%-62s %s\n", tag, cond ? "PASS" : "FAIL");
  if (!cond)
    ++gFailed;
}

struct Cand {
  long long w;
};

/// The selection Lattice used before the tree (growPhase/dissolvePhase):
///   p = w[0]/S; i = 0; while (rng > p) p += w[++i]/S;
/// clamped to the last entry if it runs off the end.
int originalScan(const std::vector<Cand> &v, double rng, double S) {
  const int n = static_cast<int>(v.size());
  double p = v[0].w / S;
  if (rng <= p)
    return 0;
  int i = 0;
  while (rng > p && i < n - 1) {
    ++i;
    p += v[i].w / S;
  }
  return i;
}

long long prefix(const std::vector<Cand> &v, int last) {
  long long s = 0;
  for (int i = 0; i <= last; i++)
    s += v[i].w;
  return s;
}

bool mirrors(const sumtree::SumTree &t, const std::vector<Cand> &v) {
  if (t.size() != static_cast<int>(v.size()))
    return false;
  long long tot = 0;
  for (int i = 0; i < t.size(); i++) {
    if (t.weight(i) != v[i].w)
      return false;
    tot += v[i].w;
  }
  return t.total() == tot;
}

struct DrawCounts {
  long long draws = 0, inexact = 0, zeroPick = 0;
  long long scanDiffers = 0, scanDiffersAwayFromBoundary = 0;
};

void draw(const sumtree::SumTree &t, const std::vector<Cand> &v, double rng,
          DrawCounts &c) {
  const double S = static_cast<double>(t.total());
  int pos = t.pick(rng, S);
  if (pos == sumtree::SumTree::PAST_END)
    pos = static_cast<int>(v.size()) - 1;
  ++c.draws;
  const double target = rng * S;
  const double lo = (pos > 0) ? static_cast<double>(prefix(v, pos - 1)) : 0.0;
  const double hi = static_cast<double>(prefix(v, pos));
  if (!(lo < target && target <= hi) && pos != static_cast<int>(v.size()) - 1)
    ++c.inexact;
  if (v[pos].w == 0)
    ++c.zeroPick;
  const int scan = originalScan(v, rng, S);
  if (scan != pos) {
    ++c.scanDiffers;
    // the boundary between the two picks
    const int b = (scan < pos) ? scan : pos;
    const double edge = static_cast<double>(prefix(v, b));
    if (std::fabs(target - edge) > 1e-9 * S)
      ++c.scanDiffersAwayFromBoundary;
  }
}

void testRandomOperations() {
  std::printf("\n--- random append / swap-remove / update / rebuild ---\n");
  std::mt19937_64 gen(20261012);
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  std::uniform_int_distribution<int> wdist(0, 40); // includes zero weights

  std::vector<Cand> v;
  sumtree::SumTree t;
  auto weightOf = [](const Cand &c) { return c.w; };
  t.rebuild(v, weightOf);

  bool alwaysMirrors = true;
  DrawCounts c;
  for (int step = 0; step < 200000; step++) {
    const double op = unif(gen);
    if (op < 0.40 || v.size() < 2) {
      Cand x{wdist(gen)};
      v.push_back(x);
      t.push(x.w);
    } else if (op < 0.65) {
      const int i = static_cast<int>(unif(gen) * v.size());
      v[i] = v.back();
      v.pop_back();
      t.swapRemove(i);
    } else if (op < 0.995) {
      const int i = static_cast<int>(unif(gen) * v.size());
      v[i].w = wdist(gen);
      t.set(i, v[i].w);
    } else {
      t.rebuild(v, weightOf);
    }
    if (step % 97 == 0 && !mirrors(t, v))
      alwaysMirrors = false;
    if (!v.empty() && t.total() > 0)
      draw(t, v, unif(gen), c);
  }
  expect(alwaysMirrors && mirrors(t, v), "tree mirrors vector after every op");
  std::printf("  draws %lld; differ from the old scan %lld\n", c.draws,
              c.scanDiffers);
  expect(c.inexact == 0, "every pick obeys the exact integer rule");
  expect(c.scanDiffersAwayFromBoundary == 0,
         "any difference from the old scan is a near-tie");
  expect(c.zeroPick == 0, "no zero-weight entry is picked");
}

void testBoundaryDraws() {
  std::printf("\n--- draws exactly on boundaries ---\n");
  std::vector<Cand> v;
  for (int i = 0; i < 5000; i++)
    v.push_back(Cand{1 + (i * 7919) % 13});
  sumtree::SumTree t;
  t.rebuild(v, [](const Cand &c) { return c.w; });
  const double S = static_cast<double>(t.total());
  DrawCounts c;
  long long run = 0;
  for (int i = 0; i < static_cast<int>(v.size()); i++) {
    run += v[i].w;
    draw(t, v, static_cast<double>(run) / S, c);          // on a boundary
    draw(t, v, std::nextafter(static_cast<double>(run) / S, 0.0), c);
    draw(t, v, std::nextafter(static_cast<double>(run) / S, 2.0), c);
  }
  std::printf("  boundary draws %lld; differ from the old scan %lld\n",
              c.draws, c.scanDiffers);
  expect(c.inexact == 0, "every boundary draw obeys the exact integer rule");
  expect(c.scanDiffersAwayFromBoundary == 0,
         "any difference from the old scan is a near-tie");
}

void testSmallCases() {
  std::printf("\n--- small and degenerate cases ---\n");
  std::vector<Cand> v{{0}, {0}, {5}, {0}};
  sumtree::SumTree t;
  t.rebuild(v, [](const Cand &c) { return c.w; });
  bool ok = true;
  for (double r : {1e-12, 0.3, 0.999999, 1.0}) {
    const int p = t.pick(r, 5.0);
    if (p != 2)
      ok = false;
  }
  expect(ok, "only the single positive weight is ever picked");

  std::vector<Cand> one{{7}};
  t.rebuild(one, [](const Cand &c) { return c.w; });
  expect(t.pick(0.5, 7.0) == 0, "single candidate is picked");

  std::vector<Cand> empty;
  t.rebuild(empty, [](const Cand &c) { return c.w; });
  expect(t.size() == 0 && t.total() == 0, "empty rebuild gives an empty tree");
  t.push(3);
  t.push(4);
  expect(t.size() == 2 && t.total() == 7, "push after empty rebuild");
  t.swapRemove(0);
  expect(t.size() == 1 && t.weight(0) == 4 && t.total() == 4,
         "swap-remove moves the last weight into the hole");
}

} // namespace

int main() {
  testRandomOperations();
  testBoundaryDraws();
  testSmallCases();
  std::printf("\n%s: %d failure(s)\n", gFailed ? "FAILED" : "ALL PASSED",
              gFailed);
  return gFailed ? 1 : 0;
}
