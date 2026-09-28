// Standalone unit test for self-consistent effective moduli.
// Compile with build_and_run.sh — no THAMES dependencies needed.

#include "Homogenization.h"

#include <cmath>
#include <cstdio>
#include <vector>

namespace {

int gFailed = 0;

void expect(bool cond, const char *tag) {
  std::printf("%-62s %s\n", tag, cond ? "PASS" : "FAIL");
  if (!cond) ++gFailed;
}

void expectNear(double actual, double expected, double tol, const char *tag) {
  const bool ok = std::fabs(actual - expected) <= tol;
  std::printf("%-62s actual=%.5f expected=%.5f %s\n", tag, actual, expected,
              ok ? "PASS" : "FAIL");
  if (!ok) ++gFailed;
}

// Hashin-Shtrikman bounds for a two-phase composite, as an independent check
// that the estimate lands where an estimate must.
double hsBulk(double f1, double K1, double G1, double f2, double K2,
              double Gref) {
  return K1 + f2 / (1.0 / (K2 - K1) + f1 / (K1 + 4.0 * Gref / 3.0));
}

// --- Test group 1: a single phase returns itself ---
void testSinglePhase() {
  std::vector<homog::Phase> p{{1.0, 20.0, 9.0}};
  const homog::Result r = homog::selfConsistent(p);
  expect(r.converged, "single phase: converged");
  expectNear(r.bulkModulus, 20.0, 1e-8, "single phase: K unchanged");
  expectNear(r.shearModulus, 9.0, 1e-8, "single phase: G unchanged");

  // Two copies of the same phase must be indistinguishable from one.
  std::vector<homog::Phase> q{{0.3, 20.0, 9.0}, {0.7, 20.0, 9.0}};
  const homog::Result s = homog::selfConsistent(q);
  expectNear(s.bulkModulus, 20.0, 1e-8, "duplicate phases: K unchanged");

  // Volume fractions need not be normalized.
  std::vector<homog::Phase> u{{3.0, 20.0, 9.0}, {7.0, 20.0, 9.0}};
  const homog::Result t = homog::selfConsistent(u);
  expectNear(t.bulkModulus, 20.0, 1e-8, "unnormalized fractions: K unchanged");
}

// --- Test group 2: porosity softens, and the estimate respects the bounds ---
void testPorosity() {
  const double Ks = 30.0, Gs = 20.0;
  double previous = 1.0e30;
  for (double porosity = 0.0; porosity <= 0.45; porosity += 0.15) {
    std::vector<homog::Phase> p{{1.0 - porosity, Ks, Gs}, {porosity, 0.0, 0.0}};
    const homog::Result r = homog::selfConsistent(p);
    char tag[80];
    std::snprintf(tag, sizeof(tag), "porosity %.2f: softer than the last",
                  porosity);
    expect(r.bulkModulus < previous, tag);
    previous = r.bulkModulus;
  }

  // At 20 % porosity the estimate must sit inside the HS bounds. The lower
  // bound uses the pore as reference and is zero for a void, so the
  // meaningful check is that it stays under the upper bound and above zero.
  const double porosity = 0.2;
  std::vector<homog::Phase> p{{1.0 - porosity, Ks, Gs}, {porosity, 0.0, 0.0}};
  const homog::Result r = homog::selfConsistent(p);
  const double upper = hsBulk(1.0 - porosity, Ks, Gs, porosity, 0.0, Gs);
  expect(r.bulkModulus > 0.0 && r.bulkModulus < upper,
         "porosity 0.20: between zero and the HS upper bound");
}

// --- Test group 3: the percolation limit of the scheme ---
void testPercolationLimit() {
  const double Ks = 30.0, Gs = 20.0;

  // Spherical phases lose all stiffness at a solid fraction of one half.
  // This is a property of the approximation, not of any material, and it is
  // checked here so that nobody mistakes it for one later.
  std::vector<homog::Phase> below{{0.35, Ks, Gs}, {0.65, 0.0, 0.0}};
  const homog::Result rb = homog::selfConsistent(below);
  expectNear(rb.bulkModulus, 0.0, 1e-6, "solid 0.35: no stiffness");

  std::vector<homog::Phase> above{{0.75, Ks, Gs}, {0.25, 0.0, 0.0}};
  const homog::Result ra = homog::selfConsistent(above);
  expect(ra.bulkModulus > 0.5, "solid 0.75: stiff");
}

// --- Test group 4: fluid in the pores, and why the caller leaves it out ---
void testFluidFilledPore() {
  const double Ks = 30.0, Gs = 20.0;
  const double Kw = 2.2; // water

  std::vector<homog::Phase> dry{{0.7, Ks, Gs}, {0.3, 0.0, 0.0}};
  std::vector<homog::Phase> wet{{0.7, Ks, Gs}, {0.3, Kw, 0.0}};
  const homog::Result rd = homog::selfConsistent(dry);
  const homog::Result rw = homog::selfConsistent(wet);

  expect(rw.bulkModulus > rd.bulkModulus,
         "water in the pores stiffens the bulk response");

  // It stiffens the SHEAR response too, which is worth asserting because it
  // catches people out: Gassmann's relation says pore fluid cannot change the
  // shear modulus, but Gassmann describes a connected fluid free to move,
  // whereas the self-consistent scheme embeds each phase as an isolated
  // inclusion. Here K and G are coupled through Feff, so raising the pore
  // fluid's K raises Geff as well. The scheme is behaving correctly; it is
  // simply not modelling a drained pore network.
  //
  // Which is exactly why the DRAINED modulus a poromechanical shrinkage
  // estimate needs is obtained by passing the porosity as empty, whatever its
  // real saturation: drained means the fluid carries no load, and that is the
  // K in (1/K - 1/K_s).
  expect(rw.shearModulus > rd.shearModulus,
         "and stiffens shear too: the scheme does not satisfy Gassmann");
}

// --- Test group 5: degenerate inputs ---
void testDegenerate() {
  {
    const std::vector<homog::Phase> none;
    const homog::Result r = homog::selfConsistent(none);
    expectNear(r.bulkModulus, 0.0, 1e-12, "no phases: zero");
  }
  {
    std::vector<homog::Phase> allVoid{{1.0, 0.0, 0.0}};
    const homog::Result r = homog::selfConsistent(allVoid);
    expect(r.converged, "all void: converged");
    expectNear(r.bulkModulus, 0.0, 1e-12, "all void: zero stiffness");
  }
  {
    std::vector<homog::Phase> zeroFractions{{0.0, 30.0, 20.0}};
    const homog::Result r = homog::selfConsistent(zeroFractions);
    expectNear(r.bulkModulus, 0.0, 1e-12, "zero fractions: zero");
  }
}

} // namespace

int main() {
  std::printf("=== Homogenization: single phase ===\n");
  testSinglePhase();
  std::printf("\n=== Homogenization: porosity ===\n");
  testPorosity();
  std::printf("\n=== Homogenization: percolation limit ===\n");
  testPercolationLimit();
  std::printf("\n=== Homogenization: fluid-filled pores ===\n");
  testFluidFilledPore();
  std::printf("\n=== Homogenization: degenerate inputs ===\n");
  testDegenerate();

  std::printf("\n%s (%d failure%s)\n", gFailed == 0 ? "ALL PASS" : "FAILURES",
              gFailed, gFailed == 1 ? "" : "s");
  return (gFailed == 0) ? 0 : 1;
}
