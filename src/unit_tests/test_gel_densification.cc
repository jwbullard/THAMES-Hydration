// Standalone unit test for GelDensificationParameters.h: the Konigsberger
// target density and the gel-envelope rule. Header only; no THAMES
// dependencies. Compile with build_and_run.sh.
//
// Expected values are computed by hand from Konigsberger et al., CCR 88
// (2016) Eq. 42 and Eq. 13 with the default coefficients, not by calling the
// code under test.

#include "GelDensificationParameters.h"

#include <cmath>
#include <cstdio>

namespace {

int gFailed = 0;

void expectNear(double actual, double expected, double tol, const char *tag) {
  const bool ok = std::fabs(actual - expected) <= tol;
  std::printf("%-60s actual=%.6f expected=%.6f %s\n", tag, actual, expected,
              ok ? "PASS" : "FAIL");
  if (!ok)
    ++gFailed;
}

const GelDensificationParameters P; // defaults
const double RHO_W = 1.0;
const double RHO_S = 2.25; // GEMS CSHQ solid, as in the cem151 runs
const double TOL = 1.0e-6;

// --- target density ---
void testTarget() {
  std::printf("\n--- gelDensityTarget ---\n");
  expectNear(gelDensityTarget(0.674, P, RHO_W), 1.624860, TOL,
             "regime II, psi 0.674 (Eq. 42)");
  expectNear(gelDensityTarget(0.942, P, RHO_W), 1.338034, TOL,
             "regime I/II boundary, psi 0.942");
  expectNear(gelDensityTarget(1.0, P, RHO_W), 1.338034, TOL,
             "above psi_I-II the onset value is held");
  expectNear(gelDensityTarget(0.426, P, RHO_W), 1.890280, TOL,
             "psi_II-III itself is regime II");
  expectNear(gelDensityTarget(0.4259999, P, RHO_W), 1.920696, 1.0e-5,
             "just below psi_II-III is regime III (Eq. 13)");
  expectNear(gelDensityTarget(0.0, P, RHO_W), 2.604, TOL,
             "psi 0: no water left, Allen solid");
}

// --- conversion to phi on the GEMS basis, no history ---
void testNoHistory() {
  std::printf("\n--- densifiedGelPorosity, first step (no history) ---\n");
  expectNear(densifiedGelPorosity(0.9, RHO_S, RHO_W, 1.0, 0.0, 0.0, P),
             0.693612, TOL, "psi 0.9 -> (2.25 - 1.382984) / 1.25");
  expectNear(densifiedGelPorosity(0.3, RHO_S, RHO_W, 1.0, 0.0, 0.0, P),
             0.101760, TOL, "psi 0.3 (regime III)");
  expectNear(densifiedGelPorosity(0.1, RHO_S, RHO_W, 1.0, 0.0, 0.0, P), 0.0,
             TOL, "target denser than GEMS solid -> floor 0");
}

// --- envelope rule, solid growing ---
void testGrowing() {
  std::printf("\n--- densifiedGelPorosity, solid growing ---\n");
  // Committed: solid 1.0 in envelope 2.0 (phi 0.5). Target phi 0.10176.
  const double phi =
      densifiedGelPorosity(0.3, RHO_S, RHO_W, 1.2, 1.0, 2.0, P);
  expectNear(phi, 0.4, TOL, "small growth: envelope held, phi 1 - 1.2/2");
  expectNear(1.2 / (1.0 - phi), 2.0, TOL, "  ... envelope unchanged at 2.0");

  // Large growth: the target is now above the hold, so the envelope grows.
  const double phi2 =
      densifiedGelPorosity(0.3, RHO_S, RHO_W, 3.0, 1.0, 2.0, P);
  expectNear(phi2, 0.101760, TOL, "large growth: target phi wins");
  const bool grew = 3.0 / (1.0 - phi2) > 2.0;
  std::printf("%-60s %s\n", "  ... envelope grew past 2.0",
              grew ? "PASS" : "FAIL");
  if (!grew)
    ++gFailed;
}

// --- envelope rule, solid shrinking ---
void testShrinking() {
  std::printf("\n--- densifiedGelPorosity, solid shrinking ---\n");
  // Committed: solid 1.0 in envelope 2.0 (phi 0.5). Solid falls to 0.8.
  // Bounds: proportional shrink phi 0.5 (envelope 1.6); same envelope
  // phi 1 - 0.8/2 = 0.6 (envelope 2.0).

  const double low =
      densifiedGelPorosity(0.3, RHO_S, RHO_W, 0.8, 1.0, 2.0, P);
  expectNear(low, 0.5, TOL, "target denser: dissolves at current porosity");
  expectNear(0.8 / (1.0 - low), 1.6, TOL,
             "  ... envelope shrinks in proportion (2.0 x 0.8)");

  const double high =
      densifiedGelPorosity(0.9, RHO_S, RHO_W, 0.8, 1.0, 2.0, P);
  expectNear(high, 0.6, TOL,
             "target more porous: capped, envelope may not grow");
  expectNear(0.8 / (1.0 - high), 2.0, TOL, "  ... envelope held at 2.0");

  // A target between the bounds is used as is: build one at phi 0.55.
  // phi = (2.25 - rho)/1.25 = 0.55 -> rho = 1.5625, which in regime II
  // (rho = (0.901 - 0.411 psi) 2.604) needs psi = 0.732267.
  const double mid =
      densifiedGelPorosity(0.732267, RHO_S, RHO_W, 0.8, 1.0, 2.0, P);
  expectNear(mid, 0.55, 1.0e-5, "target between the bounds is used");
}

// --- upper clamp ---
void testClamp() {
  std::printf("\n--- densifiedGelPorosity, clamp ---\n");
  expectNear(densifiedGelPorosity(0.3, RHO_S, RHO_W, 1.0, 1.0, 1000.0, P),
             GEL_POROSITY_MAX, TOL, "envelope hold above 0.99 is clamped");
}

} // namespace

int main() {
  testTarget();
  testNoHistory();
  testGrowing();
  testShrinking();
  testClamp();
  std::printf("\n%s: %d failure(s)\n", gFailed ? "FAILED" : "ALL PASSED",
              gFailed);
  return gFailed ? 1 : 0;
}
