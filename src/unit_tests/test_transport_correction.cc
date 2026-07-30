// Standalone unit test for TransportCorrection Phase-2 stubs
// (kineticRate, diffusionRate, solveSurfaceConcentrationLinear, pickDEff).
// Compile with build_and_run.sh — no THAMES dependencies needed.

#include "TransportCorrection.h"

#include <cmath>
#include <cstdio>

namespace {

int gFailed = 0;

void expectNear(double actual, double expected, double tol,
                const char *tag) {
  const bool ok = std::fabs(actual - expected) <= tol;
  std::printf("%-55s actual=%.6e expected=%.6e tol=%.2e %s\n",
              tag, actual, expected, tol, ok ? "PASS" : "FAIL");
  if (!ok) ++gFailed;
}

void expect(bool cond, const char *tag) {
  std::printf("%-55s %s\n", tag, cond ? "PASS" : "FAIL");
  if (!cond) ++gFailed;
}

// --- kineticRate: r = k·A·f ---
void testKineticRate() {
  std::printf("\n--- kineticRate ---\n");
  expectNear(xport::kineticRate(1e-6, 1.0, 0.5), 5e-7, 1e-18,
             "1e-6 · 1 · 0.5 = 5e-7");
  expectNear(xport::kineticRate(0.0, 100.0, 0.5), 0.0, 1e-30,
             "k=0 → 0");
  expectNear(xport::kineticRate(1e-6, 0.0, 0.5), 0.0, 1e-30,
             "A=0 → 0");
}

// --- diffusionRate: r = D·A·ΔC/δ ---
void testDiffusionRate() {
  std::printf("\n--- diffusionRate ---\n");
  expectNear(xport::diffusionRate(1e-13, 1.0, 10.0, 1e-6),
             1e-13 * 10.0 / 1e-6, 1e-18,
             "1e-13 · 1 · 10 / 1e-6 = 1e-6");
  expectNear(xport::diffusionRate(1e-13, 1.0, 10.0, 0.0), 0.0, 1e-30,
             "δ=0 guard → 0 (no divide-by-zero)");
  expectNear(xport::diffusionRate(0.0, 1.0, 10.0, 1e-6), 0.0, 1e-30,
             "D=0 → 0");
}

// --- solveSurfaceConcentrationLinear ---
void testSolveLinearAsymptotes() {
  std::printf("\n--- solveSurfaceConcentrationLinear asymptotes ---\n");
  const double C_eq = 20.0;     // mol/kg water (arbitrary units)
  const double C_bulk = 5.0;    // undersaturated (dissolution regime)

  // δ → 0 (no shell): C_surf → C_bulk (surface sees bulk directly).
  const double kInf = 1e-6;
  const double DInf = 1e-13;
  const double C0 =
      xport::solveSurfaceConcentrationLinear(kInf, C_eq, DInf, 0.0, C_bulk);
  expectNear(C0, C_bulk, 1e-18,
             "δ=0 → C_surf = C_bulk (guarded fallback)");

  // Kinetic-limited regime (k << D/δ · C_eq): C_surf → C_bulk still.
  const double C_kinLim =
      xport::solveSurfaceConcentrationLinear(1e-20, C_eq, 1.0, 1e-6, C_bulk);
  expectNear(C_kinLim, C_bulk, 1e-6,
             "kinetic limit (k → 0) → C_surf ≈ C_bulk");

  // Diffusion-limited regime (D << k · δ): C_surf → C_eq (equilibrium
  // at surface because slow diffusion can't clear the produced Ca).
  const double C_diffLim =
      xport::solveSurfaceConcentrationLinear(1.0, C_eq, 1e-30, 1e-6, C_bulk);
  expectNear(C_diffLim, C_eq, 1e-6,
             "diffusion limit (D → 0) → C_surf ≈ C_eq");
}

// --- solveSurfaceConcentrationLinear: check the closed-form algebra ---
void testSolveLinearClosedForm() {
  std::printf("\n--- solveSurfaceConcentrationLinear closed form ---\n");
  // With k=1, C_eq=10, D=2, δ=1, C_bulk=3:
  // C_surf = (1 + 2·3/1) / (1/10 + 2/1) = 7 / 2.1 = 10/3
  const double C_surf = xport::solveSurfaceConcentrationLinear(
      1.0, 10.0, 2.0, 1.0, 3.0);
  expectNear(C_surf, 10.0 / 3.0, 1e-12,
             "closed-form check: k=1, C_eq=10, D=2, δ=1, C_bulk=3");
}

// --- pickDEff (Phase 2 stub: always returns block's global dEff) ---
void testPickDEffStub() {
  std::printf("\n--- pickDEff (Phase 2 stub) ---\n");
  TransportParameters params;
  params.dEff = 3.14e-13;
  expectNear(xport::pickDEff(42, params), 3.14e-13, 1e-30,
             "Phase-2 stub returns block dEff regardless of shell id");
  expectNear(xport::pickDEff(-1, params), 3.14e-13, 1e-30,
             "even for unmapped shell id -1");
}

}  // namespace

int main() {
  std::printf("=== test_transport_correction ===\n");

  testKineticRate();
  testDiffusionRate();
  testSolveLinearAsymptotes();
  testSolveLinearClosedForm();
  testPickDEffStub();

  std::printf("\n%d probe(s) failed\n", gFailed);
  return gFailed == 0 ? 0 : 1;
}
