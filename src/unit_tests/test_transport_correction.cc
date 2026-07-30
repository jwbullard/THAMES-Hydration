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

// --- solveSurfaceConcentration (general Newton/Brent solver) ---

// Linear driving force f(Ω) = 1 - Ω (dissolution regime).
double fLinear(double omega) { return 1.0 - omega; }

// Standard's driving force with p = 2, q = 3 (nonlinear test case).
double fStandardP2Q3(double omega) {
  const double base = 1.0 - std::pow(omega, 2.0);
  if (base <= 0.0) return 0.0;
  return std::pow(base, 3.0);
}

void testGeneralSolverLinearAgreement() {
  std::printf("\n--- general solver vs closed form (linear) ---\n");
  const double k = 1.0, C_eq = 10.0, D = 2.0, delta = 1.0, C_bulk = 3.0;
  const double C_linear =
      xport::solveSurfaceConcentrationLinear(k, C_eq, D, delta, C_bulk);
  const double C_general =
      xport::solveSurfaceConcentration(k, C_eq, D, delta, C_bulk, fLinear);
  expectNear(C_general, C_linear, 1e-8,
             "general solver matches closed form for linear f");
  expectNear(C_general, 10.0 / 3.0, 1e-8,
             "general solver hits algebraic answer");
}

void testGeneralSolverNonlinear() {
  std::printf("\n--- general solver with p=2, q=3 kinetic law ---\n");
  const double k = 1e-4, C_eq = 20.0, D = 1e-6, delta = 1e-4,
               C_bulk = 5.0;
  const double C_surf = xport::solveSurfaceConcentration(
      k, C_eq, D, delta, C_bulk, fStandardP2Q3);
  // Sanity: C_surf must be strictly between C_bulk and C_eq
  // (dissolution regime, monotonic driving force).
  expect(C_surf > C_bulk && C_surf < C_eq,
         "C_surf strictly bracketed by (C_bulk, C_eq)");
  // Verify the residual is near zero at the returned root.
  const double omega = C_surf / C_eq;
  const double r_kin = k * fStandardP2Q3(omega);
  const double r_diff = D * (C_surf - C_bulk) / delta;
  expectNear(r_kin, r_diff, 1e-10,
             "flux balance holds at returned C_surf");
}

void testShellCorrectionFactor() {
  std::printf("\n--- shellCorrectionFactor bin-weighted ---\n");

  TransportParameters params;
  params.dEff = 1.0;  // pick a unit dEff for hand-math

  // Case 1: empty stats → factor = 1 (no correction).
  {
    xport::ShellStats empty;
    empty.numSitesReached = 0;
    empty.bins.emplace_back();
    empty.bins.back().siteFraction = 0.0;
    const double f = xport::shellCorrectionFactor(1.0, 10.0, empty, params);
    expectNear(f, 1.0, 1e-12, "empty stats → factor = 1");
  }

  // Case 2: single-bin monodisperse δ = 0 (no shell yet).
  //   All sites directly touching electrolyte → factor = 1.
  {
    xport::ShellStats s;
    s.numSitesReached = 100;
    xport::ShellBin bin;
    bin.deltaRep = 0.0;
    bin.siteFraction = 1.0;
    bin.dominantShellPhaseId = -1;
    s.bins.push_back(bin);
    const double f = xport::shellCorrectionFactor(1.0, 10.0, s, params);
    expectNear(f, 1.0, 1e-12, "δ=0 bin → factor = 1");
  }

  // Case 3: single-bin monodisperse δ > 0.
  //   Da = k · δ / (D · C_eq). With k=1, δ=1, D=1, C_eq=10 → Da = 0.1.
  //   Factor = 1 / (1 + 0.1) = 0.9090909...
  {
    xport::ShellStats s;
    s.numSitesReached = 100;
    xport::ShellBin bin;
    bin.deltaRep = 1.0;
    bin.siteFraction = 1.0;
    bin.dominantShellPhaseId = 42;
    s.bins.push_back(bin);
    const double f = xport::shellCorrectionFactor(1.0, 10.0, s, params);
    expectNear(f, 1.0 / 1.1, 1e-12,
               "single-bin Da=0.1 → factor = 1/1.1");
  }

  // Case 4: bimodal 50/50 with δ_thin = 0.1, δ_thick = 10.
  //   Da_thin = 1 · 0.1 / (1 · 10) = 0.01 → f_thin = 1/1.01 ≈ 0.990099
  //   Da_thick = 1 · 10  / (1 · 10) = 1.0  → f_thick = 1/2   = 0.5
  //   Weighted avg = 0.5 · 0.990099 + 0.5 · 0.5 = 0.7450495
  {
    xport::ShellStats s;
    s.numSitesReached = 200;
    xport::ShellBin thin;
    thin.deltaRep = 0.1;
    thin.siteFraction = 0.5;
    thin.dominantShellPhaseId = 1;
    xport::ShellBin thick;
    thick.deltaRep = 10.0;
    thick.siteFraction = 0.5;
    thick.dominantShellPhaseId = 2;
    s.bins.push_back(thin);
    s.bins.push_back(thick);
    const double f = xport::shellCorrectionFactor(1.0, 10.0, s, params);
    const double expected = 0.5 * (1.0 / 1.01) + 0.5 * (1.0 / 2.0);
    expectNear(f, expected, 1e-12,
               "bimodal weighted correction matches hand-math");
  }

  // Case 5: guard — C_eq = 0 → factor = 1 (no correction, defensive).
  {
    xport::ShellStats s;
    s.numSitesReached = 100;
    xport::ShellBin bin;
    bin.deltaRep = 1.0;
    bin.siteFraction = 1.0;
    s.bins.push_back(bin);
    const double f = xport::shellCorrectionFactor(1.0, 0.0, s, params);
    expectNear(f, 1.0, 1e-12, "C_eq = 0 → factor = 1 (guard)");
  }
}

void testGeneralSolverEdgeCases() {
  std::printf("\n--- general solver edge cases ---\n");
  // Equilibrium: C_bulk = C_eq → no flux, C_surf = C_bulk.
  const double C_eq = xport::solveSurfaceConcentration(
      1.0, 10.0, 1e-6, 1e-6, 10.0, fLinear);
  expectNear(C_eq, 10.0, 1e-12, "C_bulk = C_eq → C_surf = C_bulk");
  // delta = 0 → no-shell fallback returns C_bulk.
  const double Cz = xport::solveSurfaceConcentration(
      1.0, 10.0, 1e-6, 0.0, 5.0, fLinear);
  expectNear(Cz, 5.0, 1e-12, "δ = 0 → C_surf = C_bulk");
  // Null functor → fallback returns C_bulk.
  const double Cn = xport::solveSurfaceConcentration(
      1.0, 10.0, 1e-6, 1e-6, 5.0, nullptr);
  expectNear(Cn, 5.0, 1e-12, "null functor → C_surf = C_bulk");
}

}  // namespace

int main() {
  std::printf("=== test_transport_correction ===\n");

  testKineticRate();
  testDiffusionRate();
  testSolveLinearAsymptotes();
  testSolveLinearClosedForm();
  testPickDEffStub();
  testGeneralSolverLinearAgreement();
  testGeneralSolverNonlinear();
  testGeneralSolverEdgeCases();
  testShellCorrectionFactor();

  std::printf("\n%d probe(s) failed\n", gFailed);
  return gFailed == 0 ? 0 : 1;
}
