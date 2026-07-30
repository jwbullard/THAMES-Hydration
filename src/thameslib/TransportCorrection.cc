/**
@file TransportCorrection.cc
@brief Phase-2 stub implementations. Phase 3 fills in Newton solver
       and per-shell D_eff mapping.
*/

#include "TransportCorrection.h"

#include <algorithm>
#include <cmath>

#include "TransportParameters.h"

namespace xport {

double kineticRate(double k, double area, double driving_force) {
  return k * area * driving_force;
}

double diffusionRate(double dEff, double area, double deltaC,
                     double delta) {
  if (delta <= 0.0) return 0.0;
  return dEff * area * deltaC / delta;
}

double solveSurfaceConcentrationLinear(double k, double C_eq, double dEff,
                                       double delta, double C_bulk) {
  // Solve  k · (1 - C_surf / C_eq)  =  dEff · (C_surf - C_bulk) / δ
  // Rearrange:  k - k · C_surf / C_eq  =  dEff · C_surf / δ  -  dEff · C_bulk / δ
  //             k + dEff · C_bulk / δ  =  C_surf · (k / C_eq + dEff / δ)
  //             C_surf  =  (k + dEff · C_bulk / δ)  /  (k / C_eq + dEff / δ)
  if (delta <= 0.0) return C_bulk;      // no shell, C_surf = C_bulk
  if (C_eq <= 0.0) return C_bulk;       // degenerate
  const double num = k + dEff * C_bulk / delta;
  const double den = k / C_eq + dEff / delta;
  if (den <= 0.0) return C_bulk;
  return num / den;
}

double pickDEff(int /*shellPhaseId*/, const TransportParameters &params) {
  // Phase 2 stub retained through Phase 3 first cut. A future
  // refinement will introduce a per-shell-phase map so that e.g.
  // Ca+2 through C-S-H and Ca+2 through AFm get different values.
  return params.dEff;
}

double solveSurfaceConcentration(double k, double C_eq, double dEff,
                                 double delta, double C_bulk,
                                 double (*f_driving)(double omega)) {
  // No-shell fallback (matches linear closed form for these edge cases).
  if (delta <= 0.0) return C_bulk;
  if (C_eq <= 0.0) return C_bulk;
  if (f_driving == nullptr) return C_bulk;

  // Steady-state residual: F(C_surf) = k · f(Ω_surf) - dEff · (C_surf-C_bulk)/δ.
  // At C_surf = C_eq: f = 0 and the diffusion side is dEff · (C_eq - C_bulk)/δ.
  // For dissolution (C_bulk < C_eq): F(C_bulk) = k · f(Ω_bulk) > 0,
  //   F(C_eq) = -dEff · (C_eq - C_bulk)/δ < 0 → sign change ⇒ bracketed.
  // For precipitation (C_bulk > C_eq): signs flip; bracket [C_eq, C_bulk].
  //
  // Fast exit at equilibrium: no flux, C_surf = C_bulk = C_eq.
  if (std::fabs(C_bulk - C_eq) <= 0.0) return C_bulk;

  const double a = std::min(C_bulk, C_eq);
  const double b = std::max(C_bulk, C_eq);

  auto residual = [&](double C_s) -> double {
    const double omega = C_s / C_eq;
    return k * f_driving(omega) - dEff * (C_s - C_bulk) / delta;
  };

  double fa = residual(a);
  double fb = residual(b);
  // If both endpoints share sign, the bracket is degenerate — fall
  // back to the bulk value rather than iterate to a garbage answer.
  if (fa * fb > 0.0) return C_bulk;
  if (std::fabs(fa) <= 0.0) return a;
  if (std::fabs(fb) <= 0.0) return b;

  // Standard Brent's method (Numerical Recipes / Wikipedia formulation).
  double aBr = a, bBr = b, cBr = a;
  double faBr = fa, fbBr = fb, fcBr = fa;
  double d = b - a, e = d;
  const int maxIter = 60;
  const double tolRel = 1.0e-8;

  for (int iter = 0; iter < maxIter; ++iter) {
    if (fbBr * fcBr > 0.0) {
      cBr = aBr;
      fcBr = faBr;
      d = bBr - aBr;
      e = d;
    }
    if (std::fabs(fcBr) < std::fabs(fbBr)) {
      aBr = bBr; bBr = cBr; cBr = aBr;
      faBr = fbBr; fbBr = fcBr; fcBr = faBr;
    }
    const double tol1 = 2.0 * 1e-15 * std::fabs(bBr) + 0.5 * tolRel * std::fabs(b - a);
    const double xm = 0.5 * (cBr - bBr);
    if (std::fabs(xm) <= tol1 || fbBr == 0.0) return bBr;

    if (std::fabs(e) >= tol1 && std::fabs(faBr) > std::fabs(fbBr)) {
      const double s = fbBr / faBr;
      double p, q;
      if (aBr == cBr) {
        p = 2.0 * xm * s;
        q = 1.0 - s;
      } else {
        const double qq = faBr / fcBr;
        const double rr = fbBr / fcBr;
        p = s * (2.0 * xm * qq * (qq - rr) - (bBr - aBr) * (rr - 1.0));
        q = (qq - 1.0) * (rr - 1.0) * (s - 1.0);
      }
      if (p > 0.0) q = -q;
      p = std::fabs(p);
      const double min1 = 3.0 * xm * q - std::fabs(tol1 * q);
      const double min2 = std::fabs(e * q);
      if (2.0 * p < std::min(min1, min2)) {
        e = d;
        d = p / q;
      } else {
        d = xm;
        e = d;
      }
    } else {
      d = xm;
      e = d;
    }
    aBr = bBr;
    faBr = fbBr;
    if (std::fabs(d) > tol1) {
      bBr += d;
    } else {
      bBr += (xm >= 0.0 ? tol1 : -tol1);
    }
    fbBr = residual(bBr);
  }
  return bBr;
}

}  // namespace xport
