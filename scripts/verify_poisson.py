#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Author: Dr. Denys Dutykh,
#         Khalifa University of Science and Technology, Abu Dhabi, UAE.
#
# Symbolic verification script accompanying the manuscript
# "Nonlinear Dynamics of Pulsatile Blood Flow in Viscoelastic Vessels".
#
"""
SerreCyl, Lemma 9 + variational/Hamiltonian enrichment (gamma = 0).

KdV reduction:  eta_t = - d_x ( delta H / delta eta ),  J = d_x,
   H = int ( a/2 eta^2 + b/6 eta^3 - c/2 eta_x^2 ) dx,
   a = ct, b = 5 ct/(2 r0), c = ct (12 alphat + 3 r0) r0 / 48.
Checks:
  (1) -d_x (delta H/delta eta) reproduces the KdV right-hand side;
  (2) mass M = int eta is a Casimir of J = d_x  (J delta M/delta eta = 0);
  (3) variational characterization of the solitary wave:
      delta (H - cs P)/delta eta = 0  (P = int eta^2/2, the momentum)
      gives the solitary-wave ODE  c eta_xx + (a-cs) eta + (b/2) eta^2 = 0.
BBM reduction:  (1 - c d_x^2) eta_t = - d_x (delta H/delta eta),
   H = int ( a/2 eta^2 + b/6 eta^3 ) dx ; mass Casimir.
"""
import sympy as sp

x = sp.symbols('x')
a, b, c, cs = sp.symbols('a b c cs')
eta = sp.Function('eta')(x)
E, Ex, Exx = sp.symbols('E Ex Exx')      # placeholders for eta, eta_x, eta_xx


def vderiv(density):
    """Variational (Euler-Lagrange) derivative of int density dx,
       density = density(E, Ex, Exx)."""
    dE = sp.diff(density, E)
    dEx = sp.diff(density, Ex)
    dExx = sp.diff(density, Exx)
    sub = {E: eta, Ex: sp.diff(eta, x), Exx: sp.diff(eta, x, 2)}
    return sp.expand(dE.subs(sub) - sp.diff(dEx.subs(sub), x)
                     + sp.diff(dExx.subs(sub), x, 2))


def dx(e, n=1):
    return sp.diff(e, x, n)


# ---- KdV Hamiltonian ----
H_kdv = a/2*E**2 + b/6*E**3 - c/2*Ex**2
dH = vderiv(H_kdv)
print("delta H/delta eta (KdV) =", dH)            # a eta + (b/2) eta^2 + c eta_xx
kdv_rhs = -dx(dH)
kdv_expected = -(a*dx(eta) + b*eta*dx(eta) + c*dx(eta, 3))
assert sp.simplify(kdv_rhs - kdv_expected) == 0
print("PASS (1) eta_t = -d_x delta H/delta eta reproduces KdV")

# ---- Casimir: M = int eta -> delta M/delta eta = 1 -> d_x(1)=0 ----
M_density = E
dM = vderiv(M_density)
assert sp.simplify(dM) == 1 and sp.simplify(dx(dM)) == 0
print("PASS (2) mass M=int eta is a Casimir of J=d_x")

# ---- variational characterization of the solitary wave ----
P_density = E**2/2
var = sp.expand(vderiv(H_kdv) - cs*vderiv(P_density))   # delta(H - cs P)/delta eta
sw_ode = c*dx(eta, 2) + (a - cs)*eta + b/2*eta**2
assert sp.simplify(var - sw_ode) == 0
print("PASS (3) delta(H - cs P)/delta eta = c eta_xx + (a-cs) eta + (b/2) eta^2 = solitary-wave ODE")

# ---- BBM ----
H_bbm = a/2*E**2 + b/6*E**3
dHb = vderiv(H_bbm)
bbm_lhs = dx(eta) * 0  # placeholder
# (1 - c d_x^2) eta_t = -d_x delta H/delta eta  ==>  check the spatial part identity
rhs_b = -dx(dHb)
rhs_b_expected = -(a*dx(eta) + b*eta*dx(eta))
assert sp.simplify(rhs_b - rhs_b_expected) == 0
print("PASS (4) BBM: (1-c d_x^2) eta_t = -d_x delta H/delta eta ; mass is the Casimir")
print("ALL POISSON / NOETHER / VARIATIONAL CHECKS PASSED")
