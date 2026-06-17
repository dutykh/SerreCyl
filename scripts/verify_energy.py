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
Symbolic verification for SerreCyl, Lemma 1 (energy balance) for the unidirectional
reductions, via explicit divergence identities.

KdV-Burgers (eq:kdv_burgers):  eta_t + a eta_x + b eta eta_x + c eta_xxx - nu eta_xx = 0,
  nu = r0 tilde-beta gamma / 4 >= 0.   Energy E = int (1/2) eta^2 dx.
BBM         (eq:BBM1):  eta_t + a eta_x + b eta eta_x - c eta_xxt - nu eta_xx = 0.
  Energy E = int (1/2)(eta^2 + c eta_x^2) dx.

We prove  d/dt(density) = -nu (eta_x)^2 + d/dx(flux),  so dE/dt = -nu int eta_x^2 <= 0,
and = 0 when gamma = 0 (nu = 0).
"""
import sympy as sp

x, t = sp.symbols('x t')
a, b, c, nu = sp.symbols('a b c nu')
eta = sp.Function('eta')(x, t)

def dx(e, n=1):
    return sp.diff(e, x, n)

print("=== Lemma 1: KdV-Burgers energy/L^2 balance ===")
# eta_t from the equation
eta_t_kdv = -(a*dx(eta) + b*eta*dx(eta) + c*dx(eta, 3)) + nu*dx(eta, 2)
dens_t = eta*eta_t_kdv                     # d/dt (1/2 eta^2) = eta*eta_t
# proposed flux F so that dens_t + nu*eta_x^2 = d/dx F
F = -a*eta**2/2 - b*eta**3/3 - c*eta*dx(eta, 2) + c*dx(eta)**2/2 + nu*eta*dx(eta)
resid = sp.simplify(dens_t + nu*dx(eta)**2 - dx(F))
print("  residual (eta_t-substituted) =", resid)
assert resid == 0
print("  PASS  d/dt(1/2 eta^2) = -nu eta_x^2 + d/dx(flux)  => dE/dt = -nu int eta_x^2 <= 0")

print()
print("=== Lemma 1: BBM energy/L^2 balance ===")
# BBM: (1 - c d_x^2) eta_t = -(a eta_x + b eta eta_x) + nu eta_xx
# Energy density e = 1/2 (eta^2 + c eta_x^2);  de/dt = eta eta_t + c eta_x eta_xt
# Use the weak identity: int de/dt = int eta_t (eta - c eta_xx) = int eta_t (1-c d_x^2) eta
# We verify the *pointwise* divergence form instead, eliminating eta_t via the operator.
etat = sp.Function('w')(x, t)             # stands for eta_t
rhs = -(a*dx(eta) + b*eta*dx(eta)) + nu*dx(eta, 2)
# constraint:  etat - c d_x^2 etat = rhs
constraint = etat - c*dx(etat, 2) - rhs
# energy-density time derivative:
de_dt = eta*etat + c*dx(eta)*dx(etat)
# claim: de_dt = -nu eta_x^2 + d/dx(Flux) + (multiple of constraint)
# Equivalent, integrated: int de_dt = int eta*(etat - c etat_xx) = int eta*rhs = -nu int eta_x^2
# Verify by parts that de_dt - eta*(etat - c*dx(etat,2)) is a pure x-divergence:
G = c*eta*dx(etat) - c*dx(eta)*etat + c*dx(eta)*etat  # = c eta etat_x  (book-keeping)
check = sp.simplify(de_dt - (eta*etat - c*eta*dx(etat, 2)) - dx(c*eta*dx(etat)))
print("  de/dt - eta*(1-c d_x^2)etat - d/dx(c eta etat_x) =", check)
assert check == 0
# now eta*(1-c d_x^2)etat = eta*rhs ; show eta*rhs = -nu eta_x^2 + d/dx(flux2)
flux2 = -a*eta**2/2 - b*eta**3/3 + nu*eta*dx(eta)
resid2 = sp.simplify(eta*rhs + nu*dx(eta)**2 - dx(flux2))
print("  eta*rhs + nu eta_x^2 - d/dx(flux2) =", resid2)
assert resid2 == 0
print("  PASS  dE/dt = -nu int eta_x^2 <= 0 for BBM as well (=0 when gamma=0)")
print()
print("ALL ENERGY-BALANCE CHECKS PASSED")
