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
SerreCyl, Lemma 10: dispersionless (delta -> 0) limit and hodograph linearization.

The constant-radius, gamma = 0, delta = 0 reduction of (eq:bouss1, eq:bouss2) is
    eta_t  + ubar eta_x + ((r0 + eta)/2) ubar_x = 0,
    ubar_t + ubar ubar_x + betat eta_x        = 0,
a 2x2 quasilinear system U_t + M(U) U_x = 0, U = (eta, ubar)^T.  We verify that
M has real, distinct eigenvalues lambda_pm = ubar +/- sqrt( betat (r0 + eta)/2 ),
hence the system is strictly hyperbolic; it therefore possesses Riemann invariants
and (away from simple waves) is linearized by the hodograph transformation, the
exact analogue of the classical shallow-water equations (Euler-Poisson-Darboux).
"""
import sympy as sp

eta, ubar, r0, betat = sp.symbols('eta ubar r0 betat', positive=True)

M = sp.Matrix([[ubar, (r0 + eta)/2],
               [betat, ubar]])
print("characteristic matrix M =", M.tolist())

evals = list(M.eigenvals().keys())
disc = betat*(r0 + eta)/2
expected = {sp.simplify(ubar + sp.sqrt(disc)), sp.simplify(ubar - sp.sqrt(disc))}
got = {sp.simplify(e) for e in evals}
print("eigenvalues:", got)
assert got == expected, (got, expected)

# real & distinct since betat, (r0+eta) > 0  =>  disc > 0
assert sp.simplify(M.eigenvals()[ubar + sp.sqrt(disc)]) == 1  # multiplicity 1
print("PASS: lambda_pm = ubar +/- sqrt(betat (r0+eta)/2), real and distinct")
print("      => strictly hyperbolic => Riemann invariants exist => hodograph linearization (EPD).")
