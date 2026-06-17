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
SerreCyl, Consistency theorem (algebraic crux).

The velocity-potential expansion (eq:potapr)
    phi = phi0(x) + r^2 phi2 + r^4 phi4,
    phi2 = -(delta^2/4) phi0_xx,   phi4 = (delta^4/64) phi0_xxxx,
solves the cylindrical Laplace operator
    L[phi] = delta^2 r phi_xx + d_r( r phi_r )
identically through O(delta^4): the r^1 and r^3 coefficients vanish, leaving a
bulk residual exactly (delta^6/64) r^5 phi0^{(6)} = O(delta^6).  Together with the
axis condition phi_r(0)=0 (exact), this is the crux of the O(delta^4) consistency
of the cylindrical SGN system with the axisymmetric Euler equations.
"""
import sympy as sp

x, r, delta = sp.symbols('x r delta', positive=True)
phi0 = sp.Function('phi0')(x)
phi2 = -(delta**2/4)*sp.diff(phi0, x, 2)
phi4 = (delta**4/64)*sp.diff(phi0, x, 4)
phi = phi0 + r**2*phi2 + r**4*phi4

L = delta**2*r*sp.diff(phi, x, 2) + sp.diff(r*sp.diff(phi, r), r)
L = sp.expand(L)
print("L[phi] =", L)

c1 = sp.simplify(L.coeff(r, 1))
c3 = sp.simplify(L.coeff(r, 3))
print("  coeff(r^1) =", c1, "   coeff(r^3) =", c3)
assert c1 == 0 and c3 == 0, "Laplace not solved through O(delta^4)"

residual = sp.simplify(L)
expected = (delta**6/64)*r**5*sp.diff(phi0, x, 6)
print("  residual L[phi] =", residual)
assert sp.simplify(residual - expected) == 0
print("  PASS  r^1,r^3 coefficients vanish; residual = (delta^6/64) r^5 phi0^(6) = O(delta^6)")

# axis regularity:  phi_r = 2 r phi2 + 4 r^3 phi4  -> 0 at r=0
phir0 = sp.simplify(sp.diff(phi, r).subs(r, 0))
assert phir0 == 0
print("  PASS  phi_r(0) = 0 (axis condition satisfied exactly)")
print("=> Laplace + axis solved to O(delta^4); model residual O(delta^4) after averaging.")
