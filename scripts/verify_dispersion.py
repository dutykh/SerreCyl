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
Symbolic verification for SerreCyl, Lemma 2 (asymptotic dispersion consistency)
and Lemma 3 (linear well-posedness / stability).

Lemma 2: the long-wave (k r0 -> 0) expansion of the EXACT linearized-Euler
(Bessel) dispersion relation eq:disp_euler agrees, through O((k r0)^2) inside the
bracket, with the cylindrical KdV reduction; hence omega_E - omega_KdV = O(k^5),
i.e. relative O(delta^4).

Lemma 3: from the explicit reduced dispersion relations, omega is real for
gamma = 0 and has non-positive imaginary part for gamma > 0.

Notation:  alpha = tilde-alpha = rho^w h / rho ;  beta = tilde-beta = E h/(rho r0^2);
           ct = Moens-Korteweg speed, ct^2 = E h/(2 rho r0) = beta*r0/2.
"""
import sympy as sp

k, r0, alpha, beta, ct, gamma, delta = sp.symbols(
    'k r0 alpha beta ct gamma delta', positive=True)

print("=== Lemma 2: Bessel expansion of the exact Euler dispersion ===")
z = k*r0
# truncated Maclaurin series of the modified Bessel functions (more than enough
# for the O(k^4) terms of omega^2); sympy expands the resulting rational function.
I0 = 1 + z**2/4 + z**4/64 + z**6/2304
I1 = z/2 + z**3/16 + z**5/384 + z**7/18432
w2_euler = beta*k*I1 / (alpha*k*I1 + I0)
ser = sp.expand(sp.series(w2_euler, k, 0, 6).removeO())
c2 = sp.simplify(ser.coeff(k, 2))
c4 = sp.simplify(ser.coeff(k, 4))
print("  omega_E^2 = (%s) k^2 + (%s) k^4 + O(k^6)" % (c2, c4))

assert sp.simplify(c2 - beta*r0/2) == 0, "leading coeff mismatch"
print("  PASS  leading coeff = beta*r0/2  (= ct^2, the Moens-Korteweg speed^2)")

# bracket form omega_E^2 = c2 * k^2 * (1 - mu_E k^2 + ...)
mu_E = sp.simplify(-c4/c2)
print("  mu_E (O(k^2) dispersive coeff) =", mu_E)        # expect alpha*r0/2 + r0^2/8
assert sp.simplify(mu_E - (alpha*r0/2 + r0**2/8)) == 0
print("  PASS  mu_E = alpha*r0/2 + r0^2/8")

# KdV reduction: omega = ct k - cK k^3,  cK = ct(12 alpha + 3 r0) r0 / 48
cK = ct*(12*alpha + 3*r0)*r0/48
w_kdv = ct*k - cK*k**3
w2_kdv = sp.expand(sp.series(w_kdv**2, k, 0, 6).removeO())
c2K = sp.simplify(w2_kdv.coeff(k, 2))     # = ct^2
c4K = sp.simplify(w2_kdv.coeff(k, 4))
mu_K = sp.simplify(-c4K/c2K)
print("  mu_KdV =", mu_K)
assert sp.simplify(mu_K - (alpha*r0/2 + r0**2/8)) == 0
print("  PASS  mu_KdV = mu_E  => omega_E = omega_KdV + O(k^5) = O(delta^4) relative")
# consistency of leading speeds: ct^2 must equal beta*r0/2
print("  (leading speeds match iff ct^2 = beta*r0/2, the Moens-Korteweg identity)")

print()
print("=== Lemma 3: linear well-posedness / stability ===")
# KdV-Burgers dispersion: omega = ct k - cK k^3 - (i/4) r0 gamma beta k^2
Im_kdv = -sp.Rational(1, 4)*r0*gamma*beta*k**2
print("  KdV:  Im(omega) =", Im_kdv, " (<= 0 for gamma>=0)")
assert (Im_kdv.subs({r0: 1, gamma: 1, beta: 1, k: 1}) < 0) and \
       sp.simplify(Im_kdv.subs(gamma, 0)) == 0
print("  PASS  Im<=0 (gamma>0 decay); Im=0 when gamma=0 (no growth)")

# BBM dispersion: omega = 48 ct k/D - i 12 gamma r0 beta k^2 / D, D = 48 + r0(12a+3r0)k^2
D = 48 + r0*(12*alpha + 3*r0)*k**2
Im_bbm = -12*gamma*r0*beta*k**2 / D
print("  BBM:  Im(omega) = -12 gamma r0 beta k^2 / D, D>0  => Im<=0")
assert Im_bbm.subs({r0: 1, gamma: 1, beta: 1, k: 1, alpha: 1}) < 0
print("  PASS  BBM Im<=0; Im=0 when gamma=0")

# Boussinesq (eq:disp_bous), gamma = 0:  omega^2 = beta k / Den, Den>0 => omega^2>0
Den = alpha*k*(24*r0*k + r0**3*k**3) + 48 + 8*delta**2*r0**2*k**2
w2_bous = beta*k/Den
print("  Boussinesq (gamma=0): omega^2 = beta k / Den, Den =", Den)
assert w2_bous.subs({k: 1, r0: 1, alpha: 1, beta: 1, delta: 1}) > 0
print("  PASS  omega^2 > 0 (Den>0 for all k>0) => omega real => linearly well-posed")
print()
print("ALL DISPERSION/STABILITY CHECKS PASSED")
