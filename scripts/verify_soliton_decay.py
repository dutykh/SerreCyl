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
SerreCyl, Lemma 5 (adiabatic viscoelastic attenuation of the KdV solitary wave).

The KdV-Burgers reduction is  eta_t + a eta_x + b eta eta_x + c eta_xxx - nu eta_xx = 0
with nu = r0*betat*gamma/4.  The L^2 (momentum) functional P = int (1/2) eta^2
obeys  dP/dt = -nu int eta_x^2.  Inserting the slowly-varying soliton
eta = A sech^2(K zeta) with the KdV width-amplitude relation
K^2 = 10 A / (r0^2 (12 alphat + 3 r0)) yields the amplitude law

        dA/dt = - [ 8 betat gamma / (3 r0 (12 alphat + 3 r0)) ] A^2 .
"""
import sympy as sp

u = sp.symbols('u', real=True)
A, K, r0, alphat, betat, gamma, Adot = sp.symbols(
    'A K r0 alphat betat gamma Adot', positive=True)

# canonical sech integrals, via the substitution s = tanh(u), ds = sech^2(u) du,
# sech^2 = 1 - s^2,  s in (-1, 1):
s = sp.symbols('s')
I2 = sp.integrate(sp.Integer(1), (s, -1, 1))            # int sech^2 du
I4 = sp.integrate(1 - s**2, (s, -1, 1))                 # int sech^4 du
I4t2 = sp.integrate((1 - s**2)*s**2, (s, -1, 1))        # int sech^4 tanh^2 du
print("int sech^2 =", I2, " int sech^4 =", I4, " int sech^4 tanh^2 =", I4t2)
assert (I2, I4, I4t2) == (2, sp.Rational(4, 3), sp.Rational(4, 15))

# eta = A sech^2(K zeta):  int eta^2 = 4A^2/(3K),  int eta_x^2 = 16 A^2 K /15
int_eta2 = A**2/K*I4
int_etax2 = 4*A**2*K**2/K*I4t2
assert sp.simplify(int_eta2 - sp.Rational(4, 3)*A**2/K) == 0
assert sp.simplify(int_etax2 - sp.Rational(16, 15)*A**2*K) == 0
print("int eta^2 =", sp.simplify(int_eta2), "  int eta_x^2 =", sp.simplify(int_etax2))

# width-amplitude relation and dissipation coefficient
K2 = 10*A/(r0**2*(12*alphat + 3*r0))
nu = r0*betat*gamma/4
Ks = sp.sqrt(K2)

# P = (1/2) int eta^2 ; dP/dt = (dP/dA) Adot = -nu int eta_x^2  (linear in Adot)
P = sp.Rational(1, 2)*int_eta2.subs(K, Ks)
dPdA = sp.diff(P, A)
rhs = -nu*int_etax2.subs(K, Ks)
sol = sp.simplify(rhs/dPdA)
expected = -sp.Rational(8, 3)*betat*gamma*A**2/(r0*(12*alphat + 3*r0))
print("dA/dt =", sol)
print("expected =", sp.simplify(expected))
assert sp.simplify(sol - expected) == 0
print("PASS: dA/dt = -[8 betat gamma /(3 r0 (12 alphat+3 r0))] A^2  (algebraic, A(t)=A0/(1+C gamma A0 t))")
