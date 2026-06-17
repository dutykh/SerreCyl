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
SerreCyl, Lemma 7: symmetry reductions of the constant-radius Boussinesq system
(eq:bousnew1, eq:bousnew2), gamma = 0, in physical variables (eps=delta=1):

  (M)  eta_t + (r0/2) ubar_x + (1/2) eta ubar_x + eta_x ubar - (r0^3/48) ubar_xxx = 0,
  (P)  ubar_t + ubar ubar_x + beta eta_x - D ubar_xxt = 0,   D = alphat r0/2 + r0^2/6.

(a) Travelling-wave invariance under X2 + c X1: eta=H(x-ct), ubar=U(x-ct) reduces
    (M),(P) to ODEs in zeta=x-ct (every term is a function of zeta only).
(b) Self-similar invariance under the scaling generator X = ubar d_ubar - t d_t
    (admitted when beta = 0): invariants give  ubar = G(x)/t,  eta = H(x);
    substitution yields a CONSISTENT (t-independent) ODE system.
"""
import sympy as sp

x, t, c = sp.symbols('x t c')
r0, alphat, beta = sp.symbols('r0 alphat beta', positive=True)
D = alphat*r0/2 + r0**2/6


def dx(e, n=1):
    return sp.diff(e, x, n)


def dt(e, n=1):
    return sp.diff(e, t, n)


# ---------- (a) travelling-wave reduction ----------
zeta = x - c*t
H = sp.Function('H')
U = sp.Function('U')
eta_tw = H(zeta)
u_tw = U(zeta)
M_tw = dt(eta_tw) + (r0/2)*dx(u_tw) + sp.Rational(1, 2)*eta_tw*dx(u_tw) \
    + dx(eta_tw)*u_tw - (r0**3/48)*dx(u_tw, 3)
P_tw = dt(u_tw) + u_tw*dx(u_tw) + beta*dx(eta_tw) - D*dt(dx(u_tw, 2))
# every derivative is wrt the single argument (x-ct): the PDEs are ODEs in zeta.
s = sp.symbols('s')
M_tw_s = M_tw.subs(x, s + c*t)          # rewrite in terms of zeta=s
P_tw_s = P_tw.subs(x, s + c*t)
assert sp.simplify(sp.diff(M_tw_s, t)) == 0
assert sp.simplify(sp.diff(P_tw_s, t)) == 0
print("PASS (a): travelling-wave ansatz reduces (M),(P) to ODEs in zeta=x-ct")

# ---------- (b) self-similar reduction (beta = 0) ----------
G = sp.Function('G')
Hh = sp.Function('H')
u_ss = G(x)/t
eta_ss = Hh(x)
M_ss = dt(eta_ss) + (r0/2)*dx(u_ss) + sp.Rational(1, 2)*eta_ss*dx(u_ss) \
    + dx(eta_ss)*u_ss - (r0**3/48)*dx(u_ss, 3)
P_ss = dt(u_ss) + u_ss*dx(u_ss) - D*dt(dx(u_ss, 2))   # beta = 0
M_red = sp.simplify(M_ss*t)
P_red = sp.simplify(P_ss*t**2)
print("self-similar mass ODE     :", sp.expand(M_red))
print("self-similar momentum ODE :", sp.expand(P_red))
assert sp.simplify(sp.diff(M_red, t)) == 0
assert sp.simplify(sp.diff(P_red, t)) == 0
print("PASS (b): u=G(x)/t, eta=H(x) gives t-independent ODEs (a genuine self-similar reduction)")
