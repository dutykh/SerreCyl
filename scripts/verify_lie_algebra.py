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
SerreCyl, Lemma 6: Lie-algebra structure of the point symmetries in the maximal
(constant-radius, alpha=beta=r0'=0) case.  Generators
    X1 = d_x,  X2 = d_t,  X3 = t d_x + d_ubar,  X4 = ubar d_ubar - t d_t.
Verified nonzero commutators:  [X2,X3]=X1, [X2,X4]=-X2, [X3,X4]=X3 ; all others 0.
The subalgebra {X1,X2,X3} (with [X2,X3]=X1) is the (1+1)-D Galilei algebra.
"""
import sympy as sp

x, t, eta, ub = sp.symbols('x t eta ub')
f = sp.Function('f')(x, t, eta, ub)
VARS = (x, t, eta, ub)


def field(coeffs):
    """vector field as a first-order operator; coeffs: dict var->coefficient"""
    return lambda g: sum(coeffs.get(v, 0)*sp.diff(g, v) for v in VARS)


X1 = field({x: 1})
X2 = field({t: 1})
X3 = field({x: t, ub: 1})
X4 = field({ub: ub, t: -t})


def comm(A, B, g):
    return sp.expand(A(B(g)) - B(A(g)))


def check(A, B, target, name):
    lhs = comm(A, B, f)
    rhs = sp.expand(target(f))
    assert sp.simplify(lhs - rhs) == 0, (name, lhs, rhs)
    print("  PASS  %s" % name)


zero = lambda g: sp.Integer(0)
check(X2, X3, X1, "[X2,X3] = X1")
check(X2, X4, field({t: -1}), "[X2,X4] = -X2")
check(X3, X4, X3, "[X3,X4] = X3")
check(X1, X2, zero, "[X1,X2] = 0")
check(X1, X3, zero, "[X1,X3] = 0")
check(X1, X4, zero, "[X1,X4] = 0")
print("ALL COMMUTATORS VERIFIED  -> {X1,X2,X3} is the (1+1) Galilei algebra")
