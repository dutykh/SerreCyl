#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Author: Dr. Denys Dutykh,
#         Khalifa University of Science and Technology, Abu Dhabi, UAE.
#
# Figure-generating / verification script accompanying the manuscript
# "Nonlinear Dynamics of Pulsatile Blood Flow in Viscoelastic Vessels".
#
"""
Modulational-instability comparison of the cylindrical KdV and BBM reductions
(dissipation-free case gamma = 0), using the physiological parameters of the
manuscript.  The NLS amplitude equations are

  KdV :  i A_{T2} - 3 k A_zz + beta A|A|^2 = 0,   beta = b^2/(6 c_K k),
         lambda^2 = 9 k^2 q^4 + 6 beta k a0^2 q^2 ;
  BBM :  A_tau + i M A_xx + i N A|A|^2 = 0,
         lambda^2 = M^2 q^4 - 2 N M a0^2 q^2 .

Modulational instability requires lambda^2 < 0 for some q, i.e. the Lighthill
products  PQ_KdV = (-3k)(beta) = -3 k beta  and  PQ_BBM = (-M)(-N) = M N  positive.

The script (i) builds the coefficients symbolically, (ii) evaluates them for the
physiological data and a range of carrier wavenumbers k, reporting the Lighthill
sign, and (iii) plots lambda^2(q) for KdV and BBM against a focusing reference.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- physiological parameters (manuscript, Sec. dispersion / Fig. 2) ----
r0 = 0.01            # m
h = 3.0e-4           # m
E = 4.1e5            # kg m^-1 s^-2
rho_w = 1.0e3        # kg m^-3
rho = 1060.0         # kg m^-3

ct = np.sqrt(E*h/(2*rho*r0))         # Moens-Korteweg speed
alphat = rho_w*h/rho                 # tilde-alpha
betat = E*h/(rho*r0**2)              # tilde-beta
print("ct=%.4f m/s, alphat=%.4e, betat=%.4e" % (ct, alphat, betat))

# KdV / BBM equation coefficients (gamma = 0)
a = ct
b = 5*ct/(2*r0)
cK = ct*(12*alphat + 3*r0)*r0/48     # KdV dispersion coeff (eq:kdv_burgers)
cB = (12*alphat + 3*r0)*r0/48        # BBM coeff (eq:BBM1) = cK/ct


def kdv_coeffs(k):
    beta = b**2/(6*cK*k)
    PQ = -3*k*beta                   # Lighthill product (MI iff > 0)
    return beta, PQ


def bbm_coeffs(k):
    w = a*k/(1 + cB*k**2)
    p = a*(1 - cB*k**2)/(1 + cB*k**2)**2
    sigma = b*(1 + cB*k**2)/(6*a*cB*k**2)
    nu = b*(1 + cB*k**2)**2/(-3*a*cB*k**2 - a*cB**2*k**4)
    M = (w*cB + 2*k*p*cB)/(1 + cB*k**2)
    N = b*k*(nu + sigma)/(1 + cB*k**2)
    PQ = M*N                         # Lighthill product (MI iff > 0)
    return M, N, PQ


# ---- scan carrier wavenumbers (long-wave: k r0 <~ 1) ----
ks = np.linspace(1.0, 100.0, 400)
PQ_kdv = np.array([kdv_coeffs(k)[1] for k in ks])
PQ_bbm = np.array([bbm_coeffs(k)[2] for k in ks])
print("max Lighthill product over k:  KdV=%.3e  BBM=%.3e" % (PQ_kdv.max(), PQ_bbm.max()))
kdv_unstable = np.any(PQ_kdv > 0)
bbm_unstable = np.any(PQ_bbm > 0)
print("KdV modulationally unstable for some k? ", kdv_unstable)
print("BBM modulationally unstable for some k? ", bbm_unstable)

# ---- lambda^2(q) at a representative carrier wavenumber ----
k0 = 20.0          # carrier (k0 r0 = 0.2, long-wave)
a0 = 0.05          # modulation amplitude
beta, _ = kdv_coeffs(k0)
M, N, _ = bbm_coeffs(k0)
q = np.linspace(0, 3.0, 400)
lam2_kdv = 9*k0**2*q**4 + 6*beta*k0*a0**2*q**2
lam2_bbm = M**2*q**4 - 2*N*M*a0**2*q**2
# focusing reference: i A_t + A_xx + |A|^2 A  -> lambda^2 = q^4 - 2 a0^2 q^2  (unstable band)
lam2_foc = q**4 - 2*a0**2*q**2
print("k0=%.1f: beta=%.4e (KdV), M=%.4e, N=%.4e (BBM)" % (k0, beta, M, N))

fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
ax[0].plot(ks*r0, np.sign(PQ_kdv), lw=2, label="KdV: sign$(-3k\\beta)$")
ax[0].plot(ks*r0, np.sign(PQ_bbm), '--', lw=2, label="BBM: sign$(MN)$")
ax[0].axhline(0, color='k', lw=0.6)
ax[0].set_xlabel(r"$k\,r_0$ (carrier)")
ax[0].set_ylabel("Lighthill product sign")
ax[0].set_ylim(-1.4, 1.4)
ax[0].set_title("Modulational (in)stability criterion")
ax[0].legend(loc="upper right", fontsize=9)

ax[1].plot(q, lam2_kdv/lam2_kdv.max(), lw=2, label="KdV (cylindrical)")
ax[1].plot(q, lam2_bbm/lam2_bbm.max(), '--', lw=2, label="BBM (cylindrical)")
ax[1].plot(q, lam2_foc/abs(lam2_foc).max(), ':', lw=2, label="focusing NLS (ref.)")
ax[1].axhline(0, color='k', lw=0.6)
ax[1].set_xlabel(r"sideband wavenumber $q$")
ax[1].set_ylabel(r"$\lambda^2$ (normalised)")
ax[1].set_title(r"$\lambda^2>0$: stable; $\lambda^2<0$: MI band")
ax[1].legend(loc="upper left", fontsize=9)

plt.tight_layout()
out = "figs/mi_growth.pdf"
plt.savefig(out)
print("wrote", out)
