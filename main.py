#!/usr/bin/env python3
"""
alpha-free vs. traditional LLG with progress bar
------------------------------------------------
Requires: numpy, scipy, matplotlib, tqdm
Run: python llg_alpha_free_progress.py
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from tqdm import tqdm

# --------------------------------------------------------------
# 0. Physical / material constants (SI units)
# --------------------------------------------------------------
mu0   = 4 * np.pi * 1e-7         # vacuum permeability
gamma = 2.211e5 / mu0            # gyromagnetic ratio (rad·s⁻¹·T⁻¹)
Ms    = 8.0e5                    # saturation magnetisation (A/m)
Kuni  = 5.0e4                    # uniaxial anisotropy (J/m³)
Hk    = 2*Kuni/(mu0*Ms)          # anisotropy field (A/m)

# Susceptibility tensor χ (transverse part only, complex for damping)
chi0, chi1 = 1.5, 0.05
chi_perp   = chi0 - 1j*chi1
# easy axis along z => χ = χ⊥(I - zzᵀ)
chi = chi_perp * (np.eye(3) - np.outer([0,0,1],[0,0,1]))
eps = 1e-8
chi_inv = np.linalg.inv(chi + eps*np.eye(3))

# Hickey–Moodera λ prefactor (complex)
e, hbar = 1.602e-19, 1.054e-34
m0, c    = 9.109e-31, 3e8
lam      = 1j * e*hbar*mu0 / (8 * gamma * m0**2 * c**2)

# Traditional Gilbert α for comparison
alpha = 0.02

# --------------------------------------------------------------
# 1. Drive parameters
# --------------------------------------------------------------
H0    = 1e4                   # drive amplitude (A/m)
omega = 2*np.pi * 5e7         # angular frequency (rad/s)
tmax  = 4*np.pi/omega         # two full periods
steps = 2000                  # number of time‐steps
t_eval = np.linspace(0, tmax, steps)

# --------------------------------------------------------------
# 2. Helper functions
# --------------------------------------------------------------
def skew(M):
    """Return matrix C such that C·v = M × v."""
    Mx, My, Mz = M
    return np.array([[ 0, -Mz,  My],
                     [Mz,   0, -Mx],
                     [-My, Mx,   0]])

def Heff(M, t):
    """Effective field: drive along z plus uniaxial anisotropy (z‐axis)."""
    H_app = np.array([0, 0, H0 * np.sin(omega * t)])
    H_ani = (2*Kuni/(mu0*Ms**2)) * (M[2]) * np.array([0,0,1])
    return H_app + H_ani

# --------------------------------------------------------------
# 3. RHS for α-free LLG (complex)
# --------------------------------------------------------------
def rhs_alpha_free(t, M_flat):
    # M_flat is length-3 complex array
    M = M_flat
    C = skew(M)
    A = np.eye(3, dtype=complex) + gamma * lam * (np.eye(3) + chi_inv) @ C
    rhs = -gamma * (C @ Heff(M.real, t))
    # solve A·dM = rhs  =>  dM = A⁻¹ rhs
    dM = np.linalg.solve(A, rhs)
    return dM

# --------------------------------------------------------------
# 4. RHS for traditional scalar-α LLG (real)
# --------------------------------------------------------------
def rhs_traditional(t, M):
    H = Heff(M, t)
    MxH   = np.cross(M, H)
    MxMxH = np.cross(M, MxH)
    return -gamma/(1+alpha**2) * (MxH + alpha * MxMxH)

# --------------------------------------------------------------
# 5. Set up integrators
# --------------------------------------------------------------
# α-free (complex)
solver_af = ode(rhs_alpha_free).set_integrator('zvode', method='bdf',
                                               atol=1e-10, rtol=1e-8)
M0_complex = (np.array([1,0,1])/np.sqrt(2) * Ms).astype(complex)
solver_af.set_initial_value(M0_complex, 0.0)

# traditional (real)
solver_tr = ode(rhs_traditional).set_integrator('vode', method='bdf',
                                                atol=1e-10, rtol=1e-8)
M0_real = np.array([1,0,1])/np.sqrt(2) * Ms
solver_tr.set_initial_value(M0_real, 0.0)

# storage
t_af = []; M_af = []
t_tr = []; M_tr = []

# --------------------------------------------------------------
# 6. Integrate with tqdm
# --------------------------------------------------------------
print("Integrating α-free LLG...")
for t_next in tqdm(t_eval[1:], desc="α-free  "):
    solver_af.integrate(t_next)
    t_af.append(solver_af.t)
    M_af.append(solver_af.y.copy())

print("Integrating traditional LLG...")
for t_next in tqdm(t_eval[1:], desc="traditional"):
    solver_tr.integrate(t_next)
    t_tr.append(solver_tr.t)
    M_tr.append(solver_tr.y.copy())

# Prepend initial condition
t_af = np.array([0.0] + t_af)
M_af = np.vstack([M0_complex, M_af]).T  # shape (3, steps)

t_tr = np.array([0.0] + t_tr)
M_tr = np.vstack([M0_real,   M_tr]).T   # shape (3, steps)

# --------------------------------------------------------------
# 7. Build hysteresis data & areas
# --------------------------------------------------------------
Hnorm_af = H0 * np.sin(omega * t_af) / Hk
Mz_af    = M_af[2].real / Ms

Hnorm_tr = H0 * np.sin(omega * t_tr) / Hk
Mz_tr    = M_tr[2] / Ms

area_af = abs(np.trapz(Mz_af, Hnorm_af))
area_tr = abs(np.trapz(Mz_tr, Hnorm_tr))

print(f"\nLoop area α-free : {area_af:.4f}")
print(f"Loop area α={alpha:4.2f}: {area_tr:.4f}\n")

# --------------------------------------------------------------
# 8. Plot
# --------------------------------------------------------------
plt.figure(figsize=(8,5))
plt.plot(Hnorm_af, Mz_af,  'b-',  label='α-free (tensor χ)')
plt.plot(Hnorm_tr, Mz_tr,  'r--', label=f'scalar α = {alpha}')
plt.xlabel(r'$H_{\mathrm{app}}/H_k$')
plt.ylabel(r'$M_z/M_s$')
plt.title('Hysteresis Loop Comparison')
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()
