#!/usr/bin/env python3
# main_fixed.py

import warnings
warnings.filterwarnings("ignore",
    message=".*Excess work done.*",
    category=UserWarning, module="scipy.integrate._ode")

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import ode
from tqdm import tqdm

# 0. Material constants
mu0   = 4*np.pi*1e-7
gamma = 2.211e5 / mu0
Ms    = 8.0e5
Kuni  = 5.0e4
Hk    = 2*Kuni/(mu0*Ms)

# 1. Drive params
H0    = 5 * Hk
omega = 2*np.pi * 1e6
tmax  = 4*np.pi/omega
steps = 2000
t_eval = np.linspace(0, tmax, steps)
dt     = t_eval[1] - t_eval[0]

# 2. Susceptibility & λ
chi0, chi1 = 1.5, 0.05
chi_perp   = chi0 - 1j*chi1
chi        = chi_perp*(np.eye(3)-np.outer([0,0,1],[0,0,1]))
eps        = 1e-8
chi_inv    = np.linalg.inv(chi + eps*np.eye(3))

e, hbar = 1.602e-19, 1.054e-34
m0, c    = 9.109e-31, 3e8
lam      = 1j*e*hbar*mu0/(8*gamma*m0**2*c**2)
lam     *= 1e9   # moderate boost

alpha    = 0.02

def skew(M):
    x,y,z = M
    return np.array([[0,-z, y],[z,0,-x],[-y,x,0]])

def Heff(M, t):
    H_app = np.array([H0*np.sin(omega*t),0,0])
    H_ani = (2*Kuni/(mu0*Ms**2))*M[2]*np.array([0,0,1])
    return H_app + H_ani

def rhs_af(t, Mc):
    C  = skew(Mc.real)
    A  = np.eye(3, dtype=complex) + gamma*lam*(np.eye(3)+chi_inv)@C
    rhs = -gamma*(C@Heff(Mc.real,t))
    return np.linalg.solve(A, rhs)

def rhs_tr(t, M):
    H   = Heff(M,t)
    MxH = np.cross(M,H)
    return -gamma/(1+alpha**2)*(MxH + alpha*np.cross(M,MxH))

# Integrators
solver_af = ode(rhs_af).set_integrator('zvode', method='bdf',
                                       atol=1e-9, rtol=1e-7, max_step=dt)
M0c = (np.array([1,0,1])/np.sqrt(2)*Ms).astype(complex)
solver_af.set_initial_value(M0c, 0.0)

solver_tr = ode(rhs_tr).set_integrator('vode', method='bdf',
                                       atol=1e-9, rtol=1e-7, max_step=dt)
M0r = np.array([1,0,1])/np.sqrt(2)*Ms
solver_tr.set_initial_value(M0r, 0.0)

t_af, M_af = [0.0], [M0c.copy()]
t_tr, M_tr = [0.0], [M0r.copy()]

# 7. Integration with renormalization
for tn in tqdm(t_eval[1:], desc="α-free  "):
    solver_af.integrate(tn)
    Mnew = solver_af.y
    norm = np.linalg.norm(Mnew.real)
    Mnorm = Mnew*(Ms/norm)
    solver_af.set_initial_value(Mnorm, solver_af.t)
    t_af.append(solver_af.t)
    M_af.append(Mnorm.copy())

for tn in tqdm(t_eval[1:], desc="traditional"):
    solver_tr.integrate(tn)
    t_tr.append(solver_tr.t)
    M_tr.append(solver_tr.y.copy())

# Build loops
t_af = np.array(t_af);  M_af = np.vstack(M_af).T
t_tr = np.array(t_tr);  M_tr = np.vstack(M_tr).T

Haf   = H0*np.sin(omega*t_af)/Hk
Mz_af = M_af[2].real / Ms
Htr   = H0*np.sin(omega*t_tr)/Hk
Mz_tr = M_tr[2]           / Ms

area_af = abs(np.trapezoid(Mz_af, Haf))
area_tr = abs(np.trapezoid(Mz_tr, Htr))
print(f"\nLoop area α-free : {area_af:.4f}")
print(f"Loop area α={alpha:.4f}: {area_tr:.4f}")

# Plot & save
plt.figure(figsize=(8,5))
plt.plot(Haf, Mz_af, 'b-',  label='α-free')
plt.plot(Htr, Mz_tr, 'r--', label=f'α={alpha}')
plt.xlabel('$H_{app}/H_k$'); plt.ylabel('$M_z/M_s$')
plt.legend(); plt.grid(alpha=0.3); plt.tight_layout()
plt.savefig('hysteresis_fixed.png', dpi=300)
print("Saved hysteresis_fixed.png")
