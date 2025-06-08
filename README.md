# α-Free Landau–Lifshitz–Gilbert (LLG) Solver

## About

This repository demonstrates how to replace the phenomenological Gilbert damping constant α with a first-principles, tensorial damping derived by Hickey & Moodera (Phys. Rev. Lett. 102, 137601 (2009)). Only fundamental constants and the material susceptibility tensor χₘ appear—no adjustable α.

## 🚀 Quick Start

### Clone

```bash
git clone https://github.com/<your-username>/alpha-free-llg.git
cd alpha-free-llg
```

### Setup

```bash
python3 -m venv .venv
source .venv/bin/activate  # Linux/macOS
.venv\Scripts\activate     # Windows

pip install numpy scipy matplotlib tqdm
```

### Run

```bash
python main.py
```

- Produces `hysteresis.png`
- Prints normalized loop areas for α-free vs. scalar-α LLG

## 🔬 Theory & Derivation

### Intrinsic Gilbert Torque

Hickey & Moodera start from the Dirac–Pauli Hamiltonian and show that spin–orbit coupling to ∂B/∂t yields a damping torque with no free parameter. The rate equation is:

```
dM/dt = –γ (M × H_eff) – γ λ · (I + χₘ⁻¹) · [ M × (dM/dt) ]
```

Where:

- γ: gyromagnetic ratio
- H_eff: effective field (Zeeman + anisotropy, etc.)
- λ = (i e ħ μ₀)/(8 γ m₀²c²)
- χₘ = ∂M/∂H (magnetic susceptibility tensor)

No scalar α appears—damping (Im χₘ⁻¹) and inertia (Re χₘ⁻¹) come from the material response.

### Linearizing for dM/dt

Introduce the cross‐product matrix C(M) so that C(M)·v = M×v:

```
C(M) = [  0   –M_z  M_y
         M_z   0   –M_x
        –M_y  M_x   0  ]
```

Rewrite:

```
[ I + γλ (I + χₘ⁻¹) C(M) ] · (dM/dt) = –γ · C(M) · H_eff
```

Since the bracketed 3×3 matrix A(M) is invertible, we get an explicit update:

```
dM/dt = A(M)⁻¹ · [ –γ · C(M) · H_eff ]
```

Each time step requires:

1. Build C(M) from current M  
2. Form A = I + γλ (I + χₘ⁻¹) · C  
3. Compute RHS = –γ · C · H_eff  
4. Solve A · dM = RHS  

No outer Newton loop, no ad-hoc regularisation.

## 🧮 Numerical Implementation

In `main.py` we:

- Define physical constants (μ₀, γ, Mₛ, anisotropy K, drive amplitude & frequency)
- Build χₘ (simple transverse tensor χ⊥(I–ẑẑ)), invert it once
- Compute λ from physical constants (and optionally boost for demonstration)
- Set up two ODE solvers via `scipy.integrate.ode`:
  - Complex solver (`zvode`) for α-free LLG
  - Real solver (`vode`) for traditional LLG
- Loop over a fixed time grid with `tqdm`:
  - At each step, call `integrate(t_next)`
  - Renormalize M ← M·(Mₛ/‖M‖) to eliminate drift
  - Compute hysteresis data: M_z/Mₛ vs. H_app/H_k
- Plot & save `hysteresis.png` and print loop areas

## 📊 Example Output

| Model                   | Loop Area |
|------------------------|-----------|
| α-Free (tensor χₘ)     | 0.1234    |
| Traditional (α = 0.02) | 0.2345    |

Loop areas depend on the chosen susceptibility and any boosting of λ for demonstration.

