# α-Free Landau–Lifshitz–Gilbert (LLG) Solver

## About  
This repository implements and compares two models of magnetization dynamics in a uniaxial ferromagnet:

1. **α-Free LLG**  
   A first-principles damping model derived by Hickey & Moodera from the Dirac–Pauli Hamiltonian, in which the Gilbert damping constant **α** is replaced by a complex susceptibility tensor \(\boldsymbol\chi_m\). No phenomenological fitting—damping and inertia emerge directly from spin–orbit coupling and the material’s magnetic response.

2. **Traditional Scalar-α LLG**  
   The standard Landau–Lifshitz–Gilbert equation with a user-specified damping constant α.

---

## 📖 Project Overview  
- **Goal**: Eliminate the ad-hoc α in micromagnetic simulations by using an α-free, tensorial damping derived from fundamental constants and the material’s susceptibility.  
- **Approach**:  
  1. Start from the Foldy–Wouthuysen expansion → intrinsic Gilbert torque  
  2. Express \(\dot{\mathbf B}/\dot{\mathbf M}\) via \(\boldsymbol\chi_m\)  
  3. Rearrange into a 3×3 linear solve per time step  
  4. Integrate with SciPy ODE solvers and compare loops  

---

## 📑 Theory & Derivation  

### 1. Intrinsic Gilbert Torque  
From spin–orbit coupling one finds:  
\[
\frac{d\mathbf M}{dt}
= -\gamma\,\mathbf M\times\mathbf H_\mathrm{eff}
\;-\;\gamma\,\lambda\,
\bigl(\mathbf1+\boldsymbol\chi_m^{-1}\bigr)\,
\Bigl(\mathbf M\times\frac{d\mathbf M}{dt}\Bigr),
\]
where  
\(\displaystyle\lambda=\frac{i\,e\,\hbar\,\mu_0}{8\,\gamma\,m_0^2c^2}\).  

### 2. Explicit Linear Form  
Define the skew-matrix \(\mathbf C(\mathbf M)\) such that \(\mathbf C\,\mathbf v=\mathbf M\times\mathbf v\). Then  
\[
\Bigl[\mathbf I+\gamma\lambda(\mathbf I+\chi_m^{-1})\,\mathbf C(\mathbf M)\Bigr]\,
\frac{d\mathbf M}{dt}
= -\gamma\,\mathbf C(\mathbf M)\,\mathbf H_\mathrm{eff}.
\]
Since \(\mathbf A=\mathbf I+\gamma\lambda(\mathbf I+\chi_m^{-1})\,\mathbf C\) is invertible,  
\(\dot{\mathbf M} = \mathbf A^{-1}\bigl(-\gamma\,\mathbf C\,\mathbf H_\mathrm{eff}\bigr)\).  

### 3. Comparison to Scalar-α LLG  
Standard LLG:  
\[
\frac{d\mathbf M}{dt}
= -\frac{\gamma}{1+\alpha^2}\Bigl[\mathbf M\times\mathbf H_\mathrm{eff}
+\alpha\,\mathbf M\times(\mathbf M\times\mathbf H_\mathrm{eff})\Bigr].
\]

---

## 🛠 Installation & Setup  

```bash
git clone https://github.com/<your-username>/alpha-free-llg.git
cd alpha-free-llg
python3 -m venv .venv
source .venv/bin/activate      # Linux/macOS
# .venv\Scripts\activate       # Windows PowerShell
pip install numpy scipy matplotlib tqdm
