import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# ---------------------------
# Input data
# ---------------------------
# System sizes (L)
L = np.array([155.136, 177.586, 195.459])

# Critical temperatures (K) and uncertainties
Tc = np.array([6651, 6633, 6618], dtype=float)
Tc_err = np.array([20, 30, 23], dtype=float)

# ---------------------------
# Finite-size scaling function
# ---------------------------
nu = 0.63  # 3D Ising critical exponent
theta = 0.54
def Tc_scaling(L, Tc_inf, a):
    return Tc_inf + a * L**(-(1+theta)/nu)

# ---------------------------
# Fit the data
# ---------------------------
popt, pcov = curve_fit(Tc_scaling, L, Tc, sigma=Tc_err, absolute_sigma=True)
Tc_inf_fit, a_fit = popt
Tc_inf_err, a_err = np.sqrt(np.diag(pcov))

print(f"Estimated bulk critical temperature: {Tc_inf_fit:.1f} ± {Tc_inf_err:.1f} K")
print(f"Fit parameter a: {a_fit:.2e} ± {a_err:.2e}")

# ---------------------------
# Generate fit curve
# ---------------------------
# L_fit = np.linspace(min(L)*0.9, max(L)*1.1, 200)
L_fit = np.linspace(min(L)*0.8, max(L)*3, 200)
Tc_fit = Tc_scaling(L_fit, Tc_inf_fit, a_fit)

# ---------------------------
# Plot results
# ---------------------------
# plt.errorbar(L, Tc, yerr=Tc_err, fmt='o', label='Simulation data', capsize=5)
plt.scatter(L, Tc, color='red', marker='o', label='Simulation data')
# plt.plot(L_fit, Tc_fit, '-', label=f'Fit: Tc(∞)={Tc_inf_fit:.1f}±{Tc_inf_err:.1f} K')
plt.plot(L_fit, Tc_fit, '--', label=f'Fit: $T_c$(∞)={Tc_inf_fit:.1f} K')

plt.xlabel(r"System size $L_z$ (Angstrom)")
plt.ylabel("Critical Temperature $T_c$ (K)")
plt.title("Finite-size scaling of $T_c$")
plt.legend()
plt.grid(True)
plt.savefig("Tc_finite_size_scaling_corrector.png", dpi=300)
# plt.show()
