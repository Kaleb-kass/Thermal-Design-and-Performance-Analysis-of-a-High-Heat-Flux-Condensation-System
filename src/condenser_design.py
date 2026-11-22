"""
Thermal design of a horizontal shell-and-tube condenser.

Water vapour entering at 0.010627 MPa and 60 °C is fully condensed on the shell side,
rejecting 43 kW to cooling water that enters the tubes at 33 °C.

The script reproduces the second design iteration of the accompanying report
(report/Thermal-Design-High-Heat-Flux-Condensation-System.pdf):
LMTD, tube count and length for different numbers of passes, bundle and shell diameter,
heat transfer coefficients, and tube-side / shell-side pressure drops.

Usage:
    python src/condenser_design.py            # print results, save figures to figures/
    python src/condenser_design.py --show     # also open the figures in a window

Author: Kaleb Kassa Lafebo, University of Padova
"""

import argparse
import math
from pathlib import Path

import matplotlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

FIGURES_DIR = Path(__file__).resolve().parent.parent / "figures"

# =============================================================================
# 1. Design data and thermophysical properties (REFPROP)
# =============================================================================
# Vapour (shell side)
P_vapor = 0.010627      # MPa, inlet pressure
T_vapor_in = 60         # °C, inlet temperature (superheated)
T_sat = 47              # °C, saturation temperature at P_vapor
q = 43                  # kW, heat duty

# Cooling water (tube side)
T_water_in = 33         # °C
T_water_out = 42        # °C (T_sat - T_water_out = 5 K, simplified method "case 3")

# Water properties at the mean temperature 37.5 °C
Cp_water = 4.1795       # kJ/(kg·K)
rho_water = 993.11      # kg/m³
mu_water = 0.0006846    # Pa·s
k_water = 0.6251        # W/(m·K)

# Condensate at the mean film temperature, vapour at T_sat (second iteration)
h_fg = 2389.2           # kJ/kg, latent heat
rho_vapor = 0.0722      # kg/m³
mu_liquid = 0.0005926   # Pa·s
mu_vapor = 1.042e-5     # Pa·s
rho_liquid = 989.69     # kg/m³
k_liquid = 0.63612      # W/(m·K)

Ft = 1                  # LMTD correction factor (isothermal condensation)
K_overall_initial = 1000  # W/(m²·K), first guess (Bonacina et al., Trasmissione del calore)

# Tubes: copper, 1/2" OD
d_ex = 0.0127           # m, external diameter
d_in = 0.01021          # m, internal diameter (wall 1.245 mm)
k_copper = 310          # W/(m·K)
u_min, u_max = 1, 1.8   # m/s, allowable water velocity
u = 1.406               # m/s, design water velocity

# Fouling factors
fse = 0                 # m²K/W, shell side
fsi = 0.00017           # m²K/W, tube side

# Values fixed from the report's second iteration (kept so the script reproduces it exactly)
A_REQUIRED = 2.40       # m², heat transfer area from iteration 1 (Q / (K·LMTD) with K ≈ 2046)
RE_WATER = 20833.54     # tube-side Reynolds number reported (ρ·u·d_in/μ with u = 1.41 m/s)
K_OVERALL = 2049.6      # W/(m²·K), overall coefficient of iteration 2
T_WALL = 45.3           # °C, wall temperature used for the condensation coefficient
N_PASSES = 4            # selected number of tube passes
PHI_ITERATIONS = [0, 0.46, 0.34, 0.36, 0.37, 0.38]  # bundle correction factor iterations
D_SHELL = 0.173         # m, shell diameter used for shell-side pressure drops


def main(show: bool) -> None:
    if not show:
        matplotlib.use("Agg")
    FIGURES_DIR.mkdir(exist_ok=True)

    # =========================================================================
    # 2. Temperature difference and flow rates
    # =========================================================================
    LMTD = ((T_sat - T_water_out) - (T_sat - T_water_in)) / math.log((T_sat - T_water_out) / (T_sat - T_water_in))
    A_required = A_REQUIRED
    m_dot_water = q / (Cp_water * (T_water_out - T_water_in))   # kg/s
    m_dot_vapor = q / h_fg                                       # kg/s

    # =========================================================================
    # 3. Tube count, passes and length
    # =========================================================================
    Si = m_dot_water / (rho_water * u)                  # m², flow area per pass
    Nt = math.ceil(Si / (math.pi * d_in**2 / 4))        # tubes per pass

    n_passes_list = [1, 2, 4, 6, 8]
    passes_table = []
    for n_passes in n_passes_list:
        n_tubes_i = Nt * n_passes
        passes_table.append({
            "No. of passes": n_passes,
            "Total number of tubes": n_tubes_i,
            "Tube length (m)": A_required / (math.pi * d_ex * n_tubes_i),
        })
    df_tubes = pd.DataFrame(passes_table)

    selected = passes_table[n_passes_list.index(N_PASSES)]
    n_tubes = selected["Total number of tubes"]
    L_tube_4p = selected["Tube length (m)"]

    # =========================================================================
    # 4. Bundle diameter (iterative, triangular pitch 30°)
    # =========================================================================
    Pt = 1.4 * d_ex                                     # m, tube pitch
    dct_values = [round((0.866 * Pt**2 * 4 * n_tubes / (math.pi * (1 - phi)))**0.5, 3) for phi in PHI_ITERATIONS]
    Dct = dct_values[-1]

    # =========================================================================
    # 5. Heat transfer coefficients
    # =========================================================================
    reynolds_water = RE_WATER
    Pr_water = (mu_water * Cp_water * 1000) / k_water
    # Tube side: Dittus–Boelter
    alpha_water = 0.023 * (k_water / d_in) * reynolds_water**0.8 * Pr_water**0.33
    # Shell side: Nusselt film condensation on a single horizontal tube
    alpha_1 = 0.725 * (k_liquid**3 * rho_liquid * (rho_liquid - rho_vapor) * h_fg * 1000 * 9.81
                       / ((T_sat - T_WALL) * (mu_liquid * d_ex)))**0.25
    # Tube bundle: N tubes in a vertical row
    N = math.ceil((2 / 3) * (Dct / (Pt * 0.866 * 2)))
    alpha_n = alpha_1 * N**(-1 / 6)
    # Overall coefficient (formula kept for reference; iteration-2 value is used)
    # K = 1 / (d_ex/(d_in·α_w) + fsi·d_ex/d_in + fse + d_ex·ln(d_ex/d_in)/(2·k_Cu) + 1/α_N)
    K_overall = K_OVERALL

    A_required_new = q * 1000 / (Ft * K_overall * LMTD)          # m²
    L_new = A_required_new / (math.pi * d_ex * n_tubes)          # m, length of each tube
    delta_T_new = K_overall * Ft * LMTD / alpha_n                # K, check of T_sat - T_wall

    # =========================================================================
    # 6. Pressure drops
    # =========================================================================
    # Tube side (water)
    f = 0.079 * reynolds_water**-0.25
    L_path = 4 * L_new                                           # m, total flow path (4 passes)
    delta_p_friction = 2 * f * (L_path / d_in) * (rho_water * u**2)
    delta_p_conc = 8 * rho_water * u**2
    delta_p_tube = delta_p_friction + delta_p_conc
    # Momentum term of the condensing vapour (shell side)
    delta_p_momentum = (m_dot_vapor**2 / A_required_new**2) * ((1 / rho_vapor) - (1 / rho_liquid))

    # Shell side: desuperheating + two-phase condensation (Lockhart–Martinelli)
    mu_v = 1.04e-5                                               # Pa·s, vapour viscosity used below
    Cp_vapor = 1.89                                              # kJ/(kg·K)
    Q_desuperheat = m_dot_vapor * Cp_vapor * (T_vapor_in - T_sat)  # kW
    dp_desuperheat, dp_condense, dp_shell_2phase = shell_side_pressure_drop(
        m_dot_vapor, q, Q_desuperheat, d_ex, D_SHELL, L_path, mu_v, mu_liquid, rho_vapor, rho_liquid)

    # Shell side: Kern method, vapour only and liquid only
    dp_kern_vapor = kern_pressure_drop(m_dot_vapor, rho_vapor, mu_v, D_SHELL, d_ex,
                                       f_friction=None, n_baffles=4, mu_wall=0.00059124, label="vapour")
    dp_kern_liquid = kern_pressure_drop(m_dot_vapor, rho_liquid, mu_liquid, D_SHELL, d_ex,
                                        f_friction=1, n_baffles=2, mu_wall=mu_water, label="liquid")

    # =========================================================================
    # 7. Results
    # =========================================================================
    print("\nTube count and length for different numbers of passes:")
    print(df_tubes.to_string(index=False))

    rows = [
        ("Log mean temperature difference (LMTD)", f"{LMTD:.2f}", "°C"),
        ("Required heat transfer area", f"{A_required:.2f}", "m²"),
        ("Mass flow rate of water", f"{m_dot_water:.2f}", "kg/s"),
        ("Mass flow rate of vapour", f"{m_dot_vapor:.6f}", "kg/s"),
        ("Flow area per pass", f"{Si:.4f}", "m²"),
        ("Number of tubes per pass", f"{Nt}", "-"),
        ("Total number of tubes (4 passes)", f"{n_tubes:.2f}", "-"),
        ("Length of each tube (4 passes)", f"{L_tube_4p:.2f}", "m"),
        ("Reynolds number, water", f"{reynolds_water:.2f}", "-"),
        ("Prandtl number, water", f"{Pr_water:.2f}", "-"),
        ("Heat transfer coefficient, tube side (Dittus–Boelter)", f"{alpha_water:.2f}", "W/(m²·K)"),
        ("Condensation coefficient, single tube (Nusselt)", f"{alpha_1:.2f}", "W/(m²·K)"),
        ("Tubes in a vertical row, N", f"{N}", "-"),
        ("Condensation coefficient, tube bundle", f"{alpha_n:.2f}", "W/(m²·K)"),
        ("Overall heat transfer coefficient", f"{K_overall:.2f}", "W/(m²·K)"),
        ("Heat transfer area (iteration 2)", f"{A_required_new:.2f}", "m²"),
        ("Length of each tube (iteration 2)", f"{L_new:.2f}", "m"),
        ("Tube-side pressure drop, friction", f"{delta_p_friction:.2f}", "Pa"),
        ("Tube-side pressure drop, concentrated", f"{delta_p_conc:.2f}", "Pa"),
        ("Tube-side pressure drop, total", f"{delta_p_tube:.2f}", "Pa"),
        ("Shell-side momentum pressure drop", f"{delta_p_momentum:.2f}", "Pa"),
        ("T_sat − T_wall check", f"{delta_T_new:.2f}", "K"),
        ("Shell-side ΔP, desuperheating", f"{dp_desuperheat:.2f}", "Pa"),
        ("Shell-side ΔP, condensing (two-phase)", f"{dp_condense:.2f}", "Pa"),
        ("Shell-side ΔP, total (two-phase model)", f"{dp_shell_2phase:.2f}", "Pa"),
        ("Shell-side ΔP, Kern (vapour only)", f"{dp_kern_vapor:.2f}", "Pa"),
        ("Shell-side ΔP, Kern (liquid only)", f"{dp_kern_liquid:.2f}", "Pa"),
    ]
    df_summary = pd.DataFrame([(i + 1, *r) for i, r in enumerate(rows)], columns=["S.No.", "Variable", "Value", "Unit"])
    print("\nSummary of calculated results:")
    print(df_summary.to_string(index=False))

    for name, dp in (("vapour", dp_kern_vapor), ("liquid", dp_kern_liquid)):
        status = "exceeds" if dp > 2000 else "is within"
        print(f"Kern shell-side ΔP ({name}): {dp:.2f} Pa {status} the 2 kPa limit for vacuum operation.")

    # =========================================================================
    # 8. Figures
    # =========================================================================
    plot_passes(df_tubes)
    plot_dct_iterations(dct_values)
    plot_phi_extrapolation()
    print(f"\nFigures saved to {FIGURES_DIR}")
    if show:
        plt.show()


def shell_side_pressure_drop(m_dot_vapor, q, Q_desuperheat, d_ex, D_shell, L,
                             mu_vapor, mu_liquid, rho_vapor, rho_liquid):
    """Shell-side ΔP: single-phase desuperheating zone + two-phase condensing zone (Lockhart–Martinelli)."""
    P_t = 1.4 * d_ex
    B = 0.3 * D_shell                   # baffle spacing
    C = P_t - d_ex                      # clearance between tubes
    D_eq = (4 * (0.5 * P_t * 0.866 * P_t) - 0.5 * math.pi * d_ex**2) / (0.5 * math.pi * d_ex)

    # Desuperheating zone (vapour only)
    G_vapor = m_dot_vapor / (D_shell * C * B / P_t)
    Re_vapor = G_vapor * D_eq / mu_vapor
    f_vapor = 0.079 * Re_vapor**-0.25
    L_desuperheat = (Q_desuperheat / q) * L
    dp_desuperheat = 4 * f_vapor * (D_shell / D_eq) * (L_desuperheat / B) * (G_vapor**2 / (2 * rho_vapor))

    # Condensing zone (liquid-only ΔP times two-phase multiplier, average quality 0.5)
    x_avg = 0.5
    G_lo = m_dot_vapor / (D_shell * C * B / P_t)
    Re_lo = G_lo * D_eq / mu_liquid
    f_lo = 0.079 * Re_lo**-0.25
    dp_lo = 4 * f_lo * (D_shell / D_eq) * ((L - L_desuperheat) / B) * (G_lo**2 / (2 * rho_liquid))
    X_tt = ((1 - x_avg) / x_avg)**0.9 * (rho_vapor / rho_liquid)**0.5 * (mu_liquid / mu_vapor)**0.1
    dp_condense = dp_lo * (1 + 20 / X_tt + 1 / X_tt**2)

    dp_nozzle = 1.5 * (G_vapor**2 / (2 * rho_vapor))
    return dp_desuperheat, dp_condense, dp_desuperheat + dp_condense + dp_nozzle


def kern_pressure_drop(m_dot, rho, mu, D_shell, d_ex, f_friction, n_baffles, mu_wall,
                       label, pitch_ratio=1.4, B=0.104):
    """Shell-side distributed ΔP with Kern's method for a single phase.

    f_friction=None uses Kern's correlation f = exp(0.576 − 0.19·ln Re); otherwise the given
    value is used (read from the friction-factor chart for laminar flow).
    """
    P_t = pitch_ratio * d_ex
    C = P_t - d_ex
    D_eq = (4 * (0.866 * P_t**2 - math.pi * (d_ex / 2)**2)) / (math.pi * d_ex)
    A_flow = D_shell * C * B / P_t
    G = m_dot / A_flow
    Re = G * D_eq / mu
    f = math.exp(0.576 - 0.19 * math.log(Re)) if f_friction is None else f_friction

    L_k = (n_baffles + 1) * D_shell
    u_k = m_dot / (rho * A_flow)
    dp = 0.5 * f * (L_k / D_eq) * rho * u_k**2 * (mu / mu_wall)**(-0.14)
    dp_conc = 8 * rho * u_k**2
    print(f"\nKern method, {label} only: D_eq = {D_eq:.4f} m, A_flow = {A_flow:.4f} m², "
          f"u = {u_k:.4f} m/s, Re = {Re:.2f}, f = {f:.4f}")
    print(f"  distributed ΔP = {dp:.2f} Pa, concentrated ΔP = {dp_conc:.2f} Pa, total = {dp + dp_conc:.2f} Pa")
    return dp


def plot_passes(df_tubes):
    plt.figure(figsize=(8, 5))
    plt.plot(df_tubes["No. of passes"], df_tubes["Total number of tubes"], marker="o", label="Total number of tubes")
    plt.plot(df_tubes["No. of passes"], df_tubes["Tube length (m)"], marker="s", label="Tube length (m)")
    plt.xlabel("Number of passes")
    plt.ylabel("Value")
    plt.title("Number of passes vs. total tubes and tube length")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "passes_tubes_length.png", dpi=150)


def plot_dct_iterations(dct_values):
    plt.figure(figsize=(6, 4))
    plt.plot(dct_values, PHI_ITERATIONS, marker="*", linestyle="-", color="r", markersize=10)
    plt.xlabel("Dct (m)")
    plt.ylabel("phi")
    plt.title("Bundle diameter iterations: phi vs Dct")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "dct_iterations.png", dpi=150)


def plot_phi_extrapolation():
    """Correction factor phi vs Dct for 4 passes, extrapolated below the chart's 200 mm limit."""
    Dct_general = [0.2, 0.25, 0.3, 0.35, 0.4, 0.6, 0.8, 1, 1.5, 2]
    phi_general = [0.34, 0.23, 0.18, 0.155, 0.13, 0.11, 0.09, 0.07, 0.05, 0.031]
    f_phi = interp1d(Dct_general, phi_general, kind="linear", fill_value="extrapolate")
    Dct_extrap = np.linspace(0.05, 2.0, 100)

    plt.figure(figsize=(6, 4))
    plt.plot(Dct_general, phi_general, "o", color="b", label="Chart data", markersize=5)
    plt.plot(Dct_extrap, f_phi(Dct_extrap), "-", color="r", label="Extrapolated curve")
    iterations = [(0.118, 0.00), (0.16073, 0.46), (0.1454, 0.34), (0.1476, 0.36), (0.1488, 0.37), (0.14999, 0.38)]
    for i, (dct, phi) in enumerate(iterations):
        plt.scatter(dct, phi, color="g", s=10, zorder=5, label="Iterations" if i == 0 else None)
    plt.xlabel("Dct (m)")
    plt.ylabel("phi")
    plt.title("Correction factor phi vs Dct (4 passes)")
    plt.xlim(0.05, 2.5)
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.minorticks_on()
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "phi_extrapolation.png", dpi=150)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Thermal design of a shell-and-tube steam condenser.")
    parser.add_argument("--show", action="store_true", help="open the figures in a window")
    main(parser.parse_args().show)
