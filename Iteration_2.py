# Author: Kaleb Kassa
# Design a heat exchanger in which water vapor entering at a pressure 0.010627 MPa and temperature 60 °C is fully condensed, transferring the heat flow rate q = 43 kW to water, which enters at temperature t = 33°C. 
# The heat exchanger is designed to operate at a pressure of 0.010627 MPa and a temperature of 60 °C.
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Constants
# vapor inlet conditions
P_vapor = 0.010627  # MPa
T_vapor_in = 60  # °C
q = 43 # kW # heat transfer rate
T_sat= 47  # °C # saturation temperature at P_vapor
# water inlet conditions
T_water_in = 33  # °C
T_water_out = 42  # °C
# Tube side water properties at T_mean = 37.5 °C
Cp_water = 4.1795  # kJ/(kg·K) # specific heat capacity of water
rho_water = 993.11  # kg/m^3 # density of water
mu_water = 0.0006846  # Pa·s # dynamic viscosity of water
k_water = 0.6251  # W/(m·K) # thermal conductivity of water
#shell side vapor properties at T_sat = 47 °C and liquid water properties at T_mean_sat 
h_fg = 2389.2  # kJ/kg #  heat of vaporization
rho_vapor = 0.0722  # kg/m^3 # density of vapor
mu_liquid= 0.0005926  # Pa·s # dynamic viscosity of water
mu_vapor = 1.042e-5  # Pa·s # dynamic viscosity of vapor
rho_liquid = 989.69  # kg/m^3 # density of liquid water
k_liquid = 0.63612 # W/(m·K) # thermal conductivity of liquid water
# correction factror
Ft= 1 # from corection factor graph 
# Overall heat transfer coefficient from Bonancina et al. (2014)
K_overall_initial = 1000  # W/(m^2·K) # initial guess for overall heat transfer coefficient
# log mean temperature difference (LMTD) method
LMTD = ((T_sat - T_water_out) - (T_sat - T_water_in))/math.log((T_sat - T_water_out) / (T_sat - T_water_in))
# Calculate the required heat transfer area
#A_required = q*1000 / (K_overall_initial * LMTD*Ft)  # m^2
A_required = 2.40
# Calculate the mass flow rates
m_dot_water = q / (Cp_water * (T_water_out - T_water_in))  # kg/s
m_dot_vapor = q / h_fg  # kg/s
# tube material considerd copper
u_min=1 # m/s # minimum velocity
u_max=1.8 # m/s # maximum velocity
u=1.406 # m/s # velocity of water
# tube diameter
d_ex=0.0127 # m # external diameter of the tube
d_in=0.01021 # m # internal diameter of the tube
# cross-sectional area of the tube
Si=m_dot_water/(rho_water*u) # m^2
# number of tubes
Nt_1=Si/(math.pi*(d_in**2)/4) # number of tubes
Nt=math.ceil(Nt_1) # round up to the nearest whole number
n_passes_list=[1,2,4,6,8] # number of passes
# Calculate the number of tubes for each number of passes
n_tubes_list = []
for n_passes in n_passes_list:
    n_tubes = Nt*n_passes
    L = A_required/(math.pi*d_ex*n_tubes) # length of the tube
    n_tubes_list.append({
        'No. of passes': n_passes,
        'Total number of tubes': n_tubes,
        'Tube length (m)': L
    })

# Print the results as a table
df_tubes = pd.DataFrame(n_tubes_list)
#choose the number of passes
n_passes = 4
n_tubes = n_tubes_list[n_passes_list.index(n_passes)]['Total number of tubes']
L = n_tubes_list[n_passes_list.index(n_passes)]['Tube length (m)']
Pt=1.4*d_ex# m # external diameter of the tube for 30° inclination
# Iterative calculation for Dct and phi_n as in the attached table
phi_values = [0, 0.46, 0.34, 0.36, 0.37, 0.38]
dct_values = []

for phi in phi_values:
    Dct = (0.866 * (Pt**2) * 4 * n_tubes / (math.pi * (1 - phi)))**0.5
    dct_values.append(round(Dct, 3))

# Create DataFrame with extra spacing in the column names
df_phi_dct = pd.DataFrame({
    'phi      ': phi_values,
    'Dct (m)  ': dct_values
})
Dct_6 = dct_values[-1]  # Use the last value for Dct_6

# tube side (water)
#reynolds_water = (rho_water * u * d_in) / mu_water  # Reynolds number for water
reynolds_water = 20833.54
# Prdtl number for water
Pr_water = (mu_water * Cp_water*1000) / k_water  # Prandtl number for water
# Calculate the Nusselt number for water using Dittus-Boelter equation
alpha_water = 0.023 * (k_water/d_in)*(reynolds_water ** 0.8) * (Pr_water ** 0.33)  # Nusselt number for water
# Calculate the heat transfer coefficient for water
#assume T_sat - T_wall = 5K
T_wall = 45.3  # °C # wall temperature
alpha_1 = 0.725*(k_liquid**3 * rho_liquid*(rho_liquid - rho_vapor) * h_fg * 1000* 9.81/ ((T_sat - T_wall) * (mu_liquid * d_ex)))**0.25  # W/(m^2·K)
# for the bundle of tubes, the heat transfer coefficient is calculated using the following equation
N= (2/3)*(Dct_6/(Pt*0.866*2))
#rounding to the nearest highest whole number
N = math.ceil(N)  # round up to the nearest whole number
alpha_n=alpha_1*N**(-1/6)
# Calculate the overall heat transfer coefficient
fse=0    #m^2K/W # fouling factor for the shell side
fsi=0.00017 #m^2K/W # fouling factor for the tube side
k_copper = 310  # W/(m·K) # thermal conductivity of copper
#K_overall = 1/((d_ex/(d_in*alpha_water)) + (fsi*d_ex/d_in) + (d_ex* math.log(d_ex/d_in))/(2*k_copper) + (1/alpha_n))  # W/(m^2·K)
K_overall=2049.6
#new area required for the heat exchanger
A_required_new = q * 1000 / (Ft*K_overall * LMTD)  # m^2
# new length of the tube
L_new = A_required_new / (math.pi * d_ex * n_tubes)  # m
#delta_T validation
# Validate the temperature difference
delta_T_new= K_overall * Ft* LMTD / alpha_n  # K
#pressure drop calculations
# Calculate the pressure drop for the tube side (water)
f= 0.079 * (reynolds_water**-0.25)  # Darcy friction factor for turbulent flow
L=4* L_new  # total length of the tube
delta_p_friction = 2*f * (L/d_in) * (rho_water * u**2)  # pressure drop due to friction
delta_p_conc= 8 * rho_water * u**2  # pressure drop due to contraction
delta_p_total = delta_p_friction + delta_p_conc  # total pressure drop
#momentum pressure drop tube side
delta_p_shell = (m_dot_vapor**2/A_required_new **2)*((1/rho_vapor)-(1/rho_liquid))  # pressure drop fdue to momentum
# Shell-side pressure drop calculation (desuperheating + condensing)
def shell_side_pressure_drop(
    m_dot_vapor, q, Q_desuperheat, Q_condense, 
    d_ex, D_shell, L, n_tubes, mu_vapor, mu_liquid, 
    rho_vapor, rho_liquid
):
    P_t = 1.4 * d_ex  # Tube pitch
    B = 0.3 * D_shell  # Baffle spacing
    C = P_t - d_ex  # Clearance between tubes
    # Equivalent diameter for triangular pitch
    D_eq = (4 * (0.5 * P_t * 0.866 * P_t) - 0.5 * math.pi * d_ex**2) / (0.5 * math.pi * d_ex)
    # Desuperheating zone (single-phase)
    G_vapor = m_dot_vapor / (D_shell * C * B / P_t)
    Re_vapor = G_vapor * D_eq / mu_vapor
    f_vapor = 0.079 * Re_vapor**-0.25
    L_desuperheat = (Q_desuperheat / q) * L
    deltaP_desuperheat = 4 * f_vapor * (D_shell/D_eq) * (L_desuperheat/B) * (G_vapor**2/(2*rho_vapor)) #
    # Condensing zone (two-phase, Lockhart-Martinelli)
    x_avg = 0.5
    G_lo = m_dot_vapor / (D_shell * C * B / P_t)
    Re_lo = G_lo * D_eq / mu_liquid
    f_lo = 0.079 * Re_lo**-0.25
    deltaP_lo = 4 * f_lo * (D_shell/D_eq) * ((L - L_desuperheat)/B) * (G_lo**2/(2*rho_liquid))
    X_tt = ((1 - x_avg)/x_avg)**0.9 * (rho_vapor/rho_liquid)**0.5 * (mu_liquid/mu_vapor)**0.1
    phi_lo_squared = 1 + 20/X_tt + 1/X_tt**2
    deltaP_condense = deltaP_lo * phi_lo_squared
    # Nozzle losses
    deltaP_nozzle = 1.5 * (G_vapor**2 / (2*rho_vapor))
    deltaP_total = deltaP_desuperheat + deltaP_condense + deltaP_nozzle

    return deltaP_desuperheat, deltaP_condense, deltaP_total

# Example usage (add after your main calculations)
D_shell = 0.173  # m (example shell diameter)
mu_vapor = 1.04e-5  # Pa·s (vapor viscosity)
mu_liquid = 0.0005926  # Pa·s (liquid viscosity)
rho_vapor = 0.0722  # kg/m³
rho_liquid = 989.69  # kg/m³
Cp_vapor = 1.89  # kJ/(kg·K)
Q_desuperheat = m_dot_vapor * Cp_vapor * (T_vapor_in - T_sat)  # kW
Q_condense = q - Q_desuperheat  # kW

deltaP_desuperheat, deltaP_condense, deltaP_shell = shell_side_pressure_drop(
    m_dot_vapor, q, Q_desuperheat, Q_condense,
    d_ex, D_shell, L, n_tubes, mu_vapor, mu_liquid,
    rho_vapor, rho_liquid
)
# kern standalone vapor only

def kern_shell_side_pressure_drop(
    m_dot,           # Mass flow rate [kg/s]
    rho,             # Fluid density [kg/m³]
    mu,              # Dynamic viscosity [Pa·s]
    D_shell,         # Shell diameter [m]
    d_ex,            # Tube outer diameter [m]
    pitch_ratio=1.4, # Tube pitch ratio (P_t/d_ex)
   
    B=0.104
):
    # 1. Calculate geometric parameters
    P_t = pitch_ratio * d_ex          # Tube pitch [m]
    print(f"Tube pitch: {P_t:.4f} m")
    C = P_t - d_ex                    # Clearance between tubes [m]
    #B = baffle_ratio * D_shell        # Baffle spacing [m]
    
    # 2. Equivalent diameter for triangular pitch
    D_eq = (4 * ( 0.866 * P_t**2 -  math.pi * (d_ex/2)**2)) / (math.pi * d_ex)
    print(f"Equivalent diameter: {D_eq:.4f} m")
    
    # 3. Flow area and mass velocity
    A_flow = D_shell * C * B / P_t    # Flow area between baffles [m²]
    print(f"Flow area: {A_flow:.4f} m²")
    G = m_dot / A_flow                # Mass velocity [kg/(m²·s)]
    print(f"Mass velocity: {G:.4f} kg/(m²·s)")
    
    # 4. Reynolds number and friction factor
    Re = G * D_eq / mu
    f = math.exp(0.576 - 0.19 * math.log(Re))  # Kern's friction factor
    print(f"Reynolds number: {Re:.2f}")
    f1= 0.079 * Re**-0.25  # Alternative friction factor (for comparison)
    print(f"Friction factor (Kern): {f:.4f}, (Alternative): {f1:.4f}")
    print(f"Reynolds number: {Re:.2f}, Friction factor: {f:.4f}")
    # 5. Pressure drop calculation
    #L = D_shell  # Characteristic length (per baffle section)
    #deltaP = 4 * f * (D_shell/D_eq) * (L/B) * (G**2 / (2 * rho))
    n_b= 4  # Number of baffles (assumed)
    L_k= (n_b+1)*D_shell  # Total length of flow path (per baffle section)
    print(f"Total length of flow path: {L_k:.4f} m")
    mu_wall=0.00059124 # Pa·s # wall viscosity (assumed)
    u_k= m_dot / (rho * A_flow)  # Velocity in the shell side
    print(f"Velocity in the shell side: {u_k:.4f} m/s")
    deltaP = 0.5 * f * (L_k/D_eq) * rho*(u_k**2)*(mu/mu_wall)**(-0.14)  # Kern's pressure drop formula
    print(f"Pressure drop (Kern's method): {deltaP:.2f} Pa, (Alternative): {deltaP:.2f} Pa")
    delta_conc1= 8 * rho * u_k**2  # Pressure drop due to contraction
    print(f"Pressure drop due to contraction: {delta_conc1:.2f} Pa")
    delta_p_total= deltaP + delta_conc1  # Total pressure drop
    print(f"Total pressure drop (Kern's method): {delta_p_total:.2f} Pa")
    return deltaP

# --- Place this after your main calculations ---
# Example usage with your variables:
deltaP_kern = kern_shell_side_pressure_drop(
    m_dot=m_dot_vapor,
    rho=rho_vapor,
    mu=mu_vapor,
    D_shell=D_shell,
    d_ex=d_ex
)

print(f"Shell-side pressure drop (Kern method): {deltaP_kern:.2f} Pa ({deltaP_kern/1000:.3f} kPa)")

if deltaP_kern > 2000:
    print("Warning: ΔP exceeds 2 kPa limit for vacuum operation!")
else:
    print("ΔP within acceptable limits for vacuum operation.")

#kern standalone liquid only
def kern_shell_side_pressure_drop_liquid(
    m_dot,           # Mass flow rate [kg/s]
    rho,             # Fluid density [kg/m³]
    mu,              # Dynamic viscosity [Pa·s]
    D_shell,         # Shell diameter [m]
    d_ex,            # Tube outer diameter [m]
    pitch_ratio=1.4, # Tube pitch ratio (P_t/d_ex)
    B=0.104          # Baffle spacing [m]
):
    """
    Calculates shell-side pressure drop for liquid using Kern's method.
    Returns: Pressure drop [Pa]
    """
    P_t = pitch_ratio * d_ex          # Tube pitch [m]
    C = P_t - d_ex                    # Clearance between tubes [m]
    
    D_eq = (4 * ( 0.866 * P_t**2 -  math.pi * (d_ex/2)**2)) / (math.pi * d_ex)
    
    A_flow = D_shell * C * B / P_t    # Flow area between baffles [m²]
    G = m_dot / A_flow                # Mass velocity [kg/(m²·s)]
    
    Re = G * D_eq / mu
    f = 1  # Kern's friction factor from the graph (assumed for liquid) 
    print(f"Reynolds number: {Re:.2f}")
    print(f"Friction factor (Kern): {f:.4f}")
    
    n_b= 2  # Number of baffles (assumed)
    L_k= (n_b+1)*D_shell  # Total length of flow path (per baffle section)
    
    u_k= m_dot / (rho * A_flow)  # Velocity in the shell side
    deltaP = 0.5 * f * (L_k/D_eq) * rho*(u_k**2)*(mu/mu_water)**(-0.14)  # Kern's pressure drop formula
    delta_p_coc2= 8 * rho * u_k**2  # Pressure drop due to contraction
    print(f"Pressure drop (Kern's method): {deltaP:.2f} Pa")
    print(f"Pressure drop due to contraction: {delta_p_coc2:.2f} Pa")
    delta_p_total= deltaP + delta_p_coc2  # Total pressure drop
    print(f"total pressure drop liquid: {delta_p_total:.2f} Pa")
    
    return deltaP
# Example usage with your variables for liquid side:
deltaP_kern_liquid = kern_shell_side_pressure_drop_liquid(
    m_dot=m_dot_vapor,
    rho=rho_liquid,
    mu=mu_liquid,
    D_shell=D_shell,
    d_ex=d_ex
)
print(f"Shell-side pressure drop for liquid (Kern method): {deltaP_kern_liquid:.2f} Pa ({deltaP_kern_liquid/1000:.3f} kPa)")
if deltaP_kern_liquid > 2000:
    print("Warning: ΔP exceeds 2 kPa limit for vacuum operation!")
else:
    print("ΔP within acceptable limits for vacuum operation.")
# Summary of the results    
summary_data = [
    [1, "Log Mean Temperature Difference (LMTD)", f"{LMTD:.2f}", "°C"],
    [2, "Required heat transfer area", f"{A_required:.2f}", "m²"],
    [3, "Mass flow rate of water", f"{m_dot_water:.2f}", "kg/s"],
    [4, "Mass flow rate of vapor", f"{m_dot_vapor:.6f}", "kg/s"],
    [5, "Cross-sectional area of the tube", f"{Si:.4f}", "m²"],
    [6, "Number of tubes", f"{Nt}", "-"],
    [7, "Number of tubes for 4 passes", f"{n_tubes:.2f}", "-"],
    [8, "Length of each tube for 4 passes", f"{L:.2f}", "m"],
    [9, "Reynolds number for water", f"{reynolds_water:.2f}", "-"],
    [10, "Prandtl number for water", f"{Pr_water:.2f}", "-"],
    [11, "Nusselt number for water", f"{alpha_water:.2f}", "-"],
    [12, "Heat transfer coefficient for water", f"{alpha_1:.2f}", "W/(m²·K)"],
    [13, "Number of tubes in the bundle", f"{N}", "-"],
    [14, "Heat transfer coefficient for the bundle", f"{alpha_n:.2f}", "W/(m²·K)"],
    [15, "Overall heat transfer coefficient", f"{K_overall:.2f}", "W/(m²·K)"],
    [16, "New required heat transfer area", f"{A_required_new:.2f}", "m²"],
    [17, "New length of each tube", f"{L_new:.2f}", "m"],
    [18, "Pressure drop due to friction (tube side)", f"{delta_p_friction:.2f}", "Pa"],
    [19, "Pressure drop due to contraction (tube side)", f"{delta_p_conc:.2f}", "Pa"],
    [20, "Total pressure drop (tube side)", f"{delta_p_total:.2f}", "Pa"],
    [21, "Pressure drop due to momentum (tube side)", f"{delta_p_shell:.2f}", "Pa"],
    [22, "New temperature difference", f"{delta_T_new:.2f}", "K"],
    [23, "Shell-side ΔP (desuperheating)", f"{deltaP_desuperheat:.2f}", "Pa"],
    [24, "Shell-side ΔP (condensing, two-phase)", f"{deltaP_condense:.2f}", "Pa"],
    [25, "Shell-side ΔP (total)", f"{deltaP_shell:.2f}", "Pa"],
    [26, "Shell-side ΔP (Kern, vapor)", f"{deltaP_kern:.2f}", "Pa"],
    [27, "Shell-side ΔP (Kern, liquid)", f"{deltaP_kern_liquid:.2f}", "Pa"]
]

df_summary = pd.DataFrame(summary_data, columns=["S.No.", "Variable", "Value", "Unit"])

print("\nSummary Table of Calculated Results:")
print(df_summary.to_string(index=False))

# Plot phi vs Dct
plt.figure(figsize=(6,4))
plt.plot(dct_values, phi_values, marker='*', linestyle='-', color='r', markersize=10)
#plt.xticks(np.arange(0, max(dct_values) + 0.01, 0.01))
plt.ylabel('phi')
plt.xlabel('Dct (m)')
plt.title('phi vs Dct')
plt.grid(True)
plt.tight_layout()
plt.show()

# Plot to correlate number of passes, total number of tubes, and tube length
plt.figure(figsize=(8,5))
plt.plot(df_tubes['No. of passes'], df_tubes['Total number of tubes'], marker='o', label='Total number of tubes')
plt.plot(df_tubes['No. of passes'], df_tubes['Tube length (m)'], marker='s', label='Tube length (m)')
plt.xlabel('Number of Passes')
plt.ylabel('Value')
plt.title('Correlation of Number of Passes with Total Tubes and Tube Length')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# plot Dct vs phi for n=4 general case different from previous plot
from scipy.interpolate import interp1d
import numpy as np

Dct_general = [0.2, 0.25, 0.3, 0.35, 0.4, 0.6, 0.8, 1, 1.5, 2]
phi_general = [0.34, 0.23, 0.18, 0.155, 0.13, 0.11, 0.09, 0.07, 0.05, 0.031]

# Create interpolation function with extrapolation
f_phi = interp1d(Dct_general, phi_general, kind='linear', fill_value='extrapolate')

# Generate Dct values from 0.05 to 2.0
Dct_extrap = np.linspace(0.0, 2.0, 100)
phi_extrap = f_phi(Dct_extrap)

plt.figure(figsize=(6,4))
plt.plot(Dct_general, phi_general, 'o', color='b', label='Given data', markersize=5)
plt.plot(Dct_extrap, phi_extrap, '-', color='r', label='Extrapolated curve')
plt.xlabel('Dct (m)')
plt.ylabel('phi')
plt.title('Dct vs phi for n=4 General Case')
plt.xlim(0.05, 2.5)
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.minorticks_on()  # Enable minor ticks for more grid lines
plt.legend()
plt.tight_layout()
plt.show()

# plot Dct vs phi for n=4 general case different from previous plot
from scipy.interpolate import interp1d
import numpy as np

Dct_general = [0.2, 0.25, 0.3, 0.35, 0.4, 0.6, 0.8, 1, 1.5, 2]
phi_general = [0.34, 0.23, 0.18, 0.155, 0.13, 0.11, 0.09, 0.07, 0.05, 0.031]

# Create interpolation function with extrapolation
f_phi = interp1d(Dct_general, phi_general, kind='linear', fill_value='extrapolate')

# Generate Dct values from 0.05 to 2.0
Dct_extrap = np.linspace(0.05, 2.0, 100)
phi_extrap = f_phi(Dct_extrap)

plt.figure(figsize=(6,4))
plt.plot(Dct_general, phi_general, 'o', color='b', label='Given data', markersize=5)
plt.plot(Dct_extrap, phi_extrap, '-', color='r', label='Extrapolated curve')

# Mark the specified Dct points with different color
marked_points = [
    (0.118, 0.00),
    (0.16073, 0.46),
    (0.1454, 0.34),
    (0.1476, 0.36),
    (0.1488, 0.37),
    (0.14999, 0.38)
]
for dct, phi in marked_points:
    plt.scatter(dct, phi, color='g', s=10, zorder=5, label='Marked points' if dct == 0.100 else "")

plt.xlabel('Dct (m)')
plt.ylabel('phi')
plt.title('Dct vs phi for n=4 General Case')
plt.xlim(0.05, 2.5)
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.minorticks_on()
plt.legend()
plt.tight_layout()
plt.show()