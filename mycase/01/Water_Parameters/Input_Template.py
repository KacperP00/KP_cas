import cantera as ct
import numpy as np
import csv

# Create water object using the IAPWS95 equation of state
water_vap_Ideal = ct.Solution('h2o2.yaml','ohmech')
## water_vap_RK = ct.Solution('h2o2.yaml','ohmech-RK') nice to compare with ideal-gas
species = {"H2O": 1.0}

# ============================================================
# Define evenly spaced temperature & pressure ranges
# ============================================================
T_min, T_max, T_points = 300, 450.0, 150   # K
P_min, P_max, P_points = ct.one_atm, 8*ct.one_atm, 16       # Pa

T_range = np.linspace(T_min, T_max, T_points)
P_range = np.linspace(P_min, P_max, P_points)


print(f"{'P (MPa)':>8} {'rho (kg/m3)':>12} {'cp (J/kg-K)':>12} "
      f"{'mu (Pa·s)':>12} {'k (W/m-K)':>12} {'D_self (m2/s)':>14}")

# ============================================================
# CSV export
# ============================================================
csv_filename = "water_properties.csv"
with open(csv_filename, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["T [K]", "P [Pa]", "rho [kg/m3]", "cp [J/kg/K]", "mu (Pa·s)", "k (W/m-K)", "D_self (m2/s)"])

    for i, T in enumerate(T_range):
        for j, P in enumerate(P_range):
            water_vap_Ideal.TP = T, P

            rho = water_vap_Ideal.density
            cp = water_vap_Ideal.cp
            mu = water_vap_Ideal.viscosity
            k = water_vap_Ideal.thermal_conductivity

            # Heat of vaporization: h_g - h_l at same T
            water_vap_Ideal.TQ = T, 1  # saturated vapor
            h_v = water_vap_Ideal.enthalpy_mass
            water_vap_Ideal.TQ = T, 0  # saturated liquid
            h_l = water_vap_Ideal.enthalpy_mass
            Hvap = (h_v - h_l) / 1000  # convert to kJ/kg

            # Vapor pressure at this temperature
            water_vap_Ideal.TQ = T, 0
            Pvap = water_vap_Ideal.P / 1000  # kPa

            sigma = surface_tension_IAPWS(T)

            print(f"{'P (MPa)':>8} {'rho (kg/m3)':>12} {'cp (J/kg-K)':>12} "
            f"{'mu (Pa·s)':>12} {'k (W/m-K)':>12} {'D_self (m2/s)':>14}")
