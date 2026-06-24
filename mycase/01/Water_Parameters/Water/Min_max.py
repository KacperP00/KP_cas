import cantera as ct
import numpy as np
import csv

# Create water object using default the Reynolds  equation of state
water=ct.Water(backend="Reynolds")

# Create water object using alternate the IAPWS95 equation of state
#water=ct.Water(backend="IAPWS95")

# ============================================================
# Define evenly spaced temperature & pressure ranges
# ============================================================
T_min, T_max, T_points = 300, 500.0, 2   # K
P_min, P_max, P_points = 1e5, 10*1e5, 2       # bar

T_range = np.linspace(T_min, T_max, T_points)
P_range = np.linspace(P_min, P_max, P_points)


def surface_tension_IAPWS(T):
    """
    Surface tension of water using IAPWS correlation.
    Valid for 273 K < T < 647 K.
    Returns N/m.
    """
    Tc = 647.096  # K
    tau = 1 - T / Tc
    return 0.2358 * tau**1.256 * (1 - 0.625 * tau)

print(f"{'P (bar)':>8} {'T (K)':>12} {'mu (Pa·s)':>12} {'sigma (N/m)':>12}  "
      f" {'Hvap (kJ/kg)':>14} {'Pvap (kPa)':>12} "
      f" {'k (W/m-K)':>12} {'rho (kg/m3)':>12} {'cp (J/kg-K)':>12} ")

# ============================================================
# CSV export
# ============================================================
csv_filename = "water_properties.csv"
with open(csv_filename, "w", newline="") as outfile:
    writer = csv.writer(outfile)
    writer.writerow(["Pressure [bar]", "Temperature [K]", "liqViscosity [Pa·s]", "SurfaceTension [N/m]", "HeatOfVap [kJ/kg]", "vapPressure [kPa]", "liqThermalConductivity [W/m-K]", "liqDensity [kg/m3]", "liqHeatCapacity [J/kg/K]"])
    for i, P in enumerate(P_range):
        for j, T in enumerate(T_range):
            water.TP = T, P

            rho = water.density
            cp = water.cp
            mu = water.viscosity
            k = water.thermal_conductivity

            # Heat of vaporization: h_g - h_l at same T
            water.TQ = T, 1  # saturated vapor
            h_v = water.enthalpy_mass
            water.TQ = T, 0  # saturated liquid
            h_l = water.enthalpy_mass
            Hvap = (h_v - h_l) / 1000  # convert to kJ/kg

            # Vapor pressure at this temperature
            water.TQ = T, 0
            Pvap = water.P / 1000  # kPa

            sigma = surface_tension_IAPWS(T)
            writer.writerow([P/1e5, T, mu, sigma, Hvap, Pvap, k, rho, cp])
            print(f"{P/1e5:8.1f} {T:12.0f} {mu:12.2e} {sigma:12.4f} {Hvap:14.0f}  {Pvap:12.2f} {k:12.3f} {rho:12.1f} {cp:12.2f}")
