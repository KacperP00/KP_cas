import cantera as ct
import CoolProp.CoolProp as CP

import numpy as np
import dataclasses
import csv

@dataclasses.dataclass
class Result:
    T: float
    P: float
    mu: float
    k: float
    rho: float
    cp: float
    diff: float # Dodany wspolczynnik dyfuzji

MECHANISM = "h2o2.yaml"
CASES = [
    {'mode':'ohmech', 'label': 'Ideal_Gas'},
    {'mode':'ohmech-RK', 'label': 'Real_Gas'},
    {'label': 'CoolProp'},
    ]

# ============================================================
# Define evenly spaced temperature & pressure ranges
# UWAGA (Zad. 2): Dla par gazu rozproszonych w komorze (Vapour) 
# wlasciwosci w kodach CFD podaje sie pod katem cisnienia otoczenia. 
# Ustawiamy wiec sztywno 1 atm, aby badac tylko wplyw Temperatury.
# ============================================================
T_MIN, T_MAX, T_POINTS = 300, 450.0, 151   # K
P_MIN, P_MAX, P_POINTS = ct.one_atm, ct.one_atm, 1  # Sztywne cisnienie 1 atm (Pa)
T_RANGE = np.linspace(T_MIN, T_MAX, T_POINTS)
P_RANGE = np.linspace(P_MIN, P_MAX, P_POINTS)

OUTPUT_PREFIX = "Vapour_Properties"

# ============================================================
# Zadanie 1: Funkcja wyliczajaca wspolczynnik dyfuzji (Rownanie Fullera)
# ============================================================
def calc_diffusion_fuller(T, P):
    """
    Wylicza wspolczynnik dyfuzji pary wodnej (A) w powietrzu (B) [m2/s]
    w zaleznosci od T [K] i P [Pa] uzywajac metody Fullera.
    """
    M_H2O = 18.015
    M_Air = 28.97
    M_AB = 2 * (1/M_H2O + 1/M_Air)**(-1)
    
    # Obietosci dyfuzyjne wg Fullera (z tablic)
    V_H2O = 12.7  
    V_Air = 20.1  
    
    P_atm = P / 101325.0 # We wzorze Fullera cisnienie musi byc w atmosferach
    
    # Wzor Fullera [cm2/s]
    D_cm2_s = (0.00143 * T**1.75) / (P_atm * M_AB**0.5 * (V_H2O**(1/3) + V_Air**(1/3))**2)
    
    # Konwersja na [m2/s]
    D_m2_s = D_cm2_s * 1e-4
    return D_m2_s

def vapour_prop(P_range, T_range, case):

    results = []   

    if case['label'] == 'CoolProp':
        print(f'\n ---------------------- Case: CoolProp ----------------------')
        print(f"{'P (bar)':>8} {'T (K)':>12} {'mu (Pa·s)':>12} {'k (W/m-K)':>12} {'rho (kg/m3)':>12} {'cp (J/kg-K)':>12} {'diff (m2/s)':>12}")
        for i, P in enumerate(P_range):
            for j, T in enumerate(T_range):
                try:
                    mu = CP.PropsSI("VISCOSITY", "T", T, "P", P, "Water")
                    k = CP.PropsSI("CONDUCTIVITY", "T", T, "P", P, "Water")
                    rho = CP.PropsSI("D", "T", T, "P", P, "Water")
                    cp = CP.PropsSI("CPMASS", "T", T, "P", P, "Water")
                    diff = calc_diffusion_fuller(T, P)

                except Exception as e:
                    print("ERROR at T=", T, "P=", P, "→", e)
                    continue
                
                results.append(Result(T=T, P=P, mu=mu, k=k, rho=rho, cp=cp, diff=diff))
                print(f"{P/1e5:8.1f} {T:12.0f} {mu:12.2e} {k:12.3f} {rho:12.1f} {cp:12.2f} {diff:12.2e}")

    else:
        mode = case['mode']
        label = case['label']
        water_vap = ct.Solution(MECHANISM, mode)

        print(f'\n ---------------------- Case: {label} ----------------------')
        print(f"{'P (bar)':>8} {'T (K)':>12} {'mu (Pa·s)':>12} {'k (W/m-K)':>12} {'rho (kg/m3)':>12} {'cp (J/kg-K)':>12} {'diff (m2/s)':>12}")

        for i, P in enumerate(P_range):
            for j, T in enumerate(T_range):
                water_vap.TP = T, P

                mu = water_vap.viscosity
                k = water_vap.thermal_conductivity
                rho = water_vap.density
                cp = water_vap.cp
                diff = calc_diffusion_fuller(T, P)

                results.append(Result(T=T, P=P, mu=mu, k=k, rho=rho, cp=cp, diff=diff))
                print(f"{P/1e5:8.1f} {T:12.0f} {mu:12.2e} {k:12.3f} {rho:12.1f} {cp:12.2f} {diff:12.2e}")

    return results

# ============================================================
# CSV export
# ============================================================

def write_propeties_csv(results, case_label):
    with open(f"{OUTPUT_PREFIX}_{case_label}_Results.csv", "w", newline="") as outfile:
        writer = csv.writer(outfile)
        writer.writerow(
            [
                "Pressure [bar]",
                "Temperature [K]",
                "VapViscosity [Pa·s]",
                "VapThermalConductivity [W/m-K]",
                "VapDensity [kg/m3]",
                "VapHeatCapacity [J/kg-K]",
                "DiffusionCoefficient [m2/s]", # Naprawione
            ]
        )

        for result in results:
            writer.writerow(
                [
                    result.P /1e5,
                    result.T,
                    result.mu,
                    result.k,
                    result.rho,
                    result.cp,
                    result.diff,
                ]
            )

if __name__ == "__main__":
    for case in CASES:
       results = vapour_prop(P_RANGE, T_RANGE, case)
       write_propeties_csv(results, case['label'])