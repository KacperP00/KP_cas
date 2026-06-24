import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os

# ==============================================================================
# CZĘŚĆ 1: ANALIZA PORÓWNAWCZA (Wykresy i Statystyka)
# ==============================================================================
def analizuj_i_rysuj_modele():
    try:
        df_ideal = pd.read_csv("Vapour_Properties_Ideal_Gas_Results.csv", encoding='cp1252')
        df_real = pd.read_csv("Vapour_Properties_Real_Gas_Results.csv", encoding='cp1252')
        df_cool = pd.read_csv("Vapour_Properties_CoolProp_Results.csv", encoding='cp1252')
    except FileNotFoundError:
        print("Brak plików CSV z danymi oparów. Upewnij się, że są w tym samym folderze.")
        return

    # 1. Raport Statystyczny dla 1.0 bar (Warunki rozproszenia w komorze)
    p_target = 1.0
    ideal_1bar = df_ideal[df_ideal["Pressure [bar]"] == p_target].reset_index()
    cool_1bar = df_cool[df_cool["Pressure [bar]"] == p_target].reset_index()
    
    # Skupiamy się na t > 375K aby uniknąć wpływu skraplania (fazy ciekłej)
    ideal_1bar_vap = ideal_1bar[ideal_1bar["Temperature [K]"] >= 375]
    cool_1bar_vap = cool_1bar[cool_1bar["Temperature [K]"] >= 375]

    err_cp = np.mean(np.abs(ideal_1bar_vap["VapHeatCapacity [J/kg-K]"] - cool_1bar_vap["VapHeatCapacity [J/kg-K]"]) / cool_1bar_vap["VapHeatCapacity [J/kg-K]"]) * 100
    err_rho = np.mean(np.abs(ideal_1bar_vap["VapDensity [kg/m3]"] - cool_1bar_vap["VapDensity [kg/m3]"]) / cool_1bar_vap["VapDensity [kg/m3]"]) * 100

    print(f"--- RAPORT PORÓWNAWCZY DLA PARY (P = {p_target} bar, T > 375 K) ---")
    print(f"Średnia odchyłka Gęstości (Ideal Gas vs CoolProp): {err_rho:.2f}%")
    print(f"Średnia odchyłka Ciepła Właściwego Cp: {err_cp:.2f}%\n")

    # 2. Rysowanie wykresów dla dwóch skrajnych ciśnień (1.0 bar i 8.0 bar)
    p_levels = [1.0, 8.0]
    
    for p in p_levels:
        col_p = "Pressure [bar]"
        
        # Filtrujemy dane dla danego ciśnienia
        ideal_p = df_ideal[df_ideal[col_p] == p]
        real_p = df_real[df_real[col_p] == p]
        cool_p = df_cool[df_cool[col_p] == p]
        
        fig, axs = plt.subplots(2, 2, figsize=(14, 10))
        fig.suptitle(f"Zestawienie modeli termodynamicznych dla pary wodnej (Ciśnienie {p} bar)", fontsize=16)

        # Wykres 1: Ciepło właściwe Cp
        ax = axs[0, 0]
        ax.plot(ideal_p["Temperature [K]"], ideal_p["VapHeatCapacity [J/kg-K]"], label="Ideal Gas (Cantera)", linestyle='-', linewidth=2)
        ax.plot(real_p["Temperature [K]"], real_p["VapHeatCapacity [J/kg-K]"], label="Real Gas (RK Cantera)", linestyle='--', linewidth=2)
        ax.plot(cool_p["Temperature [K]"], cool_p["VapHeatCapacity [J/kg-K]"], label="CoolProp (IAPWS-95)", linestyle='-.', linewidth=2)
        ax.set_title("Ciepło właściwe ($C_p$) pary wodnej")
        ax.set_xlabel("Temperatura [K]")
        ax.set_ylabel("$C_p$ [J/(kg K)]")
        ax.grid(True, linestyle='--', alpha=0.7)
        ax.legend()

        # Wykres 2: Gęstość
        ax = axs[0, 1]
        ax.plot(ideal_p["Temperature [K]"], ideal_p["VapDensity [kg/m3]"], label="Ideal Gas", linestyle='-', linewidth=2)
        ax.plot(real_p["Temperature [K]"], real_p["VapDensity [kg/m3]"], label="Real Gas", linestyle='--', linewidth=2)
        ax.plot(cool_p["Temperature [K]"], cool_p["VapDensity [kg/m3]"], label="CoolProp", linestyle='-.', linewidth=2)
        ax.set_title("Gęstość ($\\rho$) pary wodnej")
        ax.set_xlabel("Temperatura [K]")
        ax.set_ylabel("Gęstość [kg/m³]")
        ax.grid(True, linestyle='--', alpha=0.7)
        ax.legend()
        
        # Wykres 3: Przewodność cieplna
        ax = axs[1, 0]
        ax.plot(ideal_p["Temperature [K]"], ideal_p["VapThermalConductivity [W/m-K]"], label="Ideal Gas", linestyle='-', linewidth=2)
        ax.plot(real_p["Temperature [K]"], real_p["VapThermalConductivity [W/m-K]"], label="Real Gas", linestyle='--', linewidth=2)
        ax.plot(cool_p["Temperature [K]"], cool_p["VapThermalConductivity [W/m-K]"], label="CoolProp", linestyle='-.', linewidth=2)
        ax.set_title("Przewodność cieplna ($k$) pary wodnej")
        ax.set_xlabel("Temperatura [K]")
        ax.set_ylabel("Przewodność [W/(m K)]")
        ax.grid(True, linestyle='--', alpha=0.7)
        ax.legend()
        
        # Wykres 4: Lepkość dynamiczna
        ax = axs[1, 1]
        ax.plot(ideal_p["Temperature [K]"], ideal_p["VapViscosity [Pa·s]"], label="Ideal Gas", linestyle='-', linewidth=2)
        ax.plot(real_p["Temperature [K]"], real_p["VapViscosity [Pa·s]"], label="Real Gas", linestyle='--', linewidth=2)
        ax.plot(cool_p["Temperature [K]"], cool_p["VapViscosity [Pa·s]"], label="CoolProp", linestyle='-.', linewidth=2)
        ax.set_title("Lepkość dynamiczna ($\\mu$) pary wodnej")
        ax.set_xlabel("Temperatura [K]")
        ax.set_ylabel("Lepkość [Pa·s]")
        ax.grid(True, linestyle='--', alpha=0.7)
        ax.legend()
        
        plt.tight_layout()
        filename = f'Raport_Modele_Cisnienie_{p}bar.png'
        plt.savefig(filename, dpi=300)
        plt.close()
        print(f"Wygenerowano wykres: {filename}")

# ==============================================================================
# CZĘŚĆ 2: CURVE FITTING (Wyliczanie stałych dla pc_database.f90)
# ==============================================================================

# Definicje równań z CAS solver (fortran: pc_func.f90)
def eqn100(T, A, B, C, D, E):
    return A + B*T + C*T**2 + D*T**3 + E*T**4

def eqn101(T, A, B, C, D, E):
    return np.exp(A + B/T + C*np.log(T) + D * T**E)

def eqn102(T, A, B, C, D, E):
    return (A * T**B) / (1 + C/T + D/(T**2))

def eqn106(T, A, B, C, D, E):
    Tcrit = 647.096
    Tr = T / Tcrit
    return A * (1 - Tr)**(B + C*Tr + D*Tr**2 + E*Tr**3)

def dopasuj_stale():
    try:
        df_liq = pd.read_csv("water_properties.csv", encoding='cp1252')
    except FileNotFoundError:
        print("Brak pliku 'water_properties.csv'. Uruchom najpierw skrypt dla cieczy.")
        return

    # Dane wejściowe
    T_liq = df_liq['Temperature [K]'].values
    
    # Przeliczamy jednostki dla CAS
    MW = 18.01528
    rho_mol = df_liq['liqDensity [kg/m3]'].values / MW 
    Hvap_mol = df_liq['HeatOfVap [kJ/kg]'].values * 1000.0 * MW
    cp_liq_mol = df_liq['liqHeatCapacity [J/kg/K]'].values * MW
    mu_liq = df_liq['liqViscosity [Pa·s]'].values
    pvap_Pa = df_liq['vapPressure [kPa]'].values * 1000.0
    sigma = df_liq['SurfaceTension [N/m]'].values

    print("\n--- WYNIKI CURVE-FITTINGU DLA pc_database.f90 ---")
    
    popt, _ = curve_fit(eqn100, T_liq, rho_mol, maxfev=10000)
    print("! liqDensityMol (Eqn 100)")
    print(f"pc%liqDensityMol%A = {popt[0]:.5e}_WP")
    print(f"pc%liqDensityMol%B = {popt[1]:.5e}_WP")
    print(f"pc%liqDensityMol%C = {popt[2]:.5e}_WP")
    print(f"pc%liqDensityMol%D = {popt[3]:.5e}_WP")
    print(f"pc%liqDensityMol%E = {popt[4]:.5e}_WP")
    
    popt, _ = curve_fit(eqn100, T_liq, Hvap_mol, maxfev=10000)
    print("\n! HeatOfVapMol (Eqn 100)")
    print(f"pc%HeatOfVapMol%A = {popt[0]:.5e}_WP")
    print(f"pc%HeatOfVapMol%B = {popt[1]:.5e}_WP")
    print(f"pc%HeatOfVapMol%C = {popt[2]:.5e}_WP")
    print(f"pc%HeatOfVapMol%D = {popt[3]:.5e}_WP")
    print(f"pc%HeatOfVapMol%E = {popt[4]:.5e}_WP")
    
    popt, _ = curve_fit(eqn100, T_liq, cp_liq_mol, maxfev=10000)
    print("\n! liqHeatCapacityMol (Eqn 100)")
    print(f"pc%liqHeatCapacityMol%A = {popt[0]:.5e}_WP")
    print(f"pc%liqHeatCapacityMol%B = {popt[1]:.5e}_WP")
    print(f"pc%liqHeatCapacityMol%C = {popt[2]:.5e}_WP")
    print(f"pc%liqHeatCapacityMol%D = {popt[3]:.5e}_WP")
    print(f"pc%liqHeatCapacityMol%E = {popt[4]:.5e}_WP")

    p0_101 = [70, -7000, -5, 0.0001, 1.0] 
    try:
        popt, _ = curve_fit(eqn101, T_liq, pvap_Pa, p0=p0_101, maxfev=20000)
        print("\n! vapPressure (Eqn 101)")
        print(f"pc%vapPressure%A = {popt[0]:.5e}_WP")
        print(f"pc%vapPressure%B = {popt[1]:.5e}_WP")
        print(f"pc%vapPressure%C = {popt[2]:.5e}_WP")
        print(f"pc%vapPressure%D = {popt[3]:.5e}_WP")
        print(f"pc%vapPressure%E = {popt[4]:.5e}_WP")
    except Exception as e:
        pass
        
    p0_visc = [-10, 1000, -1, 0, 1]
    try:
        popt, _ = curve_fit(eqn101, T_liq, mu_liq, p0=p0_visc, maxfev=20000)
        print("\n! liqViscosity (Eqn 101)")
        print(f"pc%liqViscosity%A = {popt[0]:.5e}_WP")
        print(f"pc%liqViscosity%B = {popt[1]:.5e}_WP")
        print(f"pc%liqViscosity%C = {popt[2]:.5e}_WP")
        print(f"pc%liqViscosity%D = {popt[3]:.5e}_WP")
        print(f"pc%liqViscosity%E = {popt[4]:.5e}_WP")
    except Exception as e:
        pass

    p0_106 = [0.2358, 1.256, 0, 0, 0]
    popt, _ = curve_fit(eqn106, T_liq, sigma, p0=p0_106, maxfev=10000)
    print("\n! SurfaceTension (Eqn 106)")
    print(f"pc%SurfaceTension%A = {popt[0]:.5e}_WP")
    print(f"pc%SurfaceTension%B = {popt[1]:.5e}_WP")
    print(f"pc%SurfaceTension%C = {popt[2]:.5e}_WP")
    print(f"pc%SurfaceTension%D = {popt[3]:.5e}_WP")
    print(f"pc%SurfaceTension%E = {popt[4]:.5e}_WP")

    # Wyliczanie parametrow transportowych pary (Eqn 102) z danych CoolProp
    try:
        df_vap = pd.read_csv("Vapour_Properties_CoolProp_Results.csv", encoding='cp1252')
        # Filtrowanie czystej fazy gazowej dla 1 bar
        df_vap_1bar = df_vap[(df_vap["Pressure [bar]"] == 1.0) & (df_vap["Temperature [K]"] >= 375)]
        
        T_vap = df_vap_1bar['Temperature [K]'].values
        mu_vap = df_vap_1bar['VapViscosity [Pa·s]'].values
        k_vap = df_vap_1bar['VapThermalConductivity [W/m-K]'].values

        # Wstepne wektory szacunkowe
        p0_102 = [1e-6, 1.5, 100, 0, 0]

        popt, _ = curve_fit(eqn102, T_vap, mu_vap, p0=p0_102, maxfev=20000)
        print("\n! vapViscosity (Eqn 102)")
        print(f"pc%vapViscosity%A = {popt[0]:.5e}_WP")
        print(f"pc%vapViscosity%B = {popt[1]:.5e}_WP")
        print(f"pc%vapViscosity%C = {popt[2]:.5e}_WP")
        print(f"pc%vapViscosity%D = {popt[3]:.5e}_WP")
        print(f"pc%vapViscosity%E = {popt[4]:.5e}_WP")

        popt, _ = curve_fit(eqn102, T_vap, k_vap, p0=p0_102, maxfev=20000)
        print("\n! vapThermalConductivity (Eqn 102)")
        print(f"pc%vapThermalConductivity%A = {popt[0]:.5e}_WP")
        print(f"pc%vapThermalConductivity%B = {popt[1]:.5e}_WP")
        print(f"pc%vapThermalConductivity%C = {popt[2]:.5e}_WP")
        print(f"pc%vapThermalConductivity%D = {popt[3]:.5e}_WP")
        print(f"pc%vapThermalConductivity%E = {popt[4]:.5e}_WP")

    except Exception as e:
        print("\n! Blad dopasowania wlasciwosci oparow:", e)

if __name__ == "__main__":
    analizuj_i_rysuj_modele()
    dopasuj_stale()