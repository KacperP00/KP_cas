import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import minimize

# Wczytanie zaktualizowanych danych eksperymentalnych (separator to srednik, przecinek jako znak dziesietny)
df_exp = pd.read_csv('article_PEN54.csv', sep=';', decimal=',', header=None, names=['Time_exp', 'LPL_exp'])
df_exp = df_exp.sort_values('Time_exp').reset_index(drop=True)

# Wczytanie danych z symulacji
cols = ['step', 'time', 'LPL', 'VPL', 'U_tip', 'D_tip', 'Chi_st', 'Extra']
df_sim = pd.read_csv('spray_PL.out', delim_whitespace=True, comment='#', names=cols)

# Funkcja optymalizacyjna (minimalizacja bledu sredniokwadratowego dla LPL > 5mm)
def objective(dt):
    f_interp = interp1d(df_sim['time'] + dt[0], df_sim['LPL'], bounds_error=False, fill_value='extrapolate')
    LPL_sim_interp = f_interp(df_exp['Time_exp'])
    mask = df_exp['LPL_exp'] > 5.0
    return np.mean((LPL_sim_interp[mask] - df_exp['LPL_exp'][mask])**2)

# Obliczenie optymalnego przesuniecia czasu
dt_opt = minimize(objective, [0.3], method='Nelder-Mead').x[0]
time_shifted = df_sim['time'] + dt_opt

# --- WYKRES 1: LPL vs czas (bez przesuniecia) ---
plt.figure(figsize=(10, 6))
plt.plot(df_sim['time'], df_sim['LPL'], label='LPL (Symulacja)', color='blue', linewidth=2)
# Rysowanie linii ze znacznikami dla danych eksperymentalnych
plt.plot(df_exp['Time_exp'], df_exp['LPL_exp'], label='LPL (Eksperyment)', color='red', marker='o', linestyle='-', zorder=5)
plt.title('Liquid penetration - porównanie - Brak przesuniecia')
plt.xlabel('Czas [ms]')
plt.ylabel('LPL [mm]')
plt.xlim(0.0, 1.2)
plt.ylim(0.0, 50.0)
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend()
plt.tight_layout()
plt.savefig('wykres1_bez_przesuniecia.png', dpi=300)
plt.close()

# --- WYKRES 2: Z przesunieciem czasowym i podwojna osia Y ---
fig, ax1 = plt.subplots(figsize=(12, 7))

# Lewa os Y (Dystans w mm)
ax1.set_xlabel('Czas [ms]', fontsize=12)
ax1.set_ylabel('Dystans [mm]', fontsize=12)
l1 = ax1.plot(time_shifted, df_sim['LPL'], label='LPL (Symulacja)', color='tab:blue', linewidth=2)
l2 = ax1.plot(time_shifted, df_sim['VPL'], label='VPL (Symulacja)', color='tab:green', linewidth=2, linestyle='--')
s1 = ax1.plot(df_exp['Time_exp'], df_exp['LPL_exp'], label='LPL (Eksperyment)', color='tab:red', marker='o', linestyle='-', zorder=5)

# Sztywne limity osi dla poprawnego wyswietlania
ax1.set_xlim(0.0, 1.2)
ax1.set_ylim(0.0, 50.0)
ax1.grid(True, linestyle='--', alpha=0.7)

# Prawa os Y (Predkosc w m/s)
ax2 = ax1.twinx()  
ax2.set_ylabel('Predkosc [m/s]', color='tab:purple', fontsize=12)  
l3 = ax2.plot(time_shifted, df_sim['U_tip'], label='U_tip (Symulacja)', color='tab:purple', linewidth=2, linestyle='-.')
ax2.set_ylim(0.0, 70.0)
ax2.tick_params(axis='y', labelcolor='tab:purple')

# Wspolna legenda
lines = l1 + l2 + l3 + s1
labels = [l.get_label() for l in lines]
ax1.legend(lines, labels, loc='upper left', framealpha=0.9)

plt.title(f'Dynamika strugi (Przesuniecie dt = {dt_opt:.3f} ms)')
fig.tight_layout()  
plt.savefig('wykres2_z_przesunieciem.png', dpi=300)
plt.close()