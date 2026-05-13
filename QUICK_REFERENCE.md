# CAS Spray Solver - Quick Reference Guide

---

## INPUTS

### Fuel & Fluid Properties
```
Fuel : npentane | n-dodecane | diesel | n-hexane | n-octane | ethanol | ...
Fuel temperature : 363.15 [K]
Fuel molecular weight : auto-loaded [kg/kmol]
Fuel molecular volume : 0.000115220 [m³/kmol]
```

### Nozzle & Injector Geometry
```
Nozzle exit diameter : 200e-6 [m]
Number of nozzles : 1 [adjustable: 1,2,4,6,8,12...]
Nozzle LD ratio : (optional) [length/diameter]
Nozzle rD ratio : (optional) [radius/diameter]
Nozzle sac diameter : (optional) [m]
Nozzle K inlet : (optional) [pressure loss coefficient]
Nozzle K expansion : (optional) [pressure loss coefficient]
Discharge coefficient of nozzle : 0.6 [-]
```

### Injection Conditions
```
Injection pressure : 150e5 [Pa]
Use nozzle flow model : .true. | .false.
In nozzle cavitation : .true. | .false.
Constant injection velocity : [m/s] (optional alternative)
Injected mass : [kg] (optional alternative)
Area of nozzle : [m²] (for injected mass option)
ROI profile : roi.txt [file path - Rate of Injection vs time]
```

### Ambient Conditions
```
Ambient pressure : 1.0e5 [Pa]
Ambient temperature : 298.15 [K]
Ambient molecular weight : 28.97e-3 [kg/kmol]
Compressibility factor : 1.0 [-]
```

### Initial Droplet Size Distribution
```
Initial DSD : Delta | Log-Normal | Rosin-Rammler | Gamma | GeneralizedGamma | MEF | Blob
Number of diameter classes : 1 [typically 1]
Initial Dm : 1.0e-5 [m] **← MUST SET**
Initial D2 : [m²] (optional; if not set, computed from Dvar)
Initial D3 : [m³] (optional)
Initial Dvar : 0.0 [m²]
Droplet grid factor : 5.0 [-]
```

### Physical Models (Activation Flags)
```
Use entrainment model : .true. | .false.
Use drag model : .true. | .false.
Use breakup model : .true. | .false.
Use evaporation model : .true. | .false.
Turbulence model : .true. | .false.
Flash boiling : .true. | .false.
```

### Spray Angle & Spreading
```
Spray half-cone angle model : Constant | Reitz-Bracco | Reitz-Bracco-Simplified | 
                              Ruiz-Chigier | Arai | Hiroyasu-Arai | Arregle | 
                              Siebers | Reddemann | Kitamura
Spray half-cone angle : [°] (if Constant model)
Tangent of spray half-cone angle : [alternative input]
C_theta : [model coefficient]
```

### Breakup Model Parameters
```
B0 : 0.61 [-] Kelvin-Helmholtz wave amplitude
B1 : 5.0 [-] KH wave damping coefficient
C3 : 0.45 [-] Breakup time constant
D0 : 1.0 [-] Rayleigh-Taylor primary
D1 : 1.0 [-] Rayleigh-Taylor secondary
Crel : 1.0 [-] Relaxation coefficient
```

### Evaporation Parameters
```
Cevap : 0.3333 [-] Evaporation rate multiplier
A0 : 1.0 [-] Heat transfer correction
A1 : 0.85 [-] Heat transfer correction
```

### Turbulence Parameters
```
Turbulence model : .true. | .false.
C_k type : Dynamic | Const | Constant
C_k : [value if Constant]
Y_O2 : 0.164215 [-] Oxygen mass fraction
```

### Combustion Model
```
Combustion model : none | MRIF
C_chi : 1.0 [-] Scalar dissipation coefficient
Flamelet interval : 1.0E-04 [s]
Number of flamelets : [integer]
Use 2D profiles for convolution : .true. | .false.
```

### Numerical Parameters - Grid & Mesh
```
Number of grid points : 1601 [axial points]
Length : 1600.0e-6 [m] domain length
Number of radial grid points : [optional, 2D]
Maximum CFL : 0.1 [-] Courant condition
```

### Solver Parameters
```
RK scheme : 1 | 2 | 3 | 4 (or Euler | TVDRK-22 | TVDRK-33 | TVDRK-44)
Conv scheme : UPWIND1 | LF | WENO3 | WENO5
Strang splitting : .true. | .false.
Maximum time : 1.12e-3 [s]
Tolerance for Newton solver : [tolerance]
Maximum iterations for Newton solver : [count]
Initial relaxation coefficient for Newton solver : [value]
Alpha for Newton solver : [value]
```

### Output Control
```
Output frequency : 10 [steps] write spray_PL.out
Datafile frequency : 500 [steps] write data.out_*, omega.out_*
Data File : data.out [prefix]
PL File : spray_PL.out [filename]
Saturation PT profile : pentane_sat.txt [file path]
```

### Optional/Advanced Inputs
```
Fixed Re : -9999 [force Reynolds number]
Fixed We : -9999 [force Weber number]
Fixed DRa : -9999 [force density ratio liquid/ambient]
Fixed DRv : -9999 [force density ratio liquid/vapor]
Fixed VRa : -9999 [force viscosity ratio]
Fixed VRv : -9999 [force viscosity ratio]
Fixed De : -9999 [force evaporation number]

f_rho_l : 0.0 [liquid density sensitivity factor]
f_mu_l : 0.0 [liquid viscosity sensitivity]
f_lambda_l : 0.0 [liquid thermal conductivity]
f_Lv : 0.0 [latent heat]
f_C_l : 0.0 [liquid heat capacity]
f_Pv : 0.0 [vapor pressure]
f_sigma : 0.0 [surface tension]
f_rho_v : 0.0 [vapor density]
f_mu_v : 0.0 [vapor viscosity]
f_lambda_v : 0.0 [vapor conductivity]
f_Cp_v : 0.0 [vapor heat capacity]
```

---

## OUTPUTS

### Primary Output: `spray_PL.out`
| Column | Variable | Units | Description |
|--------|----------|-------|-------------|
| 1 | step | - | Timestep number |
| 2 | time | s | Simulation time |
| 3 | **LPL** | m | **Liquid Penetration Length** |
| 4 | **VPL** | m | **Vapor Penetration Length** |
| 5 | U_tip | m/s | Tip velocity |
| 6 | D_tip | m | Tip diameter |
| 7 | Chi_st | 1/s | Scalar dissipation |
| 8 | Chi_st1 | 1/s | Scalar dissipation variant |

### Field Data: `data.out_XXXXXX`
| Col | Variable | Symbol | Units | Key Info |
|-----|----------|--------|-------|----------|
| 1 | z | - | m | Axial position |
| 2 | **rho** | ρ | kg/m³ | Mixture density |
| 3 | **Y_l** | Y_l | - | **LIQUID FRACTION** |
| 4 | Y_v | Y_v | - | Vapor fraction |
| 5 | Y_a | Y_a | - | Air fraction |
| 6 | Y_g | Y_g | - | Total gas (v+a) |
| 7 | u_l | u_l | m/s | Liquid velocity |
| 8 | u_g | u_g | m/s | Gas velocity |
| 9 | **dm** | d_m | m | **MEAN DROPLET DIAMETER** |
| 10 | dvar | σ_d² | m² | Droplet variance |
| 11 | d2 | - | m² | Second moment |
| 12 | d3 | - | m³ | Third moment |
| 13 | **Td** | T_d | - | **DROPLET TEMPERATURE** |
| 14 | Tg | T_g | - | Gas temperature |
| 15 | **b** | b | m | **JET WIDTH** |
| 16 | k_g | k | J/kg | Turbulent k.e. |
| 17 | eps_g | ε | J/(kg·s) | Dissipation rate |
| 18 | mu_t_g | μ_t | Pa·s | Turbulent viscosity |
| 19 | zvar_g | Z_var | - | Mixture fraction variance |
| 20 | **zmix_g** | **Z_mix** | - | **MIXTURE FRACTION** |
| 21 | **chi_g** | **χ** | 1/s | **SCALAR DISSIPATION** |
| 22 | chi_g_stl | χ_st | 1/s | Stoich. dissipation |
| 23 | dsd_type | - | - | Distribution type code |
| 24 | Pr_g | Pr | - | Prandtl number |
| 25 | Sc_g | Sc | - | Schmidt number |
| 26 | c_k | c_k | - | Turbulence model coeff |
| 27 | VRg | - | - | Viscosity ratio |
| 28 | VRtg | - | - | Turbulent visc. ratio |
| 29 | SR | - | - | Surface tension ratio |
| 30 | DRl | - | - | Liquid density ratio |
| 31 | Y_ref | - | - | Reference mass fraction |
| 32 | alpha | α | - | Profile shape parameter |
| 33 | T_v | - | K | Vapor temperature |
| 34 | DRg | - | - | Liquid/gas density ratio |

### Source Terms: `omega.out_XXXXXX`
| Col | Variable | Symbol | Units | Description |
|-----|----------|--------|-------|-------------|
| 1 | z | - | m | Axial position |
| 2 | **omega_ent** | - | kg/(m³·s) | **ENTRAINMENT SOURCE** |
| 3 | omega_drag | - | N/m³ | Drag force |
| 4 | **omega_bre1** | - | 1/s | **PRIMARY BREAKUP** |
| 5 | omega_bre2 | - | 1/s | Secondary breakup |
| 6 | **omega_vap** | - | kg/(m³·s) | **EVAPORATION SOURCE** |
| 7 | omega_T | - | K/s | Temperature change |
| 8 | omega_k_g_p | - | J/(kg·s) | k-equation production |
| 9 | omega_k_g_d | - | J/(kg·s) | k-equation dissipation |
| 10 | omega_eps_g_p | - | - | ε-equation production |
| 11 | omega_eps_g_d | - | - | ε-equation dissipation |
| 12 | omega_zvar_g_p | - | - | Z_var production |
| 13 | omega_zvar_g_d | - | - | Z_var dissipation |

### Droplet Size Distribution: `dsd.out_XXXXXX`
```
z [m]  dsdlam(1) dsdlam(2) dsdlam(3) dsdlam(4)  dm [m]  d2 [m²]  d3 [m³]
```

### Simulation Parameters: `simulation_param.out`
```
Re = Reynolds number
We = Weber number
DRa = Density ratio liquid/ambient
VRa = Viscosity ratio liquid/ambient
DRv = Density ratio liquid/vapor
VRv = Viscosity ratio liquid/vapor
WR = Molecular weight ratio
De = Dimensionless evaporation number
Spray half-cone angle = [°]
Spreading coefficient = [-]
Discharge coefficient of nozzle = [-]
Injection velocity = [m/s]
Effective jet diameter = [m]
Max. Mass flow rate = [mg/ms]
Max. Momentum flow rate = [kgm/s²]
Stoichiometric mixture fraction = [-]
```

### Optional 2D Profiles: `Zmix2D.out_*`, `Zvar2D.out_*`
```
Radial profiles of mixture fraction and variance
```

---

## POSSIBILITIES / ADJUSTABLE PARAMETERS

### BY IMPACT ON PENETRATION LENGTH (LPL):

**🔴 RANK 1: INJECTION PRESSURE** (±30% effect)
```
100 bar  → LPL reduced ~30%
150 bar  → LPL baseline (your case)
200 bar  → LPL increased ~25%
300 bar  → LPL increased ~45%
```
**How to adjust**: `Injection pressure : XXXe5`

---

**🟠 RANK 2: NOZZLE DIAMETER** (±25% effect)
```
100 µm  → LPL half (D_noz halves jet scale)
200 µm  → LPL baseline (your case)
300 µm  → LPL increased ~50%
```
**How to adjust**: `Nozzle exit diameter : XXXe-6`

---

**🟡 RANK 3: BREAKUP COEFFICIENTS** (±15% effect)
```
B1 = 2.5   → Aggressive breakup → smaller dm → LPL ↑ 10%
B1 = 5.0   → Baseline (your case)
B1 = 10.0  → Mild breakup → larger dm → LPL ↓ 10%

C3 = 0.2   → Fast breakup → LPL ↑ 8%
C3 = 0.45  → Baseline
C3 = 1.0   → Slow breakup → LPL ↓ 5%
```
**How to adjust**: 
```
B0 : 0.61
B1 : 5.0
C3 : 0.45
```

---

**🟢 RANK 4: INITIAL DROPLET SIZE** (±10% effect)
```
Dm = 5 µm    → Very fine → rapid evap → LPL ↓ 10%
Dm = 10 µm   → Baseline (your case)
Dm = 20 µm   → Coarse → slow evap → LPL ↑ 8%
```
**How to adjust**: `Initial Dm : XXXe-6`

---

**🔵 RANK 5: NUMBER OF NOZZLES** (±20% per-hole effect)
```
1 hole   → Full momentum/mass on centerline
4 holes  → 1/4 mass per hole → LPL ↓ ~20% per hole
6 holes  → Reduced circumferential jets
Note: CAS averages all holes; doesn't resolve individual jets
```
**How to adjust**: `Number of nozzles : N`

---

**🟣 RANK 6: EVAPORATION RATE** (±5% effect)
```
Cevap = 0.1  → Suppress evaporation
Cevap = 0.33 → Conservative (your case)
Cevap = 1.0  → Aggressive
Cevap = 2.0  → Very fast → VPL >> LPL (flash-boiling)
```
**How to adjust**: `Cevap : X.XXX`

---

**🟡 RANK 7: SPRAY ANGLE** (±10% effect)
```
θ = 10°  → Narrow cone → LPL ↑ 5%; concentrated
θ = 25°  → Baseline (Kitamura model, your case)
θ = 35°  → Wide cone → LPL ↓ 8%; dilute
```
**How to adjust**:
```
Spray half-cone angle model : Constant
Spray half-cone angle : 25.0
```

---

**🟠 RANK 8: FUEL TEMPERATURE** (±5% effect)
```
300 K   → Cool fuel
363 K   → Baseline (your case)
400 K   → Hot → increased evaporation
```
**How to adjust**: `Fuel temperature : XXX.XX`

---

### BY PHYSICS ASPECT:

**To Study Breakup Dynamics**:
- Vary B0, B1, C3
- Monitor **dm vs z** in `data.out_*`
- Plot **omega_bre1, omega_bre2** from `omega.out_*`

**To Study Evaporation**:
- Vary Cevap, T_fuel, Injection pressure
- Monitor **Y_v, Td, Y_l vs z** from `data.out_*`
- Track **Y_v penetration (VPL)** from `spray_PL.out`

**To Study Entrainment/Mixing**:
- Vary spray angle, turbulence model, grid resolution
- Monitor **zmix_g, zvar_g, chi_g** from `data.out_*`
- Plot **omega_ent** profile from `omega.out_*`

**To Study Flash-Boiling** (T_fuel + P_inj combinations):
- Set `Flash boiling : .true.`
- Increase T_fuel (e.g., 363 K → 400 K)
- Keep P_inj fixed or vary (e.g., 150 bar)
- Compare VPL/LPL ratio: >1 = vapor reaches ahead

**To Study Multi-Hole Injectors**:
- Set `Number of nozzles : 6` (or other count)
- Result: mass/flow per hole reduced proportionally
- Limitation: Individual jet spacing NOT modeled

**To Study Sensitivity**:
- Vary one parameter at a time
- Keep others fixed
- Report LPL, dm, zmix_g metrics

---

## MODEL SELECTION QUICK CHOICES

### Default (Recommended)
```ini
Spray half-cone angle model : Kitamura
RK scheme : 3
Conv scheme : WENO5
Use entrainment model : .true.
Use drag model : .true.
Use breakup model : .true.
Use evaporation model : .true.
Turbulence model : .false.
Flash boiling : .true.
```

### For Fast Testing
```ini
Number of grid points : 801
RK scheme : 1
Conv scheme : WENO3
Output frequency : 50
Datafile frequency : 5000
```

### For High Accuracy (Research)
```ini
Number of grid points : 3201
RK scheme : 4
Conv scheme : WENO5
Maximum CFL : 0.05
Output frequency : 5
Datafile frequency : 50
```

### For Combustion Coupling
```ini
Combustion model : MRIF
Turbulence model : .true.
C_k type : Dynamic
Y_O2 : 0.164215
Use 2D profiles for convolution : .true.
```

### For Flash-Boiling Focus
```ini
Flash boiling : .true.
Cevap : 0.5 to 2.0 (sweep parameter)
Fuel temperature : 363.15 to 400.0 K
Saturation PT profile : pentane_sat.txt
```

---

## AVAILABLE FUELS (Database)

```
Alkanes: diesel, n-pentane, n-hexane, n-heptane, n-octane, n-nonane, 
         n-decane, n-dodecane, n-tetradecane, iso-octane, iso-pentane
Alcohols: methanol, ethanol, 1-propanol, 2-propanol, 1-butanol, 2-butanol
Ethers: dimethylether, diethylether, OME1, OME4, dinpropylether, dinbutylether
Oxygenates: acetone, methyl-ethyl-ketone, γ-valerolactone
Others: methane, ethane, propane, isobutane
```

---

## FILE NAMING CONVENTION

```
Input:       cas.in
Output:      spray_PL.out                    (time-series, appended)
             data.out_000000, data.out_000100, ...  (snapshots, one per Datafile freq)
             omega.out_000000, omega.out_000100, ...
             dsd.out_XXXXXX, bpdf.out_XXXXXX        (if enabled)
             Zmix2D.out_XXXXXX, Zvar2D.out_XXXXXX  (if 2D)
             simulation_param.out                  (written once)
ROI curve:   roi.txt
P-T-sat:     pentane_sat.txt (optional)
```

---

## EXECUTION

```bash
cd cas/src && make all_opt
cd cas/example/'Aleiferis injector'
./cas_opt
```

Output files appear in current directory.

---

## EXTRACTION & PLOTTING

```python
import pandas as pd
import matplotlib.pyplot as plt

# Read spray_PL.out
data = pd.read_csv('spray_PL.out', delim_whitespace=True, comment='#', header=None)
time_ms = data[1]*1000
lpl_mm = data[2]*1000
vpl_mm = data[3]*1000

plt.plot(time_ms, lpl_mm, label='LPL')
plt.plot(time_ms, vpl_mm, label='VPL')
plt.xlabel('Time (ms)')
plt.ylabel('Penetration (mm)')
plt.legend()
plt.show()

# Read data.out_000000
fields = pd.read_csv('data.out_000000', delim_whitespace=True, comment='#', header=None, skiprows=1)
z_mm = fields[0]*1000
y_l = fields[2]
dm_um = fields[8]*1e6
td = fields[12]
zmix = fields[19]

fig, axes = plt.subplots(2,2)
axes[0,0].plot(z_mm, y_l); axes[0,0].set_ylabel('Y_l')
axes[0,1].plot(z_mm, dm_um); axes[0,1].set_ylabel('dm [µm]')
axes[1,0].plot(z_mm, td); axes[1,0].set_ylabel('Td')
axes[1,1].plot(z_mm, zmix); axes[1,1].set_ylabel('Zmix')
for ax in axes.flat: ax.set_xlabel('z [mm]')
plt.tight_layout()
plt.show()
```

---

## KEY METRICS TO REPORT

- **Liquid Penetration @ 1 ms** [mm]
- **Tip Velocity (average)** [m/s]
- **Mean Droplet Diameter @ z=50mm** [µm]
- **Mixture Fraction (centerline @ penetration tip)** [-]
- **Non-dimensional: Re, We, DRa** [-]
- **Spray Angle** [°]

---

## END
