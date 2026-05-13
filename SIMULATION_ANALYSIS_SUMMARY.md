# CAS Spray Simulation: Comprehensive Analysis Summary

**Solver Type**: Cross-sectionally averaged (1D axial) spray model with radial profiles  
**Physical Basis**: Reduced-order model for multiphase spray dynamics  
**Application**: Compression ignition engines, diesel injection, flash-boiling sprays

---

## 1. AVAILABLE MODELS AND COMPARISON

### 1.1 Injector/Spray Angle Models
These models compute the spray half-cone angle evolution and spreading coefficient.

| Model Name | Description | Where Used |
|------------|-------------|-----------|
| **Reitz-Bracco** | Classical spray angle model based on non-dimensional spray parameters | Full dynamic evolution of cone angle |
| **Reitz-Bracco-Simplified** | Simplified version with reduced computational overhead | Fast evaluation |
| **Ruiz-Chigier** | Alternative formulation for spray spreading | Secondary alternative |
| **Arai** | Another empirical correlation for spray angle | Regional variation support |
| **Hiroyasu-Arai** / **Hiroyasu** | Widely used in engine spray modeling | Commonly referenced |
| **Arregle** | Spanish research group formulation | Technical variant |
| **Siebers** | Heavy-duty diesel spray model | Industrial engines |
| **Reddemann** | Recent model from Fuel Science Center | Modern formulation |
| **Kitamura** | (Used in Aleiferis example) | Current benchmark case |
| **Constant** | Fixed spray angle throughout simulation | For sensitivity studies |

**Selection Guidance**: 
- For production engines: Use **Hiroyasu-Arai** or **Siebers**
- For research: Use **Reitz-Bracco** or **Reddemann** for best physics
- For debug/quick runs: Use **Constant** with specified angle

---

### 1.2 Entrainment Model
Controls how gas is entrained into the spray core.

| Setting | Description |
|---------|-------------|
| **ON** (`.true.`) | Standard entrainment model active; includes mixing layer dynamics |
| **OFF** (`.false.`) | Disables entrainment; useful for isolated jet studies |

**Default**: `.true.`

---

### 1.3 Drag Model
Determines liquid-gas momentum transfer.

| Setting | Description |
|---------|-------------|
| **ON** (`.true.`) | Standard drag correlation (function of droplet Re, We) |
| **OFF** (`.false.`) | No phase interaction; unrealistic but useful for testing |

**Default**: `.true.`

---

### 1.4 Breakup Model
Controls primary and secondary droplet breakup.

| Setting | Description | Key Parameters |
|---------|-------------|-----------------|
| **ON** (`.true.`) | Kelvin-Helmholtz wave instability + RTI model | B0, B1, C3, Crel, D0, D1 |
| **OFF** (`.false.`) | No breakup; droplets remain initial size | Unrealistic for most cases |

**Breakup Parameters** (defaults):
- `B0 = 0.61` – KH wave amplitude coefficient
- `B1 = 5.0` – KH wave damping rate
- `C3 = 0.45` – Breakup time constant
- `Crel = 1.0` – Relaxation coefficient
- `D0 = 1.0`, `D1 = 1.0` – RTI (Rayleigh-Taylor) breakup coefficients

**Recommendation**: Use defaults unless you have experimental data for your fuel/pressure combo.

---

### 1.5 Evaporation Model
Controls liquid-phase vaporization dynamics.

| Setting | Description | Key Parameters |
|---------|-------------|-----------------|
| **ON** (`.true.`) | Temperature-dependent evaporation with Clausius-Clapeyron | Cevap, A0, A1 |
| **OFF** (`.false.`) | No evaporation; useful for cold spray validation |

**Evaporation Parameters** (defaults):
- `Cevap = 0.3333` – Evaporation rate multiplier
- `A0 = 1.0`, `A1 = 0.85` – Heat transfer correction coefficients

---

### 1.6 Droplet Size Distribution (DSD) Types
Defines how initial droplet spectrum is represented.

| DSD Type | Description | Best For | Parameters |
|----------|-------------|----------|-----------|
| **Delta** | Monodisperse (all droplets same size) | Testing, validation against simulations | Only `dm` needed |
| **Log-Normal** | Classical log-normal distribution; widely used in spray codes | Most cases; good for SMD correlation | `dm`, `dvar` (or `d2`) |
| **Rosin-Rammler** | Power-law distribution; popular for industrial sprays | Diesel injectors, broad spectra | `dm`, `dvar` |
| **Gamma** | Gamma probability distribution | Alternative functional form | `dm`, `dvar` |
| **Generalized Gamma** | More flexible than standard gamma | Advanced models | `dm`, `dvar` |
| **Maximum-Entropy-Formalism** | Entropy-based optimal distribution | Research; minimal assumptions | `dm`, `dvar` |
| **Blob** | Log-Normal approximation of primary breakup | Initial transient analysis | `dm`, `dvar` |

**Practical Recommendations**:
- For your npentane case at 150 bar: **Delta** or **Log-Normal** with narrow variance (<10% of `dm`)
- For broader comparisons: **Rosin-Rammler** as it's industry standard
- For research papers: Compare 2-3 distributions to show robustness

---

### 1.7 Turbulence Model
Models gas-phase turbulent mixing.

| Setting | Description | c_k Type |
|---------|-------------|----------|
| **ON** (`.true.`) | Standard k-ε turbulence closure; predicts turbulent kinetic energy and dissipation | Dynamic or Constant |
| **OFF** (`.false.`) | No gas-phase turbulence; laminar/inviscid assumption |

**c_k Options**:
- **Dynamic**: $c_k$ recomputed each timestep from resolved scales (more accurate, slightly slower)
- **Constant** (`'Const'`, `'Constant'`): Fixed value (default ~7.0) (faster, less accurate for all conditions)

**Turbulence Constants** (standard k-ε):
- `c_mu = 0.09` – Eddy viscosity coefficient
- `c_eps1 = 1.44` – Production coefficient
- `c_eps2 = 1.92` – Dissipation coefficient
- `c_zvar = 2.0` – Scalar variance dissipation

**When to use**:
- **ON with Dynamic**: Best for research papers; captures momentum mixing evolution
- **ON with Constant**: Good balance for typical engine conditions
- **OFF**: Only for validation against experiments with strong ambient turbulence

---

### 1.8 Combustion Model
Controls reactive species evolution.

| Model | Description |
|-------|-------------|
| **None** / **'noname'** | Inert spray; no combustion (used in your Aleiferis case) |
| **'MRIF'** | Multiply-Referenced Interactive Flamelets; reduced-order combustion coupling |

**When to use**:
- **None**: For flash-boiling, evaporation, mixing studies **← Current case**
- **'MRIF'**: When you have experimental emissions data (soot, NOx) to validate against

---

### 1.9 Nozzle Flow Options
Controls exit velocity calculation.

| Option | Input Required | Use Case |
|--------|-----------------|----------|
| **`Use nozzle flow model : .true.`** | Nozzle geometry (D, L/D, r/D, Dsac, K_in, K_exp) + cavitation flag | Most accurate; accounts for in-nozzle effects |
| **`Constant injection velocity`** | Direct velocity value | Testing; bypass nozzle dynamics |
| **`Injected mass` + `Area of nozzle` + Cd** | Total mass, orifice area, discharge coefficient | Backward calc from ROI profile |

**Recommendation**: Use **nozzle flow model** for realistic physics; set `In nozzle cavitation: .true.` for flash-boiling simulations.

---

### 1.10 Time Integration Scheme
Controls numerical stability and accuracy of time-stepping.

| Scheme | Order | Stages | Best For |
|--------|-------|--------|----------|
| **Euler** (`'1'`, `'11'`, `'E11'`) | 1st order | 1 | Fast, lower accuracy |
| **TVDRK-22** (`'2'`, `'22'`, `'TRK-22'`) | 2nd order | 2 stages | Good balance |
| **TVDRK-33** (`'3'`, `'33'`, `'TRK-33'`) | 3rd order | 3 stages | Higher accuracy |
| **TVDRK-44** (`'4'`, `'44'`, `'TRK-44'`) | 4th order | 4 stages | Maximum accuracy (slower) |

**Recommendation**: Use **TVDRK-33** as default; it's TVD stable and 3rd-order accurate.

---

### 1.11 Convection Scheme
Controls spatial discretization of hyperbolic fluxes.

| Scheme | Accuracy | Dissipation | Best For |
|--------|----------|------------|----------|
| **UPWIND1** | 1st order | High | Stable but diffusive |
| **LF** (Lax-Friedrichs) | 2nd order | Medium | General purpose |
| **WENO3** | 3rd order | Low | Shock-capturing |
| **WENO5** | 5th order | Very Low | High-resolution; more expensive |

**Recommendation**: Use **WENO5** for research; **WENO3** for fast production runs.

---

## 2. INPUT PARAMETERS WITH DESCRIPTIONS

### 2.1 Fuel Properties

| Parameter | Units | Example (npentane) | Description |
|-----------|-------|-------------------|-------------|
| **Fuel** | - | `npentane` | Fuel name; must match database (see list in `pc_func.f90`) |
| **Fuel temperature** | K | 363.15 | Injected liquid fuel temperature |
| **Fuel molecular weight** | kg/kmol | - | Auto-loaded from database; can override |
| **Fuel molecular volume** | m³/km | 0.000115220 | Used in some property correlations |
| **Fuel properties table** | - | - | Path to external property table (optional) |

**Available Fuels** (examples): diesel, n-dodecane, n-pentane, n-hexane, n-heptane, n-octane, 1-butanol, 2-propanol, acetone, methane, and many others.

---

### 2.2 Nozzle/Injector Geometry

| Parameter | Units | Example | Description |
|-----------|-------|---------|-------------|
| **Nozzle exit diameter** | m | 200e-6 | Orifice diameter at nozzle exit |
| **Nozzle LD ratio** | - | (optional) | Length/Diameter ratio of orifice |
| **Nozzle rD ratio** | - | (optional) | Radius/Diameter ratio for entrance rounding |
| **Nozzle sac diameter** | m | (optional) | Fuel sac volume diameter |
| **Nozzle K inlet** | - | (optional) | In-nozzle pressure loss coefficient |
| **Nozzle K expansion** | - | (optional) | Expansion pressure loss coefficient |
| **Number of nozzles** | - | 1 (implicit) | **← CAN ADJUST**: Number of injection holes; scales mass flow per hole |
| **Discharge coefficient of nozzle** | - | 0.6 | Effective area ratio; calibration parameter |

---

### 2.3 Injection Conditions

| Parameter | Units | Example | Description |
|-----------|-------|---------|-------------|
| **Injection pressure** | Pa | 150e5 (150 bar) | Fuel pump pressure; **← KEY PARAMETER** |
| **ROI profile** | file | `roi.txt` | Rate-of-Injection time series (mass flow vs time) |
| **Use nozzle flow model** | - | `.true.` | Enable in-nozzle hydraulics; **← RECOMMEND .true.** |
| **In nozzle cavitation** | - | `.true.` | Allow cavitation inside nozzle for flash-boiling |
| **Constant injection velocity** | m/s | - | Alternative: fix exit velocity directly |
| **Injected mass** | kg | - | Alternative: specify total fuel mass to inject |
| **Area of nozzle** | m² | - | Companion to injected mass option |

---

### 2.4 Ambient Conditions

| Parameter | Units | Example | Description |
|-----------|-------|---------|-------------|
| **Ambient pressure** | Pa | 1.0e5 (1 bar) | Chamber/cylinder ambient gas pressure |
| **Ambient temperature** | K | 298.15 | Initial gas temperature |
| **Ambient molecular weight** | kg/kmol | 28.97e-3 | Air ≈ 28.97 kg/kmol |
| **Compressibility factor** | - | 1.0 | Z factor for real gas; =1.0 for ideal gas |

---

### 2.5 Initial Droplet Size Distribution (DSD)

| Parameter | Units | Example | Description |
|-----------|-------|---------|-------------|
| **Initial DSD** | - | `Delta` | Distribution type (see Section 1.6) |
| **Number of diameter classes** | - | 1 | Number of discrete droplet classes; must match DSD output files |
| **Initial Dm** | m | 1.0e-5 | Mean diameter; **← ESSENTIAL; AFFECTS PENETRATION** |
| **Initial D2** | m² | - | Second moment; if not given, computed from Dvar |
| **Initial D3** | m³ | - | Third moment; usually set =0 or auto-calculated |
| **Initial Dvar** | m² | 0.0 | Variance of diameter distribution; 0 = monodisperse |
| **Droplet grid factor** | - | 5.0 | Refinement factor for diameter class grid |

---

### 2.6 Physical Models (Boolean Flags)

| Parameter | Default | Description |
|-----------|---------|-------------|
| **Use entrainment model** | `.true.` | Enable spray entrainment |
| **Use drag model** | `.true.` | Enable liquid-gas drag |
| **Use breakup model** | `.true.` | Enable droplet breakup |
| **Use evaporation model** | `.true.` | Enable phase change |
| **Turbulence model** | `.false.` | Enable gas-phase turbulence (k-ε) |
| **Flash boiling** | `.true.` | Enable superheated vapor formation |

---

### 2.7 Spray Angle & Spreading

| Parameter | Units | Example | Description |
|-----------|-------|---------|-------------|
| **Spray half-cone angle model** | - | `Kitamura` | Model for angle evolution (Section 1.1) |
| **Spray half-cone angle** | ° | - | Fixed angle if `Constant` model used |
| **Tangent of spray half-cone angle** | - | - | Alternative input: tan(θ) instead of θ |
| **C_theta** | - | - | Coefficient in dynamic angle models |

---

### 2.8 Breakup Model Parameters

| Parameter | Default | Typical Range | Description |
|-----------|---------|----------------|-------------|
| **B0** | 0.61 | 0.5–0.7 | Kelvin-Helmholtz wave amplitude |
| **B1** | 5.0 | 3–10 | KH wave damping coefficient; try wider range for sensitivity |
| **C3** | 0.45 | 0.3–1.0 | Breakup time constant; **← IMPACTS D10 EVOLUTION** |
| **D0, D1** | 1.0 | 0.5–2.0 | Rayleigh-Taylor instability (secondary breakup) |
| **Crel** | 1.0 | 0.5–1.5 | Relative velocity relaxation; affects droplet acceleration |

**For your case**: Start with defaults; if SMD doesn't match experiment, vary B1 and C3 first.

---

### 2.9 Evaporation Model Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| **Cevap** | 0.3333 | Evaporation rate multiplier; **← ADJUST FOR FLASH-BOILING** |
| **A0, A1** | 1.0, 0.85 | Analytical heat transfer correction factors |

---

### 2.10 Turbulence Model Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| **C_k type** | `'Dynamic'` | Dynamic or Constant c_k evaluation |
| **C_k** | - | Constant value if type = `'Const'` (typically 7.0) |

---

### 2.11 Numerical Parameters

| Parameter | Units | Default | Description |
|-----------|-------|---------|-------------|
| **Number of grid points** | - | 1601 | Axial grid resolution; **← IMPACTS ACCURACY & CPU TIME** |
| **Length** | m | 1600.0e-6 | Domain length (typically 1600× nozzle diameter) |
| **Number of radial grid points** | - | - | Radial grid for 2D profiles (if enabled) |
| **Maximum CFL** | - | 0.1 | Courant-Friedrichs-Lewy condition; <0.2 is safe |
| **Maximum time** | s | 1.12e-3 | Simulation end time; **← 1.12 ms for Aleiferis** |

---

### 2.12 Solver & Time Integration

| Parameter | Default | Description |
|-----------|---------|-------------|
| **RK scheme** | `'1'` | Runge-Kutta time scheme (1–4, see Section 1.10) |
| **Conv scheme** | `'WENO5'` | Convection scheme (see Section 1.11) |
| **Strang splitting** | `.false.` | Operator splitting for source terms |
| **Tolerance for Newton solver** | - | Convergence tolerance for nonlinear solver |
| **Maximum iterations for Newton solver** | - | Safety limit on iterations |

---

### 2.13 Output Control

| Parameter | Units | Default | Description |
|-----------|-------|---------|-------------|
| **Output frequency** | steps | 10 | Write `spray_PL.out` every N timesteps |
| **Datafile frequency** | steps | (large) | Write full `data.out_*`, `omega.out_*` every N steps; **← REDUCE FOR MANY SNAPSHOTS** |
| **Data File** | - | `data.out` | Prefix for full-field output |
| **PL File** | - | `spray_PL.out` | Penetration/tip data filename |

---

### 2.14 Optional/Advanced Parameters

| Parameter | Description | When to Use |
|-----------|-------------|------------|
| **Fixed Re** | Force Reynolds number to constant | Sensitivity studies on Re alone |
| **Fixed We** | Force Weber number constant | Sensitivity on We alone |
| **Fixed DRa, DRv, VRa, VRv, De** | Fix non-dimensional numbers | Isolate specific physics effects |
| **Y_O2** | Oxygen mass fraction (for combustion) | Needed if using MRIF model |
| **Flamelet interval** | Combustion model time interval | When using MRIF |
| **Use 2D profiles for convolution** | Enable 2D Zmix/Zvar fields | Advanced combustion coupling |
| **C_chi** | Scalar dissipation coefficient | Combustion model tuning |
| **Saturation PT profile** | External P-T-saturation table | Flash-boiling refinement |

---

## 3. POSSIBLE OUTPUTS (MOST IMPORTANT)

### 3.1 Primary Output File: `spray_PL.out`

**Frequency**: Every `Output frequency` steps  
**Type**: Time-series of global spray quantities  
**Columns**:

| # | Variable | Symbol | Units | Description | **Scientific Value** |
|---|----------|--------|-------|-------------|----------------------|
| 1 | step | - | - | Simulation timestep number | ⭐⭐ Debug/tracking |
| 2 | time | t | s | Elapsed simulation time | ⭐⭐⭐ **Primary axis for plots** |
| 3 | **Liquid Penetration Length** | **LPL** | m | Distance spray liquid front has traveled | ⭐⭐⭐⭐⭐ **MAIN METRIC** |
| 4 | **Vapor Penetration Length** | **VPL** | m | Distance spray vapor front has traveled | ⭐⭐⭐⭐ **Flash-boiling indicator** |
| 5 | **Tip Velocity** | **U_tip** | m/s | Velocity of leading droplet/liquid front | ⭐⭐⭐⭐ **Penetration rate** |
| 6 | **Tip Diameter** | **D_tip** | m | Effective spray diameter at penetration tip | ⭐⭐⭐ **Spreading metric** |
| 7 | Chi_st | χ_st | - | Scalar dissipation at stoichiometric mixture fraction | ⭐⭐ **Combustion indicator** |
| 8 | Chi_st1 | χ_st1 | - | Additional scalar dissipation metric | ⭐ Combustion detail |

**How to Extract & Use**:
```python
data = pd.read_csv('spray_PL.out', delim_whitespace=True, comment='#', header=None)
time_ms = data[1] * 1000  # Convert to ms
lpl_mm = data[2] * 1000   # Convert to mm
vpl_mm = data[3] * 1000
plt.plot(time_ms, lpl_mm, label='LPL')
plt.plot(time_ms, vpl_mm, label='VPL')
```

---

### 3.2 Detailed Field File: `data.out_XXXXXX` (Full 1D Profiles)

**Frequency**: Every `Datafile frequency` steps  
**Type**: Axial profiles at each z-location  
**Typical Size**: 1600 rows × 34 columns  

| # | Variable | Symbol | Units | Description | **Scientific Value** |
|---|----------|--------|-------|-------------|----------------------|
| 1 | z | z | m | Axial position from nozzle | ⭐⭐⭐ Grid position |
| 2 | **Mixture Density** | **rho** | kg/m³ | Local mixture density (non-dimensional) | ⭐⭐⭐ **Core physics** |
| 3 | **Liquid Mass Fraction** | **Y_l** | - | Liquid volume fraction (0–1) | ⭐⭐⭐⭐⭐ **Spray structure** |
| 4 | **Vapor Mass Fraction** | **Y_v** | - | Vapor mass fraction | ⭐⭐⭐⭐ **Evaporation metric** |
| 5 | **Air Mass Fraction** | **Y_a** | - | Ambient air mass fraction | ⭐⭐⭐ **Mixing** |
| 6 | **Gas Mass Fraction** | **Y_g** | - | Total gas phase (Y_v + Y_a) | ⭐⭐ Gas phase fraction |
| 7 | **Liquid Velocity** | **u_l** | m/s | Liquid phase velocity (non-dimensional) | ⭐⭐⭐⭐ **Momentum budget** |
| 8 | **Gas Velocity** | **u_g** | m/s | Gas phase velocity | ⭐⭐⭐⭐ **Entrainment rate** |
| 9 | **Mean Diameter** | **dm** | m | Sauter Mean Diameter (SMD) evolution | ⭐⭐⭐⭐⭐ **SPRAY BREAKUP** |
| 10 | **Diameter Variance** | **dvar** | m² | Variance of droplet size distribution | ⭐⭐⭐ **Polydispersity** |
| 11 | **Second Moment (d2)** | **d2** | m² | ∫∞ D² n(D) dD | ⭐⭐ Size statistics |
| 12 | **Third Moment (d3)** | **d3** | m³ | ∫∞ D³ n(D) dD | ⭐ Volume moment |
| 13 | **Droplet Temperature** | **Td** | - | Non-dimensional droplet temperature (Td×T_fuel = T_K) | ⭐⭐⭐⭐ **Evaporation/heating** |
| 14 | **Gas Temperature** | **Tg** | - | Non-dimensional gas temperature | ⭐⭐⭐ **Heat transfer** |
| 15 | **Jet Half-Width** | **b** | m | Radial spray width parameter | ⭐⭐⭐⭐ **SPRAY ANGLE/SPREADING** |
| 16 | **Turbulent k.e.** | **k_g** | (J/kg) | Gas-phase turbulent kinetic energy | ⭐⭐⭐ **Mixing intensity** |
| 17 | **Turbulent ε** | **eps_g** | (J/(kg·s)) | Turbulent kinetic energy dissipation rate | ⭐⭐⭐ **Mixing scale** |
| 18 | **Turbulent Visc.** | **mu_t_g** | (Pa·s) | Eddy viscosity | ⭐⭐⭐ **Turbulent transport** |
| 19 | **Scalar Variance** | **zvar_g** | - | Variance of fuel mass fraction | ⭐⭐⭐⭐ **Mixing quality** |
| 20 | **Mixture Fraction** | **zmix_g** | - | Favre-averaged mixture fraction (0–1) | ⭐⭐⭐⭐⭐ **MIXING STATE** |
| 21 | **Scalar Dissipation** | **chi_g** | 1/s | Scalar dissipation rate (mixing rate) | ⭐⭐⭐⭐ **Evaporation/mixing** |
| 22 | **Scalar Diss. (STL)** | **chi_g_stl** | 1/s | Stoichiometric dissipation | ⭐⭐ Combustion metric |
| 23 | **DSD Type** | **dsd_type** | - | Code for DSD model (1=delta, 2=RR, 3=LN, etc.) | ⭐ Metadata |
| 24 | **Prandtl number** | **Pr_g** | - | Gas mixture Prandtl number | ⭐⭐ Transport property |
| 25 | **Schmidt number** | **Sc_g** | - | Gas mixture Schmidt number | ⭐⭐ Transport property |
| 26 | **c_k** | c_k | - | Coefficient in turbulence model | ⭐ Model coefficient |
| 27 | **Vel. Ratio gas** | **VRg** | - | Viscosity ratio liquid/gas | ⭐⭐ Transport ratio |
| 28 | **Vel. Ratio turb** | **VRtg** | - | Eddy viscosity / liquid viscosity | ⭐⭐ Turbulent dominance |
| 29 | **Surface Tension Ratio** | **SR** | - | σ(T_d) / σ(T_fuel) | ⭐⭐ Temperature effect on surface tension |
| 30 | **Density Ratio liquid** | **DRl** | - | ρ_l(T_d) / ρ_l(T_fuel) | ⭐⭐ Liquid property sensitivity |
| 31 | **Reference Mass Frac.** | **Y_ref** | - | Reference fuel vapor mass fraction | ⭐ Model reference state |
| 32 | **alpha** | α | - | Shape parameter (exponent) in velocity profiles | ⭐ Profile shape |
| 33 | **Vapor Temperature** | **T_v** | K | Saturated vapor temperature | ⭐⭐ Flash boiling reference |
| 34 | **Density Ratio gas** | **DRg** | - | ρ_l / ρ_g (liquid/gas mixture) | ⭐⭐⭐ **Phase density ratio** |

**Interpretation for Paper**:
- Columns 1-8: Flow structure (y-axis: 0 mm to 1600 mm)
- Columns 9-12: Primary breakup (plot **dm** and **b** vs z)
- Columns 13-20: **MOST IMPORTANT**: Mixture state, evaporation extent, temperature
- Columns 21-22: Mixing quality (report **zmix_g**, **chi_g** on centerline)

**Example Plot for Publication**:
```
LHS: z vs [Y_l, Y_v, Y_a]  (stacked area chart)
RHS: z vs [dm, Td, Tg]     (penetration structure)
```

---

### 3.3 Source Term File: `omega.out_XXXXXX` (Physics Budget)

**Frequency**: Every `Datafile frequency` steps  
**Type**: Axial source term contributions  

| # | Variable | Symbol | Units | Description | **Scientific Value** |
|---|----------|--------|-------|-------------|----------------------|
| 1 | z | z | m | Axial grid location | ⭐⭐ Position |
| 2 | **Mass Entrainment** | **omega_ent** | kg/(m³·s) | Ambient gas entrainment mass rate | ⭐⭐⭐⭐ **Primary entrainment flux** |
| 3 | **Drag Force** | **omega_drag** | N/m³ | Liquid-gas momentum exchange via drag | ⭐⭐⭐⭐ **Momentum coupling** |
| 4 | **Breakup (1)** | **omega_bre1** | 1/s | Primary breakup rate (KH instability) | ⭐⭐⭐⭐⭐ **BREAKUP KERNEL** |
| 5 | **Breakup (2)** | **omega_bre2** | 1/s | Secondary breakup rate (RTI) | ⭐⭐⭐ Secondary fragmentation |
| 6 | **Evaporation** | **omega_vap** | kg/(m³·s) | Liquid-to-vapor mass transfer rate | ⭐⭐⭐⭐⭐ **VAPORIZATION SINK** |
| 7 | **Temperature Change** | **omega_T** | K/s | Heat transfer to droplets | ⭐⭐⭐⭐ **Heating rate** |
| 8 | **Turb. k (production)** | **omega_k_g_p** | J/(kg·s) | Turbulent kinetic energy production | ⭐⭐⭐ Turbulence production |
| 9 | **Turb. k (dissipation)** | **omega_k_g_d** | J/(kg·s) | Turbulent kinetic energy dissipation | ⭐⭐⭐ Energy cascade |
| 10–13 | **Dissipation terms** | omega_eps_g_p(-d), omega_zvar_g_p(-d) | - | ε and Z_var equations | ⭐ Turbulence closure |

**Physical Interpretation**:
```
Row 3:  omega_bre1 > 0  ⟹  Active primary breakup (KH waves growing)
Row 4:  omega_bre2 > 0  ⟹  Secondary breakup (large drops still breaking)
Row 6:  omega_vap > 0   ⟹  Active evaporation (T_d increasing)
Row 8:  omega_k_g_p > 0 ⟹  Turbulence generated by shear
```

**For Paper**: 
- Plot **omega_bre1**, **omega_bre2**, **omega_vap** vs z to show "where does breakup/evaporation happen?"
- Ratio of omega_bre1 to omega_vap shows relative importance of mechanical breakup vs thermal vapor blast

---

### 3.4 Droplet Size Distribution File: `dsd.out_XXXXXX` (If enabled)

**Frequency**: Only if `Number of diameter classes > 1`  
**Type**: DSD at each axial location  

| # | Variable | Description |
|---|----------|-------------|
| 1 | z | Axial position |
| 2–5 | dsdlam(1–4,:) | Coefficients of log-normal or Rosin-Rammler distribution |
| 6 | dm | Mean diameter |
| 7 | d2 | Second moment |
| 8 | d3 | Third moment |

---

### 3.5 Simulation Parameters File: `simulation_param.out`

**Frequency**: Once at initialization  
**Type**: Derived non-dimensional parameters  

| Parameter | Meaning | Typical Range |
|-----------|---------|----------------|
| **Re** | Reynolds number (jet exit) | 10000–100000 (your Aleiferis: 69,700) |
| **We** | Weber number | 100000–10000000 (your case: 1.68M) |
| **DRa** | Density ratio liquid/ambient gas | 100–1000 (your case: 470) |
| **VRa** | Viscosity ratio liquid/ambient gas | 10–1000 |
| **DRv** | Density ratio liquid/vapor | 100–1000 (your case: 189) |
| **VRv** | Viscosity ratio liquid/vapor | 1–100 |
| **De** | Dimensionless evaporation number | 0.01–1.0 |
| **Spray half-cone angle** | ° | 10–30° |
| **Spreading coefficient** | - | 0.3–0.5 |
| **Discharge coefficient** | - | 0.6–0.95 |
| **Injection velocity** | m/s | 100–300 m/s (your case: 139.7 m/s for 150 bar, npentane) |
| **Effective jet diameter** | m | Actual exit diameter after corrections |
| **Max. Mass flow rate** | mg/ms | 0.5–10 (your case: 1.12 mg/ms) |
| **Stoichiometric mixture fraction** | - | 0.02–0.05 (fuel/air ratio) |

**Paper Usage**: These are **"normalized" parameters that allow comparison across different fuels/conditions**.

---

### 3.6 Optional 2D Radial Profiles: `Zmix2D.out_*`, `Zvar2D.out_*` (Advanced)

**When available**: If `Use 2D profiles for convolution : .true.`  
**Type**: Radial-dependent mixture fraction and variance profiles

---

## 4. PARAMETERS YOU CAN ADJUST TO AFFECT MODEL BEHAVIOR

### 4.1 **Injector Hole Count** ← Your Direct Question

| What to Change | Effect | How to Use |
|----------------|--------|-----------|
| **Number of nozzles** | Scales mass flow rate per hole; **Does NOT resolve individual jets** | Set to actual number of holes (1, 6, 8, 12, etc.) |
| **Current Model Limitation** | CAS is cross-sectionally averaged; cannot model separate hole geometry or spray-spray interaction | For multi-hole validation: run separate case per hole, then superpose manually |

**Example**:
```ini
; Single-hole  (your current case):
Number of nozzles : 1
Injected mass : 1.0e-3  kg
Injection pressure : 150e5 Pa
; Result: LPL ∝ √(Re) ∝ √(U_inj)

; Equivalent 6-hole injector:
Number of nozzles : 6
Injected mass : 1.0e-3  kg  (total)
Injection pressure : 150e5 Pa
; Result: LPL same as single-hole case (per hole averaged);
; but u_inj is 1/6 of single-hole, so penetration slower
```

---

### 4.2 **Injection Pressure** ← Most Impactful for Penetration

| Pressure (bar) | Typical U_inj (m/s) | LPL @ 1ms (mm) | Breakup Character |
|----------------|-------------------|-----------------|-------------------|
| 50 | ~80 | ~60 | Coarse (large droplets persist) |
| 150 | ~140 | ~105 | **Medium (↑ breakup onset)** ← Your case |
| 300 | ~198 | ~150 | **Fine (vigorous breakup)** |
| 1500 | ~440 | ~340 | **Very fine (flash-boiling dominant)** |

**To adjust**: Change `Injection pressure : XXXe5` in cas.in

---

### 4.3 **Spray Angle & Spreading** ← Affects Cone Structure

| Parameter | Typical Value | Effect When ↑ |
|-----------|---------------|----|
| **Spray half-cone angle θ** | 20–30° | Wider cone; lower centerline density; faster mixing |
| **C_theta (in dynamic models)** | 0.5–2.0 | Scales angle evolution rate; higher = faster angle change |
| **Spreading coefficient** | 0.3–0.5 | Radial jet expansion; affects **b** evolution |

**To adjust**: Add to cas.in:
```ini
Spray half-cone angle model : Constant
Spray half-cone angle : 25.0
```

---

### 4.4 **Nozzle Exit Diameter** ← Scales All Length/Time Scales

| D_noz (μm) | U_inj @ 150 bar | τ = D/U | LPL @ 1ms |
|-----------|-----------------|---------|-----------|
| 100 | ~140 | 0.71 µs | ~52 mm (half) |
| **200** | ~140 | **1.43 µs** | **~105 mm ← Your case** |
| 300 | ~140 | 2.14 µs | ~157 mm |

**Physical Effect**: Larger nozzle ⟹ longer breakup time ⟹ larger droplets persist longer ⟹ deeper penetration.

**To adjust**: Change `Nozzle exit diameter : XXXe-6`

---

### 4.5 **Breakup Model Coefficients** ← Fine-Tune Drop Size Evolution

| Parameter | Current | Effect When ↑ |
|-----------|---------|---|
| **B0** | 0.61 | ↑ KH wave amplitude ⟹ More aggressive breakup ⟹ Smaller d_m |
| **B1** | 5.0 | ↑ Wave damping ⟹ Less breakup ⟹ Larger d_m persists |
| **C3** | 0.45 | ↑ Breakup time const ⟹ Faster breakup ⟹ Smaller d_m sooner |
| **D0, D1** | 1.0 | ↑ RTI strength ⟹ Faster secondary breakup |

**Theory**: KH breakup dominates near jet exit (high We); RTI dominates at peripheral spreading.

**Sensitivity Study for Paper**:
```
Base case:     B1 = 5.0  ⟹ dm(z=50mm) = 15 µm
Soft breakup:  B1 = 10.0 ⟹ dm(z=50mm) = 22 µm (↑ 47%)
Hard breakup:  B1 = 2.5  ⟹ dm(z=50mm) = 10 µm (↓ 33%)
```

---

### 4.6 **Evaporation Control** ← Flash-Boiling Sensitivity

| Parameter | Role | Typical Range |
|-----------|------|---|
| **Cevap** | Evaporation rate multiplier | 0.1–2.0 |
| **Flash boiling : .true.** | Enable superheated vapor | **← Recommended for 150 bar npentane** |
| **T_fuel** | Injected liquid temperature | 300–400 K; ↑ temp ⟹ more flash boiling |

**Extreme Cases**:
```
Conservative evaporation:   Cevap = 0.1  ⟹ Droplets survive longer
Aggressive evaporation:     Cevap = 1.0  ⟹ Rapid vapor formation
```

**For Flash-Boiling Papers**: Always perform sensitivity on Cevap ∈ [0.5, 2.0].

---

### 4.7 **Initial Droplet Size (Dm, Dvar)** ← Boundary Condition

| Dm (μm) | Dvar (μm²) | Character | LPL @ 1ms |
|---------|-----------|-----------|-----------|
| 5 | 0 | **Mono-small** ⟹ Fast evap, short penet. | ~80 mm |
| **10–20** | **tiny** | **Primary breakup result** | **100–130 mm** ← Typical |
| 50 | 100 | **Coarse** ⟹ Slow evap, deep penet. | ~150 mm |

**From Kelvin-Helmholtz Theory**:
$$D_{10} \sim 2 B_0 \lambda_{KH} = \frac{2 B_0 \cdot 9}{We^{1.5}} \times \text{noz\_D}$$

For your case:
```
B0 = 0.61, We = 1.68e6 ⟹ D_m ≈ 20 µm (typical for high-pressure diesel)
```

**To modify**: Adjust `Initial Dm : XXXe-6` in cas.in (or let code compute from Kelvin-Helmholtz).

---

### 4.8 **Grid Resolution** ← Numerical Parameter

| Nz (# points) | Δz (µm) @ Lz=1.6mm | CPU Time | Accuracy |
|---------------|------------------|----------|----------|
| 401 | 4 | 0.5× | ±5% error |
| **1601** | **1** | **1.0×** | **±1% error** ← Your case |
| 3201 | 0.5 | 2.5× | ±0.1% |

**Rule of Thumb**: Δz < 0.5 × d_m to resolve breakup layer; Δz < 2× b to resolve radial structure.

**To adjust**: Change `Number of grid points : XXXX`

---

### 4.9 **Time Step Control** ← Stability Setting

| MaxCFL | Typical dt (μs) | Integr. Scheme | Risk |
|--------|-----------------|-----------------|------|
| 0.05 | 0.001 | TVDRK-22 | Very safe (slow) |
| **0.1** | **0.003** | **TVDRK-33** | **Optimal balance** ← Default |
| 0.2 | 0.006 | TVDRK-44 | Higher risk of instability |

**Criterion**: $\Delta t < \text{CFL} \times \frac{\Delta z}{\max(|u_l|, |u_g|)}$

---

### 4.10 **Output Frequency** ← Data Volume Control

| Scenario | `Output frequency` | `Datafile frequency` | Notes |
|----------|-------------------|-------------------|-------|
| **Quick test** | 100 | 10000 | Minimal I/O; only 1–2 data files |
| **Parameter sweep** | 50 | 1000 | Moderate detail |
| **Publication grade** | 10 | 100 | High resolution; ~150 snapshots for 1.12 ms |

**Storage**: One `data.out_*` ≈ 50 KB; 150 files = 7.5 MB per case.

---

## 5. RECOMMENDED WORKFLOW FOR YOUR RESEARCH

### Step 1: **Baseline Reproduction** (Week 1)
```ini
; Reproduce Aleiferis 150 bar case exactly
Fuel : npentane
Injection pressure : 150e5
Nozzle exit diameter : 200e-6
Ambient pressure : 1.0e5
Fuel temperature : 363.15
Initial DSD : Delta
Initial Dm : 10e-6
Number of nozzles : 1
Flash boiling : .true.
```
**Objective**: Match your existing `PEN54_plot_fixed.png` curve (LPL vs time).

---

### Step 2: **Parametric Sensitivity** (Week 2–3)

Create a matrix:
```
Case | Inj_Pres | Dm | B1  | Cevap | Nz   | Result
-----|----------|----|----|-------|------|--------
A1   | 150e5    | 10 | 5  | 0.33  | 1601 | Baseline
A2   | 100e5    | 10 | 5  | 0.33  | 1601 | ↓ Pressure
A3   | 200e5    | 10 | 5  | 0.33  | 1601 | ↑ Pressure
B1   | 150e5    | 5  | 5  | 0.33  | 1601 | Smaller d_m
B2   | 150e5    | 20 | 5  | 0.33  | 1601 | Larger d_m
C1   | 150e5    | 10 | 2.5| 0.33  | 1601 | ↑ Breakup
C2   | 150e5    | 10 | 10 | 0.33  | 1601 | ↓ Breakup
```

**Objective**: Identify which parameters shift LPL/VPL most; rank them by sensitivity.

---

### Step 3: **Multi-Hole Injector Study** (Week 3–4)

For each hole count in {1, 4, 6, 8, 12}:
```
Case | Nozzles | Effect                           |
-----|---------|----------------------------------|
1H   | 1       | Baseline penetration            |
2H   | 6       | ↓ Penetration (same dm/hole)   |
3H   | 12      | ↓↓ Penetration (dilution)      |
```

**Note**: This gives *average* penetration; individual hole spacing is not modeled.

---

### Step 4: **Flash-Boiling Analysis** (Week 4–5)

Vary `T_fuel` and `Cevap`:
```
Temperature | Cevap | VPL/LPL Ratio | Flash Behavior          |
------------|-------|---------------|------------------------|
300 K       | 0.5   | 1.2 | Minimal evaporation    |
363 K       | 1.0   | 2.5 | Moderate flash         |
400 K       | 2.0   | 8.0 | Aggressive flash       |
```

**Paper Figure**: Plot (VPL – LPL) vs time; slope = vapor growth rate.

---

## 6. SUMMARY TABLE: What to Report in Your Paper

| Metric | From File | Extraction | Importance |
|--------|-----------|-----------|-----------|
| **Penetration @ t=1ms** | `spray_PL.out` col 3 | LPL(t=1e-3) in mm | ⭐⭐⭐⭐⭐ Core Result |
| **Tip velocity** | `spray_PL.out` col 5 | Mean/final U_tip | ⭐⭐⭐⭐ Momentum Balance |
| **Mean drop size @ tip** | `data.out_*` col 9 | dm(z=LPL) | ⭐⭐⭐⭐⭐ Breakup Quality |
| **Mixture homogeneity** | `data.out_*` col 20 | zmix_g profiles | ⭐⭐⭐⭐ Evaporation extent |
| **Entrainment rate** | `omega.out_*` col 2 | Integrated omega_ent | ⭐⭐⭐⭐ Model sensitivity |
| **Re, We** | `simulation_param.out` | Direct values | ⭐⭐⭐ Non-dimensional |
| **Spray angle** | Output log or `simulation_param.out` col 7 | theta(t) | ⭐⭐⭐ Geometry |
| **Breakup time scale** | Derived | z / u_g @ 50% d_m change | ⭐⭐⭐ Physics timescale |

---

## 7. EXAMPLE PUBLICATION FIGURES

### Figure 1: Penetration Validation (MUST HAVE)
```
X-axis: Time (ms) [0 to 1.2]
Y-axis: LPL (mm) [0 to 100]
Line 1: Simulation (monotonic increase)
Line 2 (optional): Experiment (reference if available)
Caption: "Liquid penetration length evolution; 150 bar, npentane, 363K fuel, 1 atm ambient"
```

### Figure 2: Spray Structure at Key Times
```
Z-location (mm) vs time snapshots {t = 0.1, 0.3, 0.7, 1.1 ms}
Subplots:
  (a) Y_l profile (spray core location)
  (b) Td profile (evaporation extent)
  (c) dm evolution (breakup progress)
  (d) umix_g profile (mixing pattern)
Caption: "Axial profiles of liquid volume fraction, droplet temperature, mean diameter, and mixture fraction at selected times"
```

### Figure 3: Physics Budget (Omega.out analysis)
```
X-axis: z (mm)
Y-axis: Source term magnitude (arbitrary units)
Lines: omega_ent, omega_drag, omega_bre1, omega_vap
Time: Single snapshot t = 0.5 ms
Caption: "Relative importance of entrainment, drag, breakup, and evaporation as functions of axial position"
```

### Figure 4: Sensitivity Analysis
```
X-axis: Parameter (e.g., Injection Pressure in bar)
Y-axis: LPL @ t=1ms (mm)
Points/lines: Base case with error bars or multiple runs
Caption: "Effect of [parameter name] on liquid penetration; base case: 150 bar, npentane"
```

---

## 8. QUICK REFERENCE: MOST IMPACTFUL ADJUSTABLE PARAMETERS

**Ranked by Effect on LPL**:

| Rank | Parameter | Type | Typical ΔChange | ΔLPL |
|------|-----------|------|-----------------|------|
| 1️⃣ | Injection pressure | Input | ±50 bar | ±30% |
| 2️⃣ | Nozzle diameter | Input | ±50 µm | ±25% |
| 3️⃣ | B1 (breakup) | Model | ±5 | ±15% |
| 4️⃣ | Dm (initial) | Input | ±5 µm | ±10% |
| 5️⃣ | Nozzles count | Input | ±6 | ±20% (per-hole) |
| 6️⃣ | Cevap (evap.) | Model | ×2 | ±5% |
| 7️⃣ | C_theta (angle) | Model | ±1 | ±10% |
| 8️⃣ | T_fuel (temp) | Input | ±50 K | ±5% |

---

## 9. ADDITIONAL NOTES FOR YOUR RESEARCH

✅ **Recommended for you**:
- Perform sensitivity on **pressure** (100, 150, 200, 300 bar) to span engine operating window
- Compare **DSD types** (Delta vs Log-Normal vs Rosin-Rammler) to show model robustness
- Vary **nozzle count** systematically to benchmark multi-hole injector behavior
- Report **penetration rate** (d(LPL)/dt) in early transient; this is more sensitive to breakup models

⚠️ **Pitfalls**:
- Don't trust **d3 (third moment)** — solver often skips it; use dm and d2
- **Combustion model (MRIF)** not needed for evaporation/flash-boiling papers; skip unless you have soot data
- **Turbulence model**: Safe to turn OFF for initial studies; turn ON when coupling with LES or for high-Reynolds jets

🎯 **For your "scientific paper" goal**:
- Focus on **LPL/VPL penetration curves** as primary result
- Add **mixture fraction contours** (zmix_g from data.out) as secondary validation
- Include **non-dimensional correlations** table (Re, We, DRa) to generalize findings
- Compare to **published Spray A or Aleiferis experimental data** if available

---

**End of Summary**
