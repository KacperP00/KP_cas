# COMPLETE GUIDE: COMPILING AND RUNNING "CROSS-SECTIONALLY-AVERAGED-SPRAY-MODEL" IN WSL (UBUNTU)

## Initial Notes
* The model's name is quite long, so for convenience, it is highly recommended to rename the main folder to "cas".
* When working with the "Aleiferis injector" example, you should also rename its folder to something shorter and without spaces, e.g., `Al_inj`.
* The only 100% necessary modification is **Step 4**; the rest depends on your hardware.
* You can try compiling first (starting with Step 4 and 5) to see if Steps 1-3 are even necessary. Apply the first three steps only if compilation errors appear.

---

## Introduction and Problem Diagnosis
The default code for this model (written in Fortran) was designed specifically for the Intel compiler. Attempting to compile it straight away using the default Linux compiler (`gfortran`) will fail. 

During the compilation process, 3 main obstacles were encountered:
1. Lack of the Intel Fortran environment in the WSL system.
2. An outdated `Makefile` referenced the `ifort` compiler, which has been deprecated in recent Intel releases and replaced by `ifx` (resulting in Error 127).
3. A bug in the author's source code (Error #6460)—the calibration variables `A0` and `A1` were not declared in the main definition file and the parser, meaning the model did not use the declared coefficients. 

Below are the instructions to achieve a working simulation.

---

## Step 0: Clone the Repository
```bash
git clone [https://git.rwth-aachen.de/avijitsaha021/cross-sectionally-averaged-spray-model.git](https://git.rwth-aachen.de/avijitsaha021/cross-sectionally-averaged-spray-model.git)
```

---

## Step 1: Install the Intel Fortran Compiler (OneAPI)
If you are using WSL with Ubuntu, you must download the official Intel packages. Open the terminal (you can be in any directory, e.g., your home directory `~`) and execute the following system commands consecutively:

**a) Download the Intel repository security key:**
```bash
wget -O- [https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB](https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB) | gpg --dearmor | sudo tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null
```

**b) Add the Intel repository to your system's source list:**
```bash
echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] [https://apt.repos.intel.com/oneapi](https://apt.repos.intel.com/oneapi) all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list
```

**c) Update the package list and install the compiler (confirm with "y" if prompted):**
```bash
sudo apt update
sudo apt install intel-oneapi-compiler-fortran
```

---

## Step 2: Activate the Intel Environment
Simply installing the compiler is not enough. Your terminal needs to "know" where it is located. Enter the following in your terminal:
```bash
source /opt/intel/oneapi/setvars.sh
```

> **IMPORTANT NOTE:** You must run this command every single time you open a NEW terminal window and want to work with this code. 

**To avoid doing this manually:**
1. Open the `~/.bashrc` file: 
   `nano ~/.bashrc`
2. Scroll to the very bottom and add the following line: 
   `source /opt/intel/oneapi/setvars.sh`
3. Save the file and execute this in the terminal: 
   `source ~/.bashrc`

The code will now execute this command automatically upon opening.

---

## Step 3: Fix the "Makefile" (ERROR: ifort: No such file or directory)
The authors hardcoded the outdated `ifort` compiler. The newest Intel package uses the `ifx` compiler.

1. Using an editor, navigate to the `src` folder and open the file named `Makefile`.
2. Press `Ctrl + H` (Find and Replace).
3. In the "Find" field, type: `ifort`.
4. In the "Replace" field, type: `ifx`.
5. Click "Replace All" and save the file (`Ctrl + S`).

---

## Step 4: Fix the Source Code (ERROR #6460: MISSING VARIABLES A0 and A1)
During compilation, the `spray_func.f90` file will report a problem with the missing variables `spray%A0` and `spray%A1`. We must add them to the definitions file.

1. In the editor, navigate to the `src` folder and open the `spray_defs.f90` file.
2. Press `Ctrl + F` and search for the phrase: `Cevap = 1.0_WP`.
3. You will find a code snippet that looks like this:
   ```fortran
     ! Evaporation parameters
     real(WP) :: Cevap = 1.0_WP
   ```
4. Add new lines directly below it to include the missing variables. The whole section should look like this:
   ```fortran
     ! Evaporation parameters
     real(WP) :: Cevap = 1.0_WP
     real(WP) :: A0 
     real(WP) :: A1 
   ```
5. Save the file (`Ctrl + S`).
6. We also need to add the corresponding sections to the parser. Open the `src/parser.f90` file.
7. Search for the following snippet:
   ```fortran
          else if(param .eq. 'Use evaporation model') then
             read(line(idx+1:len(line)),*) spray%evap_model
          else if(param .eq. 'B0') then
             read(line(idx+1:len(line)),*) spray%B0
   ```
8. Add two conditions for `A0` and `A1` so the code looks like this:
   ```fortran
          else if(param .eq. 'Use evaporation model') then
             read(line(idx+1:len(line)),*) spray%evap_model
          else if(param .eq. 'A0') then
             read(line(idx+1:len(line)),*) spray%A0
          else if(param .eq. 'A1') then
             read(line(idx+1:len(line)),*) spray%A1
          else if(param .eq. 'B0') then
             read(line(idx+1:len(line)),*) spray%B0
   ```

---

## Step 5: Compile the Code
With the compiler installed, `ifort` replaced by `ifx`, and the author's bug patched, we can build the program.

1. In the terminal, navigate to the sources folder (depending on where you have it):
   ```bash
   cd ~/cross-sectionally-averaged-spray-model/src
   ```
2. Run the compilation process:
   ```bash
   make all_opt
   ```

If everything went well, the terminal will no longer throw any "Error", and a green executable file named `cas_opt` will be generated in the `bin` folder.

---

## Step 6: Running the First Simulation (Example)
To ensure everything works, run the built-in example provided by the authors - the Aleiferis injector.

1. From the terminal in the `src` folder, navigate to the example directory:
   ```bash
   cd ../example/"Aleiferis injector" 
   ```
   *(Or `cd ../example/Al_inj` if you renamed it as suggested at the beginning of the guide).*
2. Copy the newly created executable `cas_opt` to the current directory:
   ```bash
   cp ../../bin/cas_opt .
   ```
   > **Note:** Do not forget the dot (`.`) at the end of the command!
3. Run the simulation by passing its input file with parameters (`cas.in`):
   ```bash
   ./cas_opt cas.in
   ```

Logs and results from the model's first injector simulation should now start appearing in the terminal.

---

## Troubleshooting: Additional Source Code Changes
The following are supplementary modifications introduced during attempts to repair the code. Some of these are neutral, but if the basic guide does not yield correct results, it is probable that you should apply change (1), followed by change (3). 

**Warning:** Applying changes (2) and (4) is not recommended—they might be/are user errors.

### Change 1: `spray_func.f90` (Condition Modification)
Change the original condition from:
```fortran
        if(tauTHM(1)<tauRT(1).and.tauTHM(1)<tauKH(1))then changed by me KP
          K_bre1(k) = K_bref1(k)   
          K_bre2(k) = K_bref2(k)     
          K_bre3(k) = K_bref3(k)   
        end if
```

To the following:
```fortran
        if (spray%flash_boiling .and. tauTHM(1) < 9000.0_WP) then
          !print *, "FLASH-BOILING ACTIVE IN CELL: ", k <- no need for that anymore
          K_bre1(k) = K_bref1(k)   
          K_bre2(k) = K_bref2(k)     
          K_bre3(k) = K_bref3(k)   
        end if
```

### Change 2: `spray_func.f90` (NOT RECOMMENDED)
Changes made in the author's remark block regarding the reference vapor mass fraction. 
```fortran
    !spray%Y_ref = (-1.5_WP/spray%b)*spray%Y_v
    !spray%Y_ref = (0.5_WP/spray%b)*(spray%Y_v + 2.0_WP*Ystar_fe)/3.0_WP
    !spray%Y_ref = (0.5_WP/spray%b)*(spray%Y_v + 2.0_WP*Ystar_fe)/3.0_WP <-this one was imported 
    spray%Y_ref = (spray%Y_v + 2.0_WP*Ystar_fe)/3.0_WP !added by me KP
    !spray%Y_ref = (spray%Y_v + 2.0_WP*Ystar_fe)/3.0_WP
    !spray%Y_ref = spray%Y_v
```

### Change 3: `spray_func.f90` (Adding Variables Accumulation)
Addition of accumulation tracking for vaporization parameters:
```fortran
          spray%omega_exp(k)   = max(0.0_WP,(1.5_WP*rho(k)*Y_l(k))*sum((dsd(:,k)*K_vap_swl/di(:,k)**2)))
          spray%omega_expdm(k) = max(0.0_WP,(1.5_WP*rho(k)*Y_l(k))*sum((dsd(:,k)*K_vap_swl/di(:,k))))
          spray%omega_expd2(k) = max(0.0_WP,(1.5_WP*rho(k)*Y_l(k))*sum((dsd(:,k)*K_vap_swl)))
          spray%omega_expd3(k) = max(0.0_WP,(1.5_WP*rho(k)*Y_l(k))*sum((dsd(:,k)*K_vap_swl*di(:,k))))
         ! added by me KP
          spray%omega_vap(k)   = spray%omega_vap(k)   + spray%omega_exp(k)
          spray%omega_vapdm(k) = spray%omega_vapdm(k) + spray%omega_expdm(k)
          spray%omega_vapd2(k) = spray%omega_vapd2(k) + spray%omega_expd2(k)
          spray%omega_vapd3(k) = spray%omega_vapd3(k) + spray%omega_expd3(k)
```

### Change 4: `spray_func.f90` (NOT RECOMMENDED)
Another test change regarding reference variable formulation:
```fortran
          !Yref = spray%Y_ref(k) !spray%Y_v(k)
          !Yref = (spray%Y_v(k)/spray%Y_g(k) + 2.0_WP*Ystar_f)/3.0_WP
          !Yref = (0.5_WP/spray%b(k))*(spray%Y_v(k) + 2.0_WP*Ystar_f)/3.0_WP <- originaly imported
          Yref = (spray%Y_v(k) + 2.0_WP*Ystar_f)/3.0_WP !added by me KP
          !Yref = (spray%Y_v(k) + 2.0_WP*Ystar_f)/3.0_WP
```