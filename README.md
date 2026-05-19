##################################################################################
# Code compilation and running your first simulation                             #
##################################################################################

## NOTE! The following readme file is updated version of the original one.
Initially, the definition of A0 and A1 coefficients (for breakup model) was missing, so it had to be added
to the spray_func.f90 file and parser.f90. Other changes:
- change of ifort into ifx in Makefile 
- some Intel-related corrections (machine dependent)
- if during compilation one have additional problems, e.g. Lack of the Intel Fortran environment in the WSL system
See section after compilation.



 Originally: {Important note before you compile: 
Default Makefile is written for Intel fortran compiler. A Makefile for GNU 
fortran compiler is in ./Makefiles/ but currently does not compile due to some 
unsupported features. }

First, clone the repository:
```git clone https://github.com/KacperP00/KP_cas.git```

To compile:

1. Enter into the 'src' directory
```cd src```
2. Run ```make all_debug``` to compile with debug options
OR
2. Run ```make all_opt``` to compile with optimization options

The file is created in the /bin/ directory.

To run: (if all_opt was chosen.)

1. Enter into the /example/Al_inj directory
 ```cd example/Al_inj```
2. Copy the executable or create a shortcut link 
```cp ../../bin/cas_opt .``` 
3. Run using the input file 'cas.in'
```./cas_opt cas.in```
4. Use 'python3 result.py' and compare results of Liquid penetration with PEN54 case from author's article (pen54_article.csv or check the publication).
5. If something is not correct and one does some changes to the code, the autom_recompile_testing.sh script was designed to
delete old output files etc., recompile the code and run the simulation again along with python script for plotting. 
```./autom_recompile_testing.sh```

## Introduction and Problem Diagnosis
The default code for this model (written in Fortran) was designed specifically for the Intel compiler. Attempting to compile it straight away using the default Linux compiler (`gfortran`) will fail. 

During the compilation process, 3 main obstacles were encountered:
1. Lack of the Intel Fortran environment in the WSL system.
2. An outdated `Makefile` referenced the `ifort` compiler, which has been deprecated in recent Intel releases and replaced by `ifx` (resulting in Error 127).
3. A bug in the author's source code (Error #6460)—the calibration variables `A0` and `A1` were not declared in the main definition file and the parser, meaning the model did not use the declared coefficients. 

Below are the instructions to achieve a working simulation.

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

The code will now execute this command automatically upon opening, BUT it can influence performance of the machine (the change is done to the general file)
- One can add a conditional statement that checks if 'cd' command leads to particular directory e.g. cas, but it can in a small way influence the 'cd' command due to additional logic.
- The cleanest way is to create the environment with 'direnv' or 'conda' (or miniconda)
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





##################################################################################
# Citations                                                                      #
##################################################################################

If you use this code for your research, please cite the following publications:

1. A. Y. Deshmukh, T. Grenga, M. Davidovic, L. Schumacher, J. Palmer, 
   M.A. Reddemann, R. Kneer and H. Pitsch. A reduced-order model for multiphase 
   simulation of transient inert sprays in the context of compression ignition 
   engines. Int. J. Multiph. Flow, vol. 147, page 103872, 2022. 
   doi: http://dx.doi.org/10.1016/j.ijmultiphaseflow.2021.103872

2. A. Y. Deshmukh, M. Davidovic, T. Grenga, R. Lakshmanan, L. Cai and H. Pitsch.
   A reduced-order model for turbulent reactive sprays in compression ignition 
   engines. Combust. Flame, vol. 236, page 111751, 2022. 
   doi: http://dx.doi.org/10.1016/j.combustflame.2021.111751

3. A. Saha, A. Y. Deshmukh, T. Grenga, and H. Pitsch. Physics-based reduced-order 
   modeling of flash-boiling sprays in the context of internal combustion engines. 
   doi: https://doi.org/10.48550/arXiv.2307.03722

##################################################################################
# Funding                                                                        #
##################################################################################

This work was funded by the Deutsche Forschungsgemeinschaft (DFG, German Research
Foundation) under Germany’s Excellence Strategy—Exzellenzcluster 2186 
“The Fuel Science Center” ID: 390919832 and also supported by the 
European Research Council (ERC) under the European Union’s Horizon 2020 research 
and innovation program (grant agreement no.695747).

##################################################################################

