_Last updated: 2025-08-25_

# Project Directory
```plaintext
<project-root>/
├── README.md
├── Data_Figures_Paper/ # Data and figures for paper
├── include/
│ ├── lettuce/
│ ├── core/ # Core library code
│ ├── md/ # Molecular dynamics related code
│ ├── potential/ # Potential functions
│ └── constants.h # Global constants
├── script/
│ ├── bash/ # Bash scripts
│ └── notebooks/ # Jupyter or related notebooks
├── test/
│ ├── archive/ # Archived tests
│ ├── mc/ # Monte Carlo tests
│ ├── md/ # Molecular dynamics tests
│ └── unit/ # Unit tests
├── thirdParty/ # External dependencies
├── CMakeLists.txt # Main build configuration
├── .gitignore
├── .gitmodules
└── profile 
```

# Compiling

Lib Dependencies:
- Boost >= 1.74 (currently only boost::program_options is used)
- ADIOS2 >= v2.10.0
- Catch >= 3.9.1 (only used for testcases, which are under ~/test/unit)
lettuce_md_PARTICLE

### Normal Build
```
mkdir build && cd build
cmake -S ../. -B . -DCMAKE_CXX_FLAGS="-O2" -Dlettuce_md_PARTICLE="ParticleDot"
cmake --build . -t testMD
./test/testMD -n <NumParticles>
```

### CTest
```
mkdir build_testing && cd build_testing
cmake -S ../. -B . -DBUILD_TESTING=ON
cmake --build .
ctest --verbose --output-on-failure
```

# Benchmarks

Under Hardware/Compiler:
- AMD Ryzen 9 7950X3D 16-Core Processor
- AMD Radeon RX 7900 XT
- g++ (Ubuntu 12.3.0-1ubuntu1~22.04.2) 12.3.0
- HIP version: 6.0.32831-204d35d16
- AMD clang version 17.0.0
- Linux 6.8.0-65-generic #68~22.04.1-Ubuntu 

App-Target Runtimes

_(the default values are used, if the parameters are not specifed.)_


| Date  | Target    | ParticleType  | NumParticles  | ForceType | Thermostat    | Total_time    | ExecutionTime | Notes |
| ----  | ----      | ----          | ----          | ----      | ----          | ----          | ----          | ----  | 
| ----  | ----      | ----          | ----          | ----      | ----          | ----          | ----          | ----  | 
| 25.08 | testMD | ParticleDot | 49 | IsotropicLJ | ThermostatID::None | 10 | 1.664s | g++ ohne -O2  |
| 25.08 | testMD | ParticleDot | 49 | IsotropicLJ | ThermostatID::None | 10 | 0.11s | g++ mit -O2  |
| ----  | ----      | ----          | ----          | ----      | ----          | ----          | ----          | ----  | 
| 25.08 | testMD | ParticleDot | 256 | IsotropicLJ | ThermostatID::None | 10 | 43.008s | g++ ohne -O2  |
| 25.08 | testMD | ParticleDot | 256 | IsotropicLJ | ThermostatID::None | 10 | 2.484s | g++ mit -O2  |
| ----  | ----      | ----          | ----          | ----      | ----          | ----          | ----          | ----  | 

_(For the reason, O2 is better than O3!??)_

# Useful Parameters and Marcos

| Parameter           | Data Type    | Default             | Required     | Description                                   |
| ------------------- | ------------ | ------------------- | ------------ | --------------------------------------------- |
| areaL               | real         | 20.0                | Yes          | simulation size                               |
| collisionFr         | real         | 0.0166667 (1/60)    | Yes          | collision frequency for Andersen thermostat   |
| dt                  | real         | 0.01                | Yes          | integration step size                         |
| enableCapVelocity   | boolean         | false               | Yes          | Enable capping of velocities                  |
| exclusionRadius     | real         | 0.8                 | Yes          | exclusion radius                              |
| help, h             | boolean         | false               | Yes          | print help                                    |
| integration         | string       | "VelocityVerlet"    | Yes          | integration method (VelocityVerletStep, EulerStep, EulerSymplecticStep)  |
| mass                | real | 1         | Yes          | mass of particles |
| maxVelocity         | array of real | \[1e5, 1e3]         | Yes          | capping amount for velocity {x-y, omega}      |
| momentI             | real         | 1.0                 | Yes          | moment of inertia                             |
| neighborDistances  | array of real | \[1.2, 1.5, 2.0]    | Yes          | Distances for counting neighbors {x-y, omega} |
| particlesDensity    | real         | (none)              | No          | packing density of particles                  |
| particlesInit       | string       | "RANDOM"            | Yes          | particle initialization(RANDOM, RANDOM2, DLA, TWO) |
| particlesNum, n     | unsigned int | 49                  | Yes          | number of initial particles                   |
| particlesType       | string       | "RRUU"              | Yes          | particles type handedness & orientation (RLUU, RRUD RRUU)       |
| printOptions        | boolean         | true                | Yes          | Print all runtime options                     |
| relaxationTime      | real         | 40.0                | Yes          | relaxation time for Berendsen thermostat      |
| saveFile            | string       | "outputs/run\_0.bp" | Yes          | file path to save simulation output           |
| seed                | unsigned int | (none)              | No          | random seed                                   |
| temperature, T      | real         | 0.3                 | Yes          | temperature                                   |
| thermoInterval      | real         | 0.1                 | Yes          | interval after which to apply thermostat      |
| time, t             | real         | 10.0                | Yes          | max simulation time                           |
| writeEnergyInterval | real         | 0.005               | Yes          | measurement Energy interval                   |
| writeStateInterval  | real         | 0.05                | Yes          | measurement State interval                    |
| LJepsilon  | real         | 1               | Yes          | epsilon in Lennard-Jones force and potential                    |
| LJsigma  | real         | 1                | Yes          | sigma in Lennard-Jones force and potential                    |
| LJcutoff  | real         | 10                | Yes          | cutoff distance for Lennard-Jones interactions                   |


| Macro Name            | Data Type | Default            | Required | Description                                                                       |
| --------------------- | --------- | ------------------ | -------- | --------------------------------------------------------------------------------- |
| `DLETTUCE_THERMOSTAT` | enum      | None               | Yes      | Thermostat type: `None`, `VelocityScaling`, `Berendsen`, `Andersen` |
| `DLETTUCE_PARTICLE`   | enum      | `ParticleOriented` | Yes      | Particle model: `ParticleDot`, `ParticleOriented`                                 |
| `DLETTUCE_POTENTIAL`   | enum      | `IsotropicLJ` | Yes      | Force and Potential Model: `ChiralLJGeometric`, `IsotropicLJ`, `OrientedLJ`, `FieldCoupled`, `TabularDFT`           |


