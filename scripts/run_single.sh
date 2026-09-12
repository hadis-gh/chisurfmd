#!/usr/bin/env bash

set -eu

# ==============================================================================
SCRIPT_DIR="$(
    cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &&
    pwd
)"

PROJECT_ROOT="$(
    cd -- "$SCRIPT_DIR/.." &&
    pwd
)"
# ==============================================================================
# 1. LOAD CONFIGURATION
# ==============================================================================
config=${1:-"config.sh"}

if [ ! -f "$config" ]; then
    echo "Error: configuration file '$config' not found."
    exit 1
fi

source "$config"

# ==============================================================================
# 2. PATHS AND OUTPUT
# ==============================================================================
MD_EXE=${MD_EXE:-"$PROJECT_ROOT/build/examples/chisurfmd_md"}

outputDir=${outputDir:-"$PROJECT_ROOT/outputs"}
saveFile=${saveFile:-"${outputDir}/final_run.bp"}

mkdir -p "$outputDir"

logFile="${outputDir}/log.txt"
exec > >(tee "$logFile") 2>&1

# ==============================================================================
# 3. DEFAULT RUNTIME OPTIONS
# ==============================================================================
potentialType=${potentialType:-"OrientedLJ"}

particleInit=${particleInit:-"RANDOM"}
temperature=${temperature:-0.5}
integration=${integration:-"VelocityVerlet"}

# ==============================================================================
# 4. COMMON MD ARGUMENTS
# ==============================================================================
args=(
    --particlesInit "$particleInit"
    --particlesNum "$particlesNum"
    --particleRadius "$particleRadius"
    --momentI "$momentI"
    --chirality "$chirality"
    --alignment "$alignment"
    
    --areaL "$areaL"
    --seed "$seed"

    --integration "$integration"
    --time "$time"
    --dt "$dt"
    --writeStateInterval "$writeStateInterval"
    --writeEnergyInterval "$writeEnergyInterval"
    
    --temperature "$temperature"
    --thermoInterval "$thermoInterval"
    --collisionFr "$collisionFr"

    --saveFile "$saveFile"
)

# ==============================================================================
# 5. POTENTIAL-SPECIFIC ARGUMENTS
# ==============================================================================

if [ "$potentialType" == "OrientedLJ" ]; then
    args+=(
        --LJangularScale "$LJangularScale"
        --LJPhiOrder "$LJPhiOrder"
        --LJalpha "$LJalpha"
    )
fi

# ==============================================================================
# 6. RUN SIMULATION
# ==============================================================================

echo "=================================================="
echo "        Single Molecular Dynamics Simulation      "
echo "=================================================="
echo "Config:     $config"
echo "Potential:  $potentialType"
echo "Output:     $saveFile"
echo "--------------------------------------------------"

"$MD_EXE" "${args[@]}" || {
    echo "Error: simulation failed."
    exit 1
}

echo
echo "Simulation completed successfully."
echo "Output written to: $saveFile"
echo "=================================================="