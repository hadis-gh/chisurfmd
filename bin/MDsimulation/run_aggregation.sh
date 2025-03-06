#!/usr/bin/env bash

set -eu

start_time=$(date +%s)

config=${1:-"config.sh"}
if [ ! -f "$config" ]; then
    echo "Error: '$config' not found."
    exit 1
fi
source "$config"

MD_EXE=${MD_EXE:-"../test/testNParticleMD"}
DEPOSIT_MD_EXE=${DEPOSIT_MD_EXE:-"../test/testNParticleDepositionMD"}

outputDir=${outputDir:-"./outputs"}
mkdir -p "$outputDir"

logFile="log.txt"
exec > >(tee -a "$logFile") 2>&1
echo "Logging to: $logFile"

args=(
    --seed "$seed"
    --exclusionRadius "$exclusionRadius"
    --thermoInterval "$thermoInterval"
    --writeStateInterval "$writeStateInterval"
    --writeEnergyInterval "$writeEnergyInterval"
    --integration "$integration"
    --collisionFr "$collisionFr"
)
if [ "$simulationType" == "oriented" ]; then
    args+=(
        --LJangularScale "$LJangularScale"
        --LJPhiOrder "$LJPhiOrder"
        --LJalpha "$LJalpha"
        --momentI "$momentI"
    )
fi

echo "=================================================="
echo "    Aggregation Molecular Dynamics Simulation     "
echo "=================================================="

run_md() {
    local temperature=$1
    local equilibrationTime=$2

    printf "\nPart 1: Thermal Equilibration starting from FROZEN state: \n"
    printf "\nTemperature: %s\n" "$temperature"
    echo "------------------------------------------------------------"

    "$MD_EXE" "${args[@]}" \
        --particlesInit "initCluster/run_0.bp" \
        --temperature "$temperature" \
        --areaL 20 \
        --time "$equilibrationTime" --dt "$dt" \
        --saveFile "${outputDir}/equilibrated_${temperature}.bp" || {
        printf "Error: Thermal Equilibration simulation failed.\n"
        exit 1
    }

    printf "\nThermal Equilibration completed successfully.\n"
    echo "------------------------------------------------------------"
}

run_deposition() {
    local temperature=$1
    local depositionTime=$2

    printf "\nPart 2: Particle Deposition starting from RELAXED state: \n"
    printf "\nTemperature: %s\n" "$temperature"
    printf "Deposition Time Interval: %s\n" "$depositionTime"
    echo "------------------------------------------------------------"

    "$DEPOSIT_MD_EXE" "${args[@]}" \
        --particlesInit "${outputDir}/equilibrated_${temperature}.bp" \
        --particlesNumMax 100 \
        --areaL 50 \
        --temperature "$temperature" \
        --time "$depositionTime" --dt "$dt" \
        --saveFile "${outputDir}/aggregated_${temperature}_${depositionTime}.bp" || {
        printf "Error: Particle Deposition simulation failed.\n"
        exit 1
    }

    printf "\nParticle Deposition completed successfully.\n"
    echo "___________________________________________________________"
}

# Main Loop
temp_ranges=(0.5 0.25 0.3 0.4 0.35)
deposit_rates=(1 4 8)
equilibrationTime=10000

for T in "${temp_ranges[@]}"; do
    for rate in "${deposit_rates[@]}"; do
        run_md "$T" "$equilibrationTime"
        run_deposition "$T" "$rate"
    done
done

# Record time
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))

hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))

printf "Total execution time for Aggregation Process: %02d:%02d:%02d (hh:mm:ss)\n" "$hours" "$minutes" "$seconds"