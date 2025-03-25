#!/usr/bin/env bash

set -eu

start_time=$(date +%s)

MD_EXE=${MD_EXE:-"../build/test/testNParticleMD"}

outputDir=${outputDir:-"./outputs"}
mkdir -p "$outputDir"

logFile="benchmark_log.txt"
exec > >(tee -a "$logFile") 2>&1
echo "Logging to: $logFile"

args=(
    --seed 12345
    --areaL 100
    --exclusionRadius 0.8
    --thermoInterval 0.1
    --writeStateInterval 0.05
    --writeEnergyInterval 0.005
    --integration "VelocityVerlet"
    --collisionFr 0.01667
    --time 100
    --dt 0.01
    --temperature 0.3
)

echo "=================================================="
echo "   Running MD Benchmark for Various Particle Numbers     "
echo "=================================================="

run_md() {
    local particleNumber=$1
    local outfile="${outputDir}/run_${particleNumber}.bp"

    printf "\nRunning simulation with %s particles...\n" "$particleNumber"
    echo "------------------------------------------------------------"

    if [ -f "$outfile" ]; then
        echo "Skipping $particleNumber particles (already exists)."
        return
    fi

    # Start timing
    local start=$(date +%s)

    # Run simulation
    "$MD_EXE" "${args[@]}" --particlesNum "$particleNumber" --saveFile "$outfile" || {
        echo "Error: Simulation failed for N=$particleNumber."
        exit 1
    }

    # Stop timing
    local end=$(date +%s)
    local elapsed=$((end - start))
    printf "Simulation for %s particles completed in %d seconds.\n" "$particleNumber" "$elapsed"
    echo "------------------------------------------------------------"
}

particle_num_range=(20 40 60 80 100 150 200 300 400 500)

for num in "${particle_num_range[@]}"; do
    run_md "$num"
done

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))

hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))

printf "\nTotal execution time for benchmarking: %02d:%02d:%02d (hh:mm:ss)\n" "$hours" "$minutes" "$seconds"
