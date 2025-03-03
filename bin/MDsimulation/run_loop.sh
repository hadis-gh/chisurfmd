#!/usr/bin/env bash

# Exit on error and unset variable
set -eu

start_time=$(date +%s)

config=${1:-"config.sh"}
if [ ! -f "$config" ]; then
    echo "Error: '$config' not found."
    exit 1
fi
source "$config"

numCoolingRuns=$(echo "scale=0; ($highTemperature - $lowTemperature) / $stepTemperature" | bc)
numHeatingRuns=$(echo "scale=0; ($highTemperature - $lowTemperature) / $stepTemperature" | bc)
totalRuns=$((numCoolingRuns + numHeatingRuns))

MD_EXE=${MD_EXE:-"../test/testNParticleMD"}
outputDir=${outputDir:-"./outputs"}
mkdir -p "$outputDir"

logFile="log.txt"
exec > >(tee "$logFile") 2>&1
echo "Logging to: $logFile"


simulationType=${simulationType:-"oriented"}

run_simulation() {
    local runIndex="$1"
    local particleInit="$2"
    local temperature="$3"
    local saveFile="${outputDir}/run_${runIndex}.bp"

    printf "Starting simulation #%s | Temperature: %s\n" "$runIndex" "$temperature"

    local args=(
        --particlesInit "$particleInit"
        --temperature "$temperature"
        --time "$timeCooling" --dt "$dt"
        --saveFile "$saveFile"
        --particlesNum "$particleNum"
        --seed "$seed"
        --areaL "$areaL"
        --exclusionRadius "$exclusionRadius"
        --thermoInterval "$thermoInterval"
        --writeStateInterval "$writeStateInterval"
        --writeEnergyInterval "$writeEnergyInterval"
        --integration "$integration"
        --collisionFr "$collisionFr"
    )

    if [ "$simulationType" == "oriented" ]; then
        args+=(--LJangularScale "$LJangularScale"
               --LJPhiOrder "$LJPhiOrder"
               --LJalpha "$LJalpha"
               --momentI "$momentI")
    fi

    "$MD_EXE" "${args[@]}" || { 
        printf "Error: Simulation #%s failed.\n" "$runIndex"; 
        exit 1; 
    }

    printf "\nSimulation %s of %s completed successfully.\n\n" "$runIndex" "$totalRuns"
    echo "------------------------------------------------------------"
}

run_temperature_loop() {
    local startTemp="$1"
    local stepTemp="$2"
    local numRuns="$3"
    local startIndex="$4"

    for ((i = 1; i <= numRuns; i++)); do
        local runIndex=$((startIndex + i))
        local temperature=$(echo "$startTemp + $i * $stepTemp" | bc)
        local prevOutput="${outputDir}/run_$((runIndex - 1)).bp"

        run_simulation "$runIndex" "$prevOutput" "$temperature"
    done
}


echo "======================================================="
echo "    Temperature Loop Molecular Dynamics Simulation     "
echo "======================================================="

initialRunIndex=0
run_simulation "$initialRunIndex" "RANDOM" "$highTemperature"

run_temperature_loop "$highTemperature" "$stepTemperature" "$numCoolingRuns" 0
run_temperature_loop "$lowTemperature" "$stepTemperature" "$numHeatingRuns" "$numCoolingRuns"

# Record time
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))

hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))

printf "Total execution time: %02d:%02d:%02d (hh:mm:ss)\n" "$hours" "$minutes" "$seconds"
