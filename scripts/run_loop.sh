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

MD_EXE=${MD_EXE}
outputDir=${outputDir}
mkdir -p "$outputDir"

logFile="log.txt"
exec > >(tee "$logFile") 2>&1
echo "Logging to: $logFile"

run_simulation() {
    local runIndex="$1"
    local particlesInit="$2"
    local temperature="$3"
    local saveFile="${outputDir}/run_${runIndex}.bp"

    printf "Starting simulation #%s | Temperature: %s\n" "$runIndex" "$temperature"

    local args=(
        --particlesInit "$particlesInit"
        --particlesType "$particlesType"
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

if [ "$potentialType" == "OrientedLJ" ]; then
        args+=(--LJangularScale "$LJangularScale"
               --LJPhiOrder "$LJPhiOrder"
               --LJalpha "$LJalpha"
               --momentI "$momentI")
    elif [ "$potentialType" == "GeometricLJ" ]; then
        args+=(--patchNum "$patchNum")
    elif [ "$potentialType" == "PatchyLJ" ]; then
        args+=(--patchNums "$patchNums"
               --sigmaPatch "$sigmaPatch" 
               --patchMode "$patchMode")
    elif [ "$potentialType" == "ChiralPatchyLJ" ]; then
        args+=(--sigmaPatch "$sigmaPatch" 
               --patchNums "$patchNums"
               --epsilonSame "$epsilonSame"
               --epsilonOpp "$epsilonOpp"
               --sigmaSame "$sigmaSame"
               --sigmaOpp "$sigmaOpp"
               --chiralOffset "$chiralOffset")               
    elif [ "$potentialType" == "ChiMorse" ]; then
        args+=(--chiMorseModel "$chiMorseModel")
    else
        printf "Error: Make sure to specify a valid potentialType in config file. (Got: '%s')\n" "$potentialType"
        exit 1
    fi

    "$MD_EXE" "${args[@]}" || { 
        printf "Error: Simulation #%s failed.\n" "$runIndex"; 
        exit 1; 
    }

    printf "\nSimulation %s of %s completed successfully.\n\n" "$runIndex" "$totalRuns"
    echo "------------------------------------------------------------"
}

echo "======================================================="
echo "    Temperature Loop Molecular Dynamics Simulation     "
echo "======================================================="

initialRunIndex=0
run_simulation "$initialRunIndex" "RANDOM" "$highTemperature"

# Cooling process: 
for ((i = 1; i <= numCoolingRuns; i++)); do
    temperature=$(echo "$highTemperature - $i * $stepTemperature" | bc)
    prevOutput="${outputDir}/run_$((i - 1)).bp"

    run_simulation "$i" "$prevOutput" "$temperature"
done

# Heating process:
for ((i = 1; i <= numHeatingRuns; i++)); do
    runIndex=$((numCoolingRuns + i))
    temperature=$(echo "$lowTemperature + $i * $stepTemperature" | bc)
    prevOutput="${outputDir}/run_$((runIndex - 1)).bp"
    run_simulation "$runIndex" "$prevOutput" "$temperature"
done

echo "All simulations completed. Total runs: ${totalRuns}"

# Record time
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))

hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))

printf "Total execution time: %02d:%02d:%02d (hh:mm:ss)\n" "$hours" "$minutes" "$seconds"