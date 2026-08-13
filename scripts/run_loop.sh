#!/usr/bin/env bash

set -eu

start_time=$(date +%s)

config=${1:-"config.sh"}
if [ ! -f "$config" ]; then
    echo "Error: '$config' not found."
    exit 1
fi
source "$config"

numCoolingRuns=$(echo "scale=0; ($highTemperature - $lowTemperature) / $stepTemperature" | bc)

if (( $(echo "$timeHeating > 0" | bc -l) )); then
    numHeatingRuns="$numCoolingRuns"
else
    numHeatingRuns=0
fi

totalRuns=$((1 + numCoolingRuns + numHeatingRuns))

outputDir=${outputDir:-"outputs"}
mkdir -p "$outputDir"

logFile="log.txt"
exec > >(tee "$logFile") 2>&1
echo "Logging to: $logFile"


# ---------------------- General run function ----------------------

run_simulation() {
    local runIndex="$1"
    local particlesInit="$2"
    local temperature="$3"
    local duration="$4"
    local saveFile="${outputDir}/run_${runIndex}.bp"

    printf "Starting simulation %s/%s | Temperature: %s\n" \
        "$((runIndex + 1))" "$totalRuns" "$temperature"

    local args=(
        --particlesInit "$particlesInit"
        --temperature "$temperature"
        --time "$duration"
        --dt "$dt"
        --saveFile "$saveFile"
        --particlesNum "$particleNum"
        --chirality "$chirality"
        --alignment "$alignment"
        --seed "$((seed + runIndex))"
        --areaL "$areaL"
        --particleRadius "$particleRadius"
        --fixRadius "$fixRadius"
        --thermoInterval "$thermoInterval"
        --writeStateInterval "$writeStateInterval"
        --writeEnergyInterval "$writeEnergyInterval"
        --integration "$integration"
        --collisionFr "$collisionFr"
    )

    # ----------------------- Potential options -----------------------

    if [ "$potentialType" == "OrientedLJ" ]; then
        args+=(
            --LJangularScale "$LJangularScale"
            --LJPhiOrder "$LJPhiOrder"
            --LJalpha "$LJalpha"
            --momentI "$momentI"
        )

    elif [ "$potentialType" == "GeometricLJ" ]; then
        args+=(--patchNum "$patchNum")

    elif [ "$potentialType" == "PatchyLJ" ]; then
        args+=(
            --patchNums "$patchNums"
            --sigmaPatch "$sigmaPatch"
            --patchMode "$patchMode"
        )

    elif [ "$potentialType" == "ChiralPatchyLJ" ]; then
        args+=(
            --sigmaPatch "$sigmaPatch"
            --patchNums "$patchNums"
            --epsilonSame "$epsilonSame"
            --epsilonOpp "$epsilonOpp"
            --sigmaSame "$sigmaSame"
            --sigmaOpp "$sigmaOpp"
            --chiralOffset "$chiralOffset"
        )

    elif [ "$potentialType" == "ChiMorse" ]; then
        args+=(
            --chiMorseModel "$chiMorseModel"
            --E0 "$E0"
            --L0 "$L0"
        )

    elif [ "$potentialType" == "Tabular" ]; then
        args+=(
            --dataReferenceDir "$dataReferenceDir"
            --E0 "$E0"
            --L0 "$L0"
        )

    else
        printf "Error: invalid potentialType '%s'.\n" "$potentialType"
        exit 1
    fi

    "$MD_EXE" "${args[@]}" || {
        printf "Error: Simulation #%s failed.\n" "$runIndex"
        exit 1
    }

    printf "\nSimulation %s/%s completed successfully.\n\n" \
        "$((runIndex + 1))" "$totalRuns"

    echo "------------------------------------------------------------"
}

echo "======================================================="
echo "    Temperature Loop Molecular Dynamics Simulation     "
echo "======================================================="


# ------------------ Initial high-temperature run -----------------

run_simulation 0 "RANDOM" "$highTemperature" "$timeCooling"


# ------------------------- Cooling -------------------------------

for ((i = 1; i <= numCoolingRuns; i++)); do

    temperature=$(echo "$highTemperature - $i * $stepTemperature" | bc)
    prevOutput="${outputDir}/run_$((i - 1)).bp"

    duration="$timeCooling"

    # Longer relaxation in the ordering window
    if (( $(echo "$temperature >= $criticalLow && \
                  $temperature <= $criticalHigh" | bc -l) )); then
        duration="$timeCritical"
    fi

    run_simulation "$i" "$prevOutput" "$temperature" "$duration"
done

# ------------------------- Heating -------------------------------

for ((i = 1; i <= numHeatingRuns; i++)); do
    runIndex=$((numCoolingRuns + i))
    temperature=$(echo "$lowTemperature + $i * $stepTemperature" | bc)
    prevOutput="${outputDir}/run_$((runIndex - 1)).bp"

    run_simulation "$runIndex" "$prevOutput" "$temperature" "$timeHeating"
done


echo "All simulations completed. Total runs: $totalRuns"


# ---------------------- Execution time ---------------------------

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))

hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))

printf "Total execution time: %02d:%02d:%02d (hh:mm:ss)\n" \
    "$hours" "$minutes" "$seconds"