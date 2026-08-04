#!/usr/bin/env bash

set -euo pipefail

config=${1:-config.sh}

if [[ ! -f "$config" ]]; then
    echo "Error: config file '$config' not found."
    exit 1
fi

source "$config"

MD_EXE=${MD_EXE:-"../build/examples/testMD"}
outputDir=${outputDir:-"./outputs"}
mkdir -p "$outputDir"

logFile="${outputDir}/annealing.log"
exec > >(tee "$logFile") 2>&1

numRuns=$(awk \
    -v high="$highTemperature" \
    -v low="$lowTemperature" \
    -v step="$stepTemperature" \
    'BEGIN { printf "%d", (high - low) / step }')

totalRuns=$((2 * numRuns + 1))

run_simulation() {
    local runIndex=$1
    local particlesInit=$2
    local temperature=$3
    local runTime=$4
    local saveFile="${outputDir}/run_$(printf '%03d' "$runIndex").bp"

    echo "------------------------------------------------------------"
    printf 'Run %d/%d | T=%s\n' \
        "$((runIndex + 1))" "$totalRuns" "$temperature"

    if [[ -d "$saveFile" ]]; then
        echo "Skipping existing output: $saveFile"
        return
    fi

    args=(
        --particlesInit "$particlesInit"
        --particlesType "$particlesType"
        --temperature "$temperature"
        --time "$runTime"
        --dt "$dt"
        --saveFile "$saveFile"
        --particlesNum "$particleNum"
        --seed "$seed"
        --areaL "$areaL"
        --exclusionRadius "$exclusionRadius"
        --fixRadius "$fixRadius"
        --momentI "$momentI"
        --thermoInterval "$thermoInterval"
        --writeStateInterval "$writeStateInterval"
        --writeEnergyInterval "$writeEnergyInterval"
        --integration "$integration"
        --collisionFr "$collisionFr"
    )

    case "$potentialType" in
        ChiMorse)
            args+=(--chiMorseModel "$chiMorseModel")
            ;;
        OrientedLJ)
            args+=(
                --LJangularScale "$LJangularScale"
                --LJPhiOrder "$LJPhiOrder"
                --LJalpha "$LJalpha"
            )
            ;;
        *)
            echo "Error: unsupported potentialType '$potentialType'."
            exit 1
            ;;
    esac

    "$MD_EXE" "${args[@]}"
}

temperature_at() {
    awk \
        -v start="$1" \
        -v step="$2" \
        -v i="$3" \
        'BEGIN { printf "%.12g", start + step * i }'
}

echo "============================================================"
echo "             Temperature Annealing Simulation"
echo "============================================================"
echo "Potential: $potentialType"
echo "Temperature: $highTemperature -> $lowTemperature -> $highTemperature"
echo "Total runs: $totalRuns"
echo "Output: $outputDir"
echo "============================================================"

runIndex=0

run_simulation \
    "$runIndex" \
    "RANDOM" \
    "$highTemperature" \
    "$timeCooling"

previousOutput="${outputDir}/run_000.bp"

for ((i = 1; i <= numRuns; i++)); do
    ((runIndex += 1))

    temperature=$(temperature_at \
        "$highTemperature" \
        "-$stepTemperature" \
        "$i")

    run_simulation \
        "$runIndex" \
        "$previousOutput" \
        "$temperature" \
        "$timeCooling"

    previousOutput="${outputDir}/run_$(printf '%03d' "$runIndex").bp"
done

for ((i = 1; i <= numRuns; i++)); do
    ((runIndex += 1))

    temperature=$(temperature_at \
        "$lowTemperature" \
        "$stepTemperature" \
        "$i")

    run_simulation \
        "$runIndex" \
        "$previousOutput" \
        "$temperature" \
        "$timeHeating"

    previousOutput="${outputDir}/run_$(printf '%03d' "$runIndex").bp"
done

echo "============================================================"
echo "All $totalRuns simulations completed."
echo "Final output: $previousOutput"
echo "============================================================"