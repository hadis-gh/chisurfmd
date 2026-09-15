#!/usr/bin/env bash

set -eu


# ==============================================================================
# 1. PROJECT PATHS
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
# 2. LOAD CONFIGURATION
# ==============================================================================
config=${1:-"config.sh"}

if [ ! -f "$config" ]; then
    echo "Error: configuration file '$config' not found."
    exit 1
fi

source "$config"

# ==============================================================================
# 3. PATHS AND OUTPUT
# ==============================================================================
MD_EXE=${MD_EXE:-"$PROJECT_ROOT/build/examples/chisurfmd_md"}
DEPOSIT_MD_EXE=${DEPOSIT_MD_EXE:-"$PROJECT_ROOT/build/examples/chisurfmd_aggregation"}

outputDir=${outputDir:-"$PROJECT_ROOT/outputs/aggregation"}

mkdir -p "$outputDir"

logFile="${outputDir}/log.txt"
exec > >(tee "$logFile") 2>&1

start_time=$(date +%s)

# ==============================================================================
# 4. THERMAL EQUILIBRATION
# ==============================================================================
run_md()
{
    local temperature="$1"

    local outfile="${outputDir}/equilibrated_${temperature}.bp"

    if [ -d "$outfile" ]; then
        echo "Equilibrated state already exists: $outfile"
        return
    fi

    echo
    echo "=================================================="
    echo "Thermal equilibration"
    echo "Temperature: $temperature"
    echo "=================================================="

    local args=(
        --particlesInit "$particlesInit"
        --particlesNum "$particlesNum"

        --chirality "$chirality"
        --alignment "$alignment"

        --particleRadius "$particleRadius"
        --momentI "$momentI"

        --areaL "$areaL"
        --seed "$seed"

        --integration "$integration"
        --time "$equilibrationTime"
        --dt "$dt"

        --temperature "$temperature"
        --thermoInterval "$thermoInterval"
        --collisionFr "$collisionFr"

        --writeStateInterval "$writeStateInterval"
        --writeEnergyInterval "$writeEnergyInterval"

        --saveFile "$outfile"
    )

    if [ "$potentialType" == "OrientedLJ" ]; then
        args+=(
            --LJangularScale "$LJangularScale"
            --LJPhiOrder "$LJPhiOrder"
            --LJalpha "$LJalpha"
        )
    fi

    "$MD_EXE" "${args[@]}"
}

# ==============================================================================
# 5. PARTICLE DEPOSITION
# ==============================================================================
run_deposition()
{
    local depositionTime="$1"
    local temperature="$2"

    local inputFile="${outputDir}/equilibrated_${temperature}.bp"
    local outfile="${outputDir}/aggregated_${temperature}_${depositionTime}.bp"

    if [ ! -d "$inputFile" ]; then
        echo "Error: equilibrated state not found:"
        echo "  $inputFile"
        exit 1
    fi

    if [ -d "$outfile" ]; then
        echo "Aggregation output already exists: $outfile"
        return
    fi

    echo
    echo "=================================================="
    echo "Particle aggregation"
    echo "Temperature:        $temperature"
    echo "Deposition interval: $depositionTime"
    echo "=================================================="

    local args=(
        --particlesInit "$inputFile"
        --particlesNumMax "$particlesNumMax"

        --chirality "$chirality"
        --alignment "$alignment"

        --particleRadius "$particleRadius"
        --momentI "$momentI"

        --areaL "$areaL"
        --seed "$seed"

        --integration "$integration"
        --time "$depositionTime"
        --dt "$dt"

        --temperature "$temperature"
        --thermoInterval "$thermoInterval"
        --collisionFr "$collisionFr"

        --writeStateInterval "$writeStateInterval"
        --writeEnergyInterval "$writeEnergyInterval"

        --saveFile "$outfile"
    )

    if [ "$potentialType" == "OrientedLJ" ]; then
        args+=(
            --LJangularScale "$LJangularScale"
            --LJPhiOrder "$LJPhiOrder"
            --LJalpha "$LJalpha"
        )
    fi

    "$DEPOSIT_MD_EXE" "${args[@]}"
}

# ==============================================================================
# 6. EQUILIBRATE EACH TEMPERATURE
# ==============================================================================

echo "=================================================="
echo "        Aggregation Molecular Dynamics            "
echo "=================================================="

for temperature in "${temperatures[@]}"; do
    run_md "$temperature"
done

# ==============================================================================
# 7. PARALLEL AGGREGATION SWEEP
# ==============================================================================
export -f run_deposition

export DEPOSIT_MD_EXE
export outputDir

export particlesNumMax
export chirality
export alignment

export particleRadius
export momentI
export areaL
export seed

export integration
export dt

export thermoInterval
export collisionFr

export writeStateInterval
export writeEnergyInterval

export potentialType
export LJangularScale
export LJPhiOrder
export LJalpha


parallel \
    -j "$parallelJobs" \
    run_deposition \
    ::: "${depositionTimes[@]}" \
    ::: "${temperatures[@]}"

# ==============================================================================
# 8. EXECUTION TIME
# ==============================================================================
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))

hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))

printf \
    "\nTotal execution time: %02d:%02d:%02d\n" \
    "$hours" "$minutes" "$seconds"