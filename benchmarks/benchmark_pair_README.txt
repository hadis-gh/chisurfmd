Pair-evaluation benchmark for ChiMorse and TabularDFT
===================================================

1. Copy benchmark_pair_evaluation.cpp into the MD repository's benchmarks/ folder.
2. The two includes at the top assume these headers exist:
     chisurfmd/potential/ChiMorse.h
     chisurfmd/potential/TabularDFT.h
   Adjust those paths if your repo differs.
3. Compile as a separate executable using the same C++17 toolchain, dependency
   includes, and Release (-O3 -DNDEBUG) settings as your MD project.
   Example, *if* your existing include paths resolve all dependencies:

     g++ -std=c++17 -O3 -DNDEBUG -Iinclude \
       benchmarks/benchmark_pair_evaluation.cpp -o benchmark_pair_evaluation

   You may need your project's usual include/link settings for Boost and
   nlohmann/json; copying the MD target's settings is the easiest approach.

4. Run from a directory with access to four reference files and a JSON model:

   ./benchmark_pair_evaluation /path/to/reference_data /path/to/model.json

   Optional parameters:
   ./benchmark_pair_evaluation DATA_DIR MODEL.json N_CASES PASSES TRIALS
   Default: 2000 configurations, 10 passes, 5 trials.

5. The .dat directory must contain:
   E_all_EP.dat, E_all_EA.dat, E_all_OP.dat, E_all_OA.dat

6. The code currently samples 7-11 Angstrom and angles uniformly over
   [0, 2pi), with EP/EA/OP/OA equally represented. Adjust minR_A and
   maxR_A at the top of the file to a radial range COVERED by the reference
   table and below the cutoffs of both models. It uses E0=0.4 eV and
   L0=7.0 Angstrom, matching the defaults in the provided headers.

7. Compare the median nanoseconds per call. Ratio=Table/Chi; >1 means
   ChiMorse was faster for this specific operation and implementation.
   The first test times the scalar energy+derivatives kernels (including
   torque-on-2 conversion); the second tests the current energy-only paths.
   Neither test includes model loading, MD integration, file output, nor
   Cartesian force projection. Do NOT call this total MD performance.

8. Repeat using the full and compact JSON models, without changing other
   settings. Save the console output plus compiler/CPU details. Do not change
   library implementations before recording these baseline results.

Limitations: This file was prepared from the three user-provided headers and
mainMD.cpp. It has not been compiled against the complete project, since
its ChiSurfMD headers, data files, and dependency setup were not uploaded.
