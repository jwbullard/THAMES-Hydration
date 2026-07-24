#!/bin/bash
# Compile and run the standalone CNT sanity test.
# Zero dependencies beyond the CNT free functions.
set -e
cd "$(dirname "$0")"
clang++ -std=c++17 -O2 -Wall \
    -I../thameslib \
    -o test_nucleation_rate \
    test_nucleation_rate.cc ../thameslib/NucleationRate.cc
./test_nucleation_rate
