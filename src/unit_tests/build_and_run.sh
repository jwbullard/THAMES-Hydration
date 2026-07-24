#!/bin/bash
# Compile and run the standalone rate-law unit tests.
# Each test has zero external dependencies beyond its target sources
# and PhysicalConstants.h.
set -e
cd "$(dirname "$0")"

echo "=== test_nucleation_rate ==="
clang++ -std=c++17 -O2 -Wall \
    -I../thameslib \
    -o test_nucleation_rate \
    test_nucleation_rate.cc ../thameslib/NucleationRate.cc
./test_nucleation_rate

echo ""
echo "=== test_saturating_rate ==="
clang++ -std=c++17 -O2 -Wall \
    -I../thameslib \
    -o test_saturating_rate \
    test_saturating_rate.cc ../thameslib/SaturatingRate.cc
./test_saturating_rate
