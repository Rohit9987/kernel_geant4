#!/bin/bash
energies=$(seq 0.1 0.1 1.0)

for E in $energies
do
    energy="${E}MeV"
    echo "Running simulation for energy: $energy"
    ./exampleB4c -t 8 -m run2.mac -e $energy
done

