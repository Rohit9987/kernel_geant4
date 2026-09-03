#!/bin/bash

energies=(
  0.1 0.2 0.3 0.4 0.5 0.6 0.8
  1.0 1.25 1.5 2.0 3.0 4.0 5.0 6.0
)

for E in $energies
do
    energy="${E}MeV"
    echo "Running simulation for energy: $energy"
    ./exampleB4c -t 8 -m run2.mac -e $energy
done

