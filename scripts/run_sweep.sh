#!/bin/bash

# Parametric sweep: radius x strain rate
# All reach 20% strain
# strain = erate * dt * nsteps => nsteps = 0.20 / (erate * 0.001)

echo "Starting parametric sweep..."
echo "RAD  RATE       NSTEPS" 

# Strain rates: 0.0001, 0.0005, 0.001, 0.005
# Radii: 3, 5, 7, 10 lattice units

for RAD in 3 5 7 10; do
    for RATE in 0.0001 0.0005 0.001 0.005; do
        # Calculate steps for 20% strain: 0.20 / (RATE * 0.001)
        NSTEPS=$(python3 -c "print(int(0.20 / ($RATE * 0.001)))")
        echo "RAD=$RAD  RATE=$RATE  NSTEPS=$NSTEPS"
        lmp -in in.si_parametric \
            -var RAD $RAD \
            -var RATE $RATE \
            -var NSTEPS $NSTEPS \
            > ../outputs/log_R${RAD}_E${RATE}.txt 2>&1
        echo "  -> Done"
    done
done

echo "All runs complete."
