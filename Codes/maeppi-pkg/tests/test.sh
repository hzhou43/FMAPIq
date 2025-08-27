#!/bin/bash

if [ ! -s parms.txt ];then
    echo Please first prepare fmapb2 energy file with run.fmapb2.sh
    exit 1
fi 

# Clean up function
clean() {
    rm -f final_config.dat maeppi_sim.log maeppi_ern.log ern.log maeppi_shm_server.log 
    maeppi_shm_server cleanup parms.txt &>/dev/null || true
}

# Test setup
set -e
TEST_NAME="MAEPPI Simulation Test"

# Handle clean command
if [ "$1" = "clean" ]; then
    clean
    echo "Cleaned up generated files"
    exit 0
fi

echo "Running $TEST_NAME..."

# Clean up any previous run
clean

ncore=`grep '^cpu\scores' /proc/cpuinfo|uniq|awk '{print $4}'`
export OMP_NUM_THREADS=1 #$ncore

# Extract parameters
kbt=$(grep kbt parms.txt | awk '{print $3}')

# Start server
echo "1. Starting shared memory server..."
maeppi_shm_server start parms.txt --load-only #|tee maeppi_shm_server.log
#keyp=$(grep keyp maeppi_shm_server.log | awk '{print $3}')
keyp=`maeppi_shm_server status parms.txt --key|tail -1`

# Run simulation
echo "2. Running simulation..."
maeppi_sim 200 619.65 $kbt 4000000 1 - final_config.dat 400000 $keyp | tee maeppi_sim.log

# Calculate energy
echo "3. Calculating energy..."
maeppi_ern final_config.dat 619.65 $keyp >maeppi_ern.log

# Cleanup server
echo "4. Cleaning up server..."
maeppi_shm_server cleanup parms.txt

# Test results
echo "5. Checking results..."
grep -i "energy" maeppi_*.log | sed 's/Speed.*$//' > ern.log

if diff ern.log data/ern.log > /dev/null; then
    echo "PASS: Energy results match expected values"
else
    echo "FAIL: Energy results differ from expected"
    echo "Differences:"
    diff ern.log data/ern.log
    exit 1
fi

echo "All tests passed!"

# Clean up on successful completion
clean
