#!/bin/bash

# Clean up function
clean() {
    rm -f input.pqr sim2iq.dat
}

# Simple test script
set -e

# Handle clean command
if [ "$1" = "clean" ]; then
    clean
    echo "Cleaned up generated files"
    exit 0
fi

# Set number of threads
ncore=`grep '^cpu\scores' /proc/cpuinfo|uniq|awk '{print $4}'`
export OMP_NUM_THREADS=$ncore

echo "Running tests..."

# Test 1: Check input processing
echo "Test 1: Input processing"
sim2iqpre -f data/input.pqr -o input.pqr
diff input.pqr data/input.pqr
echo "✓ Input processing test passed"

# Test 2: Check simulation output with float comparison
echo "Test 2: Simulation output"
sim2iq input.pqr 620 619.65 sim2iq.dat --steps 4 --traj data/snap.txt

# Simple float comparison for multi-column data
compare_floats() {
    local file1=$1
    local file2=$2
    local tolerance=1e-6
    
    if awk -v tol=$tolerance '
    NR==FNR {
        for(i=1; i<=NF; i++) ref[NR,i]=$i
        refcols[NR]=NF
        next
    }
    {
        for(i=1; i<=NF; i++) {
            val1 = $i
            val2 = ref[FNR,i]
            maxval = (abs(val1) > abs(val2)) ? abs(val1) : abs(val2)
            if (maxval > 0) {
                rel_diff = abs(val1 - val2) / maxval
                if (rel_diff > tol) {
                    print "Files differ at line " FNR " column " i
                    exit 1
                }
            } else if (val1 != val2) {
                print "Files differ at line " FNR " column " i
                exit 1
            }
        }
    }
    function abs(x) {return x < 0 ? -x : x}
    ' "$file2" "$file1"; then
        echo "Float comparison passed"
    else
        return 1
    fi
}

compare_floats sim2iq.dat data/sim2iq.dat
echo "✓ Simulation output test passed"

echo "All tests passed!"

# Clean up on successful completion
clean
