#!/bin/bash

# Clean up function
clean() {
    rm -f parms.txt ang.dat fmap.ern fmap.ern.1 fmap.ern.bak subA.vdw
    rm -rf ve
}

# Test setup
set -e  # Exit on any error
TEST_NAME="FMAP Analysis Test"

# Handle clean command
if [ "$1" = "clean" ]; then
    clean
    echo "Cleaned up generated files"
    exit 0
fi

echo "Running $TEST_NAME..."

ncore=`grep '^cpu\scores' /proc/cpuinfo|uniq|awk '{print $4}'`
export OMP_NUM_THREADS=$ncore
ulimit -s unlimited

# Run the analysis
fmapb2pre data/subA.pqr 0.016 25 - - >parms.txt
#echo -1.570796 1.570796 -1.047198 >ang.dat
ln -sf ../src/maeppi/data/sspNF00072/rot.eul ang.dat
fmapb2 subA.vdw parms.txt ang.dat |tee fmap.ern

# Expected results
cat <<eof >fmap.ern.bak
vol	Sum:     7.855207e+07 Occ:  0.012010678305 Ave:      0.00715432 Net:      0.00000000
ele	Sum:     7.458214e+07 Occ:  0.012010678305 Ave:      0.03785989 Net:      0.03070557
vdw	Sum:     9.233193e+07 Occ:  0.012010678305 Ave:     -0.08854243 Net:     -0.09569676
v12	Sum:     7.853539e+07 Occ:  0.012010678305 Ave:      0.00728006 Net:      0.00012574
v+e	Sum:     7.479959e+07 Occ:  0.012010678305 Ave:      0.03613620 Net:      0.02898188

eof

head -5 fmap.ern >fmap.ern.1
echo >>fmap.ern.1

# Test 1: Compare analysis output
echo "Test 1: Checking analysis output..."
if diff fmap.ern.1 fmap.ern.bak; then
    echo "PASS: Analysis output matches expected results"
else
    echo "FAIL: Analysis output differs from expected results"
    exit 1
fi

# Test 2: Verify checksum
echo "Test 2: Checking file checksum..."
expected_md5="551c3af5055ff3f4099260b39d2d09c2"
actual_md5=$(md5sum ve/tErn.00000.mat.gz | cut -d' ' -f1)
if [ "$actual_md5" = "$expected_md5" ]; then
    echo "PASS: Checksum matches expected value"
else
    echo "FAIL: Checksum mismatch. Expected: $expected_md5, Got: $actual_md5"
    exit 1
fi

# Test 3:Check B2
echo "Test 3: Checking B2..."
fmapb2stat stats
expected_b2="14.918"
actual_b2=$(grep "FMAPb22      :" parms.txt|awk '{print $3}')
if [ "$actual_b2" = "$expected_b2" ]; then
    echo "PASS: B2 matches expected value"
else
    echo "FAIL: B2 mismatch. Expected: $expected_b2, Got: $actual_b2"
    exit 1
fi

echo "All tests passed!"

# Clean up on successful completion
#clean
