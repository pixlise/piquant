#!/bin/bash

# Set up a dir structure for the test to work off in the container
rm -rf /testenv
mkdir -p /testenv
cp -R ./test/data/ /testenv/test-data
cp -R ./test/code/* /testenv/
cp ./build/PiquantRunner /testenv/
cp ./build/Piquant /testenv/
cp CMakeLists.txt /testenv/

# Run tests in the new dir, outputs to be written to /build/test/output
cd /testenv

# Print piquant version for clarity/test logs
echo "Piquant version: "$(./Piquant version)

python3 -m unittest -v #add file names here to run specific tests, eg test_piquant_map.py
