#!/bin/bash
set -e

cmake --version
make --version
g++ --version

echo ""
echo "--------------------------------"
echo "Generating protobuf file"
echo ""

mkdir -p ./src/data-formats/
cd ./data-formats/file-formats
protoc --cpp_out="../../src/data-formats" experiment.proto
cd ../..

echo ""
echo "--------------------------------"
echo "Building PIQUANT"
echo ""

# delete last built executable
rm -f ./build/Piquant

# DEBUG:
# -g -rdynamic

mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ..
make
#docker run --rm -e "GIT_BRANCH_ENV=${CI_COMMIT_BRANCH}" -e "GIT_BRANCH_GIT=${GIT_BRANCH_FROM_GIT}" -v "$PWD":/usr/src/PIQUANT -w /usr/src/PIQUANT piquant-builder ./docker/compile/runmake.sh
