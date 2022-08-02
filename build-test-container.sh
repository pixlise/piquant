#!/bin/bash
set -e

echo "--------------------------------"
echo "Building Piquant Test Container"
echo ""

docker build -t piquant-tester -f ./Dockerfile-Tester .
