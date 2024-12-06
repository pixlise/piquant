# PIQUANT - Quantitative X-ray Fluorescence Analysis
Written for PIXL, the Planetary Instrument for X-ray Lithochemistry"

## What's here?
- `build` - Build results directory, created when building. Ignored by git.
- `data-formats` - Git sub-module containing protobuf descriptions to generate serialization code from.
- `doc` - Documentation/notes.
- `src` - Piquant C++ source code.
- `test-data` - Data required to test Piquant. Includes config files, spectrum data and expected output.
- `.github` - See description below.
- `CMakeLists.txt` - CMake file to build Piquant. See description below.
- `version.h.in` - Template for `version.h` that gets created during build process. See CMake description below.

## To compile Piquant for running in Docker
1. Install Docker.
2. Make sure git submodules are available by running `git submodule init` then `git submodule update`
3. Compile Piquant by running `./local-compile.sh`

This script starts the `compile.sh` locally on your machine by running our build container in Docker which has
all the build tools pre-installed. The compile script:
- Runs protoc to generate protobuf serialization code for reading Pixlise binary files
- Runs CMake to generate a make file for Piquant
- Finally, compiles piquant (using make) into an executable for linux to /build/Piquant

## Build Container
This is a docker container that is part of the project which contains all the required build tools to build any
of our repositories. You shouldn't have to interact with it directly, as we have "local" versions of scripts that
already run things in that docker container, but for reference, it comes from the `build-container` repository.

## Testing Piquant
1. Set up your AWS credentials so you can run tests, which need to interact with AWS S3. You'll need an AWS key and
secret, so you can put them in a credentials file, normally located in $HOME/.aws/credentials. This file is pulled
into the piquant runner docker container when running tests.
2. Run `./build-test-container.sh` once to make sure tests can run. If you change code/test data, this does not
need to be re-run. NOTE: if you change code, don't forget to run `./local-compile.sh` to update the executable used inside the test container!
3. Run `./test.sh`

This script:
- Builds a docker container to run piquant. It includes python 3 so the piquant-executing script can be run
- Runs tests, which are actually python unit tests. Each test assembles a command line and runs piquant.
Output result files are then read back and compared to expected output files in python.

### Running a specific test
Modify python3 unittest call in `./test/code/runtests.sh` to specify the test file name, eg: `python3 -m unittest -v test_piquant.py`

## Running Piquant
Piquant is a C++ program, but to execute it in docker in AWS, and using config/MSA from AWS S3, we have
a PiquantRunner executable included in its docker container. This downloads the required files for
Piquant and puts them wher ethey are needed. The runner is written in Go and also needs to be compiled.

## Compiling PiquantRunner
This is a small program written in Go (www.golang.org) and also needs to be compiled by running
`local-runner-compile.sh`.

To run its unit tests, run `local-runner-test.sh`.

## Github Actions
`.github/` contains the github actions required to build/test Piquant. It basically runs the above
described `compile.sh` and `test.sh` but with some extra commands around it to allow docker-in-docker to
work and to help speed up builds by caching the built docker container for next time.

## Versioning of Piquant
When Piquant is built, the version number in the line containing `project()` in `CMakeLists.txt` is used
to generate `build/version.h` which is built into Piquant.

The git branch the compilation is done from is also required. `CMakeLists.txt` obtains this from the
environment variable `GIT_BRANCH_ENV`. If this is empty/doesn't exist, it tries `GIT_BRANCH_GIT`.

The reason for having 2 environment variables is due to needing to build both locally and in github actions, and
there being different possibility to execute git commands to read the current branch (in github actions, the env
variable is the more dependable solution).

Running `Piquant version` will output the version and git branch it came from, for example: `3.2.1-branch`

To change the version number, edit `CMakeLists.txt`.

## CMake
CMake is used to generate make files for Piquant builds. The simple description of the source and what
compiler flags are required in `CMakeLists.txt` allow us to generate either make files for linux builds
or if required, Visual Studio solution/project files so we can compile/test Piquant in there.

More can be found at: cmake.org

## Protobuf
To help provide a more compact/compressed/reliable data format for use in Pixlise, the required overall dataset
spectrum/location CSV/housekeeping CSV files are converted into a single binary file. The description of this
format is contained in `data-formats/experiment.proto`.

More about the format can be found at https://developers.google.com/protocol-buffers

To start Piquant faster, it will be able to read from the binary format used by Pixlise, as the overhead of
downloading 1000's of individual MSA files is quite large.


## Further Reading and Help

Check out the [Piquant Wiki](https://github.com/pixlise/piquant/wiki)
