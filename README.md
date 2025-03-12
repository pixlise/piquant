[![DOI](https://zenodo.org/badge/520582671.svg)](https://zenodo.org/badge/latestdoi/520582671)

# PIQUANT - Quantitative X-ray Fluorescence Analysis
Written for PIXL, the Planetary Instrument for X-ray Lithochemistry"

## What's here?
- `build` - Build results directory, created when building. Ignored by git.
- `data-formats` - Git sub-module containing protobuf descriptions to generate serialization code from.
- `doc` - Documentation/notes.
- `src` - Piquant C++ source code.
- `test` - Test scripts and data required to test Piquant. Includes config files, spectrum data and expected output.
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

## To compile Piquant as a Windows 11 executable
1. Install Visual Studio Build Tools: winget install Microsoft.VisualStudio.2022.BuildTools --force --override "--wait --passive --add Microsoft.VisualStudio.Component.VC.Tools.x86.x64 --add Microsoft.VisualStudio.Component.Windows11SDK.22621"
2. Install CMake from https://cmake.org/download/
3. Compile Piquant by running `.\compile.bat`

This should compile PIQUANT and place the executable in .\build\Release\Piquant.exe

## Build Container
This is a docker container that is part of the project which contains all the required build tools to build any
of our repositories. You shouldn't have to interact with it directly, as we have "local" versions of scripts that
already run things in that docker container, but for reference, it comes from the `build-container` repository.

## Testing Piquant
1. Run `./build-test-container.sh` once to make sure tests can run. If you change code/test data, this does not
need to be re-run. NOTE: if you change code, don't forget to run `./local-compile.sh` to update the executable used inside the test container!
2. Run `./test.sh`

This script:
- Builds a docker container to run piquant. It includes python 3 so the piquant-executing script can be run
- Runs tests, which are actually python unit tests. Each test assembles a command line and runs piquant.
Output result files are then read back and compared to expected output files in python.

### Running a specific test
Modify python3 unittest call in `./test/code/runtests.sh` to specify the test file name, eg: `python3 -m unittest -v test_piquant.py`

## Versioning of Piquant
When Piquant is built, the version number in the line containing `project()` in `CMakeLists.txt` is used
to generate `build/version.h` which is built into Piquant.

The git branch the compilation is done from is also required. `CMakeLists.txt` obtains this from the
environment variable `GIT_BRANCH_ENV`. If this is empty/doesn't exist, it tries `GIT_BRANCH_GIT`.

The reason for having 2 environment variables is due to needing to build both locally and in gitlab CI, and
there being different possibility to execute git commands to read the current branch (in gitlab CI, the env
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
