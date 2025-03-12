rmdir /S /Q build
mkdir build
pushd build

cmake -D PIQUANT_EXCLUDE_PROTOBUF=true -G "Visual Studio 17 2022" -A x64 ..
msbuild /p:configuration=release /p:NoWarn="C4305;C4267;C4244;C4018" Piquant.vcxproj

popd
