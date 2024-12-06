echo "--------------------------------"
echo "Building Run-Container"
echo ""

if [ -n "$1" ];
then
    docker build --cache-from $1 -t piquant-runner .
else
    docker build -t piquant-runner .
fi
