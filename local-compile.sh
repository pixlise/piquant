# We pass a branch name from git, and from an env var CI_COMMIT_BRANCH. This is because the gitlab CI process doesn't give us git, it only has CI_COMMIT_BRANCH env var.
# In the test script we make a choice to use env var (for gitlab CI), if that's not set, use what git gave (this is for locally running)
if [ -n "$(which git)" ];
then
    GIT_BRANCH_FROM_GIT=$(git rev-parse --abbrev-ref HEAD)
fi
export MSYS_NO_PATHCONV=1
docker run --rm -e "GIT_BRANCH_GIT=${GIT_BRANCH_FROM_GIT}" -v "$PWD":/usr/src/PIQUANT -w /usr/src/PIQUANT ghcr.io/pixlise/build-container:golang-1.18-protoc-3.7.1-protobuf-3.11.4-angular-13.1.2-nodejs-16 /bin/bash compile.sh
