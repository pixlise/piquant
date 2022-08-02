#!/bin/bash
set -e

echo "--------------------------------"
echo "Run Piquant Tests"
echo ""

# NOTE: we pass in the git branch as GIT_BRANCH env variable - this is only required by the python test that checks the piquant version output matches
# what we're expecting from CMakeLists.txt and the git branch
# We pass a branch name from git, and from an env var CI_COMMIT_BRANCH. This is because the gitlab CI process doesn't give us git, it only has CI_COMMIT_BRANCH env var.
# In the test script we make a choice to use env var (for gitlab CI), if that's not set, use what git gave (this is for locally running)
if [ -n "$(which git)" ];
then
    GIT_BRANCH_FROM_GIT=$(git rev-parse --abbrev-ref HEAD)
fi

CONTAINER=piquant-tester

if test -z "$AWS_DEFAULT_REGION"
then
    AWS_DEFAULT_REGION="us-east-1"
fi

mkdir -p ./test/output

docker run --rm -e BUILD_VERSION=$BUILD_VERSION -e "GIT_BRANCH_ENV=${CI_COMMIT_BRANCH}" -e "GIT_BRANCH_GIT=${GIT_BRANCH_FROM_GIT}" -e AWS_ACCESS_KEY_ID=$AWS_ACCESS_KEY_ID -e AWS_SECRET_ACCESS_KEY=$AWS_SECRET_ACCESS_KEY -e AWS_DEFAULT_REGION=$AWS_DEFAULT_REGION -v "$HOME/.aws/credentials":/root/.aws/credentials -v "$PWD/":/build -w /build $CONTAINER ./test/code/runtests.sh
