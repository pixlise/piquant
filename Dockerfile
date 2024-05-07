FROM debian:bullseye-slim

MAINTAINER Peter Nemere <peter.nemere@qut.edu.au>

# Apparently this is needed now...
# https://github.com/aws/aws-sdk-go/issues/2322
#RUN apk update && apk add ca-certificates && rm -rf /var/cache/apk/*
RUN apt-get update && apt-get install -yq ca-certificates

# Copy PIQUANT executable and the runner
COPY ./build/Piquant /usr/PIQUANT/Piquant
COPY ./build/PiquantRunner /usr/PIQUANT/PiquantRunner

# Make sure they're executable
RUN chmod +x /usr/PIQUANT/Piquant /usr/PIQUANT/PiquantRunner

# Start from the root
WORKDIR /usr/PIQUANT

CMD ["./PiquantRunner"]
