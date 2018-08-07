FROM ubuntu:18.04
MAINTAINER Sylvain Laperche "sylvain.laperche@scality.com"

RUN apt-get update && apt-get -y --no-install-recommends install \
    ca-certificates \
    cmake \
    git \
    make

# Copy a standalone buildchain from the build context.
#
# To generate a standalone buildchain:
# 1. Download the Android NDK from https://developer.android.com/ndk/downloads/
# 2. Generate a toolchain using `build/tools/make_standalone_toolchain.py`
#    e.g: build/tools/make_standalone_toolchain.py --arch arm64 --install-dir /arm64
COPY arm64 /arm64

# Add the standalone toolchain to the search path.
ENV PATH $PATH:/arm64/bin

# Define what tools to use.
ENV AR aarch64-linux-android-ar
ENV AS aarch64-linux-android-clang
ENV CC aarch64-linux-android-clang
ENV CXX aarch64-linux-android-clang++
ENV LD aarch64-linux-android-ld
ENV STRIP aarch64-linux-android-strip
