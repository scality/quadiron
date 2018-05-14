FROM ubuntu:18.04
MAINTAINER Sylvain Laperche "sylvain.laperche@scality.com"

RUN apt-get update && apt-get -y --no-install-recommends install \
    ca-certificates \
    clang \
    clang-tidy \
    clang-format \
    cmake \
    doxygen \
    git \
    make

# Fix binary names.
RUN ln -s /usr/bin/clang-format-6.0 /usr/bin/clang-format &&\
    ln -s /usr/bin/clang-tidy-6.0   /usr/bin/clang-tidy
