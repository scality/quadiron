#!/bin/bash

cd `dirname $0`

mkdir -p usr/lib
mkdir -p usr/include

USR_DIR="`pwd`/usr"
echo $USR_DIR
#isa-l
if [ -d ./isa-l/.git ]; then
    pushd isa-l
    git pull
    popd
else
    git clone https://github.com/01org/isa-l.git
fi

pushd isa-l
./autogen.sh
if [[ "$OSTYPE" == "darwin"* ]]; then
    # Mac OSX
    CC=clang ./configure --target=darwin AS=yasm --prefix=$USR_DIR --libdir=$USR_DIR/lib
else
    ./configure --prefix=$USR_DIR --libdir=$USR_DIR/lib
fi
make
make install
popd
