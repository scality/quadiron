#!/bin/bash

cd `dirname $0`

#leo
if [ -d ./wirehair/.git ]; then
    pushd wirehair
    git pull
    popd
else
    git clone git@github.com:catid/wirehair.git
fi

pushd wirehair
patch -p1 < ../wirehair.patch
mkdir build
pushd build
if [[ "$OSTYPE" == "darwin"* ]]; then
    # Mac OSX
    CXX=clang++ cmake -G 'Unix Makefiles' ..
else
    cmake -G 'Unix Makefiles' ..
fi
make
popd
popd
