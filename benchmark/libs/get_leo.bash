#!/bin/bash

cd `dirname $0`

#leo
if [ -d ./leopard/.git ]; then
    pushd leopard
    git pull
    popd
else
    git clone git@github.com:catid/leopard.git
fi

pushd leopard
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
