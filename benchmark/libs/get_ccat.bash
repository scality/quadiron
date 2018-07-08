#!/bin/bash

cd `dirname $0`

#leo
if [ -d ./CauchyCaterpillar/.git ]; then
    pushd CauchyCaterpillar
    git pull
    popd
else
    git clone git@github.com:catid/CauchyCaterpillar.git
fi

pushd CauchyCaterpillar
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
