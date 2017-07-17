#!/bin/bash

for i in src test tools
do
    cpplint $i/*.cpp
    cpplint $i/*.h
done
