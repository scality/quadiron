#!/bin/bash

bin=./benchmark
samples_nb=100
chunk_size=50K

for ec_type in gf2nrsv gf2nrsc gf2nfftrs gf2nfftaddrs gfpfftrs fntrs ngff4rs; do
  for k in 5; do
    for m in 2; do
      for word_size in 1 2 4 8; do
        ${valgrind} ${bin} -e ${ec_type} -w ${word_size} -k ${k} -m ${m} -c ${chunk_size}
      done
    done
  done
done
