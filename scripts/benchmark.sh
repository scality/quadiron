#!/bin/bash

bin=$1
samples_nb=100
chunk_size=512
sce_type=enc_dec
show_type=1
threads_nb=1

if [ -z $1 ]
then
    1>&2 echo "error: missing argument"
    echo "usage: $0 PATH_TO_BENCH_DRIVER"
    exit 1;
fi

for ec_type in rs-gf2n-v rs-gf2n-c rs-gf2n-fft rs-gf2n-fft-add rs-gfp-fft rs-fnt rs-nf4; do
  for k in 5; do
    for m in 2; do
      for word_size in 1 2 4 8; do
        ${bin} -e ${ec_type} -w ${word_size} -k ${k} -m ${m} -c ${chunk_size} -s ${sce_type} -g ${threads_nb} -f ${show_type}
        show_type=0
      done
    done
  done
done

chunk_size=51200
sce_type=enc_only
for word_size in 2 8; do
  for ec_type in rs-nf4 rs-fnt rs-gfp-fft; do
    for k in 16; do
      for m in 64; do
        ${bin} -e ${ec_type} -w ${word_size} -k ${k} -m ${m} -c ${chunk_size} -s ${sce_type} -g ${threads_nb} -f ${show_type}
      done
    done
  done
done
