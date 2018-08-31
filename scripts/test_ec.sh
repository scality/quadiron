#!/bin/bash

# Copyright 2017-2018 Scality
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors
#    may be used to endorse or promote products derived from this software
#    without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

#vflag=-v
#valgrind="valgrind --leak-check=full"
#valgrind=valgrind
#valgrind="gdb --args"

bin=$1
bs=51200

if [ -z $1 ]
then
    1>&2 echo "error: missing argument"
    echo "usage: $0 PATH_TO_EC"
    exit 1;
fi

checkfail()
{
    if [ $? -ne 0 ]
    then
        echo $*: failure
        exit 1;
    fi
}

do_test()
{
    directive=$1
    fec_type=$2
    word_size=$3
    n_data=$4
    n_coding=$5
    data_loss=$6
    coding_loss=$7
    extraopts=$8

    type=`${bin} -e ${fec_type} -w ${word_size} -n ${n_data} -m ${n_coding} -p foo -c ${extraopts} -t`
    echo -n ${fec_type}_W${word_size}_N${n_data}_M${n_coding}_D${data_loss// /-}_C${coding_loss// /-}_${type},

    rm -f foo.*

    for i in `seq 0 $(expr ${n_data} - 1)`
    do
        dd if=/dev/urandom of=foo.d${i} bs=${bs} count=1 > /dev/null 2>&1
        md5sum foo.d${i} > foo.d${i}.md5sum.1
    done

    echo -n "GEN,"
    ${valgrind} ${bin} -e ${fec_type} -w ${word_size} -n ${n_data} -m ${n_coding} -p foo -c ${extraopts} ${vflag}
    checkfail "coding generation"

    for i in `seq 0 $(expr ${n_coding} - 1)`
    do
        md5sum foo.c${i} > foo.c${i}.md5sum.1
    done

    if [ "${type}" = "type_2" ]
    then
        #remove all data
        shopt -s extglob
        files=`ls foo'.'d+([0-9])`
        shopt -u extglob
        for i in $files
        do
            mv ${i} ${i}.1
        done
    fi

    j=0
    for i in $data_loss
    do
        if [ "${type}" = "type_2" ]
        then
            mv foo.c${j} foo.c${j}.1
            mv foo.c${j}.props foo.c${j}.props.1
        else
            mv foo.d${i} foo.d${i}.1
        fi
        j=`expr ${j} + 1`
    done

    for i in $coding_loss
    do
        if [ "${type}" = "type_2" ]
        then
            mv foo.c${j} foo.c${j}.1
            mv foo.c${j}.props foo.c${j}.props.1
        else
            mv foo.c${i} foo.c${i}.1
            mv foo.c${i}.props foo.c${i}.props.1
        fi
        j=`expr ${j} + 1`
    done

    if [ "${directive}" == "enconly" ]
    then
        echo
        return
    fi

    echo -n "REP,"
    ${valgrind} ${bin} -e ${fec_type} -w ${word_size} -n ${n_data} -m ${n_coding} -p foo -r ${extraopts} ${vflag}
    checkfail "repairing"

    for i in `seq 0 $(expr ${n_data} - 1)`
    do
        md5sum foo.d${i} > foo.d${i}.md5sum.2
        diff foo.d${i}.md5sum.1 foo.d${i}.md5sum.2
        checkfail "data files mismatch"
    done
    for i in `seq 0 $(expr ${n_coding} - 1)`
    do
        md5sum foo.c${i} > foo.c${i}.md5sum.2
        diff foo.c${i}.md5sum.1 foo.c${i}.md5sum.2
        checkfail "coding files mismatch"
    done

    echo
}

for i in rs-fnt_1 rs-fnt_2 rs-nf4_2 rs-nf4_4 rs-nf4_8 rs-gfp-fft_1 rs-gfp-fft_2 rs-gfp-fft_4 rs-gf2n-fft_1 rs-gf2n-fft_2 rs-gf2n-fft_4 rs-gf2n-fft_8 rs-gf2n-fft-add_1 rs-gf2n-fft-add_2 rs-gf2n-fft-add_4 rs-gf2n-fft-add_8 rs-gf2n-v_1 rs-gf2n-v_2 rs-gf2n-c_1 rs-gf2n-c_2 rs-gf2n-v_4 rs-gf2n-v_8 rs-gf2n-v_16 rs-gf2n-c_4 rs-gf2n-c_8 rs-gf2n-c_16
do
    fec_type=$(echo $i|cut -d_ -f1)
    word_size=$(echo $i|cut -d_ -f2)

    do_test enconly ${fec_type} ${word_size} 50 50 "" ""
    do_test all ${fec_type} ${word_size} 3 3 "" ""
    do_test all ${fec_type} ${word_size} 3 3 "0 1" "0"
    do_test all ${fec_type} ${word_size} 3 5 "0 1" "0"
    do_test all ${fec_type} ${word_size} 3 3 "1 2" "2"
    do_test all ${fec_type} ${word_size} 9 3 "1 2" "2"
    do_test all ${fec_type} ${word_size} 9 3 "2 3" "2"
    do_test all ${fec_type} ${word_size} 9 5 "2 3 4" "2 3"
    do_test all ${fec_type} ${word_size} 9 5 "1 3 5" "1 3"
    do_test all ${fec_type} ${word_size} 9 5 "1 3 5 7 8" ""
    do_test all ${fec_type} ${word_size} 9 5 "" "0 1 2 3 4"
done
