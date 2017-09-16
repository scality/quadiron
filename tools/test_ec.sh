#!/bin/bash

#vflag=-v
#valgrind="valgrind --leak-check=full"
#valgrind=valgrind
#valgrind="gdb --args"

bin=./ec
bs=50K
#bs=1M

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
        for i in foo.d[0-9]
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

for i in fntrs_1 fntrs_2 gf2nfftrs_1 gf2nfftrs_2 gf2nfftrs_4 gf2nfftrs_8 gf2nrsv_1 gf2nrsv_2 gf2nrsc_1 gf2nrsc_2 gf2nrsv_4 gf2nrsv_8 gf2nrsv_16 gf2nrsc_4 gf2nrsc_8 gf2nrsc_16
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
