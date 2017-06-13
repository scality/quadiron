#!/bin/sh

#vflag=-v
#valgrind="valgrind --leak-check=full"
#valgrind=valgrind
#valgrind="gdb --args"

bs=100K
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
    bin=$1
    gf_type=$2
    n_data=$3
    n_coding=$4
    data_loss=$5
    coding_loss=$6
    extraopts=$7
    echo ${bin} e=${gf_type} n=${n_data} m=${n_coding} data_loss=\"${data_loss}\" coding_loss=\"${coding_loss}\" ${extraopts}

    rm -f foo.*
    
    for i in `seq 0 $(expr ${n_data} - 1)`
    do
        dd if=/dev/urandom of=foo.d${i} bs=${bs} count=1 > /dev/null 2>&1
        md5sum foo.d${i} > foo.d${i}.md5sum.1
    done

    echo "coding generation"
    ${valgrind} ${bin} -e ${gf_type} -n ${n_data} -m ${n_coding} -p foo -c ${extraopts} ${vflag}
    checkfail "coding generation"

    for i in `seq 0 $(expr ${n_coding} - 1)`
    do
        md5sum foo.c${i} > foo.c${i}.md5sum.1
    done

    if [ "${gf_type}" = 65537 ]
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
        if [ "${gf_type}" = 65537 ]
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
        if [ "${gf_type}" = 65537 ]
        then
            mv foo.c${j} foo.c${j}.1
            mv foo.c${j}.props foo.c${j}.props.1
        else
            mv foo.c${i} foo.c${i}.1
            mv foo.c${i}.props foo.c${i}.props.1
        fi
        j=`expr ${j} + 1`
    done

    echo "repairing"
    ${valgrind} ${bin} -e ${gf_type} -n ${n_data} -m ${n_coding} -p foo -r ${extraopts} ${vflag}
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
}

#do_test ./ec 65537 3 3 "" "" $*
do_test ./ec 8 3 3 "" "" $*
do_test ./ec 16 3 3 "" "" $*

#do_test ./ec 65537 3 3 "0 1" "0" $*
do_test ./ec 8 3 3 "0 1" "0" $*
do_test ./ec 16 3 3 "0 1" "0" $*

do_test ./ec 8 3 5 "0 1" "0" $*
do_test ./ec 16 3 5 "0 1" "0" $*
 
do_test ./ec 8 3 3 "1 2" "2" $*
do_test ./ec 16 3 3 "1 2" "2" $*

do_test ./ec 8 9 3 "1 2" "2" $*
do_test ./ec 16 9 3 "1 2" "2" $*

do_test ./ec 8 9 3 "2 3" "2" $*
do_test ./ec 16 9 3 "2 3" "2" $*

do_test ./ec 8 9 5 "2 3 4" "2 3" $*
do_test ./ec 16 9 5 "2 3 4" "2 3" $*

do_test ./ec 8 9 5 "1 3 5" "1 3" $*
do_test ./ec 16 9 5 "1 3 5" "1 3" $*

do_test ./ec 8 9 5 "1 3 5 7 8" "" $*
do_test ./ec 16 9 5 "1 3 5 7 8" "" $*

do_test ./ec 8 9 5 "" "0 1 2 3 4" $*
do_test ./ec 16 9 5 "" "0 1 2 3 4" $*
