#!/bin/bash

# $? gets set to last nonzero exit status in pipe
set -o pipefail
# set -e

path="../bin" # for BOINC
if [ ".$1" = ".-p" ]; then
    path="$2"
    shift 2
fi

vdebug=3
if [ ".$1" = ".-d" ]; then
    vdebug="$2"
    shift 2
fi
echo "[`date`] $0 $*" >&2

# convert files
for f in $*; do
    if test -r "$f.txt"; then
        :
    else
        echo "[`date`] converting '$f'" >&2
        "$path/CBC1H5tool.py" "$f" | sort -g -k1,2 -k5,6 -k8 > "$f.txt"
        if [ $? -ne 0 ]; then
            echo "[`date`] $0 ERROR: CBC1H5tool.py $d '$f' failed" >&2
            rm -f "$f.txt"
            exit 1
        fi
    fi
    fs="$fs $f.txt"
done

# run validator_test on these
echo "[`date`] ../bin/validator_test_CBC1 -d $vdebug $fs" >&2
if "$path/validator_test_CBC1" -d $vdebug $fs; then
    :
else
    rm -f $fs
    exit 1
fi
