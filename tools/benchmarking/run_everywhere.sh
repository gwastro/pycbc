#!/bin/bash

ini_file=$1

source profiling_utils.sh

declare -A args
parse_args $*
verify_args

tag=${args["tag"]}
run_num=0

while [ -e ${args["pycbc-version"]}_${tag}_${run_num} ]
do
	run_num=$(( $run_num + 1 ))
done

outdir=${args["pycbc-version"]}_${tag}_${run_num}

mkdir -p ${outdir}



for node in `cat nodes`
do
	if [ "${ini_file}" == "" ]
	then
		run_file=${node}.ini
	else
		run_file=${ini_file}
	fi

	ssh ${node} "cd `pwd`; nohup ./pycbc_benchmark.sh ${run_file} outdir=${outdir} >${node}.out 2>${node}.err < /dev/null &"
done
