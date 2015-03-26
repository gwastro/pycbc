#!/bin/bash

source profiling_utils.sh

declare -A args
parse_args $*
verify_args

source ${HOME}/pycbc_installations/${args["pycbc-version"]}/etc/pycbc-user-env.sh 

schemes=${args["processing-schemes"]}
tag=${args["tag"]}
iters=${args["iterations"]}
profile=${args["profile"]}

if [ "${iters}" == "" ]
then
	iters=1
fi

iter_count=0


outdir=${args["outdir"]}

if [ "${outdir}" == "" ]
then
	run_num=0

	while [ -e ${args["pycbc-version"]}_${tag}_${run_num} ]
	do
		run_num=$(( $run_num + 1 ))
	done

	outdir=${args["pycbc-version"]}_${tag}_${run_num}
fi

mkdir -p ${outdir}

while [ $iter_count -lt $iters ]
do

	for data in clean loud grumbly
	do
		run_tests $data "${schemes}" ${outdir}/`hostname`_${data}_iter${iter_count}.dat ${profile}
	done

	iter_count=$(( $iter_count + 1 ))
done

exit 0
